# cphtitrator.tcl
#
#   This file provides the ::cphTitrator namespace, which effectively emulates
# an object containing all constant-pH information and routines relevant to
# proposing and performing Monte Carlo titrations.
#
#   If cphSystem and cphTitrator were true objects, this would essentially be
# an object composition. That is, a cphTitrator "contains" a cphSystem and
# yields information pertaining to it depending on a context (most often, the
# pH).
#
package require Tcl 8.5

source [file join [file dirname [info script]] "cphtoppar.tcl"]
source [file join [file dirname [info script]] "numtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphTitrator {
    namespace import ::cphSystem::*
    variable pH
    variable MoveSet [dict create]
    variable WeightSum 0.0
    variable maxAttempts 0

    namespace export cphSystem
    # These are defined here
    namespace export buildTitrator proposeMove getMoveSet\
           accumulateAcceptanceRate
}

# =============================================================================
# "Constructor" Routines 
# =============================================================================
proc ::cphTitrator::buildTitrator {systempH moveInfo} {
    variable ::cphTitrator::pH $systempH
    variable ::cphTitrator::maxAttempts
    variable ::cphTitrator::MoveSet
    set maxProposalAttempts [dict get $moveInfo maxProposalAttempts]
    dict unset moveInfo maxProposalAttempts
    if {$maxProposalAttempts > 0} {
        set maxAttempts [expr {int($maxProposalAttempts)}]
    } else {
        set maxAttempts [cphSystem get numresidues]
    }
    #   The default move set is to switch each residue independently (i.e. the
    # move label is just the segresidname for that residue). Settings from the
    # restart and config files are discarded as they get used so that only
    # unused options remain in the moveInfo dict when this is finished (also
    # the default options, which must also be discarded.
    #
    foreach segresidname [cphSystem get segresidnames] {
        dict set MoveSet $segresidname proposalCmd\
                "proposeResidueMove $segresidname"
        dict set MoveSet $segresidname segresidnameList $segresidname
        dict set MoveSet $segresidname weight\
                [dictPopOrDefault moveInfo $segresidname weight]
        dict set MoveSet $segresidname numsteps\
                [dictPopOrDefault moveInfo $segresidname numsteps]
        dict set MoveSet $segresidname successes 0
        dict set MoveSet $segresidname attempts 0
        if {[dict exists $moveInfo $segresidname]\
            && ![dict size [dict get $moveInfo $segresidname]]} {
            dict unset moveInfo $segresidname
        }
    }
    # Build proton transfer moves
    foreach moveLabel [dict get $moveInfo ptransfer] {
        lassign [split $moveLabel "/"] segresidname1 segresidname2
        dict set MoveSet $moveLabel proposalCmd\
                "proposeProtonTransferMove $moveLabel" 
        dict set MoveSet $moveLabel segresidnameList\
                [list $segresidname1 $segresidname2]
        dict set MoveSet $moveLabel weight\
                [dictPopOrDefault moveInfo $moveLabel weight]
        dict set MoveSet $moveLabel numsteps\
                [dictPopOrDefault moveInfo $moveLabel numsteps]
        dict set MoveSet $moveLabel successes 0
        dict set MoveSet $moveLabel attempts 0
        if {[dict exists $moveInfo $moveLabel]\
            && ![dict size [dict get $moveInfo $moveLabel]]} {
            dict unset moveInfo $moveLabel
        }
    }
    dict unset moveInfo ptransfer
    ###
    computeWeightSum
    return $moveInfo
}

# This is a kludgy shorthand for using dicts that contain setting values. If
# the dict contains the move label, then any specific settings for that move
# are assumed to be stored as a nested dict with that label as key. All values
# are removed after being returned. If no such label or setting is present,
# then a default setting is assumed to exist with just the setting as the key.
#
proc ::cphTitrator::dictPopOrDefault {settingsDict label setting} {
    upvar 1 $settingsDict SettingsDict
    if {[dict exists $SettingsDict $label $setting]} {
        set value [dict get $SettingsDict $label $setting]
        dict unset SettingsDict $label $setting
    } else {
       set value [dict get $SettingsDict default $setting] 
    }
    return $value
}

# ::cphTitrator::computeWeightSum
#
# Compute, store, and return the sum of weights for the MC move set.
#
proc ::cphTitrator::computeWeightSum {} {
    lassign [getMoveWeightsAndSum] weights tmp
    set ::cphTitrator::WeightSum [expr {lsum($weights)}]
    return
}

# =============================================================================
# Proposal Routines
# =============================================================================
# ::cphTitrator::proposeResidueMove
#
# Propose a move involving a single residue via Metropolized independence
# sampling. Return True if such a move was proposed, accepted, and stored.
#
# The core idea here is that we want to propose the most probable new state
# based on the estimated inherent pKa values. For two state systems, this is
# trivial, it's just the other state. For three or more states, this will
# naturally prefer tautomerization at extreme pKa values.
#
proc ::cphTitrator::proposeResidueMove {segresidname} {
   variable ::cphTitrator::pH
   lassign [cphSystem compute inherentWeights $pH $segresidname] weights states 
   # NB: The remaining weights normalize to pic = (1 - pi), NOT 1.
   #     If you want to see the weights as probabilities, each weight should
   #     be divided by pic (this is implicit inside weightedChoice).
   #
   set pic [expr {1. - [lindex $weights 0]}]
   set j [weightedChoice [lrange $weights 1 end] $pic]
   set pj [lindex $weights [expr {$j+1}]] ;# note the index frame due to lrange
   set du [expr {log((1. - $pj) / $pic)}]
   set accept [metropolisAcceptance $du]
   if {$accept} {
       cphSystem set trialState $segresidname [lindex $states $j]
   }
   return $accept
}

proc ::cphTitrator::proposeProtonTransferMove {moveLabel} {
    variable ::cphTitrator::pH
    lassign [split $moveLabel "/"] segresidname1 segresidname2
    set errCode\
        [::cphSystem::proposeProtonTransfer $segresidname1 $segresidname2]
    if {$errCode > 0} {
        set du1 [cphSystem compute inherent $pH $segresidname1]
        set du2 [cphSystem compute inherent $pH $segresidname2]
        set accept [metropolisAcceptance [expr {$du1 + $du2}]]
    } else {
        set accept 0
    }
    return $accept
}

proc ::cphTitrator::proposeMove {} {
    variable ::cphTitrator::maxAttempts
    variable ::cphTitrator::MoveSet
    lassign [getMoveWeightsAndSum] weights weightSum

    set accept 0
    set nattempts 0
    while {!$accept && $nattempts < $maxAttempts} {
        incr nattempts
        # Clear the previous trial if it was rejected.
        if {$nattempts > 1} {
             cphSystem update $accept $segresidnameList
        }
        set index [weightedChoice $weights $weightSum]
        set moveLabel [lindex [dict keys $MoveSet] $index]
        set proposalCmd [dict get $MoveSet $moveLabel proposalCmd]
        set accept [eval $proposalCmd]
        set numsteps [dict get $MoveSet $moveLabel numsteps]
        set segresidnameList [dict get $MoveSet $moveLabel segresidnameList]
    }
    return [list $accept $numsteps $segresidnameList $nattempts]
}

# =============================================================================
# Routines for tracking and reporting MC statistics
# =============================================================================
# ::cphTitrator::accumulateAcceptanceRate
#
# Accumulate statistics for the given move.
#
proc ::cphTitrator::accumulateAcceptanceRate {accept moveLabel} {
    variable ::cphTitrator::MoveSet
    set moveLabel [join $moveLabel "/"]
    # Alas, dict incr does not support nested keys.
    set successes [dict get $MoveSet $moveLabel successes]
    set attempts [dict get $MoveSet $moveLabel attempts]
    incr successes $accept
    incr attempts
    dict set MoveSet $moveLabel successes $successes
    dict set MoveSet $moveLabel attempts $attempts 
    return
}

# =============================================================================
# Titrator Getter Routines
#
# These all return a list of values describing each move in the titrator.
# =============================================================================
proc ::cphTitrator::getMoveSet {} {
    return $::cphTitrator::MoveSet
}

# ::cphTitrator::getMoveWeightsAndSum
#
# Return the list of MC move weights along with their (precomputed) sum.
#
proc ::cphTitrator::getMoveWeightsAndSum {} {
    variable ::cphTitrator::MoveSet
    variable ::cphTitrator::WeightSum
    set weights [list]
    dict for {moveLabel data} $MoveSet {
        lappend weights [dict get $data weight]
    }
    return [list $weights $WeightSum]
}

