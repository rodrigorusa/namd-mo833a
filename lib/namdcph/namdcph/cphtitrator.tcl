# cphtitrator.tcl
#
#   This file provides the ::cphTitrator namespace, which effectively emulates
# an object containing all constant-pH information and routines relevant to
# proposing and performing Monte Carlo titrations.
#
#   If cphSystem and cphTitrator were true objects, this would essentially be
# an object composition. That is, a cphTitrator "contains" a cphSystem and
# yields information pertaining to it depending on a context (e.g. the pH).
#
package require Tcl 8.5

source [file join [file dirname [info script]] "cphtoppar.tcl"]
source [file join [file dirname [info script]] "numtcl.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphTitrator {
    #   The core information of the cphTitrator is the moveDict, which stores
    # all possible neMD/MC moves involving one or more titratable residues. The
    # keys are of the form:
    #
    # segresidname[/segresidname..[/segresidname]]
    #
    # in analogy to the keys of ::cphSystem::residueDict. Note that a "/" is
    # used as a syntactic divider between residues. Only _specific_ move
    # information is kept in moveDict. Another dict, defaultMoveParams is used
    # to store default information (similar to ::cphSystem::residueDefDict) and
    # is used to avoid saving copies of redundant information.
    #
    namespace import ::cphSystem::*

    variable moveDict [dict create]
    variable defaultMoveParams [dict create]
    variable maxAttempts

    namespace export cphSystem
    namespace export cphTitrator
}

# =============================================================================
# cphTitrator Interface function
# =============================================================================
# ::cphTitrator::cphTitrator
#
# This is the only exported function from the cphTitrator namespace and
# provides a complete interface. The first argument is always an action (e.g.
# "get" or "set") and the rest are variable arguments.
#
# action           description
# ---------------- -----------
# get              see cphTitratorGet 
# set              see cphTitratorSet
# build            see buildTitrator
# propose          select and propose a move, return switch information
# accumulateStats  accumulate mean acceptance rate for the given move
#
# Some notes about Tcl switch statements for non-guru types like myself:
#
#   These can appear a bit delicate to those who are not familiar with their
# quirks (compared to C). For one, badly placed comments can break them. Also
# note that a switch returns the expression it first matches and does NOT use
# break statements.
#
proc ::cphTitrator::cphTitrator {action args} {
    return [switch -nocase -- $action {
        get {
            cphTitratorGet {*}$args
        }
        set {
            cphTitratorSet {*}$args
        }
        build {
            buildTitrator {*}$args 
        }
        propose { ;# sample the move set and return the selected move info.
            proposeMove {*}$args
        }
        accumulateStats { ;# accumulate MC statistics
            accumulateAcceptanceRate {*}$args
        }
        default {
            abort "Invalid cphTitrator action $action"
        }
    }]
}

# ::cphTitrator::validateMoveLabel
#
# Check that a given move label is valid. That is:
#
# 1) is it in the correct segresidname[/segresidname ...] format?
# 2) are the component segresidnames valid?
#
# A non-zero error code and msg are returned for each condition. Zero is
# returned for a valid segresidname.
#
# NB! This should only be called _after_ buildSystem.
#
proc ::cphTitrator::validateMoveLabel {moveLabel} {
    set segresidnameList [split $moveLabel "/"]
    if {![llength $segresidnameList]} {
        print "Move labels should be of the form segresidname/segresidname ...]"
        return -1
    }
    foreach segresidname $segresidnameList { 
        if {[cphSystem validate $segresidname]} {
            print "Bad segresidname $segresidname in move label $moveLabel!"
            return -2
        }
    }
    return 0
}

# Pop the value from a dict with the corresponding key. If no key exists,
# return the the default value instead.
#
proc ::cphTitrator::dictPopOrDefault {myDict key {defaultValue {}}} {
    upvar 1 $myDict MyDict
    if {[dict exists $MyDict $key]} {
        set value [dict get $MyDict $key]
        dict unset MyDict $key
    }
    if {![info exists value]} {
        set value $defaultValue
    }
    return $value
}

# =============================================================================
# Proposal Routines
# =============================================================================
# ::cphTitrator::proposeMove
#
# Propose a move from all possible moves in the set. This is done by directly
# sampling the probability mass function of the weighted moves. Once this is
# done, test the inherent probability of that move given the imposed pH. If
# rejected, try again up to the given number of "maxAttempts". Always return
# the information necessary to perform a switch for the given proposal.
#
proc ::cphTitrator::proposeMove {pH} {
    variable ::cphTitrator::maxAttempts
    variable ::cphTitrator::moveDict
    set weights [cphTitrator get weight]

    set accept 0
    set attemptsThisCycle 0
    while {!$accept && $attemptsThisCycle < $maxAttempts} {
        incr attemptsThisCycle
        # Clear the previous trial if it was rejected.
        if {$attemptsThisCycle > 1} {
             cphSystem update $accept $segresidnameList
        }

        lassign [choice [cphTitrator get moveLabels] $weights] moveLabel
        set proposalCmd [cphTitrator get proposalCmd $moveLabel]
        set accept [eval $proposalCmd $pH]
        set segresidnameList [cphTitrator get segresidnameList $moveLabel]
    }
    set numsteps [cphTitrator get numsteps $moveLabel]
    return [list $accept $numsteps $segresidnameList $attemptsThisCycle]
}

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
proc ::cphTitrator::proposeResidueMove {segresidname pH} {
   lassign [cphSystem compute inherentWeights $pH $segresidname] weights states

   set pi [lindex $weights 0]
   lassign [choice [lrange $weights 1 end] [lrange $weights 1 end]] pj j
   set du [expr {log((1. - $pj) / (1. - $pi))}]
   set accept [metropolisAcceptance $du]
   if {$accept} {
       cphSystem set trialState $segresidname [lindex $states $j]
   }
   return $accept
}

# ::cphTitrator::proposeProtonTransferMove
#
# Propose a move involving two residues in which a proton is moved from one to
# the other. Return True if such a move was proposed, accepted, and stored.
#
# Note that this automatically fails if no proton can be transferred.
#
proc ::cphTitrator::proposeProtonTransferMove {moveLabel pH} {
    lassign [split $moveLabel "/"] segresidname1 segresidname2

    if {[cphSystem propose protonTransfer $segresidname1 $segresidname2]} {
        set accept 0
    } else {
        set du1 [cphSystem compute inherent $pH $segresidname1]
        set du2 [cphSystem compute inherent $pH $segresidname2]
        set accept [metropolisAcceptance [expr {$du1 + $du2}]]
    }
    return $accept
}

# =============================================================================
# "Constructor" Routines 
# =============================================================================
# ::cphTitrator::buildTitrator
#
# Construct the titrator given a set of MC move information.
#
proc ::cphTitrator::buildTitrator {moveInfo} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams
    variable ::cphTitrator::maxAttempts

    # 1) Pop off global and default move parameters.
    #
    set maxAttempts [dictPopOrDefault moveInfo maxProposalAttempts 0]
    if {$maxAttempts <= 0} {
        # This is a reasonable default for most systems.
        set maxAttempts [cphSystem get numresidues]
    }
    set defaultMoveParams [dict merge $defaultMoveParams\
                                      [dictPopOrDefault moveInfo default]]
    if {![dict exists $defaultMoveParams default weight]} {
        cphTitrator set weight default 1.0
    }
    if {![dict exists $defaultMoveParams default protonTransfer]} {
        cphTitrator set protonTransfer default 0
    }

    # 2) Validate all remaining keys - these better be move labels!
    #
    dict for {moveLabel data} $moveInfo {
        if {[validateMoveLabel $moveLabel]} {
            abort
        }
        dict for {attr value} $data {
            cphTitrator set $attr $moveLabel $value
        }
    }
    set moveDict [dict merge $moveDict $moveInfo]

    # The default move set is to titrate each residue independently.
    foreach segresidname [cphSystem get segresidnames] {
        set moveLabel $segresidname
        cphTitrator set proposalCmd $moveLabel\
                "proposeResidueMove $segresidname"
    }
    # Build proton transfer moves, if present.
    dict for {moveLabel data} $moveInfo {
        if {![dict exists $data protonTransfer]
            || ![dict get $data protonTransfer]} {
            continue
        }
        if {[llength [split $moveLabel "/"]] != 2} {
            abort "Proton transfer moves must have exactly 2 segresidnames!"
        }
        cphTitrator set proposalCmd $moveLabel\
                "proposeProtonTransferMove $moveLabel"
    }

    # Initialize statistics for any moves that are new.
    dict for {moveLabel data} $moveDict {
        if {![dict exists $data attempts]} {
            cphTitrator set successes $moveLabel 0
            cphTitrator set attempts $moveLabel 0
        }
    }

    return
}

# =============================================================================
# Getter Routines
# =============================================================================
# ::Titrator::cphTitratorGet
#
# Getter for move attributes, called as:
#
#   <attribute> [<moveLabel>]
#
# With some specialized exceptions, <attribute> is either the name of a system
# attribute (usually a list) or else a move attribute.
#
# system attributes  description
# -----------------  -----------
# moveLabels         list of all the moveLabels
# maxAttempts        max # of attempts at move proposals
#
# move attributes  description
# ---------------  -----------
# numsteps         number of steps in a switch after successful proposal
# weight           weight for consideration in move selection
# protonTransfer   if 2 residues are involved, can they proton transfer?
# successes        the number of successful _switches_ for this move
# attempts         the number of attempted _switches_ for this move
# proposalCmd      the command needed to set up trial states
# segresidnameList list of the residues involved (NB: this is just shorthand
#                  for [split $moveLabel "/"])
#
# specialized calls  description
# -----------------  -----------
# nodefaults         Return a dictionary of all current move information, but
#                    suppress entries where a move has the default value. This
#                    gives a minimal, but full, description of the move set.
#
# For now this is _much_ simpler than the analogous cphSystemGet. In the future
# it may be useful to generalize this though.
#
proc ::cphTitrator::cphTitratorGet {attr {moveLabel {}}} {
    variable ::cphTitrator::moveDict    
    variable ::cphTitrator::maxAttempts

    set getAll [expr {![llength $moveLabel]}]
    if {!$getAll && ![dict exists $moveDict $moveLabel]
        && ![string match -nocase $moveLabel default]} {
        abort "cphTitratorGet: Invalid moveLabel $moveLabel"
    }

    return [switch -nocase -- $attr {
        moveLabels {
            dict keys $moveDict
        }
        maxAttempts {
            expr {$maxAttempts}
        }
        numsteps -
        weight -
        protonTransfer -
        successes -
        attempts -
        proposalCmd {
            if {$getAll} {
                getAllMoveAttr $attr 
            } else {
                getMoveAttr $attr $moveLabel
            }
        }
        segresidnameList {
            if {$getAll} {
                getAllSegresidnameLists
            } else {
                split $moveLabel "/"
            }
        }
        nodefaults {
            if {$getAll} {
                cannotGetAll $attr
            } else {
                getMoveNoDefaults $moveLabel
            }
        }
        default {
            abort "cphTitratorGet: Invalid attribute $attr"
        }
    }]
}

proc ::cphTitrator::cannotGetAll {attr} {
    abort "cphTitratorGet: Cannot get all $attr - must select a moveLabel"
    return -1
}

# ---------------------------
# Getters for move attributes
# ---------------------------
# Return the default value if no specific value was set.
proc ::cphTitrator::getMoveAttr {attr moveLabel} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams

    if {[dict exists $moveDict $moveLabel $attr]} {
        return [dict get $moveDict $moveLabel $attr]
    } elseif {[dict exists $defaultMoveParams $attr]} {
        return [dict get $defaultMoveParams $attr]
    }
    abort "cphTitratorGet: Error getting $attr for move $moveLabel"
    return -1
}

proc ::cphTitrator::getAllMoveAttr {attr} {
    set retList [list]
    foreach moveLabel [cphTitrator get moveLabels] {
        lappend retList [cphTitrator get $attr $moveLabel]
    }
    return $retList
}

proc ::cphTitrator::getAllSegresidnameLists {} {
    set retList [list]
    foreach moveLabel [cphTitrator get moveLabels] {
        lappend retList [split $moveLabel "/"]
    }
    return $retList
}

# -----------------------------
# Special getters for archiving
# -----------------------------
# Return the dict for a given move label, but strip out default values.
proc ::cphTitrator::getMoveNoDefaults {moveLabel} {
    set defaultNumsteps [cphTitrator get numsteps default]
    set defaultWeight [cphTitrator get weight default]
    set defaultPT [cphTitrator get protonTransfer default]
    set numsteps [cphTitrator get numsteps $moveLabel]
    set weight [cphTitrator get weight $moveLabel]
    set PT [cphTitrator get protonTransfer $moveLabel]

    set retDict [dict create]
    if {$numsteps != $defaultNumsteps} {
        dict set retDict numsteps $numsteps
    }
    if {$weight != $defaultWeight} {
        dict set retDict weight $weight
    }
    if {$PT != $defaultPT} {
        dict set retDict protonTransfer $PT
    }
    dict set retDict attempts [cphTitrator get attempts $moveLabel]
    dict set retDict successes [cphTitrator get successes $moveLabel]

    return $retDict
}

# =============================================================================
# Setter Routines
# =============================================================================
# ::cphTitrator::cphTitratorSet
#
# Setters for move attributes, called as:
#
#  <attribute> <moveLabel> <value>
#
# <attribute> is the name of a move attribute.
#
# move attributes  description
# ---------------  -----------
# numsteps         number of steps in a switch after successful proposal
# weight           weight for consideration in move selection
# protonTransfer   if 2 residues are involved, can they proton transfer?
# successes        the number of successful _switches_ for this move
# attempts         the number of attempted _switches_ for this move
# proposalCmd      the command needed to set up trial states
#
proc ::cphTitrator::cphTitratorSet {attr moveLabel value} {
    variable ::cphTitrator::moveDict
    variable ::cphTitrator::defaultMoveParams

    if {[string match -nocase $moveLabel default]} {
        set setDefault 1
    } else {
        if {[validateMoveLabel $moveLabel]} {
            abort
        }
        set setDefault 0
    }

    return [switch -nocase -- $attr {
        numsteps -
        successes -
        attempts -
        protonTransfer { ;# Require integer or boolean argument.
            set value [expr {int($value)}]
            if {$setDefault} {
                dict set defaultMoveParams $attr $value
            } else {
                dict set moveDict $moveLabel $attr $value
            }
            expr {$value}
        }
        weight { ;# Require float argument.
            set value [expr {1.*$value}]
            if {$setDefault} {
                dict set defaultMoveParams $attr $value
            } else {
                dict set moveDict $moveLabel $attr $value
            }
            expr {$value}
        }
        proposalCmd { ;# Argument is string or pure list.
            dict set moveDict $moveLabel $attr $value
            expr {$value}
        }
        default {
            abort "cphTitratorSet: Invalid attribute $attr"
            expr {-1}
        }
    }]
}

# =============================================================================
# Routines for tracking and reporting MC statistics
# =============================================================================
# ::cphTitrator::accumulateAcceptanceRate
#
# Accumulate statistics for the given move.
#
proc ::cphTitrator::accumulateAcceptanceRate {accept moveLabel} {
    # Alas, dict incr does not support nested keys.
    set successes [cphTitrator get successes $moveLabel]
    set attempts [cphTitrator get attempts $moveLabel]
    incr successes $accept
    incr attempts
    cphTitrator set successes $moveLabel $successes
    cphTitrator set attempts $moveLabel $attempts
    return
}

