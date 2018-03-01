# cphtoppar.tcl
#
#   This file provides the ::cphSystem namespace, which effectively emulates an
# object containing all constant-pH specific topology and parameter (toppar)
# information.
#
package require Tcl 8.5

source [file join [file dirname [info script]] "cphpsfgen.tcl"]
source [file join [file dirname [info script]] "namdmcmc.tcl"]

namespace eval ::cphSystem {
    #   The core information of the cphSystem is the residueDict, which stores
    # the unique and mutable information for each instance of a titratable
    # residue. The keys are "segresidnames" of the form:
    #
    # <segid>:<resid>:<resname>.
    #
    # Note the <resname> is a _titratable residue_ name and may or may not be
    # the same as a recognizable residue from the force field. For example,
    # amino acid termini are titratable residues. The residueDefDict stores
    # non-unique and immutable information for each residue definition. The
    # keys are residue names (e.g. ASP, HIS, etc.). Generally, when looking up
    # information for a specific residue, one simply needs look up the name in
    # residueDict and then look up the info in residueDefDict.
    #
    # Note that a <segid>:<resid> combination can correspond to _multiple_
    # residue names. This permits non-overlapping "sub-residues," such as amino
    # acid termini, which titrate independently of the sidechain.
    #
    variable residueDict [dict create]
    variable residueDefDict [dict create]

    namespace export cphSystem
}

# =============================================================================
# cphSystem Interface function
# =============================================================================
# ::cphSystem::cphSystem
#
# This is the only exported function from the cphSystem namespace and provides
# a complete interface to the cphSystem namespace. The first argument is always
# an action (e.g. "get" or "set") and the rest are variable arguments.
#
# action           description
# ---------------- -----------
# get              see cphSystemGet 
# set              see cphSystemSet
# build            see buildSystem
# initialize       see initializeSystem
# propose          propose a move of the given type
# compute          compute various quantities needed for MC
# update           accept/reject the proposed changes
# alchemifypsf     apply patches that make the system alchemical
# dealchemifypsf   remove patches so the system is no longer alchemical
# initializeState  set the initial state of a residue
#
proc ::cphSystem::cphSystem {action args} {
    if {[string match -nocase $action get]} {
        # Getters
        return [cphSystemGet {*}$args] 
    } elseif {[string match -nocase $action set]} {
        # Setters
        return [cphSystemSet {*}$args]
    } elseif {[string match -nocase $action build]} {
        # System building - determine system composition and definitions
        return [buildSystem {*}$args]
    } elseif {[string match -nocase $action initialize]} {
        # System initialization - establish initial states and parameters
        return [initializeSystem {*}$args]
    } elseif {[string match -nocase $action propose]} {
        # Proposal - Set new trial state - return 0 if none is found
        set type [lindex $args 0]
        set newArgs [lrange $args 1 end]
        if {[string match -nocase $type titration]} {
            return [proposeResidueTitration {*}$newArgs]
        } elseif {[string match -nocase $type tautomerization]} {
            return [proposeResidueTautomerization {*}$newArgs]
        } else {
            abort "Invalid proposal type $type"
        }
    } elseif {[string match -nocase $action compute]} {
        # Compute acceptance energy (inherent or switch correction)
        set type [lindex $args 0]
        set newArgs [lrange $args 1 end]
        if {[string match -nocase $type inherent]} {
            return [computeInherentAcceptance {*}$newArgs]
        } elseif {[string match -nocase $type inherentWeights]} {
            return [computeInherentNormedWeights {*}$newArgs]
        } elseif {[string match -nocase $type switch]} {
            return [computeSwitchAcceptance {*}$newArgs]
        } else {
            abort "Invalid energy type $type"
        }
    } elseif {[string match -nocase $action update]} {
       # Update states based on an MC result
       return [updateStates {*}$args]
    } elseif {[string match -nocase $action alchemifypsf]} {
       # Alchemify the PSF (in memory) based on trial states
       return [alchemifyResidue {*}$args]
    } elseif {[string match -nocase $action dealchemifypsf]} {
       # Dealchemify the PSF (in memory)
       return [dealchemifyResidue {*}$args]
    } elseif {[string match -nocase $action initializeState]} {
       # Assign a state
       set method [lindex $args 0]
       if {[string match -nocase $method random]} {
           return [randomizeState {*}[lrange $args 1 end]]
       } else {
           return [assignState {*}$args]
       }
    } else {
        abort "Invalid cphSystem action $action."
    }
}

# =============================================================================
# Proposal Routines
# =============================================================================
# proc ::cphSystem::updateStates
#
# Update the state of one or more residues with the given
# <segid>:<resid>:<resname> specifications based on acceptance/rejection of the
# trial states. Reset the trial states to a null value.
#
proc ::cphSystem::updateStates {accept segresidnameList} {
    if {$accept} {
        foreach segresidname $segresidnameList {
            cphSystem set state $segresidname\
                    [cphSystem get trialState $segresidname]
            cphSystem set trialState $segresidname {}
        }
    } else {
        foreach segresidname $segresidnameList {
            cphSystem set trialState $segresidname {}
        }
    }
    return 0
}

# ::cphSystem::proposeResidueTitration
#
# Propose a new trial state requiring a titration - i.e. a net change in the 
# number of protons.
#
# For consistency with tautomers (see below), this always returns true in order
# to confirm that a titration was in fact found.
#
proc ::cphSystem::proposeResidueTitration {segresidname} {
    variable ::cphSystem::resDefDict
    set possibleStates [cphSystem get trialStateList $segresidname]
    set resname [cphSystem get resname $segresidname]
    set numProtons [expr {lsum([cphSystem get occupancy $segresidname])}]
    while {true} {
        set state [choice $possibleStates]
        set occ [dict get $resDefDict $resname states $state]
        set numProtonsTrial [expr {lsum($occ)}]
        if {$numProtons != $numProtonsTrial} {
            cphSystem set trialState $segresidname $state
            return 1
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::proposeProtonTransfer
#
# Propose a new trial state requiring a proton transfer - i.e. a movement of a
# proton from residue to another.
#
# If both residues in the pair have/do not currently have a proton, then a
# transfer cannot occur and this returns false.
#
proc ::cphSystem::proposeProtonTransfer {segresidname1 segresidname2} {
    variable ::cphSystem::resDefDict
    set numProtons1 [expr {lsum([cphSystem get occupancy $segresidname1])}]
    set possibleStates1 [cphSystem get trialStateList $segresidname1]
    set resname1 [cphSystem get resname $segresidname1]
    set numProtons2 [expr {lsum([cphSystem get occupancy $segresidname2])}]
    set possibleStates2 [cphSystem get trialStateList $segresidname2]
    set resname2 [cphSystem get resname $segresidname2]
    set dn [expr {$numProtons1 - $numProtons2}]

    if {$dn == 0.0} {     
        # No proton transfer is possible.
        return 0
    } 
    # transfer from 1 to 2
    while {true} {
        set state1 [choice $possibleStates1]
        set occ1 [dict get $resDefDict $resname1 states $state1]
        set numProtonsTrial1 [expr {lsum($occ1)}]
        set state2 [choice $possibleStates2]
        set occ2 [dict get $resDefDict $resname2 states $state2]
        set numProtonsTrial2 [expr {lsum($occ2)}]
        if {$dn == [expr {$numProtonsTrial2 - $numProtonsTrial1}]} {
            cphSystem set trialState $segresidname1 $state1
            cphSystem set trialState $segresidname2 $state2
            return 1
        }
    }
    # This is an error and should never happen.
    return -1
}


# ::cphSystem::proposeResidueTautomerization
#
# Propose a new trial state requiring a tautomerization - i.e. the number of 
# protons remains unchanged. 
# 
# Unlike a state titration, a state tautomerization is not guaranteed to exist 
# for all residues. Return true if one is found, else return false.
#
# In order to ensure that a tautomer is found (if it exists), but also not
# cause an infinite loop when no tautomer exists, the number of state proposals
# is capped at maxAttempts. For a typical residue with a tautomer, there are
# three states, two of which can interconvert via tautomerization. The
# probability of selecting the tautomeric state k times in N trials is thus
# binomially distributed with p = (1 - p) = 0.5:
#
# P(k, N) = [N! / (k!(N-k)!)] 0.5^N
#
# The probability of picking that correct state at least once is thus:
#
# P(k>0, N) = 1 - P(0, N) = 1 - 0.5^N 
#
# or
#
# N = log[1 - P(k>0, N)] / log(0.5)
#
# For P(k>0, N) = 0.999, this gives N ~= 10 (the default).
#
proc ::cphSystem::proposeResidueTautomerization {segresidname {maxAttempts 10}} {
    variable ::cphSystem::resDefDict
    set possibleStates [cphSystem get trialStateList $segresidname]
    set resname [cphSystem get resname $segresidname]
    set numProtons [expr {lsum([cphSystem get occupancy $segresidname])}]
    while {true} {
        set state [choice $possibleStates]
        set occ [dict get $resDefDict $resname states $state]
        set numProtonsTrial [expr {lsum($occ)}]
        if {$numProtonsTrial == $numProtons} {
            cphSystem set trialState $segresidname $state
            return 1
        }
        incr attempts
        if {$attempts >= $maxAttempts} {
            # This probably implies that no tautomer exists. 
            return 0 
        }
    }
    # This is an error and should never happen.
    return -1
}

# ::cphSystem::computeInherentAcceptance
#
# Compute the (reduced) energy difference for a Monte Carlo move based on the 
# given segresidname, its current and trial state, and the given pH.
#
# The proposal energy is based exclusively on the "intrinsic" pKa and the
# change in the protonation vector. There are two cases: 1) tautomerizations
# and 2) titrations, the former of which is not a proper pKa as the move does
# not depend on pH (since no protons are entering or leaving the bath).
#
# case 1, n = n':
#
#     P(s --> s')
#     ----------- = 10^[sgn(s' - s) pKa_i(s, s')]
#     P(s' --> s)
#
# case 2, n != n':
#     
#     P(s --> s')
#     ----------- = 10^[sgn(n' - n) pKa_i(s, s') - (n' - n)pH]
#     P(s' --> s)
#
# where s and s' are the current and trial state indices, respectively, with
# number of protons (i.e. sum of the occupation vector elements) n and n',
# respectively. By convention, pKa_i(s, s') = pKa_i(s', s) and the antisymmetry
# of adding vs deleting protons is accounted for by the sgn function. Note
# that for tautomers, the sgn is computed by the _state index_, not the number
# of protons - this is also an arbitrary internal convention.
#
# The above ratios are correctly sampled by a Metropolis sampling of the form:
#
#     P(s --> s') = min{1, e^[-du(s, s')]},
#
# where du is either of the the exponents above times -ln(10).
#
proc ::cphSystem::computeInherentAcceptance {pH segresidname} {
    set l [cphSystem get occupancy $segresidname]
    set lp [cphSystem get trialOccupancy $segresidname]
    set dn [expr {lsum($lp) - lsum($l)}]
    set s [occupancy2Index $l]
    set sp [occupancy2Index $lp]
    set ssp [index2flatindex $s $sp]
    set pKai [lindex [cphSystem get pKaiPair $segresidname] $ssp]
    if {$dn == 0.0} {
        # tautomerization: "pKa" is positive in direction of lower index
        set sgn [expr {$sp > $s} ? 1.0 : -1.0]
    } else {
        # titration: pKa is positive in direction of fewer protons
        set sgn [expr {$dn > 0} ? 1.0 : -1.0]
    }
    return [expr {-$::LN10*($sgn*$pKai - $dn*$pH)}] 
}

# ::cphSystem::computeInherentNormedWeights
#
# Compute the normalized inherent pKa weights of all states.
#
# This is a bit tricky, bc only the relative unnormed weights are known. This
# is solved by taking all weights relative to the current state and assigning
# that state an unnormed weight of one.
#
# NB: The order of states is randomly shuffled here, because the current state
#   is omitted from the stateList and given the index 0 in the list of weights
#   (hence the weight list is one element longer!).
#
proc ::cphSystem::computeInherentNormedWeights {pH segresidname} {
    set l [cphSystem get occupancy $segresidname]
    set s [occupancy2Index $l]
    set resname [cphSystem get resname $segresidname]
    set stateList [cphSystem get trialStateList $segresidname]
    # We implicitly reference against the current state, so its unnormed weight
    # is exactly one.
    set logQs [list 1.0]
    set logQMax 1.0
    foreach state $stateList {
        set lp [state2Occupancy $resname $state]
        set dn [expr {lsum($lp) - lsum($l)}]
        set sp [occupancy2Index $lp]
        set ssp [index2flatindex $s $sp]
        set pKai [lindex [cphSystem get pKaiPair $segresidname] $ssp]
        if {$dn == 0.0} {
            # tautomerization: "pKa" is positive in direction of lower index
            set sgn [expr {$sp > $s} ? 1.0 : -1.0]
        } else {
            # titration: pKa is positive in direction of fewer protons
            set sgn [expr {$dn > 0} ? 1.0 : -1.0]
        }
        lappend logQs [expr {$::LN10*($sgn*$pKai - $dn*$pH)}]
        if {[lindex $logQs end] > $logQMax} {
            set logQMax [lindex $logQs end]
        }
    }
    return [list [normalizeLogWeights $logQs $logQMax] $stateList]
}

# ::cphSystem::computeSwitchAcceptance
#
# Compute the pairwise energy difference used to shift the work value in the
# Monte Carlo acceptance criteria.
#
#   This is only the configuration _independent_ portion of that difference 
# (i.e. it only depends on the current and trial protonation states). The
# configuration _dependent_ portion is the total work W applied to effect the
# switch. As for the proposal energy, the sign of the correction depends on
# whether this is tautomerization (1) or titration (2) move:
#
# dE(s, s') = c*{-dG(s, s') + kT ln(10) [pKa(s, s') - pKa_i(s, s')]}
# 
# case 1, n = n':  c = sgn(s' - s)
#
# case 2, n != n': c = sgn(n' - n)
#
# Here T is the _reference_ temperature at which pKa is measured and dG is
# computed, not the temperature of the simulation. See 
# computeInherentAcceptance for definition of the remaining notation.
#
proc ::cphSystem::computeSwitchAcceptance {segresidnameList} {
    set du 0.0
    foreach segresidname $segresidnameList {
        set l [cphSystem get occupancy $segresidname]
        set lp [cphSystem get trialOccupancy $segresidname]
        set dn [expr {lsum($lp) - lsum($l)}]
        set s [occupancy2Index $l]
        set sp [occupancy2Index $lp]
        set ssp [index2flatindex $s $sp]
        set dG [lindex [cphSystem get dGPair $segresidname] $ssp]
        set pKa [lindex [cphSystem get pKaPair $segresidname] $ssp]
        set pKai [lindex [cphSystem get pKaiPair $segresidname] $ssp]
        if {$dn == 0.0} {
            # tautomerization: "pKa" is positive in direction of lower index
            set sgn [expr {$sp > $s} ? 1.0 : -1.0]
        } else {
            # titration: pKa is positive in direction of fewer protons
            set sgn [expr {$dn > 0} ? 1.0 : -1.0]
        }
        set kT [expr {$::BOLTZMANN*[cphSystem get Tref $segresidname]}]
        set du [expr {$du + $sgn*($dG - $kT*$::LN10*($pKa - $pKai))}]
    }
    return $du
}

# ---------------
# Psfgen Routines
# ---------------
# ::cphSystem::alchemifyResidue
#
# Apply a trial alchemical patch to a residue.
#
proc ::cphSystem::alchemifyResidue {segresidname frcCons temp {buildH false}} {
    lassign [cphSystem get alchAtomLists $segresidname] l0atoms l1atoms
    set patch [cphSystem get hybridPatch $segresidname]
    lassign [split $segresidname ":"] segid resid
    # NB: alchPatch uses psfgen style selections!
    alchPatch $patch "$segid:$resid" $l0atoms $l1atoms $frcCons $temp $buildH
    return 0
}

# ::cphSystem::dealchemifyResidue
#
# Remove an alchemical patch from a residue.
#
proc ::cphSystem::dealchemifyResidue {segresidname} {
    lassign [cphSystem get alchAtomLists $segresidname] l0atoms l1atoms
    lassign [split $segresidname ":"] segid resid 
    # NB: alchUnpatch uses psfgen style selections!
    alchUnpatch "$segid:$resid" $l0atoms $l1atoms
    return 0
}

# =============================================================================
# "Constructor" Routines
# =============================================================================
# ::cphSystem::initializeSystem
#
# Initialize the state of the system using:
#
# 1) input data from the user (usually from a restart file)
#
# AND/OR
#
# 2) the ensemble information (i.e. the pH and temperature)
#
proc ::cphSystem::initializeSystem {pH temperature buildH stateInfo} {
    variable ::cphSystem::resDict
    dict for {segresidname resData} $resDict {
        # Assign inherent pKa values.
        if {[dict exists $stateInfo $segresidname pKai]} {
            cphSystem set pKai $segresidname\
                    [dict get $stateInfo $segresidname pKai]
        } else { ;# Default to reference pKa.
            cphSystem set pKai $segresidname [cphSystem get pKa $segresidname]
        }

        # Assign states.
        if {[dict exists $stateInfo $segresidname state]} {
            set state [dict get $stateInfo $segresidname state]
            cphSystem initializeState $segresidname $state 
        } else { ;# default randomization based on pKai and pH
            cphSystem initializeState random $segresidname $pH
        }
    }
    # Final pass - Apply the patches
    foreach segresidname [cphSystem get segresidnames] {
        lassign [split $segresidname ":"] segid resid
        patch [cphSystem get statePatch $segresidname] "$segid:$resid"
    }
    guesscoord
    foreach segresidname [cphSystem get segresidnames] {
        cphSystem alchemifypsf $segresidname 0.0 $temperature $buildH
        cphSystem dealchemifypsf $segresidname
        cphSystem update 1 $segresidname
    }
    regenerate angles dihedrals
    # NB - These changes are only reflected in _memory_ for the cphSystem and
    # psfgen. Nothing has happened to the NAMD PSF/PDB files.
    return
}

# ::cphSystem::buildSystem
#
# Build the residue definitions and residue objects for the system based on 
# the NAMD inputs.
#
proc ::cphSystem::buildSystem {resDefs resAliases segresExcls} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict [checkResDefs $resDefs]

    # Check for special sub-type residues.
    # In general, THESE WILL NOT BE MATCHED BY NAME, but rather are sub-typed
    # within one or more other residue definitions. Therefore we keep a
    # separate dict of residue names that have a sub-type (these may not even
    # be titratable otherwise, such as terminal amino acids). The keys are
    # the actual residue names and the values are a list of the possible
    # sub-types. However, even if a sub-typed residue is present, that doesn't
    # mean the sub-residue exists. There are two (optional) additional checks
    # in the sub-residue definition:
    #
    # subatoms - a list of atom names that MUST exist
    # notsubatoms - a list of atom names that MUST NOT exist
    #
    # For example, alanine is not titratable, but may contain a titratable
    # C-terminus. The terminal oxygen atoms must exist, but are only titratable
    # if they are not blocked by amination or methylation.
    #
    set subresDefDict [dict create]
    dict for {subresname resDef} $resDefDict {
        if {![dict exists $resDef subtypeof]} continue
        foreach subtypedRes [dict get $resDef subtypeof] {
            if {![dict exists $subresDefDict $subtypedRes]} {
                dict set subresDefDict $subtypedRes [list]
            }
            dict lappend subresDefDict $subtypedRes $subresname
        }
    }

    # Read in whatever files were specified to NAMD.
    set Args [list [structure] pdb [coordinates]]
    if {[isset binCoordinates]} {
        lappend Args namdbin [binCoordinates]
        if {[isset binVelocities]} {
            lappend Args velnamdbin [binVelocities]
        }
    }
    resetpsf
    readpsf {*}$Args

    set definedResidues [dict keys $resDefDict]
    set definedSubResidues [dict keys $subresDefDict]
    foreach segid [segment segids] {
        foreach resid [segment resids $segid] {
            set resname [segment residue $segid $resid]
            set segresid [format "%s:%s" $segid $resid]
            set segresidname [format "%s:%s" $segresid $resname]
            # Perform aliasing to fix residues with multiple names.
            dict for {realName altNames} $resAliases {
                if {[lsearch -nocase $altNames $resname] < 0} {
                    continue
                }
                print "aliasing $segresidname to $realName"
                psfset resname $segid $resid $realName
                set resname $realName
            }
            set resIsDefined 0
            if {[lsearch -nocase $definedResidues $resname] >= 0} {
                set resIsDefined 1
            }
            set resIsSubtyped 0
            if {[lsearch -nocase $definedSubResidues $resname] >= 0} {
                set resIsSubtyped 1
            }
            set resIsExcluded 0
            if {[lsearch -nocase $segresExcls $segresidname] >= 0} {
                set resIsExcluded 1
            }
            # Break here if nothing to do.
            if {(!$resIsDefined && !$resIsSubtyped) || $resIsExcluded} {
                continue
            }

            if {$resIsDefined} { 
                dict set resDict $segresidname [dict create]
                dict set resDict $segresidname resname $resname
            }

            if {!$resIsSubtyped} continue

            # The conditions for subtyping may still not be met...
            set atoms [segment atoms $segid $resid]
            foreach subresname [dict get $subresDefDict $resname] {
                set subatoms [list]
                if {[dict exists $resDefDict $subresname subatoms]} {
                    set subatoms [dict get $resDefDict $subresname subatoms]
                }
                set notsubatoms [list]
                if {[dict exists $resDefDict $subresname notsubatoms]} {
                    set notsubatoms\
                            [dict get $resDefDict $subresname notsubatoms]
                }

                set allSubatomsExist 1
                foreach atom $subatoms {
                    if {[lsearch -nocase $atoms $atom] < 0} {
                        set allSubatomsExist 0
                        break
                    }
                }
                set allNotSubatomsDoNotExist 1
                foreach atom $notsubatoms { 
                    if {[lsearch -nocase $atoms $atom] >= 0} {
                        set allNotSubatomsDoNotExist 0
                        break
                    }
                }
                if {$allSubatomsExist && $allNotSubatomsDoNotExist} {
                    set segresidname [format "%s:%s" $segresid $subresname]
                    dict set resDict $segresidname [dict create]
                    dict set resDict $segresidname resname $subresname
                }
            }

        }
    }
    return $resDict
}

# ::cphSystem::checkResDefs
#
# Check that a dictionary of residue definitions is valid and consistent with
# the topology (RTF) files currently loaded by psfgen.
#
proc ::cphSystem::checkResDefs {resDefs} {
    dict for {resname resDef} $resDefs {
        # Check that all required fields are present.
        foreach dtype [list dG pKa states l0atoms] {
            if {![dict exists $resDef $dtype]} {
                abort "No $dtype entry for residue $resname!"
            }
        }
        set pKa [dict get $resDef pKa]
        set dG [dict get $resDef dG]
        set states [dict get $resDef states]
        set l0atoms [dict get $resDef l0atoms]
        # If not specified, use lazy method of appending a "1". If an l1atoms
        # field is present, check that it matches l0atoms.
        if {![dict exists $resDef l1atoms]} {
            set l1atoms [list]
            foreach atom $l0atoms { 
                lappend l1atoms [format "%s1" $atom]
            }
            dict set resDefs $resname l1atoms $l1atoms
        } else {
            set l1atoms [dict get $resDefs $resname l1atoms]
        }

        if {[llength $l0atoms] != [llength $l1atoms]} { 
            abort "Mismatch in atom definitions for residue $resname"
        }
        # Check that the definitions for pKa and dG are consistent.
        if {[llength $pKa] != [llength $dG]} {
            abort "Mismatch in dG/pKa definitions for residue $resname."
        }
        # Check that the state and occupancy definitions are consistent.
        set numSites [llength [dict get $states [lindex $states 0]]]
        set maxStates [expr {int(pow(2, $numSites))}]
        if {[dict size $states] > $maxStates} {
            abort "Too many states defined for residue $resname!"
        }
        dict for {state occ} $states {
            if {[llength $occ] != $numSites} {
                abort "Bad occupancy definition for $resname state $state."
            }
            # Check that the RTFs define two patches for each state:
            # 1) state patches, with the naming convention: <resname><state>
            # 2) hybrid patches, with the naming convention: <resname>H<state>
            #
            set statePatch [format "%s%s" $resname $state]
            set hybridPatch [format "%sH%s" $resname $state]
            if {[lsearch -nocase [topology patches] $statePatch] < 0} {
                abort "No patch definition in RTFs for $statePatch!"
            }
            if {[lsearch -nocase [topology patches] $hybridPatch] < 0} {
                abort "No patch definition in RTFs for $hybridPatch!"
            }
        }
        # Build the pairwise parameters for this definition.
        dict set resDefs $resname dGPair [resDef2Matrix $resDef $dG]
        dict set resDefs $resname pKaPair [resDef2Matrix $resDef $pKa]
    }
    return $resDefs
}

# ::cphSystem::assignState
#
# Assign a protonation state by fiat.
#
#   This is a bit of a hack. The state is assigned and then an arbitrary trial
# state is chosen by rejectionless MC. The states are then swapped and the move
# is accepted as if this had happened in reverse.
#
# Arguments:
# ----------
# segresidname : string
#   residue specification as "<segid>:<resid>" - this is the same syntax as for
#   the regular psfgen patch command.
# state : string
#   state to assign 
#
# Returns:
# --------
# None
#
proc ::cphSystem::assignState {segresidname state} {
    cphSystem set state $segresidname $state
    cphSystem propose titration $segresidname
    cphSystem set state $segresidname [cphSystem get trialState $segresidname]
    cphSystem set trialState $segresidname $state
    return 0
}

# ::cphSystem::randomizeState
#
# Assign a protonation state by randomly choosing a state and performing MC 
# moves based on the pH until a new state is accepted.
#
#   This is a bit sneaky. Implementing a proper independence sampling of all of
# the states seems overkill. Rather, an originating state is chosen uniformly 
# from all possible states and the state to be assigned is chosen by a pairwise
# pH/pKai based Metropolis criterion. If this is rejected, then a new 
# originating state is chosen. This should be equivalent to independence
# sampling for many such trials. 
#
# Arguments:
# ----------
# segresidname : string
#   residue specification as "<segid>:<resid>" - this is the same syntax as for
#   the regular psfgen patch command.
# pH : float
#    pH value at which to assign protonation states
#
# Returns:
# --------
# None 
#
proc ::cphSystem::randomizeState {segresidname pH} {
    while {true} {
        set states [cphSystem get stateList $segresidname]
        cphSystem set state $segresidname [choice $states]
        cphSystem propose titration $segresidname
        set du [cphSystem compute inherent $pH $segresidname]
        if {[metropolisAcceptance $du]} { 
            return 0
        }
    }
    return -1
}

# =============================================================================
# Getter Routines
# =============================================================================
# ::cphSystem::cphSystemGet 
#
# Getters for system and residue attributes, called as:
#
#   <attribute> [<segresidname> [<args>]]
#
# <attribute> is the name of either a system attribute (segresidname selections
# are invalid) or else a residue attribute. For SOME residue attributes, not
# specifying a segresidname will return a list for all residues. Some
# attributes also require some additional arguments (see below)
#
# system attributes  description
# -----------------  -----------
# segresidnames      list of all segresidnames
# numresidues        number of residues in the system
# resdefs            list of defined resnames
# numdefs            number of defined resnames
#
# residue attributes description
# ------------------ -----------
# resname            residue name
# state              current state
# trialState         proposed trial state
# pKai               minimal pKai list for this residue
# dG                 minimal dG list for this residue 
# pKa                minimal pKa list for this residue
# Tref               reference temperature for pKa
# occupancy          occupancy vector for the current state
# trialOccupancy     occupancy vector for the trial state
# stateList          all possible states
# trialStateList     all possible trial states (not the current state)
# dGPair             pair dG for current/trial states
# pKaPair            pair pKa for current/trial states
# pKaiPair           pair pKai for current/trial states 
# statePatch         patch for the current state
# hybridPatch        alchemical patch for the trial state
# alchAtomLists*     lists of alchemical atoms at 0/1
# alchBonds*^        extraBonds entries
#
# * - segresidname selection is required
# ^ - requires additional arguments
#
proc ::cphSystem::cphSystemGet {attr {segresidname {}} args} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict

    set getAll [expr {![llength $segresidname]}]
    if {!$getAll && ![dict exists $resDict $segresidname]} {
        abort "cphSystemGet: Invalid segresidname $segresidname"
    }
    # System attributes - any selection is invalid.
    #
    if {[string match -nocase $attr segresidnames]} {
        return [dict keys $resDict]
    } elseif {[string match -nocase $attr numresidues]} {
        return [dict size $resDict]
    } elseif {[string match -nocase $attr resdefs]} {
        return [dict keys $resDefDict]
    } elseif {[string match -nocase $attr numdefs]} {
        return [dict size $resDefDict]
    }
    # Residue attributes - some of these need to be specially computed.
    #
    if {[lsearch -nocase {resname state trialState pKai} $attr] > -1} {
        if {$getAll} {
            return [getAllResAttr $attr]
        } else {
            return [getResAttr $attr $segresidname]
        }
    } 
    if {[lsearch -nocase {dG pKa Tref} $attr] > -1} {
        if {$getAll} {
            return [getAllResDefAttr $attr]
        } else {
            return [getResDefAttr $attr $segresidname]
        }
    } 
    if {[string match -nocase $attr occupancy]} {
        if {$getAll} {
            return [getAllOccupancy]
        } else {
            return [getOccupancy $segresidname]
        }
    } 
    if {[string match -nocase $attr trialOccupancy]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            return [getTrialOccupancy $segresidname]
        }
    } 
    if {[string match -nocase $attr stateList]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            set resname [getResAttr resname $segresidname]
            return [dict keys [dict get $resDefDict $resname states]] 
        }
    }
    if {[string match -nocase $attr trialStateList]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            set resname [getResAttr resname $segresidname]
            set state [getResAttr state $segresidname]
            set states [dict keys [dict get $resDefDict $resname states]]
            return [lsearch -all -inline -not -nocase $states $state]
        }
    }
    if {[lsearch -nocase {dGPair pKaPair pKaiPair} $attr] > -1} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            if {[lsearch -nocase {dGPair pKaPair} $attr] > -1} {
                return [getResDefAttr $attr $segresidname]
            } else {
                return [getResAttr $attr $segresidname]
            }
        }
    }
    if {[string match -nocase $attr statePatch]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            set resname [getResAttr resname $segresidname]
            set state [getResAttr state $segresidname]
            return [format "%s%s" $resname $state]
        }
    }
    if {[string match -nocase $attr hybridPatch]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            set resname [getResAttr resname $segresidname]
            set trialState [getResAttr trialState $segresidname]
            return [format "%sH%s" $resname $trialState]
        }
    }
    if {[string match -nocase $attr alchAtomLists]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            return [list [getResDefAttr l0atoms $segresidname]\
                         [getResDefAttr l1atoms $segresidname]]
        }
    }
    if {[string match -nocase $attr alchBonds]} {
        if {$getAll} {
            cannotGetAll $attr
        } else {
            lassign $args k
            set bondEntries [list]
            foreach l0atom [getResDefAttr l0atoms $segresidname]\
                    l1atom [getResDefAttr l1atoms $segresidname] {
                lassign [split $segresidname ":"] segid resid
                # Note that atomid indices start at one, not zero!
                set i [expr {[segment atomid $segid $resid $l0atom] - 1}]
                set j [expr {[segment atomid $segid $resid $l1atom] - 1}]
                lappend bondEntries [format "bond %d %d %f %f" $i $j $k 0]
            }
            return [join $bondEntries "\n"]
        }
    }
    abort "cphSystemGet: Invalid attribute $attr"
}

proc ::cphSystem::cannotGetAll {attr} {
    abort "cphSystemGet: Cannot get all $attr - must select a segresidname"
}

# ------------------------------
# Getters for residue attributes
# ------------------------------
proc ::cphSystem::getResAttr {attr segresidname} {
    variable ::cphSystem::resDict
    return [dict get $resDict $segresidname $attr]
}

proc ::cphSystem::getAllResAttr {attr} {
    variable ::cphSystem::resDict
    set retList [list]
    foreach resData [dict values $resDict] {
        lappend retList [dict get $resData $attr]
    }
    return $retList
}

# -----------------------------------------
# Getters for residue definition attributes
# -----------------------------------------
proc ::cphSystem::getResDefAttr {attr segresidname} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set resname [dict get $resDict $segresidname resname]
    return [dict get $resDefDict $resname $attr]
}

proc ::cphSystem::getAllResDefAttr {attr} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set retList [list]
    foreach resData [dict values $resDict] { 
        set resname [dict get $resData resname]
        lappend retList [dict get $resDefDict $resname $attr]
    }
    return $retList
}

# -----------------------------------------------
# Special getters for more complicated operations
# -----------------------------------------------
proc ::cphSystem::state2Occupancy {resname state} {
    return [dict get $::cphSystem::resDefDict $resname states $state]
}

# NB: When returning all occupancies, only the current state can be used and
#   the list is flattened.
proc ::cphSystem::getOccupancy {segresidname} {
    variable ::cphSystem::resDict
    set resname [dict get $resDict $segresidname resname]
    set state [dict get $resDict $segresidname state]
    return [state2Occupancy $resname $state]
}

proc ::cphSystem::getTrialOccupancy {segresidname} {
    variable ::cphSystem::resDict
    set resname [dict get $resDict $segresidname resname]
    set state [dict get $resDict $segresidname trialState]
    return [state2Occupancy $resname $state]
}

proc ::cphSystem::getAllOccupancy {} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    set retList [list]
    foreach resData [dict values $resDict] {
        set resname [dict get $resData resname]
        set state [dict get $resData state]
        set retList [list {*}$retList {*}[state2Occupancy $resname $state]]
    }
    return $retList
}

# ----------------------
# Other helper functions
# ----------------------
# Given an occupancy, get the state index.
# This is essentially a left to right conversion from binary to decimal.
# That is, 
#   index = l0*2**0 + l1*2**1 + l2*2**2 + ...
#
proc ::cphSystem::occupancy2Index {occupancy} {
    set index 0
    for {set i 0} {$i < [llength $occupancy]} {incr i} {
        incr index [expr {int([lindex $occupancy $i]*pow(2, $i))}]
    }
    return $index
}

# Map a pair index to the flat off-diagonal index.
#
# This maps to the _lower_ off-diagonal counting left-right, top-bottom.
#
# Ex. n = 4, 2,0 --> 1 and 3,2 --> 5
#
# |x x x x|
# |0 x x x|
# |1 2 x x|
# |3 4 5 x|
#
# The additional assumption here is that element i,j is the same as j,i.
# Therefore this internally swaps the two indices such that i > j _always_.
#
proc ::cphSystem::index2flatindex {i j} {
#    lassign [list [expr {max($i, $j)}] [expr {min($i, $j)}]] I J
    if {$i > $j} {
        set I $i
        set J $j
    } elseif {$i < $j} {
        set I $j
        set J $i
    } else {
        abort "invalid transition from state $i to state $i"
    }
    return [expr {($I*($I - 1) + 2*$J) / 2}]
}

# =============================================================================
# Setter Routines
# =============================================================================
# ::cphSystem::cphSystemSet
#
# Setters for residue attributes, called as:
#
#   <attribute> <segresidname> <value>
#
# <attribute> is the name of a residue attribute.
#
# residue attributes description
# ------------------ -----------
# state              current state
# trialState         proposed trial state
# pKai               minimal pKai list for this residue
#
proc ::cphSystem::cphSystemSet {attr segresidname value} {
    variable ::cphSystem::resDict
    variable ::cphSystem::resDefDict
    if {![dict exists $resDict $segresidname]} {
        abort "cphSystemSet: Invalid segresidname $segresidname"
    }
    if {[lsearch -nocase {state trialState} $attr] > -1} {
        set states [cphSystem get stateList $segresidname]
        if {[lsearch -nocase $states $value] < 0} {
            if {[string match -nocase $attr trialState] && ![llength $value]} {
            
            } else {
                abort "Invalid state assignment $value for residue $segresidname"
            }
        }
        dict set resDict $segresidname $attr $value
        return $value
    }
    if {[string match $attr pKai]} {
        set resname [cphSystem get resname $segresidname]
        set resDef [dict get $resDefDict $resname]
        set pKaiMatrix [resDef2Matrix $resDef $value]
        dict set resDict $segresidname pKai $value
        dict set resDict $segresidname pKaiPair $pKaiMatrix 
        return $pKaiMatrix
    }
}

# ::cphSystem::resDef2Matrix
#
# Read a residue definiition for a given thermodynamic attribute (a free energy
# or pKa) and fill in the (flat)matrix elements that are derived from this.
# 
# This permits the configuration files to have "missing" information when it is
# not used. The downside is that essentially every case needs to have its own
# conventions.
#
# Arguments:
# ----------
# matrix : list
#   
#
# Returns
# -------
#  0 - The matrix was successfully updated.
# -1 - The specific combination of state definitions is not implemented
#      Ex. A state is unexpectedly missing and leaves a gap in the cycle or
#      too many definitions are given
# -2 - There is an inconsistency in the truncated definition
#      Ex. The sum of parameters around the cycle is not 0 - THIS MIGHT SIMPLY
#      IMPLY A NEGATIVE SIGN ERROR!
#
proc ::cphSystem::resDef2Matrix {resDef data} {
    #   The logic and mathematics here is a bit dreadful. Basically, we have
    # networks of states based on the number of protonation sites (numSites).
    # A given network has, at most, maxStates = 2**numSites states, but
    # generally many fewer than this, because some of the states don't exist
    # for physical reasons (e.g. doubly protonated carboxylates).
    #   The problems begin because the residue attributes (dG and pKa) do not
    # describe states, but rather transitions between pairs of states - there
    # are formally numStates**2 pairs! Fortunately, attributes between
    # identical states will always be zero and transitions between pairs of
    # states in opposite directions just have flipped signs. There are
    # therefore only numStates*(numStates - 1)/2 pairs of states that need to
    # be stored (as a flat, off-diagonal matrix, see index2flatindex, note that
    # this conversion is symmetric!). Lookup is always done with respect to the
    # full set of states (i.e. the size is computed with maxStates, not
    # numStates). The row/column indices are computed as a left-right binary
    # conversion of the occupancy vector to a scalar (see occupancy2Index).
    #  Just to make things more complicated, many systems have equivalent
    # sites such that many fewer pairwise values are needed. For example, a
    # carboxylate formally has two sites and four states, one of which is
    # neglected. The two sites are also equivalent, so this system has, at max
    # 2**2*(2**2 - 1)/2 = 6 pairs, but only three are used in practice. Two of
    # these pairs are the same, so only a single value is specified.
    #
    # Long story short: Each case is handled specially and implemented by hand
    # (YUCK!). Fortunately, it is the state pattern that is implemented, not
    # the residue name or some overly specific nonsense. If a pattern has been
    # implemented, then any residue following that pattern works also (e.g.
    # ACE, ASP, and GLU are all one pattern).
    #
    # Final note:
    #   Because only half of the diagonal elements are stored, the signs are
    # computed on-the-fly. This assumes that all proper pKa values are positive
    # - a negative pKa should never be inserted into pKa matrix. For "improper
    # pKas" (i.e. tautomerizations), the assumption is that the value is from
    # the higher index to the lower index (this may be negative).
    #   In the reduced listing all unique transitions from the highest index
    # should be specified first before going to the next highest index for
    # which unique transitions remain. In practice, that explanation is
    # probably pretty useless, so just look at the specific code for the given
    # case when making a new definition.
    #
    set states [dict get $resDef states]
    set numSites [llength [dict get $resDef states [lindex $states 0]]]
    set numStates [dict size $states]
    set maxStates [expr {int(pow(2, $numSites))}]
    set numPairs [expr {int($maxStates*($maxStates - 1) / 2)}]
    # Are any protonation counts missing?
    set nprotonsExists [lrepeat [expr {$numSites+1}] 0]
    dict for {state occ} $states {
        set nprotons [expr {lsum($occ)}]
        for {set i 0} {$i <= $numSites} {incr i} {
            if {$i == $nprotons} {
                lset nprotonsExists $i 1
            }
        }
    }

    set Matrix [lrepeat $numPairs 0.0]
    set numValues [llength $data]
    set errorCode 0
    switch -- $numSites 1 {
        switch -- $numValues 1 {
            lset Matrix [index2flatindex 1 0] $data
        } default {
            set errorCode -1 
        }
    } 2 {
        switch -- $numStates 3 {
            if {![lindex $nprotonsExists 0]} {;# state (0,0) is deleted
                switch -- $numValues 1 {
                    lassign $data attr32
                    set attr31 $attr32
                } 2 {;# Ex. HIS
                    lassign $data attr32 attr31
                } default {
                    set errorCode -1
                }
                lset Matrix [index2flatindex 3 2] $attr32
                lset Matrix [index2flatindex 3 1] $attr31
                lset Matrix [index2flatindex 2 1] [expr {$attr31 - $attr32}]
            } elseif {![lindex $nprotonsExists 2]} {;# state (1,1) is deleted
                switch -- $numValues 1 {;# Ex. carboxylates (ASP, GLU)
                    lassign $data attr20
                    set attr10 $attr20
                } 2 {
                    lassign $data attr20 attr10
                } default {
                    set errorCode -1
                }
                lset Matrix [index2flatindex 2 1] [expr {$attr20 - $attr10}]
                lset Matrix [index2flatindex 2 0] $attr20
                lset Matrix [index2flatindex 1 0] $attr10
            } else {
                # Why would state (0,1) or (1,0) be missing?
                set errorCode -1
            }
        } 4 {
            switch -- $numValues 4 {
                lassign $data attr32 attr31 attr20 attr10
            } default {
                set errorCode -1
            }
            set attr30 [expr {$attr31 + $attr10}]
            if {$attr30 != [expr {$attr32 + $attr20}]} {
                set errorCode -2
            }
            lset Matrix [index2flatindex 3 2] $attr32
            lset Matrix [index2flatindex 3 1] $attr31
            lset Matrix [index2flatindex 3 0] $attr30
            lset Matrix [index2flatindex 2 1] [expr {$attr20 - $attr10}]
            lset Matrix [index2flatindex 2 0] $attr20
            lset Matrix [index2flatindex 1 0] $attr10
        } default {
            set errorCode -1
        }
    } 3 {
        switch -- $numStates 4 {
            if {![lindex $nprotonsExists 0] && ![lindex $nprotonsExists 1]} {
                # state (0,0,0) and 1 proton states missing
                switch -- $numValues 1 {;# Ex. LYS
                    lassign $data attr76
                    set attr75 $attr76
                    set attr73 $attr76
                } default {
                    set errorCode -1
                }
                set attr53 [expr {$attr73 - $attr75}]
                set attr63 [expr {$attr73 - $attr76}]
                set attr65 [expr {$attr75 - $attr76}]
                if {$attr53 != [expr {$attr63 - $attr65}]
                    || $attr63 != [expr {$attr53 - $attr65}]
                    || $attr65 != [expr {$attr63 - $attr53}]} {
                    set errorCode -2
                }
                lset Matrix [index2flatindex 7 6] $attr76
                lset Matrix [index2flatindex 7 5] $attr75
                lset Matrix [index2flatindex 7 3] $attr73
                lset Matrix [index2flatindex 6 5] $attr65
                lset Matrix [index2flatindex 6 3] $attr63
                lset Matrix [index2flatindex 5 3] $attr53
            } elseif {![lindex $nprotonsExists 2]
                      && ![lindex $nprotonsExists 3]} {
                # state (1, 1, 1) and 2 proton states missing
                switch -- $numValues 1 {;# Ex. phosphate monoesters
                    lassign $data attr10
                    set attr20 $attr10
                    set attr40 $attr10
                } default {
                    set errorCode -1
                }
                lset Matrix [index2flatindex 1 0] $attr10
                lset Matrix [index2flatindex 2 0] $attr20
                lset Matrix [index2flatindex 4 0] $attr40
            } else {
                set errorCode -1
            }
        } 7 {
            if {![lindex $nprotonsExists 3]} {
                switch -- $numValues 2 {;# Ex. PHOsphate w/o H3PO4
                    lassign $data attr64 attr40
                    set attr62 $attr64
                    set attr61 $attr64
                    set attr54 $attr64
                    set attr52 $attr64
                    set attr51 $attr64
                    set attr34 $attr64
                    set attr32 $attr64
                    set attr31 $attr64 
                    set attr20 $attr40
                    set attr10 $attr40

                    set attr60 [expr {$attr64 + $attr40}]
                    set attr50 [expr {$attr54 + $attr40}]
                    set attr30 [expr {$attr32 + $attr20}]
                    # All other pairs are zero

                    lset Matrix [index2flatindex 6 4] $attr64
                    lset Matrix [index2flatindex 6 2] $attr62
                    lset Matrix [index2flatindex 6 1] $attr61
                    lset Matrix [index2flatindex 5 4] $attr54
                    lset Matrix [index2flatindex 5 2] $attr52
                    lset Matrix [index2flatindex 5 1] $attr51
                    lset Matrix [index2flatindex 3 4] $attr34
                    lset Matrix [index2flatindex 3 2] $attr32
                    lset Matrix [index2flatindex 3 1] $attr31
                    lset Matrix [index2flatindex 4 0] $attr40
                    lset Matrix [index2flatindex 2 0] $attr20
                    lset Matrix [index2flatindex 1 0] $attr10
                    lset Matrix [index2flatindex 6 0] $attr60
                    lset Matrix [index2flatindex 5 0] $attr50
                    lset Matrix [index2flatindex 3 0] $attr30
                } default {
                    set errorCode -1
                }
            } else {
                set errorCode -1
            }
        # TODO: Add H3PO4 with 8 states total 
        } default {
            set errorCode -1
        }
    } default {
        set errorCode -1
    }
    switch -- $errorCode -1 {
        abort "Bad or unimplemented specification of $numValues values for\
                $numSites site residue with $numStates states"
    } -2 {
        abort "Error in thermodynamic cycle"
    }
    return $Matrix
}

