/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBALEASY_H
#define COMPUTEGLOBALEASY_H

#include "ComputeGlobalMaster.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;
class Molecule;
class SimParameters;
class SubmitReduction;

class ComputeGlobalEasy : public ComputeGlobalMaster {
protected:
  friend class ComputeGlobal;
  ComputeGlobalEasy(ComputeGlobal *, const char *);
  virtual ~ComputeGlobalEasy();

  int getAtomID(const char *segid, int resid, const char *aname);
  int getNumAtoms(const char* segid, int resid); // 0 on error
  int getAtomID(const char *segid, int resid, int index);
  double getMass(int atomid);
  int requestAtom(int atomid);
  int getPosition(int atomid, Position &position);
  int addForce(int atomid, Force force);
  void addEnergy(BigReal);

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

private:

  virtual void initialize();
  virtual void calculate();

  ComputeGlobalConfigMsg *configMsg;
  ComputeGlobalResultsMsg *resultsMsg;
  Molecule *molecule;
  SimParameters *simParams;
  SubmitReduction *reduction;

  char *configName;
  BigReal energy;

};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobalEasy.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1999/07/06 20:32:41 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalEasy.h,v $
 * Revision 1.4  1999/07/06 20:32:41  jim
 * Eliminated warnings from new generation of picky compilers.
 *
 * Revision 1.3  1999/06/17 17:05:38  jim
 * Renamed seq to step in most places.  Now has meaning only to user.
 *
 * Revision 1.2  1999/06/17 15:46:07  jim
 * Completely rewrote reduction system to eliminate need for sequence numbers.
 *
 * Revision 1.1  1999/06/03 16:50:08  jim
 * Added simplified interface to ComputeGlobal mechanism.
 *
 *
 *
 ***************************************************************************/

