/***************************************************************************/
/*      (C) Copyright 1996,1997 The Board of Trustees of the               */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/***************************************************************************/
/***************************************************************************
 * DESCRIPTION:  Forwards atoms to master node for force evaluation.
 *
 ***************************************************************************/

#ifndef COMPUTEGLOBALMASTER_H
#define COMPUTEGLOBALMASTER_H

#include "ComputeHomePatches.h"
#include "NamdTypes.h"

class ComputeGlobalConfigMsg;
class ComputeGlobalDataMsg;
class ComputeGlobalResultsMsg;
class ComputeGlobalMaster;
class ComputeMgr;

class ComputeGlobalMaster {
protected:
  friend class ComputeGlobal;
  ComputeGlobal *host;
  ComputeGlobalMaster(ComputeGlobal *);
  ~ComputeGlobalMaster();
  void recvData(ComputeGlobalDataMsg *);
  int msgcount;
  virtual void initialize();
  int initialized;
  void storedata(ComputeGlobalDataMsg *);
  void cleardata();
  AtomIDList aid;
  PositionList p;
  PositionList gcom;  // group centers of mass
  void storedefs(AtomIDList newgdef);
  AtomIDList gdef;  // group definitions
  ResizeArray<BigReal> gmass;  // group masses
  virtual void calculate();
};

#endif
/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeGlobalMaster.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 1998/02/16 00:24:37 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeGlobalMaster.h,v $
 * Revision 1.2  1998/02/16 00:24:37  jim
 * Added atom group centers of mass to Tcl interface.
 *
 * Revision 1.1  1998/02/10 05:35:04  jim
 * Split ComputeGlobal into different classes and files.
 * Switched globalForces and globalForcesTcl to tclForces and tclForcesScript.
 * Added (soon to be used) freeEnergy and freeEnergyConfig.
 *
 *
 *
 ***************************************************************************/

