/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************/
/* DESCRIPTION:							           */	
/*									   */
/***************************************************************************/

#ifndef _NODE_H
#define _NODE_H

#include "ckdefs.h"
#include "chare.h"
#include "c++interface.h"

#include "main.h"
#include "PatchMap.h"
#include "ComputeMap.h"
#include "PatchMgr.h"

class Molecule;
class Parameters;
class SimParameters;
class ConfigList;
class PDB;
class WorkDistrib;
class PatchMgr;

class Node : public groupmember
{
public:
  Molecule *molecule;
  Parameters *parameters;
  SimParameters *simParameters;
  ConfigList *configList;
  PDB *pdb;

  int workDistribGroup;
  WorkDistrib *workDistrib;
  int patchMgrGroup;
  PatchMgr *patchMgr;

  PatchMap patchMap;
  ComputeMap computeMap;

  // Charm Entry point - distributed contructor
  Node(NodeInitMsg *msg);
  ~Node(void);

  int myid(void);		   
  int numNodes(void);		   
  
  // Charm Entry point - Read in system data, get all ready to simulate
  void startup(InitMsg *initmsg);  

  void saveMolDataPointers(Molecule *, Parameters *,
			   SimParameters *, ConfigList *,
			   PDB *);

  // Run for the number of steps specified in the sim_parameters
  void run(void);                  
};

#endif /* _NODE_H */

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.9 $	$Date: 1996/10/16 08:22:39 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Node.h,v $
 * Revision 1.9  1996/10/16 08:22:39  ari
 * *** empty log message ***
 *
 * Revision 1.8  1996/09/03 22:54:25  ari
 * *** empty log message ***
 *
 * Revision 1.7  1996/08/23 22:03:52  brunner
 * Made WorkdDistrib, PatchMgr public members
 *
 * Revision 1.6  1996/08/21 23:58:25  brunner
 * *** empty log message ***
 *
 * Revision 1.5  1996/08/19 22:05:31  ari
 * *** empty log message ***
 *
 * Revision 1.4  1996/08/19 17:57:47  ari
 * *** empty log message ***
 *
 * Revision 1.3  1996/08/16 21:56:17  brunner
 * *** empty log message ***
 *
 * Revision 1.2  1996/08/16 21:42:58  brunner
 * *** empty log message ***
 *
 * Revision 1.1  1996/08/16 21:19:34  ari
 * Initial revision
 *
 *
 ***************************************************************************/
