
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Master BOC.  coordinates startup, close down of each PE
   Also owns pointers to common objects needed by system		
   Many utility static methods are owned by Node.
*/

#ifndef _SYNC_H
#define _SYNC_H

#include "charm++.h"

#include "ProcessorPrivate.h"
#include "Sync.decl.h"
#include "ResizeArrayPrimIter.h"

extern int useSync;

typedef ResizeArrayPrimIter<ComputeID> ComputeIDListIter;

class Sync : public BOCclass
{
private:
    struct _clist {
    int pid;
    int step;
    ComputeIDListIter cid;
    int doneMigration;
    } *clist;
    int capacity;

    int step;
    int counter;
    int cnum;
    int nPatcheReady;
    int numPatches;

    void triggerCompute();
public:
    Sync(void);
    ~Sync(void);
    inline static Sync *Object() { return CpvAccess(Sync_instance); }
    void openSync(); 
    void registerComp(PatchID pid, ComputeIDListIter cid, int doneMigration);
    void PatchReady(void);
};

#endif
