/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"
#include "PatchTypes.h"

class HomePatch;
class SimParameters;
class SubmitReduction;
class CollectionMgr;
class ControllerBroadcasts;
class LdbCoordinator;
class Random;

class Sequencer
{
    friend class HomePatch;
public:
    Sequencer(HomePatch *p);
    virtual ~Sequencer(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); }
    void suspend(void) { CthSuspend(); }

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void runComputeObjects(int migration = 0);

    void submitReductions(int);
    void submitCollections(int);

    void addForceToMomentum(BigReal, const int ftag = Results::normal);
    void addVelocityToPosition(BigReal);

    void rattle1(BigReal);
    void rattle2(BigReal,int);

    void maximumMove(BigReal);
    void minimizationQuenchVelocity(void);

    void rescaleVelocities(int);
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(BigReal,int);
    void berendsenPressure(int);
    void langevinPiston(int);
      int slowFreq;
    void langevinVelocities(BigReal);
    void langevinVelocitiesBBK1(BigReal);
    void langevinVelocitiesBBK2(BigReal);

    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    HomePatch *const patch;		// access methods in patch
    SubmitReduction *reduction;
    CollectionMgr *const collection;
    ControllerBroadcasts * broadcast;

    void rebalanceLoad(int timestep);

private:
    CthThread thread;
    static void threadRun(Sequencer*);

    LdbCoordinator *ldbCoordinator;
};

#endif // SEQUENCER_H

