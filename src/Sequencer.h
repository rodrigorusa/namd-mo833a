//-*-c++-*-
/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION:
 *
 ***************************************************************************/

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"

class HomePatch;
class SimParameters;
class ReductionMgr;
class CollectionMgr;

class Sequencer
{
    friend HomePatch;
public:
    Sequencer(HomePatch *p);
    ~Sequencer(void);
    void run(int numberOfCycles);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void suspend(void) { CthSuspend(); };
    void terminate(void) { CthFree(thread); CthSuspend(); };
    SimParameters *const simParams;	// for convenience
    int numberOfCycles;			// stores argument to run()
    int stepsPerCycle;			// stores info from run()
    HomePatch *const patch;		// access methods in patch
    ReductionMgr *const reduction;
    CollectionMgr *const collection;

private:
    CthThread thread;
    static void threadRun(Sequencer*);
};

#endif // SEQUENCER_H

