/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMPUTENONBONDEDPAIR_H
#define COMPUTENONBONDEDPAIR_H

#include "ComputePatchPair.h"
#include "ComputeNonbondedUtil.h"

class ComputeNonbondedPair : public ComputePatchPair, private ComputeNonbondedUtil {

public:
  ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[]);
  ~ComputeNonbondedPair();

protected :
  virtual void initialize();
  virtual int noWork();
  virtual void doForce(Position* p[2], Results* r[2], AtomProperties* a[2]);

  PositionBox<Patch> *avgPositionBox[2];

  SubmitReduction *reduction;

};

#endif

