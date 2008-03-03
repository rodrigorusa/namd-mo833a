/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef POSITIONBOX_H
#define POSITIONBOX_H

#include "NamdTypes.h"

template <class Owner> class PositionOwnerBox;

template <class Owner> class PositionBox {

  friend class PositionOwnerBox<Owner>;

  private:

    PositionBox(PositionOwnerBox<Owner>* o, int t=13) : ownerBox(o), trans(t) 
      { state = CLOSED; }

    ~PositionBox() {}
  
    enum box_state {OPEN, CLOSED} state;
    PositionOwnerBox<Owner> *ownerBox;
    int trans;

  public:

    CompAtom* open(void) { 
      if (state != OPEN) {
        state = OPEN; 
        ownerBox->openCount--;
      }
      return ownerBox->transData[trans];
    }

    CompAtom* open(int *num) { 
      *num = ownerBox->numData;
      if (state != OPEN) {
        state = OPEN; 
        ownerBox->openCount--;
      }
      return ownerBox->transData[trans];
    }

    // Closed access to the pointer
    void close(CompAtom ** const t) {
      if (state != CLOSED) {
        state = CLOSED;
        *t = 0;
    
        // Trigger callback!
        if ( ! --ownerBox->closeCount ) {
          ownerBox->close();
         }
      }
    }

    Owner *getOwner(){
        return ownerBox->getOwner();
    }

};

#endif // BOX_H
