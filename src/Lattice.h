//-*-c++-*-
#ifndef LATTICE_H
#define LATTICE_H

#include "NamdTypes.h"
#include <math.h>
#include "Tensor.h"

#define rint(X) floor((X)+0.5)

typedef Vector ScaledPosition;

class Lattice
{
public:
  Lattice(void) : p1(0), p2(0), p3(0) {};

  // maps a transformation triplet onto a single integer
  static int index(int i=0, int j=0, int k=0)
  {
    return 9 * (k+1) + 3 * (j+1) + (i+1);
  }

  // sets lattice basis vectors and origin (fixed center)
  void set(Vector A, Vector B, Vector C, Position Origin)
  {
    a1 = A; a2 = B; a3 = C; o = Origin;
    p1 = ( a1.length2() ? 1 : 0 );
    p2 = ( a2.length2() ? 1 : 0 );
    p3 = ( a3.length2() ? 1 : 0 );
    if ( ! p1 ) a1 = Vector(1.0,0.0,0.0);
    if ( ! p2 ) {
      Vector u1 = a1 / a1.length();
      Vector e_z(0.0,0.0,1.0);
      if ( fabs(e_z * u1) < 0.9 ) { a2 = cross(e_z,a1); }
      else { a2 = cross(Vector(1.0,0.0,0.0),a1); }
      a2 /= a2.length();
    }
    if ( ! p3 ) {
      a3 = cross(a1,a2);
      a3 /= a3.length();
    }
    if ( volume() < 0.0 ) a3 *= -1.0;
    recalculate();
  }

  // rescale lattice dimensions by factor, origin doesn't move
  void rescale(Tensor factor)
  {
    a1 = factor * a1;
    a2 = factor * a2;
    a3 = factor * a3;
    recalculate();
  }

  // rescale a position, keeping origin constant, assume 3D
  void rescale(Position &p, Tensor factor) const
  {
    p -= o;
    p = factor * p;
    p += o;
  }

  // transform scaled position to unscaled position
  Position unscale(ScaledPosition s) const
  {
    return (o + a1*s.x + a2*s.y + a3*s.z);
  }

  // transform unscaled position to scaled position
  ScaledPosition scale(Position p) const
  {
    p -= o;
    return Vector(b1*p,b2*p,b3*p);
  }

  // transforms a position nearest to a SCALED reference position
  Position nearest(Position data, ScaledPosition ref) const
  {
    ScaledPosition sn = scale(data);
    if ( p1 ) {
      BigReal tmp = sn.x - ref.x;
      sn.x = ref.x + tmp - rint(tmp);
    }
    if ( p2 ) {
      BigReal tmp = sn.y - ref.y;
      sn.y = ref.y + tmp - rint(tmp);
    }
    if ( p3 ) {
      BigReal tmp = sn.z - ref.z;
      sn.z = ref.z + tmp - rint(tmp);
    }
    return unscale(sn);
  }

  // transforms a position nearest to a SCALED reference position
  // adds transform for later reversal
  Position nearest(Position data, ScaledPosition ref, Transform *t) const
  {
    ScaledPosition sn = scale(data);
    if ( p1 ) {
      BigReal tmp = sn.x - ref.x;
      BigReal rit = rint(tmp);
      sn.x = ref.x + tmp - rit;
      t->i -= (int) rit;
    }
    if ( p2 ) {
      BigReal tmp = sn.y - ref.y;
      BigReal rit = rint(tmp);
      sn.y = ref.y + tmp - rit;
      t->j -= (int) rit;
    }
    if ( p3 ) {
      BigReal tmp = sn.z - ref.z;
      BigReal rit = rint(tmp);
      sn.z = ref.z + tmp - rit;
      t->k -= (int) rit;
    }
    return unscale(sn);
  }

  // reverses cumulative transformations for output
  Position reverse_transform(Position data, const Transform &t) const
  {
    return ( data - t.i*a1 - t.j*a2 - t.k*a3 );
  }

  // calculates shortest vector from p2 to p1 (equivalent to p1 - p2)
  Vector delta(Position pos1, Position pos2) const
  {
    Vector diff = pos1 - pos2;
    Vector result = diff;
    if ( p1 ) result -= a1*rint(b1*diff);
    if ( p2 ) result -= a2*rint(b2*diff);
    if ( p3 ) result -= a3*rint(b3*diff);
    return result;
  }

  // calculates shortest vector from origin to p1 (equivalent to p1 - o)
  Vector delta(Position pos1) const
  {
    Vector diff = pos1 - o;
    Vector result = diff;
    if ( p1 ) result -= a1*rint(b1*diff);
    if ( p2 ) result -= a2*rint(b2*diff);
    if ( p3 ) result -= a3*rint(b3*diff);
    return result;
  }

  Position* create(Position *d, int n, int i) const
  {
    Position *dt;
    if ( i != 13 )
    {
      dt = new Position[n];
      Vector shift = (i%3-1) * a1 + ((i/3)%3-1) * a2 + (i/9-1) * a3;
      for( int j = 0; j < n; ++j )
        dt[j] = d[j] + shift;
    }
    else
    {
      dt = d;
    }
    return dt;
  }

  void destroy(Position **d, int i) const
  {
    if ( i != 13 ) delete [] *d;
    *d = NULL;
  }

  // lattice vectors
  Vector a() const { return a1; }
  Vector b() const { return a2; }
  Vector c() const { return a3; }

  // only if along x y z axes
  int orthogonal() const {
    return ( ! ( a1.y || a1.z || a2.x || a2.z || a3.x || a3.y ) );
  }

  // origin (fixed center of cell)
  Vector origin() const
  {
    return o;
  }

  // reciprocal lattice vectors
  Vector a_r() const { return b1; }
  Vector b_r() const { return b2; }
  Vector c_r() const { return b3; }

  // periodic along this direction
  int a_p() const { return p1; }
  int b_p() const { return p2; }
  int c_p() const { return p3; }

  BigReal volume(void) const
  {
    return ( p1 && p2 && p3 ? cross(a1,a2) * a3 : 0.0 );
  }

private:
  Vector a1,a2,a3; // real lattice vectors
  Vector b1,b2,b3; // reciprocal lattice vectors (more or less)
  Vector o; // origin (fixed center of cell)
  int p1, p2, p3; // periodic along this lattice vector?

  // calculate reciprocal lattice vectors
  void recalculate(void) {
    {
      Vector c = cross(a2,a3);
      b1 = c / ( a1 * c );
    }
    {
      Vector c = cross(a3,a1);
      b2 = c / ( a2 * c );
    }
    {
      Vector c = cross(a1,a2);
      b3 = c / ( a3 * c );
    }
  }

};

#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile $
 *	$Author $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1011 $	$Date: 1999/09/03 20:46:15 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: Lattice.h,v $
 * Revision 1.1011  1999/09/03 20:46:15  jim
 * Support for non-orthogonal periodic boundary conditions.
 *
 * Revision 1.1010  1999/01/06 19:19:20  jim
 * Broadcast and Sequencers understand anisotropic volume rescaling factors.
 *
 * Revision 1.1009  1998/08/21 01:15:04  jim
 * Eliminated warnings.
 *
 * Revision 1.1008  1998/08/11 16:30:28  jim
 * Modified output from periodic boundary simulations to return atoms to
 * internally consistent coordinates.  We store the transformations which
 * were performed and undo them at the end.  It might be better to do this
 * by always keeping the original coordinates and only doing the transform
 * for the nonbonded terms but this works for now.
 *
 * Revision 1.1007  1998/04/06 16:34:08  jim
 * Added DPME (single processor only), test mode, and momenta printing.
 *
 * Revision 1.1006  1998/03/30 21:01:17  jim
 * Added nearest-image support for periodic boundary conditions to full direct.
 *
 * Revision 1.1005  1998/03/26 23:28:29  jim
 * Small changes for KCC port.  Altered use of strstream in ComputeFreeEnergy.
 *
 * Revision 1.1004  1997/03/27 08:04:18  jim
 * Reworked Lattice to keep center of cell fixed during rescaling.
 *
 * Revision 1.1003  1997/03/21 23:05:36  jim
 * Added Berendsen's pressure coupling method, won't work with MTS yet.
 *
 * Revision 1.1002  1997/03/19 11:54:24  ari
 * Add Broadcast mechanism.
 * Fixed RCS Log entries on files that did not have Log entries.
 * Added some register variables to Molecule and ComputeNonbondedExcl.C
 *
 ***************************************************************************/
