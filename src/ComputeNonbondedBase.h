/***************************************************************************/
/*                                                                         */
/*              (C) Copyright 1996 The Board of Trustees of the            */
/*                          University of Illinois                         */
/*                           All Rights Reserved                           */
/*									   */
/***************************************************************************/

/***************************************************************************
 * DESCRIPTION: ComputeNonbondedBase.h
 *
 ***************************************************************************/

// Several special cases are defined:
//   NBPAIR, NBSELF, NBEXCL switch environment (mutually exclusive)
//   MODIFY14 modified 1-4 parameters?
//   SWITCHING switching function?

#include "ComputeNonbondedHack.h"

#if defined DECLARATION || defined DEFINITION
void NODECL( CLASS :: ) NAME (Position* p PLEN, Force* f PLEN, AtomProperties* a PLEN) DECL( ; )
#endif
#ifdef DEFINITION
{
  BigReal electEnergy = 0;
  BigReal vdwEnergy = 0;

  const Real cutoff = Node::Object()->simParameters->cutoff;
  const BigReal cutoff2 = cutoff*cutoff;
  const BigReal dielectric_1 = 1/Node::Object()->simParameters->dielectric;

  const LJTable* const ljTable = LJTable::Instance();
  const Molecule* const mol = Node::Object()->molecule;

M14
(
  const BigReal scale14 = Node::Object()->simParameters->scale14;
)

SW
(
  const Real switchOn = Node::Object()->simParameters->switchingDist;
  const BigReal switchOn2 = switchOn*switchOn;
  const BigReal c0 = 1/(cutoff2-switchOn2);
  const BigReal c1 = c0*c0*c0;
  const BigReal c3 = c1 * 4;
  const BigReal c5 = 1/cutoff2;
  const BigReal c6 = -4 * c5;
)

  for ( int i = I_LOWER; i < I_UPPER; ++i )
  {
    const Position & p_i = p[I_SUB];
    const AtomProperties & a_i = a[I_SUB];

    Force & f_i = f[I_SUB];

    const BigReal M14( kq_i_u ) NOM14( kq_i ) =
    			COLOUMB * a_i.charge * dielectric_1;

    M14( const BigReal kq_i_s = kq_i_u * scale14; )

    for ( int j = J_LOWER; j < J_UPPER; ++j )
    {
      const Position & p_j = p[J_SUB];

      Vector f_vdw = p_i - p_j;

      const BigReal r2 = f_vdw.length2()

      if ( r2 > cutoff2 ) continue;

      M14 ( BigReal kq_i = kq_i_u; )

      const AtomProperties & a_j = a[J_SUB];

      const LJTable::TableEntry * NOM14( const ) lj_pars = 
		ljTable->table_val(a_i.type, a_j.type);
		
      if ( r2 <= lj_pars->exclcut2 )
      {
	if ( mol->checkexcl(a_i.id,a_j.id) ) continue;
	
M14
(
	else if ( mol->check14excl(a_i.id,a_j.id) )
	{
	  lj_pars = ljTable->table_val(a_i.type, a_j.type, 1);
	  kq_i = kq_i_s;
	}
)

      }

      Force & f_j = f[J_SUB];

      const BigReal r = sqrt(r2);
      const BigReal r_1 = 1/r;

SW
(
      BigReal switchVal;
      BigReal shiftVal;
      BigReal dSwitchVal;
      BigReal dShiftVal;

      if (r > switchOn)
      {
	const BigReal c2 = cutoff2-r2;
	const BigReal c4 = c2*(cutoff2+2*r2-3*switchOn2);
	switchVal = c2*c4*c1;
	dSwitchVal = c3*r*(c2*c2-c4);
      }
      else
      {
	switchVal = 1;
	dSwitchVal = 0;
      }

      shiftVal = 1 - r2*c5;
      dShiftVal = c6*shiftVal*r;
      shiftVal *= shiftVal;
)

      const BigReal kqq = kq_i * a_j.charge;

      BigReal f = kqq*r_1;
      
      electEnergy += f SW( * shiftVal );
      
      f *= r_1*( SW( shiftVal * ) r_1 SW( - dShiftVal ) );

      const Vector f_elec = f_vdw*f;

      f_i += f_elec;
      f_j -= f_elec;

      BigReal r_6 = r_1*r_1*r_1; r_6 *= r_6;
      const BigReal r_12 = r_6*r_6;

      const BigReal A = lj_pars->A;
      const BigReal B = lj_pars->B;

      const BigReal AmBterm = (A*r_6 - B) * r_6;

      vdwEnergy += switchVal*AmBterm;
      
      f_vdw *= ( SW( switchVal * ) 6 * (A*r_12 + AmBterm) *
			r_1 SW( - AmBterm*dSwitchVal ) )*r_1;

      f_i += f_vdw;
      f_j -= f_vdw;
    }
  }
}
#endif

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: ComputeNonbondedBase.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.4 $	$Date: 1996/11/08 02:37:12 $
 *
 ***************************************************************************
 * REVISION HISTORY:
 *
 * $Log: ComputeNonbondedBase.h,v $
 * Revision 1.4  1996/11/08 02:37:12  jim
 * split into two files to hide some preprocessor code from user
 *
 *
 ***************************************************************************/

