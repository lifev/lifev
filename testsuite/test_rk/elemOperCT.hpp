/**
 * Useful elemental operators for projection methods. 
 * __todo__ : merge with elemOper.hpp/elemOper.cpp. 
 */ 

#ifndef __ELEMOPER_CT_HH
#define __ELEMOPER_CT_HH

#include <life/lifecore/life.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/currentFE.hpp>

namespace LifeV
{

void source_pdivv(Real alpha, ElemVec& pLoc, ElemVec& elvec, 
                  const CurrentFE& fe_p, const CurrentFE& fe_u, const int iblock)
{
    int i, j, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;

    for (i=0; i < fe_u.nbNode; i++)
    {
        s = 0;
	for (iq = 0; iq < fe_u.nbQuadPt; ++iq)
	    for (j = 0; j < fe_p.nbNode; ++j)
	        s += pLoc[j]*fe_p.phi(j,iq)*fe_u.phiDer(i,iblock,iq)*fe_u.weightDet(iq);
        vec( i ) += s*alpha; 
    }
} // source_pdivv

/**
 * Interior penalty term for convective stabilization: explicit case. 
 * On input elemental velocity uLoc on fe_u = fe1. 
 * On output elemental iblock elvec on fe = fe2.
 */

void ipstab_grad_expl( const Real         coef,
                       ElemVec&           uLoc,	
                       ElemVec&           elvec,
                       const CurrentFE&   fe1,
                       const CurrentFE&   fe2,
                       const CurrentBdFE& bdfe,
                       int iblock )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemVec::vector_view vec = elvec.block(iblock);

    ElemMat::matrix_type mat_tmp( fe1.nbNode, fe2.nbNode );

    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }


    // Loop on rows
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }

    // Compute i_th block of elemental stabilization vector
    for (j=0; j<fe2.nbNode; ++j)
    {
    	sum = 0.0;
	for (i=0; i<fe1.nbNode; ++i)
	    sum += uLoc[iblock * fe1.nbNode + i ] * mat_tmp(i,j);
	vec(j) += sum;
    }

} // ipstab_grad_expl

/**
 * Interior penalty term for convective stabilization: explicit case. 
 * On input elemental velocity uLoc on fe_u = fe1. 
 * On output elemental elvec on fe = fe2.
 * __note__: better use this overload as it avoids too many calls.
 */

void ipstab_grad_expl( const Real         coef,
                       ElemVec&           uLoc,	
                       ElemVec&           elvec,
                       const CurrentFE&   fe1,
                       const CurrentFE&   fe2,
                       const CurrentBdFE& bdfe )
{
    /*
      Interior penalty stabilization: coef*\int_{face} grad u1_i . grad v1_j
    */

    ASSERT_PRE( fe_u.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );
    ASSERT_PRE( fe.hasFirstDeriv(),
                "ipstab11 needs at least the first derivatives" );

    ElemMat::matrix_type mat_tmp( fe1.nbNode, fe2.nbNode );

    Real sum, sum1, sum2;
    int i, j, ig, icoor, jcoor;
    Real x[ 3 ], rx1[ 3 ], drp1[ 3 ], rx2[ 3 ], drp2[ 3 ];
    Real phid1[ fe1.nbNode ][ fe1.nbCoor ][ bdfe.nbQuadPt ];
    Real phid2[ fe2.nbNode ][ fe2.nbCoor ][ bdfe.nbQuadPt ];
    Real b1[ 3 ], b2[ 3 ];

    fe1.coorMap( b1[ 0 ], b1[ 1 ], b1[ 2 ], 0, 0, 0 ); // translation fe1
    fe2.coorMap( b2[ 0 ], b2[ 1 ], b2[ 2 ], 0, 0, 0 ); // translation fe2

    for ( int ig = 0; ig < bdfe.nbQuadPt; ++ig )
    {  // first derivatives on quadrature points
        bdfe.coorQuadPt( x[ 0 ], x[ 1 ], x[ 2 ], ig );       // quadrature points coordinates

        // local coordinates of the quadrature point
        for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
        {
            sum1 = 0;
            sum2 = 0;
            for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
            {
                sum1 += fe1.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b1[ jcoor ] );
                sum2 += fe2.tInvJac( jcoor, icoor, 0 ) * ( x[ jcoor ] - b2[ jcoor ] );
            }
            rx1[ icoor ] = sum1;
            rx2[ icoor ] = sum2;
        }

        for ( i = 0; i < fe1.nbNode; ++i )
        {

            // first derivative on the reference element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                drp1[ icoor ] = fe1.refFE.dPhi( i, icoor, rx1[ 0 ], rx1[ 1 ], rx1[ 2 ] );
                drp2[ icoor ] = fe2.refFE.dPhi( i, icoor, rx2[ 0 ], rx2[ 1 ], rx2[ 2 ] );
            }

            // first derivative on the current element
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
            {
                sum1 = 0;
                sum2 = 0;
                for ( jcoor = 0; jcoor < fe1.nbCoor; ++jcoor )
                {
                    sum1 += fe1.tInvJac( icoor, jcoor, 0 ) * drp1[ jcoor ];
                    sum2 += fe2.tInvJac( icoor, jcoor, 0 ) * drp2[ jcoor ];
                }
                phid1[ i ][ icoor ][ ig ] = sum1;
                phid2[ i ][ icoor ][ ig ] = sum2;
            }
        }
    }


    // Loop on rows
    for ( i = 0; i < fe1.nbNode; ++i )
    {
        // Loop on columns
        for ( j = 0; j < fe2.nbNode; ++j )
        {
            sum = 0.0;
            // Loop on coordinates
            for ( icoor = 0; icoor < fe1.nbCoor; ++icoor )
                for ( ig = 0; ig < bdfe.nbQuadPt ; ++ig )
                    sum += phid1[ i ][ icoor ][ ig ] * phid2[ j ][ icoor ][ ig ] * bdfe.weightMeas( ig );
            mat_tmp( i, j ) = coef * sum;
        }
    }

    // Compute __all__ blocks of elemental stabilization vector
    for (j=0; j<fe2.nbNode; ++j)
    {
    	for (jcoor=0; jcoor<fe2.nbCoor; ++jcoor) {
	    sum = 0.0;
	    for (i=0; i<fe1.nbNode; ++i)
	        sum += uLoc[jcoor * fe1.nbNode + i ] * mat_tmp(i,j);
	    elvec.vec()[j+jcoor*fe2.nbNode] += sum;
	}
    }

} // ipstab_grad_expl

/**
 * Explicit local convective term - < v | u. \nabla u >, with a __minus__ sign. 
 * __note__: this version directly works for local vector fields (i.e. velocity) 
 * and only one currentFE as one needs no more for pure NS and it reduces calls.
 * On input : uLoc elemental velocity on fe
 * On output : elvec elemental explicit convective term on the same fe.
 *
 * __todo__: add coef argument in function for correct mass term generation.
 * __see__ : do not use this one (for debug/test purpose only) but the one below.
 */

void source_ugradu(const ElemVec& uLoc, ElemVec& elvec, const CurrentFE& fe)
{
    int iQuad, iCoor, jCoor, iNode, jNode;
    double s;
    ElemMat::matrix_type mat_tmp (nDimensions, fe.nbQuadPt);

    for (iCoor = 0; iCoor < nDimensions; ++iCoor)
    {
      for (iQuad = 0; iQuad < fe.nbQuadPt; ++iQuad)
      {
        s = 0.;
	// compute u.\nabla u at gauss point for each component
	for (iNode = 0; iNode < fe.nbNode; ++iNode)
	{
	  for (jNode = 0; jNode < fe.nbNode; ++jNode)
	  {
	    for (jCoor = 0; jCoor < nDimensions; ++jCoor)
	    {
	      s += uLoc.vec()[jNode + jCoor * fe.nbNode] * fe.phi(jNode, iQuad) * 
	            fe.phiDer(iNode, jCoor, iQuad) * uLoc.vec()[iNode+ iCoor * fe.nbNode];
	    }
	  }
	}
        
        mat_tmp(iCoor, iQuad) = s;

      }
    }

    // compute convective term
    for (iNode = 0; iNode < fe.nbNode; ++iNode)
    {
        for (iCoor = 0; iCoor < nDimensions; ++iCoor)
	{
	    s = 0.;
	    for (iQuad=0; iQuad < fe.nbQuadPt; ++iQuad)
	    {
	        // here we put the minus sign
		s -= mat_tmp(iCoor, iQuad) * fe.phi(iNode, iQuad) * fe.weightDet(iQuad);
	    }
	    elvec.vec()[iNode + iCoor * fe.nbNode] += s;
	}
    }

} // source_ugradu

/**
 * Explicit local convective term - coef * < v | u1. \nabla 2 >. 
 * __note__: this version directly works for local vector fields (i.e. velocity) 
 * and only one currentFE as one needs no more for pure NS and it reduces calls.
 * On input : uLoc1 and uLoc2  elemental velocities on fe
 * On output : elvec elemental explicit convective term on the same fe.
 *
 */

void source_u1gradu2(const Real coef, const ElemVec& uLoc1, 
                                      const ElemVec& uLoc2, 
				      ElemVec& elvec, 
				      const CurrentFE& fe)
{
    int iQuad, iCoor, jCoor, iNode, jNode;
    double s;
    ElemMat::matrix_type mat_tmp (nDimensions, fe.nbQuadPt);

    for (iCoor = 0; iCoor < nDimensions; ++iCoor)
    {
      for (iQuad = 0; iQuad < fe.nbQuadPt; ++iQuad)
      {
        s = 0.;
	// compute u1.\nabla u2 at gauss point for each component
	for (iNode = 0; iNode < fe.nbNode; ++iNode)
	{
	  for (jNode = 0; jNode < fe.nbNode; ++jNode)
	  {
	    for (jCoor = 0; jCoor < nDimensions; ++jCoor)
	    {
	      s += uLoc1.vec()[jNode + jCoor * fe.nbNode] * fe.phi(jNode, iQuad) * 
	            fe.phiDer(iNode, jCoor, iQuad) * uLoc2.vec()[iNode+ iCoor * fe.nbNode];
	    }
	  }
	}
        
        mat_tmp(iCoor, iQuad) = s;

      }
    }

    // compute convective term
    for (iNode = 0; iNode < fe.nbNode; ++iNode)
    {
        for (iCoor = 0; iCoor < nDimensions; ++iCoor)
	{
	    s = 0.;
	    for (iQuad=0; iQuad < fe.nbQuadPt; ++iQuad)
	    {
		s += mat_tmp(iCoor, iQuad) * fe.phi(iNode, iQuad) * fe.weightDet(iQuad);
	    }
	    elvec.vec()[iNode + iCoor * fe.nbNode] += coef * s;
	}
    }

} // source_u1gradu2

/**
 * Explicit stiffness strain for Crank-Nicholson step in RK2 like Chorin-Temam method
 * On input : uLoc elemental velocity vector on fe 
 * On output : all blocks elvec elemental stiff vector on fe
 * __note__: coef must be __minus one__ time the viscosity when used. 
 */

void source_stiff_strain( const Real coef, const ElemVec& uLoc, ElemVec& elvec, 
                          const CurrentFE& fe)
{
    ASSERT_PRE( fe.hasFirstDeriv(),
                "Source Stiffness Strain matrix needs at least the first derivatives" );
    double s;
    double tmp = coef * 0.5;

    ElemMat::matrix_type mat_tmp( fe.nbNode, fe.nbNode );
    ElemMat elmat(fe.nbNode, nDimensions, nDimensions);
    elmat.zero();

    for ( int i = 0; i < fe.nbNode; ++i )
    {
        for ( int j = 0; j < fe.nbNode; ++j )
        {
            s = 0;
            for ( int ig = 0; ig < fe.nbQuadPt; ++ig )
                for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
                    s += fe.phiDer( i, icoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
            mat_tmp( i, j ) = tmp * s;
        }
    }
    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        ElemMat::matrix_view mat = elmat.block( icoor, icoor );
        mat += mat_tmp;
    }

    for ( int icoor = 0; icoor < fe.nbCoor; ++icoor )
    {
        for ( int jcoor = 0; jcoor < fe.nbCoor; ++jcoor )
        {
            ElemMat::matrix_view mat = elmat.block( icoor, jcoor );
            for ( int i = 0; i < fe.nbNode; ++i )
            {
                for ( int j = 0; j < fe.nbNode; ++j )
                {
                    s = 0;
                    for ( int ig = 0; ig < fe.nbQuadPt; ++ig )
                        s += fe.phiDer( i, jcoor, ig ) * fe.phiDer( j, icoor, ig ) * fe.weightDet( ig );
                    mat( i, j ) += tmp * s;
                }
            }
        }
    }

    // get explicit term 
    for (int icoor = 0; icoor < fe.nbCoor; ++icoor)
    {
        for (int i = 0; i<fe.nbNode; ++i)
        {
  	    s = 0.0;
	    for (int jcoor = 0; jcoor < fe.nbCoor; ++jcoor)
	    {
	        for (int j = 0; j<fe.nbNode; ++j)
	        {
	              // __todo__: double check the call to the KNM from .mat() below 
		      s += uLoc[j+jcoor*fe.nbNode] * elmat.mat()(j+jcoor*fe.nbNode,i+icoor*fe.nbNode);
	        }
            }
	    elvec.vec()[i+icoor*fe.nbNode] += s;
        }
    }

} // source_stiff_strain




} // namespace LifeV

#endif // __ELEMOPER_CT_HH
