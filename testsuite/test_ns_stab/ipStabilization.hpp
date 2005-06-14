/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2004-09-22

 Copyright (C) 2004 EPFL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file ipStabilization.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-10-08
*/

#ifndef _IPSTABILIZATION_HPP
#define _IPSTABILIZATION_HPP

#define USE_OLD_PARAMETERS 0
#define WITH_DIVERGENCE 1

#include <life/lifecore/GetPot.hpp>


namespace LifeV
{

template<typename MESH, typename DOF>
class IPStabilization
{
public:

  //! Constructor
  IPStabilization( const GetPot& dataFile, 
		   const MESH&     mesh,
		   const DOF&      dof,
		   const RefFE&    refFE,
		   CurrentBdFE&    feBd,
		   const QuadRule& quadRule,
		   Real            viscosity );

    /*! compute IP stabilization terms and add them into matrix
     *  @param matrix matrix where the stabilization terms are added into
     *  @param state state vector for linearization of nonlinear stabilization
     */
    template<typename MATRIX, typename VECTOR>
    void apply( MATRIX& matrix, const VECTOR& state );

private:

    typedef ID ( *FTOP )( ID const localFace, ID const point );

    const MESH&  M_mesh;
    const DOF&   M_dof;
    CurrentFE    M_fe1;
    CurrentFE    M_fe2;
    CurrentBdFE& M_feBd;
    Real         M_gammaBeta;
    Real         M_gammaDiv;
    Real         M_gammaPress;
    Real         M_viscosity;
    ElemMat      M_elMatU;
    ElemMat      M_elMatP;
    FTOP         M_fToP;
}; // class IPStabilization

template<typename MESH, typename DOF>
 IPStabilization<MESH, DOF>::IPStabilization(  const GetPot& dataFile, 
					       const MESH&     mesh,
					       const DOF&      dof,
					       const RefFE&    refFE,
					       CurrentBdFE&    feBd,
					       const QuadRule& quadRule,                                          
					       Real            viscosity ) :
   M_mesh( mesh ),
   M_dof( dof ),
   M_fe1( refFE, getGeoMap(mesh), quadRule ),
   M_fe2( refFE, getGeoMap(mesh), quadRule ),
   M_feBd( feBd ),
   M_gammaBeta(  dataFile( "fluid/ipstab/gammaBeta", 0. ) ),
   M_gammaDiv(   dataFile( "fluid/ipstab/gammaDiv", 0. ) ),
   M_gammaPress( dataFile( "fluid/ipstab/gammaPress", 0. ) ),
   M_viscosity( viscosity ),
   M_elMatU( M_fe1.nbNode, nDimensions, nDimensions ),
   M_elMatP( M_fe1.nbNode, nDimensions+1, nDimensions+1 )
{
    switch( M_fe1.nbNode )
    {
        case 4:
            M_fToP = LinearTetra::fToP;
            break;
        case 10:
            M_fToP = QuadraticTetra::fToP;
            break;
        case 8:
            M_fToP = LinearHexa::fToP;
            break;
        case 20:
            M_fToP = QuadraticHexa::fToP;
            break;
        default:
            ERROR_MSG( "This refFE is not allowed with IP stabilisation" );
            break;
    }
}



template<typename MESH, typename DOF>
template<typename MATRIX, typename VECTOR>
void IPStabilization<MESH, DOF>::apply( MATRIX& matrix, const VECTOR& state )
{
    if ( M_gammaBeta == 0 && M_gammaDiv == 0 && M_gammaPress == 0 )
    {
        return;
    }

    Chrono chronoBeta;
    Chrono chronoUpdate;
    Chrono chronoElemComp;
    Chrono chronoAssembly;
    const UInt nDof = M_dof.numTotalDof();

    // local trace of the velocity
    ElemVec beta( M_feBd.nbNode, nDimensions );

    // loop on interior faces
    for ( UInt iFace = M_mesh.numBFaces()+1; iFace<= M_mesh.numFaces();
          ++iFace )
    {
        chronoUpdate.start();
        // update current finite elements
#if WITH_DIVERGENCE
        M_feBd.updateMeas( M_mesh.face( iFace ) );
#else
        M_feBd.updateMeasNormal( M_mesh.face( iFace ) );
        KNM<Real>& normal = M_feBd.normal;
#endif
        const Real hK2 = M_feBd.measure();
        const UInt iElAd1 = M_mesh.face( iFace ).ad_first();
        const UInt iElAd2 = M_mesh.face( iFace ).ad_second();
        M_fe1.updateFirstDeriv( M_mesh.volumeList( iElAd1 ) );
        M_fe2.updateFirstDeriv( M_mesh.volumeList( iElAd2 ) );
        chronoUpdate.stop();

        chronoBeta.start();
        // determine bmax = ||\beta||_{0,\infty,K}
        // first, get the local trace of the velocity into beta

        // local id of the face in its adjacent element
        UInt iFaEl = M_mesh.face( iFace ).pos_first();
        for ( int iNode = 0; iNode < M_feBd.nbNode; ++iNode )
        {
            UInt iloc = M_fToP( iFaEl, iNode+1 );
            for ( int iCoor = 0; iCoor < M_fe1.nbCoor; ++iCoor )
            {
                UInt ig = M_dof.localToGlobal( iElAd1, iloc+1 )-1+iCoor*nDof;
                beta.vec()[ iCoor*M_feBd.nbNode + iNode ] = state( ig );
            }
        }

        // second, calculate its max norm
        Real bmax = fabs( beta.vec()[ 0 ] );
        for ( int l = 1; l < int( M_fe1.nbCoor*M_feBd.nbNode ); ++l )
        {
            if ( bmax < fabs( beta.vec()[ l ] ) )
                bmax = fabs( beta.vec()[ l ] );
        }
        chronoBeta.stop();

        // pressure stabilization
        if ( M_gammaPress != 0.0 )
        {
#if USE_OLD_PARAMETERS
            Real coeffPress = M_gammaPress * hK2; // P1, P2 (code)
            //Real coeffPress = M_gammaPress * sqrt( hK2 ); // P1 p nonsmooth (code)
#else
            Real coeffPress = M_gammaPress * hK2 / // Pk (paper)
                std::max<Real>( bmax, M_viscosity/sqrt( hK2 ) );
#endif

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(coeffPress, M_elMatP, M_fe1, M_fe1, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly.start();
            assemb_mat(matrix, M_elMatP, M_fe1, M_dof,
                       nDimensions, nDimensions);
            chronoAssembly.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(coeffPress, M_elMatP, M_fe2, M_fe2, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly.start();
            assemb_mat(matrix, M_elMatP, M_fe2, M_dof,
                       nDimensions, nDimensions);
            chronoAssembly.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(-coeffPress, M_elMatP, M_fe1, M_fe2, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly.start();
            assemb_mat(matrix, M_elMatP, M_fe1, M_fe2, M_dof,
                       nDimensions, nDimensions);
            chronoAssembly.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(-coeffPress, M_elMatP, M_fe2, M_fe1, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly.start();
            assemb_mat(matrix, M_elMatP, M_fe2, M_fe1, M_dof,
                       nDimensions, nDimensions);
            chronoAssembly.stop();
        }

        // velocity stabilization
        if ( ( M_gammaDiv != 0 || M_gammaBeta != 0 ) && bmax > 0 )
        {
#if WITH_DIVERGENCE
#if USE_OLD_PARAMETERS
            Real coeffBeta = M_gammaBeta * hK2 / std::max<Real>(bmax, hK2); // code
#else
            Real coeffBeta = M_gammaBeta * hK2 / bmax; // paper
#endif

            Real coeffDiv = M_gammaDiv * hK2 * bmax; // (code and paper)
            //Real coeffDiv = M_gammaDiv * sqrt( hK2 ) * bmax; // ? (code)
#else
            // determine bnmax = ||\beta \cdot n||_{0,\infty,K}
            // and       bcmax = ||\beta \cross n||_{0,\infty,K}

            chronoBeta.start();
            Real bnmax = 0;
            Real bcmax = 0;
            for ( int iNode=0; iNode<M_feBd.nbNode; ++iNode ) {
                Real bn = 0;
                for ( int iCoor=0; iCoor<M_fe1.nbCoor; ++iCoor ) {
                    bn += normal(iNode, iCoor) *
                        beta.vec()[ iCoor*M_feBd.nbNode + iNode ];
                    bcmax = std::max<Real>
                        (bcmax, normal(iNode, (iCoor+1)%3) *
                         beta.vec()[ (iCoor+2)%3*M_feBd.nbNode + iNode ] -
                         normal(iNode, (iCoor+2)%3) *
                         beta.vec()[ (iCoor+1)%3*M_feBd.nbNode + iNode ]);
                }
                bnmax = std::max<Real> (bnmax, bn);
            }
            chronoBeta.stop();

            Real coeffGrad = hK2 * (M_gammaBeta*bnmax + M_gammaDiv*bcmax);
#endif
            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            ipstab_bgrad( coeffBeta, M_elMatU, M_fe1, M_fe1, beta,
                          M_feBd, 0, 0, nDimensions );
            ipstab_div( coeffDiv, M_elMatU, M_fe1, M_fe1, M_feBd );
#else
            ipstab_grad( coeffGrad, M_elMatU, M_fe1, M_fe1, M_feBd, 0, 0,
            nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    assemb_mat( matrix, M_elMatU, M_fe1, M_dof,
                                iComp, jComp );
            chronoAssembly.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            ipstab_bgrad( coeffBeta, M_elMatU, M_fe2, M_fe2, beta,
                          M_feBd, 0, 0, nDimensions );
            ipstab_div( coeffDiv, M_elMatU, M_fe2, M_fe2, M_feBd );
#else
            ipstab_grad( coeffGrad, M_elMatU, M_fe2, M_fe2, M_feBd, 0, 0,
            nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    assemb_mat( matrix, M_elMatU, M_fe2, M_dof,
                                iComp, jComp );
            chronoAssembly.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            ipstab_bgrad( -coeffBeta, M_elMatU, M_fe1, M_fe2, beta,
                          M_feBd, 0, 0, nDimensions );
            ipstab_div( -coeffDiv, M_elMatU, M_fe1, M_fe2, M_feBd );
#else
            ipstab_grad( -coeffGrad, M_elMatU, M_fe1, M_fe2, M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    assemb_mat( matrix, M_elMatU, M_fe1, M_fe2, M_dof,
                                iComp, jComp );
            chronoAssembly.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            ipstab_bgrad( -coeffBeta, M_elMatU, M_fe2, M_fe1, beta,
                         M_feBd, 0, 0, nDimensions );
            ipstab_div( -coeffDiv, M_elMatU, M_fe2, M_fe1, M_feBd );
#else
            ipstab_grad( -coeffGrad, M_elMatU, M_fe2, M_fe1, M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    assemb_mat( matrix, M_elMatU, M_fe2, M_fe1, M_dof,
                                iComp, jComp );
            chronoAssembly.stop();
        }

    } // loop on interior faces
    std::cout << std::endl;
    std::cout << "Updating of element   done in "
              << chronoUpdate.diff_cumul()   << "s." << std::endl;
    std::cout << "Determination of beta done in "
              << chronoBeta.diff_cumul()     << "s." << std::endl;
    std::cout << "Element computations  done in "
              << chronoElemComp.diff_cumul() << "s." << std::endl;
    std::cout << "Assembly              done in "
              << chronoAssembly.diff_cumul() << "s." << std::endl;

} // apply(...)


} // namespace LifeV

#endif /* _IPSTABILIZATION_HPP */
