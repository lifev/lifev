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
   \file nsipterms.hpp
   \author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
   \date 2004-10-08
*/

#ifndef _NSIPTERMS_HPP
#define _NSIPTERMS_HPP

#include <life/lifecore/chrono.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <boost/shared_ptr.hpp>
#include <elemOperCT.hpp>

#define USE_OLD_PARAMETERS 0
#define WITH_DIVERGENCE 1

namespace LifeV
{

//namespace details
//{

template<typename MESH, typename DOF>
class IPStabilization
{
public:

    typedef boost::shared_ptr<MESH> mesh_type;

    //! Constructor
    IPStabilization( const mesh_type     mesh,
                     const DOF&      dof,
                     const RefFE&    refFE,
                     CurrentBdFE&    feBd,
                     const QuadRule& quadRule,
                     Real            gammaBeta = 0,
                     Real            gammaDiv = 0,
                     Real            gammaPress = 0,
                     Real            viscosity = 1 );

    /*! compute IP stabilization terms and add them into matrix
     *  @param matrix matrix where the stabilization terms are added into
     *  @param state state vector for linearization of nonlinear stabilization
     */
    template<typename MATRIX, typename VECTOR>
    void apply( MATRIX& matrix, const VECTOR& state, bool verbose = true );

    /*! compute IP explicit stabilization and add them into rhs vector
     * @param vector vector where the stabilization terms are added to
     * @param state state vector for linearization of nonlinear stabilization
     */
    template<typename VECTOR>
    void apply_expl(VECTOR& vector, const VECTOR& beta);

    void setGammaBeta (double gammaBeta) { M_gammaBeta  = gammaBeta;}
    void setGammaDiv  (double gammaDiv)  { M_gammaDiv   = gammaDiv;}
    void setGammaPress(double gammaPress){ M_gammaPress = gammaPress;}
private:

    typedef ID ( *FTOP )( ID const localFace, ID const point );

    const mesh_type  M_mesh;

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
    ElemVec      M_elvec_u;
    FTOP         M_fToP;
}; // class IPStabilization

template<typename MESH, typename DOF>
IPStabilization<MESH, DOF>::IPStabilization( const mesh_type mesh,
                                             const DOF&      dof,
                                             const RefFE&    refFE,
                                             CurrentBdFE&    feBd,
                                             const QuadRule& quadRule,
                                             Real            gammaBeta,
                                             Real            gammaDiv,
                                             Real            gammaPress,
                                             Real            viscosity ) :
    M_mesh      ( mesh ),
    M_dof       ( dof ),
    M_fe1       ( refFE, getGeoMap(*mesh), quadRule ),
    M_fe2       ( refFE, getGeoMap(*mesh), quadRule ),
    M_feBd      ( feBd ),
    M_gammaBeta ( gammaBeta ),
    M_gammaDiv  ( gammaDiv ),
    M_gammaPress( gammaPress ),
    M_viscosity ( viscosity ),
    M_elMatU    ( M_fe1.nbNode, nDimensions    , nDimensions   ),
    M_elMatP    ( M_fe1.nbNode, nDimensions + 1, nDimensions+1 ),
    M_elvec_u   ( M_fe1.nbNode, nDimensions)
{
    switch( M_fe1.nbNode )
    {
        case 4:
            M_fToP = LinearTetra::fToP;
            break;
        case 5:
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
//            ERROR_MSG( "This refFE is not allowed with IP stabilisation" );
            break;
    }
}



template<typename MESH, typename DOF>
template<typename MATRIX, typename VECTOR>
void IPStabilization<MESH, DOF>::apply( MATRIX& matrix,  const VECTOR& state, const bool verbose )
{
    if ( M_gammaBeta == 0. && M_gammaDiv == 0. && M_gammaPress == 0. )
    {
        return;
    }

    ChronoFake chronoUpdate;
    ChronoFake chronoBeta;
    ChronoFake chronoElemComp;
    ChronoFake chronoAssembly1;
    ChronoFake chronoAssembly2;
    ChronoFake chronoAssembly3;
    ChronoFake chronoAssembly4;
    ChronoFake chronoAssembly5;
    ChronoFake chronoAssembly6;
    ChronoFake chronoAssembly7;
    ChronoFake chronoAssembly8;
    ChronoFake chronoAssembly9;
    Chrono chronoAssembly;

    const UInt nDof = M_dof.numTotalDof();

    double normInf;
    state.NormInf(&normInf);

    // local trace of the velocity
    ElemVec beta( M_feBd.nbNode, nDimensions );

    UInt myFaces(0);

    chronoAssembly.start();
    // loop on interior faces
    for ( UInt iFace = M_mesh->numBFaces() + 1; iFace<= M_mesh->numFaces();
          ++iFace )
    {
        const UInt iElAd1 = M_mesh->face( iFace ).ad_first();
        const UInt iElAd2 = M_mesh->face( iFace ).ad_second();

        if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
        {
            //std::cout << "iElAd1 = " << iElAd1 << "; iElAd2 = " << iElAd2 << std::endl;
            continue;
        }
        ++myFaces;

        chronoUpdate.start();
        // update current finite elements
#if WITH_DIVERGENCE
        M_feBd.updateMeas( M_mesh->face( iFace ) );
#else
        M_feBd.updateMeasNormal( M_mesh->face( iFace ) );
        KNM<Real>& normal = M_feBd.normal;
#endif
        const Real hK2 = M_feBd.measure();


        M_fe1.updateFirstDeriv( M_mesh->volumeList( iElAd1 ) );
        M_fe2.updateFirstDeriv( M_mesh->volumeList( iElAd2 ) );
        chronoUpdate.stop();

        Real bmax(0);
        if (normInf != 0.)
        {
            chronoBeta.start();
            // determine bmax = ||\beta||_{0,\infty,K}
            // first, get the local trace of the velocity into beta

            // local id of the face in its adjacent element
            UInt iFaEl = M_mesh->face( iFace ).pos_first();
            for ( int iNode = 0; iNode < M_feBd.nbNode; ++iNode )
            {
                UInt iloc = M_fToP( iFaEl, iNode+1 );
                for ( int iCoor = 0; iCoor < M_fe1.nbCoor; ++iCoor )
                {
                    UInt ig = M_dof.localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                    if (state.BlockMap().LID(ig + 1) >= 0)
                        beta.vec()[ iCoor*M_feBd.nbNode + iNode ] = state( ig + 1); // BASEINDEX + 1
                }
            }

            // second, calculate its max norm
            for ( int l = 0; l < int( M_fe1.nbCoor*M_feBd.nbNode ); ++l )
            {
                if ( bmax < fabs( beta.vec()[ l ] ) )
                    bmax = fabs( beta.vec()[ l ] );
            }

            chronoBeta.stop();
        }


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
            ipstab_grad( coeffPress, M_elMatP, M_fe1, M_fe1, M_feBd,
                        nDimensions, nDimensions);
//            M_elMatP.showMe();
            chronoElemComp.stop();
            chronoAssembly1.start();
//             assemb_mat(matrix, M_elMatP, M_fe1, M_dof,
//                        nDimensions, nDimensions);
            assembleMatrix(matrix, M_elMatP, M_fe1, M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly1.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad( coeffPress, M_elMatP, M_fe2, M_fe2, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly2.start();
//             assemb_mat(matrix, M_elMatP, M_fe2, M_dof,
//                        nDimensions, nDimensions);
            assembleMatrix(matrix, M_elMatP, M_fe2, M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly2.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(- coeffPress, M_elMatP, M_fe1, M_fe2, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly3.start();
//             assemb_mat(matrix, M_elMatP, M_fe1, M_fe2, M_dof,
//                        nDimensions, nDimensions);
            assembleMatrix(matrix, M_elMatP, M_fe1, M_fe2, M_dof, M_dof,
                       nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly3.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            ipstab_grad(- coeffPress, M_elMatP, M_fe2, M_fe1, M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly4.start();
//             assemb_mat(matrix, M_elMatP, M_fe2, M_fe1, M_dof,
//                        nDimensions, nDimensions);
            assembleMatrix(matrix, M_elMatP, M_fe2, M_fe1, M_dof, M_dof,
                       nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly4.stop();
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
            chronoAssembly5.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    {
//                         assemb_mat( matrix, M_elMatU, M_fe1, M_dof,
//                                     iComp, jComp );
                        assembleMatrix( matrix, M_elMatU, M_fe1, M_dof,
                                         iComp, jComp, iComp*nDof, jComp*nDof );
                    }
            chronoAssembly5.stop();

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
            chronoAssembly6.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                    {
//                         assemb_mat( matrix, M_elMatU, M_fe2, M_dof,
//                                     iComp, jComp );
                        assembleMatrix( matrix, M_elMatU, M_fe2, M_dof,
                                        iComp, jComp, iComp*nDof, jComp*nDof );
                    }
            chronoAssembly6.stop();

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
            chronoAssembly7.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                {
//                         assemb_mat( matrix, M_elMatU, M_fe1, M_fe2, M_dof,
//                                 iComp, jComp );
                    assembleMatrix( matrix, M_elMatU, M_fe1, M_fe2, M_dof, M_dof,
                                    iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly7.stop();

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
            chronoAssembly8.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                {
//                     assemb_mat( matrix, M_elMatU, M_fe2, M_fe1, M_dof,
//                                 iComp, jComp );
                    assembleMatrix( matrix, M_elMatU, M_fe2, M_fe1, M_dof, M_dof,
                                     iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly8.stop();
        }

    } // loop on interior faces
    chronoAssembly.stop();
    if (verbose)
        {
            std::cout << std::endl;
            std::cout << state.BlockMap().Comm().MyPID()
                      <<  "  .   Updating of element   done in "
                      << chronoUpdate.diff_cumul()   << " s." << std::endl;
            std::cout << "   .   Determination of beta done in "
                      << chronoBeta.diff_cumul()     << " s." << std::endl;
            std::cout << "   .   Element computations  done in "
                      << chronoElemComp.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 1              done in "
                      << chronoAssembly1.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 2              done in "
                      << chronoAssembly2.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 3              done in "
                      << chronoAssembly3.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 4              done in "
                      << chronoAssembly4.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 5              done in "
                      << chronoAssembly5.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 6              done in "
                      << chronoAssembly6.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 7              done in "
                      << chronoAssembly7.diff_cumul() << " s." << std::endl;
            std::cout << "   .   chrono 8              done in "
                      << chronoAssembly8.diff_cumul() << " s." << std::endl;
            std::cout << "   .   total                                   "
                      << chronoAssembly.diff_cumul() << " s."
                      << " myFaces = " << myFaces << std::endl;
        }

} // apply(...)

template<typename MESH, typename DOF>
template<typename VECTOR>
void IPStabilization<MESH, DOF>::apply_expl(VECTOR& vector, const VECTOR& state)
{

    // __note__: we shall not use M_gammaDiv for the moment 
    //           take it to zero in the data file
    if (M_gammaBeta == 0. && M_gammaDiv == 0.)
        return;
       
    Chrono chronoAssembly;

    const UInt nDof = M_dof.numTotalDof();
	
    double normInf;

    state.NormInf(&normInf);
    //VECTOR _vector(vector, Repeated, Zero);
    //VECTOR _vector(vector, Repeated);
    VECTOR _vector( vector.getMap(), Repeated );
    _vector *= 0.;
    bool _verb = (_vector.Comm().MyPID() == 0);

    // local trace of the velocity on a face element
    ElemVec beta(M_feBd.nbNode, nDimensions);
    // local trace of the velocity on the volumic elements sharing a face
    ElemVec uLoc1(M_fe1.nbNode, nDimensions);
    ElemVec uLoc2(M_fe2.nbNode, nDimensions);

    UInt myFaces(0);

    chronoAssembly.start();
    // loop on interior faces
    for ( UInt iFace = M_mesh->numBFaces() + 1; iFace<= M_mesh->numFaces();
          ++iFace )
    {
        const UInt iElAd1 = M_mesh->face( iFace ).ad_first();
        const UInt iElAd2 = M_mesh->face( iFace ).ad_second();

        if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
        {
            continue;
        }
        ++myFaces;

        // update current finite elements
        M_feBd.updateMeasNormal( M_mesh->face( iFace ) );
        KNM<Real>& normal = M_feBd.normal;
        const Real hK2 = M_feBd.measure();


        M_fe1.updateFirstDeriv( M_mesh->volumeList( iElAd1 ) );
        M_fe2.updateFirstDeriv( M_mesh->volumeList( iElAd2 ) );

        // __note__: in Burman and Fernandez, they take \beta . \normal on each 
	//           interior face, here we will follow Winkelman and take the
	// max of \beta . \normal on each volumic element. 
        Real bmax(0);
        if (normInf != 0.)
        {
            // determine bmax = ||\beta||_{0,\infty,K}
            
	    // first, get the local trace of the velocity into beta
            // local id of the face in its adjacent element

            UInt iFaEl = M_mesh->face( iFace ).pos_first();
            for ( int iNode = 0; iNode < M_feBd.nbNode; ++iNode )
            {
                UInt iloc = M_fToP( iFaEl, iNode+1 );
                for ( int iCoor = 0; iCoor < M_fe1.nbCoor; ++iCoor )
                {
                    UInt ig = M_dof.localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                    if (state.BlockMap().LID(ig + 1) >= 0)
                        beta.vec()[ iCoor*M_feBd.nbNode + iNode ] = state( ig + 1); // BASEINDEX + 1
                }
            }

            // second, calculate its max norm
            for ( int l = 0; l < int( M_fe1.nbCoor*M_feBd.nbNode ); ++l )
            {
                if ( bmax < fabs( beta.vec()[ l ] ) )
                    bmax = fabs( beta.vec()[ l ] );
            }

        }

        // get local velocity on elements sharing the current face
	// assume local pattern is full to avoid dereference

	for (int iNode = 0; iNode < M_fe1.nbNode; ++iNode)
	{
	    for (int iCoor = 0; iCoor < nDimensions; ++iCoor)
	    {
	    	UInt eleID1 = M_fe1.currentLocalId();
		UInt eleID2 = M_fe2.currentLocalId();
		UInt ig1 = M_dof.localToGlobal( iElAd1, iNode+1) + iCoor * nDof;
		UInt ig2 = M_dof.localToGlobal( iElAd2, iNode+1) + iCoor * nDof;
		uLoc1.vec()[ iNode + iCoor * M_fe1.nbNode ] = state(ig1);
		uLoc2.vec()[ iNode + iCoor * M_fe2.nbNode ] = state(ig2);
	    }
	}

        // velocity stabilization
        if ( ( M_gammaDiv != 0 || M_gammaBeta != 0 ) && bmax > 0 )
        {
            // determine bnmax = ||\beta \cdot n||_{0,\infty,K}
            // and       bcmax = ||\beta \cross n||_{0,\infty,K}

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

	    // temporary solution as this cross term explodes in explicit
	    bcmax = 0.0;

            // put a minus sign here since this is explicit/rhs vector
	    Real coeffGrad = -hK2 * (M_gammaBeta*bnmax + M_gammaDiv*bcmax);

            // compute <grad u+ | grad v+> term
	    M_elvec_u.zero();
	    ipstab_grad_expl(coeffGrad, uLoc1, M_elvec_u, M_fe1, M_fe1, M_feBd);

	    for (UInt iComp = 0; iComp < nDimensions; ++iComp)
		assembleVector(_vector, M_elvec_u, M_fe1, M_dof, iComp, iComp*nDof);
	    
            // compute <grad u- | grad v-> term 
            M_elvec_u.zero();
	    ipstab_grad_expl(coeffGrad, uLoc2, M_elvec_u, M_fe2, M_fe2, M_feBd);

            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
		assembleVector(_vector, M_elvec_u, M_fe2, M_dof, iComp, iComp*nDof);

            // compute -<grad u+ | grad v-> term
	    M_elvec_u.zero();
	    ipstab_grad_expl(-coeffGrad, uLoc1, M_elvec_u, M_fe1, M_fe2, M_feBd);
            
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                assembleVector(_vector, M_elvec_u, M_fe2, M_dof, iComp, iComp*nDof);

            // compute -<grad u- | grad v+> term
            M_elvec_u.zero();
            ipstab_grad_expl( -coeffGrad, uLoc2, M_elvec_u, M_fe2, M_fe1, M_feBd);
            
	    for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                assembleVector(_vector, M_elvec_u, M_fe1, M_dof, iComp, iComp*nDof);
	}

    } // loop on interior faces

    _vector.Comm().Barrier();

//SDEB
    Real state_norm;
    Real stab_norm;
    state.NormInf(&state_norm);
    _vector.NormInf(&stab_norm);
//    std::cout << "\n\n\n\n *** State norm = " << state_norm << "\n *** Stab norm  = " 
//        << stab_norm << "\n\n\n\n" << std::endl;
    if (_verb)
        std::cout << "  ip-  Relative stabilization norm: " << 
	    stab_norm/state_norm << std::endl;
//EDEB

    vector += _vector;
       
    chronoAssembly.stop();

    std::cout << "  ip-  IP terms assembling: " << 
       chronoAssembly.diff_cumul() << " sec." << std::endl;

}  // apply_expl



//} // namespace details

} // namespace LifeV

#endif /* _NSIPTERMS_HPP */
