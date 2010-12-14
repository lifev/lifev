//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Interior Penality Stabilization

    @author Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @contributor Umberto Villa <uvilla@emory.edu>
    @maintainer Umberto Villa <uvilla@emory.edu>

    @date 08-10-2004

    Implementation of I.P. stabilization for inf-sup incompatible finite elements for the Navier-Stokes equations.
 */

//TODO Change name to StabilizationIP

#ifndef _NSIPTERMS_HPP
#define _NSIPTERMS_HPP

#include <life/lifecore/chrono.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <boost/shared_ptr.hpp>

#define USE_OLD_PARAMETERS 0
#define WITH_DIVERGENCE 1

namespace LifeV
{

namespace details
{
//! IPStabilization Class
/*!
 * @brief Interior Penality Stabilization
 * @author C. Winkelmann <christoph.winkelmann@epfl.ch>
 *
 * Implementation of I.P. stabilization for inf-sup incompatible finite elements
 * for the Navier-Stokes equations. <br>
 *
 *  This function adds the following stabilization terms to the Navier-Stokes monolithic matrix:
 *  <ol>
 *  <li> Block(1,1): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_\beta [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}] + \Sigma_{f\in\mathcal{F}}\int_{f} c_d [div \mathbf{u}] [div \mathbf{v}]@f$
 *  <li> Block(1,2): 0
 *  <li> Block(2,1): 0
 *  <li> Block(2,2): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_p [\nabla p] \cdot [\nabla q]@f$
 *  </ol>
 *
 *  where
 *  <ol>
 *  <li> @f$\mathcal{F}@f$ is the set of all the internal facets in the mesh;
 *  <li> @f$[\cdot]@f$ is the jump operator: @f$[x] = x^+ - x^-@f$.
 *  </ol>
 *  and the parameter @f$ c_\beta@f$, @f$c_d@f$, @f$c_p@f$, are defined as follows:
 *  <ol>
 *  <li> @f$\displaystyle c_\beta = \frac{\gamma_\beta h}  {\|\beta\|_\infty} @f$, being @f$h@f$ the facet measure;
 *  <li> @f$\displaystyle c_d = \gamma_d h \|\beta\|_\infty @f$;
 *  <li> @f$\displaystyle c_p = \frac{\gamma_p h}{\max(\|\beta\|_\infty, \nu/h^2)}@f$.
 *  </ol>
 *  Both high Pechlet numbers and inf-sup incompatible FEM are stabilized.
 *
 */

template<typename MeshType, typename DofType>
class IPStabilization
{
public:

    //! @name Public Types
    //@{
    typedef boost::shared_ptr<MeshType> mesh_type; //deprecated
    typedef MeshType  mesh_Type;
    typedef DofType   dof_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;
    typedef boost::shared_ptr<dof_Type>  dofPtr_Type;
    //@}

    //! @name Constructor and Destructor
    //@{

    //! Default Constructor
    IPStabilization() {};

    //! Constructor
    /*!
     * @param mesh       meshPtr_Type  a pointer to the mesh
     * @param dof        dof_Type      the velocity/pressure dof (same order discretization)
     * @param refFE      RefFE         the velocity/pressure field reference finite
     * @param feBd       CurrentBdFE   the facet current fe to be used to compute the jumps
     * 								   on the interface within elements.
     * @param quadRule   QuadRule      the element quadrature rule used for the facet projection
     * @param gammaBeta  Real 		   the stabilization parameter @f$\gamma_\beta@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}]@f$
     * @param gammaDiv   Real          the stabilization parameter @f$\gamma_d@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [div \mathbf{u}] [div \mathbf{v}]@f$
     * @param gammaPress Real          the stabilization parameter @f$\gamma_p@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\nabla p] \cdot [\nabla q]@f$
     * @param viscosity  Real          the fluid viscosity @f$\nu@f$
     */
    IPStabilization( const meshPtr_Type&     mesh,
                     const dof_Type&        dof,
                     const RefFE&           refFE,
                     CurrentBdFE&           feBd,
                     const QuadRule&        quadRule,
                     Real                   gammaBeta = 0,
                     Real                   gammaDiv = 0,
                     Real                   gammaPress = 0,
                     Real                   viscosity = 1 );

    virtual ~IPStabilization() {};
    //@}

    //! @name Methods
    //@{
    //! compute IP stabilization terms and add them into matrix
    /*!
     *  This function adds the following stabilization terms to the Navier-Stokes monolithic matrix:
     *  <ol>
     *  <li> Block(1,1): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_\beta [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}] + \Sigma_{f\in\mathcal{F}}\int_{f} c_d [div \mathbf{u}] [div \mathbf{v}]@f$
     *  <li> Block(1,2): 0
     *  <li> Block(2,1): 0
     *  <li> Block(2,2): @f$\Sigma_{f\in\mathcal{F}}\int_{f} c_p [\nabla p] \cdot [\nabla q]@f$
     *  </ol>
     *  where
     *  <ol>
     *  <li> @f$\displaystyle c_\beta = \frac{\gamma_\beta h}  {\|\beta\|_\infty} @f$, being @f$h@f$ the facet measure;
     *  <li> @f$\displaystyle c_d = \gamma_\beta h \|\beta\|_\infty @f$;
     *  <li> @f$\displaystyle c_p = \frac{\gamma_p h}{\max(\|\beta\|_\infty, \nu/h^2)}@f$.
     *  </ol>
     *  Both high Pechlet numbers and inf-sup incompatible FEM are stabilized.
     *
     *	PREREQUISITE: The velocity and the pressure field should belong to the same finite element space
     *
     *  Parameters are the followings:
     *  @param dt      Real   timestep (INPUT)
     *  @param matrix  MatrixType where the stabilization terms are added into. (OUTPUT)
     *  @param state   VectorType velocity field for the linearization of the stabilization (INPUT)
     *  @param verbose bool   whenever of not to print on screen
     */
    template<typename MatrixType, typename VectorType>
    void apply( MatrixType& matrix, const VectorType& state, bool verbose = true );

    //! Display class informations
    /*!
     * Write information relative to the class on output
     * @param output ostream ostream were to write (Default cout)
     */
    void showMe(std::ostream & output = std::cout) const;
    //@}

    //! @name Set Methods
    //@{
    //! Set the stabilization parameter @f$\gamma_\beta@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}]@f$
    void setGammaBeta (const Real & gammaBeta) { M_gammaBeta  = gammaBeta;}
    //! Set the stabilization parameter @f$\gamma_d@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [div \mathbf{u}] [div \mathbf{v}]@f$
    void setGammaDiv  (const Real & gammaDiv)  { M_gammaDiv   = gammaDiv;}
    //! Set the stabilization parameter @f$\gamma_p@f$ for @f$\Sigma_{f\in\mathcal{F}}\int_{f} [\nabla p] \cdot [\nabla q]@f$
    void setGammaPress(const Real & gammaPress) { M_gammaPress = gammaPress;}
    //! Set the fluid viscosity @f$\nu@f$
    void setViscosity(const Real & viscosity) { M_viscosity = viscosity;}
    //! Set the mesh file
    void setMesh(const meshPtr_Type mesh) { M_mesh = mesh; }
    //! Set Discretization
    void setDiscretization(const dofPtr_Type& dof, const RefFE& refFE, CurrentBdFE& feBd, const QuadRule& quadRule);
    //@}
private:

    //! @name Private Types
    //@{
    //! fToP(i,j) = localId of jth point on ith local face
    typedef ID ( *FTOP )( ID const localFace, ID const point );
    //@}

    //! @name Private Constructor
    //@{
    //! Copy Constructor
    IPStabilization(const IPStabilization<mesh_Type, dof_Type> & original);
    //@}

    //! @name Private Attributes
    //@{
    //! Pointer to the mesh object
    const meshPtr_Type  M_mesh;
    //! reference to the DofType data structure
    const dofPtr_Type   M_dof;
    //! current Fe on side 1 of the current face
    boost::shared_ptr<CurrentFE>    M_feOnSide1;
    //! current Fe on side 2 of the current face
    boost::shared_ptr<CurrentFE>    M_feOnSide2;
    //! current boundary FE
    CurrentBdFE*  M_feBd;
    //! Stabilization parameter @f$\gamma_\beta@f$ for @f$\int_{face} [\beta \cdot \nabla \mathbf{u}] [\beta \cdot \nabla \mathbf{v}]@f$
    Real         M_gammaBeta;
    //! Stabilization parameter @f$\gamma_d@f$ for @f$\int_{face} [div \mathbf{u}] [div \mathbf{v}]@f$
    Real         M_gammaDiv;
    //! Stabilization parameter @f$\gamma_p@f$ for @f$\int_{face} [\nabla p] \cdot [\nabla q]@f$
    Real         M_gammaPress;
    //! Fluid viscosity @f$\nu@f$
    Real         M_viscosity;
    //! Elementary matrix
    ElemMat      M_elMatU;
    //! Elementary matrix
    ElemMat      M_elMatP;
    //! fToP(i,j) = localId of jth point on ith local face
    FTOP         M_fToP;
    //@}
}; // class IPStabilization


//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename DofType>
IPStabilization<MeshType, DofType>::IPStabilization( const meshPtr_Type & mesh,
                                                     const dof_Type&      dof,
                                                     const RefFE&    refFE,
                                                     CurrentBdFE&    feBd,
                                                     const QuadRule& quadRule,
                                                     Real            gammaBeta,
                                                     Real            gammaDiv,
                                                     Real            gammaPress,
                                                     Real            viscosity ) :
        M_mesh      ( mesh ),
        M_dof       ( new dof_Type(dof) ),
        M_feOnSide1 ( new CurrentFE(refFE, getGeoMap(*mesh), quadRule) ),
        M_feOnSide2 ( new CurrentFE(refFE, getGeoMap(*mesh), quadRule) ),
        M_feBd      ( &feBd ),
        M_gammaBeta ( gammaBeta ),
        M_gammaDiv  ( gammaDiv ),
        M_gammaPress( gammaPress ),
        M_viscosity ( viscosity ),
        M_elMatU    ( M_feOnSide1->nbFEDof(), nDimensions    , nDimensions   ),
        M_elMatP    ( M_feOnSide1->nbFEDof(), nDimensions + 1, nDimensions+1 )
{
    switch ( M_feOnSide1->nbFEDof() )
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

//=============================================================================
// Method
//=============================================================================

template<typename MeshType, typename DofType>
template<typename MatrixType, typename VectorType>
void IPStabilization<MeshType, DofType>::apply( MatrixType& matrix,  const VectorType& state, const bool verbose )
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

    const UInt nDof = M_dof->numTotalDof();

    Real normInf;
    state.NormInf(&normInf);

    // local trace of the velocity
    ElemVec beta( M_feBd->nbNode, nDimensions );

    UInt myFaces(0);

    chronoAssembly.start();
    // loop on interior faces
    for ( UInt iFace( M_mesh->numBFaces() + 1 ); iFace<= M_mesh->numFaces();
            ++iFace )
    {
        const UInt iElAd1 ( M_mesh->face( iFace ).ad_first()  );
        const UInt iElAd2 ( M_mesh->face( iFace ).ad_second() );

        if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
        {
            //std::cout << "iElAd1 = " << iElAd1 << "; iElAd2 = " << iElAd2 << std::endl;
            continue;
        }
        ++myFaces;

        chronoUpdate.start();
        // update current finite elements
#if WITH_DIVERGENCE
        M_feBd->updateMeas( M_mesh->face( iFace ) );
#else
        M_feBd->updateMeasNormal( M_mesh->face( iFace ) );
        KNM<Real>& normal = M_feBd->normal;
#endif
        const Real hK2 = M_feBd->measure();


        M_feOnSide1->updateFirstDeriv( M_mesh->volumeList( iElAd1 ) );
        M_feOnSide2->updateFirstDeriv( M_mesh->volumeList( iElAd2 ) );
        chronoUpdate.stop();

        Real bmax(0);
        if (normInf != 0.)
        {
            chronoBeta.start();
            // determine bmax = ||\beta||_{0,\infty,K}
            // first, get the local trace of the velocity into beta

            // local id of the face in its adjacent element
            UInt iFaEl ( M_mesh->face( iFace ).pos_first() );
            for ( UInt iNode ( 0 ); iNode < M_feBd->nbNode; ++iNode )
            {
                UInt iloc ( M_fToP( iFaEl, iNode+1 ) );
                for ( UInt iCoor ( 0 ); iCoor < M_feOnSide1->nbCoor(); ++iCoor )
                {
                    UInt ig ( M_dof->localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof );
                    if (state.BlockMap().LID(ig + 1) >= 0)
                        beta.vec()[ iCoor*M_feBd->nbNode + iNode ] = state( ig + 1); // BASEINDEX + 1
                }
            }

            // second, calculate its max norm
            for ( UInt l ( 0 ); l < static_cast<UInt>( M_feOnSide1->nbCoor()*M_feBd->nbNode ); ++l )
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
            Real coeffPress ( M_gammaPress * hK2 ); // P1, P2 (code)
            //Real coeffPress = M_gammaPress * sqrt( hK2 ); // P1 p nonsmooth (code)
#else
            Real coeffPress = M_gammaPress * hK2 / // Pk (paper)
                              std::max<Real>( bmax, M_viscosity/sqrt( hK2 ) );
#endif

            M_elMatP.zero();
            chronoElemComp.start();
            // coef*\int_{face} grad u1 . grad v1
            ipstab_grad( coeffPress, M_elMatP, *M_feOnSide1, *M_feOnSide1, *M_feBd,
                         nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly1.start();

            assembleMatrix(matrix, M_elMatP, *M_feOnSide1, *M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly1.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            // coef*\int_{face} grad u2 . grad v2
            ipstab_grad( coeffPress, M_elMatP, *M_feOnSide2, *M_feOnSide2, *M_feBd,
                         nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly2.start();

            assembleMatrix(matrix, M_elMatP, *M_feOnSide2, *M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly2.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            // - coef*\int_{face} grad u1 . grad v2
            ipstab_grad(- coeffPress, M_elMatP, *M_feOnSide1, *M_feOnSide2, *M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly3.start();

            assembleMatrix(matrix, M_elMatP, *M_feOnSide1, *M_feOnSide2, *M_dof, *M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly3.stop();

            M_elMatP.zero();
            chronoElemComp.start();
            // - coef*\int_{face} grad u2 . grad v1
            ipstab_grad(- coeffPress, M_elMatP, *M_feOnSide2, *M_feOnSide1, *M_feBd,
                        nDimensions, nDimensions);
            chronoElemComp.stop();
            chronoAssembly4.start();

            assembleMatrix(matrix, M_elMatP, *M_feOnSide2, *M_feOnSide1, *M_dof, *M_dof,
                           nDimensions, nDimensions, nDimensions*nDof, nDimensions*nDof);
            chronoAssembly4.stop();
        }

        // velocity stabilization
        if ( ( M_gammaDiv != 0 || M_gammaBeta != 0 ) && bmax > 0 )
        {
#if WITH_DIVERGENCE
#if USE_OLD_PARAMETERS
            Real coeffBeta ( M_gammaBeta * hK2 / std::max<Real>(bmax, hK2) ); // code
#else
            Real coeffBeta ( M_gammaBeta * hK2 / bmax ); // paper
#endif

            Real coeffDiv ( M_gammaDiv * hK2 * bmax ); // (code and paper)
            //Real coeffDiv ( M_gammaDiv * sqrt( hK2 ) * bmax ); // ? (code)
#else
            // determine bnmax = ||\beta \cdot n||_{0,\infty,K}
            // and       bcmax = ||\beta \cross n||_{0,\infty,K}

            chronoBeta.start();
            Real bnmax ( 0. );
            Real bcmax ( 0. );
            for ( UInt iNode(0); iNode<M_feBd->nbNode; ++iNode )
            {
                Real bn ( 0 );
                for ( UInt iCoor(0); iCoor<M_feOnSide1->nbCoor(); ++iCoor )
                {
                    bn += normal(iNode, iCoor) *
                          beta.vec()[ iCoor*M_feBd->nbNode + iNode ];
                    bcmax = std::max<Real>
                            (bcmax, normal(iNode, (iCoor+1)%3) *
                             beta.vec()[ (iCoor+2)%3*M_feBd->nbNode + iNode ] -
                             normal(iNode, (iCoor+2)%3) *
                             beta.vec()[ (iCoor+1)%3*M_feBd->nbNode + iNode ]);
                }
                bnmax = std::max<Real> (bnmax, bn);
            }
            chronoBeta.stop();

            Real coeffGrad = hK2 * (M_gammaBeta*bnmax + M_gammaDiv*bcmax);
#endif
            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // coef*\int_{face} (\beta1 . grad u1) (\beta2 . grad v2)
            ipstab_bgrad( coeffBeta, M_elMatU, *M_feOnSide1, *M_feOnSide1, beta,
                          *M_feBd, 0, 0, nDimensions );
            // coef*\int_{face} div u1 . div v1
            ipstab_div( coeffDiv, M_elMatU, *M_feOnSide1, *M_feOnSide1, *M_feBd );
#else
            // coef*\int_{face} grad u1 . grad v1
            ipstab_grad( coeffGrad, M_elMatU, *M_feOnSide1, *M_feOnSide1, *M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly5.start();
            for ( UInt iComp ( 0 ); iComp<nDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp<nDimensions; ++jComp )
                {
                    assembleMatrix( matrix, M_elMatU, *M_feOnSide1, *M_dof,
                                    iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly5.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // coef*\int_{face} (\beta2 . grad u2) (\beta2 . grad v2)
            ipstab_bgrad( coeffBeta, M_elMatU, *M_feOnSide2, *M_feOnSide2, beta,
                          *M_feBd, 0, 0, nDimensions );
            // coef*\int_{face} div u2 . div v2
            ipstab_div( coeffDiv, M_elMatU, *M_feOnSide2, *M_feOnSide2, *M_feBd );
#else
            // coef*\int_{face} grad u2 . grad v2
            ipstab_grad( coeffGrad, M_elMatU, *M_feOnSide2, *M_feOnSide2, *M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly6.start();
            for ( UInt iComp ( 0 ); iComp<nDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp<nDimensions; ++jComp )
                {
                    assembleMatrix( matrix, M_elMatU, *M_feOnSide2, *M_dof,
                                    iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly6.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // - coef*\int_{face} (\beta1 . grad u1) (\beta2 . grad v2)
            ipstab_bgrad( -coeffBeta, M_elMatU, *M_feOnSide1, *M_feOnSide2, beta,
                          *M_feBd, 0, 0, nDimensions );
            // - coef*\int_{face} div u1 . div v2
            ipstab_div( -coeffDiv, M_elMatU, *M_feOnSide1, *M_feOnSide2, *M_feBd );
#else
            // - coef*\int_{face} grad u1 . grad v2
            ipstab_grad( -coeffGrad, M_elMatU, *M_feOnSide1, *M_feOnSide2, *M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly7.start();
            for ( UInt iComp = 0; iComp<nDimensions; ++iComp )
                for ( UInt jComp = 0; jComp<nDimensions; ++jComp )
                {
                    assembleMatrix( matrix, M_elMatU, *M_feOnSide1, *M_feOnSide2, *M_dof, *M_dof,
                                    iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly7.stop();

            M_elMatU.zero();
            chronoElemComp.start();
#if WITH_DIVERGENCE
            // - coef*\int_{face} (\beta2 . grad u2) (\beta1 . grad v1)
            ipstab_bgrad( -coeffBeta, M_elMatU, *M_feOnSide2, *M_feOnSide1, beta,
                          *M_feBd, 0, 0, nDimensions );
            // - coef*\int_{face} div u2 . div v1
            ipstab_div( -coeffDiv, M_elMatU, *M_feOnSide2, *M_feOnSide1, *M_feBd );
#else
            // - coef*\int_{face} grad u2 . grad v1
            ipstab_grad( -coeffGrad, M_elMatU, *M_feOnSide2, *M_feOnSide1, *M_feBd, 0, 0,
                         nDimensions );
#endif
            chronoElemComp.stop();
            chronoAssembly8.start();
            for ( UInt iComp ( 0 ); iComp<nDimensions; ++iComp )
                for ( UInt jComp ( 0 ); jComp<nDimensions; ++jComp )
                {
                    assembleMatrix( matrix, M_elMatU, *M_feOnSide2, *M_feOnSide1, *M_dof, *M_dof,
                                    iComp, jComp, iComp*nDof, jComp*nDof );
                }
            chronoAssembly8.stop();
        }

    } // loop on interior faces
    chronoAssembly.stop();
    if (verbose)
    {
        Debug(7101) << "\n";
        Debug(7101) << static_cast<UInt>(state.BlockMap().Comm().MyPID())
        <<  "  .   Updating of element   done in "
        << chronoUpdate.diff_cumul()   << " s." << "\n";
        Debug(7101) << "   .   Determination of beta done in "
        << chronoBeta.diff_cumul()     << " s." << "\n";
        Debug(7101) << "   .   Element computations  done in "
        << chronoElemComp.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 1              done in "
        << chronoAssembly1.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 2              done in "
        << chronoAssembly2.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 3              done in "
        << chronoAssembly3.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 4              done in "
        << chronoAssembly4.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 5              done in "
        << chronoAssembly5.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 6              done in "
        << chronoAssembly6.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 7              done in "
        << chronoAssembly7.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   chrono 8              done in "
        << chronoAssembly8.diff_cumul() << " s." << "\n";
        Debug(7101) << "   .   total                                   "
        << chronoAssembly.diff_cumul() << " s."
        << " myFaces = " << myFaces << "\n";
    }

} // apply(...)

template<typename MeshType, typename DofType>
void IPStabilization<MeshType, DofType>::showMe(std::ostream & output) const
{
    output << "IPStabilization::showMe() " <<std::endl;
    output << "Fluid Viscosity: " << M_viscosity << std::endl;
    output << "Stabilization coefficient velocity SD jumps:         " << M_gammaBeta  << std::endl;
    output << "Stabilization coefficient velocity divergence jumps: " << M_gammaDiv   << std::endl;
    output << "Stabilization coefficient pressure gradient jumps:   " << M_gammaPress << std::endl;
    M_mesh->showMe(output);
    M_dof->showMe(output);
}

//=============================================================================
// Setters method
//=============================================================================
template<typename MeshType, typename DofType>
void IPStabilization<MeshType, DofType>::setDiscretization(const dofPtr_Type& dof, const RefFE& refFE, CurrentBdFE& feBd, const QuadRule& quadRule)
{
    M_dof = dof;
    M_feOnSide1.reset( new CurrentFE(refFE, getGeoMap(*M_mesh), quadRule) );
    M_feOnSide2.reset( new CurrentFE(refFE, getGeoMap(*M_mesh), quadRule) );
    M_feBd      ( &feBd );
    switch ( M_feOnSide1->nbFEDof() )
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

} // namespace details

} // namespace LifeV

#endif /* _NSIPTERMS_HPP */
