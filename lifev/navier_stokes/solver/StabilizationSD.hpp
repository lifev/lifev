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
   @brief Streamline diffusion and SUPG stabilization.
   @author M.A. Fernandez
   @contributor Umberto Villa <uvilla@emory.edu>
   @maintainer Umberto Villa <uvilla@emory.edu>
   @date 01/05/2005

   This file contains a c++ class implementing the streamline diffusion (SD) and the SUPG
   stabilization for the Navier-Stokes equations. Tested with P1/P1 and Q1/Q1

*/

#ifndef _SDSTABILIZATION_HPP_
#define _SDSTABILIZATION_HPP_


#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/LifeDebug.hpp>

namespace LifeV
{

//! StabilizationSD Class
/*!
 * @brief Streamline diffusion and SUPG for Navier-Stokes
 * @author M.A. Fernandez
 *
 * Implementation of streamline diffusion (SD) and SUPG for the Navier-Stokes equations. <br>
 * SUPG can be used only with equal order finite element for the discretization of velocity and pressure fields.
 */

template<typename MeshType, typename DofType>
class StabilizationSD
{
public:

    //@name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef DofType  dof_Type;
    //@}

    //! @name Constructor and Destructor
    //@{
    //! Constructor
    /*!
     * @param dataFile GetPot       dataFile where to read the stabilization parameter
     * @param mesh     mesh_Type    mesh
     * @param dof      dof_Type     velocity field degree of freedom
     * @param refFE    refFE        velocity field reference finite element
     * @param quadRule QuadratureRule     quadrature rule for the integration of the stabilization variational forms
     * @param viscosity viscosity   fluid viscosity  @f$\nu@f$
     */
    StabilizationSD( const GetPot&   dataFile,
                     const mesh_Type&     mesh,
                     const dof_Type&      dof,
                     const ReferenceFE&    refFE,
                     const QuadratureRule& quadRule,
                     Real      viscosity);
    //! ~Destructor
    virtual ~StabilizationSD() {};
    //@}

    //! @name Methods
    //@{
    //! compute the SUPG stabilization terms and add them into the monolithic N.S. matrix
    /*!
     *  This function adds the following stabilization terms to the Navier-Stokes monolithic matrix:
     *  <ol>
     *  <li> Block(1,1): @f$c_\beta(\beta \nabla \mathbf{u} , \beta \nabla \mathbf{v}) + c_\beta( - \mu L \mathbf{u} , \beta \nabla \mathbf{v} )+ c_d( div \mathbf{u} , div \nabla \mathbf{v})@f$
     *  <li> Block(1,2): @f$c_\beta(\nabla p , \beta \nabla \mathbf{v})@f$
     *  <li> Block(2,1): @f$c_\beta(\beta \nabla \mathbf{u} , \nabla q) + c_\beta( - \mu L \mathbf{u} , \nabla q )@f$
     *  <li> Block(2,2): @f$c_\beta(\nabla p , \nabla q)@f$
     *  </ol>
     *  where @f$(\cdot, \cdot)@f$ represents the @f$L^2@f$ scalar product, and
     *  <ol>
     *  <li> @f$\displaystyle c_\beta = \frac{\gamma_\beta}{\sqrt{ 4/( dt^2)
                                     + 4\|\beta\|_\infty/h^2
                                     + 16*\nu/h^4 }} @f$, where @f$h@f$ is the diameter of the element;
     *  <li> @f$\displaystyle c_d = \gamma_d h \|\beta\|_\infty @f$.
     *  </ol>
     *  Both high Pechlet numbers and inf-sup incompatible FEM are stabilized.
     *
     *	PREREQUISITE: The velocity and the pressure field should belong to the same finite element space
     *
     *  Parameters are the followings:
     *  @param dt     Real   timestep (INPUT)
     *  @param matrix MatrixType where the stabilization terms are added into. (OUTPUT)
     *  @param state  VectorType velocity field for the linearization of the stabilization (INPUT)
     */
    template<typename MatrixType, typename VectorType>
    void applySUPG(const Real dt, MatrixType& matrix, const VectorType& state );

    //! compute the SD stabilization terms and add them into the Momentum matrix
    /*!
     *  The following stabilization term is added:
     *  @f$c_\beta(\beta \nabla \mathbf{u} , \beta \nabla \mathbf{v})@f$
     *  @param dt     Real   timestep
     *  @param matrix MatrixType where the stabilization terms are added into
     *  @param state  VectorType velocity field for the linearization of the stabilization
     */
    template<typename MatrixType, typename VectorType>
    void applySD(const Real dt, MatrixType& matrix, const VectorType& state );

    //! compute the SUPG stabilization terms and add them into the right and side
    /*!
     *  Add the following stabilization terms to rhs
     *  <ol>
     *  <li> Block(1,1):  @f$(c_\beta f , \beta \nabla \mathbf{v})@f$;
     *  <li> Block(1,2):  @f$(c_\beta f , \nabla q )@f$.
     *  </ol>
     *
     *  @param dt     Real   timestep
     *  @param vector VectorType where the stabilization terms are added into
     *  @param state  VectorType velocity field for the linearization of the stabilization
     *  @param source SourceType a functor f(time, x, y, z, ic) which represents the forcing term
     *  @param time   REAL   the actual time in which source should be evaluated
     */
    template <typename VectorType, typename SourceType >
    void applyRHS(const Real dt, VectorType& vector, const VectorType& state,
                  const SourceType& source, const Real& time);

    //! Display class informations
    void showMe(std::ostream & output = std::cout) const;
    //@}

    //! @name Set Methods
    //@{
    //! Set the stabilization parameter @f$\gamma_\beta@f$ for @f$( \beta \nabla \mathbf{u} , \beta \nabla \mathbf{v})@f$
    void setGammaBeta (const Real & gammaBeta) { M_gammaBeta  = gammaBeta;}
    //! Set the stabilization parameter @f$\gamma_d@f$ for @f$( div \mathbf{u} , div \nabla \mathbf{v})@f$
    void setGammaDiv  (const Real & gammaDiv)  { M_gammaDiv   = gammaDiv;}
    //@}
private:

    //! @name Private Constructor
    //@{
    //! Default Constructor
    StabilizationSD();
    //! Copy Constructor
    StabilizationSD(const StabilizationSD<mesh_Type, dof_Type> & original);
    //@}

    //! @name Private Methods
    //@{
    // methods for elementary computations

    //! Compute the stabilization coefficients for each element
    template <typename VectorType>
    void computeParameters(const Real dt, const UInt iVol, const VectorType& state,
                           VectorElemental& beta,  Real& coeffBeta, Real& coeffDiv) const;

    //!Evaluate the varf @f$(\beta \nabla \mathbf{u}, \beta \nabla \mathbf{v})@f$
    void bgradu_bgradv(const Real& coef, VectorElemental& vel, ElemMat& elmat,const CurrentFE& fe,
                       UInt iblock, UInt jblock, UInt nb) const;

    //! Evaluate the varf @f$(\Delta \mathbf{u}, \beta \nabla \mathbf{v})@f$
    void lapu_bgradv(const Real& coef, VectorElemental& vel, MatrixElemental& elmat, const CurrentFE& fe,
                     UInt iblock, UInt jblock, UInt nb) const;

    //! Evaluate the varf @f$(\nabla p, \beta \nabla \mathbf{v})@f$
    void gradp_bgradv(const Real& coef, VectorElemental& vel, MatrixElemental& elmat, const CurrentFE& fe) const;

    //! Evaluate the varf @f$(\Delta \mathbf{u}, \nabla q)@f$
    void lapu_gradq(const Real& coef, MatrixElemental& elmat, const CurrentFE& fe) const;

    //! Evaluate the varf @f$(f, \beta \nabla \mathbf{v})@f$
    template <typename SourceType>
    void f_bgradv(const Real& coef, SourceType& source, VectorElemental& vel,
                  VectorElemental& elvec, const CurrentFE& fe, UInt iblock, const Real& time) const;

    //! Evaluate the varf @f$(f, \nabla q)@f$
    template<typename SourceType>
    void f_gradq(const Real& coef, SourceType& source, VectorElemental& elvec,
                 const CurrentFE& fe, UInt iblock, const Real& time) const;
    //@}

    //! @name Private Attributes
    //@{
    //! the mesh object
    const mesh_Type&  M_mesh;
    //! the dof object
    const dof_Type&   M_dof;
    //! current fe for the assembling of stabilization terms.
    CurrentFE    M_fe;
    //! fluid viscosity @f$\nu@f$
    Real         M_viscosity;
    //! Stabilization coefficient of @f$(c(h,dt, |\beta|, nu) \beta \nabla \mathbf{u} , \beta \nabla \mathbf{v})@f$
    Real         M_gammaBeta;
    //! Stabilization coefficient of @f$(c(h,dt) div \mathbf{u} , div \nabla \mathbf{v})@f$
    Real         M_gammaDiv;
    //! Elementary Matrix for assembling the stabilization terms
    MatrixElemental      M_elMat;
    //! Elementary Vector for assembling the stabilization terms
    VectorElemental      M_elVec;
    //@}
}; // class StabilizationSD

//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename DofType>
StabilizationSD<MeshType, DofType>::StabilizationSD( const GetPot& dataFile,
                                                     const mesh_Type&     mesh,
                                                     const dof_Type&      dof,
                                                     const ReferenceFE&    refFE,
                                                     const QuadratureRule& quadRule,
                                                     Real            viscosity):
        M_mesh( mesh ),
        M_dof( dof ),
        M_fe( refFE, getGeometricMap(mesh), quadRule ),
        M_viscosity( viscosity ),
        M_gammaBeta ( dataFile( "fluid/sdstab/gammaBeta", 0. ) ),
        M_gammaDiv  ( dataFile( "fluid/sdstab/gammaDiv", 0. ) ),
        M_elMat( M_fe.nbNode, nDimensions+1, nDimensions+1 ) ,
        M_elVec( M_fe.nbNode, nDimensions+1 ) {}

//=============================================================================
// Methods
//=============================================================================

template<typename MeshType, typename DofType>
template <typename MatrixType, typename VectorType>
void StabilizationSD<MeshType, DofType>::applySUPG(const Real dt, MatrixType& matrix, const VectorType& state )
{
    if ( M_gammaBeta == 0 && M_gammaDiv == 0)
        return;

    Chrono chronoBeta;
    Chrono chronoUpdate;
    Chrono chronoElemComp;
    Chrono chronoAssembly;

    // stabilization parameters
    Real coeffBeta, coeffDiv;

    // local velocity
    VectorElemental beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 0; iVol < M_mesh.numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstSecondDeriv( M_mesh.volumeList( iVol ) );

        // stabilization parameters computation
        chronoBeta.start();
        this->computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
        chronoBeta.stop();

        chronoElemComp.start();
        M_elMat.zero();

        // coeffBeta (\beta \nabla \mathbf{u} , \beta \nabla \mathbf{v})
        //
        this->bgradu_bgradv(coeffBeta, beta, M_elMat, M_fe, 0, 0, nDimensions);

        // coeffBeta  ( (\beta \nabla \mathbf{u} , \nabla q) + (\nabla p , \beta \nabla \mathbf{v}) )
        //
        this->gradp_bgradv(coeffBeta, beta, M_elMat, M_fe);

        // coeffBeta (\nabla p , \nabla q)
        //
        stiff( coeffBeta, M_elMat, M_fe, nDimensions, nDimensions );

        // coeffBeta ( - \mu L \mathbf{u} , \beta \nabla \mathbf{v} )
        //
        this->lapu_bgradv(-coeffBeta*M_viscosity, beta, M_elMat, M_fe, 0, 0, nDimensions);

        // coeffBeta ( - \mu L \mathbf{u} , \nabla q )
        //
        this->lapu_gradq(-coeffBeta*M_viscosity, M_elMat, M_fe);

        // coeffDiv ( div \mathbf{u} , div \nabla \mathbf{v})
        //
        stiff_div(coeffDiv, M_elMat, M_fe );

        chronoElemComp.stop();

        chronoAssembly.start();
        for ( UInt iComp = 0; iComp <= nDimensions; ++iComp )
            for ( UInt jComp = 0; jComp <= nDimensions; ++jComp )
                //assemb_mat( matrix, M_elMat, M_fe, M_dof, iComp, jComp );
                assembleMatrix( matrix, M_elMat, M_fe, M_dof,
                                iComp, jComp,
                                iComp*M_dof.numTotalDof(), jComp*M_dof.numTotalDof() );
        chronoAssembly.stop();


    }// loop on elements

    debugStream(7100) << std::endl;
    debugStream(7100) << "      Updating of element   done in "
    << chronoUpdate.diffCumul()   << "s." << std::endl;
    debugStream(7100) << "      Determination of parameters done in "
    << chronoBeta.diffCumul()     << "s." << std::endl;
    debugStream(7100) << "      Element computations  done in "
    << chronoElemComp.diffCumul() << "s." << std::endl;
    debugStream(7100) << "      Assembly              done in "
    << chronoAssembly.diffCumul() << "s." << std::endl;


} // applySUPG(...)


template<typename MeshType, typename DofType>
template <typename MatrixType, typename VectorType>
void StabilizationSD<MeshType, DofType>::applySD(const Real dt, MatrixType& matrix, const VectorType& state )
{
    if ( M_gammaBeta == 0 && M_gammaDiv == 0)
        return;

    Chrono chronoBeta;
    Chrono chronoUpdate;
    Chrono chronoElemComp;
    Chrono chronoAssembly;

    // Stabilization parameters
    Real coeffBeta/*, coeffDiv*/;

    // local velocity
    VectorElemental beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 0; iVol < M_mesh.numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstSecondDeriv( M_mesh.volumeList( iVol ) );

        // stabilization parameters computation
        chronoBeta.start();
        this->computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
        chronoBeta.stop();

        chronoElemComp.start();
        M_elMat.zero();

        // coeffBeta (\beta \nabla \mathbf{u} , \beta \nabla \mathbf{v})
        //
        this->bgradu_bgradv(coeffBeta, beta, M_elMat, M_fe, 0, 0, nDimensions);

        // coeffBeta ( - \mu L \mathbf{u} , \beta \nabla \mathbf{v} )
        //
        //this->lapu_bgradv(-coeffBeta*M_viscosity, beta, M_elMat, M_fe, 0, 0, nDimensions);

        // coeffDiv ( div \mathbf{u} , div \nabla \mathbf{v})
        //
        //stiff_div(coeffDiv, M_elMat, M_fe );

        chronoElemComp.stop();

        chronoAssembly.start();
        for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            for ( UInt jComp = 0; jComp < nDimensions; ++jComp )
                assemb_mat( matrix, M_elMat, M_fe, M_dof, iComp, jComp );
        chronoAssembly.stop();


    }// loop on elements

    debugStream(7100) << std::endl;
    debugStream(7100) << "      Updating of element   done in "
    << chronoUpdate.diffCumul()   << "s." << std::endl;
    debugStream(7100) << "      Determination of parameters done in "
    << chronoBeta.diffCumul()     << "s." << std::endl;
    debugStream(7100) << "      Element computations  done in "
    << chronoElemComp.diffCumul() << "s." << std::endl;
    debugStream(7100) << "      Assembly              done in "
    << chronoAssembly.diffCumul() << "s." << std::endl;


} // applySD(...)

template<typename MeshType, typename DofType>
template <typename VectorType, typename SourceType >
void StabilizationSD<MeshType, DofType>::applyRHS(const Real dt, VectorType& vector, const VectorType& state,
                                                  const SourceType& source, const Real& time)
{
    if ( M_gammaBeta == 0 )
        return;


    Chrono chronoBeta;
    Chrono chronoUpdate;
    Chrono chronoElemComp;
    Chrono chronoAssembly;

    // stabilization parameters
    Real coeffBeta, coeffDiv;

    // local velocity
    VectorElemental beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 0; iVol < M_mesh.numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstDeriv( M_mesh.volumeList( iVol ) );
        // local mesh parameters
        chronoUpdate.stop();

        chronoBeta.start();
        this->computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
        chronoBeta.stop();

        chronoElemComp.start();
        M_elVec.zero();

        // coeffBeta ( f , \beta \nabla \mathbf{v})
        //
        this->f_bgradv(coeffBeta, source, beta, M_elVec, M_fe, 0, time);

        // coeffBeta  ( f , \nabla q )
        //
        this->f_gradq(coeffBeta, source, M_elVec, M_fe, nDimensions, time);
        chronoElemComp.stop();

        chronoAssembly.start();
        for ( UInt iComp = 0; iComp <= nDimensions; ++iComp )
            //assemb_vec( vector, M_elVec, M_fe, M_dof, iComp);
            assembleVector( vector, M_elVec, M_fe, M_dof, iComp, iComp*M_dof.numTotalDof());
        chronoAssembly.stop();


    }// loop on elements

    debugStream(7100) << std::endl;
    debugStream(7100) << "      Updating of element   done in "
    << chronoUpdate.diffCumul()   << "s." << std::endl;
    debugStream(7100) << "      Determination of parameters done in "
    << chronoBeta.diffCumul()     << "s." << std::endl;
    debugStream(7100) << "      Element computations  done in "
    << chronoElemComp.diffCumul() << "s." << std::endl;
    debugStream(7100) << "      Assembly              done in "
    << chronoAssembly.diffCumul() << "s." << std::endl;


} // applyRHS(...)

template<typename MeshType, typename DofType>
void StabilizationSD<MeshType, DofType>::showMe(std::ostream & output) const
{
    output << "StabilizationSD::showMe() " <<std::endl;
    output << "Fluid Viscosity: " << M_viscosity << std::endl;
    output << "Stabilization coefficient SUPG/SD:  " << M_gammaBeta << std::endl;
    output << "Stabilization coefficient div grad: "<< M_gammaDiv   << std::endl;
    M_mesh.showMe(output);
    M_dof.showMe(output);
}

//=============================================================================
// Private Method
//=============================================================================
template<typename MeshType, typename DofType>
template<typename VectorType>
void StabilizationSD<MeshType, DofType>::computeParameters(const Real dt, const UInt iVol, const VectorType& state,
                                                           VectorElemental& beta, Real& coeffBeta, Real& coeffDiv) const
{

    const UInt nDof = M_dof.numTotalDof();

    // square local mesh parameter in 1-norm
    Real hK,hK2,hK4;

    hK = M_fe.diameter();
    hK2 = hK * hK;
    hK4 = hK2 * hK2;

    // determine bmax = ||\beta||_{0,\infty,K}
    // first, get the local velocity into beta
    for ( UInt iNode = 0; iNode < M_fe.nbNode; ++iNode )
    {
        for ( UInt iCoor = 0; iCoor < M_fe.nbCoor(); ++iCoor )
        {
            UInt ig = M_dof.localToGlobalMap( iVol, iNode )+iCoor*nDof;
            beta.vec()[ iCoor*M_fe.nbNode + iNode ] = state[ig];
        }
    }

    // second, calculate its max norm
    Real bmax = std::fabs( beta.vec()[ 0 ] );
    for ( UInt l = 1; l < UInt( M_fe.nbCoor()*M_fe.nbNode ); ++l )
    {
        if ( bmax < std::fabs( beta.vec()[ l ] ) )
            bmax = std::fabs( beta.vec()[ l ] );
    }

    coeffBeta = M_gammaBeta  / std::sqrt( 4./( dt * dt)
                                     + 4.*bmax*bmax/hK2
                                     + 16.*M_viscosity*M_viscosity/hK4 );
    coeffDiv = M_gammaDiv * bmax * hK;

}


template<typename MeshType, typename DofType>
void StabilizationSD<MeshType, DofType>::gradp_bgradv(const Real& coef, VectorElemental& vel,
                                                      MatrixElemental& elmat,const CurrentFE& fe)  const
{
    ASSERT_PRE(fe.hasFirstDeriv(),
               "advection_grad  matrix needs at least the first derivatives");

    MatrixElemental::matrix_type v(fe.nbCoor(),fe.nbQuadPt());
    Real s;


    // local velocity at quadrature points
    for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
    {
        for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
        {
            VectorElemental::vector_view velicoor=vel.block(icoor);
            v(icoor,ig)=0.;
            for (UInt k(0); k<fe.nbNode; ++k)
            {
                v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
            }
        }
    }

    for (UInt ic(0); ic < fe.nbCoor(); ++ic)
    {
        MatrixElemental::matrix_view mat_ic3 = elmat.block(ic,fe.nbCoor());
        MatrixElemental::matrix_view mat_3ic = elmat.block(fe.nbCoor(),ic);
        for (UInt i=0; i<fe.nbNode; ++i)
        {
            for (UInt j=0; j<fe.nbNode; ++j)
            {
                s = 0.0;
                for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
                    for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor)
                        s += fe.phiDer(j,ic,ig)*v(jcoor,ig)*fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
                mat_ic3(i,j) += coef*s;
                mat_3ic(j,i) += coef*s;
            }
        }
    }
}


template<typename MeshType, typename DofType>
void StabilizationSD<MeshType, DofType>::bgradu_bgradv(const Real& coef, VectorElemental& vel, MatrixElemental& elmat, const CurrentFE& fe,
                                                       UInt iblock, UInt jblock, UInt nb)  const
{
    ASSERT_PRE(fe.hasFirstDeriv(),
               "advection (vect) matrix needs at least the first derivatives");


    MatrixElemental::matrix_type mat_tmp(fe.nbNode,fe.nbNode);
    MatrixElemental::matrix_type v( fe.nbCoor(),fe.nbQuadPt() );
    Real s;


    // compute local vectors values
    for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
    {
        for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
        {
            VectorElemental::vector_view velicoor=vel.block(icoor);
            v(icoor,ig)=0.;
            for (UInt k(0); k<fe.nbNode; k++)
                v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
        }
    }

    for (UInt i(0); i<fe.nbNode; ++i)
    {
        for (UInt j(0); j<fe.nbNode; ++j)
        {
            s = 0.0;

            for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
                for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
                    for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor)
                        s += fe.phiDer(i,jcoor,ig)*v(jcoor,ig)*v(icoor,ig)*fe.phiDer(j,icoor,ig)*fe.weightDet(ig);
            mat_tmp(i,j) = coef*s;
        }
    }

    // copy on the components
    for (UInt icomp(0); icomp<nb; icomp++)
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block(iblock+icomp,jblock+icomp);
        for (UInt i(0); i<fe.nbDiag(); ++i)
        {
            for (UInt j(0); j<fe.nbDiag(); ++j)
            {
                mat_icomp(i,j) += mat_tmp(i,j);
            }
        }
    }
}



template<typename MeshType, typename DofType>
void StabilizationSD<MeshType, DofType>::lapu_bgradv(const Real& coef, VectorElemental& vel, MatrixElemental& elmat, const CurrentFE& fe,
                                                     UInt iblock, UInt jblock, UInt nb)  const
{


    ASSERT_PRE(fe.hasFirstDeriv(),
               "lapu_bgradv matrix needs first derivatives");

    ASSERT_PRE(fe.hasSecondDeriv(),
               "lapu_bgradv matrix needs second derivatives");

    MatrixElemental::matrix_type mat_tmp(fe.nbNode,fe.nbNode);
    MatrixElemental::matrix_type v( fe.nbCoor(),fe.nbQuadPt() );
    Real s;


    // compute local vectors values at quadrature points
    for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
    {
        for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
        {
            VectorElemental::vector_view velicoor=vel.block(icoor);
            v(icoor,ig)=0.;
            for (UInt k(0); k<fe.nbNode; ++k)
                v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
        }
    }


    // numerical integration
    for (UInt i(0); i<fe.nbNode; ++i)
    {
        for (UInt j(0); j<fe.nbNode; ++j)
        {
            s = 0.0;
            for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
                for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
                    for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor)
                        s += fe.phiDer2(j,icoor,icoor,ig)*v(jcoor,ig)*fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
            mat_tmp(i,j) = coef*s;
        }
    }

    // copy on the components
    for (UInt icomp(0); icomp<nb; ++icomp)
    {
        MatrixElemental::matrix_view mat_icomp = elmat.block(iblock+icomp,jblock+icomp);
        for (UInt i(0); i<fe.nbDiag(); ++i)
        {
            for (UInt j(0); j<fe.nbDiag(); ++j)
            {
                mat_icomp(i,j) += mat_tmp(i,j);
            }
        }
    }
}


template<typename MeshType, typename DofType>
void StabilizationSD<MeshType, DofType>::lapu_gradq(const Real& coef, MatrixElemental& elmat,const CurrentFE& fe)  const
{

    ASSERT_PRE(fe.hasFirstDeriv(),
               "lapu_gradq matrix needs first derivatives");

    ASSERT_PRE(fe.hasSecondDeriv(),
               "lapu_gradq matrix needs second derivatives");

    Real s;

    for (UInt jc(0); jc < fe.nbCoor(); ++jc) // loop on column blocks
    {
        MatrixElemental::matrix_view mat_view = elmat.block(fe.nbCoor(),jc);
        for (UInt i(0); i<fe.nbNode; ++i) // local rows
        {
            for (UInt j(0); j<fe.nbNode; ++j) // local columns
            {
                s = 0.0;
                // quadrature formula
                for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
                    for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor) // lap
                        s += fe.phiDer2(j,jcoor,jcoor,ig)*fe.phiDer(i,jc,ig)*fe.weightDet(ig);
                mat_view(i,j) += coef*s;
            }
        }
    }
}


template<typename MeshType, typename DofType>
template<typename SourceType>
void StabilizationSD<MeshType, DofType>::f_bgradv(const Real& coef, SourceType& source, VectorElemental& vel,
                                                  VectorElemental& elvec, const CurrentFE& fe, UInt iblock, const Real& time)  const
{

    ASSERT_PRE(fe.hasFirstDeriv(),
               "f_bgradv  vector needs at least the first derivatives");

    MatrixElemental::matrix_type v(fe.nbCoor(),fe.nbQuadPt());
    Real s;


    // local velocity at quadrature points
    for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
    {
        for (UInt icoor(0); icoor<fe.nbCoor(); ++icoor)
        {
            VectorElemental::vector_view velicoor=vel.block(icoor);
            v(icoor,ig)=0.;
            for (UInt k(0); k<fe.nbNode; ++k)
                v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
        }
    }

    // local vector per block
    for (UInt ic(0); ic < fe.nbCoor(); ++ic)
    {
        VectorElemental::vector_view vec_ic = elvec.block(ic+iblock);
        for (UInt i(0); i<fe.nbNode; ++i)
        {
            s = 0.0;
            for (UInt ig(0); ig<fe.nbQuadPt(); ++ig)
                for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor)
                    s += source(time,fe.quadPt(ig,0),fe.quadPt(ig,1),fe.quadPt(ig,2),ic+1)
                         *fe.phiDer(i,jcoor,ig)*v(jcoor,ig)*fe.weightDet(ig);
            vec_ic(i) += coef*s;
        }
    }
}



template<typename MeshType, typename DofType>
template<typename SourceType>
void StabilizationSD<MeshType, DofType>::f_gradq(const Real& coef, SourceType& source, VectorElemental& elvec, const CurrentFE& fe, UInt iblock, const Real& time) const
{

    ASSERT_PRE(fe.hasFirstDeriv(),
               "f_gradq  vector needs at least the first derivatives");

    Real s;

    VectorElemental::vector_view vec_ic = elvec.block(iblock);
    for (UInt i(0); i<fe.nbNode; ++i)
    {
        s = 0.0;
        for (UInt ig=0; ig<fe.nbQuadPt(); ++ig)
            for (UInt jcoor(0); jcoor<fe.nbCoor(); ++jcoor)
                s += source(time,fe.quadPt(ig,0),fe.quadPt(ig,1),fe.quadPt(ig,2),jcoor+1)
                     *fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
        vec_ic(i) += coef*s;
    }
}







} // namespace LifeV

#endif /* _SDSTABILIZATION_HPP_ */
