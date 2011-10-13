//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 16 Nov 2010

    A more detailed description of the file (if necessary)
 */

#ifndef OSEENASSEMBLER_H
#define OSEENASSEMBLER_H 1

#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/LifeV.hpp>

#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/QuadratureRuleProvider.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

namespace LifeV {

//! OseenAssembler - Assembly class for the Oseen problem
/*!
  Signes!!
  Coefficients

    @author Samuel Quinodoz

 */

template< typename mesh_type, typename matrix_type, typename vector_type>
class OseenAssembler
{
public:

    //! @name Public Types
    //@{

    typedef MapEpetra                       map_type;

    typedef FESpace<mesh_type, map_type>    fespace_type;
    typedef boost::shared_ptr<fespace_type> fespace_ptrType;

    typedef boost::shared_ptr<matrix_type>  matrix_ptrType;

    typedef LifeChrono                      chrono_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    OseenAssembler();

    //! Destructor
    virtual ~OseenAssembler(){};

    //@}


    //! @name Methods
    //@{

    //! Setup method for the FESpaces.
    /*!
      This method sets the FESpace for the assembly. With this method, the convective
      field is assumed to be the same as the velocity field. If they differ, use another
      setup method.
     */
    void setup(const fespace_ptrType& uFESpace, const fespace_ptrType& pFESpace)
    {
        setup(uFESpace,pFESpace,uFESpace);
    }

    //! Setup method for the FESpace with a different space for the convective field
    void setup(const fespace_ptrType& uFESpace, const fespace_ptrType& pFESpace, const fespace_ptrType& betaFESpace);

    //@}


    //! @name Assembly procedures
    //@{

    //! Add the viscous stress in the standard block
    void addViscousStress(matrix_ptrType matrix, const Real& viscosity)
    {
        addViscousStress(matrix,viscosity,0,0);
    };

    //! Add the viscous stress using the given offsets
    void addViscousStress(matrix_ptrType matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the stiff strain in the standard block
    void addStiffStrain(matrix_ptrType matrix, const Real& viscosity)
    {
        addStiffStrain(matrix,viscosity,0,0);
    };

    //! Add the stiff strain using the given offsets
    void addStiffStrain(matrix_ptrType matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the term involved in the gradient of the pressure term
    void addGradPressure(matrix_ptrType matrix)
    {
        addGradPressure(matrix,M_uFESpace->dof().numTotalDof()*nDimensions,0);
    };

    //! Add the term involved in the gradient of the pressure term using the given offsets
    void addGradPressure(matrix_ptrType matrix, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the term corresponding to the divergence free constraint
    /*!
     * The default choice coefficient=1.0 leads to a divergence matrix which is the transpose of the
     * pressure gradient matrix.
     */
    void addDivergence(matrix_ptrType matrix,const Real& coefficient=1.0)
    {
        addDivergence(matrix,0,M_uFESpace->dof().numTotalDof()*nDimensions,coefficient);
    };

    //! Add the divergence free constraint in a given position of the matrix.
    void addDivergence(matrix_ptrType matrix, const UInt& offsetLeft, const UInt& offsetUp,const Real& coefficient=1.0);

    //! Add the mass
    void addMass(matrix_ptrType matrix, const Real& coefficient)
    {
        addMass(matrix,coefficient,0,0);
    };

    //! Add the mass using offsets
    void addMass(matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp);

    //! Add the convective term
    void addConvection(matrix_ptrType matrix, const vector_type& beta)
    {
        addConvection(matrix,beta,0,0);
    }

    //! Add the convective term with the given offsets
    void addConvection(matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt offsetUp);

    //! Add the convective term necessary to build the Newton method
    void addNewtonConvection( matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt offsetUp );

    //! Add the convective term necessary to build the Newton method
    void addNewtonConvection( matrix_ptrType matrix, const vector_type& beta )
    {
        addNewtonConvection( matrix, beta, 0, 0 );
    }

    //! Add an explicit convection term to the right hand side
    void addConvectionRhs(vector_type& rhs, const vector_type& velocity);

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //@}

private:

    typedef CurrentFE                                    currentFE_type;
    typedef boost::scoped_ptr<currentFE_type>            currentFE_ptrType;

    typedef MatrixElemental                              localMatrix_type;
    typedef boost::scoped_ptr<localMatrix_type>          localMatrix_ptrType;

    typedef VectorElemental                              localVector_type;
    typedef boost::scoped_ptr<localVector_type>          localVector_ptrType;


    //! @name Private Methods
    //@{

    // No copy constructor
    OseenAssembler(const OseenAssembler&);

    //@}

    // Velocity FE space
    fespace_ptrType M_uFESpace;

    // Pressure FE space
    fespace_ptrType M_pFESpace;

    // Beta FE space
    fespace_ptrType M_betaFESpace;


    // CurrentFE
    currentFE_ptrType M_viscousCFE;

    currentFE_ptrType M_gradPressureUCFE;
    currentFE_ptrType M_gradPressurePCFE;

    currentFE_ptrType M_divergenceUCFE;
    currentFE_ptrType M_divergencePCFE;

    currentFE_ptrType M_massCFE;

    currentFE_ptrType M_convectionUCFE;
    currentFE_ptrType M_convectionBetaCFE;

    currentFE_ptrType M_convectionRhsUCFE;


    // Local matrix
    localMatrix_ptrType M_localViscous;

    localMatrix_ptrType M_localGradPressure;

    localMatrix_ptrType M_localDivergence;

    localMatrix_ptrType M_localMass;

    localMatrix_ptrType M_localConvection;

    localVector_ptrType M_localConvectionRhs;

};

template< typename mesh_type, typename matrix_type, typename vector_type>
OseenAssembler<mesh_type,matrix_type,vector_type>::
OseenAssembler():

    M_uFESpace(),
    M_pFESpace(),
    M_betaFESpace(),

    M_viscousCFE(),
    M_gradPressureUCFE(),
    M_gradPressurePCFE(),
    M_divergenceUCFE(),
    M_divergencePCFE(),
    M_massCFE(),
    M_convectionUCFE(),
    M_convectionBetaCFE(),
    M_convectionRhsUCFE(),

    M_localViscous(),
    M_localGradPressure(),
    M_localDivergence(),
    M_localMass(),
    M_localConvection(),
    M_localConvectionRhs()
{}


template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
setup(const fespace_ptrType& uFESpace, const fespace_ptrType& pFESpace, const fespace_ptrType& betaFESpace)
{
    ASSERT(uFESpace !=0, "Impossible to set empty FE space for the velocity. ");
    ASSERT(pFESpace !=0, "Impossible to set empty FE space for the pressure. ");
    ASSERT(betaFESpace !=0, "Impossible to set empty FE space for the convective field (beta). ");

    ASSERT(uFESpace->fieldDim() == nDimensions, "FE space for the velocity has to be vectorial");
    ASSERT(pFESpace->fieldDim() == 1, "FE space for the pressure has to be scalar");
    ASSERT(betaFESpace->fieldDim() == nDimensions, "FE space for the convective field (beta) has to be vectorial");

    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_betaFESpace = betaFESpace;

    UInt uDegree(M_uFESpace->polynomialDegree());
    UInt pDegree(M_pFESpace->polynomialDegree());
    UInt betaDegree(M_betaFESpace->polynomialDegree());

    M_viscousCFE.reset(new currentFE_type(M_uFESpace->refFE(),
                                          M_uFESpace->fe().geoMap(),
                                          QuadratureRuleProvider::provideExactnessMax(TETRA,2*uDegree-2)));

    M_gradPressureUCFE.reset(new currentFE_type(M_uFESpace->refFE(),
                                                M_uFESpace->fe().geoMap(),
                                                QuadratureRuleProvider::provideExactnessMax(TETRA,uDegree+pDegree-1)));

    M_gradPressurePCFE.reset(new currentFE_type(M_pFESpace->refFE(),
                                                M_uFESpace->fe().geoMap(),
                                                QuadratureRuleProvider::provideExactnessMax(TETRA,uDegree+pDegree-1)));

    M_divergenceUCFE.reset(new currentFE_type(M_uFESpace->refFE(),
                                              M_uFESpace->fe().geoMap(),
                                              QuadratureRuleProvider::provideExactnessMax(TETRA,uDegree+pDegree-1)));


    M_divergencePCFE.reset(new currentFE_type(M_pFESpace->refFE(),
                                              M_uFESpace->fe().geoMap(),
                                              QuadratureRuleProvider::provideExactnessMax(TETRA,uDegree+pDegree-1)));

    M_massCFE.reset(new currentFE_type(M_uFESpace->refFE(),
                                       M_uFESpace->fe().geoMap(),
                                       QuadratureRuleProvider::provideExactnessMax(TETRA,2*uDegree)));

    M_convectionUCFE.reset(new currentFE_type(M_uFESpace->refFE(),
                                              M_uFESpace->fe().geoMap(),
                                              QuadratureRuleProvider::provideExactnessMax(TETRA,2*uDegree+betaDegree-1)));

    M_convectionBetaCFE.reset(new currentFE_type(M_betaFESpace->refFE(),
                                                 M_uFESpace->fe().geoMap(),
                                                 QuadratureRuleProvider::provideExactnessMax(TETRA,2*uDegree+betaDegree-1)));

    M_convectionRhsUCFE.reset(new currentFE_type(M_betaFESpace->refFE(),
                                                 M_uFESpace->fe().geoMap(),
                                                 QuadratureRuleProvider::provideExactnessMax(TETRA,2*betaDegree+betaDegree-1)));

    M_localViscous.reset(new localMatrix_type(M_uFESpace->fe().nbFEDof(),
                                              M_uFESpace->fieldDim(),
                                              M_uFESpace->fieldDim()));
    M_localGradPressure.reset(new localMatrix_type(M_uFESpace->fe().nbFEDof(),nDimensions,0,
                                                   M_pFESpace->fe().nbFEDof(),0,1));
    M_localDivergence.reset(new localMatrix_type(M_uFESpace->fe().nbFEDof(),0,nDimensions,
                                                 M_pFESpace->fe().nbFEDof(),1,0));
    M_localMass.reset(new localMatrix_type(M_uFESpace->fe().nbFEDof(),
                                           M_uFESpace->fieldDim(),
                                           M_uFESpace->fieldDim()));
    M_localConvection.reset(new localMatrix_type(M_uFESpace->fe().nbFEDof(),
                                                 M_uFESpace->fieldDim(),
                                                 M_uFESpace->fieldDim()));

    M_localConvectionRhs.reset(new localVector_type(M_uFESpace->fe().nbFEDof(), M_uFESpace->fieldDim()));

}


template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addViscousStress(matrix_ptrType matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT(M_uFESpace != 0, "No FE space for assembling the viscous stress.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the viscous stress with no matrix.");
    ASSERT(offsetLeft + M_uFESpace->dof().numTotalDof()*(M_uFESpace->fieldDim()) <=
           UInt(matrix->matrixPtr()->NumGlobalCols()),
           " The matrix is too small (columns) for the assembly of the viscous stress");
    ASSERT(offsetUp + M_uFESpace->dof().numTotalDof()*(M_uFESpace->fieldDim()) <=
           UInt(matrix->matrixPtr()->NumGlobalRows()),
           " The matrix is too small (rows) for the assembly of the viscous stress");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbTotalDof(M_uFESpace->dof().numTotalDof());

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_viscousCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localViscous->zero();

        // local stiffness
        AssemblyElemental::stiffness(*M_localViscous,*M_viscousCFE,viscosity,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localViscous,
                            *M_viscousCFE,
                            *M_viscousCFE,
                            M_uFESpace->dof(),
                            M_uFESpace->dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*nbTotalDof + offsetUp, iFieldDim*nbTotalDof + offsetLeft);
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addStiffStrain(matrix_ptrType matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT(M_uFESpace != 0, "No FE space for assembling the stiff strain.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the stiff strain with no matrix.");
    ASSERT(offsetLeft + M_uFESpace->dof().numTotalDof()*(M_uFESpace->fieldDim()) <=
           UInt(matrix->matrixPtr()->NumGlobalCols()),
           " The matrix is too small (columns) for the assembly of the stiff strain");
    ASSERT(offsetUp + M_uFESpace->dof().numTotalDof()*(M_uFESpace->fieldDim()) <=
           UInt(matrix->matrixPtr()->NumGlobalRows()),
           " The matrix is too small (rows) for the assembly of the stiff strain");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbTotalDof(M_uFESpace->dof().numTotalDof());

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_viscousCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localViscous->zero();

        // local stiffness
        AssemblyElemental::stiffStrain(*M_localViscous,*M_viscousCFE,2.0*viscosity,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            for (UInt jFieldDim(0); jFieldDim<fieldDim; ++jFieldDim)
            {
                assembleMatrix( *matrix,
                                *M_localViscous,
                                *M_viscousCFE,
                                *M_viscousCFE,
                                M_uFESpace->dof(),
                                M_uFESpace->dof(),
                                iFieldDim, jFieldDim,
                                iFieldDim*nbTotalDof + offsetUp, jFieldDim*nbTotalDof + offsetLeft);
            }
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addGradPressure(matrix_ptrType matrix, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT(M_uFESpace != 0, "No velocity FE space for assembling the pressure gradient.");
    ASSERT(M_pFESpace != 0, "No pressure FE space for assembling the pressure gradient.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the pressure gradient with no matrix.");
    ASSERT(offsetLeft + M_pFESpace->dof().numTotalDof() <=
           UInt(matrix->matrixPtr()->NumGlobalCols()),
           "The matrix is too small (columns) for the assembly of the pressure gradient");
    ASSERT(offsetUp + M_uFESpace->dof().numTotalDof()*nDimensions <=
           UInt(matrix->matrixPtr()->NumGlobalRows()),
           " The matrix is too small (rows) for the assembly of the pressure gradient");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbUTotalDof(M_uFESpace->dof().numTotalDof());

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_gradPressureUCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_gradPressurePCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localGradPressure->zero();

        // local stiffness
        AssemblyElemental::grad(*M_localGradPressure,*M_gradPressureUCFE,*M_gradPressurePCFE,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<nDimensions; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localGradPressure,
                            *M_gradPressureUCFE,
                            *M_gradPressurePCFE,
                            M_uFESpace->dof(),
                            M_pFESpace->dof(),
                            iFieldDim, 0,
                            iFieldDim*nbUTotalDof + offsetUp, offsetLeft);
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addDivergence(matrix_ptrType matrix, const UInt& offsetLeft, const UInt& offsetUp, const Real& coefficient)
{
    ASSERT(M_uFESpace != 0, "No velocity FE space for assembling the divergence.");
    ASSERT(M_pFESpace != 0, "No pressure FE space for assembling the divergence.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the divergence with no matrix.");
    ASSERT(offsetLeft + M_uFESpace->dof().numTotalDof()*nDimensions <=
           UInt(matrix->matrixPtr()->NumGlobalCols()),
           "The matrix is too small (columns) for the assembly of the divergence");
    ASSERT(offsetUp + M_pFESpace->dof().numTotalDof() <=
           UInt( matrix->matrixPtr()->NumGlobalRows()),
           " The matrix is too small (rows) for the assembly of the divergence");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbUTotalDof(M_uFESpace->dof().numTotalDof());

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_divergenceUCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_divergencePCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localDivergence->zero();

        // local stiffness
        AssemblyElemental::divergence(*M_localDivergence,*M_divergenceUCFE,*M_divergencePCFE,fieldDim,coefficient);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<nDimensions; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localDivergence,
                            *M_divergencePCFE,
                            *M_divergenceUCFE,
                            M_pFESpace->dof(),
                            M_uFESpace->dof(),
                            0, iFieldDim,
                            offsetUp, iFieldDim*nbUTotalDof + offsetLeft);
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addConvection(matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt offsetUp)
{
    // Beta has to be repeated
    if (beta.mapType() == Unique)
    {
        addConvection(matrix,vector_type(beta,Repeated),offsetLeft,offsetUp);
        return;
    }

    ASSERT(M_uFESpace != 0, "No velocity FE space for assembling the convection.");
    ASSERT(M_betaFESpace != 0, "No convective FE space for assembling the convection.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the convection with no matrix.");
    ASSERT(offsetLeft + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalCols(),
           "The matrix is too small (columns) for the assembly of the convection");
    ASSERT(offsetUp + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalRows(),
           " The matrix is too small (rows) for the assembly of the convection");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbUTotalDof(M_uFESpace->dof().numTotalDof());
    const UInt nbQuadPt(M_convectionUCFE->nbQuadPt());

    std::vector< std::vector< Real > > localBetaValue(nbQuadPt,std::vector<Real>(nDimensions,0.0));

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionUCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_convectionBetaCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_PHI );

        // Clean the local matrix
        M_localConvection->zero();

        // Interpolate
        AssemblyElemental::interpolate(localBetaValue,*M_convectionBetaCFE,nDimensions,M_betaFESpace->dof(),iterElement,beta);

        // Local convection
        AssemblyElemental::advection(*M_localConvection,*M_convectionUCFE,localBetaValue,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<nDimensions; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localConvection,
                            *M_convectionUCFE,
                            *M_convectionUCFE,
                            M_uFESpace->dof(),
                            M_uFESpace->dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*nbUTotalDof + offsetUp, iFieldDim*nbUTotalDof + offsetLeft);
                            }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addNewtonConvection( matrix_ptrType matrix, const vector_type& beta, const UInt& offsetLeft, const UInt offsetUp )
{
    // Beta has to be repeated
    if ( beta.mapType() == Unique )
    {
        addNewtonConvection( matrix,vector_type( beta, Repeated ), offsetLeft, offsetUp );
        return;
    }

    ASSERT( M_uFESpace != 0, "No velocity FE space for assembling the convection." );
    ASSERT( M_betaFESpace != 0, "No convective FE space for assembling the convection." );
    ASSERT( matrix !=0, "Cannot perform the assembly of the convection with no matrix." );
    ASSERT( offsetLeft + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalCols(),
            "The matrix is too small (columns) for the assembly of the convection" );
    ASSERT( offsetUp + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalRows(),
            " The matrix is too small (rows) for the assembly of the convection" );

    // Some constants
    const UInt nbElements( M_uFESpace->mesh()->numElements() );
    const UInt fieldDim( M_uFESpace->fieldDim() );
    const UInt nbUTotalDof( M_uFESpace->dof().numTotalDof() );
    //const UInt nbQuadPt( M_convectionUCFE->nbQuadPt() ); //unused parameter

    // Loop over the elements
    for ( UInt iterElement( 0 ); iterElement < nbElements; ++iterElement )
    {
        // Update the diffusion current FE
        M_convectionUCFE->update( M_uFESpace->mesh()->element( iterElement ), UPDATE_DPHI | UPDATE_WDET );
        M_convectionBetaCFE->update( M_uFESpace->mesh()->element( iterElement ), UPDATE_PHI );

        // Clean the local matrix
        M_localConvection->zero();

        localVector_type betaLocal( M_uFESpace->fe().nbFEDof(), M_uFESpace->fieldDim() );

        // Create local vector
        for ( UInt iNode = 0 ; iNode < M_uFESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt iLocal = M_uFESpace->fe().patternFirst( iNode ); // iLocal = iNode

            for ( Int iComponent = 0; iComponent < fieldDim; ++iComponent )
            {
                UInt iGlobal = M_uFESpace->dof().localToGlobalMap( iterElement, iLocal ) + iComponent * nbUTotalDof;

                // un local
                betaLocal.vec()[ iLocal + iComponent*M_uFESpace->fe().nbFEDof() ] = beta( iGlobal );
            }
        }

        // Assembly
        for ( UInt iFieldDim( 0 ); iFieldDim<nDimensions; ++iFieldDim )
        {
            for ( UInt jFieldDim( 0 ); jFieldDim < nDimensions; ++jFieldDim )
            {
                AssemblyElemental::advectionNewton( 1.0, betaLocal, *M_localConvection,
                                                    *M_convectionUCFE, iFieldDim, jFieldDim );
                assembleMatrix( *matrix,
                                *M_localConvection,
                                *M_convectionUCFE,
                                *M_convectionUCFE,
                                M_uFESpace->dof(),
                                M_uFESpace->dof(),
                                iFieldDim, jFieldDim,
                                iFieldDim*nbUTotalDof + offsetUp, jFieldDim*nbUTotalDof + offsetLeft );
            }
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addConvectionRhs(vector_type& rhs, const vector_type& velocity)
{
    // velocity has to be repeated!
    if (velocity.mapType() == Unique)
    {
        addConvectionRhs(rhs, vector_type(velocity,Repeated));
        return;
    }

    // Check that the FESpace is set
    ASSERT(M_betaFESpace != 0, "No convective FE space for assembling the convection.");

    // Some constants
    const UInt nbElements(M_betaFESpace->mesh()->numElements());
    const UInt fieldDim(M_betaFESpace->fieldDim());
    const UInt dim(M_betaFESpace->dim());
    const UInt nbFEDof(M_convectionRhsUCFE->nbFEDof());

    // Temporaries
    VectorElemental localVelocity(M_betaFESpace->fe().nbFEDof(), fieldDim);

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionRhsUCFE->update( M_betaFESpace->mesh()->element(iterElement), UPDATE_PHI | UPDATE_DPHI | UPDATE_WDET );

        // Clean the local vector
        M_localConvectionRhs->zero();
        localVelocity.zero();

        for ( UInt iNode = 0 ; iNode < nbFEDof ; iNode++ )
        {
            UInt iLocal = M_betaFESpace->fe().patternFirst( iNode ); // iLocal = iNode

            for ( UInt iComponent = 0; iComponent < fieldDim; ++iComponent )
            {
                UInt iGlobal = M_betaFESpace->dof().localToGlobalMap( iterElement, iLocal ) + iComponent * dim;

                localVelocity.vec( ) [ iLocal + iComponent*nbFEDof ] = velocity( iGlobal );
            }
        }

        source_advection(localVelocity,localVelocity,*M_localConvectionRhs, *M_convectionRhsUCFE);

        // Here add in the global rhs
        for (UInt iterFDim(0); iterFDim<fieldDim; ++iterFDim)
        {
            assembleVector( rhs,
                            iterElement,
                            *M_localConvectionRhs,
                            nbFEDof,
                            M_betaFESpace->dof(),
                            iterFDim,
                            iterFDim*M_betaFESpace->dof().numTotalDof());
        }
    }
}

template< typename mesh_type, typename matrix_type, typename vector_type>
void
OseenAssembler<mesh_type,matrix_type,vector_type>::
addMass(matrix_ptrType matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp)
{

    ASSERT(M_uFESpace != 0, "No velocity FE space for assembling the mass.");
    ASSERT(matrix !=0, "Cannot perform the assembly of the mass with no matrix.");
    ASSERT(offsetLeft + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalCols(),
           "The matrix is too small (columns) for the assembly of the mass");
    ASSERT(offsetUp + M_uFESpace->dof().numTotalDof()*nDimensions <= matrix->matrixPtr()->NumGlobalRows(),
           " The matrix is too small (rows) for the assembly of the mass");

    // Some constants
    const UInt nbElements(M_uFESpace->mesh()->numElements());
    const UInt fieldDim(M_uFESpace->fieldDim());
    const UInt nbUTotalDof(M_uFESpace->dof().numTotalDof());

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massCFE->update( M_uFESpace->mesh()->element(iterElement), UPDATE_PHI | UPDATE_WDET );

        // Clean the local matrix
        M_localMass->zero();

        // local stiffness
        AssemblyElemental::mass(*M_localMass,*M_massCFE,coefficient,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<nDimensions; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localMass,
                            *M_massCFE,
                            *M_massCFE,
                            M_uFESpace->dof(),
                            M_uFESpace->dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*nbUTotalDof + offsetUp, iFieldDim*nbUTotalDof + offsetLeft);
                            }
    }
}




} // Namespace LifeV

#endif /* OSEENASSEMBLER_H */
