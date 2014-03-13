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

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/QuadratureRuleProvider.hpp>

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

namespace LifeV
{

//! OseenAssembler - Assembly class for the Oseen problem
/*!
  Signes!!
  Coefficients

    @author Samuel Quinodoz

 */

template< typename meshType, typename matrixType, typename vectorType>
class OseenAssembler
{
public:

    //! @name Public Types
    //@{

    typedef MapEpetra                        map_Type;

    typedef FESpace<meshType, map_Type>      fespace_Type;
    typedef boost::shared_ptr<fespace_Type>  fespacePtr_Type;

    typedef AssemblyElemental::function_Type function_Type;

    typedef boost::shared_ptr<matrixType>    matrixPtr_Type;

    typedef LifeChrono                       chrono_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    OseenAssembler();

    //! Destructor
    virtual ~OseenAssembler() {};

    //@}


    //! @name Methods
    //@{

    //! Setup method for the FESpaces.
    /*!
      This method sets the FESpace for the assembly. With this method, the convective
      field is assumed to be the same as the velocity field. If they differ, use another
      setup method.
     */
    void setup (const fespacePtr_Type& uFESpace, const fespacePtr_Type& pFESpace)
    {
        setup (uFESpace, pFESpace, uFESpace);
    }

    //! Setup method for the FESpace with a different space for the convective field
    void setup (const fespacePtr_Type& uFESpace, const fespacePtr_Type& pFESpace, const fespacePtr_Type& betaFESpace);

    //@}


    //! @name Assembly procedures
    //@{

    //! Add the viscous stress in the standard block
    void addViscousStress (matrixType& matrix, const Real& viscosity)
    {
        addViscousStress (matrix, viscosity, 0, 0);
    };

    //! Add the viscous stress using the given offsets
    void addViscousStress (matrixType& matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the stiff strain in the standard block
    void addStiffStrain (matrixType& matrix, const Real& viscosity)
    {
        addStiffStrain (matrix, viscosity, 0, 0);
    };

    //! Add the stiff strain using the given offsets
    void addStiffStrain (matrixType& matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the term involved in the gradient of the pressure term
    void addGradPressure (matrixType& matrix)
    {
        addGradPressure (matrix, M_uFESpace->dof().numTotalDof() *nDimensions, 0);
    };

    //! Add the term involved in the gradient of the pressure term using the given offsets
    void addGradPressure (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp);

    //! Add the term corresponding to the divergence free constraint
    /*!
     * The default choice coefficient=1.0 leads to a divergence matrix which is the transpose of the
     * pressure gradient matrix.
     */
    void addGradientTranspose (matrixType& matrix, const Real& coefficient = 1.0)
    {
        addGradientTranspose (matrix, 0, M_uFESpace->dof().numTotalDof() *nDimensions, coefficient);
    };

    //! Add the divergence free constraint in a given position of the matrix using the grad calls.
    void addGradientTranspose (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp, const Real& coefficient = 1.0);

    //! Add the term corresponding to the divergence free constraint
    /*!
     * The default choice coefficient=1.0 leads to a divergence matrix which is the transpose of the
     * pressure gradient matrix.
     */
    void addDivergence (matrixType& matrix, const Real& coefficient = 1.0)
    {
        addDivergence (matrix, 0, M_uFESpace->dof().numTotalDof() *nDimensions, coefficient);
    };

    //! Add the divergence free constraint in a given position of the matrix.
    void addDivergence (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp, const Real& coefficient = 1.0);

    //! Add the mass
    void addMass (matrixType& matrix, const Real& coefficient)
    {
        addMass (matrix, coefficient, 0, 0);
    };

    //! Add the mass using offsets
    void addMass (matrixType& matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp);

    //! Add the Pressure mass
    void addPressureMass (matrixType& matrix, const Real& coefficient)
    {
        addPressureMass (matrix, coefficient, M_uFESpace->dof().numTotalDof() *nDimensions,
                         M_uFESpace->dof().numTotalDof() *nDimensions);
    }

    //! Add the mass using offsets
    void addPressureMass (matrixType& matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp);

    //! Add a consistent stabilizing term
    void addMassDivW (matrixType& matrix, const Real& coefficient, const vectorType& beta)
    {
        addMassDivW (matrix, coefficient, beta, 0, 0);
    }

    //! Add a consistent stabilizing term with the given offsets
    void addMassDivW (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp);

    //! Add the convective term
    void addConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta)
    {
        addConvection (matrix, coefficient, beta, 0, 0);
    }

    //! Add the convective term with the given offsets
    void addConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp);

    //! Add the convective term necessary to build the Newton method
    void addNewtonConvection ( matrixType& matrix, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp );

    //! Add the convective term necessary to build the Newton method
    void addNewtonConvection ( matrixType& matrix, const vectorType& beta )
    {
        addNewtonConvection ( matrix, beta, 0, 0 );
    }

    //! Add the convective term
    void addSymmetricConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta)
    {
        addSymmetricConvection (matrix, coefficient, beta, 0, 0);
    }

    //! Add the symmetric convective term with the given offset
    void addSymmetricConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp);

    //! Add an explicit convection term to the right hand side
    void addConvectionRhs (vectorType& rhs, const Real& coefficient, const vectorType& velocity);

    void addMassRhs (vectorType& rhs, const function_Type& fun, const Real& t);

    void addFluxTerms (vectorType& vector, BCHandler const& bcHandler);
    //@}


    //! @name Set Methods
    //@{

    //! Setter for the quadrature used for the right hand side
    /*!
      Beware that calling this function might be quite heavy, so avoid using
      it when it is not necessary.
    */
    inline void setQuadRuleForMassRhs (const QuadratureRule& qr)
    {
        ASSERT (M_massRhsCFE != 0, "No Rhs currentFE for setting the quadrature rule!");
        M_massRhsCFE->setQuadRule (qr);
    }


    //@}


    //! @name Get Methods
    //@{

    //@}

private:

    typedef CurrentFE                            currentFE_Type;
    typedef boost::scoped_ptr<currentFE_Type>    currentFEPtr_Type;

    typedef MatrixElemental                              localMatrix_Type;
    typedef boost::scoped_ptr<localMatrix_Type>          localMatrixPtr_Type;

    typedef VectorElemental                              localVector_Type;
    typedef boost::scoped_ptr<localVector_Type>          localVectorPtr_Type;


    //! @name Private Methods
    //@{

    // No copy constructor
    OseenAssembler (const OseenAssembler&);

    //@}

    // Velocity FE space
    fespacePtr_Type M_uFESpace;

    // Pressure FE space
    fespacePtr_Type M_pFESpace;

    // Beta FE space
    fespacePtr_Type M_betaFESpace;

    // CurrentFE
    currentFEPtr_Type M_viscousCFE;

    currentFEPtr_Type M_gradPressureUCFE;
    currentFEPtr_Type M_gradPressurePCFE;

    currentFEPtr_Type M_divergenceUCFE;
    currentFEPtr_Type M_divergencePCFE;

    currentFEPtr_Type M_massCFE;
    currentFEPtr_Type M_massBetaCFE;
    currentFEPtr_Type M_massPressureCFE;

    currentFEPtr_Type M_convectionUCFE;
    currentFEPtr_Type M_convectionBetaCFE;

    currentFEPtr_Type M_convectionRhsUCFE;

    // CurrentFE for the mass rhs
    currentFEPtr_Type M_massRhsCFE;


    // Local matrix
    localMatrixPtr_Type M_localViscous;

    localMatrixPtr_Type M_localGradPressure;

    localMatrixPtr_Type M_localDivergence;

    localMatrixPtr_Type M_localMass;

    localMatrixPtr_Type M_localMassPressure;

    localMatrixPtr_Type M_localConvection;

    localVectorPtr_Type M_localConvectionRhs;
    // Local vector for the right hand side
    localVectorPtr_Type M_localMassRhs;

};

template< typename meshType, typename matrixType, typename vectorType>
OseenAssembler<meshType, matrixType, vectorType>::
OseenAssembler() :

    M_uFESpace(),
    M_pFESpace(),
    M_betaFESpace(),

    M_viscousCFE(),
    M_gradPressureUCFE(),
    M_gradPressurePCFE(),
    M_divergenceUCFE(),
    M_divergencePCFE(),
    M_massCFE(),
    M_massBetaCFE(),
    M_massPressureCFE(),
    M_convectionUCFE(),
    M_convectionBetaCFE(),
    M_convectionRhsUCFE(),
    M_massRhsCFE(),

    M_localViscous(),
    M_localGradPressure(),
    M_localDivergence(),
    M_localMass(),
    M_localMassPressure(),
    M_localConvection(),
    M_localConvectionRhs(),
    M_localMassRhs()

{}


template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
setup (const fespacePtr_Type& uFESpace, const fespacePtr_Type& pFESpace, const fespacePtr_Type& betaFESpace)
{
    ASSERT (uFESpace != 0, "Impossible to set empty FE space for the velocity. ");
    ASSERT (pFESpace != 0, "Impossible to set empty FE space for the pressure. ");
    ASSERT (betaFESpace != 0, "Impossible to set empty FE space for the convective field (beta). ");

    ASSERT (uFESpace->fieldDim() == nDimensions, "FE space for the velocity has to be vectorial");
    ASSERT (pFESpace->fieldDim() == 1, "FE space for the pressure has to be scalar");
    ASSERT (betaFESpace->fieldDim() == nDimensions, "FE space for the convective field (beta) has to be vectorial");

    M_uFESpace = uFESpace;
    M_pFESpace = pFESpace;
    M_betaFESpace = betaFESpace;

    UInt uDegree (M_uFESpace->polynomialDegree() );
    UInt pDegree (M_pFESpace->polynomialDegree() );
    UInt betaDegree (M_betaFESpace->polynomialDegree() );

    M_viscousCFE.reset (new currentFE_Type (M_uFESpace->refFE(),
                                            M_uFESpace->fe().geoMap(),
                                            QuadratureRuleProvider::provideExactness (TETRA, 2 * uDegree - 2) ) );

    M_gradPressureUCFE.reset (new currentFE_Type (M_uFESpace->refFE(),
                                                  M_uFESpace->fe().geoMap(),
                                                  QuadratureRuleProvider::provideExactness (TETRA, uDegree + pDegree - 1) ) );

    M_gradPressurePCFE.reset (new currentFE_Type (M_pFESpace->refFE(),
                                                  M_uFESpace->fe().geoMap(),
                                                  QuadratureRuleProvider::provideExactness (TETRA, uDegree + pDegree - 1) ) );

    M_divergenceUCFE.reset (new currentFE_Type (M_uFESpace->refFE(),
                                                M_uFESpace->fe().geoMap(),
                                                QuadratureRuleProvider::provideExactness (TETRA, uDegree + pDegree - 1) ) );


    M_divergencePCFE.reset (new currentFE_Type (M_pFESpace->refFE(),
                                                M_uFESpace->fe().geoMap(),
                                                QuadratureRuleProvider::provideExactness (TETRA, uDegree + pDegree - 1) ) );

    M_massCFE.reset (new currentFE_Type (M_uFESpace->refFE(),
                                         M_uFESpace->fe().geoMap(),
                                         QuadratureRuleProvider::provideExactness (TETRA, 2 * uDegree) ) );

    M_massBetaCFE.reset (new currentFE_Type (M_betaFESpace->refFE(),
                                             M_uFESpace->fe().geoMap(),
                                             QuadratureRuleProvider::provideExactness (TETRA, 2 * uDegree) ) );

    M_massPressureCFE.reset (new currentFE_Type (M_pFESpace->refFE(),
                                                 M_pFESpace->fe().geoMap(),
                                                 QuadratureRuleProvider::provideExactness (TETRA, 2 * pDegree) ) );

    M_convectionUCFE.reset (new currentFE_Type (M_uFESpace->refFE(),
                                                M_uFESpace->fe().geoMap(),
                                                QuadratureRuleProvider::provideExactness (TETRA, 2 * uDegree + betaDegree - 1) ) );

    M_convectionBetaCFE.reset (new currentFE_Type (M_betaFESpace->refFE(),
                                                   M_uFESpace->fe().geoMap(),
                                                   QuadratureRuleProvider::provideExactness (TETRA, 2 * uDegree + betaDegree - 1) ) );

    M_convectionRhsUCFE.reset (new currentFE_Type (M_betaFESpace->refFE(),
                                                   M_uFESpace->fe().geoMap(),
                                                   QuadratureRuleProvider::provideExactness (TETRA, 2 * betaDegree + betaDegree - 1) ) );

    M_localViscous.reset (new localMatrix_Type (M_uFESpace->fe().nbFEDof(),
                                                M_uFESpace->fieldDim(),
                                                M_uFESpace->fieldDim() ) );

    M_localGradPressure.reset (new localMatrix_Type (M_uFESpace->fe().nbFEDof(), nDimensions, 0,
                                                     M_pFESpace->fe().nbFEDof(), 0, 1) );

    M_localDivergence.reset (new localMatrix_Type (M_uFESpace->fe().nbFEDof(), 0, nDimensions,
                                                   M_pFESpace->fe().nbFEDof(), 1, 0) );

    M_localMass.reset (new localMatrix_Type (M_uFESpace->fe().nbFEDof(),
                                             M_uFESpace->fieldDim(),
                                             M_uFESpace->fieldDim() ) );

    M_localMassPressure.reset (new localMatrix_Type (M_pFESpace->fe().nbFEDof(),
                                                     M_pFESpace->fieldDim(),
                                                     M_pFESpace->fieldDim() ) );
    M_localConvection.reset (new localMatrix_Type (M_uFESpace->fe().nbFEDof(),
                                                   M_uFESpace->fieldDim(),
                                                   M_uFESpace->fieldDim() ) );

    M_localConvectionRhs.reset (new localVector_Type (M_uFESpace->fe().nbFEDof(), M_uFESpace->fieldDim() ) );

    M_massRhsCFE.reset (new currentFE_Type (M_uFESpace->refFE(), M_uFESpace->fe().geoMap(), M_uFESpace->qr() ) );

    M_localMassRhs.reset (new localVector_Type (M_uFESpace->fe().nbFEDof(), M_uFESpace->fieldDim() ) );



}


template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addViscousStress (matrixType& matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT (M_uFESpace != 0, "No FE space for assembling the viscous stress.");
    ASSERT (offsetLeft + M_uFESpace->dof().numTotalDof() * (M_uFESpace->fieldDim() ) <=
            UInt (matrix.matrixPtr()->NumGlobalCols() ),
            " The matrix is too small (columns) for the assembly of the viscous stress");
    ASSERT (offsetUp + M_uFESpace->dof().numTotalDof() * (M_uFESpace->fieldDim() ) <=
            UInt (matrix.matrixPtr()->NumGlobalRows() ),
            " The matrix is too small (rows) for the assembly of the viscous stress");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_viscousCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localViscous->zero();

        // local stiffness
        AssemblyElemental::stiffness (*M_localViscous, *M_viscousCFE, viscosity, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localViscous,
                             *M_viscousCFE,
                             *M_viscousCFE,
                             M_uFESpace->dof(),
                             M_uFESpace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbTotalDof + offsetUp, iFieldDim * nbTotalDof + offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addStiffStrain (matrixType& matrix, const Real& viscosity, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT (M_uFESpace != 0, "No FE space for assembling the stiff strain.");
    ASSERT (offsetLeft + M_uFESpace->dof().numTotalDof() * (M_uFESpace->fieldDim() ) <=
            UInt (matrix.matrixPtr()->NumGlobalCols() ),
            " The matrix is too small (columns) for the assembly of the stiff strain");
    ASSERT (offsetUp + M_uFESpace->dof().numTotalDof() * (M_uFESpace->fieldDim() ) <=
            UInt (matrix.matrixPtr()->NumGlobalRows() ),
            " The matrix is too small (rows) for the assembly of the stiff strain");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_viscousCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localViscous->zero();

        // local stiffness
        AssemblyElemental::stiffStrain (*M_localViscous, *M_viscousCFE, 2.0 * viscosity, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            for (UInt jFieldDim (0); jFieldDim < fieldDim; ++jFieldDim)
            {
                assembleMatrix ( matrix,
                                 *M_localViscous,
                                 *M_viscousCFE,
                                 *M_viscousCFE,
                                 M_uFESpace->dof(),
                                 M_uFESpace->dof(),
                                 iFieldDim, jFieldDim,
                                 iFieldDim * nbTotalDof + offsetUp, jFieldDim * nbTotalDof + offsetLeft);
            }
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addGradPressure (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp)
{
    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the pressure gradient.");
    ASSERT (M_pFESpace != 0, "No pressure FE space for assembling the pressure gradient.");
    ASSERT (offsetLeft + M_pFESpace->dof().numTotalDof() <=
            UInt (matrix.matrixPtr()->NumGlobalCols() ),
            "The matrix is too small (columns) for the assembly of the pressure gradient");
    ASSERT (offsetUp + M_uFESpace->dof().numTotalDof() *nDimensions <=
            UInt (matrix.matrixPtr()->NumGlobalRows() ),
            " The matrix is too small (rows) for the assembly of the pressure gradient");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_gradPressureUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_gradPressurePCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_WDET );

        // Clean the local matrix
        M_localGradPressure->zero();

        // local stiffness
        AssemblyElemental::grad (*M_localGradPressure, *M_gradPressureUCFE, *M_gradPressurePCFE, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localGradPressure,
                             *M_gradPressureUCFE,
                             *M_gradPressurePCFE,
                             M_uFESpace->dof(),
                             M_pFESpace->dof(),
                             iFieldDim, 0,
                             iFieldDim * nbUTotalDof + offsetUp, offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addGradientTranspose (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp, const Real& coefficient)
{
    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the pressure gradient.");
    ASSERT (M_pFESpace != 0, "No pressure FE space for assembling the pressure gradient.");
    ASSERT (offsetUp + M_pFESpace->dof().numTotalDof() <=
            UInt (matrix.matrixPtr()->NumGlobalCols() ),
            "The matrix is too small (columns) for the assembly of the pressure gradient");
    ASSERT (offsetLeft + M_uFESpace->dof().numTotalDof() *nDimensions <=
            UInt (matrix.matrixPtr()->NumGlobalRows() ),
            " The matrix is too small (rows) for the assembly of the pressure gradient");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_gradPressureUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_gradPressurePCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_WDET );

        // Clean the local matrix
        M_localGradPressure->zero();

        // local stiffness
        AssemblyElemental::grad (*M_localGradPressure, *M_gradPressureUCFE, *M_gradPressurePCFE, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleTransposeMatrix ( matrix,
                                      -coefficient,
                                      *M_localGradPressure,
                                      *M_gradPressurePCFE,
                                      *M_gradPressureUCFE,
                                      M_pFESpace->dof(),
                                      M_uFESpace->dof(),
                                      0, iFieldDim,
                                      offsetUp, iFieldDim * nbUTotalDof + offsetLeft);
        }
    }
}
template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addDivergence (matrixType& matrix, const UInt& offsetLeft, const UInt& offsetUp, const Real& coefficient)
{
    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the divergence.");
    ASSERT (M_pFESpace != 0, "No pressure FE space for assembling the divergence.");
    ASSERT (offsetLeft + M_uFESpace->dof().numTotalDof() *nDimensions <=
            UInt (matrix.matrixPtr()->NumGlobalCols() ),
            "The matrix is too small (columns) for the assembly of the divergence");
    ASSERT (offsetUp + M_pFESpace->dof().numTotalDof() <=
            UInt ( matrix.matrixPtr()->NumGlobalRows() ),
            " The matrix is too small (rows) for the assembly of the divergence");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_divergenceUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_divergencePCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_WDET );

        // Clean the local matrix
        M_localDivergence->zero();

        // local stiffness
        AssemblyElemental::divergence (*M_localDivergence, *M_divergenceUCFE, *M_divergencePCFE, fieldDim, coefficient);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localDivergence,
                             *M_divergencePCFE,
                             *M_divergenceUCFE,
                             M_pFESpace->dof(),
                             M_uFESpace->dof(),
                             0, iFieldDim,
                             offsetUp, iFieldDim * nbUTotalDof + offsetLeft);
        }
    }
}


template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp)
{
    // Beta has to be repeated
    if (beta.mapType() == Unique)
    {
        addConvection (matrix, coefficient, vectorType (beta, Repeated), offsetLeft, offsetUp);
        return;
    }

    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the convection.");
    ASSERT (M_betaFESpace != 0, "No convective FE space for assembling the convection.");
    ASSERT (static_cast<Int> (offsetLeft + M_uFESpace->dof().numTotalDof() *nDimensions) <= matrix.matrixPtr()->NumGlobalCols(),
            "The matrix is too small (columns) for the assembly of the convection");
    ASSERT (static_cast<Int> (offsetUp + M_uFESpace->dof().numTotalDof() *nDimensions) <= matrix.matrixPtr()->NumGlobalRows(),
            " The matrix is too small (rows) for the assembly of the convection");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );
    const UInt nbQuadPt (M_convectionUCFE->nbQuadPt() );

    std::vector< std::vector< Real > > localBetaValue (nbQuadPt, std::vector<Real> (nDimensions, 0.0) );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localConvection->zero();

        // Interpolate
        AssemblyElemental::interpolate (localBetaValue, *M_convectionBetaCFE, nDimensions, M_betaFESpace->dof(), iterElement, beta);

        // Local convection
        AssemblyElemental::advection (*M_localConvection, *M_convectionUCFE, coefficient, localBetaValue, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localConvection,
                             *M_convectionUCFE,
                             *M_convectionUCFE,
                             M_uFESpace->dof(),
                             M_uFESpace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbUTotalDof + offsetUp, iFieldDim * nbUTotalDof + offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addNewtonConvection ( matrixType& matrix, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp )
{
    // Beta has to be repeated
    if ( beta.mapType() == Unique )
    {
        addNewtonConvection ( matrix, vectorType ( beta, Repeated ), offsetLeft, offsetUp );
        return;
    }

    ASSERT ( M_uFESpace != 0, "No velocity FE space for assembling the convection." );
    ASSERT ( M_betaFESpace != 0, "No convective FE space for assembling the convection." );
    ASSERT ( offsetLeft + M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalCols(),
             "The matrix is too small (columns) for the assembly of the convection" );
    ASSERT ( offsetUp + M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalRows(),
             " The matrix is too small (rows) for the assembly of the convection" );

    // Some constants
    const UInt nbElements ( M_uFESpace->mesh()->numElements() );
    const UInt fieldDim ( M_uFESpace->fieldDim() );
    const UInt nbUTotalDof ( M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for ( UInt iterElement ( 0 ); iterElement < nbElements; ++iterElement )
    {
        // Update the diffusion current FE
        M_convectionUCFE->update ( M_uFESpace->mesh()->element ( iterElement ), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localConvection->zero();

        localVector_Type betaLocal ( M_uFESpace->fe().nbFEDof(), M_uFESpace->fieldDim() );

        // Create local vector
        for ( UInt iNode = 0 ; iNode < M_uFESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt iLocal = M_uFESpace->fe().patternFirst ( iNode ); // iLocal = iNode

            for ( Int iComponent = 0; iComponent < fieldDim; ++iComponent )
            {
                UInt iGlobal = M_uFESpace->dof().localToGlobalMap ( iterElement, iLocal ) + iComponent * nbUTotalDof;

                // un local
                betaLocal.vec() [ iLocal + iComponent * M_uFESpace->fe().nbFEDof() ] = beta ( iGlobal );
            }
        }

        // Assembly
        for ( UInt iFieldDim ( 0 ); iFieldDim < nDimensions; ++iFieldDim )
        {
            for ( UInt jFieldDim ( 0 ); jFieldDim < nDimensions; ++jFieldDim )
            {
                AssemblyElemental::advectionNewton ( 1.0, betaLocal, *M_localConvection,
                                                     *M_convectionUCFE, iFieldDim, jFieldDim );
                assembleMatrix ( matrix,
                                 *M_localConvection,
                                 *M_convectionUCFE,
                                 *M_convectionUCFE,
                                 M_uFESpace->dof(),
                                 M_uFESpace->dof(),
                                 iFieldDim, jFieldDim,
                                 iFieldDim * nbUTotalDof + offsetUp, jFieldDim * nbUTotalDof + offsetLeft );
            }
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addSymmetricConvection (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp)

{
    // Beta has to be repeated
    if (beta.mapType() == Unique)
    {
        addSymmetricConvection (matrix, coefficient, vectorType (beta, Repeated), offsetLeft, offsetUp);
        return;
    }

    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the convection.");
    ASSERT (M_betaFESpace != 0, "No convective FE space for assembling the convection.");
    ASSERT ( (int) offsetLeft + (int) M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalCols(),
             "The matrix is too small (columns) for the assembly of the convection");
    ASSERT ( (int) offsetUp + (int) M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalRows(),
             " The matrix is too small (rows) for the assembly of the convection");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );
    const UInt nbQuadPt (M_convectionUCFE->nbQuadPt() );

    std::vector< std::vector< Real > > localBetaValue (nbQuadPt, std::vector<Real> (nDimensions, 0.0) );
    std::vector< std::vector< std::vector< Real > > >
    localBetaGradient (nbQuadPt, std::vector<std::vector<Real> > (nDimensions, std::vector<Real> (nDimensions, 0.0) ) );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );
        M_convectionBetaCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI );

        // Clean the local matrix
        M_localConvection->zero();

        // Interpolate
        AssemblyElemental::interpolate (localBetaValue, *M_convectionBetaCFE, nDimensions, M_betaFESpace->dof(), iterElement, beta);
        // Interpolate
        AssemblyElemental::interpolateGradient (localBetaGradient, *M_convectionBetaCFE, nDimensions, M_betaFESpace->dof(), iterElement, beta);

        // Local convection
        // AssemblyElemental::advection(*M_localConvection,*M_convectionUCFE,localBetaValue,fieldDim);

        // Local convection, 1/2 \beta \grad u v + 1/2 u\grad \beta v
        AssemblyElemental::symmetrizedAdvection (*M_localConvection, *M_convectionUCFE, coefficient, localBetaGradient, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            // Assembly
            for (UInt jFieldDim (0); jFieldDim < nDimensions; ++jFieldDim)
            {
                assembleMatrix ( matrix,
                                 *M_localConvection,
                                 *M_convectionUCFE,
                                 *M_convectionUCFE,
                                 M_uFESpace->dof(),
                                 M_uFESpace->dof(),
                                 iFieldDim, jFieldDim,
                                 iFieldDim * nbUTotalDof + offsetUp, jFieldDim * nbUTotalDof + offsetLeft);
            }
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addConvectionRhs (vectorType& rhs, const Real& coefficient, const vectorType& velocity)
{
    // velocity has to be repeated!
    if (velocity.mapType() == Unique)
    {
        addConvectionRhs (rhs, coefficient, vectorType (velocity, Repeated) );
        return;
    }

    // Check that the FESpace is set
    ASSERT (M_betaFESpace != 0, "No convective FE space for assembling the convection.");

    // Some constants
    const UInt nbElements (M_betaFESpace->mesh()->numElements() );
    const UInt fieldDim (M_betaFESpace->fieldDim() );
    const UInt dim (M_betaFESpace->dim() );
    const UInt nbFEDof (M_convectionRhsUCFE->nbFEDof() );

    // Temporaries
    VectorElemental localVelocity (M_betaFESpace->fe().nbFEDof(), fieldDim);

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionRhsUCFE->update ( M_betaFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local vector
        M_localConvectionRhs->zero();
        localVelocity.zero();

        for ( UInt iNode = 0 ; iNode < nbFEDof ; iNode++ )
        {
            UInt iLocal = M_betaFESpace->fe().patternFirst ( iNode ); // iLocal = iNode

            for ( UInt iComponent = 0; iComponent < fieldDim; ++iComponent )
            {
                UInt iGlobal = M_betaFESpace->dof().localToGlobalMap ( iterElement, iLocal ) + iComponent * dim;

                localVelocity.vec( ) [ iLocal + iComponent * nbFEDof ] = velocity ( iGlobal );
            }
        }

        AssemblyElemental::source_advection (coefficient, localVelocity, localVelocity, *M_localConvectionRhs, *M_convectionRhsUCFE);

        // Here add in the global rhs
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            assembleVector ( rhs,
                             iterElement,
                             *M_localConvectionRhs,
                             nbFEDof,
                             M_betaFESpace->dof(),
                             iterFDim,
                             iterFDim * M_betaFESpace->dof().numTotalDof() );
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addMass (matrixType& matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp)
{

    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the mass.");
    ASSERT (static_cast<Int> (offsetLeft + M_uFESpace->dof().numTotalDof() *nDimensions) <= matrix.matrixPtr()->NumGlobalCols(),
            "The matrix is too small (columns) for the assembly of the mass");
    ASSERT (static_cast<Int> (offsetUp + M_uFESpace->dof().numTotalDof() *nDimensions) <= matrix.matrixPtr()->NumGlobalRows(),
            " The matrix is too small (rows) for the assembly of the mass");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_WDET );

        // Clean the local matrix
        M_localMass->zero();

        // local stiffness
        AssemblyElemental::mass (*M_localMass, *M_massCFE, coefficient, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localMass,
                             *M_massCFE,
                             *M_massCFE,
                             M_uFESpace->dof(),
                             M_uFESpace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbUTotalDof + offsetUp, iFieldDim * nbUTotalDof + offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addPressureMass (matrixType& matrix, const Real& coefficient, const UInt& offsetLeft, const UInt offsetUp)
{

    ASSERT (M_pFESpace != 0, "No pressure FE space for assembling the mass.");
    ASSERT ( (int) offsetLeft + (int) M_pFESpace->dof().numTotalDof() <= matrix.matrixPtr()->NumGlobalCols(),
             "The matrix is too small (columns) for the assembly of the mass");
    ASSERT ( (int) offsetUp + (int) M_pFESpace->dof().numTotalDof() <= matrix.matrixPtr()->NumGlobalRows(),
             " The matrix is too small (rows) for the assembly of the mass");

    // Some constants
    const UInt nbElements (M_pFESpace->mesh()->numElements() );
    const UInt fieldDim (M_pFESpace->fieldDim() );


    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massPressureCFE->update ( M_pFESpace->mesh()->element (iterElement), UPDATE_WDET );

        // Clean the local matrix
        M_localMassPressure->zero();

        // local stiffness
        AssemblyElemental::mass (*M_localMassPressure, *M_massPressureCFE, coefficient, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < fieldDim; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localMassPressure,
                             *M_massPressureCFE,
                             *M_massPressureCFE,
                             M_pFESpace->dof(),
                             M_pFESpace->dof(),
                             iFieldDim, iFieldDim,
                             offsetUp, offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addMassDivW (matrixType& matrix, const Real& coefficient, const vectorType& beta, const UInt& offsetLeft, const UInt offsetUp)
{
    // Beta has to be repeated
    if (beta.mapType() == Unique)
    {
        addMassDivW (matrix, coefficient, vectorType (beta, Repeated), offsetLeft, offsetUp);
        return;
    }

    ASSERT (M_uFESpace != 0, "No velocity FE space for assembling the mass.");
    ASSERT (M_betaFESpace != 0, "No convective FE space for assembling the divergence.");
    ASSERT ( (int) offsetLeft + (int) M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalCols(),
             "The matrix is too small (columns) for the assembly of the mass");
    ASSERT ( (int) offsetUp + (int) M_uFESpace->dof().numTotalDof() *nDimensions <= matrix.matrixPtr()->NumGlobalRows(),
             " The matrix is too small (rows) for the assembly of the mass");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbUTotalDof (M_uFESpace->dof().numTotalDof() );
    const UInt nbQuadPt (M_convectionUCFE->nbQuadPt() );

    std::vector< Real > localBetaDivergence (nbQuadPt);

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_convectionUCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_WDET );
        M_convectionBetaCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_DPHI | UPDATE_WDET );

        // Clean the local matrix
        M_localMass->zero();

        // Interpolate
        AssemblyElemental::interpolateDivergence (localBetaDivergence, *M_convectionBetaCFE, M_betaFESpace->dof(), iterElement, beta);

        // local mass, with coefficients
        AssemblyElemental::massDivW (*M_localMass, *M_convectionUCFE, coefficient, localBetaDivergence, fieldDim);

        // Assembly
        for (UInt iFieldDim (0); iFieldDim < nDimensions; ++iFieldDim)
        {
            assembleMatrix ( matrix,
                             *M_localMass,
                             *M_convectionUCFE,
                             *M_convectionUCFE,
                             M_uFESpace->dof(),
                             M_uFESpace->dof(),
                             iFieldDim, iFieldDim,
                             iFieldDim * nbUTotalDof + offsetUp, iFieldDim * nbUTotalDof + offsetLeft);
        }
    }
}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addMassRhs (vectorType& rhs, const function_Type& fun, const Real& t)
{
    // Check that the fespace is set
    ASSERT (M_uFESpace != 0, "No FE space for assembling the right hand side (mass)!");

    // Some constants
    const UInt nbElements (M_uFESpace->mesh()->numElements() );
    const UInt fieldDim (M_uFESpace->fieldDim() );
    const UInt nbFEDof (M_massRhsCFE->nbFEDof() );

    // Loop over the elements
    for (UInt iterElement (0); iterElement < nbElements; ++iterElement)
    {
        // Update the diffusion current FE
        M_massRhsCFE->update ( M_uFESpace->mesh()->element (iterElement), UPDATE_QUAD_NODES | UPDATE_WDET );

        // Clean the local matrix
        M_localMassRhs->zero();

        AssemblyElemental::bodyForces ( *M_localMassRhs,
                                        *M_massRhsCFE,
                                        fun, t,
                                        fieldDim);

        // Here add in the global rhs
        for (UInt iterFDim (0); iterFDim < fieldDim; ++iterFDim)
        {
            assembleVector ( rhs,
                             iterElement,
                             *M_localMassRhs,
                             nbFEDof,
                             M_uFESpace->dof(),
                             iterFDim,
                             iterFDim * M_uFESpace->dof().numTotalDof() );
        }
    }

}

template< typename meshType, typename matrixType, typename vectorType>
void
OseenAssembler<meshType, matrixType, vectorType>::
addFluxTerms ( vectorType&     vector,
               BCHandler const& bcHandler)
{

    for ( ID hCounter = 0; hCounter < bcHandler.size(); ++hCounter )
    {
        ASSERT ( bcHandler[ hCounter ].type()  == Flux, "Works only with Flux BC type!");
        ASSERT ( bcHandler.bcUpdateDone () , " Please call bcHandler::Update() before calling this method!");

        const BCBase&    boundaryCond (bcHandler[ hCounter ]);

        // Number of local DOF in this facet
        UInt nDofF = M_uFESpace->feBd().nbFEDof();

        // Number of total scalar Dof
        UInt totalDof = M_uFESpace->dof().numTotalDof();

        // Number of components involved in this boundary condition
        UInt nComp = boundaryCond.numberOfComponents();

        Real sum;

        const BCIdentifierNatural* pId;
        ID ibF, idDof;

        if ( !boundaryCond.isDataAVector() )
        {
            for ( ID i = 0; i < boundaryCond.list_size(); ++i )
            {
                pId = static_cast< const BCIdentifierNatural* > ( boundaryCond[ i ] );

                // Number of the current boundary facet
                ibF = pId->id();
                // Updating facet stuff
                M_uFESpace->feBd().update ( M_uFESpace->mesh()->boundaryFacet ( ibF ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );

                for ( ID idofF = 0; idofF < nDofF; ++idofF )
                {
                    for ( int ic = 0; ic < (int) nComp; ++ic)
                    {
                        idDof = pId->boundaryLocalToGlobalMap ( idofF ) + ic * totalDof;

                        sum = 0.;
                        for ( int iq = 0; iq < (int) M_uFESpace->feBd().nbQuadPt(); ++iq )
                        {
                            sum += M_uFESpace->feBd().phi ( int ( idofF ), iq ) *
                                   M_uFESpace->feBd().normal (ic , iq) *
                                   M_uFESpace->feBd().wRootDetMetric (iq);
                        }

                        vector.sumIntoGlobalValues (idDof, sum);
                    }
                }
            }
        }

    }
    vector.globalAssemble();
}





} // Namespace LifeV

#endif /* OSEENASSEMBLER_H */
