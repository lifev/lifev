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
    @brief File containing the procedures for the local assembly of the differential operators

    @author Gianmarco Mengaldo <gianmarco.mengaldo@gmail.com>
    @author Paolo Tricerri <gianmarco.mengaldo@gmail.com>
    @mantainer Paolo Tricerri <paolo.tricerri@epfl.ch>

    All the methods are described in the report StructuralSolver framework in LifeV: Description and Usage
 */


#ifndef _ELEMOPERSTRUCTURE_H_INCLUDED
#define _ELEMOPERSTRUCTURE_H_INCLUDED

#include <string>
#include <iostream>

//Trilinos includ
#include <Epetra_LAPACK.h>
#include <Epetra_BLAS.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/CurrentFEManifold.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/DOF.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

namespace LifeV
{
//! @name Public typedefs
//@{
typedef boost::numeric::ublas::matrix<Real> Matrix;
typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::zero_matrix<Real> ZeroMatrix;
//@}

/*! /namespace AssemblyElementalStructure

  This namespace is specially designed to contain the elementary
  operations (corresponding to differential operators) that build
  the local contributions to be used in the assembly procedures.

*/
namespace AssemblyElementalStructure
{

//! Save displacement according to a functor
/*!
  @param FE  the finite element space
  @param fct the functor for the saving
  @param originVector the vector from which the values are taken
  @param saveVector the vector from which the values are saved
*/
template<typename FunctorType, typename MeshType, typename MapType>
void saveVectorAccordingToFunctor ( const std::shared_ptr<FESpace<MeshType, MapType> > dispFESpace,
                                    const FunctorType functor,
                                    const VectorEpetra& originVector,
                                    const std::shared_ptr<VectorEpetra> statusVector,
                                    const std::shared_ptr<VectorEpetra> saveVector,
                                    const UInt offset)
{
    UInt dim = dispFESpace->dim();
    UInt nTotDof = dispFESpace->dof().numTotalDof();

    // We loop over the local ID on the processors of the originVector
    for( UInt i(0); i < nTotDof; ++i )
    {
        if( originVector.blockMap().LID( static_cast<EpetraInt_Type>(i) ) != -1 ) // The i-th dof is on the processor
        {
            // We look at the functor to see if the condition is satified
            bool variable = functor( i, nTotDof, offset );

            if ( variable )
            {
                // At the first time we insert the value and then we "close" the cell
                for ( UInt iComp = 0; iComp < dispFESpace->fieldDim(); ++iComp )
                {
                    Int LIDid = originVector.blockMap().LID (static_cast<EpetraInt_Type>(i + iComp * dim + offset));
                    Int GIDid = originVector.blockMap().GID (LIDid);

                    if( (*statusVector)( GIDid ) != 1.0 )
                    {
                        // Saving the value
                        (*saveVector)( GIDid ) = originVector( GIDid );

                        // Saying it has been saved at that point
                        (*statusVector) ( GIDid ) = 1.0;
                    }
                }
            }
        }
    }
}

template<typename FunctorType, typename MeshType, typename MapType>
void saveBooleanVectorAccordingToFunctor ( const std::shared_ptr<FESpace<MeshType, MapType> > dispFESpace,
                                           const FunctorType& functor,
                                           const std::shared_ptr<VectorEpetra> originVector,
                                           std::shared_ptr<VectorEpetra> saveVector,
                                           const UInt offset)
{
    UInt dim = dispFESpace->dim();
    UInt nTotDof = dispFESpace->dof().numTotalDof();

    // We loop over the local ID on the processors of the originVector
    for( UInt i(0); i < nTotDof; ++i )
    {
        if( originVector->blockMap().LID( static_cast<EpetraInt_Type>(i) ) != -1 ) // The i-th dof is on the processor
        {
            // We look at the functor to see if the condition is satified
            bool variable = functor( i, nTotDof, offset );

            if ( variable )
            {
                // At the first time we insert the value and then we "close" the cell
                for ( UInt iComp = 0; iComp < dispFESpace->fieldDim(); ++iComp )
                {
                    Int LIDid = originVector->blockMap().LID (static_cast<EpetraInt_Type>(i + iComp * dim + offset));
                    Int GIDid = originVector->blockMap().GID (LIDid);

                    (*saveVector)( GIDid ) = 1.0;
                }
            }
        }
    }
}

//! Gradient of the displacement on the local element
/*!
  This function assembles the local tensor of the gradient of the displacement field

  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param fe The current finite element
*/
void computeGradientLocalDisplacement (boost::multi_array<Real, 3>& gradientLocalDisplacement, const VectorElemental& uk_loc, const CurrentFE& fe );

//! Deformation Gradient on the local element
/*!
  This function assembles the local deformation gradient

  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param fe The current finite element
*/
void computeLocalDeformationGradient (const VectorElemental& uk_loc, std::vector<Epetra_SerialDenseMatrix>& tensorF, const CurrentFE& fe );

//! Gradient on the local element
/*!
  This function assembles the local  gradient

  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param fe The current finite element
*/
void computeLocalDeformationGradientWithoutIdentity (const VectorElemental& uk_loc, std::vector<Epetra_SerialDenseMatrix>& tensorF, const CurrentFE& fe );


//! METHODS SHARED BETWEEN LINEAR ELASTIC MODEL AND ST.VENANT-KIRCHHOFF MODEL
//! These two methods are implemented in AssemblyElemental.cpp.
//void stiff_strain( Real coef, MatrixElemental& elmat, const CurrentFE& fe );
//void stiff_div( Real coef, MatrixElemental& elmat, const CurrentFE& fe );


//! METHODS FOR ST.VENANT-KIRCHHOFF MODEL
//! Methods for the Stiffness matrix ( evaluate the RHS or the residuals )

//! Elementary first term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local first term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_derdiv ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local second term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergradbis ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local third term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_divgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fourth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the fourth term of the nonlinear stiffness matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary stiffness matrix of the current volume
  @param fe The current finite element
*/
void stiff_gradgrad ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergrad_gradbis ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth-2 term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth-2 term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergrad_gradbis_Tr ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary sixth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local sixth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_gradgradTr_gradbis ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );




//! METHODS FOR THE JACOBIAN MATRIX

//! Elementary first term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local first term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergrad ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local second term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_divgrad_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local third term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_gradgrad_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fourth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergrad_gradbis_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_dergrad_gradbis_Tr_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary sixth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local sixth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_gradgradTr_gradbis_2 ( Real coef, const boost::multi_array<Real, 3>& gradientLocalDisplacement, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary seventh term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local seventh term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_gradgradTr_gradbis_3 ( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );



//! METHODS SHARED BETWEEN NEO-HOOKEAN AND EXPONENTIAL MODELS
//! The volumetric part for Neo-Hookean and Exponential models is the same
//! Methods for the volumetric part of the stiffness vector

//! Elementary volumetric term of the nonlinear stiffness vector of the Neo-Hookean and Exponential models (see the reference)
/*!
  This function assembles the volumetric term of the nonlinear stiffness vector of the Neo-Hookean and Exponential models.

  @param coef The constant coefficient of the volumetric term of the nonlinear stiffness vector
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elvec The elementary vector of the current volume
  @param fe The current finite element
*/
void source_Pvol (Real coef, const boost::multi_array<Real, 3>&  CofFk, const std::vector<Real>& Jk, VectorElemental& elvec, const CurrentFE& fe);

//! Methods for the volumetric part of the Jacobian matrix

//! Elementary first volumetric term of the nonlinear Jacobian matrix of the Neo-Hookean and Exponential models (see the reference)
/*!
  This function assembles the local first volumetric term of the Jacobian matrix of the nonlinear volumetric stiffness vector of the Neo-Hookean and Exponential models.

  @param coef The constant coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elvec The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_Pvol_1term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second volumetric term of the nonlinear Jacobian matrix of the Neo-Hookean and Exponential models (see the reference)
/*!
  This function assembles the local second volumetric term of the Jacobian matrix of the nonlinear volumetric stiffness vector of the Neo-Hookean and Exponential models.

  @param coef The constant coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_Pvol_2term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );

//! METHODS FOR NEO-HOOKEAN MODEL
//! Methods for the isochoric part of the stiffness vector

//! Elementary nonlinear isochoric stiffness vector for Neo-Hookean model (see the reference)
/*!
  This function assembles the nonlinear isochoric part of the stiffness vector for Neo-Hookean model.

  @param coef The coefficient of the vector
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elvec The elementary vector of the current volume
  @param fe The current finite element
*/
void source_P1iso_NH (Real coef, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, VectorElemental& elvec, const CurrentFE& fe);

//! METHODS FOR TENSORIAL CALCULUS
//! In this part of the namespace, the methods to perform basics operations on tensors are defined.
//! The operations are: tensorial products, computations of invariants
/*!
  This function computes the invariants of the right Cauchy Green tensor and the cofactor of F

  @param invariants vector of invariants of C
  @param tensorF deformation gradient tensor
*/
void computeInvariantsRightCauchyGreenTensor (std::vector<LifeV::Real>& invariants,
                                              const Epetra_SerialDenseMatrix& tensorF,
                                              Epetra_SerialDenseMatrix& cofactorF);

void computeInvariantsRightCauchyGreenTensor (std::vector<LifeV::Real>& invariants,
                                              const Epetra_SerialDenseMatrix& tensorF);


/*!
  This function computes the Cauchy stress tensor given the detF, first Piola Kirchhoff and tensorF

  @param cauchy Cauchy stress tensor
  @param firstPiola first Piola-Kirchhoff tensor
  @param invariants vector of invariants of C
  @param tensorF deformation gradient tensor
*/
void computeCauchyStressTensor (Epetra_SerialDenseMatrix& cauchy,
                                Epetra_SerialDenseMatrix& firstPiola,
                                LifeV::Real det,
                                Epetra_SerialDenseMatrix& tensorF);

/*!
  This function computes the eigenvalues of \sigma
  @param cauchy Cauchy stress tensor
  @param eigenvalues vector of principal tensions
*/
void computeEigenvalues (const Epetra_SerialDenseMatrix& cauchy,
                         std::vector<LifeV::Real>& eigenvaluesR,
                         std::vector<LifeV::Real>& eigenvaluesI);

//! Methods for the isochoric part of the Jacobian matrix

//! Elementary first nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local first nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_NH_1term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary second nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local second nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_NH_2term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local third nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_NH_3term ( Real coef, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local fourth nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_NH_4term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local fifth nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_NH_5term ( Real coef, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! METHODS FOR EXPONENTIAL MODEL
//! Methods for the isochoric part of the stiffness vector

//! Elementary nonlinear isochoric stiffness vector for Exponential model (see the reference)
/*!
  This function assembles the local nonlinear isochoric stiffness vector for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elvec The elementary vector of the current volume
  @param fe The current finite element
*/
void source_P1iso_Exp ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, VectorElemental& elvec, const CurrentFE& fe );

//! Methods for the isochoric part of the Jacobian matrix

//! Elementary first nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local first nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_1term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local second nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_2term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local third nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_3term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local fourth nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_4term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary fifth nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local fifth nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_5term ( Real coef, Real coefExp, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary sixth nonlinear isochoric Jacobian matrix for Exponential model (see the reference)
/*!
  This function assembles the local sixth nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_Exp_6term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! METHOD FOR THE ST. VENANT-KIRCHHOFF LAW IN THE ISOCHORIC AND VOLUMETRIC DECOUPLED VERSION
//! This law uses, for the volumetric part, the one used by the exponential and neohookean model.
//! The isochoric part is, of course, completely new and, consequently its jacobian.

//! Elementary nonlinear isochoric stiffness vector for St. Venant-Kirchhoff Penalized model
/*!
  This function assembles the local nonlinear isochoric stiffness vector for St. Venant-Kirchhoff Penalized model

  @param lambda first Lame Coefficient
  @param mu shear modulus
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elvec The elementary vector of the current volume
  @param fe The current finite element
*/
void source_P1iso_VKPenalized ( Real lambda, Real mu, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Ic_isok, const std::vector<Real>& Ic_k, const std::vector<Real>& Jack_k, VectorElemental& elvec, const CurrentFE& fe );

//! Elementary nonlinear isochoric stiffness vector for St. Venant-Kirchhoff Penalized model
/*!
  This function assembles the local nonlinear isochoric stiffness vector for St. Venant-Kirchhoff Penalized model

  @param mu shear modulus
  @param Fk^-T
  @param FkCk The product FC
  @param Jk The determinant of F
  @param Ic_Squared The trace of C^2
  @param elvec The elementary vector of the current volume
  @param fe The current finite element
*/
void  source_P2iso_VKPenalized ( Real mu, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& FkCk, const std::vector<Real>&   Ic_Squared, const std::vector<Real>&   Jk, VectorElemental& elvec, const CurrentFE& fe );


//! Methdos for the Jacobian of the St. Venant-Kirchhoff Penalized law.
//! Elementary first nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local first nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param lambda the first Lame coefficient
  @param FkMinusTransposed The matrix F^-T
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_k The trace of C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_0term ( Real lambda, Real mu, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_k, const std::vector<Real>& IcIso_k, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary first nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local first nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param lambda the first Lame coefficient
  @param FkMinusTransposed The matrix F^-T
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_k The trace of C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_1term ( Real coeff, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_k, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary third nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local third nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda / 2 where lambda is the first Lame constant
  @param FKMinusTransposed The matrix F^-T
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_2term ( Real coef, const boost::multi_array<Real, 3 >& FkMinusTransposed, const std::vector<Real>& Jk, const std::vector<Real>& Ic_k, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local second nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef (lambda/2) where lambda is the first Lame coefficient
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_3term ( Real coef, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary fourth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local fourth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef (lambda/2) where lamba is the first Lame coefficient
  @param FkMinusTransposed The tensor F^-T
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_4term ( Real coef, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_k, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary fifth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local fifth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda first lame coefficient
  @param coef2 mu shear modulus
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_5term ( Real coef, Real secondCoef, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary sixth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local sixth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda first Lame coefficient
  @param coef2 mu sheat modulus
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_6term ( Real coef, Real secondCoef, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, const boost::multi_array<Real, 3 >& Fk, const boost::multi_array<Real, 3 >& FkMinusTransposed, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary seventh  nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local seventh nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda the first Lame coefficient
  @param coef2 mu the shear modulus
  @param FkMinusTransposed The tensor F^-T
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_7term ( Real coef, Real secondCoef, const boost::multi_array<Real, 3 >& FkMinusTransposed, const std::vector<Real>& Ic_isok, const std::vector<Real>& Ic_k, const std::vector<Real>& Jk, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary seventh  nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local seventh nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda the first Lame coefficient
  @param coef2 mu the shear modulus
  @param FkMinusTransposed The tensor F^-T
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_8term ( Real coef, const std::vector<Real>& Jack_k, const boost::multi_array<Real, 3 >& FkMinusTransposed, const boost::multi_array<Real, 3 >& FkCk, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary seventh  nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local seventh nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef lambda the first Lame coefficient
  @param coef2 mu the shear modulus
  @param FkMinusTransposed The tensor F^-T
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param Ic_k The first invariant of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_9term ( Real coef, const std::vector<Real>& Jack_k, const std::vector<Real>& Ic_kSquared, const boost::multi_array<Real, 3 >& FkMinusTransposed, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary eigth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local eigth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef mu the shear modulus
  @param Ck The right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_10term ( Real coef, const std::vector<Real>& Jack_k, const boost::multi_array<Real, 3 >& Ck, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary sixth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local sixth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef mu the shear modulus
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation tensor
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_11term ( Real coef, const std::vector<Real>& Jk, const boost::multi_array<Real, 3 >& Fk, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary tenth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local tenth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef mu the shear modulus
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation tensor
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_12term ( Real coef, const std::vector<Real>& Jk, const boost::multi_array<Real, 3 >& Fk, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary eleventh nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local eleventh nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef mu the shear modulus
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the right Cauchy-Green tensor C^2
  @param FkMinusTransposed The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_13term ( Real coef, const std::vector<Real>& Jk, const std::vector<Real>& Ic_kSquared,  const boost::multi_array<Real, 3 >& FkMinusTransposed, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary twelveth nonlinear isochoric Jacobian matrix for VK-Penalized model (see the reference)
/*!
  This function assembles the local twelveth nonlinear isochoric Jacobian matrix for VK-Penalized model.

  @param coef mu the shear modulus
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ck The first invariant of the right Cauchy-Green tensor C^2
  @param Fk The deformation gradient matrix
  @param FkMinusTransposed The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_VKPenalized_14term ( Real coef, const std::vector<Real>& Jk, const boost::multi_array<Real, 3 >& FkCk, const boost::multi_array<Real, 3 >& FkMinusTransposed, MatrixElemental& elmat, const CurrentFE& fe );

//! Methods for second order exponential law
//! Methods for the first Piola-Kirchhoff tensor
void source_P1iso_SecondOrderExponential ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& trCisok, VectorElemental& elvec, const CurrentFE& fe );

//! Methods for the Jacobian matrix
//! Elementary first nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local first nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_1term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local second nonlinear isochoric Jacobian matrix for Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_2term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary third nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local third nonlinear isochoric Jacobian matrix for Second Order Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_3term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary fourth nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local fourth nonlinear isochoric Jacobian matrix for Second Order Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Fk The deformation gradient that depends on the local displacement uk_loc
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param Ic_k The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_4term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const boost::multi_array<Real, 3 >& Fk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, const std::vector<Real>& Ic_k, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary fifth nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local fifth nonlinear isochoric Jacobian matrix for Second Order Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_5term ( Real coef, Real coefExp, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

//! Elementary sixth nonlinear isochoric Jacobian matrix for Second Order Exponential model (see the reference)
/*!
  This function assembles the local sixth nonlinear isochoric Jacobian matrix for Second Order Exponential model.

  @param coef The pre-exponential coefficient
  @param coefExp The expoenential coefficient
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Ic_isok The first invariant of the isochoric part of the right Cauchy-Green tensor C
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
*/
void stiff_Jac_P1iso_SecondOrderExp_6term ( Real coef, Real coefExp, const boost::multi_array<Real, 3 >& CofFk, const std::vector<Real>& Jk, const std::vector<Real>& Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

} //! End namespace AssemblyElementalStructure

} //! End namespace LifeV
#endif
