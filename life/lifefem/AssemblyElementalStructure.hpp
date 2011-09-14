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

    @contributor Gianmarco Mengaldo <gianmarco.mengaldo@gmail.ch>
    @mantainer Paolo Tricerri <paolo.tricerri@epfl.ch>

    All the methods are described in Mengaldo G. Master Thesis "Nonlinear Fluid-Structure Interaction with applications in computational haemodynamics" Politecnico di Milano 2011
 */


#ifndef _ELEMOPERSTRUCTURE_H_INCLUDED
#define _ELEMOPERSTRUCTURE_H_INCLUDED

#include <life/lifearray/MatrixElemental.hpp>
#include <life/lifearray/VectorElemental.hpp>

#include <life/lifecore/LifeV.hpp>

#include <life/lifefem/CurrentBoundaryFE.hpp>
#include <life/lifefem/CurrentFE.hpp>
#include <life/lifefem/DOF.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

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

//! Elementary mass for constant mass coefficient
/*!
  This function assembles the local mass matrix when the mass coefficient is constant.

  @param coef The constant coefficient of the matrix
  @param elmat The elementary mass matrix of the current volume
  @param fe The current finite element
  @param iblock The component of v that is concerned
  @param jblock The component of v that is concerned
  @param nb The number of the dimensions of the problem
 */
void mass( Real coef, MatrixElemental& elmat, const CurrentFE& fe, int iblock, int jblock, UInt nb );


//! Elementary mass for not-constant mass coefficient
/*!
  This function assembles the local mass matrix when the mass coefficient is variable in time and in space.

  @param qpt_coef The not-constant coefficient of the matrix
  @param elmat The elementary mass matrix of the current volume
  @param fe The current finite element
  @param iblock The component of v that is concerned
  @param jblock The component of v that is concerned
  @param nb The number of the dimensions of the problem
 */
void mass( const std::vector<Real>& qpt_coef, MatrixElemental& elmat, const CurrentFE& fe, int iblock, int jblock, UInt nb );


  
//! METHODS SHARED BETWEEN LINEAR ELASTIC MODEL AND ST.VENANT-KIRCHHOFF MODEL

//! Elementary first term of the linear stiffness matrix for linear elasticity (see the reference)
/*!
  This function assembles the local first term of the linear stiffness matrix for linear elasticity.

  @param coef The constant coefficient of the matrix 
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_strain( Real coef, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second term of the linear stiffness matrix for linear elasticity (see the reference)
/*!
  This function assembles the local second term of the linear stiffness matrix for linear elasticity.

  @param coef The constant coefficient of the matrix 
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_div( Real coef, MatrixElemental& elmat, const CurrentFE& fe );



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
void stiff_derdiv( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local second term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local third term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_divgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fourth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the fourth term of the nonlinear stiffness matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary stiffness matrix of the current volume
  @param fe The current finite element
 */
void stiff_gradgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergrad_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth-2 term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth-2 term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergrad_gradbis_Tr( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary sixth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local sixth term of the nonlinear stiffness matrix for St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_gradgradTr_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );




//! METHODS FOR THE JACOBIAN MATRIX

//! Elementary first term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local first term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local second term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_divgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local third term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_gradgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fourth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fourth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergrad_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary fifth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local fifth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_dergrad_gradbis_Tr_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary sixth term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local sixth term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_gradgradTr_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary seventh term of the Jacobian matrix for the nonlinear stiffness matrix of the St.Venant-Kirchhoff model (see the reference)
/*!
  This function assembles the local seventh term of the Jacobian matrix of the nonlinear stiffness matrix of the St.Venant-Kirchhoff model.

  @param coef The constant coefficient of the matrix 
  @param uk_loc The local displacement (remark: the nonlinear matrix depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_gradgradTr_gradbis_3( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );



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
void source_Pvol(Real coef, const  KNMK<Real> CofFk, const KN<Real> Jk, VectorElemental& elvec, const CurrentFE& fe);

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
void stiff_Jac_Pvol_1term( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary second volumetric term of the nonlinear Jacobian matrix of the Neo-Hookean and Exponential models (see the reference)
/*!
  This function assembles the local second volumetric term of the Jacobian matrix of the nonlinear volumetric stiffness vector of the Neo-Hookean and Exponential models.

  @param coef The constant coefficient of the matrix
  @param CofFk The cofactor of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_Jac_Pvol_2term( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );



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
void source_P1iso_NH(Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, VectorElemental& elvec, const CurrentFE& fe); 



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
void stiff_Jac_P1iso_NH_1term( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_NH_2term( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


//! Elementary third nonlinear isochoric Jacobian matrix for Neo-Hookean model (see the reference)
/*!
  This function assembles the local third nonlinear isochooric Jacobian matrix for Neo-Hookean model.

  @param coef The coefficient of the matrix
  @param Jk The determinant of the deformation gradient F that depends on the local displacement uk_loc (remark: the nonlinear vector depends on current displacement)
  @param elmat The elementary matrix of the current volume
  @param fe The current finite element
 */
void stiff_Jac_P1iso_NH_3term( Real coef, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_NH_4term( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_NH_5term( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );



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
void source_P1iso_Exp( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, VectorElemental& elvec, const CurrentFE& fe );

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
void stiff_Jac_P1iso_Exp_1term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_Exp_2term( Real coef, Real coefExp, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_Exp_3term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_Exp_4term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk,const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_Exp_5term( Real coef, Real coefExp, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );


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
void stiff_Jac_P1iso_Exp_6term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

} //! End namespace AssemblyElementalStructure

} //! End namespace LifeV
#endif
