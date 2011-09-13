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
 *  @file
 *  @brief  This class has been created only to split the methods for structural problems from the others methods of AssemblyElemental. It contains all the methods are  private. All the methods have been developed to implement the materials which are in LifeV at the moment of writing: linear elastic, Venant-Kirchhoff, neohookean and exponential.
 *
 *  @version 1.0
 *  @date 28-08-2011
 *  @author Paolo Tricerri
 *  @author Gianmarco Mengaldo
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _STRUCTURALASSEMBLER_H_
#define _STRUCTURALASSEMBLER_H_

#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/LifeV.hpp>

#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/AssemblyElemental.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{
//! StructuralAssembler - This class handle the methods to compute the First Piola Kirchhoff tensor of different materials. In particular, linear elastic, St. Venant-Kirchhoff Neohookean and Exponential. Update the documentation!!!

class StructuralAssembler
{
public:

  //! @name Constructors & Deconstructor
  //@{
  StructuralAssembler();

  ~StructuralAssembler();
  //@}


  //! @name Public Methods
  //@{
  
  //! METHODS FOR LINEAR ELASTIC MODEL
  void stiff_strain( Real coef, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_div( Real coef, MatrixElemental& elmat, const CurrentFE& fe );



  //! METHODS FOR ST.VENANT-KIRCHHOFF MODEL
  //! Methods for the Stiffness matrix ( evaluate the RHS or the residuals )
  void stiff_derdiv( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_dergradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_divgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_gradgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_dergrad_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_dergrad_gradbis_Tr( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_gradgradTr_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //! Methods for the jacobian stiffness matrix ( called in StructuralMaterial::updateJacobian )
  // These first methods are  already implemented for the stiffness matrix
  //void stiff_derdiv( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //void stiff_divgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //void stiff_gradgrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //void stiff_dergrad_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //void stiff_dergrad_gradbis_Tr( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  //void stiff_gradgradTr_gradbis( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_dergrad( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_divgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_gradgrad_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_dergrad_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_dergrad_gradbis_Tr_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_gradgradTr_gradbis_2( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );
  void stiff_gradgradTr_gradbis_3( Real coef, const VectorElemental& uk_loc, MatrixElemental& elmat, const CurrentFE& fe );



  //! METHODS SHARED BETWEEN NEO-HOOKEAN AND EXPONENTIAL MODELS
  //! The volumetric part for Neo-Hookean and Exponential models is the same
  //! Methods for the volumetric part of the stiffness vector
  void source_Pvol(Real coef, const  KNMK<Real> CofFk, const KN<Real> Jk, VectorElemental& elvec, 
		    const CurrentFE& fe);


  //! Methods for the volumetric part of the Jacobian matrix
  void stiff_Jac_Pvol_1term( Real coef,  const KNMK<Real> CofFk,  const KN<Real> Jk, 
			      MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_Pvol_2term( Real coef, const KNMK<Real> CofFk,  const KN<Real> Jk, 
  			      MatrixElemental& elmat, const CurrentFE& fe );



  //! METHODS FOR NEO-HOOKEAN MODEL
  //! Methods for the isochoric part of the stiffness vector
  void source_P1iso_NH(Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok , VectorElemental& elvec, const CurrentFE& fe); 

  //! Methods for the isochoric part of the Jacobian matrix
  void stiff_Jac_P1iso_NH_1term( Real coef,  const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk ,  MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_NH_2term( Real coef,  const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_NH_3term( Real coef,  const KN<Real> Jk, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_NH_4term( Real coef, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk , MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_NH_5term( Real coef, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );



  //! METHODS FOR EXPONENTIAL MODEL
  //! Methods for the isochoric part of the stiffness vector
  void source_P1iso_Exp( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, VectorElemental& elvec, const CurrentFE& fe );

   //! Methods for the isochoric part of the Jacobian matrix
  void stiff_Jac_P1iso_Exp_1term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_Exp_2term( Real coef, Real coefExp, const KNMK<Real> Fk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_Exp_3term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_Exp_4term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KNMK<Real> Fk,const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_Exp_5term( Real coef, Real coefExp, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  void stiff_Jac_P1iso_Exp_6term( Real coef, Real coefExp, const KNMK<Real> CofFk, const KN<Real> Jk, const KN<Real> Ic_isok, MatrixElemental& elmat, const CurrentFE& fe );

  //@}


};
}
#endif /* _STRUCTURALASSEMBLER_H_ */
