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
    @file elemOperHdiv.hpp
    @brief functions for Hdiv that were in elemOper.hpp

    @author Simone Deparis <simone.deparis@epfl.ch>
    @date 09 Mar 2010

 */

#ifndef ELEMOPERHDIV_H
#define ELEMOPERHDIV_H 1

#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/currentHdivFE.hpp>
#include <life/lifefem/refHybridFE.hpp>

namespace LifeV {

  //----------------------------------------------------------------------
  //
  //!@name  Operators for H(div) finite elements
  //!@{
  //
  //----------------------------------------------------------------------
  //!Transpose of Elementary divergence matrix
  void grad_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
		  const CurrentFE& fe_p, int iblock, int jblock );

  //!Elementary divergence matrix
  void div_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
		 const CurrentFE& fe_p, int iblock, int jblock );

  //!Elementary tp v dot n matrix
  void TP_VdotN_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe, const RefHybridFE& vdotnfe, int iblock, int jblock );

  //!Elementary tp^2  matrix
  void TP_TP_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe, int iblock, int jblock );

//!@}
//!@name Mass matrix
//!@{
  //-------------Mass matrix---------------------------------------
  /*!
    Weighted Mass matrix with a permeability tensor which is a constant scalar matrix
    (i.e. \f$K^{-1} = coef \cdot Id\f$, coef being the inverse of the permeability).
  */
  void mass_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe,
		  int iblock = 0, int jblock = 0 );

  //-------------Mass matrix---------------------------------------
  /*!
    Weighted Mass matrix with permeability matrix which is a constant
    per element symmetric positive definite matrix (non diagonal a
    priori) and ALREADY INVERTED (with Lapack LU or Choleski for
    instance).
  */
  void mass_Hdiv( Matrix const& Invperm, ElemMat& elmat, const CurrentHdivFE& fe,
		  int iblock = 0, int jblock = 0 );


  /*!
    Weighted Mass matrix with a permeability that is a scalar function.
    The inverse function of the permeability should be provided.
  */
  void mass_Hdiv( Real ( *Invperm ) ( const Real&, const Real&, const Real& ),
		  ElemMat& elmat, const CurrentHdivFE& fe, int iblock, int jblock );

  /*!
    Weighted Mass matrix with a permeability which is a constant scalar
    (i.e. K^{-1} = coef, coef is the inverse of the permeability).
  */
  void mass_Mixed_Hdiv( Real coef, ElemMat& elmat, const CurrentFE& fe,
			const CurrentHdivFE& hdivfe, int iblock, int jblock );

//!@}

} // Namespace LifeV

#endif /* ELEMOPERHDIV_H */
