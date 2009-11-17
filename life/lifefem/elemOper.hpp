/*-*- mode: c++ -*-
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/*!@file elemOper.hpp
   Element matrix operations
*/

#ifndef _ELEMOPER_H_INCLUDED
#define _ELEMOPER_H_INCLUDED
#include <life/lifecore/life.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/currentBdFE.hpp>
#include <life/lifefem/currentHdivFE.hpp>
#include <life/lifefem/refHybridFE.hpp>
#include <life/lifefem/dof.hpp>

namespace LifeV
{
  //----------------------------------------------------------------------
  //
  //!@name               Operators for classical finite elements
  //@{
  //----------------------------------------------------------------------
  //!coef(t,x,y,z,u)
  void mass( Real (*coef)(Real,Real,Real,Real,Real),
	     ElemMat& elmat, const CurrentFE& fe,
	     const Dof& dof,
	     const ScalUnknown<Vector>& U,Real t);

  void stiff( Real (*coef)(Real,Real,Real,Real,Real),
	      ElemMat& elmat, const CurrentFE& fe,
	      const Dof& dof,
	      const ScalUnknown<Vector>& U,Real t);
  void source( Real (*fct)(Real,Real,Real,Real,Real),
	       ElemVec& elvec, const CurrentFE& fe,
	       const Dof& dof,
	       const ScalUnknown<Vector>& U,Real t);

  void mass( Real coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock = 0, int jblock = 0 );
  void mass( Real coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock, int jblock, int nb );
  //! Mass term with coefficients given for each quadrature point
  void mass( const std::vector<Real>& qpt_coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock, int jblock, int nb );

  void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
	      int iblock = 0, int jblock = 0 );
  void stiff( Real coef, Real ( *fct ) ( Real, Real, Real ), ElemMat& elmat,
	      const CurrentFE& fe, int iblock, int jblock );
  void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
	      int iblock, int jblock, int nb );
  //! Stiff term with coefficient given for each quadrature node
  void stiff( const std::vector<Real>& qpt_coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock, int jblock, int nb );
  void stiff_curl( Real coef, ElemMat& elmat, const CurrentFE& fe,
		   int iblock, int jblock, int nb );





  //! \f$ coef \cdot ( e(u) , e(v) )\f$
  void stiff_strain( Real coef, ElemMat& elmat, const CurrentFE& fe );

  //! \f$ coef \cdot ( div u , div v )\f$
  void stiff_div( Real coef, ElemMat& elmat, const CurrentFE& fe );

  //! \f$ coef \cdot ( [\nabla u^k]^T \nabla u : \nabla v  )\f$
  void stiff_dergradbis( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe );

  //! \f$ coef \cdot ( [\nabla u]^T \nabla u^k + [\nabla u^k]^T \nabla u : \nabla v  )\f$ for Newton on St-Venant
  void stiff_dergrad( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe );

  //! \f$ coef \cdot ( \tr { [\nabla u^k]^T \nabla u }, \nabla\cdot  v  ) \f$ for Newton on St-Venant
  void stiff_derdiv( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe );




  void grad( const int icoor, Real coef, ElemMat& elmat,
	     const CurrentFE& fe_u, const CurrentFE& fe_p,
	     int iblock = 0, int jblock = 0 );
  void div( const int icoor, Real coef, ElemMat& elmat,
	    const CurrentFE& fe_u, const CurrentFE& fe_p,
	    int iblock = 0, int jblock = 0 );
  void grad_div( Real coef_grad, Real coef_div, ElemMat& elmat,
		 const CurrentFE& fe_u, const CurrentFE& fe_p,
		 int block_pres );

  void stab_stokes( Real visc, Real coef_stab, ElemMat& elmat,
		    const CurrentFE& fe, int block_pres );
  void advection( Real coef, ElemVec& vel, ElemMat& elmat,
		  const CurrentFE& fe, int iblock, int jblock, int nb );

  void source( Real constant, ElemVec& elvec, const CurrentFE& fe, int iblock );

  void source( Real constant, ElemVec& elvec, const CurrentFE& fe, Real t, int iblock );


  // right-hand sides for Chorin-Teman projection scheme
  void source_divuq(Real alpha, ElemVec& uLoc,  ElemVec& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock = 0 );
  void source_gradpv(Real alpha, ElemVec& pLoc,  ElemVec& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock );


  //!@}
  //!@name Elementary operations for the interior penalty stabilization
  //!@{
  //

  //! \f$ coef < \nabla p1, \nabla q2 >\f$
  void ipstab_grad( const Real coef, ElemMat& elmat,
		    const CurrentFE& fe1, const CurrentFE& fe2,
		    const CurrentBdFE& bdfe, int iblock = 0, int jblock = 0 );

//! \f$ coef < \nabla u1, \nabla v2 >\f$
  void ipstab_grad( const Real coef, ElemMat& elmat,
		    const CurrentFE& fe1, const CurrentFE& fe2,
		    const CurrentBdFE& bdfe, int iblock, int jblock, int nb );

//! \f$ coef < \nabla\cdot  u1, \nabla\cdot  v2 >\f$
  void ipstab_div( const Real coef, ElemMat& elmat,
		   const CurrentFE& fe1, const CurrentFE& fe2,
		   const CurrentBdFE& bdfe, int iblock = 0, int jblock = 0 );
//! \f$ coef < \beta1 . \nabla u1, \beta2 . \nabla v2 >\f$
  void ipstab_bgrad( const Real coef, ElemMat& elmat,
		     const CurrentFE& fe1, const CurrentFE& fe2,
		     const ElemVec& beta, const CurrentBdFE& bdfe,
		     int iblock, int jblock, int nb );
//! \f$ coef < |\beta . n|^2 / |\beta| \nabla p1, \nabla q2 >\f$
  void ipstab_bagrad( const Real coef, ElemMat& elmat,
		      const CurrentFE& fe1, const CurrentFE& fe2,
		      const ElemVec& beta, const CurrentBdFE& bdfe,
		      int iblock = 0, int jblock = 0 );

//!\f$ coef < |\beta\cdot n| \nabla p1, \nabla q2 >\f$
//! p1 lives in fe1
//! q2 lives in fe2
//! beta lives in fe3

void ipstab_bagrad( const Real           coef,
                    ElemMat&             elmat,
                    const CurrentFE&     fe1,
                    const CurrentFE&     fe2,
                    const CurrentFE&     fe3,
                    const ElemVec&       beta,
                    const CurrentBdFE&   bdfe,
                    int iblock = 0, int jblock = 0 );

//!@}
///////////////////////////////////////


  ////////////////////////////////////////
  //! Convective term with a local vector coefficient (useful for Navier-Stokes problem)
  void grad( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
	     const CurrentFE& fe1, const CurrentFE& fe2,
	     int iblock, int jblock );

  //! Convective term with a local vector coefficient (useful for Navier-Stokes problem+adv-diff)
  void grad( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
	     const CurrentFE& fe1, const CurrentFE& fe2, const CurrentFE& fe3,
	     int iblock = 0, int jblock = 0 );

  //! Convective term with a local vector coefficient for Navier-Stokes problem in Skew-Symmetric form
  void grad_ss( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
		const CurrentFE& fe1, const CurrentFE& fe2,
		int iblock = 0, int jblock = 0 );

  //! StreamLine Diffusion
  void stiff_sd( Real coef, const ElemVec& vec_loc, ElemMat& elmat, const CurrentFE& fe,
		 const CurrentFE& fe2, int iblock = 0, int jblock = 0, int nb = 1 );

  /////////////////////////////////////////////
  //
  //! source  \f$ \int fct \phi_i \f$
  //
  template <typename UsrFct>
  void source( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, int iblock = 0 )
  {
    ASSERT_PRE( fe.hasQuadPtCoor(), "Source with space dependent fonction need updated quadrature point coordinates. Call for example updateFirstDerivQuadPt() instead of updateFirstDeriv()." );
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    for ( i = 0;i < fe.nbNode;i++ )
      {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
	  {
            s += fe.phi( i, ig ) * fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ),
                                        iblock ) * fe.weightDet( ig );
	  }
        vec( i ) += s;
      }
  }

  /////////////////////////////////////////////
  //
  //! source  \f$ \int fct(t) \phi_i\f$
  //
  template <typename UsrFct>
  void source( const UsrFct& fct, ElemVec& elvec, const CurrentFE& fe, Real t, int iblock = 0 )
  {
    ASSERT_PRE( fe.hasQuadPtCoor(), "Source with space dependent fonction need updated quadrature point coordinates. Call for example updateFirstDerivQuadPt() instead of updateFirstDeriv()." );
    int i, ig;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    for ( i = 0;i < fe.nbNode;i++ )
      {
        s = 0;
        for ( ig = 0;ig < fe.nbQuadPt;ig++ )
	  {
            s += fe.phi( i, ig ) * fct( fe.quadPt( ig, 0 ), fe.quadPt( ig, 1 ), fe.quadPt( ig, 2 ), t,
                                        iblock ) * fe.weightDet( ig );
	  }
        vec( i ) += s;
      }
  }


  void source( Real coef, ElemVec& f, ElemVec& elvec, const CurrentFE& fe,
	       int fblock = 0, int eblock = 0 );

  void source_fhn( Real coef_f, Real coef_a, ElemVec& u, ElemVec& elvec, const CurrentFE& fe,
		   int fblock = 0, int eblock = 0 );

//!@name Shape derivative terms for Newton FSI
//!@{
  //
  //! \f$ coef \cdot ( \nabla (-w^k):[I\nabla\cdot  d - (\nabla d)^T] u^k + convect^T[I\nabla\cdot  d - (\nabla d)^T] (\nabla u^k)^T , v  ) \f$ for Newton FSI
  //
  //! Remark: convect = \f$u^n-w^k\f$
  //
  void source_mass1( Real coef,
                     const ElemVec& uk_loc,
                     const ElemVec& wk_loc,
                     const ElemVec& convect_loc,
                     const ElemVec& d_loc,
                     ElemVec& elvec,
                     const CurrentFE& fe );

  //
  //! \f$coef \cdot ( \nabla u^k dw, v  )\f$ for Newton FSI
  //
  //
  void source_mass2( Real coef, const ElemVec& uk_loc, const ElemVec& dw_loc,
		     ElemVec& elvec, const CurrentFE& fe );



  void  source_mass3( Real coef, const ElemVec& un_loc, const ElemVec& uk_loc, const ElemVec& d_loc,
		      ElemVec& elvec, const CurrentFE& fe );





  //
  //! \f$coef \cdot ( [-p^k I + 2\mu e(u^k)] [I\nabla\cdot d - (\nabla d)^T] , \nabla v  )\f$ for Newton FSI
  //
  void source_stress( Real coef, Real mu, const ElemVec& uk_loc, const ElemVec& pk_loc,
		      const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u,
		      const CurrentFE& fe_p );

  //
  //! \f$+ \mu ( \nabla u^k \nabla d + [\nabla d]^T[\nabla u^k]^T : \nabla v )\f$
  //
  void source_stress2( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u );



/**
 *Shape terms for the CE system in FSI Newton. It is the sum of the terms \c source_mass1\c, \c source_stress\c and \c source_stress2\c.
 *the difference between this method and the previous ones is that here the shape terms are assembled in a matrix, instead
 *of a vector. This implies some extra loop and the explicit construction of the tensor \f$[I\nabla\cdot  d - (\nabla d)^T]\f$
 *instead of it's multiplication times a vector. However the computation of these terms need to be done once per Newton
 *iterations, instead of once per Jacobian-vector multiplication.

 *Note that the term \c source_mass2\c is not considered here because the fluid domain velocity \f w\f is trated explicitly
 */
void shape_terms_vel( Real coef,
                       Real mu,
                       const ElemVec& uk_loc,
                       const ElemVec& wk_loc,
                       const ElemVec& convect_loc,
                       const ElemVec& pk_loc,
                       ElemMat& elmat,
                       const CurrentFE& fe,
                       const CurrentFE& fe_p,
                      ElemMat& /*elmatP*/,
                      int iblock=0);
  //
  //! \f$coef * (  (\nabla u^k):[I\nabla\cdot  d - (\nabla d)^T] , q  )\f$ for Newton FSI
  //
void source_press( Real coef, const ElemVec& uk_loc, ElemMat& elmat,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock=0 );



  //
  //! \f$coef * (  (\nabla u^k):[I\nabla\cdot  d - (\nabla d)^T] , q  )\f$ for Newton FSI
  //
  void source_press( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock=0 );


  void source_press2( Real coef, const ElemVec& p_loc, const ElemVec& d_loc, ElemVec& elvec,
		      const CurrentFE& fe, int iblock=0 );
//@}
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


  void mass_divw( Real coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
		  int iblock, int jblock, int nb );

  //! Idem \c mass_divw \c, but with coefficient given by quadrature node
  void mass_divw(const std::vector<Real>& coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
		  int iblock, int jblock, int nb );


  void mass_gradu( Real coef, const ElemVec& u0_loc, ElemMat& elmat, const CurrentFE& fe );



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

//!@name Cholesky
//!@{
  //-------------Cholesky---------------------------------------
  /*!
    Cholesky decomposition and solution for a KNM matrix.
  */
  void choldc( KNM<Real> &a, KN<Real> &p );
  void cholsl( KNM<Real> &a, KN<Real> &p, KN<Real> &b, KN<Real> &x );
//!@}


#ifdef UNDEFINED
//these methods will become useful in a future fully implicit version of Newton FSI
void  source_mass3( Real coef, const ElemVec& un_loc, const ElemVec& uk_loc,
                    ElemMat& elmat, const CurrentFE& fe );

void source_mass2( Real coef, const ElemVec& uk_loc,
		     ElemMat& elmat, const CurrentFE& fe , double& alpha);

void source_mass22( Real coef, const ElemVec& uk_loc,
                      ElemMat& elmat, const CurrentFE& fe , double& alpha);

#endif

}
#endif
