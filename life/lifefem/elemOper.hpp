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
  //               Operators for classical finite elements
  //
  //----------------------------------------------------------------------
  //coef(t,x,y,z,u)
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
  // Mass term with coefficients given for each quadrature point
  void mass( const std::vector<Real>& qpt_coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock, int jblock, int nb );

  void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
	      int iblock = 0, int jblock = 0 );
  void stiff( Real coef, Real ( *fct ) ( Real, Real, Real ), ElemMat& elmat,
	      const CurrentFE& fe, int iblock, int jblock );
  void stiff( Real coef, ElemMat& elmat, const CurrentFE& fe,
	      int iblock, int jblock, int nb );
  // Stiff term with coefficient given for each quadrature node
  void stiff( const std::vector<Real>& qpt_coef, ElemMat& elmat, const CurrentFE& fe,
	     int iblock, int jblock, int nb );  
  void stiff_curl( Real coef, ElemMat& elmat, const CurrentFE& fe,
		   int iblock, int jblock, int nb );





  // coef * ( e(u) , e(v) )
  void stiff_strain( Real coef, ElemMat& elmat, const CurrentFE& fe );

  // coef * ( div u , div v )
  void stiff_div( Real coef, ElemMat& elmat, const CurrentFE& fe );

  // coef * ( [\grad u^k]^T \grad u : \grad v  )
  void stiff_dergradbis( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe );

  // coef * ( [\grad u]^T \grad u^k + [\grad u^k]^T \grad u : \grad v  ) for Newton on St-Venant
  void stiff_dergrad( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe );

  // coef * ( \tr { [\grad u^k]^T \grad u }, \div v  ) for Newton on St-Venant
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


  //
  // Elementary operations for the interior penalty stabilization
  //

  // coef < \grad p1, \grad q2 >
  void ipstab_grad( const Real coef, ElemMat& elmat,
		    const CurrentFE& fe1, const CurrentFE& fe2,
		    const CurrentBdFE& bdfe, int iblock = 0, int jblock = 0 );

  // coef < \grad u1, \grad v2 >
  void ipstab_grad( const Real coef, ElemMat& elmat,
		    const CurrentFE& fe1, const CurrentFE& fe2,
		    const CurrentBdFE& bdfe, int iblock, int jblock, int nb );

  // coef < \div u1, \div v2 >
  void ipstab_div( const Real coef, ElemMat& elmat,
		   const CurrentFE& fe1, const CurrentFE& fe2,
		   const CurrentBdFE& bdfe, int iblock = 0, int jblock = 0 );
  // coef < \beta1 . \grad u1, \beta2 . \grad v2 >
  void ipstab_bgrad( const Real coef, ElemMat& elmat,
		     const CurrentFE& fe1, const CurrentFE& fe2,
		     const ElemVec& beta, const CurrentBdFE& bdfe,
		     int iblock, int jblock, int nb );
  // coef < |\beta . n|^2 / |\beta| \grad p1, \grad q2 >
  void ipstab_bagrad( const Real coef, ElemMat& elmat,
		      const CurrentFE& fe1, const CurrentFE& fe2,
		      const ElemVec& beta, const CurrentBdFE& bdfe,
		      int iblock = 0, int jblock = 0 );

// coef < |\beta.n| \grad p1, \grad q2 >
// p1 lives in fe1
// q2 lives in fe2
// beta lives in fe3

void ipstab_bagrad( const Real           coef,
                    ElemMat&             elmat,
                    const CurrentFE&     fe1,
                    const CurrentFE&     fe2,
                    const CurrentFE&     fe3,
                    const ElemVec&       beta,
                    const CurrentBdFE&   bdfe,
                    int iblock = 0, int jblock = 0 );


///////////////////////////////////////


  ////////////////////////////////////////
  // Convective term with a local vector coefficient (useful for Navier-Stokes problem)
  void grad( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
	     const CurrentFE& fe1, const CurrentFE& fe2,
	     int iblock, int jblock );

  // Convective term with a local vector coefficient (useful for Navier-Stokes problem+adv-diff)
  void grad( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
	     const CurrentFE& fe1, const CurrentFE& fe2, const CurrentFE& fe3,
	     int iblock = 0, int jblock = 0 );

  // Convective term with a local vector coefficient for Navier-Stokes problem in Skew-Symmetric form
  void grad_ss( const int icoor, const ElemVec& vec_loc, ElemMat& elmat,
		const CurrentFE& fe1, const CurrentFE& fe2,
		int iblock = 0, int jblock = 0 );

  // StreamLine Diffusion
  void stiff_sd( Real coef, const ElemVec& vec_loc, ElemMat& elmat, const CurrentFE& fe,
		 const CurrentFE& fe2, int iblock = 0, int jblock = 0, int nb = 1 );

  /////////////////////////////////////////////
  //
  // source  \int fct phi_i
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
  // source  \int fct(t) phi_i
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
  //
  // coef * ( \grad (-w^k):[I\div d - (\grad d)^T] u^k + convect^T[I\div d - (\grad d)^T] (\grad u^k)^T , v  ) for Newton FSI
  //
  // Remark: convect = u^n-w^k
  //
  void source_mass1( Real coef,
                     const ElemVec& uk_loc,
                     const ElemVec& wk_loc,
                     const ElemVec& convect_loc,
                     const ElemVec& d_loc,
                     ElemVec& elvec,
                     const CurrentFE& fe );

  void source_mass1( Real coef,
                     const ElemVec& uk_loc,
                     const ElemVec& wk_loc,
                     const ElemVec& convect_loc,
                     ElemMat& elmat,
                     const CurrentFE& fe );

  void source_mass11( Real coef,
                     const ElemVec& uk_loc,
                     const ElemVec& wk_loc,
                     const ElemVec& convect_loc,
                     ElemMat& elmat,
                     const CurrentFE& fe );


  //
  // coef * ( \grad u^k dw, v  ) for Newton FSI
  //
  //
  void source_mass2( Real coef, const ElemVec& uk_loc, const ElemVec& dw_loc,
		     ElemVec& elvec, const CurrentFE& fe );

  void source_mass2( Real coef, const ElemVec& uk_loc,
		     ElemMat& elmat, const CurrentFE& fe , double& alpha);
  void source_mass22( Real coef, const ElemVec& uk_loc,
                      ElemMat& elmat, const CurrentFE& fe , double& alpha);


  void  source_mass3( Real coef, const ElemVec& un_loc, const ElemVec& uk_loc, const ElemVec& d_loc,
		      ElemVec& elvec, const CurrentFE& fe );

void  source_mass3( Real coef, const ElemVec& un_loc, const ElemVec& uk_loc,
                    ElemMat& elmat, const CurrentFE& fe );





  //
  // coef * ( [-p^k I + 2*mu e(u^k)] [I\div d - (\grad d)^T] , \grad v  ) for Newton FSI
  //
  void source_stress( Real coef, Real mu, const ElemVec& uk_loc, const ElemVec& pk_loc,
		      const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u,
		      const CurrentFE& fe_p );

  void source_stress( Real coef, Real mu, const ElemVec& uk_loc, const ElemVec& pk_loc,
		      ElemMat& elmat, const CurrentFE& fe_u,
		      const CurrentFE& fe_p );

  //
  // + \mu ( \grad u^k \grad d + [\grad d]^T[\grad u^k]^T : \grad v )
  //
  void source_stress2( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec, const CurrentFE& fe_u );

  void source_stress2( Real coef, const ElemVec& uk_loc, ElemMat& elmat, const CurrentFE& fe_u );


  //
  // coef * (  (\grad u^k):[I\div d - (\grad d)^T] , q  ) for Newton FSI
  //
  void source_press( Real coef, const ElemVec& uk_loc, const ElemVec& d_loc, ElemVec& elvec,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock=0 );

  void source_press( Real coef, const ElemVec& uk_loc, ElemMat& elmat,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock=0 );

  void source_presss( Real coef, const ElemVec& uk_loc, ElemMat& elmat,
                   const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock=0 );

  void source_press2( Real coef, const ElemVec& p_loc, const ElemVec& d_loc, ElemVec& elvec,
		      const CurrentFE& fe, int iblock=0 );
  //----------------------------------------------------------------------
  //
  //                  Operators for H(div) finite elements
  //
  //----------------------------------------------------------------------
  //------------Transpose of Elementary divergence matrix--------------
  void grad_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
		  const CurrentFE& fe_p, int iblock, int jblock );

  //----------------Elementary divergence matrix-----------------------
  void div_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe_u,
		 const CurrentFE& fe_p, int iblock, int jblock );

  //----------------Elementary tp v dot n matrix-----------------------
  void TP_VdotN_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe, const RefHybridFE& vdotnfe, int iblock, int jblock );

  //----------------Elementary tp^2  matrix-----------------------
  void TP_TP_Hdiv( Real coef, ElemMat& elmat, const RefHybridFE& tpfe, int iblock, int jblock );

  //-------------Mass matrix---------------------------------------
  /*!
    Weighted Mass matrix with a permeability tensor which is a constant scalar matrix
    (i.e. K^{-1} = coef * Id, coef being the inverse of the permeability).
  */
  void mass_Hdiv( Real coef, ElemMat& elmat, const CurrentHdivFE& fe,
		  int iblock = 0, int jblock = 0 );


  void mass_divw( Real coef, const ElemVec& w_loc, ElemMat& elmat, const CurrentFE& fe,
		  int iblock, int jblock, int nb );

  // Idem mass_divw, but with coefficient given by quadrature node
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


  //-------------Cholesky---------------------------------------
  /*!
    Cholesky decomposition and solution for a KNM matrix.
  */
  void choldc( KNM<Real> &a, KN<Real> &p );
  void cholsl( KNM<Real> &a, KN<Real> &p, KN<Real> &b, KN<Real> &x );
}
#endif
