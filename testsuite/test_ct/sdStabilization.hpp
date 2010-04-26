/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
      Date: 2005-04

 Copyright (C) 2005 INRIA

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
/**
   \file sdStabilization.hpp
   \author M.A. Fernandez
   \date 01/05/2005

   \brief This file contains a c++ class implementing the Stream-line Diffusion
          stabilization for the Navier-Stokes equations. Tested with P1/P1 and Q1/Q1

*/

#ifndef _SDSTABILIZATION_HPP_
#define _SDSTABILIZATION_HPP_


#include <life/lifecore/GetPot.hpp>

namespace LifeV
{

  template<typename MESH, typename DOF>
  class SDStabilization
  {
  public:

    typedef boost::shared_ptr<MESH> mesh_type;
    //! Constructor
    SDStabilization(
                     const mesh_type     mesh,
                     const DOF&      dof,
                     const RefFE&    refFE,
                     const QuadRule& quadRule,
                           Real	     gammaBeta = 0.,
                           Real      gammaDiv = 0.,
                           Real      viscosity = 1.);


    /*! compute SD stabilization terms and add them into matrix
     *  @param matrix matrix where the stabilization terms are added into
     *  @param state state vector for linearization of nonlinear stabilization
     */
    template<typename MATRIX, typename VECTOR>
    void apply(const Real dt, MATRIX& matrix, const VECTOR& state );

    template<typename MATRIX, typename VECTOR>
    void applyCT(const Real dt, MATRIX& matrix, const VECTOR& state );

    template <typename VECTOR, typename SOURCE >
    void apply(const Real dt, VECTOR& vector, const VECTOR& state, const SOURCE& source, const Real& time);
    void setGammaBeta (double gammaBeta) {M_gammaBeta = gammaBeta;}
    void setGammaDiv  (double gammaDiv)  {M_gammaDiv  = gammaDiv;}
  private:

    const mesh_type  M_mesh;
    const DOF&   M_dof;
    CurrentFE    M_fe;
    Real         M_viscosity;
    Real         M_gammaBeta;
    Real         M_gammaDiv;
    ElemMat      M_elMat;
    ElemVec      M_elVec;



    // methods for elementary computations
    template <typename VECTOR>
    void M_computeParameters(const Real dt, const ID iVol, const VECTOR& state, ElemVec& beta,  Real& coeffBeta, Real& coeffDiv);

    void bgradu_bgradv(const Real& coef, ElemVec& vel, ElemMat& elmat,const CurrentFE& fe,
		       int iblock, int jblock, int nb);

    void lapu_bgradv(const Real& coef, ElemVec& vel, ElemMat& elmat, const CurrentFE& fe,
		     int iblock, int jblock, int nb);

    void gradp_bgradv(const Real& coef, ElemVec& vel, ElemMat& elmat, const CurrentFE& fe);

    void lapu_gradq(const Real& coef, ElemMat& elmat, const CurrentFE& fe);

    template <typename SOURCE>
    void f_bgradv(const Real& coef, SOURCE& source, ElemVec& vel,
		  ElemVec& elvec, const CurrentFE& fe, int iblock, const Real& time);

    template<typename SOURCE>
    void f_gradq(const Real& coef, SOURCE& source, ElemVec& elvec, const CurrentFE& fe, int iblock, const Real& time);


  }; // class SDStabilization

  template<typename MESH, typename DOF>
  SDStabilization<MESH, DOF>::SDStabilization(
					       const mesh_type     mesh,
					       const DOF&      dof,
					       const RefFE&    refFE,
					       const QuadRule& quadRule,
                                               Real 	       gammaBeta,
					       Real            gammaDiv,
					       Real            viscosity):
    M_mesh( mesh ),
    M_dof( dof ),
    M_fe( refFE, getGeoMap(*mesh), quadRule ),
    M_viscosity( viscosity ),
    M_gammaBeta ( gammaBeta ),
    M_gammaDiv  ( gammaDiv ),
    M_elMat( M_fe.nbNode, nDimensions+1, nDimensions+1 ) ,
    M_elVec( M_fe.nbNode, nDimensions+1 ) {}


  template<typename MESH, typename DOF>
  template <typename MATRIX, typename VECTOR>
  void SDStabilization<MESH, DOF>::apply(const Real dt, MATRIX& matrix, const VECTOR& state )
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
    ElemVec beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 1; iVol <= M_mesh->numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstSecondDeriv( M_mesh->volumeList( iVol ) );

	// stabilization parameters computation
	chronoBeta.start();
	this->M_computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
	chronoBeta.stop();

	chronoElemComp.start();
	M_elMat.zero();

	// coeffBeta (beta \nabla u , \beta \nabla v)
	//
	this->bgradu_bgradv(coeffBeta, beta, M_elMat, M_fe, 0, 0, nDimensions);

	// coeffBeta  ( (beta \nabla u , \nabla q) + (\nabla p , \beta \nabla v) )
	//
	this->gradp_bgradv(coeffBeta, beta, M_elMat, M_fe);
	
	// coeffBeta (\nabla p , \nabla q)
	//
	stiff( coeffBeta, M_elMat, M_fe, nDimensions, nDimensions );

	// coeffBeta ( - \mu L u , \beta \nabla v )
	//
	this->lapu_bgradv(-coeffBeta*M_viscosity, beta, M_elMat, M_fe, 0, 0, nDimensions);

	// coeffBeta ( - \mu L u , \nabla q )
	//	
        this->lapu_gradq(-coeffBeta*M_viscosity, M_elMat, M_fe);

	// coeffDiv ( \div u , \div \nabla v)
	//
	stiff_div(coeffDiv, M_elMat, M_fe );

	chronoElemComp.stop();

	chronoAssembly.start();
	for ( UInt iComp = 0; iComp <= nDimensions; ++iComp )
	  for ( UInt jComp = 0; jComp <= nDimensions; ++jComp )
	    assemb_mat( matrix, M_elMat, M_fe, M_dof, iComp, jComp );
	chronoAssembly.stop();


    }// loop on elements

    std::cout << std::endl;
    std::cout << "      Updating of element   done in "
              << chronoUpdate.diff_cumul()   << "s." << std::endl;
    std::cout << "      Determination of parameters done in "
              << chronoBeta.diff_cumul()     << "s." << std::endl;
    std::cout << "      Element computations  done in "
              << chronoElemComp.diff_cumul() << "s." << std::endl;
    std::cout << "      Assembly              done in "
              << chronoAssembly.diff_cumul() << "s." << std::endl;
	

  } // apply(...)


 template<typename MESH, typename DOF>
  template <typename MATRIX, typename VECTOR>
  void SDStabilization<MESH, DOF>::applyCT(const Real dt, MATRIX& matrix, const VECTOR& state )
  {
    if ( M_gammaBeta == 0 && M_gammaDiv == 0)
      return;

    const UInt nDof = M_dof.numTotalDof();

    Chrono chronoBeta;
    Chrono chronoUpdate;
    Chrono chronoElemComp;
    Chrono chronoAssembly;

    // stabailization parameters
    Real coeffBeta, coeffDiv;

    // local velocity
    ElemVec beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 1; iVol <= M_mesh->numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstSecondDeriv( M_mesh->volumeList( iVol ) );

	// stabilization paramteres computation
	chronoBeta.start();
	this->M_computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
	chronoBeta.stop();

	chronoElemComp.start();
	M_elMat.zero();

	// coeffBeta (beta \nabla u , \beta \nabla v)
	//
	this->bgradu_bgradv(coeffBeta, beta, M_elMat, M_fe, 0, 0, nDimensions);

	// coeffBeta ( - \mu L u , \beta \nabla v )
	//
	this->lapu_bgradv(-coeffBeta*M_viscosity, beta, M_elMat, M_fe, 0, 0, nDimensions);

	// coeffDiv ( \div u , \div \nabla v)
	//
	//stiff_div(coeffDiv, M_elMat, M_fe );

	chronoElemComp.stop();

	chronoAssembly.start();
	for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
	  for ( UInt jComp = 0; jComp < nDimensions; ++jComp )
	  {
            // assemb_mat( matrix, M_elMat, M_fe, M_dof, iComp, jComp );
            assembleMatrix( matrix, M_elMat, M_fe, M_dof, iComp, jComp, 
                            iComp*nDof, jComp*nDof);
          }
	chronoAssembly.stop();


    }// loop on elements

    std::cout << std::endl;
    std::cout << "      Updating of element   done in "
              << chronoUpdate.diff_cumul()   << "s." << std::endl;
    std::cout << "      Determination of parameters done in "
              << chronoBeta.diff_cumul()     << "s." << std::endl;
    std::cout << "      Element computations  done in "
              << chronoElemComp.diff_cumul() << "s." << std::endl;
    std::cout << "      Assembly              done in "
              << chronoAssembly.diff_cumul() << "s." << std::endl;
	

  } // applyCT(...)

  template<typename MESH, typename DOF>
  template <typename VECTOR, typename SOURCE >
  void SDStabilization<MESH, DOF>::apply(const Real dt, VECTOR& vector, const VECTOR& state, const SOURCE& source, const Real& time)
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
    ElemVec beta( M_fe.nbNode, nDimensions );

    // loop on elements
    for ( UInt iVol = 1; iVol <= M_mesh->numVolumes(); iVol++ )
    {
        chronoUpdate.start();
        // update current finite elements
        M_fe.updateFirstDeriv( M_mesh->volumeList( iVol ) );
	// local mesh parameters
	chronoUpdate.stop();

	chronoBeta.start();
	this->M_computeParameters(dt, iVol, state, beta, coeffBeta, coeffDiv);
	chronoBeta.stop();

	chronoElemComp.start();
	M_elVec.zero();

	// coeffBeta ( f , \beta \nabla v)
	//
	this->f_bgradv(coeffBeta, source, beta, M_elVec, M_fe, 0, time);

	// coeffBeta  ( f , \nabla q )
	//
	this->f_gradq(coeffBeta, source, M_elVec, M_fe, nDimensions, time);
	chronoElemComp.stop();

	chronoAssembly.start();
	for ( UInt iComp = 0; iComp <= nDimensions; ++iComp )
	  assemb_vec( vector, M_elVec, M_fe, M_dof, iComp);
	chronoAssembly.stop();


    }// loop on elements

    std::cout << std::endl;
    std::cout << "      Updating of element   done in "
              << chronoUpdate.diff_cumul()   << "s." << std::endl;
    std::cout << "      Determination of parameters done in "
              << chronoBeta.diff_cumul()     << "s." << std::endl;
    std::cout << "      Element computations  done in "
              << chronoElemComp.diff_cumul() << "s." << std::endl;
    std::cout << "      Assembly              done in "
              << chronoAssembly.diff_cumul() << "s." << std::endl;
	

  } // apply(...)


  template<typename MESH, typename DOF>
  template<typename VECTOR>
  void SDStabilization<MESH, DOF>::M_computeParameters(const Real dt, const ID iVol, const VECTOR& state,
						       ElemVec& beta, Real& coeffBeta, Real& coeffDiv) {

    const UInt nDof = M_dof.numTotalDof();

    // square local mesh parameter in 1-norm
    Real hK,hK2,hK4;

    hK = M_fe.diameter();
    hK2 = hK * hK;
    hK4 = hK2 * hK2;

    // determine bmax = ||\beta||_{0,\infty,K}
    // first, get the local velocity into beta
    for ( int iNode = 0; iNode < M_fe.nbNode; ++iNode )
      {
	for ( int iCoor = 0; iCoor < M_fe.nbCoor(); ++iCoor )
	  {
	    //UInt ig = M_dof.localToGlobal( iVol, iNode+1 )-1+iCoor*nDof;
	    UInt ig = M_dof.localToGlobal( iVol, iNode+1) + iCoor*nDof; 
            beta.vec()[ iCoor*M_fe.nbNode + iNode ] = state[ig];
	  }
      }
	
    // second, calculate its max norm
    Real bmax = fabs( beta.vec()[ 0 ] );
    for ( int l = 1; l < int( M_fe.nbCoor()*M_fe.nbNode ); ++l )
      {
	if ( bmax < fabs( beta.vec()[ l ] ) )
	  bmax = fabs( beta.vec()[ l ] );
      }

    coeffBeta = M_gammaBeta  / sqrt( 4/( dt * dt)
				     + 4*bmax*bmax/hK2
				     + 16*M_viscosity*M_viscosity/hK4 );
    coeffDiv = M_gammaDiv * bmax * hK;

  }


  template<typename MESH, typename DOF>
  void SDStabilization<MESH, DOF>::gradp_bgradv(const Real& coef, ElemVec& vel,
					       ElemMat& elmat,const CurrentFE& fe) {
    ASSERT_PRE(fe.hasFirstDeriv(),
	       "advection_grad  matrix needs at least the first derivatives");

    ElemMat::matrix_type v(fe.nbCoor(),fe.nbQuadPt());
    Real s;


    // local velocity at quadrature points
    for(int ig=0;ig<fe.nbQuadPt();ig++){
      for(int icoor=0;icoor<fe.nbCoor();icoor++){
	ElemVec::vector_view velicoor=vel.block(icoor);
	v(icoor,ig)=0.;
	for(int k=0;k<fe.nbNode;k++){
	  v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
	}
      }
    }

    for (int ic=0; ic < fe.nbCoor(); ++ic) {
      ElemMat::matrix_view mat_ic3 = elmat.block(ic,fe.nbCoor());
      ElemMat::matrix_view mat_3ic = elmat.block(fe.nbCoor(),ic);
      for(int i=0;i<fe.nbNode;i++){
	for(int j=0;j<fe.nbNode;j++){
	  s = 0.0;
	  for(int ig=0;ig<fe.nbQuadPt();ig++)
	    for(int jcoor=0;jcoor<fe.nbCoor();jcoor++)
	      s += fe.phiDer(j,ic,ig)*v(jcoor,ig)*fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
	  mat_ic3(i,j) += coef*s;
	  mat_3ic(j,i) += coef*s;
	}
      }
    }
  }





  template<typename MESH, typename DOF>
  void SDStabilization<MESH, DOF>::bgradu_bgradv(const Real& coef, ElemVec& vel, ElemMat& elmat, const CurrentFE& fe,
		     int iblock, int jblock, int nb)
  {
    ASSERT_PRE(fe.hasFirstDeriv(),
	       "advection (vect) matrix needs at least the first derivatives");


    ElemMat::matrix_type mat_tmp(fe.nbNode,fe.nbNode);
    ElemMat::matrix_type v( fe.nbCoor(),fe.nbQuadPt() );
    Real s;


    // compute local vectors values
    for(int ig=0; ig<fe.nbQuadPt(); ig++)
      {
	for(int icoor=0; icoor<fe.nbCoor(); icoor++)
	  {
	    ElemVec::vector_view velicoor=vel.block(icoor);
	    v(icoor,ig)=0.;
	    for(int k=0;k<fe.nbNode;k++)
	      v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
	  }
      }

    for(int i=0;i<fe.nbNode;i++){
      for(int j=0;j<fe.nbNode;j++){
	s = 0.0;
	
	for(int ig=0;ig<fe.nbQuadPt();ig++)
	  for(int icoor=0;icoor<fe.nbCoor();icoor++)
	    for(int jcoor=0;jcoor<fe.nbCoor();jcoor++)
	      s += fe.phiDer(i,jcoor,ig)*v(jcoor,ig)*v(icoor,ig)*fe.phiDer(j,icoor,ig)*fe.weightDet(ig);
	mat_tmp(i,j) = coef*s;
      }
    }

    // copy on the components
    for(int icomp=0;icomp<nb;icomp++){
      ElemMat::matrix_view mat_icomp = elmat.block(iblock+icomp,jblock+icomp);
      for(int i=0;i<fe.nbDiag();i++){
	for(int j=0;j<fe.nbDiag();j++){
	  mat_icomp(i,j) += mat_tmp(i,j);
	}
      }
    }
  }



  template<typename MESH, typename DOF>
  void SDStabilization<MESH, DOF>::lapu_bgradv(const Real& coef, ElemVec& vel, ElemMat& elmat, const CurrentFE& fe,
		     int iblock, int jblock, int nb)
  {


    ASSERT_PRE(fe.hasFirstDeriv(),
	       "lapu_bgradv matrix needs first derivatives");

    ASSERT_PRE(fe.hasSecondDeriv(),
	       "lapu_bgradv matrix needs second derivatives");

    ElemMat::matrix_type mat_tmp(fe.nbNode,fe.nbNode);
    ElemMat::matrix_type v( fe.nbCoor(),fe.nbQuadPt() );
    Real s;


    // compute local vectors values at quadrature points
    for(int ig=0; ig<fe.nbQuadPt(); ig++)
      {
	for(int icoor=0; icoor<fe.nbCoor(); icoor++)
	  {
	    ElemVec::vector_view velicoor=vel.block(icoor);
	    v(icoor,ig)=0.;
	    for(int k=0;k<fe.nbNode;k++)
	      v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
	  }
      }


    // numerical integration
    for(int i=0;i<fe.nbNode;i++)
      {
	for(int j=0;j<fe.nbNode;j++)
	  {
	    s = 0.0;
	    for(int ig=0;ig<fe.nbQuadPt();ig++)
	      for(int icoor=0;icoor<fe.nbCoor();icoor++)
		for(int jcoor=0;jcoor<fe.nbCoor();jcoor++)
		  s += fe.phiDer2(j,icoor,icoor,ig)*v(jcoor,ig)*fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
	    mat_tmp(i,j) = coef*s;
	  }
      }

    // copy on the components
    for(int icomp=0;icomp<nb;icomp++){
      ElemMat::matrix_view mat_icomp = elmat.block(iblock+icomp,jblock+icomp);
      for(int i=0;i<fe.nbDiag();i++){
	for(int j=0;j<fe.nbDiag();j++){
	  mat_icomp(i,j) += mat_tmp(i,j);
	}
      }
    }
  }




  template<typename MESH, typename DOF>
  void SDStabilization<MESH, DOF>::lapu_gradq(const Real& coef, ElemMat& elmat,const CurrentFE& fe)
  {

    ASSERT_PRE(fe.hasFirstDeriv(),
	       "lapu_gradq matrix needs first derivatives");

    ASSERT_PRE(fe.hasSecondDeriv(),
	       "lapu_gradq matrix needs second derivatives");

    Real s;

    for (int jc=0; jc < fe.nbCoor(); ++jc) // loop on column blocks
      {
	ElemMat::matrix_view mat_view = elmat.block(fe.nbCoor(),jc);
	for(int i=0;i<fe.nbNode;++i) // local rows
	  {
	    for(int j=0;j<fe.nbNode;++j) // local columns
	      {
		s = 0.0;
		// quadrature formula
		for(int ig=0;ig<fe.nbQuadPt();++ig)
		  for(int jcoor=0;jcoor<fe.nbCoor();++jcoor) // lap
		    s += fe.phiDer2(j,jcoor,jcoor,ig)*fe.phiDer(i,jc,ig)*fe.weightDet(ig);
		mat_view(i,j) += coef*s;
	      }
	  }
      }
  }












  template<typename MESH, typename DOF>
  template<typename SOURCE>
  void SDStabilization<MESH, DOF>::f_bgradv(const Real& coef, SOURCE& source, ElemVec& vel,
					    ElemVec& elvec, const CurrentFE& fe, int iblock, const Real& time)
  {

    ASSERT_PRE(fe.hasFirstDeriv(),
	       "f_bgradv  vector needs at least the first derivatives");

    ElemMat::matrix_type v(fe.nbCoor(),fe.nbQuadPt());
    Real s;


    // local velocity at quadrature points
    for(int ig=0;ig<fe.nbQuadPt();ig++)
      {
	for(int icoor=0;icoor<fe.nbCoor();icoor++)
	  {
	    ElemVec::vector_view velicoor=vel.block(icoor);
	    v(icoor,ig)=0.;
	    for(int k=0;k<fe.nbNode;k++)
	      v(icoor,ig) += velicoor(k)*fe.phi(k,ig); // velocity on the intgt point
	  }
      }

    // local vector per block
    for (int ic=0; ic < fe.nbCoor(); ++ic)
      {
	ElemVec::vector_view vec_ic = elvec.block(ic+iblock);
	for(int i=0;i<fe.nbNode;i++)
	  {
	    s = 0.0;
	    for(int ig=0;ig<fe.nbQuadPt();ig++)
	      for(int jcoor=0;jcoor<fe.nbCoor();jcoor++)
		s += source(time,fe.quadPt(ig,0),fe.quadPt(ig,1),fe.quadPt(ig,2),ic+1)
		  *fe.phiDer(i,jcoor,ig)*v(jcoor,ig)*fe.weightDet(ig);
	    vec_ic(i) += coef*s;
	  }
      }
  }



  template<typename MESH, typename DOF>
  template<typename SOURCE>
  void SDStabilization<MESH, DOF>::f_gradq(const Real& coef, SOURCE& source, ElemVec& elvec, const CurrentFE& fe, int iblock, const Real& time)
  {

    ASSERT_PRE(fe.hasFirstDeriv(),
	       "f_gradq  vector needs at least the first derivatives");

    Real s;

    ElemVec::vector_view vec_ic = elvec.block(iblock);
    for(int i=0;i<fe.nbNode;i++)
      {
	s = 0.0;
	for(int ig=0;ig<fe.nbQuadPt();ig++)
	  for(int jcoor=0;jcoor<fe.nbCoor();jcoor++)
	    s += source(time,fe.quadPt(ig,0),fe.quadPt(ig,1),fe.quadPt(ig,2),jcoor+1)
	      *fe.phiDer(i,jcoor,ig)*fe.weightDet(ig);
	vec_ic(i) += coef*s;
      }
  }







} // namespace LifeV

#endif /* _SDSTABILIZATION_HPP_ */
