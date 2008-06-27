/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file oneDModelSolver.cpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file implements a Taylor-Galerkin solver for 1D model.
*/

#include <life/lifesolver/oneDModelSolver.hpp>


namespace LifeV
{
  OneDModelSolver::OneDModelSolver(const GetPot& data_file,
				   const OneDNonLinModelParam& onedparam):
    // const LinearSimpleParam& onedparam):
    OneDModelHandler(data_file),
    _M_oneDParam(onedparam),
    _M_fluxFun(_M_oneDParam),
    _M_sourceFun(_M_oneDParam),
    //! id of left and right bc nodes
    _M_leftNodeId( 0 ),
    _M_leftInternalNodeId( _M_leftNodeId + 1 ),
    _M_rightNodeId( _M_dimDof - 1 ),
    _M_rightInternalNodeId( _M_rightNodeId  - 1 ),
    //! boundary edges
    _M_leftEdge( _M_mesh.edgeList( 1 ) ),
    _M_rightEdge( _M_mesh.edgeList(  _M_nb_elem ) ),
    // elementary matrices
    _M_elmatMass (_M_fe.nbNode,1,1),
    _M_elmatStiff(_M_fe.nbNode,1,1),
    _M_elmatGrad (_M_fe.nbNode,1,1),
    _M_elmatDiv  (_M_fe.nbNode,1,1),
    // vectorial unknowns and rhs
    _M_U1_thistime(_M_dimDof),
    _M_U2_thistime(_M_dimDof),
    _M_W1_thistime(_M_dimDof),
    _M_W2_thistime(_M_dimDof),
    _M_U1_timestep(_M_dimDof),
    _M_U2_timestep(_M_dimDof),
    _M_rhs1(_M_dimDof),
    _M_rhs2(_M_dimDof),
    // vectors and matrices of the non-linear function
    _M_Flux1(_M_dimDof),
    _M_Flux2(_M_dimDof),
    _M_diffFlux11(_M_nb_elem),
    _M_diffFlux12(_M_nb_elem),
    _M_diffFlux21(_M_nb_elem),
    _M_diffFlux22(_M_nb_elem),
    _M_Source1(_M_dimDof),
    _M_Source2(_M_dimDof),
    _M_diffSrc11(_M_nb_elem),
    _M_diffSrc12(_M_nb_elem),
    _M_diffSrc21(_M_nb_elem),
    _M_diffSrc22(_M_nb_elem),
    // mass matrix (to be inverted)
    _M_massMatrix(_M_dimDof),
    _M_factorMassMatrix(_M_dimDof),
    _M_tridiagSlv(_M_dimDof),
    // matrices used to build the rhs
    _M_massMatrixDiffSrc11(_M_dimDof),
    _M_massMatrixDiffSrc12(_M_dimDof),
    _M_massMatrixDiffSrc21(_M_dimDof),
    _M_massMatrixDiffSrc22(_M_dimDof),
    _M_stiffMatrixDiffFlux11(_M_dimDof),
    _M_stiffMatrixDiffFlux12(_M_dimDof),
    _M_stiffMatrixDiffFlux21(_M_dimDof),
    _M_stiffMatrixDiffFlux22(_M_dimDof),
    _M_gradMatrix(_M_dimDof),
    _M_gradMatrixDiffFlux11(_M_dimDof),
    _M_gradMatrixDiffFlux12(_M_dimDof),
    _M_gradMatrixDiffFlux21(_M_dimDof),
    _M_gradMatrixDiffFlux22(_M_dimDof),
    _M_divMatrixDiffSrc11(_M_dimDof),
    _M_divMatrixDiffSrc12(_M_dimDof),
    _M_divMatrixDiffSrc21(_M_dimDof),
    _M_divMatrixDiffSrc22(_M_dimDof),
    // Handle boundary conditions
    _M_bcH(_M_U1_thistime, _M_U2_thistime, _M_W1_thistime, _M_W2_thistime, _M_fluxFun),
	_M_CFL(data_file("miscellaneous/showCFL",0)),
	_M_UW(data_file("miscellaneous/alternate_solver",0))
  {
	_M_oneDstring2initializeVarMap["P"] = OneDInitPressure;
	_M_oneDstring2initializeVarMap["A"] = OneDInitArea;
	_M_oneDstring2initializeVarMap["Q"] = OneDInitFlux;
	_M_oneDstring2initializeVarMap["W1"] = OneDInitReimann1;
	_M_oneDstring2initializeVarMap["W2"] = OneDInitReimann2;

    Debug( 6030 ) << "[OneDModelSolver::OneDModelSolver] O-  Nb of unknowns: " << _M_dimDof     << "\n";
    Debug( 6030 ) << "[OneDModelSolver::OneDModelSolver] O-  Computing the constant matrices... \n";

    Chrono chrono;
    chrono.start();

    //! Matrices initialization
    _M_massMatrix.zero();
    // _M_factorMassMatrix.zero();
    _M_factorMassMatrix.zero();
    _M_gradMatrix.zero();

    //-------------------------------------------
    //! update first the constant matrices (cst w.r. to time iter)

    //! set the coeff to 1.
    _M_coeffMass = 1.;
    _M_coeffGrad = 1.;

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

      //! set the elementary matrices to 0.
      _M_elmatMass.zero();
      _M_elmatGrad.zero();

      //! update the current element
      _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
      //std::cout << _M_fe.currentId() << std::endl;

      //! update the mass and grad matrices
      mass( _M_coeffMass, _M_elmatMass, _M_fe,0, 0 );
      grad( 0 , - _M_coeffGrad, _M_elmatGrad, _M_fe, _M_fe, 0, 0 );

      //! assemble the mass and grad matrices
      assemb_mat( _M_massMatrix, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );
      assemb_mat( _M_gradMatrix, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

    } //! end loop on elements

    _M_factorMassMatrix = _M_massMatrix;

    //! Dirichlet boundary conditions set in the mass matrix
    _updateBCDirichletMatrix( _M_factorMassMatrix );

    //! cholesky or lapack lu factorization of the mass matrix
    _M_tridiagSlv.Factor( _M_factorMassMatrix );

    for( UInt i=0; i<_M_dimDof; i++)
      {
	_M_rhs1(i)=0.;
	_M_rhs2(i)=0.;
      }

    chrono.stop();

    Debug( 6030 ) << "[OneDModelSolver::OneDModelSolver] \tdone in " << chrono.diff() << " s.\n";

  }

  /*! Update the coefficients
    (from the flux, source functions and their derivatives)
  */
  void OneDModelSolver::
  _updateMatrixCoefficients(const UInt& ii, const UInt& jj , const UInt& iedge)
  {
    Real dFluxdUelem(0), dSrcdUelem(0);

    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    if( ii == 1 && jj == 1 ) {
      dFluxdUelem = _M_diffFlux11( iedge - 1 ); //! iedge starts from 1...
      dSrcdUelem  = _M_diffSrc11( iedge - 1 );
    }
    if( ii == 1 && jj == 2 ) {
      dFluxdUelem = _M_diffFlux12( iedge - 1 );
      dSrcdUelem  = _M_diffSrc12( iedge - 1 );
    }
    if( ii == 2 && jj == 1 ) {
      dFluxdUelem = _M_diffFlux21( iedge - 1 );
      dSrcdUelem  = _M_diffSrc21( iedge - 1 );
    }
    if( ii == 2 && jj == 2 ) {
      dFluxdUelem = _M_diffFlux22( iedge - 1 );
      dSrcdUelem  = _M_diffSrc22( iedge - 1 );
    }

    _M_coeffGrad  = dFluxdUelem; //! term gradDiffFlux(U) [* S(U)]
    _M_coeffDiv   = dSrcdUelem;  //! term  divDiffSrc(U) [* F(U)]
    _M_coeffStiff = dFluxdUelem; //! term stiffDiffFlux(U) [* F(U)]
    _M_coeffMass  = dSrcdUelem;  //! term  massDiffSrc(U) [* S(U)]
  }

  //! Update the element matrices with the updated
  //! current element and updated coefficients
  void OneDModelSolver::_updateElemMatrices()
  {
    //! set the elementary matrices to 0.
    _M_elmatMass.zero();
    _M_elmatStiff.zero();
    _M_elmatGrad.zero();
    _M_elmatDiv.zero();

    //! update the mass matrix
    mass( _M_coeffMass, _M_elmatMass, _M_fe,0, 0 );
    //  std::cout << "Elem Mass matrix :" << std::endl;
    //  _M_elmatMass.showMe( std::cout );

    //! update the stiffness matrix
    stiff( _M_coeffStiff, _M_elmatStiff, _M_fe,0 ,0 );
    // std::cout << "Elem Stiff matrix :" << std::endl;
    // _M_elmatStiff.showMe( std::cout );

    /*! update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - _M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0 , - _M_coeffGrad, _M_elmatGrad, _M_fe, _M_fe, 0, 0 );
    //  std::cout << "Elem Grad matrix :" << std::endl;
    //  _M_elmatGrad.showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this _M_elmatDiv is not necessary
      as it is the transpose of the _M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0 , - _M_coeffDiv, _M_elmatDiv, _M_fe, _M_fe, 0, 0 );
    //  std::cout << "Elem Div matrix :" << std::endl;
    //  _M_elmatDiv.showMe( std::cout );
  }


  //! assemble the matrices
  int OneDModelSolver::
  _assemble_matrices(const UInt& ii, const UInt& jj )
  {
    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    if( ii == 1 && jj == 1 ) {
      //! assemble the mass matrix
      assemb_mat( _M_massMatrixDiffSrc11, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the stiffness matrix
      assemb_mat( _M_stiffMatrixDiffFlux11, _M_elmatStiff, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the gradient matrix
      assemb_mat( _M_gradMatrixDiffFlux11, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the divergence matrix
      assemb_mat( _M_divMatrixDiffSrc11, _M_elmatDiv, _M_fe, _M_dof1D , 0, 0 );

      return 0;
    }
    if( ii == 1 && jj == 2 ) {
      //! assemble the mass matrix
      assemb_mat( _M_massMatrixDiffSrc12, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the stiffness matrix
      assemb_mat( _M_stiffMatrixDiffFlux12, _M_elmatStiff, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the gradient matrix
      assemb_mat( _M_gradMatrixDiffFlux12, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the divergence matrix
      assemb_mat( _M_divMatrixDiffSrc12, _M_elmatDiv, _M_fe, _M_dof1D , 0, 0 );

      return 0;
    }
    if( ii == 2 && jj == 1 ) {
      //! assemble the mass matrix
      assemb_mat( _M_massMatrixDiffSrc21, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the stiffness matrix
      assemb_mat( _M_stiffMatrixDiffFlux21, _M_elmatStiff, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the gradient matrix
      assemb_mat( _M_gradMatrixDiffFlux21, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the divergence matrix
      assemb_mat( _M_divMatrixDiffSrc21, _M_elmatDiv, _M_fe, _M_dof1D , 0, 0 );

      return 0;
    }
    if( ii == 2 && jj == 2 ) {
      //! assemble the mass matrix
      assemb_mat( _M_massMatrixDiffSrc22, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the stiffness matrix
      assemb_mat( _M_stiffMatrixDiffFlux22, _M_elmatStiff, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the gradient matrix
      assemb_mat( _M_gradMatrixDiffFlux22, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

      //! assemble the divergence matrix
      assemb_mat( _M_divMatrixDiffSrc22, _M_elmatDiv, _M_fe, _M_dof1D , 0, 0 );

      return 0;
    }

    ERROR_MSG("Invalid values for the _assemble_matrices method.");
    return 1;
  }

  /*! update the matrices
    _M_massMatrixDiffSrcij, _M_stiffMatrixDiffFluxij
    _M_gradMatrixDiffFluxij, and _M_divMatrixDiffSrcij (i,j=1,2)

    from the values of diffFlux(Un) and diffSrc(Un)
    that are computed with _updateMatrixCoefficients.
  */
  void OneDModelSolver::_updateMatrices()
  {
    //--------------------------------------------------------
    // Chrono chrono;
    // chrono.start();
    // std::cout << "o-loop over the matrices INIT... ";

    //! Matrices initialization
    _M_massMatrixDiffSrc11.zero();
    _M_massMatrixDiffSrc12.zero();
    _M_massMatrixDiffSrc21.zero();
    _M_massMatrixDiffSrc22.zero();

    _M_stiffMatrixDiffFlux11.zero();
    _M_stiffMatrixDiffFlux12.zero();
    _M_stiffMatrixDiffFlux21.zero();
    _M_stiffMatrixDiffFlux22.zero();

    _M_gradMatrixDiffFlux11.zero();
    _M_gradMatrixDiffFlux12.zero();
    _M_gradMatrixDiffFlux21.zero();
    _M_gradMatrixDiffFlux22.zero();

    _M_divMatrixDiffSrc11.zero();
    _M_divMatrixDiffSrc12.zero();
    _M_divMatrixDiffSrc21.zero();
    _M_divMatrixDiffSrc22.zero();

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices... ";
    */

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

      //! update the current element
      _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
      //  std::cout << _M_fe.currentId() << std::endl;

      for(UInt ii = 1; ii <= 2; ii ++) {
	for(UInt jj = 1; jj <= 2; jj ++) {

	  //! update the _M_coeff*
	  _updateMatrixCoefficients( ii , jj, iedge);

	  //! update the _M_elmat*
	  _updateElemMatrices();

	  //! assemble the global matrices
	  _assemble_matrices( ii, jj );
	}
      }
    } //! end loop on elements

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices BC DIR... ";
    */

    /*
    //! useless ??????
    //! taking into account the dirichlet bc
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc11 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc12 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc21 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc22 );

    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux11 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux12 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux21 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux22 );

    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux11 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux12 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux21 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux22 );

    _updateBCDirichletMatrix( _M_divMatrixDiffSrc11 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc12 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc21 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc22 );
    // chrono.stop();
    // std::cout << "done in " << chrono.diff() << " s." << std::endl;
    */
  }


  /*! modify the matrix to take into account
    the Dirichlet boundary conditions
    (works for P1Seg and canonic numbering!)
  */
  void OneDModelSolver::
  _updateBCDirichletMatrix( TriDiagMatrix<Real>& mat )
  {
    UInt firstDof = 0;
    UInt lastDof  = mat.OrderMatrix()-1;

    /* 
    //! unsymmetric treatment (LU must be used!)
    //! modify the first row
    mat.Diag()( firstDof )   = 1.;
    mat.UpDiag()( firstDof ) = 0.;

    //! modify the last row
    mat.Diag()( lastDof )      = 1.;
    mat.LowDiag()( lastDof-1 ) = 0.;
    */

    //! symmetric treatment (cholesky can be used)
    //! modify the first row
    mat.Diag()( firstDof )    = 1.;
    mat.UpDiag()( firstDof )  = 0.;
    mat.LowDiag()( firstDof ) = 0.; //!and second row

    //! modify the last row
    mat.Diag()( lastDof )      = 1.;
    mat.UpDiag()( lastDof-1 )  = 0.; //!and penultimate row
    mat.LowDiag()( lastDof-1 ) = 0.;
  }

  /*! modify the vector to take into account
    the Dirichlet boundary conditions
    (works for P1Seg and canonic numbering!)

    \param val_left  : Dirichlet value inserted to the left
    \param val_right : Dirichlet value inserted to the right
  */
  void  OneDModelSolver::
  _updateBCDirichletVector()
  {
    UInt firstDof = 0;
    UInt lastDof  = _M_rhs1.size()-1;

    /* 
    //! unsymmetric treatment (LU must be used!)
    //! first row modified
    _M_rhs1( firstDof ) = _M_bcDirLeft.first;
    _M_rhs2( firstDof ) = _M_bcDirLeft.second;

    //! last row modified
    _M_rhs1( lastDof ) = _M_bcDirRight.first;
    _M_rhs2( lastDof ) = _M_bcDirRight.second;
    */

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _M_rhs1( firstDof ) = _M_bcDirLeft.first;
    _M_rhs2( firstDof ) = _M_bcDirLeft.second;
    //! second row modified (for symmetry)
    _M_rhs1( firstDof + 1 ) += - _M_massMatrix.LowDiag()( firstDof ) * _M_bcDirLeft.first;
    _M_rhs2( firstDof + 1 ) += - _M_massMatrix.LowDiag()( firstDof ) * _M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _M_rhs1( lastDof ) = _M_bcDirRight.first;
    _M_rhs2( lastDof ) = _M_bcDirRight.second;
    //! penultimate row modified (for symmetry)
    _M_rhs1( lastDof - 1 ) += - _M_massMatrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.first;
    _M_rhs2( lastDof - 1 ) += - _M_massMatrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.second;

  }


  //! compute the _M_bcDirLeft and _M_bcDirRight and set them to the new values
  void OneDModelSolver::_computeBC( const Real& time_val )
  {
	_M_bcH.applyBC(time_val, _M_bcDirLeft, _M_bcDirRight );
  }


  //! get the flux function
  const NonLinearFluxFun1D& OneDModelSolver::FluxFun() const
  {
    return _M_fluxFun;
  }


  //! get the source function
  const NonLinearSourceFun1D& OneDModelSolver::SourceFun() const
  {
    return _M_sourceFun;
  }


  //! get the left edge
  Edge1D OneDModelSolver::LeftEdge() const
  {
    return _M_leftEdge;
  }


  //! get the right edge
  Edge1D OneDModelSolver::RightEdge() const
  {
    return _M_rightEdge;
  }


  //! get the left node
  UInt OneDModelSolver::LeftNodeId() const
  {
    return _M_leftNodeId;
  }


  //! get the left internal node (neighboring node)
  UInt OneDModelSolver::LeftInternalNodeId() const
  {
    return _M_leftInternalNodeId;
  }


  //! get the right node
  UInt OneDModelSolver::RightNodeId() const
  {
    return _M_rightNodeId;
  }


  //! get the right internal node (neighboring node)
  UInt OneDModelSolver::RightInternalNodeId() const
  {
    return _M_rightInternalNodeId;
  }


  //! get the Dirichlet boundary conditions (left)
  OneDModelSolver::Vec2D OneDModelSolver::BCValuesLeft() const
  {
    return Vec2D( _M_U1_thistime( LeftNodeId() ),
                  _M_U2_thistime( LeftNodeId() ) );
  }


  //! get the value at neighboring node (left)
  OneDModelSolver::Vec2D OneDModelSolver::BCValuesInternalLeft() const
  {
    return Vec2D( _M_U1_thistime( LeftInternalNodeId() ),
                  _M_U2_thistime( LeftInternalNodeId() ) );
  }


  //! get the Dirichlet boundary conditions (right)
  OneDModelSolver::Vec2D OneDModelSolver::BCValuesRight() const 
  {
    return Vec2D( _M_U1_thistime( RightNodeId() ),
                  _M_U2_thistime( RightNodeId() ) );
  }


  //! get the value at neighboring node (right)
  OneDModelSolver::Vec2D OneDModelSolver::BCValuesInternalRight() const 
  {
    return Vec2D( _M_U1_thistime( RightInternalNodeId() ),
                  _M_U2_thistime( RightInternalNodeId() ) );
  }


  //! set the Dirichlet boundary conditions (right)
  void OneDModelSolver::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
  {
    _M_bcDirRight.first  = bcR1;
    _M_bcDirRight.second = bcR2;
  }


  //! set the Dirichlet boundary conditions (left)
  void OneDModelSolver::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
  {
    _M_bcDirLeft.first  = bcL1;
    _M_bcDirLeft.second = bcL2;
  }


  void OneDModelSolver::setBCLeft_internalnode()
  {
    _M_bcH.setBCLeft_internalnode();
  }    


  void OneDModelSolver::setBCRight_internalnode()
  {
    _M_bcH.setBCRight_internalnode();
  }


  //! simple cfl computation (correct for constant mesh)
  void OneDModelSolver::CheckCFL() const
  {
    Real CFL = 0.;

    //! length of the first edge (arbitrary as they are all supposed equal).
    Real deltaX,
      deltaX_min = _M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x();
    Real Ainode, Qinode;

    Real lambda1_max = 0.;
    Real lambda2_max = 0.;
    Real eigval1, eigval2;
    Real tmp11, tmp12, tmp21, tmp22;

    for ( UInt inode=0; inode < _M_dimDof ; inode++ ) {

      Ainode = _M_U1_thistime( inode );
      Qinode = _M_U2_thistime( inode );

      //! compute the eigenvalues at node
      _M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
					       eigval1, eigval2,
					       tmp11, tmp12,
					       tmp21, tmp22,
					       inode );


      lambda1_max = std::max<Real>( std::fabs(eigval1), lambda1_max );
      lambda2_max = std::max<Real>( std::fabs(eigval2), lambda2_max );
    }

    for ( UInt inode=1; inode < _M_dimDof ; inode++ ) {
      
      deltaX = _M_mesh.edgeList( inode ).pt2().x() - _M_mesh.edgeList( inode ).pt1().x();
    
      deltaX_min = std::min<Real>( std::fabs(deltaX), deltaX_min );
    }

    CFL = _M_time_step / deltaX_min * std::max<Real>( lambda1_max , lambda2_max );

    if ( _M_CFL ) 
    	std::cout << "CFL = " << CFL << std::endl;
	
    /*
      std::cout << "Old CFL = " << _M_time_step / 
      (_M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x())
      * std::max<Real>( lambda1_max , lambda2_max )
      << std::endl;
    */	
    //    }
    if( CFL > 0.5774 )
      std::cout << "\n[OneDModelSolver::CheckCFL] CFL not respected in " << _M_post_file
		<< ": CFL = " << CFL << std::endl;
    //      ASSERT( CFL < 0.5774 , "CFL not respected" );
  }


  /*! Axpy product for 2D vectors (pairs)
    Axpy(alpha, x, beta, y) -> y = a*A*x + beta*y

    A is given by two pairs corresponding to the 2 lines.
    A = [line1;
    line2 ]
  */
  void OneDModelSolver::Axpy(const Vec2D& line1, const Vec2D& line2,
			     const Real& alpha,  const Vec2D& x,
			     const Real& beta,   Vec2D& y) const
  {
    y.first  = alpha * ( line1.first * x.first + line1.second * x.second ) + beta * y.first;
    y.second = alpha * ( line2.first * x.first + line2.second * x.second ) + beta * y.second;
  }


  //! update the P1 flux vector from U: _M_Fluxi = F_h(Un) i=1,2
  //! BEWARE: works only for P1Seg elements
  void OneDModelSolver::_updateFlux()
  {
    Real Aii, Qii;

    for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
      Aii = _M_U1_thistime( ii );
      Qii = _M_U2_thistime( ii );
      _M_Flux1( ii ) = _M_fluxFun( Aii, Qii, 1, ii );
      _M_Flux2( ii ) = _M_fluxFun( Aii, Qii, 2, ii );
    }
  }


  /*! call _updateFlux and update the P0 derivative of flux vector from U:
    _M_diffFluxij = dF_h/dU(Un) i,j=1,2

    _M_diffFluxij(elem) = 1/2 [ dF/dU(U(node1(elem))) + dF/dU(U(node2(elem))) ]

    (mean value of the two extremal values of dF/dU)

    BEWARE: works only for P1Seg elements
  */
  void OneDModelSolver::_updateFluxDer()
  {
    //! first update the Flux vector
    _updateFlux();

    //! then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < _M_nb_elem ; ielem++ ) {
      //! for P1Seg and appropriate mesh only!
      ii = ielem;        //! left node of current element
      iip1 = ielem + 1;  //! right node of current element

      Aii = _M_U1_thistime( ielem );
      Qii = _M_U2_thistime( ielem );
      Aiip1 = _M_U1_thistime( ielem + 1 );
      Qiip1 = _M_U2_thistime( ielem + 1 );

      //! diffFlux11
      tmp =  _M_fluxFun.diff(   Aii,   Qii, 1, 1, ii );
      tmp += _M_fluxFun.diff( Aiip1, Qiip1, 1, 1, iip1 );
      _M_diffFlux11( ielem ) = 0.5 * tmp;
      //! diffFlux12
      tmp =  _M_fluxFun.diff(   Aii,   Qii, 1, 2, ii );
      tmp += _M_fluxFun.diff( Aiip1, Qiip1, 1, 2, iip1 );
      _M_diffFlux12( ielem ) = 0.5 * tmp;
      //! diffFlux21
      tmp =  _M_fluxFun.diff(   Aii,   Qii, 2, 1, ii );
      tmp += _M_fluxFun.diff( Aiip1, Qiip1, 2, 1, iip1 );
      _M_diffFlux21( ielem ) = 0.5 * tmp;
      //! diffFlux22
      tmp =  _M_fluxFun.diff(   Aii,   Qii, 2, 2, ii );
      tmp += _M_fluxFun.diff( Aiip1, Qiip1, 2, 2, iip1 );
      _M_diffFlux22( ielem ) = 0.5 * tmp;

    }
  }

  //! update the P1 source vector from U: _M_Sourcei = S_h(Un) i=1,2
  //! BEWARE: works only for P1Seg elements
  void OneDModelSolver::_updateSource()
  {
    Real Aii, Qii;

    for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
      Aii = _M_U1_thistime( ii );
      Qii = _M_U2_thistime( ii );
      _M_Source1( ii ) = _M_sourceFun( Aii, Qii, 1, ii );
      _M_Source2( ii ) = _M_sourceFun( Aii, Qii, 2, ii );
    }
  }

  /*! call _updateSource and update the P0 derivative of source vector from U:
    _M_diffSrcij = dS_h/dU(Un) i,j=1,2

    _M_diffSrcij(elem) = 1/2 [ dS/dU(U(node1(elem))) + dS/dU(U(node2(elem))) ]

    (mean value of the two extremal values of dS/dU)

    BEWARE: works only for P1Seg elements
  */
  void OneDModelSolver::_updateSourceDer()
  {
    //! first update the Source vector
    _updateSource();

    //! then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < _M_nb_elem ; ielem++ ) {
      //! for P1Seg and appropriate mesh only!
      ii = ielem;        //! left node of current element
      iip1 = ielem + 1;  //! right node of current element

      Aii = _M_U1_thistime( ielem );
      Qii = _M_U2_thistime( ielem );
      Aiip1 = _M_U1_thistime( ielem + 1 );
      Qiip1 = _M_U2_thistime( ielem + 1 );

      //! diffSrc11
      tmp =  _M_sourceFun.diff(   Aii,   Qii, 1, 1, ii );
      tmp += _M_sourceFun.diff( Aiip1, Qiip1, 1, 1, iip1 );
      _M_diffSrc11( ielem ) = 0.5 * tmp;
      //! diffSrc12
      tmp =  _M_sourceFun.diff(   Aii,   Qii, 1, 2, ii );
      tmp += _M_sourceFun.diff( Aiip1, Qiip1, 1, 2, iip1 );
      _M_diffSrc12( ielem ) = 0.5 * tmp;
      //! diffSrc21
      tmp =  _M_sourceFun.diff(   Aii,   Qii, 2, 1, ii );
      tmp += _M_sourceFun.diff( Aiip1, Qiip1, 2, 1, iip1 );
      _M_diffSrc21( ielem ) = 0.5 * tmp;
      //! diffSrc22
      tmp =  _M_sourceFun.diff(   Aii,   Qii, 2, 2, ii );
      tmp += _M_sourceFun.diff( Aiip1, Qiip1, 2, 2, iip1 );
      _M_diffSrc22( ielem ) = 0.5 * tmp;

    }
  }


  //! Initialize from Reimann invariants
  //! Initialize with constant initial conditions
  void OneDModelSolver::initialize(const Real& u10, const Real& u20, const std::string& var )
  {
  	if( var == "physical")
  	{
    	_M_U1_thistime = ScalarVector( _M_U1_thistime.size(), u10 );
    	_M_U2_thistime = ScalarVector( _M_U2_thistime.size(), u20 );

      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );
	}
	else if( var == "reimann" )
	{    
		_M_W1_thistime = ScalarVector( _M_U1_thistime.size(), u10 );
		_M_W2_thistime = ScalarVector( _M_U2_thistime.size(), u20 );
	
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.U_from_W( _M_U1_thistime[ielem], _M_U2_thistime[ielem],
				_M_W1_thistime[ielem], _M_W2_thistime[ielem], ielem );
    }
	else
	{
		std::cout << "[OneDModelSolver::initialize] trying to initialize " << var << " variables!" << std::endl;
		abort();
	}
	
    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step ); 

    //! create matlab scripts  
    create_movie_file();
	
    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

  }

  //! Initialize when initial conditions concentration
  void OneDModelSolver::initialize(const Function& /*c0*/, Real /*t0*/, Real /*dt*/)
  {
    ERROR_MSG("Not yet implemented");
  }

  // ! Initialize when initial values for the concentration are read from file
  void OneDModelSolver::initialize(const std::string & vname)
  {
    ERROR_MSG("Not yet implemented");
    /* 
    //! create matlab scripts  
    create_movie_file();
    */
    std::fstream Resfile(vname.c_str(),std::ios::in | std::ios::binary);
    if (Resfile.fail()) {
      std::cerr<<" Error in initialize: File not found or locked"<<std::endl;
      abort();
    }
    Resfile.read((char*)&_M_U1_thistime(1),_M_U1_thistime.size()*sizeof(Real));
    Resfile.close();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

  }

  //! Initialize only Flux (Area read from OneDNonLinParam)
  void OneDModelSolver::initialize(const Real& u20)
  {
    _M_U2_thistime = ScalarVector( _M_U2_thistime.size(), u20 );

    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
    	_M_U1_thistime[ielem]=_M_oneDParam.Area0(ielem);

		_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );
    }
    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step ); 

    //! create matlab scripts  
    create_movie_file();
	
    openFileBuffers();
	
    output2FileBuffers( 0. );
	
    postProcess( 0. );
  
  }


  /* Initialize with non constant data:
     parameter var = non constant variable
     1. for Area
     2. for Flux
     3. for W1
     4. for W2
     parameter value1 = outside the step
     parameter value2 = inside the step
     parameter value3 = constant variable
     if var = 1, constant flux = value3 is assumed
     if var = 2, constant area = value3 is assumed
     if var = 3, constant W2 = value3 is assumed
     if var = 4, constant W1 = value3 is assumed
     parameter firstnode, lastnode = edges of the step
     hypothesis: constant parameters along the vessel
  */
  void OneDModelSolver::initialize(GetPot& data_file)
  //const std::string& var,
  //					const Real& value1, const Real& value2, const Real& value3,
//				   	const UInt& firstnode, const UInt& lastnode)
  {
  	// the discontinuity is comprised between firstnode and lastnode
   	UInt firstnode( data_file("initialize/firstnode",1) );
   	UInt lastnode( data_file("initialize/lastnode",2) );
      
    if( firstnode > lastnode )
      {
        std::cerr<<" Error in initialize: incorrect step limits" <<std::endl;
        abort();
      }
    if( lastnode > _M_dimDof )
      {
        std::cerr<<" Error in initialize: step outside tube boundaries" <<std::endl;
        abort();
      }

      // read initialization type from data file (see OneDModelSolver::initialize)
      std::string init_var( data_file("initialize/var","pressure") );
      Real multiplier( data_file("initialize/multiplier",1.) );
      
      // tell me what I am doing
      Debug( 6030 ) << "[OneDModelSolver::initialize] 0- Initializing with values:\n";
      Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tinitialize var = " << init_var << "\n";
	  Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tfirstnode = " << firstnode << "\n";
	  Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tlastnode = " << lastnode << "\n";
	  Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tmultiplier = " << multiplier << "\n";

      Real value1, value2;

      switch( _M_oneDstring2initializeVarMap[init_var] ){
	// case 1, 2: initialize physical variables to desired value	
      case OneDInitPressure:
		value1 = _M_oneDParam.A_from_P( data_file("initialize/value",0.) );
		value2 = 0;
    	_M_U1_thistime = ScalarVector( _M_U1_thistime.size(), value1 );
    	_M_U2_thistime = ScalarVector( _M_U2_thistime.size(), value2 );

		for (UInt inode=firstnode; inode <= lastnode ; ++inode ) {
	  		_M_U1_thistime[inode - 1] = value1 * multiplier;
		}
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
			_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );
      	}
	break;

      case OneDInitArea:
		value1 = data_file("initialize/value",0.);
		value2 = 0;
		
		_M_U1_thistime = ScalarVector( _M_U1_thistime.size(), value1 );
    	_M_U2_thistime = ScalarVector( _M_U2_thistime.size(), value2 );

		for (UInt inode=firstnode; inode <= lastnode ; ++inode )
	  		_M_U1_thistime[inode - 1] = value1 * multiplier;
	  		
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );

	break;

      case OneDInitFlux:
		value1 = 0;
		value2 = data_file("initialize/value",0.);
    	_M_U1_thistime = ScalarVector( _M_U1_thistime.size(), value2 );
    	_M_U2_thistime = ScalarVector( _M_U2_thistime.size(), value1 );

		for (UInt inode=firstnode; inode <= lastnode ; ++inode )
	  		_M_U2_thistime[inode - 1] = value1 * multiplier;
	  		
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );

	break;

      case OneDInitReimann1:
		value1 = data_file("initialize/value",0.);
		value2 = -value1;
		std::cout << "[OneDModelSolver::initialize] WARNING! Initializing W2 = - W1"
				<< " (assuming Q = 0)" << std::endl;
    	_M_W1_thistime = ScalarVector( _M_W1_thistime.size(), value1 );
    	_M_W2_thistime = ScalarVector( _M_W2_thistime.size(), value2 );

		for (UInt inode=firstnode; inode <= lastnode ; ++inode )
		{
	  		_M_W1_thistime[inode - 1] = value1 * multiplier;
	  		_M_W2_thistime[inode - 1] = value2 * multiplier;
      	}	  		
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.U_from_W( _M_U1_thistime[ielem], _M_U2_thistime[ielem],
				_M_W1_thistime[ielem], _M_W2_thistime[ielem], ielem );

	break;

      case OneDInitReimann2:
		value1 = -data_file("initialize/value",0.);
		value2 = -value1;
		std::cout << "[OneDModelSolver::initialize] WARNING! Initializing W1 = - W2"
				<< " (assuming Q = 0)" << std::endl;
    	_M_W1_thistime = ScalarVector( _M_W1_thistime.size(), value2 );
    	_M_W2_thistime = ScalarVector( _M_W2_thistime.size(), value1 );

		for (UInt inode=firstnode; inode <= lastnode ; ++inode )
		{
	  		_M_W2_thistime[inode - 1] = value1 * multiplier;
	  		_M_W1_thistime[inode - 1] = value2 * multiplier;
		}
      	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
			_M_oneDParam.U_from_W( _M_U1_thistime[ielem], _M_U2_thistime[ielem],
				_M_W1_thistime[ielem], _M_W2_thistime[ielem], ielem );

	break;

      default:	
        ERROR_MSG("No such initializing option.");

      }

	Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tvalue1 = " << value1 << "\n";
	Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tvalue1_step = " << value1 * multiplier << "\n";
    Debug( 6030 ) << "[OneDModelSolver::initialize]\t\tvalue2 = " << value2 << "\n";

    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step ); 

    //! create matlab scripts  
    create_movie_file();
	
    openFileBuffers();
    output2FileBuffers( 0. );
	
    postProcess( 0. );
  
  }


  void OneDModelSolver::savesol()
  {
    _M_U1_timestep = _M_U1_thistime;
    _M_U2_timestep = _M_U2_thistime;
  }

  // here i need to keep the boundary values!!
  void OneDModelSolver::loadsol()
  {
    Vec2D U_leftbd   = Vec2D ( _M_U1_thistime( 0 ) ,
			       _M_U2_thistime( 0 ) );
    Vec2D U_rightbd   = Vec2D ( _M_U1_thistime( _M_rhs1.size()-1 ) ,
				_M_U2_thistime( _M_rhs1.size()-1 ) );

    _M_U1_thistime = _M_U1_timestep;
    _M_U2_thistime = _M_U2_timestep;

    _M_U1_thistime( 0 ) = U_leftbd.first;
    _M_U2_thistime( 0 ) = U_leftbd.second;
    _M_U1_thistime(_M_rhs1.size()-1) = U_rightbd.first;
    _M_U2_thistime(_M_rhs1.size()-1) = U_rightbd.second;

  	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
		_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );

  }


  //! Update the right hand side for time advancing
  void OneDModelSolver::timeAdvance( const Real& time_val )
  {
    Chrono chrono;

    Real dt2over2 = _M_time_step * _M_time_step * 0.5;

    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] o-  updates of flux and sources... " << "\n";
    chrono.start();

    //! output cfl
    CheckCFL();

    //! update the vector containing the values of the flux at the nodes
    //! and its jacobian
    _updateFluxDer();
    //! update the vector containing the values of the source term at the nodes
    //! and its jacobian
    _updateSourceDer();
    chrono.stop();
    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] \tdone in " << chrono.diff() << " s.\n";
    

    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] o-  updates of matrices... " << "\n";
    chrono.start();
    //! update the matrices for the non-linear terms
    _updateMatrices();
    chrono.stop();
    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] \tdone in " << chrono.diff() << " s.\n";

    /*!
      ---------------------------------------------------
      Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T )
      ---------------------------------------------------

      (Un+1, phi) =                             //! massFactor^{-1} * Un+1
      (Un, phi)                         //!            mass * U
      + dt     * (       Fh(Un), dphi/dz )         //!            grad * F(U)
      - dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     //! gradDiffFlux(U) * S(U)
      + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //!   divDiffSrc(U) * F(U)
      - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //!stiffDiffFlux(U) * F(U)
      - dt     * (       Sh(Un), phi )             //!            mass * S(U)
      + dt^2/2 * (diffSh(Un) Sh(Un), phi )         //!  massDiffSrc(U) * S(U)
      ---------------------------------------------------
    */

    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] o-  Matrix vector products... " << "\n";
    chrono.start();

    //! Reminder of the function Axpy:
    //! Axpy(alpha, x, beta, y) -> y = alpha*A*x + beta*y

    //!---------------------------------------------------
    //! 1/ compute _M_rhs1 (system in U1=A)
    //!---------------------------------------------------
    //! rhs1 = mass * Un1
    _M_massMatrix.Axpy( 1., _M_U1_thistime , 0., _M_rhs1 );

    //! rhs1 = rhs1 + dt * grad * F1(Un)
    _M_gradMatrix.Axpy( _M_time_step, _M_Flux1 , 1., _M_rhs1 );

    //! rhs1 = rhs1 - dt^2/2 * gradDiffFlux11 * S1(Un)
    _M_gradMatrixDiffFlux11.Axpy( -dt2over2, _M_Source1 , 1., _M_rhs1 );
    //! rhs1 = rhs1 - dt^2/2 * gradDiffFlux12 * S2(Un)
    _M_gradMatrixDiffFlux12.Axpy( -dt2over2, _M_Source2 , 1., _M_rhs1 );

    //! rhs1 = rhs1 + dt^2/2 * divDiffSrc11 * F1(Un)
    _M_divMatrixDiffSrc11.Axpy( dt2over2, _M_Flux1 , 1., _M_rhs1 );
    //! rhs1 = rhs1 + dt^2/2 * divDiffSrc12 * F2(Un)
    _M_divMatrixDiffSrc12.Axpy( dt2over2, _M_Flux2 , 1., _M_rhs1 );

    //! rhs1 = rhs1 - dt^2/2 * stiffDiffFlux11 * F1(Un)
    _M_stiffMatrixDiffFlux11.Axpy( -dt2over2, _M_Flux1 , 1., _M_rhs1 );
    //! rhs1 = rhs1 - dt^2/2 * stiffDiffFlux12 * F2(Un)
    _M_stiffMatrixDiffFlux12.Axpy( -dt2over2, _M_Flux2 , 1., _M_rhs1 );

    //! rhs1 = rhs1 - dt * mass * S1(Un)
    _M_massMatrix.Axpy( - _M_time_step, _M_Source1 , 1., _M_rhs1 );

    //! rhs1 = rhs1 + dt^2/2 * massDiffSrc11 * S1(Un)
    _M_massMatrixDiffSrc11.Axpy( dt2over2, _M_Source1 , 1., _M_rhs1 );
    //! rhs1 = rhs1 + dt^2/2 * divDiffSrc12 * S2(Un)
    _M_massMatrixDiffSrc12.Axpy( dt2over2, _M_Source2 , 1., _M_rhs1 );


    //!---------------------------------------------------
    //! 2/ compute _M_rhs2 (system in U2=Q)
    //!---------------------------------------------------
    //! rhs2 = mass * Un2
    _M_massMatrix.Axpy( 1., _M_U2_thistime , 0., _M_rhs2 );

    //! rhs2 = rhs2 + dt * grad * F2(Un)
    _M_gradMatrix.Axpy( _M_time_step, _M_Flux2 , 1., _M_rhs2 );

    //! rhs2 = rhs2 - dt^2/2 * gradDiffFlux21 * S1(Un)
    _M_gradMatrixDiffFlux21.Axpy( -dt2over2, _M_Source1 , 1., _M_rhs2 );
    //! rhs2 = rhs2 - dt^2/2 * gradDiffFlux22 * S2(Un)
    _M_gradMatrixDiffFlux22.Axpy( -dt2over2, _M_Source2 , 1., _M_rhs2 );

    //! rhs2 = rhs2 + dt^2/2 * divDiffSrc21 * F1(Un)
    _M_divMatrixDiffSrc21.Axpy( dt2over2, _M_Flux1 , 1., _M_rhs2 );
    //! rhs2 = rhs2 + dt^2/2 * divDiffSrc22 * F2(Un)
    _M_divMatrixDiffSrc22.Axpy( dt2over2, _M_Flux2 , 1., _M_rhs2 );

    //! rhs2 = rhs2 - dt^2/2 * stiffDiffFlux21 * F1(Un)
    _M_stiffMatrixDiffFlux21.Axpy( -dt2over2, _M_Flux1 , 1., _M_rhs2 );
    //! rhs2 = rhs2 - dt^2/2 * stiffDiffFlux22 * F2(Un)
    _M_stiffMatrixDiffFlux22.Axpy( -dt2over2, _M_Flux2 , 1., _M_rhs2 );

    //! rhs2 = rhs2 - dt * mass * S2(Un)
    _M_massMatrix.Axpy( - _M_time_step, _M_Source2 , 1., _M_rhs2 );

    //! rhs2 = rhs2 + dt^2/2 * massDiffSrc21 * S1(Un)
    _M_massMatrixDiffSrc21.Axpy( dt2over2, _M_Source1 , 1., _M_rhs2 );
    //! rhs2 = rhs2 + dt^2/2 * divDiffSrc22 * S2(Un)
    _M_massMatrixDiffSrc22.Axpy( dt2over2, _M_Source2 , 1., _M_rhs2 );

    //!---------------------------------------------------
    //! 3/ take into account the BOUNDARY CONDITIONS
    //!---------------------------------------------------
    //! compute the values for the boundary conditions
    _computeBC( time_val );
    //! take into account the bc
    _updateBCDirichletVector();

    // *******************************************************
    chrono.stop();
    Debug( 6030 ) << "[OneDModelSolver::timeAdvance] \trhs computed in " << chrono.diff() << " s.\n";


  }



  void OneDModelSolver::iterate( const Real& time_val , const int& count)
  {
    Debug( 6030 ) << "[OneDModelSolver::iterate] o-  Solving the system... t = " << time_val
		  << ", iter = " << count  << "... \n";
    Chrono chrono;
    chrono.start();

  	if( _M_UW ){
	    Real Ainode, Qinode;
	
	    Real lambda1_plus = 0.;
	    Real lambda2_plus = 0.;
	    Real lambda1_minus = 0.;
	    Real lambda2_minus = 0.;
	    Real eigval1, eigval2;
	    Real tmp11, tmp12, tmp21, tmp22;

  		Real deltaX = _M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x();

   		//! working on riemann invariants
   		ScalUnknown<Vector> W1_UW(_M_dimDof);
   		ScalUnknown<Vector> W2_UW(_M_dimDof);

		//! converting boundary conditions on physical variables
		//! in boundary conditions on characteristic variables
		_M_oneDParam.W_from_U( W1_UW(0), W2_UW(0),
				_M_rhs1[0], _M_rhs2[0], 0 );
		_M_oneDParam.W_from_U( W1_UW(_M_dimDof-1), W2_UW(_M_dimDof-1),
				_M_rhs1[_M_dimDof-1], _M_rhs2[_M_dimDof-1], _M_dimDof-1 );
  	
    	for ( UInt ii=1; ii < (_M_dimDof-1) ; ii++ ) {
    		//! compute the eigenvalues at node
		    Ainode = _M_U1_thistime( ii );
      		Qinode = _M_U2_thistime( ii );
      		_M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
					       eigval1, eigval2,
					       tmp11, tmp12,
					       tmp21, tmp22,
					       ii );

		    lambda1_plus = std::max<Real>( eigval1, 0. );
		    lambda1_minus = std::min<Real>( eigval1, 0. );
		    lambda2_plus = std::max<Real>( eigval2, 0. );
		    lambda2_minus = std::min<Real>( eigval2, 0. );
   		//! update the solution for the next time step
    		W1_UW[ii] = _M_W1_thistime[ii]
    			- (_M_time_step / deltaX) * lambda1_plus * ( _M_W1_thistime[ii] - _M_W1_thistime[ii-1]) 
    			- (_M_time_step / deltaX) * lambda1_minus * ( _M_W1_thistime[ii+1] - _M_W1_thistime[ii]) 
    			- _M_time_step * ( tmp11 * _M_Source1[ii] + tmp12 * _M_Source2[ii] );
    		W2_UW[ii] = _M_W2_thistime[ii]
    			- (_M_time_step / deltaX) * lambda2_plus * ( _M_W2_thistime[ii] - _M_W2_thistime[ii-1]) 
    			- (_M_time_step / deltaX) * lambda2_minus * ( _M_W2_thistime[ii+1] - _M_W2_thistime[ii]) 
    			- _M_time_step * ( tmp21 * _M_Source1[ii] + tmp22 * _M_Source2[ii] );
    	}
    	
    _M_W1_thistime = W1_UW;
    _M_W2_thistime = W2_UW;
    
   	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
		_M_oneDParam.U_from_W( _M_U1_thistime[ielem], _M_U2_thistime[ielem],
				_M_W1_thistime[ielem], _M_W2_thistime[ielem],
				ielem );

  	}
  	else{
    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs1 );

    //! solve the system: rhs2 = massFactor^{-1} * rhs2
    _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs2 );

    //! update the solution for the next time step
    _M_U1_thistime = _M_rhs1;
    _M_U2_thistime = _M_rhs2;

   	for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
		_M_oneDParam.W_from_U( _M_W1_thistime[ielem], _M_W2_thistime[ielem],
				_M_U1_thistime[ielem], _M_U2_thistime[ielem], ielem );
  	}
  	
    chrono.stop();
    Debug( 6030 ) << "[OneDModelSolver::iterate] \tdone in " << chrono.diff() << " s.\n";

    output2FileBuffers( time_val );
	
  }



  void OneDModelSolver::openFileBuffers()
  {
    std::string file_output;

    boost::shared_ptr<std::ostringstream> buf;

    file_output = _M_post_dir + "/" + "Area" + _M_post_file + ".m";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + "Portata" + _M_post_file + ".m";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + "Reimann1_" + _M_post_file + ".m";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + "Reimann2_" + _M_post_file + ".m";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );


  
    file_output = _M_post_dir + "/" + _M_post_file + "A.mtv";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + _M_post_file + "Q.mtv";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + _M_post_file + "W1.mtv";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );

    file_output = _M_post_dir + "/" + _M_post_file + "W2.mtv";

    buf.reset( new std::ostringstream("") );
    _M_post_process_buffer.insert( std::map<std::string, boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
    _M_post_process_buffer_offset.insert( std::map<std::string, long>::value_type( file_output, buf->tellp() ) );


  }



  void OneDModelSolver::output2FileBuffers( const Real& time_val )//, const std::string str )
  {
  
    if( !( static_cast<int>( std::floor( (time_val/_M_time_step) + .5 ) ) %  _M_verbose ) ){

   
      std::string file_output;
    
      file_output = _M_post_dir + "/" + "Area" + _M_post_file + ".m";
    
      output_to_matlab( file_output, time_val, _M_U1_thistime, "A" );
	      
      file_output = _M_post_dir + "/" + "Portata" + _M_post_file + ".m";
    
      output_to_matlab( file_output, time_val, _M_U2_thistime, "Q" );
    
      file_output = _M_post_dir + "/" + "Reimann1_" + _M_post_file + ".m";
    
      output_to_matlab( file_output, time_val, _M_W1_thistime, "W1" );
    
      file_output = _M_post_dir + "/" + "Reimann2_" + _M_post_file + ".m";
    
      output_to_matlab( file_output, time_val, _M_W2_thistime, "W2" );
    
    
    
      file_output = _M_post_dir + "/" + _M_post_file + "A.mtv";
    
      output_to_plotmtv( file_output, time_val, _M_mesh.pointList(), _M_U1_thistime, "A" );
    
      file_output = _M_post_dir + "/" + _M_post_file + "Q.mtv";
    
      output_to_plotmtv( file_output, time_val, _M_mesh.pointList(), _M_U2_thistime, "Q" );
    
      file_output = _M_post_dir + "/" + _M_post_file + "W1.mtv";
    
      output_to_plotmtv( file_output, time_val, _M_mesh.pointList(), _M_W1_thistime, "W1" );
    
      file_output = _M_post_dir + "/" + _M_post_file + "W2.mtv";
    
      output_to_plotmtv( file_output, time_val, _M_mesh.pointList(), _M_W2_thistime, "W2" );


    } 
  
  }



  void OneDModelSolver::closeFileBuffers()
  {
    // as I have a boost::shared_ptr, I expect the objects to be deallocated now that the pointers are destroyed
    _M_post_process_buffer.erase( _M_post_process_buffer.begin(), _M_post_process_buffer.end() );
    _M_post_process_buffer_offset.erase( _M_post_process_buffer_offset.begin(), _M_post_process_buffer_offset.end() );
  }



  void OneDModelSolver::postProcess( const Real& time_val )
  {

    std::string str;

    std::ofstream outfile;
    
    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;
    
//    UInt time_index( static_cast<UInt>( std::floor( time_val / _M_postProcessTimeStep + 0.5 ) ) );
    
    if (time_val==0.){
      // the code is entering this for sure, as
      // initialize invokes postProcess( 0. )
      
      std::string file_output;
      
      Real deltax( _M_oneDParam.Length() / static_cast<Real>(_M_oneDParam.ParamSize() - 1) );
      
      file_output = "dati.m";
      outfile.open(file_output.c_str(), std::ios::app );
      outfile << "z" << _M_post_file
	      << " = (" << _M_x_left
	      << ":" << deltax
	      << ":" << _M_x_right << ");\n"
	      << std::endl;
      outfile.close();

    }
    
    
    Debug( 6030 ) << "[OneDModelSolver::postProcess] o- Dumping solutions on files (1d)!" << "\n";
    // dump solutions on files (buffers must be active!)
    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end(); iter++ ){
      outfile.open( (*iter).first.c_str(), std::ios::app );
      outfile << (*(*iter).second).str();
      outfile.close();
    }
    
    
    resetFileBuffers();
    
    
  };



  void OneDModelSolver::resetFileBuffers( )
  {
    closeFileBuffers();
    openFileBuffers();
  }



  void OneDModelSolver::seekpFileBuffers( )
  {

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;
  
    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end(); iter++ ){
      (*iter).second->seekp( _M_post_process_buffer_offset[(*iter).first] );

    }

  }



  void OneDModelSolver::tellpFileBuffers( )
  {

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;
  
    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end(); iter++ ){
      _M_post_process_buffer_offset[(*iter).first] = (*iter).second->tellp();

    }

  }



  void OneDModelSolver::gplot( )
  {
    _M_GracePlot.Plot( _M_mesh.pointList(), _M_U1_thistime );
  }

  //! output for Plotmtv.
  void OneDModelSolver::output_to_plotmtv(std::string fname, Real time_val,
					  const std::vector< Point1D >& ptlist,
					  const ScalUnknown<Vector>& U,
					  std::string vname )
  {
 
    boost::shared_ptr<std::ostringstream> buf;

    buf = _M_post_process_buffer[fname];

 
    (*buf) << "$ DATA = CURVE2D\n % xlabel='z'\n"
	   << "% toplabel='Section,time=" << time_val << "'\n % ylabel='" << vname << "'\n"
	   << std::flush;

    for(UInt ii = 0; ii < U.size(); ii++){
      (*buf) << ptlist[ii].x() << " " << U[ii] << "\n" << std::flush;
    }


  }

  //! output for Matlab.
  void OneDModelSolver::output_to_matlab( std::string fname, Real time_val,
					  const ScalUnknown<Vector>& U,
					  std::string vname )
  {

    boost::shared_ptr<std::ostringstream> buf;

    buf = _M_post_process_buffer[fname];

    (*buf) << vname << _M_post_file
	   << "( " << (static_cast<int>( std::floor( time_val/_M_time_step + 0.5 ) ) /  _M_verbose )+1
	   << ", : ) = [ " << std::flush;
  
    for ( UInt ii=LeftNodeId(); ii <= RightNodeId() ; ++ii ) {
      (*buf) << U[ii] << "; ";
    }
    (*buf) << "]';\n" << std::endl;


  }

  //! Create a Matlab script to visualize output matlab files
  void OneDModelSolver::create_movie_file(int /*i*/)
  {
    std::ofstream outfile;
    std::string file_output;
    file_output = _M_post_dir + "/" + "areamovie"+_M_post_file+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Area"<<_M_post_file<<";\n"
	    <<"[n,m]=size(A" << _M_post_file << ");\n"
	    <<"Max=max(max(A" << _M_post_file << "));\n"
	    <<"Min=min(min(A" << _M_post_file << "));\n"
	    <<"for i=1:n\n"
	    <<"  plot(A" << _M_post_file << "(i,:));\n"
	    <<"  title((i-1)*"<<_M_time_step<<"*"<<_M_verbose<<");\n"
	    <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
	    <<"  %pause;\n"
	    <<"  F(i) = getframe;\n"
	    <<"end\n"
	    <<"movie(F)";
    
    outfile.close();
    file_output = _M_post_dir + "/" + "portatamovie"+_M_post_file+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Portata"<<_M_post_file<<";\n"
	    <<"[n,m]=size(Q" << _M_post_file << ");\n"
	    <<"Max=max(max(Q" << _M_post_file << "));\n"
	    <<"Min=min(min(Q" << _M_post_file << "));\n"
	    <<"for i=1:n\n"
	    <<"  plot(Q" << _M_post_file << "(i,:));\n"
	    <<"  title((i-1)*"<<_M_time_step<<"*"<<_M_verbose<<");\n"
	    <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
	    <<"  %pause;\n"
	    <<"  F(i) = getframe;\n"
	    <<"end\n"
	    <<"movie(F)";
    outfile.close();

  }

}
