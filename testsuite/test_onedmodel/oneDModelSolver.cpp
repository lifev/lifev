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

  \brief This file implements a solver for 1D model.
*/

#include "oneDModelSolver.hpp"


namespace LifeV
{
OneDModelSolver::OneDModelSolver(const GetPot& data_file):
  OneDModelHandler(data_file),
  _M_elmatMass (_M_fe.nbNode,1,1),
  _M_elmatStiff(_M_fe.nbNode,1,1),
  _M_elmatGrad (_M_fe.nbNode,1,1),
  _M_elmatDiv  (_M_fe.nbNode,1,1),
  _M_elvec     (_M_fe.nbNode,1),
  _M_U_thistime(_M_dimDof),
  _M_U_nexttime(_M_dimDof),
  _M_FluxU(_M_dimDof),
  _M_rhs(_M_dimDof),
  _M_massMatrix(_M_dimDof),
  _M_stiffMatrix(_M_dimDof),
  _M_gradMatrix(_M_dimDof),
  _M_divMatrix(_M_dimDof),
  _M_factorMassMatrix(_M_dimDof),
  _M_massupdiag2( _M_dimDof - 2 ),
  _M_massipiv( _M_dimDof )
{

  cout << endl;
  cout << "O-  Nb of unknowns: " << _M_dimDof     << endl;
  cout << "O-  Computing mass matrix... \n";

  Chrono chrono;
  chrono.start();

  //! Matrices initialization
  _M_massMatrix.zero();
  _M_stiffMatrix.zero();
  _M_gradMatrix.zero();
  _M_divMatrix.zero();
  _M_factorMassMatrix.zero();

  //  _M_massMatrix.showMe(std::cout, _M_verbose);

  _M_rhs = 1.;
  _M_bcDirLeft  = 1.;
  _M_bcDirRight = 0;

  _M_coeffMass  = 1.;
  _M_coeffStiff = 1.;
  _M_coeffGrad  = 1.;
  _M_coeffDiv   = 1.;

  //! Elementary computation and matrix assembling
  //! Loop on elements
  for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

    //! update _M_fe and _M_elmat*
    _updateElemMatrices( iedge );

    //! assemble the mass matrix
    assemb_mat( _M_massMatrix, _M_elmatMass, _M_fe, _M_dof , 0, 0 );

    //! assemble the stiffness matrix
    assemb_mat( _M_stiffMatrix, _M_elmatStiff, _M_fe, _M_dof , 0, 0 );

    //! assemble the gradient matrix
    assemb_mat( _M_gradMatrix, _M_elmatGrad, _M_fe, _M_dof , 0, 0 );

    //! assemble the divergence matrix
    assemb_mat( _M_divMatrix, _M_elmatDiv, _M_fe, _M_dof , 0, 0 );
  } //! end loop on elements

  //! Dirichlet boundary conditions set in matrices
  _updateBCDirichletMatrix( _M_massMatrix );
  _updateBCDirichletMatrix( _M_stiffMatrix );
  _updateBCDirichletMatrix( _M_gradMatrix );
  _updateBCDirichletMatrix( _M_divMatrix );

  _M_factorMassMatrix = _M_massMatrix;
  //! factorization of the mass matrix
  _factorizeMassMatrix();

  /*
  cout << "\n\n\tMass matrix " << endl;
  _M_massMatrix.showMe( std::cout , _M_verbose );
  cout << "\n\n\tStiffness matrix " << endl;
  _M_stiffMatrix.showMe( std::cout , _M_verbose );
  cout << "\n\n\tGradient matrix " << endl;
  _M_gradMatrix.showMe( std::cout , _M_verbose );
  cout << "\n\n\tDivergence matrix " << endl;
  _M_divMatrix.showMe( std::cout , _M_verbose );
  cout << "\n\n\tFACTORIZED Mass matrix " << endl;
  _M_factorMassMatrix.showMe( std::cout , _M_verbose );

  ScalUnknown<Vector> vec( _M_dimDof );
  for (int ii=0; ii < _M_stiffMatrix.OrderMatrix() ; ii++ ) {
    vec( ii ) = ii;
    cout <<  ii << " " << vec(ii) << endl;;
  }

  _solveMassMatrix( vec );
  cout << "solve " << endl;
  for (int ii=0; ii < _M_stiffMatrix.OrderMatrix() ; ii++ ) {
    cout <<  ii << " " << vec(ii) << endl;;
  }

  vec = 1000.;
  _M_stiffMatrix.Axpy( 1., _M_rhs , 1., vec );
  cout << "matvec " << endl;
  for (int ii=0; ii < _M_stiffMatrix.OrderMatrix() ; ii++ ) {
    cout <<  ii << " " << vec(ii) << endl;;
  }
  */

  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;

}

//! Update the element matrices with the current element (from 1)
void OneDModelSolver::_updateElemMatrices( const UInt& iedge )
{
  //! set the elementary matrices to 0.
  _M_elmatMass.zero();
  _M_elmatStiff.zero();
  _M_elmatGrad.zero();
  _M_elmatDiv.zero();

  //! update the current element
  _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
  //  std::cout << _M_fe.currentId() << std::endl;

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

/*! modify the matrix to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)
*/
void OneDModelSolver::
_updateBCDirichletMatrix( TriDiagMatrix<double>& mat )
{
  UInt firstDof = 0;
  UInt lastDof  = mat.OrderMatrix()-1;
  //! modify the first row
  mat.Diag()( firstDof )   = 1.;
  mat.UpDiag()( firstDof ) = 0.;
  //! modify the last row

  //! FOR NEUMANN!!!!!!
  //  mat.Diag()( lastDof )      = 1.;
  //  mat.LowDiag()( lastDof-1 ) = 0.;


}

/*! modify the vector to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)

  \param vec : the rhs vector
  \param val_left  : Dirichlet value inserted to the left
  \param val_right : Dirichlet value inserted to the right
*/
void  OneDModelSolver::
_updateBCDirichletVector( ScalUnknown<Vector>& vec,
			  const double& val_left,
			  const double& val_right )
{
  UInt firstDof = 0;
  UInt lastDof  = vec.size()-1;
  //! first row and last row are modified
  vec( firstDof ) = val_left;
  //! FOR NEUMANN!!!!!!
  //  vec( lastDof  ) = val_right;
}

//! update the flux from the current unknown: FfluxU = F_h(U_h^n)
void OneDModelSolver::_updateFlux()
{
  double celerity = 2.;
  for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
    _M_FluxU( ii ) = celerity * _M_U_thistime( ii );
  }
}

/*! LU factorize with lapack _M_factorMassMatrix
  SUBROUTINE DGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
  BEWARE: it modifies _M_factorMassMatrix!
*/
void OneDModelSolver::_factorizeMassMatrix()
{
  int INFO = 0;
  int OrderMat =_M_factorMassMatrix.OrderMatrix();

  //! solve with lapack (for tridiagonal matrices)
  dgttrf_( &OrderMat, _M_factorMassMatrix.LowDiag(), _M_factorMassMatrix.Diag(),
	   _M_factorMassMatrix.UpDiag(), _M_massupdiag2, _M_massipiv, &INFO);
  ASSERT_PRE(!INFO,"Lapack factorization of tridiagonal matrix not achieved.");

}
/*! lapack LU solve AFTER FACTORIZATION of _M_factorMassMatrix
  SUBROUTINE DGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, INFO )

*  TRANS   (input) CHARACTER
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)

*/
void OneDModelSolver::_solveMassMatrix( ScalUnknown<Vector>& vec )
{
  int INFO = 0;
  int NBRHS = 1;//  nb columns of the vec := 1.
  int OrderMat = _M_factorMassMatrix.OrderMatrix();

  ASSERT_PRE( OrderMat == (int)vec.size() ,
	      "The right-hand side must have the same dimensions as the tridiag matrix.");

  //! solve with lapack (for tridiagonal matrices)
  dgttrs_( "N", &OrderMat, &NBRHS, _M_factorMassMatrix.LowDiag(), _M_factorMassMatrix.Diag(),
	   _M_factorMassMatrix.UpDiag(), _M_massupdiag2, _M_massipiv,
	  vec.giveVec(), &OrderMat, &INFO);
  ASSERT_PRE(!INFO,"Lapack solve of tridiagonal matrix not achieved.");

}

/*! direct LU solve with lapack _M_factorMassMatrix
  (use it once as it changes the matrix !)
  SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
*/
void OneDModelSolver::_directsolveMassMatrix( ScalUnknown<Vector>& vec )
{
  int INFO = 0;
  int NBRHS = 1;//  nb columns of the vec := 1.
  int OrderMat = _M_factorMassMatrix.OrderMatrix();

  ASSERT_PRE( OrderMat == (int)vec.size() ,
	      "The right-hand side must have the same dimensions as the tridiag matrix.");

  //! solve with lapack (for tridiagonal matrices)
  dgtsv_( &OrderMat, &NBRHS, _M_factorMassMatrix.LowDiag(), _M_factorMassMatrix.Diag(),
	   _M_factorMassMatrix.UpDiag(),
	  vec.giveVec(), &OrderMat, &INFO);
  ASSERT_PRE(!INFO,"Lapack solve of tridiagonal matrix not achieved.");

}

//! Update the right hand side for time advancing
void OneDModelSolver::timeAdvance()
{
  cout << "  o-  Updating right hand side... ";

  Chrono chrono;
  chrono.start();

  double dt2over2 = _M_time_step * _M_time_step * 0.5;

  //! _M_FluxU = F_h( U_h^n )
  _updateFlux();

  //! Reminder of the function Axpy:
  //! Axpy(alpha, x, beta, y) -> y = alpha*A*x + beta*y

  //! rhs = mass * Un
  _M_massMatrix.Axpy( 1., _M_U_thistime , 0., _M_rhs );

  //! rhs = rhs + dt * grad * F_h(Un)
  _M_gradMatrix.Axpy( _M_time_step, _M_FluxU , 1., _M_rhs );

  //! rhs = rhs - dt^2/2 * stiff * F_h(Un)
  _M_stiffMatrix.Axpy( -dt2over2, _M_FluxU , 1., _M_rhs );

  //! take into account the bc
  _updateBCDirichletVector( _M_rhs, _M_bcDirLeft, _M_bcDirRight );

  // *******************************************************
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;
}

void OneDModelSolver::iterate()
{
  cout << "  o-  Solving the system... ";
  Chrono chrono;
  chrono.start();

  //! solve the mass matrix and return the result in _M_rhs
  _solveMassMatrix( _M_rhs );

  /*
  cout << "\n\tsolution at time n+1 " << endl;
  for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
    cout <<  ii << " " << _M_rhs(ii) << endl;;
  }
  */

  //! solution for the next time step
  _M_U_thistime = _M_rhs;
   // *******************************************************
  chrono.stop();
  cout << "done in " << chrono.diff() << " s." << endl;

}


void OneDModelSolver::gplot( )
{
  _M_GracePlot.Plot( _M_mesh.pointList(), _M_U_thistime );
}
}
