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

#include "oneDModelSolver.hpp"


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
    //! elementary matrices
    _M_elmatMass (_M_fe.nbNode,1,1),
    _M_elmatStiff(_M_fe.nbNode,1,1),
    _M_elmatGrad (_M_fe.nbNode,1,1),
    _M_elmatDiv  (_M_fe.nbNode,1,1),
    //! vectorial unknowns and rhs
    _M_U1_thistime(_M_dimDof),
    _M_U2_thistime(_M_dimDof),
    _M_rhs1(_M_dimDof),
    _M_rhs2(_M_dimDof),
    //! vectors and matrices of the non-linear function
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
    //! mass matrix (to be inverted)
    _M_massMatrix(_M_dimDof),
    _M_factorMassMatrix(_M_dimDof),
    _M_tridiagSlv(_M_dimDof),
    //! matrices used to build the rhs
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
    _M_divMatrixDiffSrc22(_M_dimDof)
{

    std::cout << "\n";
    std::cout << "O-  Nb of unknowns: " << _M_dimDof     << "\n";
    std::cout << "O-  Computing the constant matrices... \n" << std::endl;

    std::string fname1 = _M_post_dir + "/" + _M_post_file + "A.mtv";
    std::string fname2 = _M_post_dir + "/" + _M_post_file + "Q.mtv";
    output_to_plotmtv( fname1, 0., _M_mesh.pointList(), _M_U1_thistime, 0);
    output_to_plotmtv( fname2, 0., _M_mesh.pointList(), _M_U2_thistime, 0);

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

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}

/*! Update the coefficients
  (from the flux, source functions and their derivatives)
*/
void OneDModelSolver::
_updateMatrixCoefficients(const UInt& ii, const UInt& jj , const UInt& iedge)
{
    Real dFluxdUelem = 0;
    Real dSrcdUelem = 0;

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


//! compute the _M_bcDirLeft and _M_bcDirRight
void OneDModelSolver::_computeBCValues( const Real& time_val )
{
    //! node on the left boundary
    UInt leftDof = 0;
    //! node on the right boundary
    UInt rightDof  = _M_rhs1.size()-1;

    Edge1D leftedge  = _M_mesh.edgeList( 1 ); //! from 1 to nb_elem
    Edge1D rightedge = _M_mesh.edgeList( _M_nb_elem ); //! from 1 to nb_elem

    //! eigen values of the jacobian diffFlux (= dF/dU)
    Real  eigval1, eigval2;
    //! left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D left_eigvec1, left_eigvec2;

    /*! interpolate the eigenvectors
      Vec2D left_eigvec1_internalBd, left_eigvec2_internalBd;
      Vec2D left_eigvec1_charact_pt, left_eigvec2_charact_pt;
      Real eigval1_internalBd, eigval2_internalBd;
    */

    //! first and second line of the identity matrix
    Vec2D line1_id( 1. , 0. );
    Vec2D line2_id( 0. , 1. );

    //! right hand side for the 2x2 linear system to be solved at each side
    Vec2D rhsBC;
    //! result of the 2x2 linear system to be solved at each side
    Vec2D resBC;

    //! bc values of U at one side
    Vec2D UBC_data;

    //! quasi linear source term
    Vec2D qlSource;

    //! value of U at the boundary, at the neighboring internal node
    //!            and at the foot of the characteristics
    Vec2D U_boundary, U_internalBd, U_charact_pt;
    Real  boundaryPoint, internalBdPoint;

    switch ( _M_test_case ) {
    case 1: //! both eigenvalues are > 0
        // *******************************************************
        //! left BC.
        // *******************************************************
        rhsBC.first  = std::sin( 2 * M_PI * time_val ); //! Dirichlet bc.
        rhsBC.second = 1.; //! Dirichlet bc.

        //! solve the linear system Id * U(tn+1, z=z_right) = rhsBC (a bit artificial...)
        resBC = _solveLinearSyst2x2( line1_id, line2_id, rhsBC );

        _M_bcDirLeft = resBC;

        // *******************************************************
        //! right BC.
        // *******************************************************
        //-------------------------------------
        //!compatibility conditions

        boundaryPoint   = rightedge.pt2().x(); //!< point on the boundary (on the right of the edge!)
        internalBdPoint = rightedge.pt1().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( rightDof ) ,
                               _M_U2_thistime( rightDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( rightDof - 1 ) ,
                               _M_U2_thistime( rightDof - 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                rightDof);

        ASSERT( eigval1 > 0. && eigval2 > 0. ,
                "The eigenvalues do not have the expected signs.");

        //-------------------------------------
        //! first characteristics
        //! interpolation of U at the foot of the 1rst characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval1,
                                       U_boundary, U_internalBd);

        //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
        rhsBC.first = dot( left_eigvec1 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, rightDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, rightDof);

        //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
        rhsBC.first -= _M_time_step * dot( left_eigvec1 , qlSource );

        //-------------------------------------
        //! second characteristics
        //! interpolation of U at the foot of the 2nd characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval2,
                                       U_boundary, U_internalBd);

        //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
        rhsBC.second = dot( left_eigvec2 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, rightDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, rightDof);

        //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
        rhsBC.second -= _M_time_step * dot( left_eigvec2 , qlSource );


        /*!-------------------------------------
          solve the linear system L * U(tn+1, z=z_right) = rhsBC
          i.e.
          [ left_eigvec1^T;  [U1(tn+1,z_right);   = [ left_eigvec1^T * U(tn, charact_pt1);
          left_eigvec2^T ]  U2(tn+1,z_right) ]      left_eigvec2^T * U(tn, charact_pt2) ]

          - dt * [ left_eigvec1^T * qlSource(tn, charact_pt1);
          left_eigvec2^T * qlSource(tn, charact_pt2) ]


          with L such that L * dF/dU * L^{-1} = diag(eigval1, eigval2)
          -------------------------------------
        */
        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirRight = resBC;
        break;
    case 2: //! both eigenvalues are < 0
        // *******************************************************
        //! left BC.
        // *******************************************************
        //-------------------------------------
        //!compatibility conditions

        boundaryPoint   = leftedge.pt1().x(); //!< point on the boundary (on the left of the edge!)
        internalBdPoint = leftedge.pt2().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( leftDof ) ,
                               _M_U2_thistime( leftDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( leftDof + 1 ) ,
                               _M_U2_thistime( leftDof + 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                leftDof);


        ASSERT( eigval1 < 0. &&  eigval2 < 0.  ,
                "The eigenvalues do not have the expected signs.");


        //-------------------------------------
        //! first characteristics
        //! interpolation of U at the foot of the 1rst characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval1,
                                       U_boundary, U_internalBd);

        //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
        rhsBC.first = dot( left_eigvec1 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, leftDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, leftDof);

        //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
        rhsBC.first -= _M_time_step * dot( left_eigvec1 , qlSource );

        //-------------------------------------
        //! second characteristics
        //! interpolation of U at the foot of the 2nd characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval2,
                                       U_boundary, U_internalBd);

        //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
        rhsBC.second = dot( left_eigvec2 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, leftDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, leftDof);

        //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
        rhsBC.second -= _M_time_step * dot( left_eigvec2 , qlSource );

        /*!-------------------------------------
          solve the linear system L * U(tn+1, z=z_left) = rhsBC
          i.e.
          [ left_eigvec1^T;  [U1(tn+1,z_left);   = [ left_eigvec1^T * U(tn, charact_pt1);
          left_eigvec2^T ]  U2(tn+1,z_left) ]      left_eigvec2^T * U(tn, charact_pt2) ]


          with L such that L * dF/dU * L^{-1} = diag(eigval1, eigval2)
          -------------------------------------
        */
        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirLeft = resBC;


        // *******************************************************
        //! right BC.
        // *******************************************************
        rhsBC.first  = std::sin( 2 * M_PI * time_val ); //! Dirichlet bc.
        rhsBC.second = 1.; //! Dirichlet bc.

        //! solve the linear system Id * U(tn+1, z=z_right) = rhsBC (a bit artificial...)
        resBC = _solveLinearSyst2x2( line1_id, line2_id, rhsBC );

        _M_bcDirRight = resBC;
        break;

    case 3: //! eig1 > 0 and eig2 < 0.
        // *******************************************************
        //! left BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        // rhsBC.first = _M_oneDParam.Area0(0) * ( 1 + 0.1 * std::sin( 2 * M_PI * time_val ) );
        if( time_val <= 5e-2 ){
            rhsBC.first = _M_oneDParam.Area0( leftDof ) *
                ( 1 + 20000 * _M_oneDParam.DensityRho() / _M_oneDParam.Beta0( leftDof ) );
        } else {
            rhsBC.first = _M_oneDParam.Area0( leftDof );
        }

        //-------------------------------------
        //!compatibility condition

        boundaryPoint   = leftedge.pt1().x(); //!< point on the boundary (on the left of the edge!)
        internalBdPoint = leftedge.pt2().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( leftDof ) ,
                               _M_U2_thistime( leftDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( leftDof + 1 ) ,
                               _M_U2_thistime( leftDof + 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                leftDof);

        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        leftDof + 1);

        left_eigvec2_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval2,
        left_eigvec2, left_eigvec2_internalBd);

        std::cout << "\nEigenVector2 " << left_eigvec2.first << " " << left_eigvec2.second << std::endl;
        std::cout << "EigenVector2 charac " << left_eigvec2_charact_pt.first
        << " " << left_eigvec2_charact_pt.second << std::endl;

        left_eigvec2.first  = left_eigvec2_charact_pt.first;
        left_eigvec2.second = left_eigvec2_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        //------------------------------------------------------------

        if ( _M_verbose > 1 )
            std::cout << "\nEigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;

        ASSERT( eigval1 > 0. &&  eigval2 < 0.  ,
                "The eigenvalues do not have the expected signs.");


        //-------------------------------------
        //! second characteristics
        //! interpolation of U at the foot of the 2nd characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval2,
                                       U_boundary, U_internalBd);

        //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
        rhsBC.second = dot( left_eigvec2 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, leftDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, leftDof);

        //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
        rhsBC.second -= _M_time_step * dot( left_eigvec2 , qlSource );

        /*!-------------------------------------
          solve the linear system Theta_left * U(tn+1, z=z_left) = rhsBC
          i.e.
          [  1  0         ;  [U1(tn+1,z_left);   = [ dirichlet_value1;
          left_eigvec2^T ]  U2(tn+1,z_left) ]      left_eigvec2^T * U(tn, charact_pt2) ]
          -------------------------------------
        */
        resBC = _solveLinearSyst2x2( line1_id, left_eigvec2, rhsBC );

        _M_bcDirLeft = resBC;


        // *******************************************************
        //! right BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        rhsBC.second = 1.;

        //-------------------------------------
        //!compatibility condition
        boundaryPoint   = rightedge.pt2().x(); //!< point on the boundary (on the right of the edge!)
        internalBdPoint = rightedge.pt1().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( rightDof ) ,
                               _M_U2_thistime( rightDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( rightDof - 1 ) ,
                               _M_U2_thistime( rightDof - 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                rightDof);
        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        rightDof - 1);

        left_eigvec1_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval1,
        left_eigvec1, left_eigvec1_internalBd);

        std::cout << "\nEigenVector1 " << left_eigvec1.first << " " << left_eigvec1.second << std::endl;
        std::cout << "EigenVector1 charac " << left_eigvec1_charact_pt.first
        << " " << left_eigvec1_charact_pt.second << std::endl;

        left_eigvec1.first  = left_eigvec1_charact_pt.first;
        left_eigvec1.second = left_eigvec1_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        if ( _M_verbose > 1 )
            std::cout << "EigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;


        ASSERT( eigval1 > 0. && eigval2 < 0. ,
                "The eigenvalues do not have the expected signs.");

        //-------------------------------------
        //! first characteristics
        //! interpolation of U at the foot of the 1rst characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval1,
                                       U_boundary, U_internalBd);

        //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
        rhsBC.first = dot( left_eigvec1 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, rightDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, rightDof);

        //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
        rhsBC.first -= _M_time_step * dot( left_eigvec1 , qlSource );

        resBC = _solveLinearSyst2x2( left_eigvec1, line2_id, rhsBC );

        _M_bcDirRight = resBC;
        break;

    case 4: //! eig1 > 0 and eig2 < 0.
        /*! In this test case, we assume that the 2 values (U1 and U2) are given
          at each boundary. We use these 2 values to compute the incoming characteristics
          and set it as the boundary condition.
          The second boundary condition remains the compatibility condition.
          Once again, we treat the problem explicitly after linearisation.
        */
        // *******************************************************
        //! left BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        // rhsBC.first = _M_oneDParam.Area0(0) * ( 1 + 0.1 * std::sin( 2 * M_PI * time_val ) );
        if( time_val <= 5e-2 ){
            UBC_data.first = _M_oneDParam.Area0( leftDof ) *
                ( 1 + 20000 * _M_oneDParam.DensityRho() / _M_oneDParam.Beta0( leftDof ) );
        } else {
            UBC_data.first = _M_oneDParam.Area0( leftDof );
        }
        UBC_data.second = 0;


        //-------------------------------------
        //! Compute the eigen values and vectors
        boundaryPoint   = leftedge.pt1().x(); //!< point on the boundary (on the left of the edge!)
        internalBdPoint = leftedge.pt2().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( leftDof ) ,
                               _M_U2_thistime( leftDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( leftDof + 1 ) ,
                               _M_U2_thistime( leftDof + 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                leftDof);

        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        leftDof + 1);

        left_eigvec2_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval2,
        left_eigvec2, left_eigvec2_internalBd);

        std::cout << "\nEigenVector2 " << left_eigvec2.first << " " << left_eigvec2.second << std::endl;
        std::cout << "EigenVector2 charac " << left_eigvec2_charact_pt.first
        << " " << left_eigvec2_charact_pt.second << std::endl;

        left_eigvec2.first  = left_eigvec2_charact_pt.first;
        left_eigvec2.second = left_eigvec2_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        //------------------------------------------------------------

        if ( _M_verbose > 1 )
            std::cout << "\nEigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;

        ASSERT( eigval1 > 0. &&  eigval2 < 0.  ,
                "The eigenvalues do not have the expected signs.");


        //-------------------------------------
        //! Dirichlet bc. (treated through the incoming characteristics)
        //! rhsBC1 = left_eigvec1 dot UBC_data
        rhsBC.first = dot( left_eigvec1 , UBC_data );


        //-------------------------------------
        //!compatibility condition
        //-------------------------------------
        //! second characteristics
        //! interpolation of U at the foot of the 2nd characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval2,
                                       U_boundary, U_internalBd);

        //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
        rhsBC.second = dot( left_eigvec2 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, leftDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, leftDof);

        //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
        rhsBC.second -= _M_time_step * dot( left_eigvec2 , qlSource );

        /*!-------------------------------------
          solve the linear system Theta_left * U(tn+1, z=z_left) = rhsBC
          i.e.
          [ left_eigvec1^T ;  [U1(tn+1,z_left);   = [ left_eigvec1^T * UBC_data;
          left_eigvec2^T ]   U2(tn+1,z_left) ]      left_eigvec2^T * U(tn, charact_pt2) - deltaT * left_eigvec2^T * qlSource(tn, charact_pt2)]
          -------------------------------------
        */
        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirLeft = resBC;

        // *******************************************************
        //! right BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        UBC_data.first = _M_oneDParam.Area0( leftDof );
        UBC_data.second = 0.;

        //-------------------------------------
        //! Compute the eigen values and vectors
        boundaryPoint   = rightedge.pt2().x(); //!< point on the boundary (on the right of the edge!)
        internalBdPoint = rightedge.pt1().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( rightDof ) ,
                               _M_U2_thistime( rightDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( rightDof - 1 ) ,
                               _M_U2_thistime( rightDof - 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                rightDof);
        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        rightDof - 1);

        left_eigvec1_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval1,
        left_eigvec1, left_eigvec1_internalBd);

        std::cout << "\nEigenVector1 " << left_eigvec1.first << " " << left_eigvec1.second << std::endl;
        std::cout << "EigenVector1 charac " << left_eigvec1_charact_pt.first
        << " " << left_eigvec1_charact_pt.second << std::endl;

        left_eigvec1.first  = left_eigvec1_charact_pt.first;
        left_eigvec1.second = left_eigvec1_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        if ( _M_verbose > 1 )
            std::cout << "EigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;


        ASSERT( eigval1 > 0. && eigval2 < 0. ,
                "The eigenvalues do not have the expected signs.");

        //-------------------------------------
        //! Dirichlet bc. (treated through the incoming characteristics)
        //! rhsBC2 = left_eigvec2 dot UBC_data
        rhsBC.second = dot( left_eigvec2 , UBC_data );

        //-------------------------------------
        //!compatibility condition
        //-------------------------------------
        //! first characteristics
        //! interpolation of U at the foot of the 1rst characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval1,
                                       U_boundary, U_internalBd);

        //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
        rhsBC.first = dot( left_eigvec1 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, rightDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, rightDof);

        //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
        rhsBC.first -= _M_time_step * dot( left_eigvec1 , qlSource );

        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirRight = resBC;
        break;


    case 999: //! eig1 > 0 and eig2 < 0. ONLY LEFT!!!!!!!!!!!
        /*! In this test case, we assume that the 2 values (U1 and U2) are given
          at each boundary. We use these 2 values to compute the incoming characteristics
          and set it as the boundary condition.
          The second boundary condition remains the compatibility condition.
          Once again, we treat the problem explicitly after linearisation.
        */
        // *******************************************************
        //! left BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        // rhsBC.first = _M_oneDParam.Area0(0) * ( 1 + 0.1 * std::sin( 2 * M_PI * time_val ) );
        if( time_val <= 5e-2 ){
            UBC_data.first = _M_oneDParam.Area0( leftDof ) *
                ( 1 + 20000 * _M_oneDParam.DensityRho() / _M_oneDParam.Beta0( leftDof ) );
        } else {
            UBC_data.first = _M_oneDParam.Area0( leftDof );
        }
        UBC_data.second = 0;


        //-------------------------------------
        //! Compute the eigen values and vectors
        boundaryPoint   = leftedge.pt1().x(); //!< point on the boundary (on the left of the edge!)
        internalBdPoint = leftedge.pt2().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( leftDof ) ,
                               _M_U2_thistime( leftDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( leftDof + 1 ) ,
                               _M_U2_thistime( leftDof + 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                leftDof);

        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        leftDof + 1);

        left_eigvec2_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval2,
        left_eigvec2, left_eigvec2_internalBd);

        std::cout << "\nEigenVector2 " << left_eigvec2.first << " " << left_eigvec2.second << std::endl;
        std::cout << "EigenVector2 charac " << left_eigvec2_charact_pt.first
        << " " << left_eigvec2_charact_pt.second << std::endl;

        left_eigvec2.first  = left_eigvec2_charact_pt.first;
        left_eigvec2.second = left_eigvec2_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        //------------------------------------------------------------

        if ( _M_verbose > 1 )
            std::cout << "\nEigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;

        ASSERT( eigval1 > 0. &&  eigval2 < 0.  ,
                "The eigenvalues do not have the expected signs.");


        //-------------------------------------
        //! Dirichlet bc. (treated through the incoming characteristics)
        //! rhsBC1 = left_eigvec1 dot UBC_data
        rhsBC.first = dot( left_eigvec1 , UBC_data );


        //-------------------------------------
        //!compatibility condition
        //-------------------------------------
        //! second characteristics
        //! interpolation of U at the foot of the 2nd characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval2,
                                       U_boundary, U_internalBd);

        //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
        rhsBC.second = dot( left_eigvec2 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, leftDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, leftDof);

        //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
        rhsBC.second -= _M_time_step * dot( left_eigvec2 , qlSource );

        /*!-------------------------------------
          solve the linear system Theta_left * U(tn+1, z=z_left) = rhsBC
          i.e.
          [ left_eigvec1^T ;  [U1(tn+1,z_left);   = [ left_eigvec1^T * UBC_data;
          left_eigvec2^T ]   U2(tn+1,z_left) ]      left_eigvec2^T * U(tn, charact_pt2) - deltaT * left_eigvec2^T * qlSource(tn, charact_pt2)]
          -------------------------------------
        */
        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirLeft = resBC;

        break;

    case 9999: //! eig1 > 0 and eig2 < 0. ONLY RIGHT!!!!!!!!!!!
        /*! In this test case, we assume that the 2 values (U1 and U2) are given
          at each boundary. We use these 2 values to compute the incoming characteristics
          and set it as the boundary condition.
          The second boundary condition remains the compatibility condition.
          Once again, we treat the problem explicitly after linearisation.
        */

        // *******************************************************
        //! right BC.
        // *******************************************************
        //-------------------------------------
        //! Dirichlet bc.
        UBC_data.first = _M_oneDParam.Area0( leftDof );
        UBC_data.second = 0.;

        //-------------------------------------
        //! Compute the eigen values and vectors
        boundaryPoint   = rightedge.pt2().x(); //!< point on the boundary (on the right of the edge!)
        internalBdPoint = rightedge.pt1().x(); //!< neighboring point (internal)

        //! values of U on the boundary
        U_boundary   = Vec2D ( _M_U1_thistime( rightDof ) ,
                               _M_U2_thistime( rightDof ) );
        //! values of U on the neighboring node of the boundary point
        U_internalBd = Vec2D ( _M_U1_thistime( rightDof - 1 ) ,
                               _M_U2_thistime( rightDof - 1 ) );

        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                eigval1, eigval2,
                                                left_eigvec1.first, left_eigvec1.second,
                                                left_eigvec2.first, left_eigvec2.second,
                                                rightDof);
        /*! interpolate the eigenvectors
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        //! compute the eigenvalues/eigenvectors of the flux jacobian (computed at the internal bd point)
        _M_fluxFun.jacobian_EigenValues_Vectors(U_internalBd.first, U_internalBd.second,
        eigval1_internalBd, eigval2_internalBd,
        left_eigvec1_internalBd.first, left_eigvec1_internalBd.second,
        left_eigvec2_internalBd.first, left_eigvec2_internalBd.second,
        rightDof - 1);

        left_eigvec1_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
        _M_time_step, eigval1,
        left_eigvec1, left_eigvec1_internalBd);

        std::cout << "\nEigenVector1 " << left_eigvec1.first << " " << left_eigvec1.second << std::endl;
        std::cout << "EigenVector1 charac " << left_eigvec1_charact_pt.first
        << " " << left_eigvec1_charact_pt.second << std::endl;

        left_eigvec1.first  = left_eigvec1_charact_pt.first;
        left_eigvec1.second = left_eigvec1_charact_pt.second;
        //wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
        */

        if ( _M_verbose > 1 )
            std::cout << "EigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;


        ASSERT( eigval1 > 0. && eigval2 < 0. ,
                "The eigenvalues do not have the expected signs.");

        //-------------------------------------
        //! Dirichlet bc. (treated through the incoming characteristics)
        //! rhsBC2 = left_eigvec2 dot UBC_data
        rhsBC.second = dot( left_eigvec2 , UBC_data );

        //-------------------------------------
        //!compatibility condition
        //-------------------------------------
        //! first characteristics
        //! interpolation of U at the foot of the 1rst characteristics
        U_charact_pt = _interpolLinear(boundaryPoint, internalBdPoint,
                                       _M_time_step, eigval1,
                                       U_boundary, U_internalBd);

        //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
        rhsBC.first = dot( left_eigvec1 , U_charact_pt );

        //! take into account the (quasi linear) source term
        /*!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
          THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
        */
        qlSource.first  = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         1, rightDof);
        qlSource.second = _M_sourceFun.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                         2, rightDof);

        //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
        rhsBC.first -= _M_time_step * dot( left_eigvec1 , qlSource );

        resBC = _solveLinearSyst2x2( left_eigvec1, left_eigvec2, rhsBC );

        _M_bcDirRight = resBC;
        break;

    default:
        ERROR_MSG("No such test case.");
    }
}

//! linear interpolation at the foot of the characteristic
//! determined by the given eigenvalue
OneDModelSolver::Vec2D
OneDModelSolver::_interpolLinear(const Real& point_bound, const Real& point_internal,
                                 const Real& deltaT, const Real& eigenvalue,
                                 const Vec2D& U_bound, const Vec2D& U_intern) const
{
    Real deltaX = std::abs(point_bound - point_internal);

    Real cfl =  eigenvalue * deltaT / deltaX;

    Real weight;   //!< weight in the linear approximation

    if ( point_bound < point_internal ) { //! the edge is on the left of the domain
        ASSERT( -1. < cfl && cfl < 0. ,
                "This characteristics is wrong!\nEither it is not outcoming (eigenvalue>0 at the left of the domain),\n or CFL is too high.");

        weight = - cfl;
    }
    else {  //! the edge is on the right of the domain
        ASSERT( 0. < cfl && cfl < 1. ,
                "This characteristics is wrong!\nEither it is not outcoming (eigenvalue<0 at the right of the domain),\n or CFL is too high.");

        weight = cfl;
    }

    Vec2D u_interp( ( 1 - weight ) * U_bound.first  + weight * U_intern.first ,
                    ( 1 - weight ) * U_bound.second + weight * U_intern.second );
    return u_interp;

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

//! simple cfl computation (correct for constant mesh)
void OneDModelSolver::CheckCFL() const
{
    Real CFL = 0.;

    //! length of the first edge (arbitrary as they are all supposed equal).
    Real deltaX = _M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x();
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

    CFL = _M_time_step / deltaX * std::max<Real>( lambda1_max , lambda2_max );
    if ( _M_verbose > 1) {
        std::cout << "CFL = " << CFL << std::endl;
    }

    ASSERT( CFL < 0.5774 , "CFL not respected" );
}

/*! solve a 2x2 linear system by the Cramer method (for the boundary systems)
  (beware of ill-conditioning!...).
  Matrix A is given by two pairs corresponding to the 2 lines.
  A = [line1;
  line2 ]
  return A^{-1} * rhs2d
*/
OneDModelSolver::Vec2D
OneDModelSolver::_solveLinearSyst2x2(const Vec2D& line1, const Vec2D& line2,
                                     const Vec2D& rhs2d) const
{
    Real aa11, aa12, aa21, aa22;

    aa11 = line1.first;  aa12 = line1.second;
    aa21 = line2.first;  aa22 = line2.second;

    Real determinant = aa11 * aa22 - aa12 * aa21;

    ASSERT( determinant != 0,
            "Error: the 2x2 system on the boundary is not invertible. \nCheck the boundary conditions.");
    Vec2D res( (   aa22 * rhs2d.first - aa12 * rhs2d.second ) / determinant,
               ( - aa21 * rhs2d.first + aa11 * rhs2d.second ) / determinant );

    return res;
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

//! 2D dot product
Real OneDModelSolver::dot(const Vec2D& vec1, const Vec2D& vec2) const
{
    return vec1.first * vec2.first + vec1.second * vec2.second;
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

// ! Initialize with constant initial conditions concentration
void OneDModelSolver::initialize(const Real& u10, const Real& u20)
{
    _M_U1_thistime = ScalarVector( _M_U1_thistime.size(), u10 );
    _M_U2_thistime = ScalarVector( _M_U1_thistime.size(), u20 );
}

// ! Initialize when  initial conditions concentration
void OneDModelSolver::initialize(const Function& /* c0 */,
                                 Real /* t0 */,
                                 Real /* dt */)
{
    ERROR_MSG("Not yet implemented");
}

// ! Initialize when initial values for the concentration are read from file
void OneDModelSolver::initialize(const std::string & vname)
{
    ERROR_MSG("Not yet implemented");

    std::fstream Resfile(vname.c_str(),std::ios::in | std::ios::binary);
    if (Resfile.fail()) {
        std::cerr<<" Error in initialize: File not found or locked"<<std::endl;
        abort();
    }
    Resfile.read((char*)&_M_U1_thistime(0),_M_U1_thistime.size()*sizeof(Real));
    Resfile.close();
}


//! Update the right hand side for time advancing
void OneDModelSolver::timeAdvance( const Real& time_val )
{
    Chrono chrono;

    Real dt2over2 = _M_time_step * _M_time_step * 0.5;

    std::cout << "  o-  updates of flux and sources... ";
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
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    std::cout << "  o-  updates of matrices... ";
    chrono.start();
    //! update the matrices for the non-linear terms
    _updateMatrices();
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

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

    std::cout << "  o-  Matrix vector products... ";
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
    _computeBCValues( time_val );
    //! take into account the bc
    _updateBCDirichletVector();

    // *******************************************************
    chrono.stop();
    std::cout << "rhs computed in " << chrono.diff() << " s." << std::endl;

}

void OneDModelSolver::iterate( const Real& time_val , const int& count)
{
    std::cout << "  o-  Solving the system... t = " << time_val
              << ", iter = " << count  << "... ";
    Chrono chrono;
    chrono.start();

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs1 );

    //! solve the system: rhs2 = massFactor^{-1} * rhs2
    _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs2 );

    //! update the solution for the next time step
    _M_U1_thistime = _M_rhs1;
    _M_U2_thistime = _M_rhs2;

    if( !(count % 5) ){
        std::string fname1 = _M_post_dir + "/" + _M_post_file + "A.mtv";
        std::string fname2 = _M_post_dir + "/" + _M_post_file + "Q.mtv";
        output_to_plotmtv( fname1, time_val, _M_mesh.pointList(), _M_U1_thistime, count);
        output_to_plotmtv( fname2, time_val, _M_mesh.pointList(), _M_U2_thistime, count);
    }
    // *******************************************************
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}


void OneDModelSolver::gplot( )
{
    _M_GracePlot.Plot( _M_mesh.pointList(), _M_U1_thistime );
}

//! output for Plotmtv.
void OneDModelSolver::output_to_plotmtv(std::string fname, Real time_val,
                                        const std::vector< Point1D >& ptlist,
                                        const ScalUnknown<Vector>& U,
                                        const int& count )
{

    FILE *fp;
    if (count==0){
        fp = fopen(fname.c_str(),"w");
    } else {
        fp = fopen(fname.c_str(),"a");
    }
    fprintf(fp,"$ DATA = CURVE2D\n %% xlabel='z'\n");
    fprintf(fp,"%% toplabel='Section,time=%f'\n %% ylabel='A'\n",time_val);

    for(UInt ii = 0; ii < U.size(); ii++){
        fprintf(fp,"%10.6f %10.6f\n", ptlist[ii].x(), U[ii]);
    }
    fclose(fp);

}

}
