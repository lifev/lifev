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
  \file oneDModelSolver.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a solver class for the 1D model.

  ---------------------------------------------------
  1 dimensional hyperbolic equation:
  ---------------------------------------------------
      dU/dt  +  dF(U)/dz  +  S(U) = 0
  ---------------------------------------------------

  with U = [U1,U2]^T in R^2.

  The non linear flux function F(U) and source function S(U)
  are quite independant of this solver : they are taken into 
  account only via two classes that define a vectorial function
  and its derivatives.

  More precisely:
  two functions (_M_fluxFun and _M_sourceFun) have to be defined
  separately to allow the update (_updateFluxDer, _updateSourceDer) 
  of the corresponding vectors (_M_Fluxi, _M_diffFluxij,
  _M_Sourcei, _M_diffSrcij). I also separated the treatment of the
  parameters that exist in the functions.
  Normally, one should be able to create a template to allow
  the user to select between different problems (linear, non-linear,
  etc.). 


    -----------------------------------------------------
    Solver based on the 2nd order Taylor-Galerkin scheme.
    -----------------------------------------------------
    (see for instance Formaggia-Veneziani, Mox report no 21 june 2003)

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

    Approximation of the unknowns and non-linearities:

    Let's define:
    a) (phi_i)_{i in nodes} is the basis of P1 (the "hat" functions)
    b) (1_{i+1/2})_{i+1/2 in elements} is the basis of P0 (constant per
        element). The vertices of the element "i+1/2" are the nodes "i" and "i+1".

   Then:

    c) Uh    is in P1 : U = sum_{i in nodes} U_i phi_i

    d) Fh(U) is in P1 : F(U) = sum_{i in nodes} F(U_i) phi_i 
    e) diffFh(U) is in P0 : 
diffFlux(U) = sum_{i+1/2 in elements} 1/2 { dF/dU(U_i) + dF/dU(U_i+1) } 1_{i+1/2}
    (means of the two extremal values of the cell)

    We note that d) allows us to define easily
    f) dF/dz(U) = sum_{i in nodes} F(U_i) d(phi_i)/dz 

    g) Sh(U) is in P1 : S(U) = sum_{i in nodes} S(U_i) phi_i
    h) diffSh(U) is in P0 : 
diffSrc(U) = sum_{i+1/2 in elements} 1/2 { dS/dU(U_i) + dS/dU(U_i+1) } 1_{i+1/2}
    (means of the two extremal values of the cell)


    -----------------------------------------------------
    The option taken here is to define the different tridiagonal matrix
    operators (div, grad, mass, stiff) and reconstruct them at each time
    step (as they depend on diffFlux and diffSrc). They are thus rebuilt
    at the element level and reassembled.
    Afterwards, there remains to do only some tridiagonal matrix vector
    products to obtain the right hand side.

    This procedure might appear a bit memory consuming (there are 18
    tridiagonal matrices stored), but it has the advantage of being
    very clear. If it is too costly, it should be quite easy to improve 
    it.

    -----------------------------------------------------
*/

#ifndef _ONEDMODELSOLVER_H_
#define _ONEDMODELSOLVER_H_

#include <string>

#include <oneDModelHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>

#include <life/lifefem/assemb.hpp>
#include <life/lifecore/chrono.hpp>

#include <life/lifearray/tridiagMatrix.hpp>
#include <triDiagCholesky.hpp>
#include <triDiagLU.hpp>
#include <oneDNonLinModelParam.hpp>
#include <vectorFunction1D.hpp>


namespace LifeV
{
/*!
  \class OneDModelSolver

  This class contains a solver class for the 1D model.

*/
class OneDModelSolver:
        public OneDModelHandler
{

public:

    //! Constructor
    /*!
      \param data_file GetPot data file
    */
    OneDModelSolver(const GetPot& data_file,
                    const OneDNonLinModelParam& onedparam);
    // const LinearSimpleParam& onedparam);

    //! return the solution at current time step
    const ScalUnknown<Vector>& U1_thistime() const { return _M_U1_thistime;}

    //! return the solution at current time step
    const ScalUnknown<Vector>& U2_thistime() const { return _M_U2_thistime;}

    //! Sets initial condition for the concentration
    void initialize(const Real& u10, const Real& u20);

    //! Sets initial condition for the concentration
    //! (incremental approach): the initial time is t0, the time step dt
    void initialize(const Function& c0, Real t0, Real dt);

    //! Sets initial condition for the concentration from file
    void initialize(const std::string & vname);

    //! Update the right  hand side  for time advancing
    void timeAdvance( const Real& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate( const Real& time , const int& count);

    //! get the Dirichlet boundary conditions (left)
    Vec2D BCValuesLeft() const;
    //! get the value at neighboring node (left)
    Vec2D BCValuesInternalLeft() const;
    //! get the Dirichlet boundary conditions (right)
    Vec2D BCValuesRight() const;
    //! get the value at neighboring node (right)
    Vec2D BCValuesInternalRight() const;

    //! set the Dirichlet boundary conditions (left)
    void setBCValuesLeft( const Real& bcL1, const Real& bcL2 );
    //! set the Dirichlet boundary conditions (right)
    void setBCValuesRight( const Real& bcR1, const Real& bcR2 );

    //! get the flux function
    NonLinearFluxFun1D const& FluxFun() const;
    //! get the source function
    NonLinearSourceFun1D const& SourceFun() const;

    //! get the left edge
    Edge1D LeftEdge() const;
    //! get the right edge
    Edge1D RightEdge() const;

    //! get the left node
    UInt LeftNodeId() const;
    //! get the left internal node (neighboring node)
    UInt LeftInternalNodeId() const;
    //! get the right node
    UInt RightNodeId() const;
    //! get the right internal node (neighboring node)
    UInt RightInternalNodeId() const;

    //! simple cfl computation (correct for constant mesh)
    void CheckCFL() const;

    //! plotting
    void gplot();

    //! writer
    void output_to_plotmtv(std::string fname, Real time_val, 
                           const std::vector< Point1D >& ptlist, 
                           const ScalUnknown<Vector>& U,
                           const int& count);


private:

    //!@{  TO BE Templatized !!
    //! the parameters
    const OneDNonLinModelParam& _M_oneDParam;
    // const LinearSimpleParam& _M_oneDParam;

    //! the flux function
    NonLinearFluxFun1D _M_fluxFun;
    //! the source function 
    NonLinearSourceFun1D _M_sourceFun ;
    /*
    //! the flux function
    LinearSimpleFluxFun1D _M_fluxFun;
    //! the source function
    LinearSimpleSourceFun1D _M_sourceFun ;
    */
    //!@} end of TO BE Templatized !!

    const UInt _M_leftNodeId;
    const UInt _M_leftInternalNodeId;
    const UInt _M_rightNodeId;
    const UInt _M_rightInternalNodeId;

    //! boundary edges
    const Edge1D _M_leftEdge;
    const Edge1D _M_rightEdge;

    //! coefficient in front of the corresponding _M_elmat*
    Real _M_coeffMass;
    Real _M_coeffStiff;
    Real _M_coeffGrad;
    Real _M_coeffDiv;


    ElemMat _M_elmatMass;  //!< element mass matrix
    ElemMat _M_elmatStiff; //!< element stiffness matrix
    ElemMat _M_elmatGrad;  //!< element gradient matrix
    ElemMat _M_elmatDiv;   //!< element divergence matrix

    //  ElemVec _M_elvec; // Elementary right hand side

    //! Unknown at present time step
    ScalUnknown<Vector> _M_U1_thistime;
    ScalUnknown<Vector> _M_U2_thistime;

    //! Right hand sides of the linear system i: "mass * _M_Ui = _M_rhsi"
    ScalUnknown<Vector> _M_rhs1;
    ScalUnknown<Vector> _M_rhs2;

    //! Flux F(U) (in P1)
    ScalUnknown<Vector> _M_Flux1;
    ScalUnknown<Vector> _M_Flux2;
    //! diffFlux = dF(U)/dU (in P0) 
    ScalUnknown<Vector> _M_diffFlux11;
    ScalUnknown<Vector> _M_diffFlux12;
    ScalUnknown<Vector> _M_diffFlux21;
    ScalUnknown<Vector> _M_diffFlux22;

    //! Source term S (in P1)
    ScalUnknown<Vector> _M_Source1;
    ScalUnknown<Vector> _M_Source2;
    //! diffSrc = dSource(U)/dU (in P0)
    ScalUnknown<Vector> _M_diffSrc11;
    ScalUnknown<Vector> _M_diffSrc12;
    ScalUnknown<Vector> _M_diffSrc21;
    ScalUnknown<Vector> _M_diffSrc22;



    //! tridiagonal mass matrix
    TriDiagMatrix<Real> _M_massMatrix;

    //! factorized tridiagonal mass matrix
    TriDiagMatrix<Real> _M_factorMassMatrix;

    //!@{  TO BE Templatized !!
    //! cholesky factorization
    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _M_tridiagSlv;
    //! lapack LU factorization
    // TriDiagLU< Real, TriDiagMatrix<Real>, Vector > _M_tridiagSlv;
    //!@} end of  TO BE Templatized !!
   

    //! tridiagonal mass matrices multiplied by diffSrcij
    TriDiagMatrix<Real> _M_massMatrixDiffSrc11;
    TriDiagMatrix<Real> _M_massMatrixDiffSrc12;
    TriDiagMatrix<Real> _M_massMatrixDiffSrc21;
    TriDiagMatrix<Real> _M_massMatrixDiffSrc22;

    //! tridiagonal stiffness matrices multiplied by diffFluxij
    TriDiagMatrix<Real> _M_stiffMatrixDiffFlux11;
    TriDiagMatrix<Real> _M_stiffMatrixDiffFlux12;
    TriDiagMatrix<Real> _M_stiffMatrixDiffFlux21;
    TriDiagMatrix<Real> _M_stiffMatrixDiffFlux22;

    //! tridiagonal gradient matrix
    TriDiagMatrix<Real> _M_gradMatrix;
    //! tridiagonal gradient matrices multiplied by diffFluxij
    TriDiagMatrix<Real> _M_gradMatrixDiffFlux11;
    TriDiagMatrix<Real> _M_gradMatrixDiffFlux12;
    TriDiagMatrix<Real> _M_gradMatrixDiffFlux21;
    TriDiagMatrix<Real> _M_gradMatrixDiffFlux22;

    //! tridiagonal divergence matrices multiplied by diffSrcij
    TriDiagMatrix<Real> _M_divMatrixDiffSrc11;
    TriDiagMatrix<Real> _M_divMatrixDiffSrc12;
    TriDiagMatrix<Real> _M_divMatrixDiffSrc21;
    TriDiagMatrix<Real> _M_divMatrixDiffSrc22;


    //! Update the coefficients
    //! (from the flux, source functions and their derivatives)
    void _updateMatrixCoefficients(const UInt& ii, const UInt& jj,
                                   const UInt& iedge);

    //! Update the element matrices with the current element
    void _updateElemMatrices();

    //! assemble the matrices
    int _assemble_matrices(const UInt& ii, const UInt& jj );

    /*! update the matrices
      _M_massMatrixDiffSrcij, _M_stiffMatrixDiffFluxij
      _M_gradMatrixDiffFluxij, and _M_divMatrixDiffSrcij (i,j=1,2)

      from the values of diffFlux(Un) and diffSrc(Un)
      that are computed with _updateMatrixCoefficients.

      call of  _updateMatrixCoefficients,
      _updateElemMatrices and _assemble_matrices.
    */
    void _updateMatrices();

    /*! modify the matrix to take into account
      the Dirichlet boundary conditions
      (works for P1Seg and canonic numbering!)
    */
    void _updateBCDirichletMatrix( TriDiagMatrix<Real>& mat );

    /*! modify the vector to take into account
      the Dirichlet boundary conditions
      (works for P1Seg and canonic numbering!)
    */
    void _updateBCDirichletVector();

    /*! compute the _M_bcDirLeft and _M_bcDirRight
      from the external boundary condition
      and the compatibility condition.

      use the extrapolation of the pseudo-characteristics 

      (from the solution at time n, obtain the 
      value of the linearized characteristics
      at time n+1.)
      Used as Compatibility condition.
    */
    void _computeBCValues( const Real& time_val );

    //! linear interpolation at the foot of the characteristic
    //! determined by the given eigenvalue
    Vec2D _interpolLinear(const Real& point_bound, const Real& point_internal,
                          const Real& deltaT, const Real& eigenvalue, 
                          const Vec2D& U_bound, const Vec2D& U_intern) const;

    //! solve a 2x2 linear system by the Cramer method (for the boundary systems)
    //! (beware of ill-conditioning!...).
    Vec2D _solveLinearSyst2x2(const Vec2D& line1, const Vec2D& line2,
                              const Vec2D& rhs2d) const;

    //! Axpy product for 2D vectors (pairs)
    //! Axpy(alpha, x, beta, y) -> y = a*A*x + beta*y
    void Axpy(const Vec2D& line1, const Vec2D& line2,
              const Real& alpha,  const Vec2D& x,
              const Real& beta,   Vec2D& y) const;

    //! 2D dot product
    Real dot(const Vec2D& vec1, const Vec2D& vec2) const;

    //! update the P1 flux vector from U: _M_Fluxi = F_h(Un) i=1,2
    void _updateFlux();
    //! update the P1 source vector from U: _M_Sourcei = S_h(Un) i=1,2
    void _updateSource();

    //! call _updateFlux and update the P0 derivative of flux vector from U: 
    //! _M_diffFluxij = dF_h/dU(Un) i,j=1,2
    void _updateFluxDer();
    //! call _updateSource and update the P0 derivative of source vector from U: 
    //! _M_diffSrcij = dS_h/dU(Un) i,j=1,2
    void _updateSourceDer();

};
}

#endif
