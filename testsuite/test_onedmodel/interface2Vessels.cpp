/* -*- mode: c++ -*-

This file is part of the LifeV library

Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
Date: 2004-11-02

Copyright (C) 2004 EPFL, INRIA, Politecnico di Milano

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
   \file interface2Vessels.cpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-11-02

   \brief This file implements a Taylor-Galerkin solver for 1D model.
*/
#include <life/lifealg/clapack.hpp>

#include "interface2Vessels.hpp"


namespace LifeV
{

Interface2Vessels::Interface2Vessels( OneDModelSolver const& tube_alpha,
                                      OneDModelSolver const& tube_beta ) :
    _M_Un_alpha_bd( tube_alpha.BCValuesRight() ),
    _M_Un_alpha_int( tube_alpha.BCValuesInternalRight() ),
    _M_edge_alpha( tube_alpha.RightEdge() ),
    _M_dof_alpha( tube_alpha.RightNodeId() ),
    _M_fluxFun_alpha( tube_alpha.FluxFun() ),
    _M_sourceFun_alpha( tube_alpha.SourceFun() ),
    _M_bcDir_alpha( 2 ),
    _M_Un_beta_bd( tube_beta.BCValuesLeft() ),
    _M_Un_beta_int( tube_beta.BCValuesInternalLeft() ),
    _M_edge_beta( tube_beta.LeftEdge() ),
    _M_dof_beta( tube_beta.LeftNodeId() ),
    _M_fluxFun_beta( tube_beta.FluxFun() ),
    _M_sourceFun_beta( tube_beta.SourceFun() ),
    _M_bcDir_beta( 2 ),
    _M_time_step( tube_alpha.timestep() )
{
}

void Interface2Vessels::
updateInterface2Vessels( OneDModelSolver const& tube_alpha,
                         OneDModelSolver const& tube_beta )
{
    _M_Un_alpha_bd  = tube_alpha.BCValuesRight();
    _M_Un_alpha_int = tube_alpha.BCValuesInternalRight();
    _M_Un_beta_bd   = tube_beta.BCValuesLeft();
    _M_Un_beta_int  = tube_beta.BCValuesInternalLeft();
}

Vector Interface2Vessels::BcDir_alpha() const
{
    return _M_bcDir_alpha;
}
Vector Interface2Vessels::BcDir_beta() const
{
    return _M_bcDir_beta;
}

int Interface2Vessels::computeInterface2TubesValues()
{


    //----------------------------------------------
    const UInt f_size(4);

    //! unknown of non linear equation f(x) = 0
    Vector x(f_size);
    //! non linear function f
    Vector f(f_size);
    //! jacobian of the non linear function
    Matrix jac(f_size,f_size);
    //! transpose of the jacobian of the non linear function
    Matrix jac_trans(f_size,f_size);

    //! tmp matrix for lapack lu inversion
    boost::numeric::ublas::vector<Int> ipiv(f_size);

    x[0] = _M_Un_alpha_bd.first; //!< A_alpha
    x[1] = _M_Un_alpha_bd.second;//!< Q_alpha
    x[2] = _M_Un_beta_bd.first;  //!< A_beta
    x[3] = _M_Un_beta_bd.second; //!< Q_beta

    //! lapack variable
    int INFO[1] = {0};
    int NBRHS[1] = {1};//  nb columns of the rhs := 1.
    int NBU[1] = {f_size};

    //! newton raphson iteration
    //! (use the newton class???)
    for ( UInt iter = 0 ; iter < 10 ; iter ++) {

        // std::cout << "=======\n\tNewton iter = " << iter << std::endl;

        //! compute f(x) and its jacobian df(x)
        f_jac( x, f, jac);

    //  std::cout << "---After call of f_jac:\nx : " << x << "\nf : " << f << "\njac : " << jac << std::endl;

        //! transpose to pass to fortran storage (lapack!)
        jac_trans = trans(jac);

        //! Compute f <-  ( df(x)^{-1} f(x) ) (lu dcmp)
        dgesv_(NBU, NBRHS, &jac_trans(0,0), NBU , &ipiv(0), &f(0), NBU, INFO);
        ASSERT_PRE(!INFO[0],"Lapack LU resolution of y = df(x)^{-1} f(x) is not achieved.");
        x += - f;

        // std::cout << "---After lapack inversion:\nx : " << x << "\ndf(x)^{-1}f(x) : " << f << std::endl;

    /* //Write a correct test here!
        //! convergence if Q_alpha == Q_beta
        if ( std::fabs( x[1] - x[3] ) < 1e-12 ) {
            _M_bcDir_alpha( 0 ) = x[0];
            _M_bcDir_alpha( 1 ) = x[1];
            _M_bcDir_beta( 0 )  = x[2];
            _M_bcDir_beta( 1 )  = x[3];

            std::cout << "\n\tNewton finished : iter =" << iter << "\n========" << std::endl;
            return 0;
        }
    */
    }
    //! dummy convergence test ( if Q_alpha == Q_beta )
    if ( std::fabs( x[1] - x[3] ) < 1e-12 ) {
      _M_bcDir_alpha( 0 ) = x[0];
      _M_bcDir_alpha( 1 ) = x[1];
      _M_bcDir_beta( 0 )  = x[2];
      _M_bcDir_beta( 1 )  = x[3];

      // std::cout << "\n\tNewton finished\n========" << std::endl;
      return 0;
    }
    //! no convergence
    return 1;
}

/*!
    Function and its gradient
    f : R4 -> R4
    x=[A_alpha, Q_alpha, A_beta, Q_beta]

    @ Continuity of the flux
    f1(x) = x[3] - x[1] = Q_beta - Q_alpha
    @ Continuity of the total pressure
    f2(x) = Pt_beta(x[2],x[3]) - Pt_alpha(x[0],x[1])
      with Pt(A,Q) = P(A) + rho/2 * (Q/A)^2
                   = beta0 ( (A/A0)^beta1 - 1 ) + rho/2 * (Q/A)^2
    @ Compatibility condition for tube alpha (explicit and linearized)
    f3(x) = left_eigenvector1_alpha dot {   (A_alpha, Q_alpha)
                                          - (U_alpha + dt G_alpha)_extrapol }
    @ Compatibility condition for tube beta (explicit and linearized)
    f4(x) = left_eigenvector2_beta dot {   (A_beta, Q_beta)
                                          - (U_beta + dt G_beta)_extrapol }

    @@ jacobian of f in x:
    jac(x) = [ [ 0 , -1 , 0 , 1 ] ;
           [ - dP_alpha/dA_alpha + rho Q_alpha^2/A_alpha^3, - rho Q_alpha/A_alpha,
               dP_beta/dA_beta   - rho Q_beta^2/A_beta^3  ,   rho Q_beta/A_beta   ] ;
           [ left_eigenvector1_alpha_1 , left_eigenvector1_alpha_2, 0 , 0 ] ;
           [ 0 , 0, left_eigenvector2_beta_1 , left_eigenvector2_beta_2 ]   ]

 */
void Interface2Vessels::
f_jac( const Vector& x, Vector& f, Matrix& jac ) const
{
    //! eigen values of the jacobian diffFlux (= dF/dU)
    Real  eigval1, eigval2;
    //! left eigen vectors for the eigen values eigval1 and eigval2
    Vec2D left_eigvec1, left_eigvec2;

    //! right hand side for the 2x2 linear system to be solved at each side
    Real rhsBC1, rhsBC2;

    //! quasi linear source term
    Vec2D qlSource;

    //! value of U at the boundary, at the neighboring internal node
    //!            and at the foot of the characteristics
    Vec2D U_boundary, U_internalBd, U_charact_pt;
    Real  boundaryPoint, internalBdPoint;

    int verbose = 2;

    Real A_alpha = x[0];
    Real Q_alpha = x[1];
    Real A_beta  = x[2];
    Real Q_beta  = x[3];

    //-------------
    //! Get the dof of the interfaces (to extract the paramaters values)
    //-------------
    //! identity of the point on the boundary of domain alpha
    //! (on the right of the domain alpha)
    UInt rightDof = _M_dof_alpha;
    //! identity of the point on the boundary of domain beta
    //! (on the left of the domain beta)
    UInt leftDof  = _M_dof_beta;
    //-------------

    // *******************************************************
    //! 0/ Continuity of the flux
    // *******************************************************
    f(0) = Q_beta - Q_alpha;
    //! Jacobian
    jac( 0, 0 ) =  0.; //!< df0/dA_alpha
    jac( 0, 1 ) = -1.; //!< df0/dQ_alpha
    jac( 0, 2 ) =  0.; //!< df0/dA_beta
    jac( 0, 3 ) =  1.; //!< df0/dQ_beta

    // *******************************************************
    //! 1/ Continuity of the total pressure
    // *******************************************************
    f(1) = (   _M_fluxFun_beta.totalPressure(  A_beta,  Q_beta,  leftDof )
             - _M_fluxFun_alpha.totalPressure( A_alpha, Q_alpha, rightDof ) );

    //! Jacobian
    //!df1/dA_alpha:
    jac( 1, 0 ) = - _M_fluxFun_alpha.totalPressureDiff(A_alpha, Q_alpha, 1, rightDof );
    //!df1/dQ_alpha:
    jac( 1, 1 ) = - _M_fluxFun_alpha.totalPressureDiff(A_alpha, Q_alpha, 2, rightDof );
    //!df1/dA_beta:
    jac( 1, 2 ) = + _M_fluxFun_beta.totalPressureDiff( A_beta,  Q_beta,  1, leftDof );
    //!df1/dQ_beta:
    jac( 1, 3 ) = + _M_fluxFun_beta.totalPressureDiff( A_beta,  Q_beta,  2, leftDof );

    // *******************************************************
    //! 2/ ALPHA : right compatibility condition (for tube alpha)
    // *******************************************************
    //-------------------------------------
    //! Compute the eigen values and vectors
    boundaryPoint   = _M_edge_alpha.pt2().x(); //!< point on the boundary (on the right of the edge!)
    internalBdPoint = _M_edge_alpha.pt1().x(); //!< neighboring point (internal)

    //! values of U on the boundary
    U_boundary   = Vec2D ( _M_Un_alpha_bd );
    //! values of U on the neighboring node of the boundary point
    U_internalBd = Vec2D ( _M_Un_alpha_int );

    //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
    _M_fluxFun_alpha.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                  eigval1, eigval2,
                                                  left_eigvec1.first, left_eigvec1.second,
                                                  left_eigvec2.first, left_eigvec2.second,
                                                  rightDof);

    if ( verbose > 3 )
        std::cout << "EigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;


    ASSERT( eigval1 > 0. && eigval2 < 0. ,
            "The eigenvalues do not have the expected signs.");

    //-------------------------------------
    //!compatibility condition
    //-------------------------------------
    //! first characteristics
    //! interpolation of U at the foot of the 1rst characteristics
    U_charact_pt = interpolLinear(boundaryPoint, internalBdPoint,
                                  _M_time_step, eigval1,
                                  U_boundary, U_internalBd);

    //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
    rhsBC1 = dot( left_eigvec1 , U_charact_pt );

    //! take into account the (quasi linear) source term
    //!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
    //! THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
    qlSource.first  = _M_sourceFun_alpha.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                           1, rightDof);
    qlSource.second = _M_sourceFun_alpha.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                           2, rightDof);

    //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
    rhsBC1 -= _M_time_step * dot( left_eigvec1 , qlSource );

    //! return f(2): left_eigvec1 dot (A_alpha_n+1, Q_alpha_n+1)
    f(2) = left_eigvec1.first * A_alpha + left_eigvec1.second * Q_alpha - rhsBC1;
    //! Jacobian
    jac( 2, 0 ) =  left_eigvec1.first; //!< df2/dA_alpha
    jac( 2, 1 ) =  left_eigvec1.second;//!< df2/dQ_alpha
    jac( 2, 2 ) =  0.; //!< df2/dA_beta
    jac( 2, 3 ) =  0.; //!< df2/dQ_beta

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // *******************************************************
    //! 3/ BETA : left compatibility condition (for tube beta).
    // *******************************************************
    //-------------------------------------
    //! Compute the eigen values and vectors
    boundaryPoint   = _M_edge_beta.pt1().x(); //!< point on the boundary (on the left of the edge!)
    internalBdPoint = _M_edge_beta.pt2().x(); //!< neighboring point (internal)

    //! values of U on the boundary
    U_boundary   = Vec2D ( _M_Un_beta_bd );
    //! values of U on the neighboring node of the boundary point
    U_internalBd = Vec2D ( _M_Un_beta_int );

    //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
    _M_fluxFun_beta.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                                 eigval1, eigval2,
                                                 left_eigvec1.first, left_eigvec1.second,
                                                 left_eigvec2.first, left_eigvec2.second,
                                                 leftDof);

    //------------------------------------------------------------

    if ( verbose > 3 )
        std::cout << "\nEigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;

    ASSERT( eigval1 > 0. &&  eigval2 < 0.  ,
            "The eigenvalues do not have the expected signs.");

    //-------------------------------------
    //!compatibility condition
    //-------------------------------------
    //! second characteristics
    //! interpolation of U at the foot of the 2nd characteristics
    U_charact_pt = interpolLinear(boundaryPoint, internalBdPoint,
                                  _M_time_step, eigval2,
                                  U_boundary, U_internalBd);

    //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
    rhsBC2 = dot( left_eigvec2 , U_charact_pt );

    //! take into account the (quasi linear) source term
    //!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
    //! THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
    qlSource.first  = _M_sourceFun_beta.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                          1, leftDof);
    qlSource.second = _M_sourceFun_beta.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                          2, leftDof);

    //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
    rhsBC2 -= _M_time_step * dot( left_eigvec2 , qlSource );

    //! return f(3): left_eigvec2 dot (A_beta_n+1, Q_beta_n+1)
    f(3) = left_eigvec2.first * A_beta + left_eigvec2.second * Q_beta - rhsBC2;
    //! Jacobian
    jac( 3, 0 ) =  0.; //!< df3/dA_alpha
    jac( 3, 1 ) =  0.; //!< df3/dQ_alpha
    jac( 3, 2 ) =  left_eigvec2.first; //!< df3/dA_beta
    jac( 3, 3 ) =  left_eigvec2.second;//!< df3/dQ_beta

}

Real Interface2Vessels::dot(const Vec2D& vec1, const Vec2D& vec2) const
{
    return vec1.first * vec2.first + vec1.second * vec2.second;
}

Interface2Vessels::Vec2D
Interface2Vessels::interpolLinear(const Real& point_bound, const Real& point_internal,
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

}
