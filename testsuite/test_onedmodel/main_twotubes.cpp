/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
       Date: 2004-10-26

  Copyright (C) 2004 Politecnico di Milano

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
   \file main_twotubes.cpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-10-26
 */

#include "lifeV.hpp"
#include "chrono.hpp"
#include "dataOneDModel.hpp"
#include "oneDModelSolver.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"
#include <sstream>

using namespace LifeV;

typedef pair< Real, Real > Vec2D;
void computeInterface2TubesValues( const Vec2D& Un_alpha_bd,
                                   const Vec2D& Un_alpha_int,
                                   const Edge1D& edge_alpha,
                                   const NonLinearFluxFun1D&   flux_alpha,
                                   const NonLinearSourceFun1D& source_alpha,
                                   const Vec2D& Un_beta_bd,
                                   const Vec2D& Un_beta_int,
                                   const Edge1D& edge_beta,
                                   const NonLinearFluxFun1D&   flux_beta,
                                   const NonLinearSourceFun1D& source_beta,
                                   const Real& time_step,
                                   Vec2D& bcDir_alpha,
                                   Vec2D& bcDir_beta
                                   );

void f_jac( const Vec2D& Un_alpha_bd,
            const Vec2D& Un_alpha_int,
            const Edge1D& edge_alpha,
            const NonLinearFluxFun1D&   fluxFun_alpha,
            const NonLinearSourceFun1D& sourceFun_alpha,
            const Vec2D& Un_beta_bd,
            const Vec2D& Un_beta_int,
            const Edge1D& edge_beta,
            const NonLinearFluxFun1D&   fluxFun_beta,
            const NonLinearSourceFun1D& sourceFun_beta,
            const Real& time_step,
            const Vector<Real>& x,
            Vector<Real>& f,
            Matrix<Real>& jac
            );

//! 2D dot product
Real dot(const Vec2D& vec1, const Vec2D& vec2)
{
    return vec1.first * vec2.first + vec1.second * vec2.second;
}

Vec2D interpolLinear(const Real& point_bound, const Real& point_internal,
                     const Real& deltaT, const Real& eigenvalue,
                     const Vec2D& U_bound, const Vec2D& U_intern) 
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


int main(int argc, char** argv)
{

    //! ********** Reading from data file ******************************************
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    /*
      GetPot data_file(data_file_name);
    */

    GetPot data_file( "datanl" );

    OneDNonLinModelParam onedparamNL( data_file );
    onedparamNL.initParam(4e6);

    GetPot data_file_t2( "datanl2" );
    OneDNonLinModelParam onedparamNL_t2( data_file_t2 );
    onedparamNL_t2.initParam(8e6);


    LinearSimpleParam onedparamLin(data_file);
    LinearSimpleParam onedparamLin_t2(data_file_t2);

    std::cout << "======\n\tNon Linear model tube 1" << std::endl;
    onedparamNL.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "======\n\tNon Linear model tube 2" << std::endl;
    onedparamNL_t2.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;

    std::cout << "======\n\tLinear model tube 1" << std::endl;
    onedparamLin.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "======\n\tLinear model tube 2" << std::endl;
    onedparamLin_t2.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;

    //  OneDModelSolver onedm(data_file, onedparamLin);
    OneDModelSolver onedm   (data_file, onedparamNL);
    OneDModelSolver onedm_t2(data_file_t2, onedparamNL_t2);

    std::cout << "======\n\ttube 1" << std::endl;
    onedm.showMeData();
    std::cout << "======\n\ttube 2" << std::endl;
    onedm_t2.showMeData();
    //  onedm.showMeHandler(cout, 6);



    // Initialization
    //
    Real dt     = onedm.timestep();
    Real startT = onedm.inittime();
    Real T      = onedm.endtime();

    Real dt_t2     = onedm_t2.timestep();
    Real startT_t2 = onedm_t2.inittime();
    Real T_t2      = onedm_t2.endtime();

    ASSERT_PRE( dt == dt_t2, "Same time step required!");
    ASSERT_PRE( startT == startT_t2, "Same initial time required!");
    ASSERT_PRE( T == T_t2, "Same final time required!");

    ASSERT( onedm.xRight() == onedm_t2.xLeft(),
            "Same interface point required!");

    /*
      Real u1_0 = 0.; //! constant initial condition
      Real u2_0 = 0.; //! constant initial condition
    */
    //! tube 1
    Real u1_0 = onedparamNL.Area0(0); //! constant initial condition
    Real u2_0 = 0.; //! constant initial condition

    std::cout << "initialize tube 1 with constant (u1_0, u2_0)" << std::endl;
    onedm.initialize(u1_0, u2_0);

    //! tube 2
    u1_0 = onedparamNL_t2.Area0(0); //! constant initial condition
    u2_0 = 0.; //! constant initial condition

    std::cout << "initialize tube 2 with constant (u1_0, u2_0)" << std::endl;
    onedm_t2.initialize(u1_0, u2_0);

    std::cout << "startT T dt " << startT << " " <<  T << " " << dt << std::endl;

    char ch;
    std::cout << "Hit return to continue" << std::endl;
    std::cin.get(ch);


    Real bcR1, bcR2, bcL1, bcL2;
    bcR1 = 0.;
    bcR2 = 0.;
    bcL1 = 0.;
    bcL2 = 0.;


    // Temporal loop
    //
    Chrono chrono;
    int count = 0;
    for ( Real time = startT + dt ; time <= T ; time += dt ) {
        count++;
        std::cout << "Iteration " <<  count  << ", t = " << time
                  << "s... \n";
        chrono.start();

        //! compute the interface values

        //! set the interface values
        onedm.setBCValuesRight( bcR1, bcR2 );
        onedm_t2.setBCValuesLeft( bcL1, bcL2 );

        //! tube 1
        onedm.timeAdvance( time );
        onedm.iterate( time , count );

        //! tube 2
        onedm_t2.timeAdvance( time );
        onedm_t2.iterate( time , count );

        chrono.stop();
        std::cout << "Iteration " <<  count  << " computed in " << chrono.diff() << " s.\n" << std::endl;

        if ( data_file( "miscellaneous/show_graceplot", 0 ) )
            onedm.gplot();

    }

    std::cout << "Hit return to close" << std::endl;
    std::cin.get(ch);

    return 0;

}

void computeInterface2TubesValues( const Vec2D& Un_alpha_bd,
                                   const Vec2D& Un_alpha_int,
                                   const Edge1D& edge_alpha,
                                   const NonLinearFluxFun1D&   fluxFun_alpha,
                                   const NonLinearSourceFun1D& sourceFun_alpha,
                                   const Vec2D& Un_beta_bd,
                                   const Vec2D& Un_beta_int,
                                   const Edge1D& edge_beta,
                                   const NonLinearFluxFun1D&   fluxFun_beta,
                                   const NonLinearSourceFun1D& sourceFun_beta,
                                   const Real& time_step,
                                   Vec2D& bcDir_alpha,
                                   Vec2D& bcDir_beta
                                   )
{
    //! unknown of non linear equation f(x) = 0
    Vector<Real> x(4);
    //! non linear function f
    Vector<Real> f(4);
    //! jacobian of the non linear function
    Matrix<Real> jac(4,4);

    x[0] = Un_alpha_bd.first;
    x[1] = Un_alpha_bd.second;
    x[2] = Un_beta_bd.first;
    x[3] = Un_beta_bd.second;
 
    //! lapack variable
    int INFO[1] = {0};
    int NBRHS[1] = {1};//  nb columns of the rhs := 1.
    int NBU[1] = {4};

    //! newton raphson iteration
    for ( UInt iter = 0 ; iter < 100 ; iter ++) {
        //! compute f(x) and its jacobian df(x)
        f_jac( Un_alpha_bd, Un_alpha_int, edge_alpha, fluxFun_alpha, sourceFun_alpha,
	       Un_beta_bd, Un_beta_int, edge_beta, fluxFun_beta, sourceFun_beta,
	       time_step, x, f, jac);
	
	//   jac <- L and Lt where L Lt is the Cholesky factorization of df(x)
        dpotrf_("L", NBU, jac , NBU , INFO );
        ASSERT_PRE(!INFO[0],"Lapack factorization of jacobian matrix is not achieved.");

        // Compute f <-  ( df(x)^{-1} f(x) ) (solve triangular system)
        dtrtrs_("L", "N", "N", NBU, NBRHS, jac, NBU, f, NBU, INFO);
        ASSERT_PRE(!INFO[0],"Lapack Computation y = df(x)^{-1} f(x)  is not achieved.");

	x += - f;
    }

}

void f_jac( const Vec2D& Un_alpha_bd,
            const Vec2D& Un_alpha_int,
            const Edge1D& edge_alpha,
            const NonLinearFluxFun1D&   fluxFun_alpha,
            const NonLinearSourceFun1D& sourceFun_alpha,
            const Vec2D& Un_beta_bd,
            const Vec2D& Un_beta_int,
            const Edge1D& edge_beta,
            const NonLinearFluxFun1D&   fluxFun_beta,
            const NonLinearSourceFun1D& sourceFun_beta,
            const Real& time_step,
            const Vector<Real>& x,
            Vector<Real>& f,
            Matrix<Real>& jac
            )
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

    // *******************************************************
    //! Continuity of the flux
    // *******************************************************
    f(0) = Q_beta - Q_alpha;
    //! Jacobian
    jac( 0, 0 ) =  0.; //!< df0/dA_alpha
    jac( 0, 1 ) = -1.; //!< df0/dQ_alpha
    jac( 0, 2 ) =  0.; //!< df0/dA_beta
    jac( 0, 3 ) =  1.; //!< df0/dQ_beta

    // *******************************************************
    //! Continuity of the total pressure
    // *******************************************************
    f(1) = (   fluxFun_beta.totalPressure(  A_beta,  Q_beta,  leftDof ) 
	     - fluxFun_alpha.totalPressure( A_alpha, Q_alpha, rightDof ) );
    //! Jacobian
    //!df1/dA_alpha:
    jac( 1, 0 ) = fluxFun_alpha.totalPressureDiff(A_alpha, Q_alpha, 1, rightDof );
    //!df1/dQ_alpha:
    jac( 1, 1 ) = fluxFun_alpha.totalPressureDiff(A_alpha, Q_alpha, 2, rightDof );
    //!df1/dA_beta:
    jac( 1, 2 ) = fluxFun_beta.totalPressureDiff( A_beta,  Q_beta,  1, leftDof );
    //!df1/dQ_beta:
    jac( 1, 3 ) = fluxFun_beta.totalPressureDiff( A_beta,  Q_beta,  2, leftDof );
    

    // *******************************************************
    //! ALPHA : right compatibility condition (for tube alpha)
    // *******************************************************
    //-------------------------------------
    //! Compute the eigen values and vectors
    boundaryPoint   = edge_alpha.pt2().x(); //!< point on the boundary (on the right of the edge!)
    internalBdPoint = edge_alpha.pt1().x(); //!< neighboring point (internal)

    //! values of U on the boundary
    U_boundary   = Vec2D ( Un_alpha_bd );
    //! values of U on the neighboring node of the boundary point
    U_internalBd = Vec2D ( Un_alpha_int );

    //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
    fluxFun_alpha.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                               eigval1, eigval2,
                                               left_eigvec1.first, left_eigvec1.second,
                                               left_eigvec2.first, left_eigvec2.second,
                                               rightDof);

    if ( verbose > 1 )
        std::cout << "EigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;


    ASSERT( eigval1 > 0. && eigval2 < 0. ,
            "The eigenvalues do not have the expected signs.");

    //-------------------------------------
    //!compatibility condition
    //-------------------------------------
    //! first characteristics
    //! interpolation of U at the foot of the 1rst characteristics
    U_charact_pt = interpolLinear(boundaryPoint, internalBdPoint,
                                  time_step, eigval1,
                                  U_boundary, U_internalBd);

    //! rhsBC1 = left_eigvec1 dot U(tn, z = charact_pt1)
    rhsBC1 = dot( left_eigvec1 , U_charact_pt );

    //! take into account the (quasi linear) source term
    //!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
    //! THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
    qlSource.first  = sourceFun_alpha.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                        1, rightDof);
    qlSource.second = sourceFun_alpha.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                        2, rightDof);

    //! rhsBC1 = rhsBC1 - deltaT * left_eigvec1 dot qlSource(tn, z = charact_pt1)
    rhsBC1 -= time_step * dot( left_eigvec1 , qlSource );

    //! return f(2): left_eigvec1 dot (A_alpha_n+1, Q_alpha_n+1)
    f(2) = left_eigvec1.first * A_alpha + left_eigvec1.second * Q_alpha - rhsBC1;
    //! Jacobian
    jac( 2, 0 ) =  left_eigvec1.first; //!< df2/dA_alpha
    jac( 2, 1 ) =  left_eigvec1.second;//!< df2/dQ_alpha
    jac( 2, 2 ) =  0.; //!< df2/dA_beta
    jac( 2, 3 ) =  0.; //!< df2/dQ_beta

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
    // *******************************************************
    //! BETA : left compatibility condition (for tube beta).
    // *******************************************************
    //-------------------------------------
    //! Compute the eigen values and vectors
    boundaryPoint   = edge_beta.pt1().x(); //!< point on the boundary (on the left of the edge!)
    internalBdPoint = edge_beta.pt2().x(); //!< neighboring point (internal)

    //! values of U on the boundary
    U_boundary   = Vec2D ( Un_beta_bd );
    //! values of U on the neighboring node of the boundary point
    U_internalBd = Vec2D ( Un_beta_int );

    //! compute the eigenvalues/eigenvectors of the flux jacobian (computed on the boundary point)
    fluxFun_beta.jacobian_EigenValues_Vectors(U_boundary.first, U_boundary.second,
                                              eigval1, eigval2,
                                              left_eigvec1.first, left_eigvec1.second,
                                              left_eigvec2.first, left_eigvec2.second,
                                              leftDof);

    //------------------------------------------------------------

    if ( verbose > 1 )
        std::cout << "\nEigenValue 1 " << eigval1 << " EigenValue 2 " << eigval2 << std::endl;

    ASSERT( eigval1 > 0. &&  eigval2 < 0.  ,
            "The eigenvalues do not have the expected signs.");

    //-------------------------------------
    //!compatibility condition
    //-------------------------------------
    //! second characteristics
    //! interpolation of U at the foot of the 2nd characteristics
    U_charact_pt = interpolLinear(boundaryPoint, internalBdPoint,
                                  time_step, eigval2,
                                  U_boundary, U_internalBd);

    //! rhsBC2 = left_eigvec2 dot U(tn, z = charact_pt2)
    rhsBC2 = dot( left_eigvec2 , U_charact_pt );

    //! take into account the (quasi linear) source term
    //!BEWARE: HERE THE PARAMETERS ARE TAKEN AT rightDof
    //! THEY SHOULD BE TAKEN AT THE CHARACTERISTICS!!
    qlSource.first  = sourceFun_beta.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                       1, leftDof);
    qlSource.second = sourceFun_beta.QuasiLinearSource(U_charact_pt.first, U_charact_pt.second,
                                                       2, leftDof);

    //! rhsBC2 = rhsBC2 - deltaT * left_eigvec2 dot qlSource(tn, z = charact_pt2)
    rhsBC2 -= time_step * dot( left_eigvec2 , qlSource );

    //! return f(3): left_eigvec2 dot (A_beta_n+1, Q_beta_n+1)
    f(3) = left_eigvec2.first * A_beta + left_eigvec2.second * Q_beta - rhsBC2;
    //! Jacobian
    jac( 3, 0 ) =  0.; //!< df3/dA_alpha
    jac( 3, 1 ) =  0.; //!< df3/dQ_alpha
    jac( 3, 2 ) =  left_eigvec2.first; //!< df3/dA_beta
    jac( 3, 3 ) =  left_eigvec2.second;//!< df3/dQ_beta


}

