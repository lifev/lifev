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
int computeInterface2TubesValues( const Vec2D& Un_alpha_bd,
				  const Vec2D& Un_alpha_int,
				  const Edge1D& edge_alpha, const UInt& dof_alpha,
				  const NonLinearFluxFun1D&   flux_alpha,
				  const NonLinearSourceFun1D& source_alpha,
				  const Vec2D& Un_beta_bd,
				  const Vec2D& Un_beta_int,
				  const Edge1D& edge_beta, const UInt& dof_beta,
				  const NonLinearFluxFun1D&   flux_beta,
				  const NonLinearSourceFun1D& source_beta,
				  const Real& time_step,
				  Vec2D& bcDir_alpha,
				  Vec2D& bcDir_beta
				  );
int test_newton( );

void f_jac( const Vec2D& Un_alpha_bd,
	    const Vec2D& Un_alpha_int,
	    const Edge1D& edge_alpha, const UInt& dof_alpha,
	    const NonLinearFluxFun1D&   fluxFun_alpha,
	    const NonLinearSourceFun1D& sourceFun_alpha,
	    const Vec2D& Un_beta_bd,
	    const Vec2D& Un_beta_int,
	    const Edge1D& edge_beta, const UInt& dof_beta,
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
    std::cout << "cfl " << cfl << std::endl;
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
    /*const char* data_file_name = command_line.follow("data", 2, "-f","--file");
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
    /*
    std::cout << "======\n\tLinear model tube 1" << std::endl;
    onedparamLin.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    std::cout << "======\n\tLinear model tube 2" << std::endl;
    onedparamLin_t2.showMeData(std::cout);
    std::cout << "-----------------------------" << std::endl;
    */

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

    //! interface values
    Vec2D bcDir_t1, bcDir_t2;

    std::cout << "++++++++++++++++++++++++++++++\n\tTemporal loop starting ... " << std::endl;
    // Temporal loop
    //
    Chrono chrono;
    int count = 0;
    for ( Real time = startT + dt ; time <= T ; time += dt ) {
        count++;
        std::cout << "Iteration " <<  count  << ", t = " << time
                  << "s... \n";
        chrono.start();

	Edge1D foo = onedm.RightEdge();
	

 	std::cout << "right edge " << foo.pt1().x() <<  " " << foo.pt2().x() << std::endl;

       //! compute the interface values
	Vec2D bcDir_t1, bcDir_t2;
	//! compute interface values
	int cvg_newton = computeInterface2TubesValues
	  ( onedm.BCValuesRight(),
	    onedm.BCValuesInternalRight(),
	    onedm.RightEdge(), 
	    onedm.RightNodeId(),
	    onedm.FluxFun(), 
	    onedm.SourceFun(),
	    onedm_t2.BCValuesLeft(),
	    onedm_t2.BCValuesInternalLeft(),
	    onedm_t2.LeftEdge(), 
	    onedm_t2.LeftNodeId(),
	    onedm_t2.FluxFun(), 
	    onedm_t2.SourceFun(),
	    onedm.timestep(),
	    bcDir_t1, bcDir_t2
	    );

	std::cout << "bc dir tube 1 " << bcDir_t1.first <<  " " << bcDir_t1.second 
		  << "\nbc dir tube 2 " << bcDir_t2.first <<  " " << bcDir_t2.second << std::endl;

        //! set the interface values
        onedm.setBCValuesRight( bcDir_t1.first, bcDir_t1.second );
        onedm_t2.setBCValuesLeft( bcDir_t2.first, bcDir_t2.second );

        ASSERT_PRE( !cvg_newton,"Newton iteration for interface values computation not achieved.");

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


int computeInterface2TubesValues( const Vec2D& Un_alpha_bd,
				  const Vec2D& Un_alpha_int,
				  const Edge1D& edge_alpha, const UInt& dof_alpha,
				  const NonLinearFluxFun1D&   fluxFun_alpha,
				  const NonLinearSourceFun1D& sourceFun_alpha,
				  const Vec2D& Un_beta_bd,
				  const Vec2D& Un_beta_int,
				  const Edge1D& edge_beta, const UInt& dof_beta,
				  const NonLinearFluxFun1D&   fluxFun_beta,
				  const NonLinearSourceFun1D& sourceFun_beta,
				  const Real& time_step,
				  Vec2D& bcDir_alpha,
				  Vec2D& bcDir_beta
				  )

//int test_newton()
{
    const UInt f_size(4);
    
    //! unknown of non linear equation f(x) = 0
    Vector<Real> x(f_size);
    //! non linear function f
    Vector<Real> f(f_size);
    //! jacobian of the non linear function
    Matrix<Real> jac(f_size,f_size);
    //! transpose of the jacobian of the non linear function
    Matrix<Real> jac_trans(f_size,f_size);

    //! tmp matrix for lapack lu inversion
    Vector<Int> ipiv(f_size);

    x[0] = Un_alpha_bd.first; //!< A_alpha
    x[1] = Un_alpha_bd.second;//!< Q_alpha
    x[2] = Un_beta_bd.first;  //!< A_beta
    x[3] = Un_beta_bd.second; //!< Q_beta

    //! lapack variable
    int INFO[1] = {0};
    int NBRHS[1] = {1};//  nb columns of the rhs := 1.
    int NBU[1] = {f_size};

    //! newton raphson iteration
    //! (use the newton class???)
    for ( UInt iter = 0 ; iter < 100 ; iter ++) {

        std::cout << "\titer = " << iter << std::endl;

        //! compute f(x) and its jacobian df(x)
        f_jac( Un_alpha_bd, Un_alpha_int, edge_alpha, dof_alpha, fluxFun_alpha, sourceFun_alpha,
	       Un_beta_bd, Un_beta_int, edge_beta, dof_beta, fluxFun_beta, sourceFun_beta,
               time_step, x, f, jac);
      

	std::cout << "x : " << x << "\nf : " << f << "\njac : " << jac << std::endl;

	//! transpose to pass to fortran storage (lapack!)
	jac_trans = trans(jac);
	
	//! Compute f <-  ( df(x)^{-1} f(x) ) (lu dcmp)
	dgesv_(NBU, NBRHS, &jac_trans(0,0), NBU , &ipiv(0), &f(0), NBU, INFO);
        ASSERT_PRE(!INFO[0],"Lapack LU resolution of y = df(x)^{-1} f(x) is not achieved.");
	x += - f;

        std::cout << "x : " << x << "\nf : " << f << std::endl;

	//! convergence if Q_alpha == Q_beta
	if ( std::fabs( x[1] - x[3] ) < 1e-8 ) {
	    bcDir_alpha = Vec2D( x[0], x[1] );
	    bcDir_beta  = Vec2D( x[2], x[3] );
	    return 0;
	}
	  
    }
    //! no convergence
    return 1;

}

void f_jac( const Vec2D& Un_alpha_bd,
            const Vec2D& Un_alpha_int,
            const Edge1D& edge_alpha, const UInt& dof_alpha,
            const NonLinearFluxFun1D&   fluxFun_alpha,
            const NonLinearSourceFun1D& sourceFun_alpha,
            const Vec2D& Un_beta_bd,
            const Vec2D& Un_beta_int,
            const Edge1D& edge_beta, const UInt& dof_beta,
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

    std::cerr << "x = " << x << std::endl;
   

    Real A_alpha = x[0];
    Real Q_alpha = x[1];
    Real A_beta  = x[2];
    Real Q_beta  = x[3];

    //-------------
    //! Get the dof of the interfaces (to extract the paramaters values)
    //-------------
    //! identity of the point on the boundary of domain alpha
    //! (on the right of the domain alpha)
    UInt rightDof = dof_alpha; 
    //! identity of the point on the boundary of domain beta
    //! (on the left of the domain beta)
    UInt leftDof  = dof_beta;
    //-------------

    std::cout << "rightDof = " << rightDof 
	      << "\tleftDof = " << leftDof << std::endl;


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

    std::cout << "boundaryPoint = " <<  boundaryPoint << " internalBdPoint = " << internalBdPoint
	      << "U_boundary = " <<  U_boundary.first << " " << U_boundary.second
	      << std::endl;


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

#if 0
//! test function for newton iter
int f_jac( const Vector<Real>& x,
	   Vector<Real>& f,
	   Matrix<Real>& jac
	   )
{
    // *******************************************************
    f(0) = 2*x(0) + 3*x[1] - 11;
    //! Jacobian
    jac( 0, 0 ) =  2.;
    jac( 0, 1 ) =  3.;

    // *******************************************************
    f(1) = 2*x[0]*x[0] + x[0] + x[1]*x[1] - 12;
    //! Jacobian
    jac( 1, 0 ) = 4*x[0] + 1;
    jac( 1, 1 ) = 2*x[1];
}
#endif

