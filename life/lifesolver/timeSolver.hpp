/*

  $Id: timeSolver.hpp,v 1.2 2004-12-14 15:12:33 winkelma Exp $
  TimeSolver class
  the purpouse of this class is to handle more easily
  the solution of M\partial u/\partial t + Au=b
  I've generalized the method introducing a parameter \theta,
  so:
  theta=0	Euler esplicit
  theta=1	Euler implicit
  theta=0.5	Crank Nicolson (best method)

*/

#ifndef _TIME_SOLVER_HPP_
#define _TIME_SOLVER_HPP_

#include "SolverAztec.hpp"
#include "vecUnknown.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include <vector>
#include "bcHandler.hpp"
#include "currentFE.hpp"
#include "currentBdFE.hpp"
#include "dof.hpp"
#include "pattern.hpp"
#include <stdio.h>	//for NULL


#include <boost/shared_ptr.hpp>

//#include <boost/numeric/ublas/matrix_sparse.hpp>  //gives only problems!
//#include <boost/numeric/ublas/io.hpp>
//TODO think to a possible convertion between
//using MSRMatr<TYPE> as default and switching to sparse_matrix<TYPE>

/* for types of vecors here I use ScalUnknown<Vector>
   for types of matrices MSRMatr<double>
   I liked the idea of templates of types, but
   Aztec allows only MSR or VBR, anyway I've seen
   MSRMatr and ScalUnknown are used a lot
 */

namespace LifeV
{

//template <typename VectorType,typename MatrixType>
class TimeSolver: public SolverAztec
{
public:
  typedef ScalUnknown<Vector> vector_type;
  typedef MSRMatr<Real> matrix_type;
  typedef const ScalUnknown<Vector> & const_vector_type;
  typedef const MSRMatr<Real> & const_matrix_type;

  TimeSolver():SolverAztec("data"),P(),q(),u_k(),
  		b_k(),b_k_1(),u_k_1(),M_k(),M_k_1(),A_k(),
                A_k_1(),theta(0.5),dt(1.0)

  {
    CONSTRUCTOR("TimeSolver 2");
  }
  ~TimeSolver()
  {
    DESTRUCTOR("TimeSolver");
  }
  TimeSolver(UInt _dim):SolverAztec("data"),
  		P(),q(_dim),u_k(_dim),
  		b_k(_dim),b_k_1(_dim),u_k_1(_dim),M_k(),M_k_1(),A_k(),
                A_k_1(),theta(0.5),dt(1.0)
  {
    CONSTRUCTOR("TimeSolver 3");
  }

  	
  /* solves u_k on a generic equation of the type:
     (M/dt+A_k*theta)u_k=theta*b_k+(1-theta)b_k_1-A_k_1(1-theta)u_k_1+M/dt*u_k_1
   */
  TimeSolver(Real _dt,const_matrix_type _M_k,const_matrix_type _M_k_1
		,const_matrix_type _A_k,const_matrix_type _A_k_1,
             	const_vector_type _b_k,const_vector_type _b_k_1,
		const_vector_type _u_k_1,Real _theta);
  void setM_k(const_matrix_type _M_k)
  {
    M_k=_M_k;
    initPq(M_k);	
  }
  void setM_k_1(const_matrix_type _M_k_1)
  {
    M_k_1=_M_k_1;
  }
  void setA_k(const_matrix_type _A_k)
  {
    A_k=_A_k;
  }
  void setb_k(const_vector_type _b_k)
  {
    b_k=_b_k;
  }
  void setb_k_1(const_vector_type _b_k_1)
  {
    b_k_1=_b_k_1;
  }
  void setu_k_1(const_vector_type _u_k_1)
  {
    u_k_1=_u_k_1;
    u_k=_u_k_1;		//to initialize size
  }
  void initPq(const_matrix_type _M)	
  {
    //just in case of empty constructor
    P=matrix_type(*_M.Patt());
    q=vector_type(_M.Patt()->nRows());
  }

  void setA_k_1(const_matrix_type _A_k_1)
  {
    A_k_1=_A_k_1;
  }
  void settheta(Real _theta)
  {
    theta=_theta;
  }
  void setdt(Real _dt)
  {
    dt=_dt;
  }

  /* I need flexibility using this class, so I don't update automatically */
  void timeIncrement()
  {
    u_k_1=u_k;
    A_k_1=A_k;
    b_k_1=b_k;
    M_k_1=M_k;
  }

  void update()	//update the matrices/vectors of the equation to solve
  {
Debug(6013)<<"entering call to TimeSolver::update()\n";
    dti=1/dt;
    /* MSRMatr::operator + is still not defined, I implement it at least for
       2matrix with same pattern.
       MSRMatr::flop (4 args) is not commmited stil, I added it for better
       efficiency */
    //P=(M*dti+A_k*theta);
    matrix_type M=M_k;
    M.flop(theta,M_k,(1-theta),M_k_1);

    P.flop(dti,M,theta,A_k);

/* MM: I wan't able to generate some sparse_matrix matrices instead
   of MSRMatr, becouse in no debug mode multiplication
   ublas...vectors * MSRMatr doesn't work.
   This is time comsuming anyway but there is a lot of
   code that uses MSRMatr, so it's not easy to switch all to
   sparse_matrix type
*/

    q=b_k*theta+b_k_1*(1.0-theta)-(A_k_1.operator*(u_k_1*(1.0-theta)))+ (M.operator*(u_k_1*(1/dt)));


//  this works in debug mode ublasvectors*MSRMatr allowed
//    q=b_k*theta+b_k_1*(1.0-theta)-(A_k_1*(u_k_1*(1.0-theta)))+ (M*(u_k_1*(1/dt)));
    setMatrix(P);
  }

  /* call update before solve if solution could change */
  vector_type solve()
  {
    update();     //check if update at every step
    SolverAztec::solve(u_k,q);
    return u_k;
  }



private:
  //I construct a system: P*u_k=q to get the solution
  matrix_type P;
  vector_type q;

  vector_type u_k;
  vector_type b_k,b_k_1,u_k_1;
  matrix_type M_k,M_k_1,A_k,A_k_1;
  Real theta,dt;

  //constructed vars to improve a bit the efficiency
  Real dti;



};






}
#endif
