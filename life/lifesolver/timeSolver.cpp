/*

  TimeSolver class

*/

#include "timeSolver.hpp"




/* for types of vecors here I use ScalUnknown<Vector>
   for types of matrices MSRMatr<double>
   I liked the idea of templates of types, but
   Aztec allows only MSR or VBR, anyway I've seen
   MSRMatr and ScalUnknown are used a lot
 */

namespace LifeV
{


  /* solves u_k on a generic equation of the type:
     M=theta M_k + (1-theta)M_k_1
     (M/dt+A_k*theta)u_k=theta*b_k+(1-theta)b_k_1-A_k_1(1-theta)u_k_1+M/dt*u_k_1
     note that matrix M and A (and so P) must have the same pattern
     for theta in [0,1] we have different methods:
     theta=0	Eulero esplicit
     theta=1	Eulero implicit
     theta=0.5	Crank Nicolson (better)
   */
TimeSolver::TimeSolver(Real _dt,const_matrix_type _M_k,const_matrix_type _M_k_1,
		const_matrix_type _A_k,const_matrix_type _A_k_1,
             	const_vector_type _b_k,const_vector_type _b_k_1,
		const_vector_type _u_k_1,Real _theta):
	     SolverAztec("data"),
	     P(matrix_type(*_M_k.Patt())),q(vector_type(_M_k.Patt()->nRows())),
	     u_k(_u_k_1), b_k(_b_k),b_k_1(_b_k_1),
	     u_k_1(_u_k_1),M_k(_M_k),M_k_1(_M_k_1),A_k(_A_k),
	     A_k_1(_A_k_1),theta(_theta),dt(_dt)
{
  CONSTRUCTOR("TimeSolver 1");
}



}

