#ifndef _INEXACTLU_H
#define _INEXACTLU_H

/*----------------------------------------------------------------------

  A class to manage the three steps of the resolution of the problem
  
  [ A  B^T ] [ x ]    [ f ]
  [        ] [   ]  = [   ]
  [ B  0   ] [ y ]    [ g ]
  
  by inexact LU factorization.

  We denote by A1, \tilde A and A2 three approximations of A

  Step 1:                   A1 x0 = f
  Step 2: B (\tilde A)^{-1} B^T y = B x0 - g
  Step 3:                       x = x0 +  A2^{-1} B^T y 

  The operator of Step 2 is a class "SaddlePoint" (see saddlePointCG.h)

  (
  ! Example:
  ! if A = 1/dt M + K (M: mass, K: stiffness + advection)
  !
  ! If A1 = \tilde A = A2 = A           : Uzawa algorithm (exact factorisation)
  ! If A1 = A2 = A and \tilde A = 1/dt M: Yosida algorithm
  ! If A1 = A and A2 = \tilde A = 1/dt M: algebraic Chorin-Temam algorithm
  )
  Remark: the Uzawa algorithm can be performed directly on the class
          SaddlePoint.
	  
	                                          J.-F. Gerbeau, 01/2000
  ----------------------------------------------------------------------*/
/*
  Methods required for the template classes:
  
  *** OpA1 and OpA2:
  
      1) int solve(Vector& x,const Vector& f) const;
         --> solve Ax=f with initial guess x. Return 0 if ok, 1 if ko.

  *** SaddlePoint

      1)   int solveDual(Vector& y,const Vector& x,const Vector& g) const;
         --->  B (\tilde A)^{-1} B^T y = B x0 - g

      2)    int solveDual(Vector& y,const Vector& x) const;
         --->  B (\tilde A)^{-1} B^T y = B x0
	 
  *** OpB: 
  
      1) Vector trans_mul(const Vector& y) const ;
         --> compute and return B^{T}y

*/
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
class InexactLU
{
  const OpA1* _A1;
  const OpA2* _A2;
  const SaddlePoint* _saddle;
  const OpB*  _B;
public:
  InexactLU(const OpA1& A1,const OpA2& A2,const SaddlePoint& saddle,
	    const OpB& B);
  int step1(Vector& x0,const Vector& f ) const;
  int step2(Vector& y, const Vector& x0) const;
  int step2(Vector& y, const Vector& x0, const Vector& g) const;
  int step3(Vector& x, const Vector& x0, const Vector y) const;
  int allSteps(Vector& x,Vector& y,const Vector& f) const;
  int allSteps(Vector& x,Vector& y,const Vector& f,const Vector& g) const;
};

//----------------------------------------------------------------------
//
//                  IMPLEMENTATIONS
//
//----------------------------------------------------------------------
//
//____________________
// constructor
//____________________
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::
InexactLU(const OpA1& A1,const OpA2& A2,const SaddlePoint& saddle,
	  const OpB& B):_A1(&A1),_A2(&A2),_saddle(&saddle),_B(&B)
{}
//____________________
//  step 1 : A1 x0 = f
//____________________
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::step1(Vector& x0,const Vector& f )
  const
{
  return(  _A1->solve(x0,f)  );
}
//___________________________________________
//  step 2 : B \tilde A^{-1} B^T y = B x0 - g
//___________________________________________
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::step2(Vector& y,
						 const Vector& x0,
						 const Vector& g) const
{
  return(  _saddle->solveDual( y , x0 ,g ) );
}
//
// idem when g=0
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::step2(Vector& y,
						 const Vector& x0) const
{
  return(  _saddle->solveDual( y , x0 ) );
}
//__________________________________
//  step 3 : x = x0 +  A2^{-1} B^T y 
//__________________________________
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::
step3(Vector& x, const Vector& x0,const Vector y) const
{
  int i;
  x = 0.;
  i = _A2->solve( x , _B->trans_mul(y) );
  x = x0 - x;
  return(i);
}
//_______________________________
// allSteps : step1, step2, step3 
//_______________________________
//
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::
allSteps(Vector& x,Vector& y,const Vector& f,const Vector& g) const
{
  int i=0;
  Vector x0 = x;
  i += step1(x0,f);
  i += 2*step2(y,x0,g);
  i += 4*step3(x,x0,y);
  return(i);
  /*
    i=0 : everything ok !,
    i=1 : step1 failed,
    i=2 : step2 failed,
    i=4 : step3 failed
    i=3 : step1 & step2 failed,
    i=5 : step1 & step3 failed,
    and so on ... 
  */
}
//__________________
// allSteps when g=0
//------------------
template<class OpA1,class OpA2,class SaddlePoint,class OpB,class Vector>
int InexactLU<OpA1,OpA2,SaddlePoint,OpB,Vector>::
allSteps(Vector& x,Vector& y,const Vector& f) const
{
  int i=0;
  Vector x0 = x;
  i += step1(x0,f);
  i += 2*step2(y,x0);
  i += 4*step3(x,x0,y);
  return(i);
}
//
#endif

