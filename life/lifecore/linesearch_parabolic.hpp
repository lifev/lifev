#ifndef _LINESEARCH_PARAB_H_
#define _LINESEARCH_PARAB_H_

/*
  This line search algorithm comes from C.T. Kelley, Iterative methods for linear
  and nonlinear equations, SIAM 1995 (Chap. 8).
  
    (i)   lambda given (usually 1 when Newton method is used)
    (ii)  sol_test = sol + lambda step_test
    (iii) if residu_test < (1 - alpha lambda) residu
          then sol = sol_test
	  else choose lambda via a three point parabolic interpolation
	       (that does not net the derivative) and apply a safeguarding
	       step: if lambda < sigma0 lambda_cur then lambda = sigma0 lambda_cur
                     if lambda > sigma1 lambda_cur then lambda = sigma1 lamnda_cur

    Constant parameters:
    
    sigma0, sigma1: safeguarding bounds (default values 0.1 and 0.5)
    alpha         : parameter to measure sufficient decrease (default 1e-4)
    max_linesearch: maximum number of steplength reductions before
                    failure is reported (default 50)
    
    
  */

template<class Fct,class Vector,class Real, class Norm>
void lineSearch_parab(Fct& f,Norm& norm, Vector& residual,Vector& sol,Vector& step,Real& normRes,
		 Real& lambda,int iter)
{  
  //----------------------------------------------------------------------
  const Real sigma0         = 0.1;
  const Real sigma1         = 0.5;
  const Real alpha          = 1.e-4;
  const int  max_linesearch = 50;
  //----------------------------------------------------------------------
  static Vector sol_cur = sol; // static to avoid too many alloc/de-alloc
  int iter_linesearch;
  Real lambda_cur,lambda_old,normRes_test,res2,res_test2,res_test_old2,c1,c2;
  //
  res2=normRes*normRes;
  lambda_old=lambda;
  lambda_cur=lambda;
  sol_cur = sol;
  sol += lambda * step;
  f.evalResidual(residual,sol,iter);
  normRes_test = norm(residual);
  res_test2 = normRes_test * normRes_test;
  res_test_old2 = res_test2;
  iter_linesearch = 0;
  while( normRes_test >= (1-alpha*lambda)*normRes ){
    iter_linesearch++;
    // parabolic interpolation of lambda:
    c2 = lambda_old*(res_test2 - res2) - lambda_cur*(res_test_old2 - res2);
    if(c2 >= 0) lambda = sigma1*lambda_cur;
    else{
      c1 = lambda_cur*lambda_cur*(res_test_old2 - res2)
	- lambda_old*lambda_old*(res_test2 - res2);
      lambda = -c1*.5/c2;
      if(lambda < sigma0*lambda_cur) lambda = sigma0*lambda_cur;
      if(lambda > sigma1*lambda_cur) lambda = sigma1*lambda_cur;
    }
    cout << "--- line search " << iter_linesearch << " : residual test = "
	 << normRes_test << ", reduction = " << lambda << endl;
    // update sol_test
    sol = sol_cur + lambda * step;
    lambda_old = lambda_cur;
    lambda_cur = lambda;
    // eval norms
    f.evalResidual(residual,sol,iter);
    normRes_test = norm(residual);
    res_test_old2 = res_test2;
    res_test2 = normRes_test*normRes_test;
    if(iter_linesearch>max_linesearch){
      cout << "!!! Too many iterations in the line search algorithm" << endl;
      exit(1);
    }
  }
  normRes = normRes_test;
}
#endif
