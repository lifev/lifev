/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#ifndef _SADDLEPOINTCG_H
#define _SADDLEPOINTCG_H
#include <cg.h>

namespace LifeV
{
/*----------------------------------------------------------------------

  A class to deal with the following saddle point problem:

  [ A  B^T ] [ x ]    [ f ]
  [        ] [   ]  = [   ]
  [ B  0   ] [ y ]    [ g ]


  The Schur complement (dual problem) B A^{-1} B^T y = B A^{-1} f  -  g
  is solved by preconditionned conjugate gradient.

	                                          J.-F. Gerbeau, 01/2000
  ----------------------------------------------------------------------*/
/*
  Methods required for the template classes:

  *** OpA:

      1) Vector operator*(const Vector& u) const ;
         --> compute and return A*x

      2) int solve(Vector& x,const Vector& f) const;
         --> solve Ax=f with initial guess x. Return 0 if ok, 1 if ko.

      3) Vector solve(const Vector& f) const;
         --> solve Ax=f and return x

  *** OpB:

      1) Vector operator*(const Vector& x) const ;
         --> compute and return Bx
      2) Vector trans_mul(const Vector& y) const ;
         --> compute and return B^{T}y

  *** Prec: (it is the preconditionner class for the Schur complement problem)

      1) Vector solve(const Vector& f) const;

*/

template<class OpA,class OpB,class Prec,class Vector>
class SaddlePointCG
{
  const OpA* _A;
  const OpB* _B;
  const Prec* _P; // preconditionner for the dual problem
  double _tol;    // tolerance for the resolution of the dual problem
  int _max_iter;  // max iter for the resolution of the dual problem
public:
  SaddlePointCG(const OpA& A,const OpB& B,const Prec& P,double tol,
		int max_iter);
  inline const OpA* operatorA() {return _A;}
  inline const OpB* operatorB() {return _B;}
  Vector operator*(const Vector& y) const; // B A^{-1} B^T y
  Vector rhsDual(const Vector& f, const Vector& g) const;// B A^{-1} f - g
  Vector rhsDual(const Vector& f                 ) const;// B A^{-1} f
  int solvePrimal(Vector& x,const Vector& y,
		  const Vector& f) const;        // A x = f - B^T y
  int solveDualRHS(Vector& y,const Vector& b) const;// B A^{-1} B^T y = b
  int solveDual(Vector& y,const Vector& x) const;// B A^{-1} B^T y = Bx
  int solveDual(Vector& y,const Vector& x,
		const Vector& g) const; // B A^{-1} B^T y = Bx - g
};
//----------------------------------------------------------------------
//
//                  IMPLEMENTATIONS
//
//----------------------------------------------------------------------
//----------------------------------------
// constructor
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
SaddlePointCG<OpA,OpB,Prec,Vector>::SaddlePointCG(const OpA& A,const OpB& B,
					      const Prec& P,double tol,
					      int max_iter):
  _A(&A),_B(&B),_P(&P),_tol(tol),_max_iter(max_iter)
{}

//----------------------------------------
//  the product B A^{-1} B^T y
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
Vector SaddlePointCG<OpA,OpB,Prec,Vector>::operator*(const Vector& y) const
{
  return( (*_B) * ( _A->solve( _B->trans_mul(y) ) ) );
}

//----------------------------------------
// compute B A^{-1} f - g
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
Vector SaddlePointCG<OpA,OpB,Prec,Vector>::rhsDual(const Vector& f,
						   const Vector& g) const
{
  return( (*_B) * ( _A->solve( f ) ) - g);

}

//----------------------------------------
// compute B A^{-1} f
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
Vector SaddlePointCG<OpA,OpB,Prec,Vector>::rhsDual(const Vector& f) const
{
  return( (*_B) * ( _A->solve( f ) ));
}

//----------------------------------------
//  solve Ax = f - B^T y
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
int SaddlePointCG<OpA,OpB,Prec,Vector>::solvePrimal(Vector& x,const Vector& y,
						    const Vector& f) const
{
  return( _A->solve( x, f - _B->trans_mul(y) ) );
}

//----------------------------------------
// solve  B A^{-1} B^T y = f
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
int SaddlePointCG<OpA,OpB,Prec,Vector>::solveDualRHS(Vector& y,
						     const Vector& b) const
{
  int i, maxit = _max_iter;
  double tol = _tol;
  i = CG( (*this) , y , b , (*_P) , maxit , tol);
  if(i)
    cout << "SaddlePointCG::solveDualRHS --> Convergence failed ! max iter = "
	 << maxit << ", tolerance = " << tol << endl;
  else
    cout << "SaddlePointCG::solveDualRHS --> iter = "
	 << maxit << ", tolerance = " << tol << endl;
  return i;
}
//----------------------------------------
// solve  B A^{-1} B^T y = B x - g
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
int SaddlePointCG<OpA,OpB,Prec,Vector>::solveDual(Vector& y,
						  const Vector& x,
						  const Vector& g) const
{
  int i, maxit = _max_iter;
  double tol = _tol;
  i = CG( (*this) , y , (*_B) * x - g , (*_P) , maxit , tol);
  if(i)
    cout << "SaddlePointCG::solveDual --> Convergence failed ! max iter = "
	 << maxit << ", tolerance = " << tol << endl;
  else
    cout << "SaddlePointCG::solveDual --> iter = " << maxit << ", tolerance = " << tol << endl;
  return i;
}


//----------------------------------------
// solve  B A^{-1} B^T y = B x
//----------------------------------------
template<class OpA,class OpB,class Prec,class Vector>
int SaddlePointCG<OpA,OpB,Prec,Vector>::solveDual(Vector& y,
						  const Vector& x) const
{
  int i, maxit = _max_iter;
  double tol = _tol;
  i = CG( (*this) , y , (*_B) * x, (*_P) , maxit , tol);
  if(i)
    cout << "SaddlePointCG::solveDual --> Convergence failed ! max iter = "
	 << maxit <<", tolerance = " << tol << endl;
  else
    cout << "SaddlePointCG::solveDual --> iter = " << maxit
	 << ", tolerance = " << tol << endl;
  return i;
}
}
#endif


