
/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2004-09-23

  Copyright (C) 2004 EPFL

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
#ifndef _GENERALIZEDAITKEN_H
#define _GENERALIZEDAITKEN_H

namespace LifeV
{

template<class Vector, class Real>
class generalizedAitken
{
  /*
  compute the acceleration with the vector variant of Aitken.
  lk,lk1.  are column vectors (lambda_k lambda_{k-1})
  Rk, Rk-1  are  column vectors of residuals(possibly more than one) .
  lk: lambda at the iterate k
  lk1: lambda at the iterate k-1
  Rk: residual(s) associated to lk. in some cases = lk-f(lk)
  Rk1: residual(s) associated to lk1.
  defaultOmega: defaultOmega can only be a scalar.
     If defaultOmega < 0, then do not use Aitken, but fixed relaxation
      parameter=abs(defaultOmega)
     If you want to assign a fixed parameter different from 
     the default one, you can use 
       Vector deltaLambda(lam, res, omega)
     be sure that omega has the size sizeOmega!
  If length(Rk1)==0, then Aitken acceleration is not applicable,
     hence we need a default value for omega.
     !: should be of the size 1,size(Rk,2)
  sizeOmega is the number of relaxation parameters to be computed.
    usually sizeOmega=1, but for example fo r a N-N preconditioner, =2.
  */
private:
  Vector* lk, lk1;
  Vector* Rk, Rk1;
  Real defaultOmega;
  int sizeOmega;
public:
  generalizedAitken(const Real& defOmega=0.1, const int size=1);
  ~generalizedAitken();
  void restart();
  Vector deltaLambda(const Vector& lam, const Vector* res);
  Vector deltaLambda(const Vector& lam, const Vector* res, Real* omega);
  // in this case, omega is taken as the default value

};




template<class Vector, class Real>
generalizedAitken<class Vector, class Real>::generalizedAitken(const Real& defOmega=0.1, const int size=1):
  lk(nil),lk1(nil),Rk(nil),Rk1(nil), sizeOmega(size),
  defaultOmega(defOmega)
{}

template<class Vector, class Real>
generalizedAitken<class Vector, class Real>::~generalizedAitken()
{
  restart();
}

template<class Vector, class Real>
generalizedAitken<class Vector, class Real>::restart()
{
  delete[] lk;
  delete[] lk1;
  delete[] Rk;
  delete[] Rk1;
}


template<class Vector, class Real>
generalizedAitken<class Vector, class Real>::Vector deltaLambda(const Vector& lam, const Vector* residue, Real* omega)
{
  Vector res(lam.size());
  res=0;
  int i;

  if (defaultOmega<=0){
    for (i=0; i<sizeOmega; i++){
      res+= omega[i] * residue[i];
    }
    return(res);
  }
  

  if (lk1 == nil){ // back to default omega
    lk1 = lk;
    Rk1 = Rk;
    lk = new Vector;
    Rk = new Vector[sizeOmega];

    *lk = lam;
    for (i=0; i<sizeOmega; i++){
      *Rk = residue[i];
      res+= omega[i] * Rk[i];
    }
    return(res);
  }

  Vector* temp;
  temp = lk1;  lk1 = lk;  lk = temp;
  temp = Rk1;  Rk1 = lk;  Rk = temp;
  *lk = lam;
  for (i=0; i<sizeOmega; i++){
    *Rk = residue[i];
  }
  
  Vector *resDiff;
  resDiff = new Vector[sizeOmega];
  for (i=0; i<sizeOmega; i++){
    resDiff[i] = Rk[i] - Rk1[i];
  }
  
  if (sizeOmega == 1) {
    omega[0]= - resDiff[0] * (*lk[0] - lk1[0]) / resDiff[0].norm();
  } else {
    cout << "not yet implemented" << endl;
    //  *omega=-pinv(Rk-Rk1)*(lk-lk1);
  }

  for (i=0; i<sizeOmega; i++){
    res+= omega[i] * Rk[i];
  }
  return(res);
}


template<class Vector, class Real>
generalizedAitken<class Vector, class Real>::Vector deltaLambda(const Vector& lam, const Vector* residue)
{
  Real* omega;
  omega = new Real[sizeOmega];
  int i;
  for(i=0; i<sizeOmega; i++){
    omega[i] = abs(defaultOmega);
  }
  return(deltaLambda(lam,residue,omega) );

}

}
#endif
