/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2009-03-25

  Copyright (C) 2009 EPFL

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
   \file ComposedPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2009-03-25
 */

#include "ComposedPreconditioner.hpp"
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

namespace LifeV
{

ComposedPreconditioner::ComposedPreconditioner()
  :
  M_P(),
  M_Inverse(0),
  M_Transpose(0),
  M_Setted(0),
  M_meanIter(0),
  M_numCalled(0)
{
}

ComposedPreconditioner::~ComposedPreconditioner()
{

}

// we expecte an external matrix pointer, we will not destroy the matrix!
int ComposedPreconditioner::
push_back( prec_type& P,
           const bool useInverse,
           const bool useTranspose)
{
    if (P.get() ==  0)   return(M_Setted);

    if (P->Comm().MyPID() == 0 && M_numCalled > 0)
        cout << "Previous number of call: "
             << M_numCalled << ", mean iters: " << M_meanIter << endl;

    M_meanIter=0;
    M_numCalled=0;

    // push back the nth operator
    M_P.        push_back(P);
    M_Inverse.  push_back(useInverse);
    M_Transpose.push_back(useTranspose);

    if (useInverse && useTranspose) assert(false); // I am not sure I have coded this in the right way

    M_Setted++;

    return(M_Setted);

}

// we expecte an external matrix pointer, we will not destroy the matrix!
int ComposedPreconditioner::
replace( prec_type& P,
         UInt const& index,
         const bool useInverse,
         const bool useTranspose)
{
    if (P.get() ==  0)   return(M_Setted);

    ASSERT(index <= M_Setted, "ComposedPreconditioner::replace: index too large");

    if (P->Comm().MyPID() == 0 && M_numCalled > 0)
        cout << "Previous number of call: "
             << M_numCalled << ", mean iters: " << M_meanIter << endl;

    M_meanIter=0;
    M_numCalled=0;

    // replace the nth operator
    M_P        [index] = P;
    M_Inverse  [index] = useInverse;
    M_Transpose[index] = useTranspose;

    if (useInverse && useTranspose) assert(false); // I am not sure I have coded this in the right way

    return(M_Setted);

}

int ComposedPreconditioner::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return Apply_DirectOperator(X,Y);

}

int ComposedPreconditioner::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return Apply_InverseOperator(X,Y);
}


int ComposedPreconditioner::
Apply_DirectOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (M_Setted > 0 );

  bool useTranspose;
  Epetra_MultiVector Z(X);

  // P = M_P(0) * ... * M_P(Qp-1)

  for (int q(M_Setted-1); q>=0 ; q--) {
    if (M_Inverse[q]) {
        useTranspose = M_P[q]->UseTranspose (); // storing original status
        assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
        M_P[q]->ApplyInverse(Z,Y);
        M_P[q]->SetUseTranspose(useTranspose); // put back to original status
    } else {
        useTranspose = M_P[q]->UseTranspose (); // storing original status
        M_P[q]->Apply(Z,Y);
        M_P[q]->SetUseTranspose(useTranspose); // put back to original status
    }

    Z = Y;
  }

  return(EXIT_SUCCESS);

}

int ComposedPreconditioner::
Apply_InverseOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (M_Setted > 0);

  // P = M_P(0) * ... * M_P(Qp-1)

  bool useTranspose;
  Epetra_MultiVector Z(X);

  for (int q(0); q<M_Setted ; q++) {
    if (M_Inverse[q]) {
        useTranspose = M_P[q]->UseTranspose (); // storing original status
        M_P[q]->Apply(Z,Y);
        M_P[q]->SetUseTranspose(useTranspose); // put back to original status
    } else {
        useTranspose = M_P[q]->UseTranspose (); // storing original status
        assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
        M_P[q]->ApplyInverse(Z,Y);
        M_P[q]->SetUseTranspose(useTranspose); // put back to original status
    }

    Z = Y;
  }

  return(EXIT_SUCCESS);

}


const char * ComposedPreconditioner::
Label() const
{
  return("Composed Operator P1*P2* ... *Pn");
}

bool ComposedPreconditioner::
HasNormInf() const
{
  return(false);
}

const Epetra_Comm & ComposedPreconditioner::
Comm() const
{
  assert (M_Setted > 0);
  //cout << "Prec: M_setted = " << M_Setted << endl;
  return( M_P[0]->Comm());
}

const Epetra_Map & ComposedPreconditioner::
OperatorDomainMap() const
{
    // P = M_P(0) * ... * M_P(Qp-1)
  assert (M_Setted > 0 );
  return( M_P[M_Setted-1]->OperatorDomainMap());
}

const Epetra_Map & ComposedPreconditioner::
OperatorRangeMap() const
{
    // P = M_P(0) * ... * M_P(Qp-1)
  assert(M_Setted > 0);
  return( M_P[0]->OperatorRangeMap() );

}

int ComposedPreconditioner::SetUseTranspose(bool /*UseTranspose*/)
{
  assert(false);
  return(EXIT_FAILURE);
}

double ComposedPreconditioner::NormInf()  const
{
  assert(false);
  return(EXIT_FAILURE);
}


double ComposedPreconditioner::Condest()  const
{
    double cond(1);

    for (int q(0); q<M_Setted ; q++) {
        if (M_Inverse[q])
            cond = cond /  M_P[q]->Condest();
        else
            cond *= M_P[q]->Condest();
    }
    return(cond);
}

bool ComposedPreconditioner::UseTranspose()  const
{
  return(false);
}


} // namespace LifeV








