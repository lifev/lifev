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
   \file ComposedPreconditioner.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2009-03-25
 */

/**
 * \file ComposedPreconditioner.hpp
 * \author Simone Deparis, Paolo Crosetto
 * \date 2009-03-25
 * Class handling a preconditioner defined as a composition of operators. The class is templated
 * so that the operators that define the preconditioner can be general
 * (matrices, other preconditioners, other composed preconditioners, etc.).
 * It is similar to the thyra package in Trilinos although it doesn't depend on it.
 *\end{array}\right)$
 */

#ifndef COMPOSEDPRECONDITIONER_HPP
#define COMPOSEDPRECONDITIONER_HPP

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Epetra_Operator.h>

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/displayer.hpp>

namespace LifeV
{

  template<typename Operator>
class  ComposedPreconditioner
    : public Epetra_Operator {

public:

    /** @name Typedefs
     */
    //@{
    //typedef Ifpack_Preconditioner                prec_raw_type;
    typedef Operator                    prec_raw_type;
    typedef typename boost::shared_ptr<prec_raw_type>     prec_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{
    /**
       The constructor builds an empty chain of composed preconditioner.
       The resulting preconditioner shall be equal to the identity.
    */
      ComposedPreconditioner(const Epetra_Comm* comm=0 );

    /**
       Copy constructor: it doesn't make a real copy of the operators, it just copy the list of shared pointer
       (note that the implementation is mandatory in order to avoid a bit-copy of the shared_ptr vector)
    */
      ComposedPreconditioner( ComposedPreconditioner<Operator>& P);

      ~ComposedPreconditioner( );

      /**
         Method to add an operator at the end of the operators list.
         we expect an external Operator shared pointer, we will not destroy the matrix!
         remember: P = M_P(0) * ... * M_P(Qp-1)
         \param P: the operator
         \param useTranspose: flag specifying if we want to applly the transposed operator
         \param summed: flag specifying if we want the operator to be summed. Note that the operator would be "out
         of the vector", summed in a second step after all the multiplications.
       */
      int push_back(prec_type  P,
                    const bool useInverse=false,
                    const bool useTranspose=false,
                    const bool summed=false);

      /**
         Method to add an operator at the end of the operators list.
         Deprecated. The input parameter P is a standard pointer.
         \param P: the operator
         \param useTranspose: flag specifying if we want to applly the transposed operator
         \param summed: flag specifying if we want the operator to be summed. Note that the operator would be "out
         of the vector", summed in a second step after all the multiplications.
       */

  int push_back(prec_raw_type*  P,
                  const bool useInverse=false,
                const bool useTranspose=false,
                const bool summed=false);

    //! we expecte an external matrix shared pointer, we will not destroy the matrix!
    //! remember: P = M_P(0) * ... * M_P(Qp-1)
      /**
         Method to replace an operator in a specified position of the operators vector.
         \param P: input operator
         \param index: position (starting from 0)
         \param useInverse: flag specifying if the operator is already inverted (if AllpyInverse() is implemented)
       */
    int replace(prec_type&  P,
                UInt const& index,
                const bool useInverse=false,
                const bool useTranspose=false);

    /** @name Epetra_Operator Implementation
     */
    //@{
      /**
         Method to set the flag M_useTranspose, which specifies if the composed preconditioner is to be transposed.
       */
    int SetUseTranspose(bool UseTranspose);

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    // Simone:here use the ApplyInverse from the local preconditioners.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      /** Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\le i\le m} \sum_{j=1}^n |a_{ij}| \f].

       \warning This method must not be called unless HasNormInf() returns true.
     */
    double NormInf() const;

    double Condest() const;

    const char * Label() const;

    bool UseTranspose() const;

    bool HasNormInf() const;

    const Epetra_Comm & Comm() const;

    const Epetra_Map & OperatorDomainMap() const;

    const Epetra_Map & OperatorRangeMap() const;
    //@}

      const std::vector<prec_type>& getP() const {return M_P;}

      //bool getInverse(){return M_allInverse;}

      const bool getAllTranspose() const {return M_allTranspose;}

      const std::vector<bool>  getTranspose() const {return M_Transpose;}

      const std::vector<ID>  getSummed() const {return M_Summed;}

      const std::vector<bool> getInverse() const {return M_Inverse;}

      const int getNumber() const {return M_Setted;}

      const double getMeanIter() const {return M_meanIter;}

      const int getNumCalled() const {return M_numCalled;}

      void reset();

      void swap(UInt first, UInt second);

      void setComm(const Epetra_Comm& comm);

      const Displayer& displayer(){return M_displayer;}

      void resetP(){M_P.clear();}

protected:

      // P = M_P(0) * ... * M_P(n-1)
      mutable std::vector<prec_type> M_P;
      std::vector<bool> M_Inverse;
      std::vector<bool> M_Transpose;
      std::vector<ID> M_Summed;
      bool M_allTranspose;

    mutable int M_Setted;

    mutable double M_meanIter;
    mutable int M_numCalled;


      int Apply_DirectOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      int Apply_InverseOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      Displayer M_displayer;

};

  template <typename Operator>
ComposedPreconditioner<Operator>::ComposedPreconditioner(const Epetra_Comm* comm)
  :
  M_P(),
  M_Inverse(),
  M_Transpose(),
  M_Summed(),
  M_Setted(0),
  M_meanIter(0),
  M_numCalled(0),
  M_displayer(comm)
{
}

  template <typename Operator>
ComposedPreconditioner<Operator>::ComposedPreconditioner( ComposedPreconditioner<Operator>& P)
  :
    M_P(),
    M_Inverse(P.getInverse()),
    M_Transpose(P.getTranspose()),
    M_Summed(P.getSummed()),
    M_allTranspose(P.getAllTranspose()),
    M_Setted(P.getNumber()),
    M_meanIter(P.getMeanIter()),
    M_numCalled(P.getNumCalled()),
    M_displayer(&P.getP()[0]->Comm())
{
    for(int i=0; i<M_Setted; ++i)
        {
            M_P.push_back(P.getP()[i]);
        }
}

  template <typename Operator>
ComposedPreconditioner<Operator>::~ComposedPreconditioner()
{

}

  template <typename Operator>
int ComposedPreconditioner<Operator>::
push_back( prec_raw_type* P,
           const bool useInverse,
           const bool useTranspose,
           const bool summed)
{
  prec_type prec(P);
  push_back(P, useInverse, useTranspose, summed);
}


// we expect an external matrix pointer, we will not destroy the matrix!
  template <typename Operator>
int ComposedPreconditioner<Operator>::
push_back( prec_type P,
           const bool useInverse,
           const bool useTranspose,
           const bool summed)
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
    if(summed)
        M_Summed.push_back(M_Setted);

    //    if (useInverse && useTranspose) assert(false); // I am not sure I have coded this in the right way

    M_Setted++;

    return(M_Setted);

}

// we expect an external matrix pointer, we will not destroy the matrix!
  template <typename Operator>
int ComposedPreconditioner<Operator>::
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

  template <typename Operator>
int ComposedPreconditioner<Operator>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return Apply_DirectOperator(X,Y);

}

  template <typename Operator>
int ComposedPreconditioner<Operator>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return Apply_InverseOperator(X,Y);
}

  template <typename Operator>
int ComposedPreconditioner<Operator>::
Apply_DirectOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (M_Setted > 0 );

  bool useTranspose;
  Epetra_MultiVector Z(X);

  // P = M_P(0) * ... * M_P(Qp-1)

  ID k=0;
  for (int q(M_Setted-1); q>=0 ; q--)
  {
      if(M_Summed.size() && q==M_Summed[k])
      {
          ++k;
      }
      else
      {
          if (M_Inverse[q]) {
              useTranspose = M_P[q]->UseTranspose (); // storing original status
              assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
              M_P[q]->ApplyInverse(Z,Y);
              M_P[q]->SetUseTranspose(useTranspose); // put back to original status
          } else {
              useTranspose = M_P[q]->UseTranspose (); // storing original status
              assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
              M_P[q]->Apply(Z,Y);
              M_P[q]->SetUseTranspose(useTranspose); // put back to original status
          }
          Z = Y;
      }
  }

  for(UInt k=0; k<M_Summed.size(); ++k)
  {
      if(M_Inverse[M_Summed[k]])
          M_P[M_Summed[k]]->ApplyInverse(X,Y);
      else
          M_P[M_Summed[k]]->Apply(X,Y);
      Z.Update(1., Y, 1.);
  }

  Y=Z;
  return(EXIT_SUCCESS);

}

  template <typename Operator>
int ComposedPreconditioner<Operator>::
Apply_InverseOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (M_Setted > 0);

  // P = M_P(0) * ... * M_P(Qp-1)

  bool useTranspose;
  Epetra_MultiVector Z(X);
  ID k=0;

  for (int q(0); q<M_Setted ; q++, Z = Y)
  {
      if(M_Summed.size() && q==M_Summed[k])
      {
          ++k;
      }
      else
      {
          if (M_Inverse[q]) {
              useTranspose = M_P[q]->UseTranspose (); // storing original status
              assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
              M_P[q]->Apply(Z,Y);//y+=Pi*z
              M_P[q]->SetUseTranspose(useTranspose); // put back to original status
          } else {
              useTranspose = M_P[q]->UseTranspose (); // storing original status
              assert (M_P[q]->SetUseTranspose(M_Transpose[q]) != -1);
              M_P[q]->ApplyInverse(Z,Y);//y+=Pi^(-1)*z
              M_P[q]->SetUseTranspose(useTranspose); // put back to original status
          }
          Z=Y;
      }
  }
  Y.Scale(0.);
  for(UInt k=0; k<M_Summed.size(); ++k)
  {
      if(M_Inverse[M_Summed[k]])
      {
          useTranspose = M_P[M_Summed[k]]->UseTranspose (); // storing original status
          assert (M_P[M_Summed[k]]->SetUseTranspose(M_Transpose[M_Summed[k]]) != -1);
          M_P[M_Summed[k]]->Apply(X,Y);
          M_P[M_Summed[k]]->SetUseTranspose(useTranspose); // put back to original status
      }
      else
      {
          useTranspose = M_P[M_Summed[k]]->UseTranspose (); // storing original status
          assert (M_P[M_Summed[k]]->SetUseTranspose(M_Transpose[M_Summed[k]]) != -1);
          M_P[M_Summed[k]]->ApplyInverse(X,Y);
          M_P[M_Summed[k]]->SetUseTranspose(useTranspose); // put back to original status
      }
      Z.Update(1., Y, 1.);
  }

  Y=Z;
  return(EXIT_SUCCESS);
}


  template <typename Operator>
const char * ComposedPreconditioner<Operator>::
Label() const
{
  return("Composed Operator P1*P2* ... *Pn");
}

  template <typename Operator>
bool ComposedPreconditioner<Operator>::
HasNormInf() const
{
  return(false);
}

  template <typename Operator>
const Epetra_Comm & ComposedPreconditioner<Operator>::
Comm() const
{
  assert (M_Setted > 0);
  //cout << "Prec: M_setted = " << M_Setted << endl;
  return( M_P[0]->Comm());
}

  template <typename Operator>
const Epetra_Map & ComposedPreconditioner<Operator>::
OperatorDomainMap() const
{
    // P = M_P(0) * ... * M_P(Qp-1)
  assert (M_Setted > 0 );
  return( M_P[M_Setted-1]->OperatorDomainMap());
}

  template <typename Operator>
const Epetra_Map & ComposedPreconditioner<Operator>::
OperatorRangeMap() const
{
    // P = M_P(0) * ... * M_P(Qp-1)
  assert(M_Setted > 0);
  return( M_P[0]->OperatorRangeMap() );

}

  template <typename Operator>
int ComposedPreconditioner<Operator>::SetUseTranspose(bool useTranspose)
{
    if(useTranspose!=UseTranspose())
        {
            M_allTranspose=useTranspose;
            int size=M_P.size()-1;
            if(size)
            for(int i=0; i<size/2; ++i)
                {
                    prec_type swp;
                    swp=M_P[i];
                    M_P[i]=M_P[size-i];
                    M_P[size-i]=swp;
                    M_P[i]->SetUseTranspose(!M_P[i]->UseTranspose());
                    assert(M_P[size-i]->SetUseTranspose(!M_P[size-i]->UseTranspose())!=-1);
                }
            else
                assert(M_P[0]->SetUseTranspose(!M_P[0]->UseTranspose())!=-1);
            M_P[(int)(size/2)]->SetUseTranspose(!M_P[(int)(size/2)]->UseTranspose());
        }
    return 0;
}

  template <typename Operator>
double ComposedPreconditioner<Operator>::NormInf()  const
{
  assert(false);
  return(EXIT_FAILURE);
}

  template <typename Operator>
bool ComposedPreconditioner<Operator>::UseTranspose()  const
{
    return(M_allTranspose);
}

  template <typename Operator>
double ComposedPreconditioner<Operator>::Condest()  const
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

  template <typename Operator>
  void ComposedPreconditioner<Operator>::reset()
  {
  for (int q(0); q<M_Setted ; q++)
    {
      M_P[q].reset();
    }
  }

template <typename Operator>
void ComposedPreconditioner<Operator>::swap(UInt first, UInt second)
{
    prec_type tmpPrec;
    tmpPrec = M_P[first];
    M_P[first] = M_P[second];
    M_P[second] = tmpPrec;
}

  template <typename Operator> void
  ComposedPreconditioner<Operator>::setComm(const Epetra_Comm& comm)
  {
      //this->M_comm=comm;//copy of the communicator
      M_displayer.SetCommunicator(comm);
  }

} // namespace LifeV
#endif
