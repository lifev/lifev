/* -*- mode: c++ -*- */
//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/**
 * \file ComposedOperator.hpp
 * \author Simone Deparis, Paolo Crosetto
 * \mantainer Paolo Crosetto
 * \date 2009-03-25
 * Class handling a preconditioner defined as a composition of operators. The class is templated
 * so that the operators that define the preconditioner can be general
 * (matrices, other preconditioners, other composed preconditioners, etc.).
 * It is similar to the thyra package in Trilinos although it doesn't depend on it.
 *
 */

#ifndef COMPOSEDPRECONDITIONER_HPP
#define COMPOSEDPRECONDITIONER_HPP


#include <Epetra_Operator.h>
#include <Epetra_MultiVector.h>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

namespace LifeV
{
//! ComposedOperator -
/*!
  @author Simone Deparis, Paolo Crosetto
  @mantainer Paolo Crosetto
  * Class handling a preconditioner defined as a composition of operators. The class is templated
  * so that the operators that define the preconditioner can be general
  * (matrices, other preconditioners, other composed preconditioners, etc.).
  */

template<typename OperatorType>
class  ComposedOperator
    : public Epetra_Operator
{

public:

    //! @name Typedefs
    //@{
    //typedef Ifpack_Preconditioner                       prec_Type;
    typedef OperatorType                                  operator_Type;
    typedef typename std::shared_ptr<operator_Type>     operatorPtr_Type;
    //@}

    //! @name Constructors, destructor
    //@{
    /**
       The constructor builds an empty chain of composed operators.
       The resulting operator shall be equal to the identity.
    */
    ComposedOperator ( const  std::shared_ptr<Epetra_Comm>& comm );

    /**
       Copy constructor: it doesn't make a real copy of the operators, it just copies the vector of shared pointer
       (note that the implementation is mandatory in order to avoid a bit-copy of the shared_ptr vector)
    */
    ComposedOperator ( const ComposedOperator<operator_Type>& P);

    virtual ~ComposedOperator();
    //@}

    //!@name Public Methods
    //@{
    /**
       Method to add an operator at the end of the operators list.
       we expect an external Operator shared pointer, we will not destroy the matrix!
       remember: P = M_operator(0) * ... * M_operator(Qp-1)
       \param P: the operator
       \param useTranspose: flag specifying if we want to applly the transposed operator
       \param summed: flag specifying if we want the operator to be summed. Note that the operator would be "out
       of the vector", summed in a second step after all the multiplications.
    */
    UInt push_back (operatorPtr_Type  P,
                    const bool useInverse = false,
                    const bool useTranspose = false,
                    const bool summed = false);

    /**
       Method to add an operator at the end of the operators list.
       Deprecated. The input parameter P is a standard pointer.
       \param P: the operator
       \param useTranspose: flag specifying if we want to applly the transposed operator
       \param summed: flag specifying if we want the operator to be summed. Note that the operator would be "out
       of the vector", summed in a second step after all the multiplications.
    */

    UInt push_back (operator_Type*  P,
                    const bool useInverse = false,
                    const bool useTranspose = false,
                    const bool summed = false);

    //! we expecte an external matrix shared pointer, we will not destroy the matrix!
    //! remember: P = M_operator(0) * ... * M_operator(Qp-1)
    /**
       Method to replace an operator in a specified position of the operators vector.
       \param P: input operator
       \param index: position (starting from 0)
       \param useInverse: flag specifying if the operator is already inverted (if AllpyInverse() is implemented)
    */
    UInt replace (operatorPtr_Type  P,
                  const UInt index,
                  const bool useInverse = false,
                  const bool useTranspose = false);

    //! Resets all the factors in the operator
    /**
       Calls reset() on all the factors
     */
    void reset();

    //! swaps two factors
    /**
       Swaps the factors and all the related flags (useTranspose, etc...)
       \param first: first factor
       \param second: second factor
     */
    void swap (UInt first, UInt second);

    //! resets the sandard vector holding all the factors
    /**
       calls clear on the std::vector
       \todo change name in e.g. resetOperator
     */
    void resetOperator()
    {
        M_operator.clear();
    }
    //@}

    /** @name Implementation of Epetra_Operator Methods
     */
    //@{
    /**
       Method to set the flag M_useTranspose, which specifies if the composed operator is to be transposed.
    */
    int SetUseTranspose (bool UseTranspose);

    /**
       Method returning the result of applying the operator to a vector
       \param X: input vector
       \param Y: output vector
     */
    int Apply (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    // Simone:here use the ApplyInverse from the local operators.
    int ApplyInverse (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /** Returns the quantity \f$ \| A \|_\infty\f$ such that
        \f[\| A \|_\infty = \max_{1\le i\le m} \sum_{j=1}^n |a_{ij}| \f].

        \warning This method must not be called unless HasNormInf() returns true.
    */
    double NormInf() const;

    //! Returns an estimation of the condition number
    /**
       the estimation is devised from an estimation on each factor with the rule \f$ K_{12}\leq K_1K_2 \f$
    */
    double Condest() const;

    /**
       Method returning a string identifying the type of operator
     */
    const char* Label() const;

    /**
       Returns true if all the operators composed have to be considered transpose (i.e. if M_allTranspose==true)
     */
    bool UseTranspose() const;

    /**
       Never called: it is implemented for compatibility with Epetra_Operator
     */
    bool HasNormInf() const;

    /**
       Returns the communicator
     */
    const Epetra_Comm& Comm() const;

    /**
       Returns a shared pointer to the communicator
     */
    const std::shared_ptr<Epetra_Comm> commPtr() const;

    /**
       returns the OperatorDomainMap of the operators contained
     */
    const Epetra_Map& OperatorDomainMap() const;

    /**
       returns the OperatorRangeMap of the operators contained
     */
    const Epetra_Map& OperatorRangeMap() const;
    //@}

    //! Public Get Methods
    //@{

    //!returns a const reference tor the std::vector of factors
    /**
     */
    const std::vector<operatorPtr_Type>& Operator() const
    {
        return M_operator;
    }

    //!returns the std::vector of factors by (non const) reference
    /**
     */
    std::vector<operatorPtr_Type>& OperatorView() const
    {
        return M_operator;
    }

    //! Returns true if the operator is transposed
    /**
       The operator transposed means that all its factors have to be considered transposed
     */
    bool allTranspose() const
    {
        return M_allTranspose;
    }

    //! Returns a vector saying for each factor if it is considered transposed or not
    const std::vector<bool>&  transpose() const
    {
        return M_transpose;
    }

    //! Returns a vector of UInt specifying the form of the operator
    /**
       The size of the vector corresponds to the number of '+' operation in the operator.
       Each element of the output vector says how many factors are there between two '+' operation. For instance for
       the operator AB+CDE the vector M_summed returned would be [2,3]
     */
    const std::vector<ID>&  summed() const
    {
        return M_summed;
    }

    //! Returns a vector saying which factors have to be considered as inverse
    const std::vector<bool>& inverse() const
    {
        return M_inverse;
    }
    //! returns the number of factors present in the operator
    UInt number() const
    {
        return M_set;
    }

    const double& meanIter() const
    {
        return M_meanIter;
    }

    int numCalled() const
    {
        return M_numCalled;
    }

    const Displayer& displayer()
    {
        return M_displayer;
    }
    //@}

    //!@name Setter Methods
    //!{
    //! Sets the MPI communicator
    void setComm (const std::shared_ptr<Epetra_Comm> comm);
    //!}

protected:

    //!@name Protected Methods
    //!{

    int Apply_DirectOperator (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int Apply_InverseOperator (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
    //@}

    //!@name Protected Members
    //!{
    // P = M_operator(0) * ... * M_operator(n-1)
    mutable std::vector<operatorPtr_Type> M_operator;
    std::vector<bool> M_inverse;
    std::vector<bool> M_transpose;
    std::vector<ID> M_summed;
    bool M_allTranspose;

    UInt M_set;

    mutable double M_meanIter;
    mutable int M_numCalled;

    Displayer M_displayer;
    //@}
};

// ===================================================
//! Constructors and Destructor
// ===================================================

template <typename operator_Type>
ComposedOperator<operator_Type>::ComposedOperator ( const std::shared_ptr<Epetra_Comm>& comm)
    :
    M_operator(),
    M_inverse(),
    M_transpose(),
    M_summed(),
    M_allTranspose (false),
    M_set (0),
    M_meanIter (0),
    M_numCalled (0),
    M_displayer (comm)
{
}

template <typename operator_Type>
ComposedOperator<operator_Type>::ComposedOperator ( const ComposedOperator<operator_Type>& P)
    :
    M_operator(),
    M_inverse (P.inverse() ),
    M_transpose (P.transpose() ),
    M_summed (P.summed() ),
    M_allTranspose (P.allTranspose() ),
    M_set (P.number() ),
    M_meanIter (P.meanIter() ),
    M_numCalled (P.numCalled() ),
    M_displayer (P.commPtr() )
{
    for (UInt i = 0; i < M_set; ++i)
    {
        M_operator.push_back (P.Operator() [i]);
    }
}

template <typename operator_Type>
ComposedOperator<operator_Type>::~ComposedOperator()
{
    M_operator.clear();
}

// ===================================================
//! Public Methods
// ===================================================

template <typename operator_Type>
UInt ComposedOperator<operator_Type>::
push_back ( operator_Type* P,
            const bool useInverse,
            const bool useTranspose,
            const bool summed)
{
    operatorPtr_Type prec (P);
    return push_back (P, useInverse, useTranspose, summed);
}


// we expect an external matrix pointer, we will not destroy the matrix!
template <typename operator_Type>
UInt ComposedOperator<operator_Type>::
push_back ( operatorPtr_Type  P,
            const bool useInverse,
            const bool useTranspose,
            const bool summed)
{
    if (P.get() ==  0)
    {
        M_displayer.leaderPrint (" CP-  Precondtioner not set:                 ", M_set, "\n");
        return (M_set);
    }

    // M_displayer.leaderPrint(" CP-  Previous number of call:                 ", M_numCalled, "\n");
    // M_displayer.leaderPrint(" CP-  Mean iters:                              ", M_meanIter, "\n" );

    M_meanIter = 0;
    M_numCalled = 0;

    // push back the nth operator

    M_operator.        push_back (P);
    M_inverse.  push_back (useInverse);
    M_transpose.push_back (useTranspose);
    if (summed)
    {
        M_summed.push_back (M_set);
    }

    //    if (useInverse && useTranspose) assert(false); // I am not sure I have coded this in the right way

    M_set++;

    return (M_set);
}

// we expect an external matrix pointer, we will not destroy the matrix!
template <typename operator_Type>
UInt ComposedOperator<operator_Type>::
replace ( operatorPtr_Type  P,
          const UInt index,
          const bool useInverse,
          const bool useTranspose)
{
    if (P.get() ==  0)
    {
        return (M_set);
    }

    ASSERT (index <= M_set, "ComposedOperator::replace: index too large");

    M_displayer.leaderPrint (" CP-  Previous number of call:                 ", M_numCalled, "\n");
    M_displayer.leaderPrint (" CP-  Mean iters:                              ", M_meanIter, "\n" );

    M_meanIter = 0;
    M_numCalled = 0;

    // replace the nth operator
    M_operator        [index] = P;
    M_inverse  [index] = useInverse;
    M_transpose[index] = useTranspose;

    if (useInverse && useTranspose)
    {
        assert (false);    // I am not sure I have coded this in the right way
    }

    return (M_set);

}


template <typename operator_Type>
void ComposedOperator<operator_Type>::reset()
{
    for (UInt q (0); q < M_set ; q++)
    {
        M_operator[q].reset();
    }
}


template <typename operator_Type>
void ComposedOperator<operator_Type>::swap (UInt first, UInt second)
{
    operatorPtr_Type tmpPrec;
    tmpPrec = M_operator[first];
    M_operator[first] = M_operator[second];
    M_operator[second] = tmpPrec;
}

// ===================================================
//! Implementation of Epetra_Operator methods
// ===================================================

template <typename operator_Type>
int ComposedOperator<operator_Type>::SetUseTranspose (bool useTranspose)
{
    if (useTranspose != UseTranspose() )
    {
        M_allTranspose = useTranspose;
        int size = M_operator.size() - 1;
        if (size)
            for (int i = 0; i < size / 2; ++i)
            {
                operatorPtr_Type swp;
                swp = M_operator[i];
                M_operator[i] = M_operator[size - i];
                M_operator[size - i] = swp;
                M_operator[i]->SetUseTranspose (!M_operator[i]->UseTranspose() );
                assert (M_operator[size - i]->SetUseTranspose (!M_operator[size - i]->UseTranspose() ) != -1);
            }
        else
        {
            assert (M_operator[0]->SetUseTranspose (!M_operator[0]->UseTranspose() ) != -1);
        }
        M_operator[ (int) (size / 2)]->SetUseTranspose (!M_operator[ (int) (size / 2)]->UseTranspose() );
    }
    return 0;
}

template <typename operator_Type>
int ComposedOperator<operator_Type>::
Apply (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return Apply_DirectOperator (X, Y);

}

template <typename operator_Type>
int ComposedOperator<operator_Type>::
ApplyInverse (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    return Apply_InverseOperator (X, Y);
}

template <typename operator_Type>
double ComposedOperator<operator_Type>::NormInf()  const
{
    assert (false);
    return (EXIT_FAILURE);
}

template <typename operator_Type>
double ComposedOperator<operator_Type>::Condest()  const
{
    double cond (1);

    for (UInt q (0); q < M_set ; q++)
    {
        if (M_inverse[q])
        {
            cond = cond /  M_operator[q]->condest();
        }
        else
        {
            cond *= M_operator[q]->condest();
        }
    }
    return (cond);
}

template <typename operator_Type>
const char* ComposedOperator<operator_Type>::
Label() const
{
    return ("Composed operator_Type P1*P2* ... *Pn");
}

// ===================================================
//! Get Methods
// ===================================================

template <typename operator_Type>
bool ComposedOperator<operator_Type>::UseTranspose()  const
{
    return (M_allTranspose);
}

template <typename operator_Type>
bool ComposedOperator<operator_Type>::
HasNormInf() const
{
    return (false);
}


template <typename operator_Type>
const Epetra_Comm& ComposedOperator<operator_Type>::
Comm() const
{
    assert (M_set > 0);
    //cout << "Prec: M_setted = " << M_set << endl;
    return ( *M_displayer.comm() );
}

template <typename operator_Type>
const std::shared_ptr<Epetra_Comm> ComposedOperator<operator_Type>::
commPtr() const
{
    return ( M_displayer.comm() );
}


template <typename operator_Type>
const Epetra_Map& ComposedOperator<operator_Type>::
OperatorDomainMap() const
{
    // P = M_operator(0) * ... * M_operator(Qp-1)
    assert (M_set > 0 );
    return ( M_operator[M_set - 1]->OperatorDomainMap() );
}

template <typename operator_Type>
const Epetra_Map& ComposedOperator<operator_Type>::
OperatorRangeMap() const
{
    // P = M_operator(0) * ... * M_operator(Qp-1)
    assert (M_set > 0);
    return ( M_operator[0]->OperatorRangeMap() );

}

// ===================================================
//! Set Methods
// ===================================================

template <typename operator_Type> void
ComposedOperator<operator_Type>::setComm (const std::shared_ptr<Epetra_Comm> comm)
{
    //this->M_comm=comm;//copy of the communicator
    M_displayer.setCommunicator (comm);
}

// ===================================================
//! Protected Methods
// ===================================================

template <typename operator_Type>
int ComposedOperator<operator_Type>::
Apply_DirectOperator (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    assert (M_set > 0 );

    bool useTranspose;
    Epetra_MultiVector Z (X);

    // P = M_operator(0) * ... * M_operator(Qp-1)

    ID k = 0;
    for (Int q (M_set - 1); q >= 0 ; q--)
    {
        if (M_summed.size() && q == (Int) M_summed[k])
        {
            ++k;
        }
        else
        {
            if (M_inverse[q])
            {
                useTranspose = M_operator[q]->UseTranspose (); // storing original status
                assert (M_operator[q]->SetUseTranspose (M_transpose[q]) != -1);
                M_operator[q]->ApplyInverse (Z, Y);
                M_operator[q]->SetUseTranspose (useTranspose); // put back to original status
            }
            else
            {
                useTranspose = M_operator[q]->UseTranspose (); // storing original status
                assert (M_operator[q]->SetUseTranspose (M_transpose[q]) != -1);
                M_operator[q]->Apply (Z, Y);
                M_operator[q]->SetUseTranspose (useTranspose); // put back to original status
            }
            Z = Y;
        }
    }

    for (UInt k = 0; k < M_summed.size(); ++k)
    {
        if (M_inverse[M_summed[k]])
        {
            M_operator[M_summed[k]]->ApplyInverse (X, Y);
        }
        else
        {
            M_operator[M_summed[k]]->Apply (X, Y);
        }
        Z.Update (1., Y, 1.);
    }

    Y = Z;
    return (EXIT_SUCCESS);

}

template <typename operator_Type>
int ComposedOperator<operator_Type>::
Apply_InverseOperator (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
    assert (M_set > 0);

    // P = M_operator(0) * ... * M_operator(Qp-1)

    bool useTranspose;
    Epetra_MultiVector Z (X);
    ID k = 0;

    for (UInt q (0); q < M_set ; q++, Z = Y)
    {
        if (M_summed.size() && q == M_summed[k])
        {
            ++k;
        }
        else
        {
            if (M_inverse[q])
            {
                useTranspose = M_operator[q]->UseTranspose (); // storing original status
                ASSERT (M_operator[q]->SetUseTranspose (M_transpose[q]) != -1, "SetUseTranspose failed in the preconditioner of your choice");
                M_operator[q]->Apply (Z, Y); //y+=Pi*z
                M_operator[q]->SetUseTranspose (useTranspose); // put back to original status
            }
            else
            {
                useTranspose = M_operator[q]->UseTranspose (); // storing original status
                //ifdef DEBUG
                ASSERT (!M_transpose[q] || M_operator[q]->SetUseTranspose (M_transpose[q]) == -1, "Something went wrong in SetUseTranspose for your preconditioner.\n");
                //#endif
                M_operator[q]->ApplyInverse (Z, Y); //y+=Pi^(-1)*z
                M_operator[q]->SetUseTranspose (useTranspose); // put back to original status
            }
            Z = Y;
        }
    }
    Y.Scale (0.);
    for (UInt k = 0; k < M_summed.size(); ++k)
    {
        if (M_inverse[M_summed[k]])
        {
            useTranspose = M_operator[M_summed[k]]->UseTranspose (); // storing original status
            assert (M_operator[M_summed[k]]->SetUseTranspose (M_transpose[M_summed[k]]) != -1);
            M_operator[M_summed[k]]->Apply (X, Y);
            M_operator[M_summed[k]]->SetUseTranspose (useTranspose); // put back to original status
        }
        else
        {
            useTranspose = M_operator[M_summed[k]]->UseTranspose (); // storing original status
            assert (M_operator[M_summed[k]]->SetUseTranspose (M_transpose[M_summed[k]]) != -1);
            M_operator[M_summed[k]]->ApplyInverse (X, Y);
            M_operator[M_summed[k]]->SetUseTranspose (useTranspose); // put back to original status
        }
        Z.Update (1., Y, 1.);
    }

    Y = Z;
    return (EXIT_SUCCESS);
}

} // namespace LifeV
#endif
