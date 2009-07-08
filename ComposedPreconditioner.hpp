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


#ifndef COMPOSEDPRECONDITIONER_HPP
#define COMPOSEDPRECONDITIONER_HPP

#include <boost/shared_ptr.hpp>
#include <vector>
#include <Epetra_Operator.h>

#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>

#include <Ifpack_Preconditioner.h>

namespace LifeV
{


class  ComposedPreconditioner
    : public Epetra_Operator {

public:

    /** @name Typedefs
     */
    //@{
    typedef Ifpack_Preconditioner                prec_raw_type;
    //typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{
    /**
       The constructor builds an empty chain of composed preconditioner.
       The resulting preconditioner shall be equal to the identity.
    */
    ComposedPreconditioner( );

    ~ComposedPreconditioner( );

    //! we expecte an external matrix shared pointer, we will not destroy the matrix!
    //! remember: P = M_P(0) * ... * M_P(Qp-1)
    int push_back(prec_type&  P,
                  const bool useInverse=false,
                  const bool useTranspose=false);

    //! we expecte an external matrix shared pointer, we will not destroy the matrix!
    //! remember: P = M_P(0) * ... * M_P(Qp-1)
    int replace(prec_type&  P,
                UInt const& index,
                const bool useInverse=false,
                const bool useTranspose=false);

    /** @name Epetra_Operator Implementation
     */
    //@{
    int SetUseTranspose(bool UseTranspose);

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    // Simone:here use the ApplyInverse from the local preconditioners.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /* Returns the quantity \f$ \| A \|_\infty\f$ such that
       \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

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



protected:

    // P = M_P(0) * ... * M_P(n-1)
    mutable std::vector<prec_type> M_P;
    std::vector<bool> M_Inverse;
    std::vector<bool> M_Transpose;

    mutable int M_Setted;

    mutable double M_meanIter;
    mutable int M_numCalled;


    int Apply_DirectOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    int Apply_InverseOperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


};

} // namespace LifeV
#endif
