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

/*!
 *  @file
 *  @brief File containing the Composed Preconditioner Class
 *
 *  @date 07-05-2009
 *  @author Simone Deparis <simone.deparis@epfl.ch>, Paolo Crosetto <paolo.crosetto@epfl.ch>
 *
 *  @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @maintainer Paolo Crosetto <paolo.crosetto@epfl.ch>
 * Class handling a preconditioner defined as a multiplication of several preconditioners. The class inherits
 * from the base Preconditioner class and implements some methds which are used in the EpetraOperator, so that it can be
 * passed as a template argument to the ComposedOperator object.
 */

#ifndef PreconditionerComposed_HPP
#define PreconditionerComposed_HPP

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/ComposedOperator.hpp>

namespace LifeV
{

//! PreconditionerComposed -
/*!
  @author Simone Deparis, Paolo Crosetto
  @mantainer Paolo Crosetto
  * Class handling a preconditioner defined as a multiplication of several preconditioners. The class inherits
  * from the base Preconditioner class and implements some methds which are used in the EpetraOperator, so that it can be
  * passed as a template argument to the ComposedOperator object. It contains an instance of ComposedOperator, which will be filled with
  * shared pointers to trilinos preconditioners (more generally Epetra_Operators)
  * and a vector of shared pointers to LifeV matrices (of MatrixEpetra type, in M_operVector). Notice that the matrices are not copied when the
  * push_back method is called, and not even if the
  * copy constructor is called. Only the shared pointers are copied. The same happens to the preconditioners through the constructor of
  * the ComposedOperator class, which behave in the same way.
  */
class PreconditionerComposed:
    public Preconditioner
{
public:

    /** @name Typedefs
     */
    //@{

    typedef Preconditioner                                             super_Type;
    typedef ComposedOperator<Preconditioner>                           prec_Type;
    typedef std::shared_ptr<prec_Type>                               precPtr_Type;
    typedef std::shared_ptr<Preconditioner>                          epetraPrecPtr_Type;
    typedef super_Type::operator_raw_type                              operator_Type;
    typedef std::shared_ptr<operator_Type>                             operatorPtr_Type;
    typedef super_Type::list_Type                                      list_Type;
    //typedef std::shared_ptr<Preconditioner>                          epetraPrec_Type;


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    PreconditionerComposed ( std::shared_ptr<Epetra_Comm> comm = std::shared_ptr<Epetra_Comm>() );

private: //set it private to be sure that it's never called
    //!Copy Constructor
    /**
       This copy constructor does not copy the matrices, but just the shared_ptrs. It calls the copy constructor of ComposedOperator.
     */
    PreconditionerComposed ( PreconditionerComposed& P );

public:
    //! constructor from matrix A.
    //! @param A MatrixEpetra<double> matrix upon which construct the preconditioner
    //    PreconditionerComposed(operatorPtr_Type& A);

    //! default destructor

    virtual ~PreconditionerComposed();

    //@}


    //!@name  Public Methods
    //@{
    //! Sets the data from GetPot
    void                   setDataFromGetPot ( const GetPot&      dataFile,
                                               const std::string& section);
    //! Creates the Trilinos Teuchos parameter list reading from data file
    void                   createParametersList ( list_Type& /*list*/, const GetPot& dataFile, const std::string& section, const std::string& subSection );

    //! Returns an estimate of the condition number
    double                 condest ();

    //! same as push_back
    int                    buildPreconditioner (operatorPtr_Type& A);

    //! same as push_back
    int                    buildPreconditioner (operatorPtr_Type& A,
                                                const bool useInverse,
                                                const bool useTranspose = false);
    //! Builds a preconditioner based on A and pushes it back in the composedPreconditioner.
    int                    push_back          (operatorPtr_Type& A,
                                               const bool useInverse = false,
                                               const bool useTranspose = false
                                              );

    //! Builds a preconditioner based on A and replaces it in the composedPreconditioner.
    int                    replace            (operatorPtr_Type& A,
                                               const UInt index,
                                               const bool useInverse = false,
                                               const bool useTranspose = false);

    //! resets the pointer to the preconditioner M_prec
    void                   resetPreconditioner();

    //! returns the operator vectir
    const std::vector<operatorPtr_Type>& operVector() const
    {
        return M_operVector;
    }
    //@}

    //!@name Implementation of Methods from Epetra_Operator
    //@{
    //! returns the communicator
    const Epetra_Comm& Comm()
    {
        return preconditioner()->Comm();
    }

    //! sets the M_useTranspose flag
    int            SetUseTranspose ( bool useTranspose = false )
    {
        return M_prec->SetUseTranspose (useTranspose);
    }
    //! returns the M_useTranspose flag
    bool            UseTranspose(  )
    {
        return M_prec->UseTranspose();
    }

    //! Applies the inverse operator to an input vector
    /**
       \param X: input vector
       \param Y: output vector
     */
    virtual int ApplyInverse (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_prec->ApplyInverse (X, Y);
    }

    //! Applies the operator to an input vector
    /**
       \param X: input vector
       \param Y: output vector
     */
    virtual int Apply (const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
        return M_prec->Apply (X, Y);
    }

    //! returns the range map
    const Epetra_Map& OperatorRangeMap() const
    {
        return M_prec->OperatorRangeMap();
    }

    //! returns the domain map
    const Epetra_Map& OperatorDomainMap() const
    {
        return M_prec->OperatorDomainMap();
    }
    //@}

    //!@name Get Methods
    //@{
    //! returns a raw pointer to the preconditioner base class
    super_Type::prec_raw_type*  preconditioner ();

    //! returms the number of factors in the preconditioner
    UInt number() const
    {
        return M_prec->number();
    }

    //! returns a shared pointer to the preconditioner
    super_Type::prec_type              preconditionerPtr()
    {
        return M_prec;
    }

    //! returns a shared pointer to the preconditioner
    const precPtr_Type                 composedPreconditionerPtr()
    {
        return M_prec;
    }

    //! returns a string identifying the preconditioner type
    std::string            preconditionerType()
    {
        return "PreconditionerComposed";
    }
    //@}

    //!@name Static Methods
    //@{

    //! Factory method
    static Preconditioner* createComposedPreconditioner()
    {
        return new PreconditionerComposed();
    }
    //@}

private:

    //!@name Private Methods
    //@{
    void myCreateParametersList (const GetPot& dataFile, const std::string& section, const std::string& subSection);

    Int createPrec (operatorPtr_Type& oper,
                    std::shared_ptr<Preconditioner>& prec);
    //@}


    //!@name Private Members
    //@{
    std::vector<operatorPtr_Type> M_operVector; // we need to keep track of all the operators.
    precPtr_Type                  M_prec;
    static bool registerComposed;
    //@}
};

} // namespace LifeV

#endif // PreconditionerComposed_HPP
