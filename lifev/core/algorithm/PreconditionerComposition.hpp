//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the PreconditionerComposition class.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 30-11-2010

    The class provides an efficient way of dealing with preconditioner
    composed as a multiplication of matrices.
 */

#ifndef PRECONDITIONERCOMPOSITION_HPP
#define PRECONDITIONERCOMPOSITION_HPP 1

#include <vector>
#include <boost/shared_ptr.hpp>

#include <Teuchos_ParameterList.hpp>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/ComposedOperator.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorBlockStructure.hpp>

namespace LifeV
{

//! PreconditionerComposition - Class to manage preconditioners composed by matrices multiplication.
/*!
    @author Gwenol Grandperrin

    This class makes use of ComposedOperator to handle matrices composition.
 */
class PreconditionerComposition:
    public Preconditioner
{
public:

    //! @name Public Types
    //@{
    typedef Preconditioner                         super_Type;
    typedef std::shared_ptr<super_Type>            superPtr_Type;
    typedef Epetra_Operator                        operator_Type;
    typedef std::shared_ptr<operator_Type>       operatorPtr_Type;
    typedef ComposedOperator<operator_Type>        prec_Type;
    typedef std::shared_ptr<prec_Type>           precPtr_Type;

    typedef MatrixEpetra<Real>                     matrix_Type;
    typedef std::shared_ptr<matrix_Type>           matrixPtr_Type;

    typedef Teuchos::ParameterList                 list_Type;
    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    /*!
        @param comm Communicator (std::shared_ptr<Epetra_Comm>() by default)
     */
#ifdef HAVE_MPI
    PreconditionerComposition ( std::shared_ptr<Epetra_Comm> comm = std::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerComposition ( std::shared_ptr<Epetra_Comm> comm = std::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

private:

    //! Copy constructor
    /*!
        @param precComp PreconditionerComposition
     */
    PreconditionerComposition ( const PreconditionerComposition& precComp );

public:
    //! Destructor
    ~PreconditionerComposition();

    //@}

    //! @name Methods
    //@{

    virtual void createParametersList ( list_Type& list,
                                        const GetPot& dataFile,
                                        const std::string& section,
                                        const std::string& subSection ) = 0;

    //! Build the preconditioner
    /*!
      @param A the base matrix for computing the preconditioner
    */
    virtual int buildPreconditioner ( matrixPtr_Type& A ) = 0;

    //! Reset the preconditioner
    void resetPreconditioner();

    //! Return an estimation of the conditionement number of the preconditioner
    Real condest();

    //@}

    //! @name Epetra Operator Interface Methods
    //@{

    Int SetUseTranspose ( const bool useTranspose = false );

    Int Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    Int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    bool UseTranspose();

    const Epetra_Map& OperatorRangeMap() const;

    const Epetra_Map& OperatorDomainMap() const;

    //@}

    //! @name Set Methods
    //@{

    //! Setter using GetPot
    /*!
        This method use GetPot to load data from a file and then set
        the preconditioner.
        @param dataFile is a GetPot dataFile
        @param section is the section containing the data
     */
    virtual void setDataFromGetPot ( const GetPot& dataFile,
                                     const std::string& section ) = 0;

    //! Method to setup the solver using Teuchos::ParameterList
    /*!
        @param list Teuchos::ParameterList object
     */
    virtual void setParameters ( Teuchos::ParameterList& list ) = 0;

    /*!
        copies the shared_ptr to the communicator in the member M_comm and builds a new instance
    */
    void setComm ( std::shared_ptr<Epetra_Comm> comm );

    //@}

    //! @name Get Methods
    //@{

    //! Preconditioner is set?
    /*!
     *  @return true
     */
    bool isPreconditionerSet() const;

    /** Get a standard pointer to the preconditioner. In most of the cases is more safe to use getPrecPtr(), which
     returns a std::shared_ptr*/
    operator_Type* preconditioner();

    /** get a std::shared_ptr to the preconditioner. The only requirement on the preconditioner is that
     it must derive from the Epetra_Operator object*/
    operatorPtr_Type preconditionerPtr();

    //! Return the type name of the preconditioner.
    /*!
     *  @return type of the preconditioner
     */
    virtual std::string preconditionerType();

    //! Return the number of operators in the composition
    UInt numOperators() const;

    //@}

protected:

    //! @name Private Methods
    //@{

    //! Add A to the right of the composition
    int pushBack ( matrixPtr_Type A,
                   const bool useInverse   = false,
                   const bool useTranspose = false );

    int pushBack ( operatorPtr_Type oper,
                   const bool useInverse     = false,
                   const bool useTranspose   = false,
                   matrixPtr_Type baseMatrix = matrixPtr_Type() );

    //! Use a preconditioner to build the inverse of A and add it to the right of the composition
    int pushBack ( matrixPtr_Type A,
                   superPtr_Type preconditioner,
                   const bool useInverse   = false,
                   const bool useTranspose = false );

    //! Use a preconditioner to build the inverse of A and add it to the right of the composition
    int pushBack ( matrixPtr_Type embeddedA,
                   superPtr_Type preconditioner,
                   const VectorBlockStructure& blockStructure,
                   const UInt& blockIndex,
                   const MapEpetra& fullMap,
                   const bool useInverse   = false,
                   const bool useTranspose = false,
                   const bool buildPreconditioner = true );

    int pushBack ( operatorPtr_Type embeddedOperator,
                   const VectorBlockStructure& blockStructure,
                   const UInt& blockIndex,
                   const MapEpetra& fullMap,
                   const bool useInverse,
                   const bool useTranspose );

    //@}

    std::shared_ptr<Epetra_Comm> M_comm;

private:
    precPtr_Type                   M_prec;
    std::vector< matrixPtr_Type >  M_precBaseOperators;
};

} // Namespace LifeV

#endif /* PRECONDITIONERCOMPOSITION_HPP */
