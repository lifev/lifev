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
#include <life/lifecore/Life.hpp>
#include <life/lifealg/Preconditioner.hpp>
#include "ComposedOperator.hpp"

namespace LifeV {

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

    typedef Preconditioner                         super;
    typedef ComposedOperator<Preconditioner>       prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>       prec_type;
    typedef super::operator_raw_type               operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type>   operator_type;
    typedef Teuchos::ParameterList                 list_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    /*!
        @param comm Communicator (boost::shared_ptr<Epetra_Comm>() by default)
     */
    PreconditionerComposition( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm>() );

    //! Copy constructor
    /*!
        @param precComp PreconditionerComposition
     */
    PreconditionerComposition( const PreconditionerComposition& precComp );

    //! Destructor
    ~PreconditionerComposition();

    //@}



    //! @name Methods
    //@{

    void createList( list_type& list,
                     const GetPot& dataFile,
                     const std::string& section,
                     const std::string& subSection ) = 0;

    //! Build the preconditioner
    /*!
      @param A the base matrix for computing the preconditioner
    */
    int buildPreconditioner(operator_type& A) = 0;

    //! Reset the preconditioner
    void precReset();

    //! Return an estimation of the conditionement number of the preconditioner
    double Condest();

    //! Add A to the right of the composition
    int push_back( operator_type& A,
                   const bool useInverse=false,
                   const bool useTranspose=false );

    //! Replace the operators at position i by A
    int replace( operator_type& A,
                 const UInt index,
                 const bool useInverse=false,
                 const bool useTranspose=false );

    //@}



    //! @name Epetra Operator Interface Methods
    //@{

    int SetUseTranspose( const bool useTranspose = false );

    int Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    int ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    bool UseTranspose();

    Epetra_Map& OperatorRangeMap() const;

    Epetra_Map& OperatorDomainMap() const;

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
    void setDataFromGetPot ( const GetPot& dataFile,
                             const std::string& section ) = 0;

    //@}


    //! @name Get Methods
    //@{

    //! Preconditioner is set?
    /*!
     *  @return true
     */
    bool set() const;

    /** Get a standard pointer to the preconditioner. In most of the cases is more safe to use getPrecPtr(), which
     returns a boost::shared_ptr*/
    prec_raw_type* getPrec();

    /** get a boost::shared_ptr to the preconditioner. The only requirement on the preconditioner is that
     it must derive from the Epetra_Operator object*/
    prec_type getPrecPtr();

    //! Return the type name of the preconditioner.
    /*!
     *  @return type of the preconditioner
     */
    std::string precType();

    //! Return the number of operators in the composition
    UInt numOperators() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Short description of this method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
     */
    //void privateMethodOne();

    //@}

    prec_type                  M_prec;
    std::vector<operator_type> M_operators;
};

} // Namespace LifeV

#endif /* PRECONDITIONERCOMPOSITION_HPP */
