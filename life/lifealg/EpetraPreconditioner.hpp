//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
 *  @file
 *  @brief Epetra preconditioner
 *
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @date 09-11-2006
 */

#ifndef _EPETRAPRECONDITIONER_HPP_
#define _EPETRAPRECONDITIONER_HPP_

#include <boost/shared_ptr.hpp>

#include <Teuchos_ParameterList.hpp>

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/displayer.hpp>

#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV {

// Forward declaration
class SolverTrilinos;

class EpetraPreconditioner
{
public:

    //! @name Typedefs
    //@{

    typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef EpetraMatrix<double>                 operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type> operator_type;

    typedef Displayer::comm_Type                 comm_Type;
    typedef Displayer::comm_PtrType              comm_PtrType;

    typedef Teuchos::ParameterList               list_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor.
    EpetraPreconditioner( const comm_PtrType& comm = comm_PtrType() );

    /** Copy constructor*/
    EpetraPreconditioner( const EpetraPreconditioner& P, const comm_PtrType& comm = comm_PtrType() );

    //! Default virtual destructor
    virtual ~EpetraPreconditioner();

    //@}


    //! @name Methods
    //@{

    virtual void createList( const GetPot& dataFile, const std::string& section, list_Type& list);

    //! Build the preconditioner
    /*!
     *  @param A the base matrix for computing the preconditioner
     */
    virtual int buildPreconditioner( operator_type& A ) = 0;

    virtual void precReset() = 0;

    //! Compute the condition number of the preconditioner
    /*!
     *  @return Condition number of the preconditioner
     */
    virtual double Condest() = 0;

    //@}


    //! @name Epetra Operator Interface Methods
    //@{

    virtual int SetUseTranspose( const bool useTranspose = false );

    virtual int Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    virtual int ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;

    virtual bool UseTranspose();

    virtual const Epetra_Map& OperatorRangeMap() const;

    virtual const Epetra_Map& OperatorDomainMap() const;

    //@}


    //! @name Set Methods
    //@{

    void setList( const list_Type& list );

    virtual void setDataFromGetPot ( const GetPot& dataFile, const std::string& section, const UInt& listNumber=1 ) = 0;

    virtual void setSolver( SolverTrilinos& /*solver*/ );

    //@}


    //! @name Get Methods
    //@{

    //! Return if the preconditioner has been created
    /*!
     *  @return true if the preconditioner has been created.
     */
    const bool& preconditionerCreated();

    //! Preconditioner is set?
    /*!
     *  @return true
     */
    virtual bool set() const = 0;

    /** Get a standard pointer to the preconditioner. In most of the cases is more safe to use getPrecPtr(), which
     returns a boost::shared_ptr*/
    virtual prec_raw_type* getPrec() = 0;

    /** get a boost::shared_ptr to the preconditioner. The only requirement on the preconditioner is that
     it must derive from the Epetra_Operator object*/
    virtual prec_type getPrecPtr() = 0;

    //! Return the type name of the preconditioner.
    /*!
     *  @return type of the preconditioner
     */
    virtual std::string precType() = 0;

    const list_Type& getList() const;

    list_Type& list();

    //@}

protected:

    std::string                         M_precType;
    Displayer                           M_displayer;
    list_Type                           M_List;
    bool                                M_preconditionerCreated;

private:

    static void createPreconditionerList( const GetPot& dataFile, const std::string& section, list_Type& list) {};

};

typedef singleton<factory<EpetraPreconditioner, std::string> > PRECFactory;

} // namespace LifeV

#endif
