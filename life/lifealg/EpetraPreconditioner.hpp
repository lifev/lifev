/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

  Copyright (C) 2006 EPFL

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
   \file EpetraPreconditioner.hpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
 */


#ifndef _EPETRAPRECONDITIONER_HPP_
#define _EPETRAPRECONDITIONER_HPP_


#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>


#include <boost/shared_ptr.hpp>

#include <Ifpack_config.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/displayer.hpp>

#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV
{

// Forward declaration
class SolverTrilinos;

class EpetraPreconditioner
{
public:

    /** @name Typedefs
     */
    //@{
    typedef Epetra_Operator                      prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>     prec_type;

    typedef EpetraMatrix<double>                 operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type> operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    EpetraPreconditioner(const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>() );

    /** Copy constructor*/
    EpetraPreconditioner(EpetraPreconditioner& P, const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>() );

    //! default virtual destructor
    virtual ~EpetraPreconditioner();

    //@}


    /** @name  Methods
     */

    virtual void            setDataFromGetPot ( const GetPot& dataFile, const std::string& section ) = 0;

    virtual double          Condest() = 0;

    /** Get a standard pointer to the preconditioner. In most of the cases is more safe to use getPrecPtr(), which
     returns a boost::shared_ptr*/
    virtual prec_raw_type*  getPrec() = 0;

    /** get a boost::shared_ptr to the preconditioner. The only requirement on the preconditioner is that
     it must derive from the Epetra_Operator object*/
    virtual prec_type       getPrecPtr()=0;

    //! Return the type name of the preconditioner.
    /*!
     *  @return type of the preconditioner
     */
    virtual std::string     precType() = 0;

    virtual int             buildPreconditioner(operator_type& A) = 0;

    virtual void            precReset() = 0;

    //! returns true if prec exists
    /*const*/
    virtual bool            set() const = 0;

    virtual void            setSolver( SolverTrilinos& /*solver*/ ) {}

    // Teuchos list management
    void                    setList(Teuchos::ParameterList list);
    const Teuchos::ParameterList& getList() const;
    const int& getOverlapLevel() const;

    //! Return if the preconditioner has been created
    /*!
     *  @return true if the preconditioner has been created.
     */
    bool preconditionerCreated();

protected:

    Displayer                           M_displayer;
    int                                 M_overlapLevel;
    operator_raw_type::matrix_ptrtype   M_Oper;
    Teuchos::ParameterList              M_List;
    bool                                M_preconditionerCreated;

};

typedef singleton<factory<EpetraPreconditioner,  std::string> > PRECFactory;

} // namespace LifeV

#endif
