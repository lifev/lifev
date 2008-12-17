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
#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV
{
// namespace Epetra
// {

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
    EpetraPreconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
//     EpetraPreconditioner(operator_type& A);

    //! default virtual destructor

    virtual ~EpetraPreconditioner();

    //@}


    /** @name  Methods
     */

    virtual void           setDataFromGetPot ( const GetPot& dataFile, const std::string& section ) = 0;

    virtual double         Condest() = 0;

    virtual prec_raw_type* getPrec() = 0;

    //std::string            precType() { return M_precType; }

    virtual int            buildPreconditioner(operator_type& A) = 0;

    virtual void           precReset() = 0;

    //! returns true if prec exists
    /*const*/
    virtual bool           set() const = 0;

protected:

    int                    M_overlapLevel;

    operator_type          M_Oper;

    Teuchos::ParameterList M_List;
    //    std::string            M_precType;


private:



};



// } // namespace Epetra

typedef boost::shared_ptr<EpetraPreconditioner>                 prec_ptr;
typedef singleton<factory<EpetraPreconditioner,  std::string> > PRECFactory;


} // namespace LifeV
#endif
