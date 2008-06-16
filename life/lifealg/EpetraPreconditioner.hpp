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
namespace Epetra
{

class Preconditioner
{
public:

    /** @name Typedefs
     */
    //@{
    typedef Ifpack_Preconditioner prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type> prec_type;
    typedef EpetraMatrix<double>  operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type> operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    Preconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    Preconditioner(operator_type& A);

    //@}


    /** @name  Methods
     */

    void        setDataFromGetPot ( const GetPot& dataFile, const std::string& section );

    double      Condest ();

    prec_raw_type*   getPrec();

    std::string precType(){return M_precType;}

    int         buildPreconditioner(operator_type& A);

    void        precReset();

    //! returns true if prec exists
    const bool  set() const {return M_Prec;}

private:

    void        createList( const GetPot& dataFile );

    prec_type      M_Prec;
    operator_type  M_Oper;

    Teuchos::ParameterList M_List;
    std::string            M_precType;
    int                    M_overlapLevel;

};

} // namespace Epetra
} // namespace LifeV
#endif
