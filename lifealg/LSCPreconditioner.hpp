/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-08

  Copyright (C) 2010 EPFL

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
   \file LSCPreconditioner.hpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-08
 */


#ifndef _LSCPRECONDITIONER_HPP_
#define _LSCPRECONDITIONER_HPP_

#include <boost/shared_ptr.hpp>

#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <lifemc/lifealg/TekoPreconditioner.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

namespace LifeV {

//! LSCPreconditioner
/*!
 *  @author Gwenol Grandperrin
 *
 *  The LSCPreconditioner class provides the LSC block preconditioner
 *  available in the Teko package of Trilinos
 */
class LSCPreconditioner:
        public TekoPreconditioner
{
public:

    /** @name Public Types
     */
    //@{
    typedef EpetraPreconditioner                    super;

    typedef Teko::Epetra::EpetraBlockPreconditioner prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>        prec_type;

    typedef super::operator_raw_type                operator_raw_type;
    typedef super::operator_type                    operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    LSCPreconditioner();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner(operator_type& A);

    //! default destructor
    ~LSCPreconditioner();

    //@}

    /** @name  Methods
     */

    //! Setter using GetPot
    /*!
        This method use GetPot to load data from a file and then set
        the preconditioner.
        @param dataFile is a GetPot dataFile
        @param section is the section containing the data
     */
    void        setDataFromGetPot ( const GetPot&      dataFile,
                                               const std::string& section );

    //! Return an estimation of the conditionement number of the preconditioner
    double      Condest ();

    //! Return the name of the preconditioner to be used in the factory
    std::string precType(){return M_precType;}

    //! Build the preconditioner
    int         buildPreconditioner(operator_type& A);

    int         numBlockRow() const;
    int         numBlockCol() const;
protected:

    std::string M_precType;
    int         M_velocityBlockSize;
    int         M_pressureBlockSize;

};

void createLSCList( const GetPot&              dataFile,
                    const std::string&         section,
                    Teuchos::ParameterList&    list);

inline EpetraPreconditioner* createLSC(){ return new LSCPreconditioner(); }
namespace
{
	static bool registerLSC = PRECFactory::instance().registerProduct( "LSC", &createLSC );
}

} // namespace LifeV

#endif
