/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-10-14

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
   \file PreconditionerTeko.hpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-10-14
 */


#ifndef _PRECONDITIONERTEKO_HPP_
#define _PRECONDITIONERTEKO_HPP_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifealg/Preconditioner.hpp>
#include <lifemc/lifealg/PreconditionerBlock.hpp>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_EpetraBlockPreconditioner.hpp"
#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

namespace LifeV {

//! PreconditionerTeko
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerTeko class provides a wrapper for Trilinos
 *  block preconditioners implimented with the Teko package
 */
class PreconditionerTeko:
        public PreconditionerBlock
{
public:

    /** @name Typedefs
     */
    //@{
    typedef Preconditioner                          super;

    typedef Teko::Epetra::EpetraBlockPreconditioner prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>        prec_type;

    typedef super::operator_raw_type                operator_raw_type;
    typedef super::operator_type                    operator_type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    PreconditionerTeko(const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>());

    /** Copy constructor*/
    PreconditionerTeko(PreconditionerTeko& P, const boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>() );

    //! default virtual destructor
    virtual ~PreconditionerTeko();

    //@}


    /** @name  Methods
     */

    //! Return a pointer on the preconditioner
    super::prec_raw_type* getPrec();

    //! Return a pointer on the preconditioner
    super::prec_type      getPrecPtr();

    //! Reset the preconditioner
    void                  precReset();

    //! returns true if prec exists
    /*const*/
    bool                  set() const;

    virtual int           SetUseTranspose( const bool useTranspose=false );
    virtual bool          UseTranspose();
    virtual int           ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
    virtual int           Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
    virtual const Epetra_Map & OperatorRangeMap() const;
    virtual const Epetra_Map & OperatorDomainMap() const;

protected:
    prec_type             M_prec;

    void buildBlockGIDs(std::vector<std::vector<int> > & gids,
                        const MapEpetra & map,
                        const std::vector<int>& blockSizes);
    void buildPreconditionerTeko(RCP<Teko::BlockPreconditionerFactory> precFact,
                                 operator_type& oper,
                                 const std::vector<int>& blockSizes);

private:
    operator_raw_type::matrix_ptrtype M_oper;

};

} // namespace LifeV

#endif
