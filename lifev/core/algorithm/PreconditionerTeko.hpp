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
    @file
    @brief PreconditionerTeko

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 14-10-2010
 */

#ifndef PRECONDITIONERTEKO_HPP
#define PRECONDITIONERTEKO_HPP 1

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAVE_TEKO

#include <vector>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerBlock.hpp>

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

namespace LifeV
{

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
    typedef Preconditioner                          super_Type;

    typedef Epetra_Operator                         operator_Type;
    typedef std::shared_ptr<Epetra_Operator>      operatorPtr_Type;

    typedef Teko::Epetra::EpetraBlockPreconditioner preconditioner_Type;
    typedef std::shared_ptr<preconditioner_Type>  preconditionerPtr_Type;

    typedef MatrixEpetra<Real>                      matrix_Type;
    typedef std::shared_ptr<matrix_Type>          matrixPtr_Type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    PreconditionerTeko ( const std::shared_ptr<Epetra_Comm>& comm = std::shared_ptr<Epetra_Comm>() );

    /** Copy constructor*/
    PreconditionerTeko ( PreconditionerTeko& P, const std::shared_ptr<Epetra_Comm>& comm = std::shared_ptr<Epetra_Comm>() );

    //! default virtual destructor
    virtual ~PreconditionerTeko();

    //@}


    /** @name  Methods
     */

    //! Return a pointer on the preconditioner
    operator_Type* preconditioner();

    //! Return a pointer on the preconditioner
    operatorPtr_Type preconditionerPtr();

    //! Reset the preconditioner
    void resetPreconditioner();

    //! returns true if prec exists
    bool isPreconditionerSet() const;

    virtual int SetUseTranspose ( const bool useTranspose = false );
    virtual bool UseTranspose();
    virtual int ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
    virtual int Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const;
    virtual const Epetra_Map& OperatorRangeMap() const;
    virtual const Epetra_Map& OperatorDomainMap() const;

protected:

    void buildBlockGIDs ( std::vector<std::vector<int> >& gids,
                          const MapEpetra& map,
                          const std::vector<int>& blockSizes);
    void buildPreconditionerTeko ( RCP<Teko::BlockPreconditionerFactory> precFact,
                                   matrixPtr_Type& oper,
                                   const std::vector<int>& blockSizes );

    preconditionerPtr_Type M_prec;

private:
    matrix_Type::matrix_ptrtype M_oper;

};

} // namespace LifeV

#endif // HAVE_LIFEV_TEKO

#endif /* PRECONDITIONERTEKO_HPP */
