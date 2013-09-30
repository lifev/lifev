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
    @brief ConfinedOperator

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>

    @date 21-08-2012
 */

#ifndef _CONFINEDOPERATOR_HPP_
#define _CONFINEDOPERATOR_HPP_

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include <lifev/core/array/VectorBlockStructure.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/operator/LinearOperator.hpp>

namespace LifeV
{
namespace Operators
{

//! @class ConfinedOperator
/*! @brief Class which wrap an operator to apply it only on a restriction.
 *
 */
class ConfinedOperator : public LinearOperator
{
public:

    //! @name Public Typedefs and Enumerators
    //@{
    typedef Epetra_Operator                        operator_Type;
    typedef boost::shared_ptr<operator_Type>       operatorPtr_Type;
    typedef VectorBlockStructure                   blockStructure_Type;
    typedef Epetra_MultiVector                     vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;
    typedef Epetra_Comm                            comm_Type;
    typedef Epetra_Map                             map_Type;
    //@}

    //! null constructor and destructor
    //@{
#ifdef HAVE_MPI
    ConfinedOperator ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    ConfinedOperator ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif
    ~ConfinedOperator();
    //@}

    //! @name Attribute set methods
    //@{

    //! If set true, transpose of this operator will be applied.
    virtual int SetUseTranspose ( bool useTranspose );

    void setOperator ( operatorPtr_Type oper );

    void setFullMap ( const MapEpetra& map );

    void setBlockStructure ( const blockStructure_Type& blockStructure );

    void setBlockIndex ( UInt index );

    //@}

    //! @name Mathematical methods
    //@{

    //! Returns the result of a Epetra_Operator applied to a vector_Type X in Y.
    virtual int Apply ( const vector_Type& X, vector_Type& Y ) const;

    //! Returns the result of a Epetra_Operator inverse applied to an vector_Type X in Y.
    virtual int ApplyInverse ( const vector_Type& X, vector_Type& Y ) const;

    //! Returns the infinity norm of the global matrix.
    double NormInf() const;

    //@}

    //! @name Attribute access methods
    //@{

    //! Returns a character string describing the operator
    virtual const char* Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const comm_Type& Comm() const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const map_Type& OperatorDomainMap() const;

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const map_Type& OperatorRangeMap() const;

    //@}

protected:

    operatorPtr_Type               M_oper;
    blockStructure_Type            M_blockStructure;
    UInt                           M_blockIndex;
    boost::shared_ptr<Epetra_Comm> M_comm;
    boost::shared_ptr<map_Type>    M_map;

};

} /*end namespace Operators */
} /*end namespace LifeV */
#endif /* _CONFINEDOPERATOR_HPP_ */
