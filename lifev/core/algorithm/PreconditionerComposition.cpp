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

#include <lifev/core/algorithm/PreconditionerComposition.hpp>
#include <lifev/core/operator/ConfinedOperator.hpp>

// <--Check for necessity
#include <Ifpack_config.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

PreconditionerComposition::PreconditionerComposition ( std::shared_ptr<Epetra_Comm> comm ) :
    super_Type ( comm ),
    M_comm ( comm ),
    M_prec ( new prec_Type ( comm ) )
{

}

PreconditionerComposition::PreconditionerComposition ( const PreconditionerComposition& precComp ) :
    super_Type ( precComp, precComp.M_comm ),
    M_comm ( precComp.M_comm ),
    M_prec ( new prec_Type ( * ( precComp.M_prec.get() ) ) )
{

}

PreconditionerComposition::~PreconditionerComposition()
{
    M_prec.reset();
    M_precBaseOperators.clear();
}

// ===================================================
// Methods
// ===================================================

void
PreconditionerComposition::resetPreconditioner()
{
    M_prec.reset ( new prec_Type ( M_comm ) );
    M_precBaseOperators.clear();
    this->M_preconditionerCreated = false;
}

Real
PreconditionerComposition::condest()
{
    // todo Is there a way to obtain condest?
    //return M_prec->Condest();
    return 0.0;
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
int
PreconditionerComposition::SetUseTranspose ( const bool useTranspose )
{
    return M_prec->SetUseTranspose ( useTranspose );
}

int
PreconditionerComposition::Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->Apply ( X, Y );
}

int
PreconditionerComposition::ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->ApplyInverse ( X, Y );
}

bool
PreconditionerComposition::UseTranspose()
{
    return M_prec->UseTranspose();
}

const Epetra_Map&
PreconditionerComposition::OperatorRangeMap() const
{
    return M_prec->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerComposition::OperatorDomainMap() const
{
    return M_prec->OperatorDomainMap();
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerComposition::setComm ( std::shared_ptr<Epetra_Comm> comm )
{
    M_comm = comm;
    M_prec->setComm ( comm );
}

// ===================================================
// Get Methods
// ===================================================
bool
PreconditionerComposition::isPreconditionerSet() const
{
    return M_prec != nullptr ? true : false;
}

PreconditionerComposition::operator_Type*
PreconditionerComposition::preconditioner()
{
    return M_prec.get();
}

PreconditionerComposition::operatorPtr_Type
PreconditionerComposition::preconditionerPtr()
{
    return M_prec;
}

std::string
PreconditionerComposition::preconditionerType()
{
    return M_precType;
}

UInt
PreconditionerComposition::numOperators() const
{
    return M_prec->number();
}

// ===================================================
// Protected Methods
// ===================================================
int
PreconditionerComposition::pushBack ( matrixPtr_Type A,
                                      const bool useInverse,
                                      const bool useTranspose )
{
    //std::cout << "[DEBUG] pushBack() matrix version" << std::endl;
    M_prec->push_back ( std::dynamic_pointer_cast<operator_Type> ( A->matrixPtr() ), useInverse, useTranspose );

    return EXIT_SUCCESS;
}

int
PreconditionerComposition::pushBack ( operatorPtr_Type oper,
                                      const bool useInverse,
                                      const bool useTranspose,
                                      matrixPtr_Type baseMatrix )
{
    if ( baseMatrix.get() != 0 )
    {
        M_precBaseOperators.push_back ( baseMatrix );
    }
    M_prec->push_back ( oper, useInverse, useTranspose );

    return EXIT_SUCCESS;
}

int
PreconditionerComposition::pushBack ( matrixPtr_Type A,
                                      superPtr_Type preconditioner,
                                      const bool useInverse,
                                      const bool useTranspose )
{
    //std::cout << "[DEBUG] pushBack() preconditioner version" << std::endl;
    M_precBaseOperators.push_back ( A );
    preconditioner->buildPreconditioner ( A );
    operatorPtr_Type oper ( preconditioner->preconditionerPtr() );
    M_prec->push_back ( oper, useInverse, useTranspose );

    return EXIT_SUCCESS;
}

int
PreconditionerComposition::pushBack ( matrixPtr_Type embeddedA,
                                      superPtr_Type preconditioner,
                                      const VectorBlockStructure& blockStructure,
                                      const UInt& blockIndex,
                                      const MapEpetra& fullMap,
                                      const bool useInverse,
                                      const bool useTranspose,
                                      const bool buildPreconditioner )
{
    // Add the operator
    M_precBaseOperators.push_back ( embeddedA );

    // Build the preconditioner
    if ( buildPreconditioner )
    {
        preconditioner->buildPreconditioner ( embeddedA );
    }
    operatorPtr_Type precOper ( preconditioner->preconditionerPtr() );

    // Wrap the preconditioner in a ConfinedOperator
    Operators::ConfinedOperator* confinedOperator = new Operators::ConfinedOperator ( M_comm );
    confinedOperator->setOperator ( precOper );
    confinedOperator->setFullMap ( fullMap );
    confinedOperator->setBlockStructure ( blockStructure );
    confinedOperator->setBlockIndex ( blockIndex );
    operatorPtr_Type oper ( confinedOperator );

    // Add the operator
    M_prec->push_back ( oper, useInverse, useTranspose );

    return EXIT_SUCCESS;
}

int
PreconditionerComposition::pushBack ( operatorPtr_Type embeddedOperator,
                                      const VectorBlockStructure& blockStructure,
                                      const UInt& blockIndex,
                                      const MapEpetra& fullMap,
                                      const bool useInverse,
                                      const bool useTranspose )
{
    // Wrap the preconditioner in a ConfinedOperator
    Operators::ConfinedOperator* confinedOperator = new Operators::ConfinedOperator ( M_comm );
    confinedOperator->setOperator ( embeddedOperator );
    confinedOperator->setFullMap ( fullMap );
    confinedOperator->setBlockStructure ( blockStructure );
    confinedOperator->setBlockIndex ( blockIndex );
    operatorPtr_Type oper ( confinedOperator );

    // Add the operator
    M_prec->push_back ( oper, useInverse, useTranspose );

    return EXIT_SUCCESS;
}

} // Namespace LifeV
