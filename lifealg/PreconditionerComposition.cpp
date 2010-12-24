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

#include <PreconditionerComposition.hpp>

// <--Check for necessity
#include <Ifpack_config.h>
#include <Ifpack.h>
#include <Ifpack_Preconditioner.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Ifpack_AdditiveSchwarz.h>
#include <Ifpack_Amesos.h>
#include <Ifpack_ILU.h>
// -->

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================

PreconditionerComposition::PreconditionerComposition( boost::shared_ptr<Epetra_Comm> comm ):
    super (comm ),
    M_prec(new prec_raw_type(comm))
{

}

PreconditionerComposition::PreconditionerComposition( const PreconditionerComposition& precComp )
{

}

PreconditionerComposition::~PreconditionerComposition()
{

}

// ===================================================
// Methods
// ===================================================
/*
void PreconditionerComposition::createList( list_type& list,
                     const GetPot& dataFile,
                     const std::string& section,
                     const std::string& subSection )
{

}

int PreconditionerComposition::buildPreconditioner(operator_type& A)
{

}
*/

void PreconditionerComposition::precReset()
{
    M_prec.reset();
    this->M_preconditionerCreated = false;
}

double PreconditionerComposition::Condest()
{
    return M_prec->Condest();
}

int PreconditionerComposition::push_back( operator_type& A,
                                          const bool useInverse,
                                          const bool useTranspose )
{
    M_operators.push_back(A);

    return EXIT_SUCCESS;
}

int PreconditionerComposition::replace( operator_type& A,
                                        const UInt index,
                                        const bool useInverse,
                                        const bool useTranspose )
{
    ASSERT(index <= M_operators.size(), "ComposedPreconditioner::replace: index too large");
    M_operators[index] = A;
    //M_prec->replace(index,useInverse,useTranspose);

    return EXIT_SUCCESS;
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
int PreconditionerComposition::SetUseTranspose( const bool useTranspose )
{
    return M_prec->SetUseTranspose(useTranspose);
}

int PreconditionerComposition::Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->Apply(X,Y);
}

int PreconditionerComposition::ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->ApplyInverse(X,Y);
}

bool PreconditionerComposition::UseTranspose()
{
    return M_prec->UseTranspose();
}

const Epetra_Map& PreconditionerComposition::OperatorRangeMap() const
{
    return M_prec->OperatorRangeMap();
}

const Epetra_Map& PreconditionerComposition::OperatorDomainMap() const
{
    return M_prec->OperatorDomainMap();
}

// ===================================================
// Set Methods
// ===================================================
/*
void PreconditionerComposition::setDataFromGetPot ( const GetPot& dataFile,
                                                    const std::string& section )
{

}
*/

// ===================================================
// Get Methods
// ===================================================
bool PreconditionerComposition::set() const
{
    return M_prec;
}

PreconditionerComposition::prec_raw_type* PreconditionerComposition::getPrec()
{
    return M_prec.get();
}

PreconditionerComposition::operator_type PreconditionerComposition::preconditionerPtr()
{
    return M_prec;
}

std::string PreconditionerComposition::preconditionerType()
{
    return "PreconditionerComposition";
}

UInt PreconditionerComposition::numOperators() const
{
    return M_prec->number();
}

// ===================================================
// Private Methods
// ===================================================


} // Namespace LifeV
