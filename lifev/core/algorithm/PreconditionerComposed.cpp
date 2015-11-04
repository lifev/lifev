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

#include <lifev/core/algorithm/PreconditionerComposed.hpp>

namespace LifeV
{

// ===================================================
// Constructurs and Destructor
// ===================================================
PreconditionerComposed::PreconditionerComposed ( std::shared_ptr<Epetra_Comm> comm) :
    super_Type (comm ),
    M_operVector (0),
    M_prec (new prec_Type (comm) )
    //M_precType()
{
}

/*
PreconditionerComposed::PreconditionerComposed(PreconditionerComposed& P):
    super_Type(P, std::dynamic_pointer_cast<ComposedOperator<Ifpack_Preconditioner> >(P.preconditionerPtr())->commPtr()),
    M_operVector(P.operVector()),
    M_prec(new prec_Type(*std::dynamic_pointer_cast<prec_Type>(P.preconditionerPtr())))
    //M_precType(P.preconditionerType())
{
    //    *M_prec=*P.preconditioner();
}
*/

PreconditionerComposed::~PreconditionerComposed()
{
    M_prec.reset();
    M_operVector.clear();
}


// ===================================================
// Public Methods
// ===================================================
void
PreconditionerComposed::setDataFromGetPot ( const GetPot&      dataFile,
                                            const std::string& section )
{
    myCreateParametersList (dataFile, section, "Composed" );
}

void
PreconditionerComposed::createParametersList (list_Type& /*list*/,
                                              const GetPot&      dataFile,
                                              const std::string& section,
                                              const std::string& subSection )
{
    myCreateParametersList (dataFile,section,subSection);
}

void PreconditionerComposed::myCreateParametersList (const GetPot&      dataFile,
                                                     const std::string& section,
                                                     const std::string& subSection)
{
    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters
    ASSERT ( !M_prec->number(), "Error, when initializing the preconditioner, it must be empty" );
    for ( UInt i (0); i < dataFile.vector_variable_size ( ( section + "/" + subSection + "/list" ).data() ); ++i )
    {
        epetraPrecPtr_Type tmp ( PRECFactory::instance().createObject ( dataFile ( ( section + "/" + subSection + "/list" ).data(), "ML", i ) ) );
        M_prec->push_back (tmp);
        M_prec->OperatorView() [i]->createParametersList (M_prec->OperatorView() [i]->parametersList(), dataFile, section, dataFile ( ( section + "/" + subSection + "/sections" ).data(), "ML", i ) );
    }
}

double
PreconditionerComposed::condest()
{
    return M_prec->Condest();
}

Preconditioner::prec_raw_type*
PreconditionerComposed::preconditioner()
{
    return M_prec.get();
}

int
PreconditionerComposed::buildPreconditioner (operatorPtr_Type& oper)
{
    //M_prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back (oper);
}

int
PreconditionerComposed::buildPreconditioner (operatorPtr_Type& oper,
                                             const bool useInverse,
                                             const bool useTranspose)
{
    //M_prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back (oper, useInverse, useTranspose);
}

int
PreconditionerComposed::push_back (operatorPtr_Type& oper,
                                   const bool useInverse,
                                   const bool useTranspose
                                  )
{
    if (!M_prec.get() )
    {
        M_prec.reset (new prec_Type (M_displayer.comm() ) );
    }
    M_operVector.push_back (oper);
    LifeChrono chrono;
    epetraPrecPtr_Type prec;

    this->M_displayer.leaderPrint ( std::string ("ICP-  Preconditioner type:                     ") + M_prec->Operator() [M_operVector.size() - 1]->preconditionerType() + std::string ("\n") );
    this->M_displayer.leaderPrint ( "ICP-  Computing preconditioner ...             " );
    chrono.start();
    createPrec (oper, M_prec->OperatorView() [M_operVector.size() - 1]);
    chrono.stop();
    this->M_displayer.leaderPrintMax ("done in ", chrono.diff() );
    M_prec->replace (prec, useInverse, useTranspose); // \TODO to reset as push_back
    if ( M_prec->Operator().size() == M_operVector.size() )
    {
        this->M_preconditionerCreated = true;
    }
    return EXIT_SUCCESS;
}

int
PreconditionerComposed::replace (operatorPtr_Type& oper,
                                 const UInt index,
                                 const bool useInverse,
                                 const bool useTranspose)
{
    ASSERT (index <= M_operVector.size(), "PreconditionerComposed::replace: index too large");

    M_operVector[index] = oper;
    LifeChrono chrono;
    //ifpack_prec_type prec;
    this->M_displayer.leaderPrint ( std::string ("ICP-  Preconditioner type:                     ") + M_prec->Operator() [index]->preconditionerType() + std::string ("\n") );
    this->M_displayer.leaderPrint ( "ICP-  Computing preconditioner ...             " );
    chrono.start();
    createPrec (oper, M_prec->OperatorView() [index]);
    chrono.stop();
    this->M_displayer.leaderPrintMax ("done in ", chrono.diff() );

    M_prec->replace (M_prec->Operator() [index], index, useInverse, useTranspose);

    return EXIT_SUCCESS;
}

void
PreconditionerComposed::resetPreconditioner()
{
    //M_operVector.reset();
    M_prec.reset (new prec_Type (M_displayer.comm() ) );

    this->M_preconditionerCreated = false;
}


// ===================================================
// Private Methods
// ===================================================
Int
PreconditionerComposed::createPrec (operator_type& oper,
                                    std::shared_ptr<Preconditioner>& prec )
{
    return prec->buildPreconditioner ( oper );
}



bool PreconditionerComposed::registerComposed = PRECFactory::instance().registerProduct ( "Composed", &PreconditionerComposed::createComposedPreconditioner );

} // namespace LifeV
