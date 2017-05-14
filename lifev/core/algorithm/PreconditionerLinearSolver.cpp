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
    @brief LinearSolver preconditioner

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-09-2011
 */

#include <lifev/core/algorithm/PreconditionerLinearSolver.hpp>
#include <lifev/core/util/LifeDebug.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/filter/GetPot.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerLinearSolver::PreconditionerLinearSolver ( std::shared_ptr<Epetra_Comm> comm ) :
    Preconditioner          ( comm ),
    M_printSubiterationCount ( false ),
    M_precName              ( "" ),
    M_precDataSection       ( "" )
{

}

PreconditionerLinearSolver::~PreconditionerLinearSolver()
{
    M_solver.reset();
    M_preconditioner.reset();
}

// ===================================================
// Methods
// ===================================================

void
PreconditionerLinearSolver::createParametersList ( list_Type&         list,
                                                   const GetPot&      dataFile,
                                                   const std::string& section,
                                                   const std::string& subSection )
{
    createLinearSolverList ( list, dataFile, section, subSection, M_displayer.comm()->MyPID() == 0 );
}

void
PreconditionerLinearSolver::createLinearSolverList ( list_Type&         list,
                                                     const GetPot&      dataFile,
                                                     const std::string& section,
                                                     const std::string& subsection,
                                                     const bool&        verbose )
{
    bool displayList = dataFile ( ( section + "/displayList" ).data(), false);

    // If this option is true, the solver will print the iteration count
    const std::string solverParamFile = dataFile ( ( section + "/" + subsection + "/parameters_file" ).data(), "none" );
    list = * ( Teuchos::getParametersFromXmlFile ( solverParamFile ) );

    if ( displayList && verbose )
    {
        std::cout << "PreconditionerLinearSolver parameters list:" << std::endl;
        std::cout << "-----------------------------" << std::endl;
        list.print ( std::cout );
        std::cout << "-----------------------------" << std::endl;
    }
}

Int
PreconditionerLinearSolver::buildPreconditioner ( operator_type& matrix )
{
    // Setup the solver
    M_solver.reset ( new solver_Type ( this->M_displayer.comm() ) );
    M_solver->setParameters ( M_list.sublist ( "LinearSolver" ) );
    M_solver->setOperator ( matrix );

    // Setup the preconditioner for the solver
    M_preconditioner.reset ( PRECFactory::instance().createObject ( M_precName ) );
    ASSERT ( M_preconditioner.get() != 0, " Preconditioner not set" );
    M_preconditioner->setDataFromGetPot ( M_dataFile, M_precDataSection );
    M_solver->setPreconditioner ( M_preconditioner );
    M_solver->buildPreconditioner();
    M_solver->setupSolverOperator();
    M_solver->solver()->setUsedForPreconditioning ( M_printSubiterationCount );

    this->M_preconditionerCreated = true;

    return 0;
}

void
PreconditionerLinearSolver::resetPreconditioner()
{
    M_solver.reset();
    this->M_preconditionerCreated = false;
}

Real
PreconditionerLinearSolver::condest()
{
    return 0.0;
}

void
PreconditionerLinearSolver::showMe ( std::ostream& output ) const
{
    M_solver->showMe (output);
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
PreconditionerLinearSolver::SetUseTranspose ( const bool useTranspose )
{
    return M_solver->solver()->SetUseTranspose (useTranspose);
}

bool
PreconditionerLinearSolver::UseTranspose()
{
    return M_solver->solver()->UseTranspose();
}

Int
PreconditionerLinearSolver::Apply ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_solver->solver()->Apply ( X, Y );
}

Int
PreconditionerLinearSolver::ApplyInverse ( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    if ( M_solver )
    {
        M_solver->solver()->ApplyInverse ( X, Y );
        if ( M_printSubiterationCount )
        {
            M_displayer.leaderPrint ( "> ", M_solver->numIterations(), " subiterations\n" );
        }
    }
    return 0;
}

const Epetra_Map&
PreconditionerLinearSolver::OperatorRangeMap() const
{
    return M_solver->solver()->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerLinearSolver::OperatorDomainMap() const
{
    return M_solver->solver()->OperatorDomainMap();
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerLinearSolver::setDataFromGetPot ( const GetPot& dataFile, const std::string& section )
{
    createLinearSolverList ( M_list, dataFile, section, "LinearSolver", M_displayer.comm()->MyPID() == 0 );
    M_printSubiterationCount = this->M_list.get ( "Print Subiteration Count", false );
    M_precName               = this->M_list.get ( "Preconditioner", "ML" );
    M_precDataSection        = this->M_list.get ( "Preconditioner Data Section", "" );
    M_dataFile               = dataFile;
}

void
PreconditionerLinearSolver::setSolver ( SolverAztecOO& /*solver*/ )
{

}

// ===================================================
// Get Methods
// ===================================================
//! Return true if the preconditioner is set
bool
PreconditionerLinearSolver::isPreconditionerSet() const
{
    return M_solver != nullptr ? true : false;
}

PreconditionerLinearSolver::prec_raw_type*
PreconditionerLinearSolver::preconditioner()
{
    return M_solver->solver().get();
}

PreconditionerLinearSolver::prec_type
PreconditionerLinearSolver::preconditionerPtr()
{
    return M_solver->solver();
}

PreconditionerLinearSolver::solverPtr_Type
PreconditionerLinearSolver::solverPtr()
{
    return M_solver;
}

std::string
PreconditionerLinearSolver::preconditionerType()
{
    return "LinearSolver preconditioner";
}

// ===================================================
// Private Methods
// ===================================================


} // namespace LifeV

