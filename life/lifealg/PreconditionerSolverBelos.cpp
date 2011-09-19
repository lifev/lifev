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
    @brief Solver Belos Operator

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-09-2011
 */

#include <life/lifealg/PreconditionerSolverBelos.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifefilters/GetPot.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerSolverBelos::PreconditionerSolverBelos( boost::shared_ptr<Epetra_Comm> comm ) :
        Preconditioner( comm )
{

}

PreconditionerSolverBelos::~PreconditionerSolverBelos()
{
    M_prec.reset();
}

// ===================================================
// Methods
// ===================================================

void
PreconditionerSolverBelos::createParametersList( list_Type&         list,
                                                  const GetPot&      dataFile,
                                                  const std::string& section,
                                                  const std::string& subSection )
{
    createSolverBelosList( list, dataFile, section, subSection );
}

void
PreconditionerSolverBelos::createSolverBelosList( list_Type&         list,
                                                  const GetPot&      dataFile,
                                                  const std::string& section,
                                                  const std::string& subsection )
{
    bool displayList = dataFile( ( section + "/displayList" ).data(), false);

    bool found;

    // Flexible Gmres will be used to solve this problem
    bool flexibleGmres = dataFile( ( section + "/" + subsection + "/flexible_gmres" ).data(), false, found );
    if ( found ) list.set( "Flexible Gmres", flexibleGmres );

    // Relative convergence tolerance requested
    Real tolerance = dataFile( ( section + "/" + subsection + "/tol" ).data(), 1.e-6, found );
    if ( found ) list.set( "Convergence Tolerance", tolerance );

    // Maximum number of iterations allowed
    Int maxIter = dataFile( ( section + "/" + subsection + "/max_iter"      ).data(), 200, found );
    if ( found ) list.set( "Maximum Iterations", maxIter );

    // Output Frequency
    Int outputFrequency = dataFile( ( section + "/" + subsection + "/output_frequency" ).data(), 1, found );
    if ( found ) list.set( "Output Frequency", outputFrequency );

    // Blocksize to be used by iterative solver
    Int blockSize = dataFile( ( section + "/" + subsection + "/block_size" ).data(), 10, found );
    if ( found ) list.set( "Block Size", blockSize );

    // Maximum number of blocks in Krylov factorization
    Int numBlocks = dataFile( ( section + "/" + subsection + "/num_blocks" ).data(), 10, found );
    if ( found ) list.set( "Num Blocks", numBlocks );

    // Maximum number of restarts allowed
    Int maximumRestarts = dataFile( ( section + "/" + subsection + "/maximum_restarts" ).data(), 0, found );
    if ( found ) list.set( "Maximum Restarts", maximumRestarts );

    // Set the output style (General Brief)
    std::string outputStyle = dataFile( ( section + "/" + subsection + "/output_style" ).data(), "brief", found );
    if ( found )
    {
        if ( outputStyle == "brief" )   list.set( "Output Style", Belos::Brief );
        if ( outputStyle == "general" ) list.set( "Output Style", Belos::General );
    }

    // Setting the desired output informations
    bool msgEnable( false );
    int msg = Belos::Errors;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_warnings" ).data(), true, found );
    if ( found && msgEnable ) msg += Belos::Warnings;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_iterations_details" ).data(), true, found );
    if ( found && msgEnable ) msg += Belos::IterationDetails;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_ortho_details" ).data(), false, found );
    if ( found && msgEnable ) msg += Belos::OrthoDetails;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_final_summary" ).data(), false, found );
    if ( found && msgEnable ) msg += Belos::FinalSummary;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_timing_details" ).data(), false, found );
    if ( found && msgEnable ) msg += Belos::TimingDetails;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_status_test_details" ).data(), false, found );
    if ( found && msgEnable ) msg += Belos::StatusTestDetails;
    msgEnable = dataFile( ( section + "/" + subsection + "/enable_debug" ).data(), false, found );
    if ( found && msgEnable ) msg += Belos::Debug;
    list.set( "Verbosity", msg );

    // LifeV features

    // Reuse the preconditioner from one to another call
    bool reusePreconditioner = dataFile( ( section + "/" + subsection + "/reuse_preconditioner" ).data(), true, found );
    if ( found ) list.set( "Reuse preconditioner", reusePreconditioner );

    // Max iterations allowed to reuse the preconditioner
    Int maxItersForReuse = dataFile( ( section + "/" + subsection + "/max_iters_for_reuse" ).data(), static_cast<Int> ( maxIter*8./10. ), found );
    if ( found ) list.set( "max iteration for reuse", maxItersForReuse );

    // If quitOnFailure is enabled and if some problems occur
    // the simulation is stopped
    bool quitOnFailure = dataFile( ( section + "/" + subsection + "/quit_on_failure").data(), false, found );
    if ( found ) list.set( "Quit on failure", quitOnFailure );

    // All the information different from warnings and errors are
    // not displayed
    bool silent = dataFile( ( section + "/" + subsection + "/silent").data(), false, found );
    if ( found ) list.set( "Silent", silent );

    std::string prec = dataFile( ( section + "/" + subsection + "/prec" ).data(), "ML" );
    list.set( "prec", prec );
    std::string precDataSection = dataFile( ( section + "/" + subsection + "/prec_data_section" ).data(), "" );
    list.set( "prec data section", ( section + "/" + subsection+"/"+precDataSection ).data() );

    if ( displayList ) list.print( std::cout );
}

Int
PreconditionerSolverBelos::buildPreconditioner( operator_type& matrix )
{
    M_prec.reset( new precOperator_Type( this->M_displayer.comm() ) );
    M_prec->buildSolver( matrix, M_list );
    M_prec->buildPreconditioner( M_dataFile, M_precDataSection );

    this->M_preconditionerCreated = true;

    return 0;
}

void
PreconditionerSolverBelos::resetPreconditioner()
{
    M_prec.reset();
    this->M_preconditionerCreated = false;
}

Real
PreconditionerSolverBelos::condest()
{
    return 0.0;
}

void
PreconditionerSolverBelos::showMe( std::ostream& output ) const
{
    M_prec->showMe(output);
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
PreconditionerSolverBelos::SetUseTranspose( const bool useTranspose )
{
    return M_prec->SetUseTranspose(useTranspose);
}

bool
PreconditionerSolverBelos::UseTranspose()
{
    return M_prec->UseTranspose();
}

Int
PreconditionerSolverBelos::Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->Apply( X, Y );
}

Int
PreconditionerSolverBelos::ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    if( M_prec )
    {
        M_prec->ApplyInverse( X, Y );
    }
    return 0;
}

const Epetra_Map&
PreconditionerSolverBelos::OperatorRangeMap() const
{
    return M_prec->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerSolverBelos::OperatorDomainMap() const
{
    return M_prec->OperatorDomainMap();
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerSolverBelos::setDataFromGetPot ( const GetPot& dataFile, const std::string& section )
{
    createSolverBelosList( M_list, dataFile, section, "SolverBelos" );
    M_solverPrecName    = this->M_list.get( "prec", "ML" );
    M_precDataSection   = this->M_list.get( "prec data section", "" );
    M_dataFile = dataFile;
}

void
PreconditionerSolverBelos::setSolver( SolverAztecOO& /*solver*/ )
{

}

// ===================================================
// Get Methods
// ===================================================
//! Return true if the preconditioner is set
bool
PreconditionerSolverBelos::isPreconditionerSet() const
{
    return M_prec;
}

PreconditionerSolverBelos::prec_raw_type*
PreconditionerSolverBelos::preconditioner()
{
    return M_prec.get();
}

PreconditionerSolverBelos::prec_type
PreconditionerSolverBelos::preconditionerPtr()
{
    return M_prec;
}

std::string
PreconditionerSolverBelos::preconditionerType()
{
    return "SolverBelos preconditioner";
}

// ===================================================
// Private Methods
// ===================================================


} // namespace LifeV

