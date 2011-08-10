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
    @brief SolverBelos

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 03-08-2011
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifealg/SolverBelos.hpp>
#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
SolverBelos::SolverBelos() :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManagerType    (BlockGmres),
        M_solverManager        (),
        M_problem              ( new LinearProblem_type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer() ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false),
        M_quitOnFailure        (false),
        M_lossOfPrecision      (false),
        M_maxNumItersReached   (false)
{
    M_problem->setLabel("SolverBelos");
}

SolverBelos::SolverBelos( const boost::shared_ptr<Epetra_Comm>& comm ) :
        M_leftPreconditioner   (),
        M_rightPreconditioner  (),
        M_solverManagerType    (BlockGmres),
        M_solverManager        (),
        M_problem              ( new LinearProblem_type() ),
        M_parameterList        (),
        M_displayer            ( new Displayer(comm) ),
        M_maxIterForReuse      ( 0 ),
        M_reusePreconditioner  (false),
        M_quitOnFailure        (false),
        M_lossOfPrecision      (false),
        M_maxNumItersReached   (false)
{
    M_problem->setLabel("SolverBelos");
}

SolverBelos::~SolverBelos()
{

}

// ===================================================
// Methods
// ===================================================
Int
SolverBelos::solve( vector_type& solution )
{
    // Reset status informations
    M_lossOfPrecision = false;
    M_maxNumItersReached = false;
    bool failure = false;

    // Setting the unknown in the system
    Teuchos::RCP<vector_type::vector_type> solutionPtr(&(solution.epetraVector()),false);
    M_problem->setLHS( solutionPtr );

    // Build preconditioners if needed
    bool retry( true );
    if ( !isPreconditionerSet() || !M_reusePreconditioner  )
    {
        if( M_baseMatrixForPreconditioner.get()==0 )
        {
            M_displayer->leaderPrint( "SLV-  Build preconditioner using the problem matrix\n");
            buildPreconditioner( M_matrix );
        }
        else
        {
            M_displayer->leaderPrint( "SLV-  Build preconditioner using the base matrix provided\n");
            buildPreconditioner( M_baseMatrixForPreconditioner );
        }
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_displayer->leaderPrint( "SLV-  Reusing precond ...                 \n" );
    }

    bool set = M_problem->setProblem();
    if (set == false) {
        M_displayer->leaderPrint( "SLV-  ERROR: SolverBelos failed to set up correctly!\n");
        return -1;
    }

    // Setup the Solver Manager
    setupSolverManager();

    // Solve the linear system
    Belos::ReturnType ret = M_solverManager->solve();

    // Getting informations post-solve
    Int numIters = M_solverManager->getNumIters();
    bool M_lossOfPrecision = M_solverManager->isLOADetected();

    // todo missing the automatic second run

    if(M_lossOfPrecision)
    {
        M_displayer->leaderPrint("SLV-  WARNING: Loss of accuracy detected!\n");
        failure = true;
    }
    if(ret == Belos::Converged)
    {
        M_displayer->leaderPrint( "SLV-  Convergence in " , numIters, " iterations\n" );
    }
    else
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos failed to converged to the desired precision!\n");
        M_maxNumItersReached = true;
        failure = true;
    }

    // If quitOnFailure is enabled and if some problems occur
    // the simulation is stopped
    if(M_quitOnFailure && failure)
        exit(-1);

    return numIters;
}

Real
SolverBelos::computeResidual( vector_type& solution )
{
    if( M_problem->getOperator()==null || M_problem->getRHS()==null )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos can not compute the residual if the operator or the rhs is not set!\n");
        return -1;
    }
    vector_type res( solution.map() );
    M_problem->computeCurrResVec(&res.epetraVector(),&solution.epetraVector(),M_problem->getRHS().get());
    Real residual;
    res.norm2( &residual );
    return residual;
}

std::string
SolverBelos::printStatus()
{
    std::ostringstream stat;
    std::string str;

    /*
    // AztecOO informations
    if ( status[AZ_why] == AZ_normal         ) stat << "Normal Convergence    ";
    else if ( status[AZ_why] == AZ_maxits    ) stat << "Maximum iters reached ";
    else if ( status[AZ_why] == AZ_loss      ) stat << "Accuracy loss         ";
    else if ( status[AZ_why] == AZ_ill_cond  ) stat << "Ill-conditioned       ";
    else if ( status[AZ_why] == AZ_breakdown ) stat << "Breakdown             ";

    stat << setw(12) << "res = " << status[AZ_scaled_r];
    stat << setw(4)  << " " << (Int)status[AZ_its] << " iters. ";
    stat << std::endl;
    */

    if (M_lossOfPrecision) stat    << "Accuracy loss ";
    if (M_maxNumItersReached) stat << "Maximum number of iterations reached ";
    str = stat.str();
    return str;
}

Int
SolverBelos::solveSystem( const vector_type& rhsFull,
                              vector_type&       solution,
                              matrix_ptrtype&    baseMatrixForPreconditioner )

{
    // todo redo the implementation or delete the method
    // todo deal with preconditioner
    bool retry( true );

    LifeChrono chrono;

    M_displayer->leaderPrint( "SLV-  Setting up the solver ...                \n" );

    if ( baseMatrixForPreconditioner.get() == 0 )
        M_displayer->leaderPrint( "SLV-  WARNING: BaseMatrixForPreconditioner is empty     \n" );

    if ( !isPreconditionerSet() || !M_reusePreconditioner  )
    {
        buildPreconditioner( baseMatrixForPreconditioner );
        // do not retry if I am recomputing the preconditioner
        retry = false;
    }
    else
    {
        M_displayer->leaderPrint( "SLV-  Reusing precond ...                 \n" );
    }

    Int numIter = solveSystem( rhsFull, solution, M_rightPreconditioner );

    // If we do not want to retry, return now.
    // otherwise rebuild the preconditioner and solve again:
    if ( numIter < 0  && retry )
    {
        chrono.start();

        M_displayer->leaderPrint( "SLV-  Iterative solver failed, numiter = " , - numIter, "\n" );
        //M_displayer->leaderPrint( "SLV-  maxIterSolver = " , M_maxIter );
        M_displayer->leaderPrint( "SLV-  retrying:\n" );

        buildPreconditioner( baseMatrixForPreconditioner );

        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff() );
        // Solving again, but only once (retry = false)
        numIter = solveSystem( rhsFull, solution, M_rightPreconditioner );

        if ( numIter < 0 )
            M_displayer->leaderPrint( " ERROR: Iterative solver failed again.\n" );
    }

    if ( std::abs(numIter) > M_maxIterForReuse )
        resetPreconditioner();

    return numIter;
}

void
SolverBelos::setPreconditionerFromGetPot( const GetPot& dataFile,  const std::string& section )
{
    std::string precType = dataFile( (section + "/prectype").data(), "Ifpack" );
    M_rightPreconditioner.reset( PRECFactory::instance().createObject( precType ) );

    ASSERT( M_rightPreconditioner.get() != 0, " Preconditioner not set" );

    M_rightPreconditioner->setDataFromGetPot( dataFile, section );
}

void
SolverBelos::buildPreconditioner( matrix_ptrtype& preconditioner )
{
    LifeChrono chrono;
    Real condest(-1);

    if(M_leftPreconditioner)
    {
        chrono.start();
        M_displayer->leaderPrint( "SLV-  Computing the left preconditioner...    " );
        M_leftPreconditioner->buildPreconditioner( preconditioner );
        condest = M_leftPreconditioner->condest();
        Teuchos::RCP<OP> leftPrec(M_leftPreconditioner->preconditioner(),false);
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( leftPrec ) );
        M_problem->setLeftPrec(belosPrec);
        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff());
        M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
    }
    if(M_rightPreconditioner)
    {
        chrono.start();
        M_displayer->leaderPrint( "SLV-  Computing the right preconditioner...    " );
        M_rightPreconditioner->buildPreconditioner( preconditioner );
        condest = M_rightPreconditioner->condest();
        Teuchos::RCP<OP> rightPrec(M_rightPreconditioner->preconditioner(),false);
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( rightPrec ) );
        M_problem->setRightPrec(belosPrec);
        chrono.stop();
        M_displayer->leaderPrintMax( "done in " , chrono.diff());
        M_displayer->leaderPrint( "SLV-  Estimated condition number               " , condest, "\n" );
    }
}

void
SolverBelos::resetPreconditioner()
{
    M_leftPreconditioner->resetPreconditioner();
    M_rightPreconditioner->resetPreconditioner();
    M_problem->setLeftPrec(Teuchos::RCP<OP>(null));
    M_problem->setRightPrec(Teuchos::RCP<OP>(null));
}

bool
SolverBelos::isPreconditionerSet() const
{
    return ( M_leftPreconditioner.get() !=0 && M_leftPreconditioner->preconditionerCreated() )||
           ( M_rightPreconditioner.get() !=0 && M_rightPreconditioner->preconditionerCreated() )||
             M_problem->isLeftPrec() || M_problem->isRightPrec();
}

void
SolverBelos::showMe( std::ostream& output ) const
{
    if(M_displayer->isLeader())
    {
        output << "SolverBelos, parameters list:" << std::endl;
        output << "-----------------------------" << std::endl;
        output << M_parameterList << endl;
        output << "-----------------------------" << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
SolverBelos::setSolverManager( const SolverManagerType& solverManager )
{
    M_solverManagerType = solverManager;
}

void
SolverBelos::setCommunicator( const boost::shared_ptr<Epetra_Comm>& comm )
{
    M_displayer->setCommunicator( comm );
}

void SolverBelos::setMatrix( matrix_ptrtype& matrix )
{
    M_matrix = matrix;
    Teuchos::RCP<Epetra_FECrsMatrix> A = Teuchos::rcp(M_matrix->matrixPtr());
    M_problem->setOperator( A );
}

void
SolverBelos::setOperator( Epetra_Operator& oper )
{
    M_problem->setOperator(rcp(&oper));
}

void
SolverBelos::setRightHandSide(const vector_type& rhs)
{
    Teuchos::RCP<const vector_type::vector_type> rhsPtr(&(rhs.epetraVector()),false);
    M_problem->setRHS(rhsPtr);
}

void
SolverBelos::setPreconditioner( prec_type& preconditioner, PrecApplicationType precType )
{
    if(precType == RightPreconditioner)
    {
        // If a right Epetra_Operator exists it must be deleted
        if(M_problem->isRightPrec())
        {
            M_problem->setRightPrec(Teuchos::RCP<OP>(null));
        }

        M_rightPreconditioner = preconditioner;
    }
    else
    {
        // If a left Epetra_Operator exists it must be deleted
        if(M_problem->isLeftPrec())
        {
            M_problem->setLeftPrec(Teuchos::RCP<OP>(null));
        }
        M_leftPreconditioner  = preconditioner;
    }
}

void
SolverBelos::setPreconditioner( comp_prec_type& preconditioner, PrecApplicationType precType )
{
    if(precType == RightPreconditioner)
    {
        // If a right LifeV::Preconditioner exists it must be deleted
        M_rightPreconditioner.reset();

        Teuchos::RCP<OP> rightPrec=Teuchos::rcp(preconditioner);
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( rightPrec ) );
        M_problem->setRightPrec(belosPrec);
    }
    else
    {
        // If a left LifeV::Preconditioner exists it must be deleted
        M_leftPreconditioner.reset();

        Teuchos::RCP<OP> leftPrec=Teuchos::rcp(preconditioner);
        Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = rcp( new Belos::EpetraPrecOp( leftPrec ) );
        M_problem->setLeftPrec(belosPrec);
    }
}

void
SolverBelos::setParameters( const GetPot& dataFile, const std::string& section )
{
    bool found;

    // Flexible Gmres will be used to solve this problem
    bool flexibleGmres = dataFile( ( section + "/flexible_gmres" ).data(), false, found );
    if ( found ) M_parameterList.set( "Flexible Gmres", flexibleGmres );

    // Relative convergence tolerance requested
    Real tolerance = dataFile( ( section + "/tol" ).data(), 1.e-6, found );
    if ( found ) M_parameterList.set( "Convergence Tolerance", tolerance );

    // Maximum number of iterations allowed
    Int maxIter = dataFile( ( section + "/max_iter"      ).data(), 200, found );
    if ( found ) M_parameterList.set( "Maximum Iterations", maxIter );

    // Output Frequency
    Int outputFrequency = dataFile( ( section + "/max_iter" ).data(), 1, found );
    if ( found ) M_parameterList.set( "Output Frequency", 1 );

    Int blockSize = dataFile( ( section + "/block_size" ).data(), 10, found );
    if ( found ) M_parameterList.set( "Block Size", blockSize );               // Blocksize to be used by iterative solver

    // Maximum number of blocks in Krylov factorization
    Int numBlocks = dataFile( ( section + "/num_blocks" ).data(), 10, found );
    if ( found ) M_parameterList.set( "Num Blocks", numBlocks );

    // Maximum number of restarts allowed
    Int maximumRestarts = dataFile( ( section + "/maximum_restarts" ).data(), 0, found );
    if ( found ) M_parameterList.set( "Maximum Restarts", maximumRestarts );

    // Set the output style (General Brief)
    std::string outputStyle = dataFile( ( section + "/output_style" ).data(), "brief", found );
    if ( found )
    {
        if (outputStyle == "brief") M_parameterList.set( "Output Style", Belos::Brief );
        if (outputStyle == "general") M_parameterList.set( "Output Style", Belos::General );
    }

    // Setting the desired output informations
    int msg = Belos::Errors;
    dataFile( ( section + "/enable_warnings" ).data(), true, found );
    if ( found ) msg += Belos::Warnings;
    dataFile( ( section + "/enable_iterations_details" ).data(), true, found );
    if ( found ) msg += Belos::IterationDetails;
    dataFile( ( section + "/enable_ortho_details" ).data(), false, found );
    if ( found ) msg += Belos::OrthoDetails;
    dataFile( ( section + "/enable_final_summary" ).data(), false, found );
    if ( found ) msg += Belos::FinalSummary;
    dataFile( ( section + "/enable_timing_details" ).data(), false, found );
    if ( found ) msg += Belos::TimingDetails;
    dataFile( ( section + "/enable_status_test_details" ).data(), false, found );
    if ( found ) msg += Belos::StatusTestDetails;
    dataFile( ( section + "/enable_debug" ).data(), false, found );
    if ( found ) msg += Belos::Debug;
    M_parameterList.set("Verbosity", msg );

    // LifeV features

    // Reuse the preconditioner from one to another call
    bool reusePreconditioner = dataFile( (section + "/reuse_preconditioner").data(), true, found );
    if ( found ) M_reusePreconditioner = reusePreconditioner;
    //if ( found ) M_parameterList.set( "Reuse preconditioner", true );

    // Max iterations allowed to reuse the preconditioner
    Int maxIterForReuse = dataFile( ( section + "/max_iter_reuse").data(), static_cast<Int> ( maxIter*8./10.), found );
    if ( found ) M_maxIterForReuse=maxIterForReuse;

    // If quitOnFailure is enabled and if some problems occur
    // the simulation is stopped
    bool quitOnFailure = dataFile( ( section + "/quit_on_failure").data(), false, found );
    if ( found ) M_quitOnFailure=quitOnFailure;
}

void
SolverBelos::setParameters( const Teuchos::ParameterList& list )
{
    M_parameterList.setParameters( list );
}

void
SolverBelos::resetParameters()
{
    M_parameterList = Teuchos::ParameterList();
}

void
SolverBelos::setReusePreconditioner( const bool reusePreconditioner )
{
    M_reusePreconditioner = reusePreconditioner;
}

void setQuitOnFailure(const bool enable)
{
    M_quitOnFailure = enable;
}

// ===================================================
// Get Methods
// ===================================================

Int
SolverBelos::numIterations() const
{
    return M_solverManager->getNumIters();
}


Real
SolverBelos::trueResidual()
{
    if( !M_problem->isProblemSet() )
    {
        M_displayer->leaderPrint( "SLV-  WARNING: SolverBelos can not compute the residual if the linear system is not set!\n");
        return -1;
    }
    MV res( *(M_problem->getRHS().get()) );
    M_problem->computeCurrResVec(&res,M_problem->getLHS().get(),M_problem->getRHS().get());
    Real residual;
    res.Norm2( &residual );
    return residual;
}

SolverBelos::prec_type&
SolverBelos::preconditioner( PrecApplicationType precType )
{
    if(precType == RightPreconditioner)
    {
        return M_rightPreconditioner;
    }
    return M_leftPreconditioner;
}

Teuchos::ParameterList&
SolverBelos::getParametersList()
{
    return M_parameterList;
}

SolverBelos::SolverManager_ptrtype
SolverBelos::solver()
{
    return M_solverManager;
}

boost::shared_ptr<Displayer>
SolverBelos::displayer()
{
    return M_displayer;
}

void
SolverBelos::setupSolverManager()
{
    // If a SolverManager already exists we simply clean it!
    if(!M_solverManager.is_null())
    {
        M_solverManager.reset();
    }

    switch(M_solverManagerType)
    {
        case BlockCG:
            // Create the block CG iteration
            M_solverManager = rcp( new Belos::BlockCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
            break;
        case BlockGmres:
            // Create the block GMRes iteration
            // Create the flexible, block GMRes iteration
            M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );
            break;
        case GCRODR:
            M_solverManager = rcp( new Belos::GCRODRSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );
            break;
        case GmresPoly:
            M_solverManager = rcp( new Belos::GmresPolySolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );
            break;
        case PCPG:
            M_solverManager = rcp( new Belos::PCPGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
            break;
        case PseudoBlockCG:
            // Create the pseudo block CG iteration
            M_solverManager = rcp( new Belos::PseudoBlockCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
            break;
        case PseudoBlockGmres:
            // Create the pseudo block GMRes iteration
            M_solverManager = rcp( new Belos::PseudoBlockGmresSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false) ) );
            break;
        case RCG:
            M_solverManager = rcp( new Belos::RCGSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
            break;
        case TFQMR:
            // Create TFQMR iteration
            M_solverManager = rcp( new Belos::TFQMRSolMgr<Real,MV,OP>( M_problem, rcp(&M_parameterList,false)) );
            break;
    }
}

} // namespace LifeV
