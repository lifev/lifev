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
    @brief BelosOperator

    @author Umberto Villa <umberto.villa@gmail.com>

    @date 28-09-2010
 */

#include<lifev/core/operator/BelosOperator.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"

#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosGmresPolySolMgr.hpp>
#include <BelosPCPGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#ifdef HAVE_TRILINOS_GT_10_6
#include <BelosMinresSolMgr.hpp>
#endif
#include <BelosRCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"
#pragma GCC diagnostic warning "-Wunused-but-set-variable"

namespace LifeV
{
namespace Operators
{

BelosOperator::BelosOperator():
        SolverOperator(),
        M_linProblem( Teuchos::rcp( new LinearProblem ) )
{
    M_name = "BelosOperator";
}

BelosOperator::~BelosOperator()
{
    M_belosPrec = Teuchos::null;
    M_linProblem = Teuchos::null;
    M_solverManager = Teuchos::null;
}

int BelosOperator::doApplyInverse( const vector_Type& X, vector_Type& Y ) const
{

    Teuchos::RCP<vector_Type> Xcopy( new vector_Type( X ) );
    Y.PutScalar( 0.0 );
    bool set = M_linProblem->setProblem( Teuchos::rcp( &Y, false ), Xcopy );
    if ( set == false )
    {
        std::cout << std::endl << "SLV-  ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
        return -12;
    }

    M_solverManager->setProblem( M_linProblem );


    // Solving the system
    Belos::ReturnType ret = M_solverManager->solve();

    // Update the number of performed iterations
    M_numIterations = M_solverManager->getNumIters();

    // Update of the status
    if( M_solverManager->isLOADetected() )
    {
        M_lossOfAccuracy = yes;
    }
    else
    {
        M_lossOfAccuracy = no;
    }

    if( ret == Belos::Converged )
    {
        M_converged = yes;
        return 0;
    }
    else
    {
        M_converged = no;
        return -1;
    }

}

void BelosOperator::doSetOperator()
{
    Teuchos::RCP<OP> tmpPtr( M_oper.get(), false );
    M_linProblem->setOperator( tmpPtr );
}

void BelosOperator::doSetPreconditioner()
{
    Teuchos::RCP<OP> tmpPtr( M_prec.get(), false );
    M_belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( tmpPtr ), true );

    // The line below produces a memory leak; It has been kept as an example to illustrate
    // why it has been changed.
    // M_belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( Teuchos::rcp( M_prec ) ), false );
}

void BelosOperator::doSetParameterList()
{
    if( !M_pList->sublist( "Trilinos: Belos List" ).isParameter( "Verbosity" ) )
         M_pList->sublist( "Trilinos: Belos List" ).set( "Verbosity", Belos::Errors + Belos::Warnings +
                                                         Belos::TimingDetails + Belos::StatusTestDetails );

    if( M_tolerance > 0 )
        M_pList->sublist( "Trilinos: Belos List" ).set( "Convergence Tolerance", M_tolerance );

    std::string solverType( M_pList->get<std::string>( "Solver Manager Type" ) );
    allocateSolver( getSolverManagerTypeFromString( solverType ) );
    M_solverManager->setParameters( sublist( M_pList, "Trilinos: Belos List", true ) );

    std::string precSideStr( M_pList->get<std::string>( "Preconditioner Side" ) );
    PreconditionerSide precSide( getPreconditionerSideFromString( precSideStr ) );

    switch(precSide)
    {
    case None:
        break;
    case Left:
        M_linProblem->setLeftPrec( M_belosPrec );
        break;
    case Right:
        M_linProblem->setRightPrec( M_belosPrec );
        break;
    default:
        exit(1);
    }

}

//============================================================================//
//                     Protected or Private Methods                           //
//============================================================================//
void BelosOperator::allocateSolver( const SolverManagerType & solverManagerType )
{
       // If a SolverManager already exists we simply clean it!
        if ( !M_solverManager.is_null() )
        {
            M_solverManager = Teuchos::null;
        }

        switch ( solverManagerType )
        {
            case NotAValidSolverManager:
                std::cout<<"SLV-  ERROR: Not a Valid Solver Manager \n";
                exit(1);
                break;
            case BlockCG:
                // Create the block CG iteration
                M_solverManager = rcp( new Belos::BlockCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case PseudoBlockCG:
                // Create the pseudo block CG iteration
                M_solverManager = rcp( new Belos::PseudoBlockCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case RCG:
                M_solverManager = rcp( new Belos::RCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case BlockGmres:
                M_solverManager = rcp( new Belos::BlockGmresSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case PseudoBlockGmres:
                M_solverManager = rcp( new Belos::PseudoBlockGmresSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case GmresPoly:
                M_solverManager = rcp( new Belos::GmresPolySolMgr<Real,vector_Type,operator_Type>() );
                break;
            case GCRODR:
                M_solverManager = rcp( new Belos::GCRODRSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case PCPG:
                M_solverManager = rcp( new Belos::PCPGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case TFQMR:
                // Create TFQMR iteration
                M_solverManager = rcp( new Belos::TFQMRSolMgr<Real,vector_Type,operator_Type>() );
                break;
#ifdef HAVE_TRILINOS_GT_10_6
            case MINRES:
                // Create MINRES iteration
                M_solverManager = rcp( new Belos::MinresSolMgr<Real,vector_Type,operator_Type>() );
                break;
#endif
            default:
                ERROR_MSG("Belos solver not found!");
         }

}

BelosOperator::SolverManagerType
BelosOperator::getSolverManagerTypeFromString ( const std::string& str )
{
    if( str == "BlockCG" )
        return BlockCG;
    else if( str == "PseudoBlockCG" )
        return PseudoBlockCG;
    else if( str == "RCG" )
        return RCG;
    else if( str == "BlockGmres" )
        return BlockGmres;
    else if( str == "PseudoBlockGmres" )
        return PseudoBlockGmres;
    else if( str == "GmresPoly" )
        return GmresPoly;
    else if( str == "GCRODR" )
        return GCRODR;
    else if( str == "PCPG" )
        return PCPG;
    else if( str == "TFQMR" )
        return TFQMR;
    else if( str == "MINRES" )
        return MINRES;
    else
        return NotAValidSolverManager;
}

BelosOperator::PreconditionerSide
BelosOperator::getPreconditionerSideFromString( const std::string& str )
{
    if( str == "Right" || str == "right" )
        return Right;
    else if( str == "Left" || str == "left" )
        return Left;
    else
        return None;
}

void
BelosOperator::doResetSolver()
{
    M_solverManager->reset(Belos::Problem);
    M_belosPrec = Teuchos::null;
}

} // Namespace Operators

} // Namespace LifeV
