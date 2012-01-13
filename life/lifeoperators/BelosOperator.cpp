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

#include<life/lifeoperators/BelosOperator.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wextra"

#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosGmresPolySolMgr.hpp>
#include <BelosPCPGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosRCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Wextra"

namespace LifeV
{
namespace Operators
{

std::auto_ptr<BelosOperator::solverManagerMap_Type> BelosOperator::S_solverManagerMap(BelosOperator::singletonSolverManagerMap());
std::auto_ptr<BelosOperator::precSideMap_Type> BelosOperator::S_precSideMap(BelosOperator::singletonPrecSideMap());

BelosOperator::BelosOperator():
		SolverOperator(),
		M_linProblem(Teuchos::rcp(new LinearProblem))
{
	M_name = "BelosOperator";
}

int BelosOperator::doApplyInverse(const vector_Type& X, vector_Type& Y)
{

	Teuchos::RCP<vector_Type> Xcopy(new vector_Type(X) );
	Y.PutScalar(0.0);
	bool set = M_linProblem->setProblem(Teuchos::rcp(&Y, false), Xcopy);
	if (set == false)
	{
		std::cout << std::endl << "SLV-  ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
		return -12;
	}

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

	if(ret == Belos::Converged)
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
	M_linProblem->setOperator(M_oper);
}

void BelosOperator::doSetPreconditioner()
{
	M_belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( M_prec ) );
}

void BelosOperator::doSetParameterList()
{
	if( !M_pList->sublist( "Belos List" ).isParameter( "Verbosity" ) )
		 M_pList->sublist( "Belos List" ).set( "Verbosity", Belos::Errors + Belos::Warnings +
											             Belos::TimingDetails + Belos::StatusTestDetails );

	std::string solverType(M_pList->get<std::string>("Solver Manager Type"));
	allocateSolver( (*S_solverManagerMap)[solverType]);
	M_solverManager->setParameters(sublist(M_pList, "Belos List", true));

	std::string precSideStr( M_pList->get<std::string>("Preconditioner Side"));
	PreconditionerSide precSide((*S_precSideMap)[precSideStr]);

	switch(precSide)
	{
	case None:
		break;
	case Left:
		M_linProblem->setLeftPrec(M_belosPrec);
		break;
	case Right:
		M_linProblem->setRightPrec(M_belosPrec);
		break;
	default:
		exit(1);
	}

	M_solverManager->setProblem(M_linProblem);

}

//============================================================================//
//                     Protected or Private Methods                           //
//============================================================================//
void BelosOperator::allocateSolver(const SolverManagerType & solverManagerType)
{
	   // If a SolverManager already exists we simply clean it!
	    if ( !M_solverManager.is_null() )
	    {
	        M_solverManager.reset();
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
	     }

}

BelosOperator::solverManagerMap_Type * BelosOperator::singletonSolverManagerMap()
{
	solverManagerMap_Type * map(new solverManagerMap_Type);
    (*map)["BlockCG"] = BlockCG;
    (*map)["PseudoBlockCG"] = PseudoBlockCG;
    (*map)["RCG"] = RCG;
    (*map)["BlockGmres"] = BlockGmres;
    (*map)["PseudoBlockGmres"] = PseudoBlockGmres;
    (*map)["GmresPoly"] = GmresPoly;
    (*map)["GCRODR"] = GCRODR;
    (*map)["PCPG"] = PCPG;
    (*map)["TFQMR"] = TFQMR;

    return map;
}

BelosOperator::precSideMap_Type * BelosOperator::singletonPrecSideMap()
{
	precSideMap_Type * map(new precSideMap_Type);
    (*map)["None"] = None;
    (*map)["Right"] = Right;
    (*map)["Left"] = Left;

    (*map)["none"] = None;
    (*map)["right"] = Right;
    (*map)["left"] = Left;

    return map;
}


} /* end namespace Operators */
} /* end namespace */
