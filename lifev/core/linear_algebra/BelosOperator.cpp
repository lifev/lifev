/*
 * BelosOperator.cpp
 *
 *  Created on: Sep 28, 2010
 *      Author: uvilla
 */

#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosGmresPolySolMgr.hpp>
#include <BelosPCPGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosRCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>

#include<lifev/operator/linear_algebra/BelosOperator.hpp>

namespace LifeV
{
namespace Operators
{

std::auto_ptr<BelosOperator::solverManagerMap_Type> BelosOperator::S_solverManagerMap(BelosOperator::singletonSolverManagerMap());
std::auto_ptr<BelosOperator::precSideMap_Type> BelosOperator::S_precSideMap(BelosOperator::singletonPrecSideMap());

BelosOperator::BelosOperator():
        InvertibleOperator(),
        M_linProblem(Teuchos::rcp(new LinearProblem))
{
    M_name = "BelosOperator";
}

int BelosOperator::doApplyInverse(const vector_Type& X, vector_Type& Y) const
{

    Teuchos::RCP<vector_Type> Xcopy(new vector_Type(X) );
    Y.PutScalar(0.0);
    bool set = M_linProblem->setProblem(Teuchos::rcp(&Y, false), Xcopy);
    if (set == false)
    {
        std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
        return -12;
    }

    Belos::ReturnType ret = M_solverManager->solve();

    if(ret == Belos::Converged)
        return 0;
    else
        return -1;

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
    if(! M_pList->sublist("options").isParameter("Verbosity"))
        M_pList->sublist("options").set( "Verbosity", Belos::Errors + Belos::Warnings +
                Belos::TimingDetails + Belos::StatusTestDetails );

    std::string solverType(M_pList->get<std::string>("Solver Type"));
    allocateSolver( (*S_solverManagerMap)[solverType]);
    M_solverManager->setParameters(sublist(M_pList, "options", true));

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
                std::cout<<"Not a Valid Solver Manager \n";
                exit(1);
                break;
            case BlockCG:
                // Create the block CG iteration
                M_solverManager = Teuchos::rcp( new Belos::BlockCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case PseudoBlockCG:
                // Create the pseudo block CG iteration
                M_solverManager = Teuchos::rcp( new Belos::PseudoBlockCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case RCG:
                M_solverManager = Teuchos::rcp( new Belos::RCGSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case BlockGmres:
                M_solverManager = Teuchos::rcp( new Belos::BlockGmresSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case PseudoBlockGmres:
                M_solverManager = Teuchos::rcp( new Belos::PseudoBlockGmresSolMgr<Real,vector_Type,operator_Type>() );
                break;
            case GmresPoly:
                M_solverManager = Teuchos::rcp( new Belos::GmresPolySolMgr<Real,vector_Type,operator_Type>() );
                 break;
             case GCRODR:
                 M_solverManager = Teuchos::rcp( new Belos::GCRODRSolMgr<Real,vector_Type,operator_Type>() );
                 break;
             case PCPG:
                 M_solverManager = Teuchos::rcp( new Belos::PCPGSolMgr<Real,vector_Type,operator_Type>() );
                 break;
             case TFQMR:
                 // Create TFQMR iteration
                 M_solverManager = Teuchos::rcp( new Belos::TFQMRSolMgr<Real,vector_Type,operator_Type>() );
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
