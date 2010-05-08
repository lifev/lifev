//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @file
 *  @brief MultiScale Solver
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 28-09-2009
 */

#include <lifemc/lifesolver/MS_Solver.hpp>

namespace LifeV {

std::string MS_ProblemFolder = "";
UInt        MS_ProblemStep   = 0;
bool        MS_ExitFlag      = EXIT_SUCCESS;

// ===================================================
// Constructors
// ===================================================
MS_Solver::MS_Solver() :
    M_multiscale        ( new MS_Model_MultiScale() ),
    M_algorithm         (),
    M_dataPhysics       ( new MS_PhysicalData() ),
    M_comm              (),
    M_displayer         (),
    M_chrono            ()
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::MS_Solver() \n";
#endif

    //Define the maps of MS objects
    MS_MapsDefinition();

    //Register the objects
    FactoryModels::instance().registerProduct   ( MultiScale,        &createMultiScale );
    FactoryModels::instance().registerProduct   ( Fluid3D,           &createFluid3D );
    FactoryModels::instance().registerProduct   ( FSI1D,             &createFSI1D );
    FactoryCouplings::instance().registerProduct( Stress,            &createStress );
    FactoryCouplings::instance().registerProduct( FluxStress,        &createFluxStress );
    FactoryCouplings::instance().registerProduct( BoundaryCondition, &createBoundaryCondition );
    FactoryAlgorithms::instance().registerProduct( Aitken,           &createAitken );
    FactoryAlgorithms::instance().registerProduct( Newton,           &createNewton );
}

MS_Solver::MS_Solver( const MS_Solver& solver ) :
    M_multiscale        ( solver.M_multiscale ),
    M_algorithm         ( solver.M_algorithm ),
    M_dataPhysics       ( solver.M_dataPhysics ),
    M_comm              ( solver.M_comm ),
    M_displayer         ( solver.M_displayer ),
    M_chrono            ( solver.M_chrono )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::MS_Solver( solver ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Solver&
MS_Solver::operator=( const MS_Solver& solver )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::operator=( solver ) \n";
#endif

    if ( this != &solver )
    {
        M_multiscale        = solver.M_multiscale;
        M_algorithm         = solver.M_algorithm;
        M_dataPhysics       = solver.M_dataPhysics;
        M_comm              = solver.M_comm;
        M_displayer         = solver.M_displayer;
        M_chrono            = solver.M_chrono;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================
void
MS_Solver::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_multiscale->SetCommunicator( M_comm );
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

void
MS_Solver::SetupProblem( const std::string& FileName, const std::string& problemFolder )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SetupData( FileName, problemFolder ) \n";
#endif

    // Load data file
    GetPot DataFile( FileName );

    // Define the folder containing the problem
    MS_ProblemFolder = problemFolder;

    // Define the step of the problem
    if ( DataFile( "Solver/Restart/Restart", false ) )
        MS_ProblemStep = DataFile( "Solver/Restart/RestartFromStepNumber", 0 ) + 1;

    // Time & Physics containers
    M_dataPhysics->ReadData( DataFile );
    M_multiscale->SetGlobalData( M_dataPhysics );

    // Setup data from data file
    M_multiscale->SetupData( DataFile( "Problem/MS_problem", "./MultiScaleData/Models/Model.dat" ) );

    // Setup Models
    M_multiscale->SetupModel();

    // Algorithm parameters
    M_algorithm = Algorithm_ptrType( FactoryAlgorithms::instance().createObject( algorithmMap[ DataFile( "Solver/Algorithm/AlgorithmType", "Aitken" ) ] ) );
    M_algorithm->SetCommunicator( M_comm );
    M_algorithm->SetMultiScaleProblem( M_multiscale );
    M_algorithm->SetupData( FileName );
}

bool
MS_Solver::SolveProblem()
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SolveProblem() \n";
#endif

    // Move to the "true" first time step when restarting a simulation
    if ( MS_ProblemStep > 0 )
    {
        M_dataPhysics->GetDataTime()->updateTime();
        M_dataPhysics->GetDataTime()->setInitialTime( M_dataPhysics->GetDataTime()->getTime() );
    }

    for ( ; M_dataPhysics->GetDataTime()->canAdvance(); M_dataPhysics->GetDataTime()->updateTime() )
    {
        M_chrono.start();

        if ( M_displayer->isLeader() )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE SIMULATION" << std::endl;
            std::cout << "             time = " << M_dataPhysics->GetDataTime()->getTime() << " s; "  <<
                          "time step number = " << M_dataPhysics->GetDataTime()->getTimeStepNumber()  << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        if ( M_dataPhysics->GetDataTime()->isFirstTimeStep() )
            M_multiscale->BuildSystem();
        else
            M_multiscale->UpdateSystem();

        // SolveSystem
        M_multiscale->SolveSystem();

        // Algorithm - SubIterate
        M_algorithm->SubIterate();

        // SaveSolution
        M_multiscale->SaveSolution();

        if ( M_algorithm->GetSubiterationsMaximumNumber() == 0 ) // If we use an explicit coupling algorithm
            M_multiscale->InitializeCouplingVariables();         // we need to manually update the coupling variables

        M_chrono.stop();

        if ( M_displayer->isLeader() )
            std::cout << " MS-  Total iteration time:                    " << M_chrono.diff() << " s" << std::endl;
    }

    return MS_ExitFlag;
}

void
MS_Solver::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        std::cout << std::endl << std::endl
                  << "=============== MultiScale Solver Information ===============" << std::endl << std::endl;

        std::cout << "Problem folder      = " << MS_ProblemFolder << std::endl
                  << "Problem step        = " << MS_ProblemStep << std::endl << std::endl;

        M_dataPhysics->ShowMe();

        std::cout << std::endl << std::endl;
    }

    M_multiscale->ShowMe();
    M_algorithm->ShowMe();

    if ( M_displayer->isLeader() )
        std::cout << "=============================================================" << std::endl << std::endl;
}

} // Namespace LifeV
