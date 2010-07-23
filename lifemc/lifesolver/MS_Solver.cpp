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

std::string MS_ProblemFolder = "./";
UInt        MS_ProblemStep   = 0;
bool        MS_ExitFlag      = EXIT_SUCCESS;

// ===================================================
// Constructors
// ===================================================
MS_Solver::MS_Solver() :
    M_model             (),
    M_algorithm         (),
    M_globalData        ( new MS_GlobalDataContainer_Type() ),
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
    MS_Model_Factory::instance().registerProduct   (  MultiScale,          &MS_createMultiScale );
    MS_Model_Factory::instance().registerProduct   (  Fluid3D,             &MS_createFluid3D );
    MS_Model_Factory::instance().registerProduct   (  OneDimensional,      &MS_createOneDimensional );
    MS_Model_Factory::instance().registerProduct   (  FSI3D,               &MS_createModelFSI3D );

    MS_Coupling_Factory::instance().registerProduct(  Stress,              &MS_createStress );
    MS_Coupling_Factory::instance().registerProduct(  FluxStress,          &MS_createFluxStress );
    MS_Coupling_Factory::instance().registerProduct(  BoundaryCondition,   &MS_createBoundaryCondition );

    MS_Algorithm_Factory::instance().registerProduct( Aitken,              &MS_createAitken );
    MS_Algorithm_Factory::instance().registerProduct( Explicit,            &MS_createExplicit );
    MS_Algorithm_Factory::instance().registerProduct( Newton,              &MS_createNewton );
}

MS_Solver::MS_Solver( const MS_Solver& solver ) :
    M_model             ( solver.M_model ),
    M_algorithm         ( solver.M_algorithm ),
    M_globalData        ( solver.M_globalData ),
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
        M_model             = solver.M_model;
        M_algorithm         = solver.M_algorithm;
        M_globalData        = solver.M_globalData;
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

    // Create the main model and set the communicator
    M_model = MS_Model_PtrType( MS_Model_Factory::instance().createObject( MS_modelsMap[ DataFile( "Problem/ProblemType", "MultiScale" ) ] ) );
    M_model->SetCommunicator( M_comm );

    // Setup data
    M_globalData->ReadData( DataFile );
    M_model->SetGlobalData( M_globalData );
    M_model->SetupData( DataFile( "Problem/ProblemFile", "./MultiScaleData/Models/Model.dat" ) );

    // Setup Models
    M_model->SetupModel();

    // Algorithm parameters
    if ( M_model->GetType() == MultiScale )
    {
        M_algorithm = MS_Algorithm_PtrType( MS_Algorithm_Factory::instance().createObject( MS_algorithmsMap[ DataFile( "Solver/Algorithm/AlgorithmType", "Newton" ) ] ) );
        M_algorithm->SetCommunicator( M_comm );
        M_algorithm->SetModel( M_model );
        M_algorithm->SetupData( FileName );
    }
}

bool
MS_Solver::SolveProblem()
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SolveProblem() \n";
#endif

    // Save initial solution if it is the very first time step
    if ( !MS_ProblemStep )
        M_model->SaveSolution();

    // Move to the "true" first time-step
    M_globalData->GetDataTime()->updateTime();
    M_globalData->GetDataTime()->setInitialTime( M_globalData->GetDataTime()->getTime() );

    for ( ; M_globalData->GetDataTime()->canAdvance(); M_globalData->GetDataTime()->updateTime() )
    {
        M_chrono.start();

        if ( M_displayer->isLeader() )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "                    MULTISCALE SIMULATION" << std::endl;
            std::cout << "             time = " << M_globalData->GetDataTime()->getTime() << " s; "  <<
                          "time step number = " << M_globalData->GetDataTime()->getTimeStepNumber()  << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        if ( M_globalData->GetDataTime()->isFirstTimeStep() )
            M_model->BuildSystem();
        else
            M_model->UpdateSystem();

        // SolveSystem
        M_model->SolveSystem();

        // If it is a MultiScale model, call algorithms for subiterations
        if ( M_model->GetType() == MultiScale )
            M_algorithm->SubIterate();

        // SaveSolution
        M_model->SaveSolution();

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

        std::cout << "Problem folder                = " << MS_ProblemFolder << std::endl
                  << "Problem step                  = " << MS_ProblemStep << std::endl << std::endl;

        M_globalData->ShowMe();

        std::cout << std::endl << std::endl;
    }

    M_model->ShowMe();
    if ( M_model->GetType() == MultiScale )
        M_algorithm->ShowMe();

    if ( M_displayer->isLeader() )
        std::cout << "=============================================================" << std::endl << std::endl;
}

} // Namespace LifeV
