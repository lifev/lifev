/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-28

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Solver.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-28
 */

#include <lifemc/lifesolver/MS_Solver.hpp>

namespace LifeV {

// ===================================================
//! Constructors
// ===================================================
MS_Solver::MS_Solver() :
    M_multiscale        (),
    M_dataPhysics       ( new MS_PhysicalData() ),
    M_dataTime          (),
    M_comm              (),
    M_displayer         (),
    M_chrono            (),
    M_couplingVariables (),
    M_couplingResiduals (),
    M_generalizedAitken (),
    M_subITMax          (),
    M_tolerance         ()
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::MS_Solver() \n";
#endif

    //Models & Couplings mapCreation
    modelsMap["MultiScale"]           = MultiScale;
    modelsMap["Fluid3D"]              = Fluid3D;
    couplingsMap["BoundaryCondition"] = BoundaryCondition;
    couplingsMap["FluxStress"]        = FluxStress;

    //Models & Couplings Factory registration
    FactoryModels::instance().registerProduct   ( MultiScale,        &createMultiScale );
    FactoryModels::instance().registerProduct   ( Fluid3D,           &createFluid3D );
    FactoryCouplings::instance().registerProduct( FluxStress,        &createFluxStress );
    FactoryCouplings::instance().registerProduct( BoundaryCondition, &createBoundaryCondition );
}

MS_Solver::MS_Solver( const MS_Solver& solver ) :
    M_multiscale        ( solver.M_multiscale ),
    M_dataPhysics       ( solver.M_dataPhysics ),
    M_dataTime          ( solver.M_dataTime ),
    M_comm              ( solver.M_comm ),
    M_displayer         ( solver.M_displayer ),
    M_chrono            ( solver.M_chrono ),
    M_couplingVariables ( solver.M_couplingVariables ),
    M_couplingResiduals ( solver.M_couplingResiduals ),
    M_generalizedAitken ( solver.M_generalizedAitken ),
    M_subITMax          ( solver.M_subITMax ),
    M_tolerance         ( solver.M_tolerance )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::MS_Solver( solver ) \n";
#endif

}

// ===================================================
//! Methods
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
        M_dataPhysics       = solver.M_dataPhysics;
        M_dataTime          = solver.M_dataTime;
        M_comm              = solver.M_comm;
        M_displayer         = solver.M_displayer;
        M_chrono            = solver.M_chrono;
        M_couplingVariables = solver.M_couplingVariables;
        M_couplingResiduals = solver.M_couplingResiduals;
        M_generalizedAitken = solver.M_generalizedAitken;
        M_subITMax          = solver.M_subITMax;
        M_tolerance         = solver.M_tolerance;
    }

    return *this;
}

void
MS_Solver::SetCommunicator( const boost::shared_ptr< Epetra_Comm >& comm )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SetCommunicator( comm ) \n";
#endif

    M_comm = comm;
    M_multiscale.SetCommunicator( M_comm );
    M_displayer.reset( new Displayer( M_comm.get() ) );
}

void
MS_Solver::SetupProblem( const std::string& dataFile )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SetupData( dataFile ) \n";
#endif

    GetPot DataFile( dataFile );

    // Main MultiScale problem
    M_multiscale.SetDataFile( DataFile( "Problem/MS_problem", "./MultiScaleData/Models/Model.dat" ) );

    // Time & Physics containers
    M_dataTime.reset( new DataTime( DataFile, "Algorithm/time_discretization" ) ); //Add here Aitken??
    M_dataPhysics->ReadData( DataFile );

    // Sub-Iterations parameters
    M_generalizedAitken.setDefault( DataFile( "Algorithm/Aitken_method/omega", 1.e-3 ) );
    M_subITMax  = DataFile( "Algorithm/Aitken_method/subITMax", 100 );
    M_tolerance = DataFile( "Algorithm/Aitken_method/tolerance", 1.e-10 );

    // Setup MultiScale problem
    M_multiscale.SetData( M_dataPhysics, M_dataTime );
    M_multiscale.SetupData();
    M_multiscale.SetupModel();

    // Build coupling variables and residuals vectors
    M_multiscale.SetupImplicitCoupling( M_couplingVariables, M_couplingResiduals );
}

void
MS_Solver::SolveProblem( void )
{

#ifdef DEBUG
    Debug( 8000 ) << "MS_Solver::SolveProblem() \n";
#endif

    Real residual;

    for ( ; M_dataTime->canAdvance(); M_dataTime->updateTime() )
    {
        M_chrono.start();

        if ( M_displayer->isLeader() )
        {
            std::cout << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
            std::cout << "            MULTISCALE SIMULATION TIME: " << M_dataTime->getTime() << " s" << std::endl;
            std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl << std::endl;
        }

        // Build or Update System
        if ( M_dataTime->isFirstTimeStep() )
            M_multiscale.BuildSystem();
        else
            M_multiscale.UpdateSystem();

        // SolveSystem
        M_multiscale.SolveSystem();

        // Temporary Computation of a Block Vector
        //VectorType blocksVector( M_couplingVariables ); blocksVector = 0.0;
        //for ( UInt i = 1 ; i < blocksVector.size() ; i = i+2)
        //    blocksVector[i] = 1.0;
        //std::cout << "blocksVector: " << std::endl;
        //blocksVector.ShowMe();

        // Sub-Iterations with Aitken method
        M_generalizedAitken.restart( true );
        for ( UInt subIT = 0; subIT < M_subITMax; ++subIT )
        {
            residual = M_couplingResiduals.WeightNorm2();

            if ( M_displayer->isLeader() )
            {
                std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;
                std::cout << " MS-  Residual:                                " << residual << std::endl;
            }

            // To be moved in a post-processing class
            std::cout << " MS-  CouplingVariables:\n" << std::endl;
            M_couplingVariables.ShowMe();
            std::cout << " MS-  CouplingResiduals:\n" << std::endl;
            M_couplingResiduals.ShowMe();

            // Verify tolerance
            if ( residual <= M_tolerance )
                break;

            // Update Coupling Variables
            //M_generalizedAitken.restart(); // To have fixed omega (if omega = 1, Fixed Point)
            //M_couplingVariables += M_generalizedAitken.computeDeltaLambdaScalar( M_couplingVariables, M_couplingResiduals, true );
            M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVector( M_couplingVariables, M_couplingResiduals, true, true );
            //M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVectorBlock( M_couplingVariables, M_couplingResiduals, blocksVector, 2, true );

            // SolveSystem
            M_multiscale.SolveSystem();
        }

        // SaveSolution
        M_multiscale.SaveSolution();

        M_chrono.stop();

        if ( M_displayer->isLeader() )
            std::cout << " MS-  Total iteration time:                    " << M_chrono.diff() << " s" << std::endl;
    }
}

void
MS_Solver::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        std::cout << std::endl << std::endl
                  << "=============== MultiScale Solver Information ===============" << std::endl << std::endl;

        std::cout << "Initial time       = " << M_dataTime->getInitialTime() << std::endl
                  << "End time           = " << M_dataTime->getEndTime() << std::endl
                  << "TimeStep           = " << M_dataTime->getTimeStep() << std::endl << std::endl;
        std::cout << std::endl << std::endl;

        std::cout << "Max Sub-iterations = " << M_subITMax << std::endl
                  << "Tolerance          = " << M_tolerance << std::endl << std::endl;
        std::cout << std::endl << std::endl;

        M_multiscale.ShowMe();

        std::cout << "=============================================================" << std::endl << std::endl;
    }
}

} // Namespace LifeV
