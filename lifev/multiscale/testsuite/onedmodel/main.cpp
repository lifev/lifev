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
 *  @file
 *  @brief File containing the One Dimensional Test
 *
 *  @date 01-09-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 01-01-2010
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  This is a test to verify that the One Dimensional Model works correctly.
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MapEpetra.hpp>

// Mathcard includes
#include <lifev/one_d_fsi/fem/OneDFSIBCHandler.hpp>
#include <lifev/one_d_fsi/solver/OneDFSISolver.hpp>
#include <lifev/multiscale/solver/MultiscaleModelFSI1D.hpp>

#include "ud_functions.hpp"

using namespace LifeV;
using namespace Multiscale;

bool checkValue(const Real val, const Real test, const Real tol = 1.e-5, const bool verbose = true)
{
    Real norm = abs(val - test);

    if ( verbose )
        std::cout << "value = " << val << " computed value = " << test << " diff = " << norm << std::endl;

    return (norm < tol);
}

Int main(Int argc, char** argv)
{
    //Setup main communicator
    boost::shared_ptr<Epetra_Comm>  comm;

#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init( &argc, &argv );
#endif

    //MPI Preprocessing
#ifdef EPETRA_MPI
    Int nprocs;
    Int rank;

    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if ( rank == 0 )
    {
        std::cout << "MPI Processes: " << nprocs << std::endl;
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

    comm->Barrier();
#else
    std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_SerialComm() );
#endif

    // *********************************
    // Useful typedefs
    // *********************************
    typedef MultiscaleModelFSI1D::physics_Type       physics_Type;
    typedef MultiscaleModelFSI1D::flux_Type          flux_Type;
    typedef MultiscaleModelFSI1D::source_Type        source_Type;

    typedef MultiscaleModelFSI1D::bc_Type            bc_Type;
    typedef bc_Type::bcFunction_Type                 bcFunction_Type;

    // *********************************
    // Reading from data file
    // *********************************
    GetPot command_line(argc,argv);

    // checking if we are checking for the nightly build
    const bool check = command_line.search(2, "-c", "--check");

    string fileName = command_line.follow("data", 2, "-f","--file");
    GetPot dataFile( fileName );

    // *********************************
    // Build the 1D model
    // *********************************
    MultiscaleModelFSI1D oneDModel;
    oneDModel.setCommunicator( comm );

    // Scale, Rotate, Translate 1D (if necessary)
    boost::array< Real, NDIM >    geometryScale;
    boost::array< Real, NDIM >    geometryRotate;
    boost::array< Real, NDIM >    geometryTranslate;

    geometryScale[0] = dataFile( "1D_Model/space_discretization/transform", 1., 0);
    geometryScale[1] = dataFile( "1D_Model/space_discretization/transform", 1., 1);
    geometryScale[2] = dataFile( "1D_Model/space_discretization/transform", 1., 2);

    geometryRotate[0] = dataFile( "1D_Model/space_discretization/transform", 0., 3) * M_PI / 180;
    geometryRotate[1] = dataFile( "1D_Model/space_discretization/transform", 0., 4) * M_PI / 180;
    geometryRotate[2] = dataFile( "1D_Model/space_discretization/transform", 0., 5) * M_PI / 180;

    geometryTranslate[0] = dataFile( "1D_Model/space_discretization/transform", 0., 6);
    geometryTranslate[1] = dataFile( "1D_Model/space_discretization/transform", 0., 7);
    geometryTranslate[2] = dataFile( "1D_Model/space_discretization/transform", 0., 8);

    oneDModel.setGeometry( geometryScale, geometryRotate, geometryTranslate );

    oneDModel.setupData( fileName );

    // *********************************
    // BC for the 1D model
    // *********************************

    // Set BC using BCInterface
    //oneDModel.GetBCInterface().FillHandler( fileName, "1D_Model" );

    // Set BC using standard approach
    Sin sinus( 0, 10, .01, 0.);
    bcFunction_Type sinusoidalFunction( boost::bind( &Sin::operator(), &sinus, _1 ) );

    // Absorbing
    bc_Type::bcFunctionSolverDefinedPtr_Type absorbing ( new OneDFSIFunctionSolverDefinedAbsorbing( OneDFSI::right, OneDFSI::W2 ) );
    bcFunction_Type absorbingFunction ( boost::bind( &OneDFSIFunctionSolverDefinedAbsorbing::operator(),
                                                     dynamic_cast<OneDFSIFunctionSolverDefinedAbsorbing *> ( &( *absorbing ) ), _1, _2 ) );

    // BC to test A_from_P conversion
    //Constant constantArea( 1.00 );
    //bcFunction_Type constantAreaFunction( boost::bind( &Constant::operator(), &constantArea, _1 ) );

    //Constant constantPressure( 24695.0765959599 );
    //bcFunction_Type constantPressureFunction( boost::bind( &Constant::operator(), &constantPressure, _1 ) );

    oneDModel.bc().setBC( OneDFSI::left,  OneDFSI::first, OneDFSI::Q,  sinusoidalFunction  );
    oneDModel.bc().setBC( OneDFSI::right, OneDFSI::first, OneDFSI::W2, absorbingFunction );

    //oneDModel.bc().setBC( OneDFSI::right, OneDFSI::first, OneDFSI::A, constantAreaFunction );
    //oneDModel.bc().setBC( OneDFSI::right, OneDFSI::first, OneDFSI::P,   constantPressureFunction );

    oneDModel.setupModel();

    absorbing->setSolution( oneDModel.solution() );
    absorbing->setFluxSource( oneDModel.flux(), oneDModel.source() );

    // *********************************
    // Tempolar loop
    // *********************************
    std::cout << "\nTemporal loop:" << std::endl;

    LifeChrono chronoTotal;
    LifeChrono chronoSystem;
    LifeChrono chronoIteration;

    Int count = 0;
    chronoTotal.start();
    for ( ; oneDModel.data().dataTime()->canAdvance() ; oneDModel.data().dataTime()->updateTime(), ++count )
    {
        std::cout << std::endl << "--------- Iteration " << count << " time = " << oneDModel.data().dataTime()->time() << std::endl;

        chronoIteration.start();

        if ( oneDModel.data().dataTime()->isFirstTimeStep() )
            oneDModel.buildModel();
        else
            oneDModel.updateModel();

        chronoSystem.start();

        oneDModel.solveModel();

        chronoSystem.stop();

        oneDModel.updateSolution();

        //Save solution
        if ( count%50 == 0 || oneDModel.data().dataTime()->isLastTimeStep() )
            oneDModel.saveSolution();

        chronoIteration.stop();

        std::cout << " System solved in " << chronoSystem.diff() << " s, (total time " << chronoIteration.diff() << " s)." << std::endl;
    }

    chronoTotal.stop();
    std::cout << std::endl << " Simulation ended successfully in " << chronoTotal.diff()  << " s" << std::endl;

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    if ( check )
    {
        bool ok = true;
        Int rightNodeID = oneDModel.solver()->boundaryDOF( OneDFSI::right );

        ok = ok && checkValue( 0.999998  , (*oneDModel.solution("A"))[rightNodeID]);
        ok = ok && checkValue(-0.00138076, (*oneDModel.solution("Q"))[rightNodeID]);
        ok = ok && checkValue(-0.00276153, (*oneDModel.solution("W1"))[rightNodeID]);
        ok = ok && checkValue( 0.00000000, (*oneDModel.solution("W2"))[rightNodeID]);

        ok = ok && checkValue( 0.999999  , (*oneDModel.solution("A"))[rightNodeID - 1]);
        ok = ok && checkValue(-0.00040393, (*oneDModel.solution("Q"))[rightNodeID - 1]);
        ok = ok && checkValue(-0.00080833, (*oneDModel.solution("W1"))[rightNodeID - 1]);
        ok = ok && checkValue( 0.00000045, (*oneDModel.solution("W2"))[rightNodeID - 1]);

        if (ok)
        {
            std::cout << "Test successful" << std::endl;
            return 0;
        }
        else
        {
            std::cout << "Test unsuccessful" << std::endl;
            return -1;
        }
    }
    else
        return 0;
}
