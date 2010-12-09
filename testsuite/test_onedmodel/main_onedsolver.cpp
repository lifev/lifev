//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
 *  @brief OneDimensionalModel Test
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-09-2004
 *
 *  @version 2.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 01-01-2010
 */

#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifealg/EpetraMap.hpp>

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

#include <lifemc/lifesolver/MS_Model_1D.hpp>

#include "ud_functions.hpp"

using namespace LifeV;

bool checkValue(const double val, const double test, const double tol = 1.e-5, const bool verbose = true)
{
    Real norm = abs(val - test);

    if ( verbose )
        std::cout << "value = " << val << " computed value = " << test << " diff = " << norm << std::endl;

    return (norm < tol);
}

int main(int argc, char** argv)
{
    //Setup main communicator
    boost::shared_ptr<Epetra_Comm>  comm;

#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init( &argc, &argv );
#endif

    //MPI Preprocessing
#ifdef EPETRA_MPI
    int nprocs;
    int rank;

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
    typedef MS_Model_1D::Physics_Type                Physics_Type;
    typedef MS_Model_1D::Flux_Type                   Flux_Type;
    typedef MS_Model_1D::Source_Type                 Source_Type;

    typedef MS_Model_1D::BC_Type                     BC_Type;
    typedef BC_Type::BCFunction_Type                 BCFunction_Type;

    // *********************************
    // Reading from data file
    // *********************************
    GetPot command_line(argc,argv);

    // checking if we are checking for the nightly build
    const bool check = command_line.search(2, "-c", "--check");

    string FileName = command_line.follow("data", 2, "-f","--file");
    GetPot DataFile(FileName);

    // *********************************
    // Build the 1D model
    // *********************************
    MS_Model_1D OneDModel;
    OneDModel.SetCommunicator(comm);

    // Scale, Rotate, Translate 1D (if necessary)
    boost::array< Real, NDIM >    geometryScale;
    boost::array< Real, NDIM >    geometryRotate;
    boost::array< Real, NDIM >    geometryTranslate;

    geometryScale[0] = DataFile( "1D_Model/space_discretization/transform", 1., 0);
    geometryScale[1] = DataFile( "1D_Model/space_discretization/transform", 1., 1);
    geometryScale[2] = DataFile( "1D_Model/space_discretization/transform", 1., 2);

    geometryRotate[0] = DataFile( "1D_Model/space_discretization/transform", 0., 3) * Pi / 180;
    geometryRotate[1] = DataFile( "1D_Model/space_discretization/transform", 0., 4) * Pi / 180;
    geometryRotate[2] = DataFile( "1D_Model/space_discretization/transform", 0., 5) * Pi / 180;

    geometryTranslate[0] = DataFile( "1D_Model/space_discretization/transform", 0., 6);
    geometryTranslate[1] = DataFile( "1D_Model/space_discretization/transform", 0., 7);
    geometryTranslate[2] = DataFile( "1D_Model/space_discretization/transform", 0., 8);

    OneDModel.SetGeometry( geometryScale, geometryRotate, geometryTranslate );

    OneDModel.SetupData( FileName );

    // *********************************
    // BC for the 1D model
    // *********************************

    // Set BC using BCInterface
    //OneDModel.GetBCInterface().FillHandler( FileName, "1D_Model" );

    // Set BC using standard approach
    Sin sinus;
    BCFunction_Type sinusoidalFunction( boost::bind( &Sin::operator(), &sinus, _1 ) );

    // Absorbing
    BC_Type::BCFunction_Default_PtrType absorbing ( new OneDimensionalModel_BCFunction_Absorbing( OneD_right, OneD_W2 ) );
    absorbing->setSolution( OneDModel.GetSolution() );
    absorbing->setFluxSource( OneDModel.GetFlux(), OneDModel.GetSource() );

    BCFunction_Type absorbingFunction ( boost::bind( &OneDimensionalModel_BCFunction_Absorbing::operator(),
                                                     dynamic_cast<OneDimensionalModel_BCFunction_Absorbing *> ( &( *absorbing ) ), _1, _2 ) );

    // BC to test A_from_P conversion
    //Constant constantArea( 1.05 );
    //BCFunction_Type constantAreaFunction( boost::bind( &Constant::operator(), &constantArea, _1 ) );

    //Constant constantPressure( 24695.0765959599 );
    //BCFunction_Type constantPressureFunction( boost::bind( &Constant::operator(), &constantPressure, _1 ) );

    OneDModel.GetBC().setBC( OneD_left,  OneD_first, OneD_Q,  sinusoidalFunction  );
    OneDModel.GetBC().setBC( OneD_right, OneD_first, OneD_W2, absorbingFunction );

    //OneDModel.GetBC().setBC( OneD_right, OneD_first, OneD_A,   constantAreaFunction );
    //OneDModel.GetBC().setBC( OneD_right, OneD_first, OneD_P,   constantPressureFunction );

    OneDModel.SetupModel();

    // *********************************
    // Tempolar loop
    // *********************************
    printf("\nTemporal loop:\n");

    Chrono chronoTotal;
    Chrono chronoSystem;
    Chrono chronoIteration;

    int count = 0;
    chronoTotal.start();
    for ( ; OneDModel.GetData().dataTime()->canAdvance() ; OneDModel.GetData().dataTime()->updateTime(), ++count )
    {
        std::cout << std::endl << "--------- Iteration " << count << " time = " << OneDModel.GetData().dataTime()->getTime() << std::endl;

        chronoIteration.start();

        if ( OneDModel.GetData().dataTime()->isFirstTimeStep() )
            OneDModel.BuildSystem();
        else
            OneDModel.UpdateSystem();

        chronoSystem.start();

        OneDModel.SolveSystem();

        chronoSystem.stop();

        //Save solution
        if ( count%50 == 0 || OneDModel.GetData().dataTime()->isLastTimeStep() )
            OneDModel.SaveSolution();

        chronoIteration.stop();

        std::cout << " System solved in " << chronoSystem.diff()
                  << " s, (total time " << chronoIteration.diff() << " s)." << std::endl;
    }

    chronoTotal.stop();
    std::cout << std::endl << " Simulation ended successfully in " << chronoTotal.diff()
              << " s" << std::endl;

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    if ( check )
    {
        bool ok = true;
        int rightnodeid = OneDModel.GetSolver()->rightNodeId();

        ok = ok && checkValue( 0.999998  , (*OneDModel.GetSolution("A"))[rightnodeid - 0]);
        ok = ok && checkValue(-0.00138076, (*OneDModel.GetSolution("Q"))[rightnodeid - 0]);
        ok = ok && checkValue(-0.00276153, (*OneDModel.GetSolution("W1"))[rightnodeid - 0]);
        ok = ok && checkValue( 0.00000000, (*OneDModel.GetSolution("W2"))[rightnodeid - 0]);

        ok = ok && checkValue( 0.999999  , (*OneDModel.GetSolution("A"))[rightnodeid - 1]);
        ok = ok && checkValue(-0.00040393, (*OneDModel.GetSolution("Q"))[rightnodeid - 1]);
        ok = ok && checkValue(-0.00080833, (*OneDModel.GetSolution("W1"))[rightnodeid - 1]);
        ok = ok && checkValue( 0.00000045, (*OneDModel.GetSolution("W2"))[rightnodeid - 1]);

        if (ok)
        {
            std::cout << "Test succesful" << std::endl;
            return 0;
        }
        else
        {
            std::cout << "Test unseccesful" << std::endl;
            return -1;
        }
    }
    else
        return 0;
}
