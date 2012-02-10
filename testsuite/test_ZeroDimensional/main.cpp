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
 *  @brief File containing the Multiscale Test
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  This is a very general main file to run a Multiscale simulation.
 *
 *  Models available:
 *  <ol>
 *      <li> Fluid3D
 *      <li> FSI3D
 *      <li> OneDimensional
 *      <li> Multiscale
 *      <li> Windkessel0D
 *  </ol>
 *
 *  Couplings available:
 *  <ol>
 *      <li> BoundaryCondition
 *      <li> Stress
 *      <li> FlowRateStress
 *      <li> FlowRate
 *      <li> FlowRateValve
 *  </ol>
 *
 *  Algorithms available:
 *  <ol>
 *      <li> Aitken
 *      <li> Explicit
 *      <li> Newton
 *      <li> Broyden
 *  </ol>
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <sys/stat.h>
#include <sys/types.h>

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
#include <life/lifecore/LifeV.hpp>

// Mathcard includes
#include <lifemc/lifesolver/MultiscaleModel0D.hpp>

using namespace LifeV;
using namespace Multiscale;

bool checkValue(const Real val, const Real test, const Real tol = 1.e-5, const bool verbose = true)
{
  Real norm = abs(val - test);

  if ( verbose )
    std::cout << "value = " << val << " computed value = " << test << " diff = " << norm << std::endl;

  return (norm < tol);
}

Int
main( Int argc, char** argv )
{
    //Setup main communicator
    boost::shared_ptr< Epetra_Comm > comm;

    //Setup MPI variables
    Int numberOfProcesses(1);
    Int rank(0);

#ifdef HAVE_MPI
    MPI_Init( &argc, &argv );

    MPI_Comm_size( MPI_COMM_WORLD, &numberOfProcesses );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif

    if ( rank == 0 )
        std::cout << "MPI Processes: " << numberOfProcesses << std::endl;

#ifdef EPETRA_MPI
    if ( rank == 0 )
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    std::cout << "SERIAL Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_SerialComm() );
#endif


// *********************************
// Reading from data file
// *********************************
    GetPot command_line(argc,argv);

    // checking if we are checking for the nightly build
    const bool check = command_line.search(2, "-c", "--check");
    string fileName = command_line.follow("data", 2, "-f","--file");
    //-----------------------------------Test 3 : make Data Container---------------
    MultiscaleModel0D zeroDmodel;
    zeroDmodel.setCommunicator( comm );
    zeroDmodel.setupData(fileName + ".dat");
    zeroDmodel.showMe();

    // *********************************
    // Tempolar loop
    // *********************************
    std::cout << "\nTemporal loop:" << std::endl;

    LifeChrono chronoTotal;
    LifeChrono chronoSystem;
    LifeChrono chronoIteration;

    Int count = 0;
    chronoTotal.start();
    zeroDmodel.data().dataTime()->setInitialTime(0);
    zeroDmodel.setupModel();
    zeroDmodel.saveSolution();

    zeroDmodel.data().dataTime()->updateTime();
    zeroDmodel.data().dataTime()->setInitialTime(zeroDmodel.data().dataTime()->time());

   // Move to the "true" first time-step

    for ( ; zeroDmodel.data().dataTime()->canAdvance() ; zeroDmodel.data().dataTime()->updateTime(), ++count )
    {
        std::cout << std::endl << "--------- Iteration " << count << " time = " << zeroDmodel.data().dataTime()->time() << std::endl;

        chronoIteration.start();

        if ( zeroDmodel.data().dataTime()->isFirstTimeStep() )
            zeroDmodel.buildModel();
        else
            zeroDmodel.updateModel();

        chronoSystem.start();

        zeroDmodel.solveModel();

        chronoSystem.stop();

        //Save solution
        zeroDmodel.saveSolution();

        chronoIteration.stop();

        std::cout << " System solved in " << chronoSystem.diff() << " s, (total time " << chronoIteration.diff() << " s)." << std::endl;
    }

    chronoTotal.stop();
    std::cout << std::endl << " Simulation ended successfully in " << chronoTotal.diff()  << " s" << std::endl;
    //-------------------------------------------------------------------------------------

    if ( check )
      {
        bool ok = true;
    
        ok = ok && checkValue( 0.001329039627  , zeroDmodel.data().circuitData() ->Nodes() ->nodeListAt(1) ->voltage());
        ok = ok && checkValue( 0.000787475119  , zeroDmodel.data().circuitData() ->Elements() ->elementListAt(1) ->current());
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


    bool exitFlag = EXIT_SUCCESS;


#ifdef HAVE_MPI
    if ( rank == 0 )
        std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif
    std::cout << "Finished" << std::endl;

    return exitFlag;
}
