/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*!
  \file main_net.cpp
  \author Tiziano Passerini <passerini@mate.polimi.it>
  \date 12/2006

  This program implements the test case presented in the OneDNet tutorial
  available online at
  http://www.lifev.org/lifev/documentation/how-to

*/

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/GetPot.hpp>

#include <life/lifesolver/oneDModelSolver.hpp>
#include <life/lifesolver/oneDNet.hpp>

#include <life/lifecore/debug.hpp>

#include "ud_functions.hpp"

#include "mpi.h"

LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "lifev_oned_net" ,
                            "lifev_oned_net" ,
                            "0.1",
                            "1D net test case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2005 EPFL");

    about.addAuthor("Tiziano Passerini", "developer",
        "tiziano@mathcs.emory.edu", "");
    return about;

}


using namespace LifeV;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
        cout << "% using MPI" << endl;
#else
    Epetra_SerialComm Comm;
    cout << "% using serial Version" << endl;
#endif

    // *********************************
    // Useful typedefs
    // *********************************
    typedef NonLinearFluxFun1D Flux1D;
    typedef NonLinearSourceFun1D Source1D;
    typedef OneDModelSolver<OneDNonLinModelParam, Flux1D, Source1D> onedsolver_type;

    // *********************************
    // ***** Reading from data file
    // *********************************
    GetPot command_line(argc,argv);
    string data_file_name = command_line.follow("datanet", 2, "-f","--file");
    GetPot data_file(data_file_name);

    // *********************************
    // Set up some useful variables
    // *********************************
    std::ofstream outfile;
    std::string file_output, str, prefix;

    // Create a OneDNet variable
    OneDNet< onedsolver_type > onednet( data_file );

    // for debugging purposes
    for( UInt i = 0; i < onednet.getNumVessels(); ++i )
      {
        std::cout << "[main] process " << Comm.MyPID() << " owns tube "
        << i+1 << "? " << onednet.tube_map().MyGID(i+1) << std::endl;
      }

    // data for inflow signal
    Real signal_rest_value; //( data_file("1dnetwork/tube1/inflow/signal_rest_value", 0.) );
    Real signal_amplitude; //( data_file("1dnetwork/tube1/inflow/signal_amplitude", 1.) );
    // when imposing a sin
    Real sin_period; //( data_file("1dnetwork/tube1/inflow/sin_period", 8.*std::atan(1.) ) );
    Real sin_phase; //( data_file("1dnetwork/tube1/inflow/sin_phase", 0. ) );

    /*
     We have a total of 6 tubes in this test case.
     3 of them are "inflow" tubes, the others are "outflow" tubes
     */
    UInt nTubes(3);
    ID inflow_tubes_ID[nTubes], outflow_tubes_ID[nTubes];
    inflow_tubes_ID[0] = 1; inflow_tubes_ID[1] = 3; inflow_tubes_ID[2] = 5;
    outflow_tubes_ID[0] = 2; outflow_tubes_ID[1] = 4; outflow_tubes_ID[2] = 6;

    for( UInt inTube = 0; inTube < nTubes; ++inTube ) {
      ID tubeID( inflow_tubes_ID[inTube] );
      prefix = "1dnetwork/tube1/inflow/";
      std::string::iterator prefix_end = prefix.end(), section_end;
      // Conversion using lexical_cast
      try {
          // store interface index in string variable
          str=boost::lexical_cast<std::string>(tubeID);
      }
      catch (boost::bad_lexical_cast &) {
          // conversion failed
          std::cout << "\n[OneDNet::OneDNet]My PID = " << Comm.MyPID()
          << ": lexical conversion failed!" << std::endl;
      }
      // update prefix
      prefix.replace( prefix_end - 3, prefix_end - 1, str );
      // point to the correct section in data file
      data_file.set_prefix( prefix.c_str() );

      // data for inflow signal
      signal_rest_value = data_file("signal_rest_value", 0.);
      signal_amplitude = data_file("signal_amplitude", 1.);
      // when imposing a sin
      sin_period = data_file("sin_period", 8.*std::atan(1.) );
      sin_phase = data_file("sin_phase", 0. );

      if( onednet.tube_map().MyGID(tubeID) )
        {
          Debug( 6030 ) << "[main] process " << Comm.MyPID()
          << " setting BC for tube " << tubeID << "!\n";

          OneDBCFunctionPointer mysin( new Sin( signal_rest_value,
              signal_amplitude, sin_period, sin_phase ) );

          onednet.setBC( mysin, tubeID, "left", "first", "Q" );
        }
    }

    // Initialize 1D solvers in the network (need to access data_file)
    onednet.initialize( data_file );

    // Network time parameters
    Real dt     = onednet.timestep();
    Real startT = onednet.inittime();
    Real T      = onednet.endtime();
    Real postprocess_dt = data_file("1dnetwork/postprocess_timestep",0.01);
    Real preloadT =  data_file("1dnetwork/preload_length",0.);

    Debug(6030) << "[main] process " << Comm.MyPID() << " startT T dt postprocess_dt "
                << startT << " " <<  T << " " << dt << " " << postprocess_dt << "\n";

    // Initialize 1D solvers in the network (need to access data_file)
    onednet.initialize( data_file );

    Chrono chrono;

    // Iteration counter, needed by iterate()
    int count = 0;

    Real tin, tend;
    for( int loop = 0; loop < 2; ++loop )
    {
    	if( !loop ){ // Preload loop
    		tin = startT+dt;
    		tend = preloadT;
    		if ( Comm.MyPID() == 0 ) printf("\nPreload loop:\n");
    	}
    	else{ // Time loop
    		tin = preloadT+dt;
    		tend = T;
    		if ( Comm.MyPID() == 0 ) printf("\nTime loop:\n");
    	}

    	for (Real time=tin ; dt*floor( time/dt + 1/2 ) <= tend; time+=dt)
    	{
        count++;
        chrono.start();
        // Update rhs for 1D solvers
        onednet.timeAdvance( time );
        // Solve 1D problems
        onednet.iterate( time , count );
        // Print results to file (but not at each time step - default every 0.01 seconds)
        if( !( static_cast<int>( std::floor( time/dt + 0.5 ) ) %
               static_cast<int>( std::floor( postprocess_dt/dt  + 0.5 ) ) ) )
            onednet.postProcess( time );

        chrono.stop();

        // Some echoes on screen
        printf("\033[0GIteration %d", count);
        printf("\033[15G, t = %f", time);
        printf("\033[25Gs... computed in %f", chrono.diff());
        printf("\033[50Gs.");

    }

	}
	printf("\nSimulation ended successfully.\n");

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}
