/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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

    about.addAuthor("Tiziano Passerini", "developer", "christophe.prudhomme@epfl.ch", "");
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

    std::ofstream outfile;
    std::string file_output,str;

    // Reading from data file
    GetPot command_line(argc,argv);
    string data_file_name = command_line.follow("datanet", 2, "-f","--file");
    GetPot data_file(data_file_name);

    // Create a OneDNet variable
    OneDNet< OneDModelSolver, OneDNonLinModelParam > onednet( data_file );

  	std::cout << "[main] process " << Comm.MyPID() << " owns tube 1? " << onednet.tube_map().MyGID(1) << std::flush;
  	std::cout << "[main] process " << Comm.MyPID() << " owns tube 2? " << onednet.tube_map().MyGID(2) << std::flush;

  	if( onednet.tube_map().MyGID(1) )
    {
    	Debug( 6030 ) << "[main] process " << Comm.MyPID() << " setting BC for tube 1!\n";
//    	OneDBCFunctionPointer heart_pressure ( new Heart(  data_file,
//                                                       onednet.Param(1),
//                                                       onednet.Solver(1).Mesh(),
//                                                       onednet.Solver(1).FluxFun(), onednet.Solver(1).SourceFun(),
//                                                       onednet.Solver(1).U1_thistime(), onednet.Solver(1).U2_thistime(),
//                                                       onednet.Solver(1).W1_thistime(), onednet.Solver(1).W2_thistime(),
//                                                       onednet.Solver(1).timestep(),
//                                                       "left" /*border*/,
//                                                       "W1"  /*var*/,
//                                                      false /*type*/ ) );
//        onednet.setBC( heart_pressure, 1, "left", "first", "Q" );
        // � necessario che sia OneDBCFunctionPointer
        OneDBCFunctionPointer mysin( new Sin( 0, 1, 0.5 ) );
        //OneDBCFunctionPointer mysinneg( new Sin( 0, -1, 0.5 ) );
        onednet.setBC( mysin, 1, "left", "first", "Q" );
    }
    if( onednet.tube_map().MyGID(2) )
    {
    	Debug( 6030 ) << "[main] process " << Comm.MyPID() << " setting BC for tube 2!\n";
//    	OneDBCFunctionPointer resistence ( new Resistence( data_file("1dnetwork/tube2/parameters/R",0.),
//                                                       onednet.Param(2),
//                                                       onednet.Solver(2).Mesh(),
//                                                       onednet.Solver(2).FluxFun(), onednet.Solver(2).SourceFun(),
//                                                       onednet.Solver(2).U1_thistime(), onednet.Solver(2).U2_thistime(),
//                                                       onednet.Solver(2).W1_thistime(), onednet.Solver(2).W2_thistime(),
//                                                       onednet.Solver(2).timestep(),
//                                                       "right" /*border*/,
//                                                       "W2"  /*var*/) );
//    	onednet.setBC( resistence, 2, "right", "first", "W2" );
    }
//	onednet.setBC( heart_pressure, 3, "left", "first", "Q" );
//	onednet.setBC( heart_pressure, 5, "left", "first", "Q" );
//	onednet.setBC( resistence, 4, "right", "first", "W2" );
//	onednet.setBC( resistence, 6, "right", "first", "W2" );

    // Initialize 1D solvers in the network (need to access data_file)
    onednet.initialize( data_file );

    // Network time parameters
    Real dt     = onednet.timestep();
    Real startT = onednet.inittime();
    Real T      = onednet.endtime();
    Real postprocess_dt = data_file("1dnetwork/postprocess_timestep",0.01);
    Debug(6030) << "[main] startT T dt postprocess_dt "
                << startT << " " <<  T << " " << dt << " " << postprocess_dt << "\n";

    // Temporal loop
    printf("\nTemporal loop:\n");

    Chrono chrono;
    // Iteration counter, needed by iterate()
    int count = 0;

    for ( Real time = startT + dt ; time <= T ; time += dt ) {
        count++;
        chrono.start();
        // Update rhs for 1D solvers
        onednet.timeAdvance( time );
        // Solve 1D problems
        onednet.iterate( time , count );
        // Print results to file (but not at each time step - default every 0.01 seconds)
        if( !( static_cast<int>( std::floor( time/dt + 0.5 ) ) %
               static_cast<int>( std::floor( postprocess_dt/dt + 0.5 ) ) ) )
        {
//          std::cout << "[main] going to postprocess? " << std::flush;
        	onednet.postProcess( time );
        }

        chrono.stop();

        // Some echoes on screen
        printf("\033[0GIteration %d", count);
        printf("\033[20G, t = %f", time);
        printf("\033[40Gs... computed in %f", chrono.diff());
        printf("\033[70Gs.");

    }

    printf("\nSimulation ended successfully.\n");

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}
