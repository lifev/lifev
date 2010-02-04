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
//! \author:Vincent Martin 09/04
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

//#include <lifemc/lifesolver/dataOneDModel.hpp>
#include <lifemc/lifesolver/oneDModelSolver.hpp>
#include <lifemc/lifefem/oneDBCFunctions.hpp>
#include "ud_functions.hpp"

#include <sstream>

using namespace LifeV;


int main(int argc, char** argv)
{



#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::flush;

    MPI_Init(&argc,&argv);

    Epetra_Comm* comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
//    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    std::cout << "ok" << std::endl;
#else
    Epetra_Comm* comm = new Epetra_SerialComm();
#endif;


  // *********************************
  // Useful typedefs
  // *********************************
  typedef NonLinearFluxFun1D                          Flux1D;
  typedef NonLinearSourceFun1D                        Source1D;
  typedef OneDNonLinModelParam                        Params1D;
  typedef OneDModelSolver<Params1D, Flux1D, Source1D> onedsolver_type;

  // *********************************
  // ***** Reading from data file
  // *********************************
  GetPot command_line(argc,argv);
  string data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  // *********************************
  // Build the 1D model
  // *********************************
  std::cout << "====== Building 1D model " << std::endl;

  std::cout << "    Building the mesh ... " << std::flush;

  const std::string section = "onedmodel";

  DataOneDModel data  (data_file, section);
  Params1D      params(data_file, section);
  std::cout << "ok." << std::endl;

  std::cout << data.mesh() << std::endl;

  const RefFE*    refFE = &feSegP1;
  const QuadRule* qR    = &quadRuleSeg3pt;
  const QuadRule* bdQr  = &quadRuleSeg1pt;


  //  boost::shared_ptr<BasicOneDMesh> mesh = &data.mesh();

  //   RegionMesh* pmesh = &data.mesh();
  //   boost::shared_ptr<RegionMesh> mesh(mesh);

  std::cout << "    Building FE Space ... " << std::flush;

  //std::cout << feSegP1.nbCoor << std::endl;
  FESpace<RegionMesh, EpetraMap> odFESpace(data.mesh(), *refFE, *qR, *bdQr, 1, *comm);
  std::cout << "nbCoorFE " << odFESpace.fe().nbCoor << std::endl;
  std::cout << "ok." << std::endl;

  //odFESpace.dof().showMe(std::cout, true);
  std::cout << "    Building Solver   ... " << std::flush;
  onedsolver_type onedm(data, params, odFESpace, *comm);


  onedm.setup(data_file, section);
  std::cout << "ok." << std::endl;

  onedm.showMe(std::cout);


  //

  //OneDBCFunctionPointer sinusoidal_flux ( new Sin() );
//   OneDBCFunctionPointer pressure( new PressureRamp<Flux1D, Source1D, OneDNonLinModelParam>
//                                   (odFESpace,
//                                    onedm.FluxFun(),
//                                    onedm.SourceFun(),
//                                    onedm.U_thistime(),
//                                    data.timestep(),
//                                    "left" /*border*/,
//                                    "W1"  /*var*/,
//                                    onedm.oneDParams()));

  OneDBCFunctionPointer resistence ( new Resi<Flux1D, Source1D, OneDNonLinModelParam>
                                    ( data_file("parameters/R",0.),
                                      onedm.oneDParams(),
                                      odFESpace,
                                      onedm.FluxFun(),
                                      onedm.SourceFun(),
                                      onedm.U_thistime(),
                                      data.timestep(),
                                      "right" /*border*/,
                                      "W2"  /*var*/, true));

  //   onedm.bcH().setBC( sinusoidal_flux, "left",  "first", "Q"  );

  OneDBCFunctionPointer sinusoidal_flux ( new Sin() );
  onedm.bcH().setBC( sinusoidal_flux, "left", "first", "Q" );

  //onedm.bcH().setBC( pressure,         "left",  "first", "W1"  );
  onedm.bcH().setBC( resistence,      "right",  "first", "W2" );

  // Initialization
  //
  Real dt             = data.timestep();
  Real startT         = data.inittime();
  Real T              = data.endtime();
  UInt postprocess_it = data_file("miscellaneous/postprocess_timestep", 10);

  Debug(6030) << "[main] startT T dt postprocess_it "
              << startT << " " <<  T << " " << dt
              << " " << postprocess_it << "\n";


  onedm.initialize(data_file, section);
//   Real u1_0 = 3.14; //! constant initial condition
//   Real u2_0 = 0.;    //! constant initial condition

//   //     if (fsi->isSolid())
//   //     {
//   //        std::cout << "     1d- initialize tube with constant A_0 = " << u1_0 << " and Q_0 = " << u2_0 << std::endl;
//   onedm.initialize(u1_0, u2_0);
//   std::cout << "    1d- initialize tube with constant A_0 = " << onedm.BCValuesLeft()[0]
//             << " and Q_0 = " << onedm.BCValuesLeft()[1] << std::endl;
//   //     }
  // Temporal loop
  printf("\nTemporal loop:\n");

  Chrono chronota;
  Chrono chronoit;

  Chrono chrono;

  int count = 0;

  for (Real time=startT + dt ; time <= T; time+=dt)
  {
      std::cout << std::endl;
      std::cout << "--------- Iteration " << count << " time = " << time << std::endl;

    count++;

    chrono.start();
    Debug(6030) << "[main] 1d model time advance\n";
    chronota.start();
    onedm.timeAdvance( time );

    chronota.stop();

    Debug(6030) << "[main] 1d model iterate\n";
    chronoit.start();
    onedm.iterate( time , count );
    chronoit.stop();

    //    std::cout <<
    if( !( static_cast<int>( std::floor( count%postprocess_it))))
//     if( !( static_cast<int>( std::floor( time/dt + 0.5 ) ) %
//            static_cast<int>( std::floor( postprocess_dt/dt  + 0.5 ) ) ) )
		{
            std::cout << "PostProcessing at time ... " << time << std::endl;
            onedm.postProcess( time );
        }

//     if ( data_file( "miscellaneous/show_graceplot", 0 ) )
//         onedm.gplot();


    chrono.stop();

    std::cout << " time adv. computed in " << chronota.diff()
              << " iter computed in " << chronoit.diff()
              << " total " << chrono.diff() << " s." << std::endl;
//     printf("\033\n");
//     printf("\033[0GIteration %d", count);
//     printf("\033[25G, t = %f", time);
//     printf("\033[25Gs, time adv. computed in %f", chronota.diff());
//     printf("\033[60Gs, iter. computed in %f", chronoit.diff());
//     printf("\033[85Gs. total %f", chrono.diff());
//     printf("\033[100Gs.\n");
//     printf("\033\n");

  }

  printf("\nSimulation ended successfully.\n");

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif


  return 0;



}
