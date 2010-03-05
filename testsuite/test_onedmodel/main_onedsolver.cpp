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
#include <lifemc/lifesolver/oneDBCHandler.hpp>

#include <life/lifealg/SolverAmesos.hpp>

#include <lifemc/lifesolver/MS_Model_1D.hpp>

#include "ud_functions.hpp"

#include <sstream>

using namespace LifeV;


bool checkValue(const double val, const double test, const double tol = 1.e-5, const bool verbose = true)
{
    double norm = abs(val - test);
    std::cout << "value = " << val << " computed value = " << test << " diff = " << norm << std::endl;

    return (norm < tol);
}


int main(int argc, char** argv)
{


#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::flush;

    MPI_Init(&argc,&argv);

    boost::shared_ptr<Epetra_MpiComm> comm;
    comm.reset(new Epetra_MpiComm( MPI_COMM_WORLD ));
    int ntasks;
//    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    std::cout << "ok" << std::endl;
#else
    boost::shared_ptr<Epetra_SerialComm> comm;
    comm.reset(new Epetra_SerialComm());
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

    // checking if we are checking for the nightly build
    const bool check = command_line.search(2, "-c", "--check");

  string data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  // *********************************
  // Build the 1D model
  // *********************************
  std::cout << "====== Building 1D model " << std::endl;

  const std::string section = "onedmodel";

  std::cout << "    1d- Building the data ... " << std::flush;
  DataOneDModel data  (data_file, section);
  std::cout << "ok" << std::endl;

  std::cout << "    1d- Building the params ... " << std::flush;
  Params1D      params(data_file, section);
  std::cout << "ok" << std::endl;

  std::cout << "    1d- Building the flux function ... " << std::flush;
  Flux1D        flux(params);
  std::cout << "ok" << std::endl;


  std::cout << "    1d- Building the source function ... " << std::flush;
  Source1D      source(params);
  std::cout << "ok" << std::endl;


  MS_Model_1D od;
  od.SetCommunicator(comm);
  od.SetupData(data_file, section);
  od.SetupModel();

  //

  OneDBCFunctionPointer resistence ( new Resi<Flux1D, Source1D, OneDNonLinModelParam>
                                    ( data_file("parameters/R",0.),
                                      params,
                                      od.GetFESpace(),
                                      flux,
                                      source,
                                      od.GetSolver().U_thistime(),
                                      data.timestep(),
                                      "right" /*border*/,
                                      "W2"  /*var*/, true));

  OneDBCFunctionPointer sinusoidal_flux ( new Sin() );

  od.GetBC().setBC( sinusoidal_flux, "left",  "first", "Q"  );
  od.GetBC().setBC( resistence,      "right", "first", "W2" );
  od.GetBC().setDefaultBC(od.GetFESpace(), source, data.timestep());

  // Initialization
  //
  Real dt             = data.timestep();
  Real startT         = data.inittime();
  Real T              = data.endtime();
  UInt postprocess_it = data_file("miscellaneous/postprocess_timestep", 10);

  Debug(6030) << "[main] startT T dt postprocess_it "
              << startT << " " <<  T << " " << dt
              << " " << postprocess_it << "\n";


  //onedm.initialize();
  od.BuildSystem();

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
    od.UpdateSystem();
    chronota.stop();



    Debug(6030) << "[main] 1d model iterate\n";
    chronoit.start();
    od.SolveSystem();
    chronoit.stop();

    int leftnodeid = od.GetSolver().LeftNodeId();

    if( !( static_cast<int>( std::floor( count%postprocess_it))))
		{
            std::cout << "PostProcessing at time ... " << time << std::endl;
            //onedm.postProcess( time );
        }


    chrono.stop();

    std::cout << " time adv. computed in " << chronota.diff()
              << " iter computed in " << chronoit.diff()
              << " total " << chrono.diff() << " s." << std::endl;
  }


  printf("\nSimulation ended successfully.\n");

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif


  if (check)
  {
      bool ok = true;
      int rightnodeid = od.GetSolver().RightNodeId();


      ok = ok && checkValue( 0.999998  , od.GetSolver().U1_thistime()[rightnodeid - 0]);
      ok = ok && checkValue(-0.00138076, od.GetSolver().U2_thistime()[rightnodeid - 0]);
      ok = ok && checkValue(-0.00276153, od.GetSolver().W1_thistime()[rightnodeid - 0]);
      ok = ok && checkValue( 0.00000000, od.GetSolver().W2_thistime()[rightnodeid - 0]);

      ok = ok && checkValue( 0.999999  , od.GetSolver().U1_thistime()[rightnodeid - 1]);
      ok = ok && checkValue(-0.00040393, od.GetSolver().U2_thistime()[rightnodeid - 1]);
      ok = ok && checkValue(-0.00080833, od.GetSolver().W1_thistime()[rightnodeid - 1]);
      ok = ok && checkValue( 0.00000045, od.GetSolver().W2_thistime()[rightnodeid - 1]);

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
//       bool = bool &&  onedm.U2_thistime()[rightnodeid - 0] << std::endl;
//       std::cout << onedm.W1_thistime()[rightnodeid - 0] << std::endl;
//       std::cout << onedm.W2_thistime()[rightnodeid - 0] << std::endl;

//       std::cout << onedm.U1_thistime()[rightnodeid - 1] << std::endl;
//       std::cout << onedm.U2_thistime()[rightnodeid - 1] << std::endl;
//       std::cout << onedm.W1_thistime()[rightnodeid - 1] << std::endl;
//       std::cout << onedm.W2_thistime()[rightnodeid - 1] << std::endl;

  }
  else
      return 0;



}
