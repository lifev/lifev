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

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>

#include <lifemc/lifesolver/MS_Model_1D.hpp>

#include "ud_functions.hpp"

#include <sstream>

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

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::flush;

    MPI_Init(&argc,&argv);

    boost::shared_ptr<Epetra_MpiComm> comm;
    comm.reset(new Epetra_MpiComm( MPI_COMM_WORLD ));

    std::cout << "ok" << std::endl;
#else
    boost::shared_ptr<Epetra_SerialComm> comm;
    comm.reset(new Epetra_SerialComm());
#endif

  // *********************************
  // Useful typedefs
  // *********************************
  typedef MS_Model_1D::Physics_Type                Physics_Type;
  typedef MS_Model_1D::Flux_Type                   Flux_Type;
  typedef MS_Model_1D::Source_Type                 Source_Type;
  typedef MS_Model_1D::BC_Type                     BC_Type;

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
  MS_Model_1D OneDModel;
  OneDModel.SetCommunicator(comm);
  OneDModel.SetupData( data_file_name );
  OneDModel.SetupModel();

  // *********************************
  // Boundary Conditions of the 1D model
  // *********************************
  BC_Type::OneDBCFunction_PtrType resistence ( new Resi( data_file("PhysicalParameters/R",0.),
                                                         OneDModel.GetPhysics(),
                                                         OneDModel.GetFESpace(),
                                                         OneDModel.GetFlux(),
                                                         OneDModel.GetSource(),
                                                         OneDModel.GetSolver().U_thistime(),
                                                         OneDModel.GetData().dataTime()->getTimeStep(),
                                                         "right" /*border*/,
                                                         "W2"  /*var*/, true )
                                             );

  BC_Type::OneDBCFunction_PtrType sinusoidal_flux ( new Sin() );

  OneDModel.GetBC().setBC( sinusoidal_flux, "left",  "first", "Q"  );
  OneDModel.GetBC().setBC( resistence,      "right", "first", "W2" );
  OneDModel.GetBC().setDefaultBC(OneDModel.GetFESpace(), OneDModel.GetSource(), OneDModel.GetData().dataTime()->getTimeStep());

  // *********************************
  // Initialization
  // *********************************
  Real dt             = OneDModel.GetData().dataTime()->getTimeStep();
  Real startT         = OneDModel.GetData().dataTime()->getInitialTime();
  Real T              = OneDModel.GetData().dataTime()->getEndTime();
  OneDModel.BuildSystem();

  // *********************************
  // Tempolar loop
  // *********************************
  printf("\nTemporal loop:\n");

  Chrono chronota;
  Chrono chronoit;

  Chrono chrono;

  int count = 0;

  for ( Real time=startT + dt ; time <= T; time+=dt )
  {
      std::cout << std::endl;
      std::cout << "--------- Iteration " << count << " time = " << time << std::endl;

    count++;

    chrono.start();
    Debug(6030) << "[main] 1d model time advance\n";
    chronota.start();
    OneDModel.GetData().dataTime()->updateTime();
    OneDModel.UpdateSystem();
    chronota.stop();

    Debug(6030) << "[main] 1d model iterate\n";
    chronoit.start();
    OneDModel.SolveSystem();
    chronoit.stop();

    //if ( !( static_cast<int>( std::floor( count%postprocess_it))) )
    {
        std::cout << "PostProcessing at time ... " << time << std::endl;
        OneDModel.GetSolver().postProcess( time );
    }

    chrono.stop();

    std::cout << " time adv. computed in " << chronota.diff() << " iter computed in " << chronoit.diff()
              << " total " << chrono.diff() << " s." << std::endl;
  }

  printf("\nSimulation ended successfully.\n");

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  if ( check )
  {
      bool ok = true;
      int rightnodeid = OneDModel.GetSolver().RightNodeId();


      ok = ok && checkValue( 0.999998  , OneDModel.GetSolver().U1_thistime()[rightnodeid - 0]);
      ok = ok && checkValue(-0.00138076, OneDModel.GetSolver().U2_thistime()[rightnodeid - 0]);
      ok = ok && checkValue(-0.00276153, OneDModel.GetSolver().W1_thistime()[rightnodeid - 0]);
      ok = ok && checkValue( 0.00000000, OneDModel.GetSolver().W2_thistime()[rightnodeid - 0]);

      ok = ok && checkValue( 0.999999  , OneDModel.GetSolver().U1_thistime()[rightnodeid - 1]);
      ok = ok && checkValue(-0.00040393, OneDModel.GetSolver().U2_thistime()[rightnodeid - 1]);
      ok = ok && checkValue(-0.00080833, OneDModel.GetSolver().W1_thistime()[rightnodeid - 1]);
      ok = ok && checkValue( 0.00000045, OneDModel.GetSolver().W2_thistime()[rightnodeid - 1]);

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
