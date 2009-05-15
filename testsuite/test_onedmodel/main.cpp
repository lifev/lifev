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
#include <life/lifesolver/dataOneDModel.hpp>
#include <life/lifesolver/oneDModelSolver.hpp>
#include <life/lifecore/GetPot.hpp>
#include "ud_functions.hpp"

#include <sstream>


int main(int argc, char** argv)
{
  using namespace LifeV;

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
  string data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  // *********************************
  // Build the 1D model
  // *********************************
  std::cout << "======\n\tBuilding 1D model " << std::endl;
  onedsolver_type onedm(data_file);
  onedm.showMe( std::cout );


  OneDBCFunctionPointer sinusoidal_flux ( new Sin() );

  OneDBCFunctionPointer resistence (
      new Resistance<Flux1D, Source1D, OneDNonLinModelParam>
      ( data_file("parameters/R",0.),
          onedm.oneDParam(),
          onedm.Mesh(),
          onedm.FluxFun(), onedm.SourceFun(),
          onedm.U_thistime(),
          onedm.timestep(),
          "right" /*border*/,
          "W2"  /*var*/) );

  onedm.bcH().setBC( sinusoidal_flux, "left", "first", "Q" );
  onedm.bcH().setBC( resistence, "right", "first", "W2" );

  // Initialization
  //
  Real dt     = onedm.timestep();
  Real startT = onedm.inittime();
  Real T      = onedm.endtime();
  Real postprocess_dt = data_file("miscellaneous/postprocess_timestep",0.01);
  Debug(6030) << "[main] startT T dt postprocess_dt "
              << startT << " " <<  T << " " << dt
              << " " << postprocess_dt << "\n";

  onedm.initialize(data_file);

  // Temporal loop
  printf("\nTemporal loop:\n");

  Chrono chrono;
  int count = 0;
  for (Real time=startT+dt ; time <= T; time+=dt) {
    count++;
    chrono.start();

    onedm.timeAdvance( time );
    onedm.iterate( time , count );

    if( !( static_cast<int>( std::floor( time/dt + 0.5 ) ) %
     static_cast<int>( std::floor( postprocess_dt/dt  + 0.5 ) ) ) )
		onedm.postProcess( time );

    if ( data_file( "miscellaneous/show_graceplot", 0 ) )
      onedm.gplot();

    chrono.stop();

    printf("\033[0GIteration %d", count);
    printf("\033[15G, t = %f", time);
    printf("\033[25Gs... computed in %f", chrono.diff());
    printf("\033[50Gs.");

  }

  printf("\nSimulation ended successfully.\n");

  return 0;
}
