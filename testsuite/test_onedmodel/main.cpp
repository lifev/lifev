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
  
  // ********** Reading from data file ******************************************
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  OneDNonLinModelParam onedparamNL(data_file);

  //LinearSimpleParam onedparamLin(data_file);

  std::cout << "======\n\tNon Linear model " << std::endl;
  onedparamNL.showMeData(std::cout);
  std::cout << "-----------------------------" << std::endl;
//  std::cout << "======\n\tLinear model " << std::endl;
//  onedparamLin.showMeData(std::cout);
//  std::cout << "-----------------------------" << std::endl;

    //OneDModelSolver onedm(data_file, onedparamLin);
  OneDModelSolver onedm(data_file, onedparamNL );

  OneDBCFunctionPointer sinusoidal_flux ( new Sin() );						

  OneDBCFunctionPointer resistence ( new Resistence( data_file("parameters/R",0.),
								onedparamNL,
								onedm.Mesh(), 
  								onedm.FluxFun(), onedm.SourceFun(),	
								onedm.U1_thistime(), onedm.U2_thistime(),  
								onedm.W1_thistime(), onedm.W2_thistime(),  
								onedm.timestep(),
								"right" /*border*/,
								"W2"  /*var*/) );						
  
  //onedm.bcH().setBC( sinusoidal_flux, "left", "first", "Q" );
  //onedm.bcH().setBC( resistence, "right", "first", "W2" );
	
  onedm.showMeData();
  //  onedm.showMeHandler(cout, 6);

  // Initialization
  //
  Real dt     = onedm.timestep();
  Real startT = onedm.inittime();
  Real T      = onedm.endtime();
  Real postprocess_dt = data_file("miscellaneous/postprocess_timestep",0.01);
  Debug(6030) << "[main] startT T dt postprocess_dt "
  			<< startT << " " <<  T << " " << dt << " " << postprocess_dt << "\n";

  Real u1_0 = onedparamNL.Area0(0);
  Debug(6030)<< "[main] initializing tube with area " << u1_0 <<"\n";
  Real u2_0 = 0.; //! constant initial condition
  Debug(6030)<< "[main] initializing tube with flux " << u2_0 <<"\n";

  //onedm.initialize(u1_0, u2_0);
  onedm.initialize(data_file);

  /*
    if(startT > 0.0){
    cout << "initialize velocity and pressure with data from file" << std::endl;
    std::ostringstream indexin;
    std::string vinname, cinname;
    indexin << (startT*100);
    vinname = "fluid.res"+indexin.str();
    onedm.initialize(vinname);}
    else{
    std::cout << "initialize velocity and pressure with u0 and p0" << std::endl;
    onedm.initialize(u0,p0,0.0,dt);
    }
  */

//  char ch;
//  std::cout << "Hit return to continue" << std::endl;
//  std::cin.get(ch);

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
   	printf("\033[20G, t = %f", time);
   	printf("\033[40Gs... computed in %f", chrono.diff());
   	printf("\033[70Gs.");

  }

  printf("\nSimulation ended successfully.\n");

//  std::cout << "Hit return to close" << std::endl;
//  std::cin.get(ch);

  return 0;
}
