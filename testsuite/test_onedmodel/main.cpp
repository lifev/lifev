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

#include "lifeV.hpp"
#include "chrono.hpp"
#include "dataOneDModel.hpp"
#include "oneDModelSolver.hpp"
#include "ud_functions.hpp"
#include "GetPot.hpp"
#include <sstream>


int main(int argc, char** argv)
{
  using namespace LifeV;
  //! ********** Reading from data file ******************************************

  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  OneDNonLinModelParam onedparamNL(data_file_name);
  LinearSimpleParam onedparamLin(data_file_name);

  std::cout << "======\n\tNon Linear model " << std::endl;
  onedparamNL.showMeData(std::cout);
  std::cout << "-----------------------------" << std::endl;
  std::cout << "======\n\tLinear model " << std::endl;
  onedparamLin.showMeData(std::cout);
  std::cout << "-----------------------------" << std::endl;

  OneDModelSolver onedm(data_file_name, onedparamLin);

  onedm.showMeData();
  //  onedm.showMeHandler(cout, 6);



  // Initialization
  //
  Real dt     = onedm.timestep();
  Real startT = onedm.inittime();
  Real T      = onedm.endtime();

  Real u1_0 = 0.; //! constant initial condition
  Real u2_0 = 0.; //! constant initial condition

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

  std::cout << "initialize with constant (u1_0, u2_0)" << std::endl;
  onedm.initialize(u1_0, u2_0);

  std::cout << "startT T dt " << startT << " " <<  T << " " << dt << std::endl;

  char ch;
  std::cout << "Hit return to continue" << std::endl;
  std::cin.get(ch);

  // Temporal loop
  //
  Chrono chrono;
  int count = 0;
  for (Real time=startT+dt ; time <= T; time+=dt) {
    count++;
    std::cout << "Iteration " <<  count  << ", t = " << time
	      << "s... \n";
    chrono.start();

    onedm.timeAdvance( time );
    onedm.iterate( time , count );

    chrono.stop();
    std::cout << "Iteration " <<  count  << " computed in " << chrono.diff() << " s.\n" << std::endl;

    if ( data_file( "miscellaneous/show_graceplot", 0 ) )
      onedm.gplot();

    // ************* saving result on file *****************************************
    std::ostringstream indexout;
    indexout << (time*100);
    std::string voutname;
    voutname = onedm.PostDirectory() + "/res.res" + indexout.str();
    // fstream Resfile(voutname.c_str(),ios::out | ios::binary);
    std::fstream Resfile(voutname.c_str(),std::ios::out );
    // Resfile.write((char*)&onedm.u()(1),onedm.u().size()*sizeof(Real));
    Resfile.write((char*)&onedm.U1_thistime()(1),
		  onedm.U1_thistime().size()*sizeof(Real));
    Resfile.close();

    //onedm.postProcess();
  }

  std::cout << "Hit return to close" << std::endl;
  std::cin.get(ch);

  return 0;
}
