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
#include <iostream>
#include "GetPot.hpp"
#include "fhnSolver.hpp"
#include "chrono.hpp"

/*

  Fitzhugh-Nagumo solver

  Author: J.F. Gerbeau
  
  usage: test_fhn              : read the data file "data" and run the solver
         test_fhn -f otherdata : read the data file "otherdata" and run the solver
	 test_fhn -h           : read the data file, print help and exit
	 test_fhn -i           : read the data file, print the read values and exit

*/
int main(int argc, char** argv)
{
    using namespace LifeV;
    GetPot command_line(argc,argv);
    const char* data_file_name = command_line.follow("data", 2, "-f","--file");
    GetPot data_file(data_file_name);
    if( command_line.search(2, "-i","--info") ) {
        data_file.print();
        exit(0);
    }
    Chrono chrono;
    //
    std::cout << "*** Initialisation --->" << std::endl;
    chrono.start();
    FhNSolver pb(data_file);
    pb.dataFhNShowMe();
    pb.dataAztecShowMe();
    pb.dataTransientShowMe(std::cout);
    if(command_line.search(2,"-h","--help")){
        std::cout << std::endl << std::endl;
        std::cout <<"usage: test_fhn              : read the data file 'data' \n";
        std::cout <<"       test_fhn -f otherdata : read the data file 'otherdata' \n";
        std::cout <<"       test_fhn -h           : help and exit\n";
        std::cout <<"       test_fhn -i           : read the data file, print the read values and exit\n";
        std::cout << std::endl;
        std::cout << "Help for the data file:\n";
        pb.dataAztecHelp();
        pb.dataAztecShowMe();
        pb.dataFhNHelp();
        pb.dataFhNShowMe();
        exit(0);
    }
    //
    chrono.stop();
    std::cout << "<--- Initialisation done in " << chrono.diff() << "s." << std::endl;

    if(pb.verbose) std::cout << "*** Resolution of the problem --->\n";
    chrono.start();
    pb.solve();
    chrono.stop();
    if(pb.verbose)std::cout << "<--- Solved in " << chrono.diff()
                       << "s." << std::endl << std::endl;

    std::cout << "FhN normal end." << std::endl;
    return 0;
}


