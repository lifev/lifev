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
#include "darcySolver.hpp"
#include "chrono.hpp"

using namespace std;

/*

  Darcy soler using Mixed Hybrid finite element


  usage: darcy              : read the data file "data" and run the solver
         darcy -f otherdata : read the data file "otherdata" and run the solver
	 darcy -h           : read the data file, print help and exit
	 darcy -i           : read the data file, print the read values and exit

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
    cout << "*** Initialisation --->" << endl;
    chrono.start();
    DarcySolver pb(data_file);
    if(command_line.search(2,".hpp","--help")){
        cout << endl << endl;
        cout <<"usage: darcy              : read the data file 'data' \n";
        cout <<"       darcy -f otherdata : read the data file 'otherdata' \n";
        cout <<"	   darcy -h           : help and exit\n";
        cout <<"       darcy -i           : : read the data file, print the read values and exit\n";
        cout << endl;
        cout << "Help for the data file:\n";
        pb.dataAztecHelp();
        pb.dataAztecShowMe();
        pb.dataDarcyHelp();
        pb.dataDarcyShowMe();
        exit(0);
    }
    //
    chrono.stop();
    cout << "<--- Initialisation done in " << chrono.diff() << "s." << endl;

    if(pb.verbose)
        cout << "*** Compute the matrix --->" << endl;
    chrono.start();
    pb.computeHybridMatrixAndSourceRHS();
    chrono.stop();
    if(pb.verbose)
        cout << "<--- matrix computation done in "<< chrono.diff() << "s." << endl;
    pb.applyBC();
    //
    if(pb.verbose) cout << "*** Resolution of the hybrid system --->\n";
    chrono.start();
    pb.solveDarcy();
    chrono.stop();
    if(pb.verbose)cout << "<--- Linear system solved in " << chrono.diff()
                       << "s." << endl << endl;

    if(pb.verbose) cout << "*** Compute pressure and flux --->" << endl;
    chrono.start();
    pb.computePresFlux();
    chrono.stop();
    if(pb.verbose)
        cout << "<---  done in " << chrono.diff() << "s." << endl << endl;

    if(pb.verbose) cout << "*** Postproc --->" << endl;
    chrono.start();
    pb.postProcessTraceOfPressureRT0();
    pb.postProcessVelocityRT0();
    pb.postProcessPressureQ0();
    pb.postProcessPressureQ1();
    pb.postProcessVelocityQ1();
    chrono.stop();
    if(pb.verbose)
        cout << "<---  done in " << chrono.diff() << "s." << endl << endl;
    //

    return 0;
}


