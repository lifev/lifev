/* -*- mode: c++ -*-
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <iostream>
#include <string>
using namespace std;

#include "GetPot.hpp"

using getpot::StringVector;

int main(int argc, char** argv)
{
    GetPot cl(argc, argv);
    int num_files=cl.vector_variable_size("--input-file");
    if (num_files == 0)
    {
        cerr<< "AT LEAST ONE INPUT FILE MUST BE GIVEN" <<endl;
        cerr<< "input-file='file 1 file 2 ...'"<<endl;
        return 1;
    }
    StringVector arg;
    for (int i =0; i != num_files; ++i) arg.push_back(string(cl("--input-file","NONE",i)));
    GetPot inputdata(arg);
    inputdata.print();
    //GetPot inputdata(cl("input-file","NONE",0));

    //for (int i=1; i<num_files; ++i) meshdata(cl("input-file",i));
    inputdata.set_prefix("mesh/");
    int check=inputdata("check",0);
    int facemarkers=inputdata("facemarkers",0);
    int edgemarkers=inputdata("edgemarkers",0);
    string logfile=string(inputdata("logfile","NONE"));
    cout << check <<" " <<facemarkers<<" " <<edgemarkers<< " " << logfile<<endl;


}

