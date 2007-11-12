/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-16

  Copyright (C) 2005 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file main.cpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-16
 */

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include "ethiersteinman.hpp"
#include "mpi.h"


LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "life_turek_cylinder" ,
                            "life_turek_cylinder" ,
                            "0.1",
                            "3D cylinder test case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2005 EPFL");

    about.addAuthor("Christophe Prud'homme", "developer", "christophe.prudhomme@epfl.ch", "");
    about.addAuthor("Christoph Winkelmann", "developer", "christoph.winkelmann@epfl.ch", "");
    return about;

}


using namespace LifeV;


int
main( int argc, char** argv )
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

//**************** cylinder
//    MPI_Init(&argc,&argv);

    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");
    Ethiersteinman es( argc, argv, makeAbout(), desc );
    es.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  return( EXIT_SUCCESS );
}


