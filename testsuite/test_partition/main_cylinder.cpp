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
   \file main_cylinder.cpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
   \author Radu Popescu <radu.popescu@epfl.ch>
   \date 2010-07-02
 */

#ifdef TWODIM
#error test_cylinder cannot be compiled in 2D
#endif

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include "cylinder.hpp"
#include "EqualSolutions.hpp"

#include <mpi.h>


LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "life_cylinder" ,
                            "life_cylinder" ,
                            "0.1",
                            "3D cylinder test case w/ prepartitioned mesh in HDF5",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2010 EPFL");
    about.addAuthor("Gilles Fourestey", "developer", "gilles.fourestey@epfl.ch", "");
    about.addAuthor("Radu Popescu", "developer", "radu.popescu@epfl.ch", "");
    return about;
}

using namespace LifeV;

std::set<UInt> parseList( const std::string& list )
{
    std::string stringList = list;
    std::set<UInt> setList;
    if ( list == "" )
    {
        return setList;
    }
    UInt commaPos = 0;
    while ( commaPos != std::string::npos )
    {
        commaPos = stringList.find( "," );
        setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
        stringList = stringList.substr( commaPos+1 );
    }
    setList.insert( atoi( stringList.c_str() ) );
    return setList;
}

namespace LifeV
{
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}
}

int
main( int argc, char** argv )
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int rank;
  int return_value = EXIT_SUCCESS;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//**************** cylinder
//    MPI_Init(&argc,&argv);

    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");

    Cylinder cyl( argc, argv, makeAbout(), desc );
    cyl.run();

// Test validity of solution

    if (rank == 0) {
        if (equalSolutions("cylinder_ref.h5", "cylinder.h5", 2, 1e-6)) {
            std::cout << "TEST PARTITION WAS SUCCESSFUL.\n";
        } else {
            std::cout << "TEST PARTITION FAILED.\n";
            return_value = EXIT_FAILURE;
        }
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

  return return_value;
}

