/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  F. Nobile     <fabio.nobile@polimi.it>
              M. Pozzoli    <matteo1.pozzoli@mail.polimi.it>
              C. Vergara    <christian.vergara@polimi.it>
       Date: 2010-02-15

  Copyright (C) 2009 EPFL

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
/* ========================================================

/* ========================================================

Solve the problem

              \frac{\partial u}{\partial t} - \Delta u = f

               u = u0 on the boundary

linear_function.hpp:

  uexact = exp(-sin(Pi/2*t))*(x+y+z);

  f = (3 \pi^2 + 1 )  exp(-t)  sin( \pi x) sin(\pi y) sin ( \pi z) on a cube

nonlinear_function.hpp:

   uexact = exp(-sin(Pi/2*t))*cos(x *Pi)*cos(y*Pi)*cos(z*Pi);

   f = Pi2/4*( sin(Pi/2*t)+cos(Pi/2*t)*cos(Pi/2*t) )*exp(-sin(Pi/2*t))*cos(x *Pi)*cos(y*Pi)*cos(z*Pi);
*/
/**
  \file timeAdvance.hpp
  \author F. Nobile, M. Pozzoli, C. Vergara
  \date 2010-02-15
*/


// ===================================================
//! Includes
// ===================================================
#include "Epetra_config.h"
#ifdef HAVE_MPI
	#include "mpi.h"
	#include "Epetra_MpiComm.h"
#else
	#include "Epetra_SerialComm.h"
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include "timeAdvance.hpp"

// ===================================================
//! Program information
// ===================================================
LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "Test TimeAdvance Order I" ,
                            "LifeV Test TimeAdvance Order I" ,
                            "1.0",
                            "Time Advance test case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2010 MOX");

    about.addAuthor("Fabio Nobile", "Developer", "fabio.nobile@polimi.it", ""); 
    about.addAuthor("Matteo Pozzoli", "Developer", "matteo1.pozzoli@mail.polimi.it", "");
    about.addAuthor("Christian Vergara", "Developer", "christian.vergara@polimi.it", "");
    return about;
}





// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;

namespace
{
	static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
	static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}


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



// ===================================================
//! Main
// ===================================================
int
main( int argc, char** argv )
{   
   #ifdef HAVE_MPI
      MPI_Init(&argc, &argv);

      boost::shared_ptr<Epetra_MpiComm> Comm(new Epetra_MpiComm( MPI_COMM_WORLD ) );
      if ( Comm->MyPID() == 0 )
	cout << "% using MPI" << endl;
    #else
      boost::shared_ptr<Epetra_SerialComm> Comm( new Epetra_SerialComm() );
      cout << "% using serial Version" << endl;
   #endif


    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");

    problem ProblemOrderI( argc, argv, Comm, makeAbout(), desc );
    ProblemOrderI.run();


	#ifdef HAVE_MPI
		MPI_Finalize();
		std::cout << "MPI Finalization" << std::endl;
	#endif

    return( EXIT_SUCCESS );
}
