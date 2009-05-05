/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  C. Malossi    <cristiano.malossi@epfl.ch>
       Date: 2009-05-05

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





/**
   \file main.hpp
   \authorC. Malossi <cristiano.malossi@epfl.ch>
   \date 2009-05-05
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

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>





// ===================================================
//! Program information
// ===================================================
LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "Test GetPot" ,
                            "LifeV Test GetPot" ,
                            "1.0",
                            "LifeV Test GetPot",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2009 EPFL");

    about.addAuthor("Cristiano Malossi", "Developer", "cristiano.malossi@epfl.ch", "");

    return about;
}





// ===================================================
//! Namespaces
// ===================================================
using namespace LifeV;





// ===================================================
//! Main
// ===================================================
int main(int argc, char** argv)
{

	#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
		std::cout << "MPI Initialization" << std::endl;
	#endif

    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow( "data", 2, "-f", "--file" );

    GetPot dataFile( data_file_name );

    int  A = dataFile( "variables/A", 0 );
    std::cout << "A: " << A << std::endl;

    Real B = dataFile( "variables/B", 0.0 );
    std::cout << "B: " << B << std::endl;

    bool C = dataFile( "variables/C", false );
    std::cout << "C: " << C << std::endl;

    std::string D = dataFile( "variables/D", "empty" );
    std::cout << "D: " << D << std::endl;

    std::string E = dataFile( "variables/E", "empty empty" );
    std::cout << "E: " << E << std::endl;


	#ifdef HAVE_MPI
		MPI_Finalize();
		std::cout << "MPI Finalization" << std::endl;
	#endif

    return( EXIT_SUCCESS );
}
