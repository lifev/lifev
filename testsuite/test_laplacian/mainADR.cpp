/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  L. Iapichino  <laura.iapichino@epfl.ch>
              C. Malossi    <cristiano.malossi@epfl.ch>
              A. Manzoni    <andrea.manzoni@epfl.ch>
       Date: 2009-03-24

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

Simple Laplacian test with Dirichlet Boundary condition

Solve the problem

               - \Delta u = f

               u = 0 on the boundary


 3D: with the source term f = 12 \pi^2 sin(2 \pi x) sin(2 \pi y) sin (2 \pi z) on a cube
 2D: with the source term f = 8 \pi^2 sin(2 \pi x) sin(2 \pi y) on a square

 the rhs is computed as rhs = Mass_Matrix * f_iterpolated


 More generally this test can solve the problem:

               - \nu \Delta u + \beta \nabla u + \sigma u = f

               u = g on the boundary

 being \nu and \sigma constants defined in the data file and \beta interpolated.

*/


/**
   \file main.hpp
   \author L. Iapichino <laura.iapichino@epfl.ch>, C. Malossi <cristiano.malossi@epfl.ch>, A. Manzoni <andrea.manzoni@epfl.ch>
   \date 2009-03-24
 */




// ===================================================
//! Includes
// ===================================================
#include <Epetra_Comm.h>

#include <boost/program_options.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include "adrProblem.hpp"





// ===================================================
//! Program information
// ===================================================
LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "Test Laplacian" ,
                            "LifeV Test Laplacian" ,
                            "1.0",
                            "Laplacian test case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2009 EPFL");

    about.addAuthor("Laura Iapichino", "Developer", "laura.iapichino@epfl.ch", "");
    about.addAuthor("Cristiano Malossi", "Developer", "cristiano.malossi@epfl.ch", "");
    about.addAuthor("Andrea Manzoni", "Developer", "andrea.manzoni@epfl.ch", "");

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





// ===================================================
//! Main
// ===================================================
int main(int argc, char** argv)
{

	#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
		std::cout << "MPI Initialization" << std::endl;
	#endif


    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");

    ADRProblem problem( argc, argv, makeAbout(), desc );
    problem.run();


	#ifdef HAVE_MPI
		MPI_Finalize();
		std::cout << "MPI Finalization" << std::endl;
	#endif

    return( EXIT_SUCCESS );
}
