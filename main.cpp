/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-04-22

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
   \file main.cpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-04-22
 */





// ===================================================
//! Includes
// ===================================================
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
	#include <mpi.h>
#endif

#include <string>

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <lifemc/lifecore/SpiritParser.hpp>





// ===================================================
//! Program information
// ===================================================
LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "SpiritParser Test" ,
                            "A spirit parser tester" ,
                            "1.0",
                            "Some practical examples about how the spirit parser is working!",
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
int
main( int argc, char** argv )
{
	#ifdef HAVE_MPI
		std::cout << "MPI Initialization" << std::endl;
		MPI_Init( &argc, &argv );
	#endif

	std::string expression;
	Real result;

	std::cout << std::setprecision(15) << std::endl;

	// Initialization of the parser
	SpiritParser parser;

	// TEST 1:
	expression = "1 + 1";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	// TEST 2:
	expression = "(0.8 > 0.9)";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	// TEST 3:
	expression = "(0.8 < 0.9)";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	// TEST 4:
	expression = "x=1, y=2, 5 * (x^2 + y^2)";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	// TEST 5:
	expression = "sin(3/4*pi) * sin(3/4*pi) + cos(3/4*pi)^2";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	// TEST 6:
	expression = "144^0.5 * sqrt(144)";
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << expression << " = " << result << std::endl;

	#ifdef HAVE_MPI
		std::cout << std::endl << "MPI Finalization" << std::endl;
		MPI_Finalize();
	#endif

	return( EXIT_SUCCESS );
}
