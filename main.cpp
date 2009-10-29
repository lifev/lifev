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

#include <iomanip>
#include <string>

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/chrono.hpp>

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
	Real tolerance = 1e-15;

	std::cout << std::setprecision(30) << std::endl;

	// Initialization of the parser
	SpiritParser parser;

	std::cout << "READY TO TEST WITH 10+ EXPRESSIONS:" << std::endl << std::endl;

	// TEST 0:
	expression = "-sqrt(4)+1*2"; // = 0
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  0:  " << expression << " = " << result << std::endl;
	if ( std::fabs(result - 0) > tolerance )
		return( EXIT_FAILURE );


	// TEST 1:
	expression = "showme; (1+1)/(2+2)"; // = 0.5
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  1:  " << expression << " = " << result << std::endl;
	if ( std::abs(result - 0.5) > tolerance )
		return( EXIT_FAILURE );



	// TEST 2:
	expression = "1-2-3+4*5-6-7+8*9+1"; // = 76
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  2:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 76) > tolerance)
		return( EXIT_FAILURE );



	// TEST 3:
	expression = "-(1+1)+(250+250+(-2))"; // = 496
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  3:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 496) > tolerance)
		return( EXIT_FAILURE );



	// TEST 4:
	expression = "(0.8 > 0.9)"; // = 0
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  4:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 0) > tolerance)
		return( EXIT_FAILURE );



	// TEST 5:
	expression = "(0.8 < 0.9)"; // = 1
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  5:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 1) > tolerance)
		return( EXIT_FAILURE );



	// TEST 6:
	expression = "sin(3/4*pi) * sin(3/4*pi) + cos(3/4*pi)^2"; // = 1
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  6:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 1) > tolerance)
		return( EXIT_FAILURE );



	// TEST 7:
	expression = "144^0.5 * sqrt(144)"; // = 144
	parser.setString(expression);
	result = parser.evaluate();
	std::cout << "TEST  7:  " << expression << " = " << result << std::endl;
	if (std::abs(result - 144) > tolerance)
		return( EXIT_FAILURE );



	// TEST 8:
	expression = "c=2; (0., c, c*c, c*c*c)"; // (0, 2, 4, 8)
	parser.setString(expression);
	std::cout << "TEST  8:  " << expression << " = (" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << ", " << parser.evaluate(4) << ")" << std::endl;
	if ( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 2 ) > tolerance || std::abs( parser.evaluate(3) - 4 ) > tolerance || std::abs( parser.evaluate(4) - 8 ) > tolerance )
		return( EXIT_FAILURE );



	// TEST 9:
	expression = "(0, 0, -x^2+y^2)))";
	parser.setString(expression);
	parser.setVariable("x", 1);
	parser.setVariable("y", 2); // (0, 0, -5)
	std::cout << "TEST  9a: " << "x = " << 1 << ", y = " << 2 << " ==> ";
	std::cout << expression << " = (" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << ")" << std::endl;
	if ( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 0 ) > tolerance || std::abs( parser.evaluate(3) - 3 ) > tolerance )
		return( EXIT_FAILURE );



	parser.setString(expression);
	parser.setVariable("x", 4);
	parser.setVariable("y", 5); // (0, 0, -41)
	std::cout << "TEST  9b: " << "x = " << 4 << ", y = " << 5 << " ==> ";
	std::cout << expression << " = (" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << ")" << std::endl;
	if ( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 0 ) > tolerance || std::abs( parser.evaluate(3) - 9 ) > tolerance )
		return( EXIT_FAILURE );

	std::cout << std::endl << "TEST ENDS SUCCESFULLY -> NOW TESTING PERFORMANCES" << std::endl;

	// PERFORMANCE TEST
	Chrono chrono;

	expression = "sqrt(((1+pi)*2)^3)"; //We test ONE expression containing different operations
	parser.setString(expression);

	chrono.start();
	UInt nEvaluations = 10000000;
	for (UInt i = 0 ; i < nEvaluations ; ++i)
		parser.evaluate();
	chrono.stop();

	std::cout << std::endl << "Total time for " << nEvaluations << " evaluations of expression f=" << expression << " --> " << chrono.diff() << " s" << std::endl;

	#ifdef HAVE_MPI
		std::cout << std::endl << "MPI Finalization" << std::endl;
		MPI_Finalize();
	#endif

	return( EXIT_SUCCESS );
}
