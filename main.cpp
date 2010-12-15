//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing the Spirit Parser Test
 *
 *  @date 22-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  This is a test to verify that spiritParser performs correct computations.
 *  Note that the parser works only with boost v1.41+.
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <iomanip>
#include <string>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#endif

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

// LifeV includes
#include <life/lifecore/life.hpp>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/chrono.hpp>

// Mathcard includes
#include <lifemc/lifecore/Parser.hpp>

using namespace LifeV;

bool EXIT_FLAG = EXIT_SUCCESS;
std::string success = "OK     ";
std::string failed  = "FAILED ";

inline std::string
check( const bool& expression )
{
    if ( expression )
    {
        EXIT_FLAG = EXIT_FAILURE;
        return failed;
    }

    return success;
}

Int
main( Int argc, char** argv )
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
    Parser parser;

    std::cout << "READY TO TEST WITH 10+ EXPRESSIONS:" << std::endl << std::endl;

    // TEST 0:
    expression = "-sqrt(4)+1*2"; // = 0
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  0:  " << check( std::fabs(result - 0) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 1:
    expression = "(1+1)/(2+2)"; // = 0.5
    parser.setString(expression);
    result = parser.evaluate();

    std::cout << "TEST  1:  " << check( std::abs(result - 0.5) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 2:
    expression = "1-2-3+4*5-6-7+8*9+1"; // = 76
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  2:  " << check( std::abs(result - 76) > tolerance )
              << expression << " = " << result << std::endl;;

    // TEST 3:
    expression = "-(1+1)+(250+250+(-2))"; // = 496
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  3:  " << check( std::abs(result - 496) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 4:
    expression = "(0.8 > 0.9)"; // = 0
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  4:  " << check( std::abs(result - 0) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 5:
    expression = "(0.8 < 0.9)"; // = 1
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  5:  " << check( std::abs(result - 1) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 6:
    expression = "sin(3/4*pi) * -sin(3/4*pi) + -(cos(3/4*pi))^2"; // = -1
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  6:  " << check( std::abs(result - -1) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 7:
    expression = "144^0.5 * sqrt(144)"; // = 144
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  7:  " << check( std::abs(result - 144) > tolerance )
              << expression << " = " << result << std::endl;

    // TEST 8:
    expression = "abc = -2^3^-3; abc"; // = -0.001953125
    parser.setString(expression);
    result = parser.evaluate();
    std::cout << "TEST  8:  " << check( std::abs(result - -0.001953125) > tolerance )
              << expression << " = " << result << std::endl;
    // TEST 9:
    expression = "c=2; [0., c, c*c, c*c*c]"; // (0, 2, 4, 8)
    parser.setString(expression);
    std::cout << "TEST  9:  " << check( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 2 ) > tolerance || std::abs( parser.evaluate(3) - 4 ) > tolerance || std::abs( parser.evaluate(4) - 8 ) > tolerance )
              << expression << " = [" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << ", " << parser.evaluate(4) << "]" << std::endl;

    // TEST 10:
    expression = "[0, 0, -(x^2)+y^2]";
    parser.setString(expression);
    parser.setVariable("x", 1);
    parser.setVariable("y", 2); // (0, 0, -5)
    std::cout << "TEST 10a: " << check( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 0 ) > tolerance || std::abs( parser.evaluate(3) - 3 ) > tolerance )
              << "x = " << 1 << ", y = " << 2 << " ==> "
              << expression << " = [" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << "]" << std::endl;


    parser.setString(expression);
    parser.setVariable("x", 4);
    parser.setVariable("y", 5); // (0, 0, -41)
    std::cout << "TEST 10b: " << check( std::abs( parser.evaluate(1) - 0 ) > tolerance || std::abs( parser.evaluate(2) - 0 ) > tolerance || std::abs( parser.evaluate(3) - 9 ) > tolerance )
              << "x = " << 4 << ", y = " << 5 << " ==> "
              << expression << " = [" << parser.evaluate(1) << ", " << parser.evaluate(2) << ", " << parser.evaluate(3) << "]" << std::endl;

    std::cout << std::endl << "TEST ENDS SUCCESFULLY" << std::endl;

//    // PERFORMANCE TEST
//	Chrono chrono;
//
//	expression = "sqrt(((index+pi)*2)^3)"; //We test ONE expression containing different operations
//	parser.SetString(expression);
//
//	chrono.start();
//	UInt nEvaluations = 10000000; // 10 Milions
//	for (UInt i = 0 ; i < nEvaluations ; ++i)
//	{
//	    parser.setVariable("index", i);
//		parser.evaluate();
//	}
//	chrono.stop();
//
//	std::cout << std::endl << "Total time for " << nEvaluations << " evaluations of expression f=" << expression << " --> " << chrono.diff() << " s" << std::endl;

#ifdef HAVE_MPI
    std::cout << std::endl << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    return( EXIT_FLAG );
}
