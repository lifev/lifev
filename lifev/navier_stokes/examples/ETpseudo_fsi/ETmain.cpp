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
    @file
    @brief

    @author Matteo Pozzoli <matteo1.pozzoli@mail.polimi.it>
    @author Claudia Colciago <claudia.colciago@epfl.ch>
    @date 08-03-2011
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
// #pragma GCC diagnostic warning "-Wunused-variable"
// #pragma GCC diagnostic warning "-Wunused-parameter"

#include "ETRobinMembraneSolver.hpp"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>

#include "ud_functions.hpp"


using namespace LifeV;

std::set<UInt> parseList( const std::string& list )
{
     std::string stringList = list;
     std::set<UInt> setList;
     if ( list == "" )
     {
         return setList;
     }
     size_t commaPos = 0;
     while ( commaPos != std::string::npos )
     {
         commaPos = stringList.find( "," );
         setList.insert( atoi( stringList.substr( 0, commaPos ).c_str() ) );
         stringList = stringList.substr( commaPos+1 );
     }
     setList.insert( atoi( stringList.c_str() ) );
     return setList;
}

int
main( int argc, char** argv )
{
 #ifdef HAVE_MPI
     MPI_Init(&argc, &argv);
 #endif

 //**************** cylinder
 //    MPI_Init(&argc,&argv);

     ETRobinMembraneSolver cyl( argc, argv );
     cyl.run();

 #ifdef HAVE_MPI
     MPI_Finalize();
 #endif
       return( EXIT_SUCCESS );
}
