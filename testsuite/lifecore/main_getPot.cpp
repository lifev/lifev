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
    @brief Test GetPot type conversion

    @author C. Malossi <cristiano.malossi@epfl.ch>
    @contributor
    @maintainer

    @date 2009-05-05

	Check if GetPot is able to read correctly an integer, real, boolean, string value
 */

// ===================================================
//! Includes
// ===================================================
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/life.hpp>



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
    const string data_file_name = command_line.follow( "data_getPot", 2, "-f", "--file" );

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
