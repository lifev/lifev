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
    @brief GetPot2Teuchos converter

    @author Umberto Villa <uvilla@emory.edu>
    @contributor
    @maintainer

    @date 04-10-2010

	This utility converts a GetPot data file into a Teuchos XML data file.

 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include<boost/shared_ptr.hpp>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/LifeV.hpp>
#include <life/lifefilters/GetPot.hpp>

#include "converter.hpp"

// Do not edit
int main(int argc, char **argv)
{
    using namespace LifeV;
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    std::cout<< "MPI Initialization\n";
#endif

#ifdef EPETRA_MPI
    boost::shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    boost::shared_ptr<Epetra_Comm> comm(new Epetra_SerialComm());
#endif

    std::string dataFileName;
    GetPot command_line(argc, argv);
    dataFileName = command_line.follow("data", 2, "-f", "--file");

    Teuchos::ParameterList pList;
    bool check;
    check = fillFromGetPot(dataFileName, pList);

    writeParameterListToXmlFile(pList, "data.xml");

    comm.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout<< "MPI Finalization \n";
#endif

    return check;
}
