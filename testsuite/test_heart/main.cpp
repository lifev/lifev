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
    @brief main for the test_heart

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#include <Epetra_ConfigDefs.h>
#include <Epetra_MpiComm.h>

#include <life/lifecore/life.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>
#include <heart.hpp>


namespace LifeV
{
namespace
{
static bool regIF = PRECFactory::instance().registerProduct( "Ifpack", &createIfpack );
static bool regML = PRECFactory::instance().registerProduct( "ML", &createML );
}
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

Int main( Int argc, char** argv )
{
    //! Initializing Epetra communicator
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
        cout << "% using MPI" << endl;
    Heart heart( argc, argv );
    heart.run();

    //! Finalizing Epetra communicator
    MPI_Finalize();
    return( EXIT_SUCCESS );
}
