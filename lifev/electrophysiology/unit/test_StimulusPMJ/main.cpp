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

    @date
    @author
    @contributor
    @mantainer
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
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <fstream>
#include <string>


#include <lifev/electrophysiology/solver/ElectroETAMonodomainSolver.hpp>


#include <lifev/electrophysiology/stimulus/StimulusPMJ.hpp>
#include <lifev/core/LifeV.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"
#include <Teuchos_ScalarTraitsDecl.hpp>

#include <sys/stat.h>
#include <time.h>       /* time */

using namespace LifeV;

Int main ( Int argc, char** argv )
{

//    bool verbose = false;

    //! Initializing Epetra communicator
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm>  Comm ( new Epetra_MpiComm (MPI_COMM_WORLD) );

    //*********************************************//
    // creating output folder
    //*********************************************//
    GetPot command_line (argc, argv);


    const string stimulus_datafile_name = command_line.follow ("StimulationParameters.xml", 2, "-s", "--stimulus");
    Teuchos::ParameterList stimulusList = * ( Teuchos::getParametersFromXmlFile ( stimulus_datafile_name ) );

    StimulusPMJ stimulus;
     stimulus.setParameters ( stimulusList );
    if ( Comm->MyPID() == 0 )
    {
        stimulus.showMe();
    }

    Real DeltaX,DeltaY,DeltaZ,DeltaT;
    DeltaX = 1000;
    DeltaY = 1000;
    DeltaZ = 1000;
    DeltaT = 1000;

   if ( Comm->MyPID() == 0 )
    {
        std::cerr<<"Starting test"<<std::endl;
    }

 Teuchos::ScalarTraits<Real>::seedrandom (time (NULL) );
    for(ID i=0;i<100;i++)
    {
	Real x,y,z,t;
	x = (Teuchos::ScalarTraits<Real>::random() + stimulusList.get ( "pacing_site_X", 1.0 ) )*DeltaX;
	y = (Teuchos::ScalarTraits<Real>::random() + stimulusList.get ( "pacing_site_Y", 1.0 ) )*DeltaY;
	z = (Teuchos::ScalarTraits<Real>::random() + stimulusList.get ( "pacing_site_Z", 1.0 ) )*DeltaZ;
	t = (Teuchos::ScalarTraits<Real>::random() + 1.0 )*DeltaX;
        stimulus.appliedCurrent (  t,  x,  y,  z,  i );
    }
   if ( Comm->MyPID() == 0 )
    {
        std::cerr<<"Ending test"<<std::endl;
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize();

    return ( EXIT_SUCCESS );
}
