/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s):  <simone.deparis@epfl.ch>
       Date: 2008-08-08

  Copyright (C) 2008 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file ensightToHdf5.cpp
   \author  <simone.deparis@epfl.ch>
   \date 2008-08-08
 */

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
//#include "life/lifesolver/NavierStokesSolver.hpp"
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>

#include <life/lifesolver/Oseen.hpp>

#include <iostream>

#include "ensightToHdf5.hpp"


using namespace LifeV;


struct EnsightToHdf5::Private
{
    Private() {}

    std::string    data_file_name;

    Epetra_Comm*   comm;
};




EnsightToHdf5::EnsightToHdf5( int argc,
                              char** argv,
                              LifeV::AboutData const& /*ad*/,
                              LifeV::po::options_description const& /*od*/ )
    :
    d( new Private )
{
    GetPot command_line(argc, argv);
    const char* data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    d->data_file_name = data_file_name;

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    //    MPI_Init(&argc,&argv);

    d->comm = new Epetra_MpiComm( MPI_COMM_WORLD );
    int ntasks;
    int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm = new Epetra_SerialComm();
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}

void
EnsightToHdf5::run()
{
    typedef Oseen< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;
    // Reading from data file
    //
    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);


    // fluid solver



    DataNavierStokes<RegionMesh3D<LinearTetra> > dataNavierStokes( dataFile );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(*dataNavierStokes.mesh(), *d->comm);

    std::string uOrder =  dataFile( "fluid/discretization/vel_order", "P1");

    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_press;
    const QuadRule* qR_press;
    const QuadRule* bdQr_press;

    if ( uOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 velocity " << std::flush;
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt; // DoE 5
        bdQr_vel  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( uOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 velocity ";
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;  // DoE 2
            bdQr_vel  = &quadRuleTria3pt;   // DoE 2
        }
        else
            if ( uOrder.compare("P1Bubble") == 0 )
            {
                if (verbose) std::cout << "P1-bubble velocity " << std::flush;
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;  // DoE 2
                bdQr_vel  = &quadRuleTria3pt;   // DoE 2
            }


    std::string pOrder =  dataFile( "fluid/discretization/press_order", "P1");
    if ( pOrder.compare("P2") == 0 )
    {
        if (verbose) std::cout << "P2 pressure " << std::flush;
        refFE_press = &feTetraP2;
        qR_press    = &quadRuleTetra15pt; // DoE 5
        bdQr_press  = &quadRuleTria3pt;   // DoE 2
    }
    else
        if ( pOrder.compare("P1") == 0 )
        {
            if (verbose) std::cout << "P1 pressure";
            refFE_press = &feTetraP1;
            qR_press    = &quadRuleTetra4pt;  // DoE 2
            bdQr_press  = &quadRuleTria3pt;   // DoE 2
        }

    dataNavierStokes.setMesh(meshPart.mesh());

    if (verbose)
        std::cout << "Building the velocity FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,
                                                             *refFE_vel,
                                                             *qR_vel,
                                                             *bdQr_vel,
                                                             3,
                                                             *d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;

    if (verbose)
        std::cout << "Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,
                                                             *refFE_press,
                                                             *qR_press,
                                                             *bdQr_press,
                                                             1,
                                                             *d->comm);

    if (verbose)
        std::cout << "ok." << std::endl;



    UInt totalVelDof   = uFESpace.map().getMap(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpace.map().getMap(Unique)->NumGlobalElements();


    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    Oseen< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              *d->comm);
    EpetraMap fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);
    fluid.buildSystem();

    // Initialization

    Real dt     = dataNavierStokes.timestep();
    Real t0     = dataNavierStokes.inittime();
    Real tFinal = dataNavierStokes.endtime ();

    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;
    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > importer;


    // hdf5 exporter, still under development
    exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );

    vector_ptrtype velAndPressureUnique ( new vector_type(fluid.solution(), Unique ) );

    importer.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, meshPart.mesh(), "ethiersteinman", d->comm->MyPID()) );
    vector_ptrtype velAndPressureRepeated ( new vector_type(fluid.solution(), Repeated ) );


    importer->addVariable( ExporterData::Vector, "velocity", velAndPressureRepeated,
                           UInt(0), uFESpace.dof().numTotalDof() );

    importer->addVariable( ExporterData::Scalar, "pressure", velAndPressureRepeated,
                           UInt(3*uFESpace.dof().numTotalDof()),
                           UInt(pFESpace.dof().numTotalDof()) );
    importer->import( t0 );

    *velAndPressureUnique = *velAndPressureRepeated;

    exporter->addVariable( ExporterData::Vector, "velocity", velAndPressureUnique,
                           UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "pressure", velAndPressureUnique,
                           UInt(3*uFESpace.dof().numTotalDof()),
                           UInt(pFESpace.dof().numTotalDof()) );
    exporter->postProcess( t0 );


    // Temporal loop

    Chrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {

        chrono.stop();
        importer->import( time );

        *velAndPressureUnique = *velAndPressureRepeated;

        exporter->postProcess( time );

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }

}


//////////////////////



