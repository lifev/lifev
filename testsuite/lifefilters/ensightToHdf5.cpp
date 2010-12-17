//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief Convert from Ensight to HDF5
 *
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @date 08-08-2008
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>
#include <life/lifesolver/Oseen.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/quadRule.hpp>

#include <iostream>

#include "ensightToHdf5.hpp"

using namespace LifeV;


typedef Oseen< RegionMesh3D<LinearTetra> >::vector_type  vector_type;
typedef boost::shared_ptr<vector_type> vector_ptrtype;


template <class Mesh, class Map>
void
computeP0pressure(FESpace< Mesh, Map >& pFESpace,
                  FESpace< Mesh, Map >& p0FESpace,
                  const FESpace<Mesh, Map >& uFESpace,
                  const vector_type& velAndPressureExport,
                  vector_type& P0pres, Real /*time*/);

struct EnsightToHdf5::Private
{
    Private() {}

    std::string    data_file_name;

    boost::shared_ptr<Epetra_Comm>   comm;
};

EnsightToHdf5::EnsightToHdf5( int argc,
                              char** argv )
        :
        d( new Private )
{
    GetPot command_line(argc, argv);
    string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    d->data_file_name = data_file_name;

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    // MPI_Init(&argc,&argv);
    d->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
    int ntasks;
    // int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset( new Epetra_SerialComm() );
#endif

    if (!d->comm->MyPID())
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << ntasks << " running." << std::endl;
    }
}

void
EnsightToHdf5::run()
{

    // Reading from data file
    GetPot dataFile( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);

    // Fluid solver
    boost::shared_ptr<DataNavierStokes> dataNavierStokes(new DataNavierStokes());
    dataNavierStokes->setup( dataFile );

    DataMesh dataMesh;
    dataMesh.setup(dataFile, "fluid/space_discretization");

    boost::shared_ptr<RegionMesh3D<LinearTetra> > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    readMesh(*fullMeshPtr, dataMesh);

    // writeMesh("test.mesh", *fullMeshPtr);
    // Scale, Rotate, Translate (if necessary)
    boost::array< Real, NDIM >    geometryScale;
    boost::array< Real, NDIM >    geometryRotate;
    boost::array< Real, NDIM >    geometryTranslate;

    geometryScale[0] = dataFile( "fluid/space_discretization/transform", 1., 0);
    geometryScale[1] = dataFile( "fluid/space_discretization/transform", 1., 1);
    geometryScale[2] = dataFile( "fluid/space_discretization/transform", 1., 2);

    geometryRotate[0] = dataFile( "fluid/space_discretization/transform", 0., 3) * Pi / 180;
    geometryRotate[1] = dataFile( "fluid/space_discretization/transform", 0., 4) * Pi / 180;
    geometryRotate[2] = dataFile( "fluid/space_discretization/transform", 0., 5) * Pi / 180;

    geometryTranslate[0] = dataFile( "fluid/space_discretization/transform", 0., 6);
    geometryTranslate[1] = dataFile( "fluid/space_discretization/transform", 0., 7);
    geometryTranslate[2] = dataFile( "fluid/space_discretization/transform", 0., 8);

    fullMeshPtr->transformMesh( geometryScale, geometryRotate, geometryTranslate );

    partitionMesh< RegionMesh3D<LinearTetra> >   meshPart(fullMeshPtr, d->comm);

    std::string uOrder =  dataFile( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile( "fluid/space_discretization/press_order", "P1");

    //dataNavierStokes->dataMesh()->setMesh(meshPart.mesh());

    if (verbose) std::cout << "Building the velocity FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > uFESpace(meshPart,uOrder,3,d->comm);

    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > pFESpace(meshPart,pOrder,1,d->comm);

    if (verbose) std::cout << "ok." << std::endl;

    if (verbose) std::cout << "Building the P0 pressure FE space ... " << std::flush;

    FESpace< RegionMesh3D<LinearTetra>, EpetraMap > p0FESpace(meshPart, feTetraP0, quadRuleTetra1pt,
                                                              quadRuleTria1pt, 1,d->comm);

    if (verbose) std::cout << "ok." << std::endl;



    UInt totalVelDof   = uFESpace.map().map(Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpace.map().map(Unique)->NumGlobalElements();
    UInt totalP0PresDof = p0FESpace.map().map(Unique)->NumGlobalElements();

    if (verbose) std::cout << "Total Velocity Dof = " << totalVelDof << std::endl;
    if (verbose) std::cout << "Total Pressure Dof = " << totalPressDof << std::endl;
    if (verbose) std::cout << "Total P0Press  Dof = " << totalP0PresDof << std::endl;

    if (verbose) std::cout << "Calling the fluid constructor ... ";

    Oseen< RegionMesh3D<LinearTetra> > fluid (dataNavierStokes,
                                              uFESpace,
                                              pFESpace,
                                              d->comm);
    EpetraMap fullMap(fluid.getMap());

    if (verbose) std::cout << "ok." << std::endl;

    fluid.setUp(dataFile);

    // Initialization
    Real dt     = dataNavierStokes->dataTime()->getTimeStep();
    Real t0     = dataNavierStokes->dataTime()->getInitialTime();
    Real tFinal = dataNavierStokes->dataTime()->getEndTime();

    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > exporter;
    boost::shared_ptr< Exporter<RegionMesh3D<LinearTetra> > > importer;

    std::string const exporterType =  dataFile( "exporter/type", "hdf5");
    std::string const exporterName =  dataFile( "exporter/filename", "ethiersteinman");
    std::string const exportDir    =  dataFile( "exporter/exportDir", "./");

    std::string const importerType =  dataFile( "importer/type", "ensight");
    std::string const importerName =  dataFile( "importer/filename", "ethiersteinman");
    std::string const importDir    =  dataFile( "importer/importDir", "importDir");

#ifdef HAVE_HDF5
    if (exporterType.compare("hdf5") == 0)
        exporter.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, exporterName ) );
    else
#endif
        exporter.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, exporterName ) );

    exporter->setPostDir( exportDir ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId( meshPart.meshPartition(), d->comm->MyPID() );

#ifdef HAVE_HDF5
    if (importerType.compare("hdf5") == 0)
        importer.reset( new Hdf5exporter<RegionMesh3D<LinearTetra> > ( dataFile, importerName ) );
    else
#endif
        importer.reset( new Ensight<RegionMesh3D<LinearTetra> > ( dataFile, importerName ) );

    // todo this will not work with the Ensight filter (it uses M_importDir, a private variable)
    importer->setPostDir( importDir ); // This is a test to see if M_post_dir is working
    importer->setMeshProcId( meshPart.meshPartition(), d->comm->MyPID() );

    vector_ptrtype velAndPressureExport ( new vector_type(*fluid.solution(), exporter->mapType() ) );
    vector_ptrtype velAndPressureImport ( new vector_type(*fluid.solution(), importer->mapType() ) );

    if ( exporter->mapType() == Unique )
        velAndPressureExport->setCombineMode(Zero);

    importer->addVariable( ExporterData::Vector, "velocity", velAndPressureImport,
                           UInt(0), uFESpace.dof().numTotalDof() );
    importer->addVariable( ExporterData::Scalar, "pressure", velAndPressureImport,
                           UInt(3*uFESpace.dof().numTotalDof() ),
                           UInt(  pFESpace.dof().numTotalDof() ) );
    importer->import( t0 );

    *velAndPressureExport = *velAndPressureImport;


    vector_ptrtype P0pres ( new vector_type(p0FESpace.map()) );
    MPI_Barrier(MPI_COMM_WORLD);
    computeP0pressure(pFESpace, p0FESpace, uFESpace, *velAndPressureImport, *P0pres, t0);

    exporter->addVariable( ExporterData::Vector, "velocity", velAndPressureExport,
                           UInt(0), uFESpace.dof().numTotalDof() );

    exporter->addVariable( ExporterData::Scalar, "pressure", velAndPressureExport,
                           UInt(3*uFESpace.dof().numTotalDof() ),
                           UInt(  pFESpace.dof().numTotalDof() ) );

    exporter->addVariable(ExporterData::Scalar, "P0pressure", P0pres,
                          UInt(0),
                          UInt(p0FESpace.dof().numTotalDof()),
                          UInt(0), ExporterData::Cell );
    exporter->postProcess( t0 );

    // Temporal loop
    Chrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt/2.; time += dt, iter++)
    {
        std::cout << "Doing " << time << std::endl;
        chrono.stop();
        importer->import( time );

        *velAndPressureExport = *velAndPressureImport;
        MPI_Barrier(MPI_COMM_WORLD);
        computeP0pressure(pFESpace, p0FESpace, uFESpace, *velAndPressureImport, *P0pres, time);

        exporter->postProcess( time );

        chrono.stop();
        if (verbose) std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
    }
}


template<class Mesh, class Map>
void
computeP0pressure(FESpace< Mesh, Map >& pFESpace,
                  FESpace< Mesh, Map >& p0FESpace,
                  const FESpace< Mesh, Map >& uFESpace,
                  const vector_type& velAndPressure,
                  vector_type& P0pres, Real time)
{

    int MyPID;
    MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);
    UInt offset = 3*uFESpace.dof().numTotalDof();

    std::vector<int> gid0Vec(0);
    gid0Vec.reserve(p0FESpace.mesh()->numVolumes());
    std::vector<Real> val0Vec(0);
    val0Vec.reserve(p0FESpace.mesh()->numVolumes());

    for (UInt ivol=1; ivol<= pFESpace.mesh()->numVolumes(); ++ivol)
    {

        pFESpace.fe().update( pFESpace.mesh()->volumeList( ivol ), UPDATE_DPHI );
        p0FESpace.fe().update( p0FESpace.mesh()->volumeList( ivol ) );

        UInt eleID = pFESpace.fe().currentLocalId();

        double tmpsum=0.;
        for (UInt iNode=0; iNode < (UInt) pFESpace.fe().nbFEDof(); iNode++)
        {
            int ig = pFESpace.dof().localToGlobal( eleID, iNode + 1 );
            tmpsum += velAndPressure(ig+offset);
            gid0Vec.push_back( p0FESpace.fe().currentId() );
            val0Vec.push_back( tmpsum / (double) pFESpace.fe().nbFEDof() );
        }
    }

    P0pres.setCoefficients(gid0Vec, val0Vec);
    P0pres.globalAssemble();
    MPI_Barrier(MPI_COMM_WORLD);

}
