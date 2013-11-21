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
 *  @brief Convert from Ensight to HDF5
 *
 *  @author Simone Deparis <simone.deparis@epfl.ch>
 *  @date 08-08-2008
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/navier_stokes/solver/OseenSolver.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>

#include <iostream>

#include "ensightToHdf5.hpp"

using namespace LifeV;

typedef RegionMesh<LinearTetra>                  mesh_Type;
typedef OseenSolver< mesh_Type >::vector_Type    vector_Type;
typedef OseenSolver< mesh_Type >::vectorPtr_Type vectorPtr_Type;
typedef FESpace< mesh_Type, MapEpetra >          feSpace_Type;
typedef boost::shared_ptr< feSpace_Type >        feSpacePtr_Type;


void
computeP0pressure (const feSpacePtr_Type& pFESpacePtr,
                   const feSpacePtr_Type& p0FESpacePtr,
                   const feSpacePtr_Type& uFESpacePtr,
                   const vector_Type& velAndPressureExport,
                   vector_Type& P0pres, Real /*time*/ );

struct EnsightToHdf5::Private
{
    Private() {}

    std::string    data_file_name;

    boost::shared_ptr<Epetra_Comm>   comm;
};

EnsightToHdf5::EnsightToHdf5 ( int argc,
                               char** argv )
    :
    d ( new Private )
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );

    d->data_file_name = data_file_name;

#ifdef EPETRA_MPI
    std::cout << "mpi initialization ... " << std::endl;

    // MPI_Init(&argc,&argv);
    d->comm.reset ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
    // int err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#else
    d->comm.reset ( new Epetra_SerialComm() );
#endif

    if (!d->comm->MyPID() )
    {
        std::cout << "My PID = " << d->comm->MyPID() << " out of " << d->comm->NumProc ()  << " running." << std::endl;
    }
}

void
EnsightToHdf5::run()
{

    // Reading from data file
    GetPot dataFile ( d->data_file_name.c_str() );

    bool verbose = (d->comm->MyPID() == 0);

    // Fluid solver
    boost::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( dataFile );

    MeshData meshData;
    meshData.setup (dataFile, "fluid/space_discretization");

    boost::shared_ptr<mesh_Type > fullMeshPtr (new mesh_Type ( d->comm ) );
    readMesh (*fullMeshPtr, meshData);

    // writeMesh("test.mesh", *fullMeshPtr);
    // Scale, Rotate, Translate (if necessary)
    boost::array< Real, NDIM >    geometryScale;
    boost::array< Real, NDIM >    geometryRotate;
    boost::array< Real, NDIM >    geometryTranslate;

    geometryScale[0] = dataFile ( "fluid/space_discretization/transform", 1., 0);
    geometryScale[1] = dataFile ( "fluid/space_discretization/transform", 1., 1);
    geometryScale[2] = dataFile ( "fluid/space_discretization/transform", 1., 2);

    geometryRotate[0] = dataFile ( "fluid/space_discretization/transform", 0., 3) * M_PI / 180;
    geometryRotate[1] = dataFile ( "fluid/space_discretization/transform", 0., 4) * M_PI / 180;
    geometryRotate[2] = dataFile ( "fluid/space_discretization/transform", 0., 5) * M_PI / 180;

    geometryTranslate[0] = dataFile ( "fluid/space_discretization/transform", 0., 6);
    geometryTranslate[1] = dataFile ( "fluid/space_discretization/transform", 0., 7);
    geometryTranslate[2] = dataFile ( "fluid/space_discretization/transform", 0., 8);

    MeshUtility::MeshTransformer<mesh_Type, mesh_Type::markerCommon_Type > _transformMesh (*fullMeshPtr);
    _transformMesh.transformMesh ( geometryScale, geometryRotate, geometryTranslate );

    boost::shared_ptr<mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, d->comm);
        meshPtr = meshPart.meshPartition();
    }

    std::string uOrder =  dataFile ( "fluid/space_discretization/vel_order", "P1");
    std::string pOrder =  dataFile ( "fluid/space_discretization/press_order", "P1");

    //oseenData->meshData()->setMesh(meshPtr);

    if (verbose)
    {
        std::cout << "Building the velocity FE space ... " << std::flush;
    }

    feSpacePtr_Type uFESpacePtr ( new feSpace_Type (meshPtr, uOrder, 3, d->comm) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the pressure FE space ... " << std::flush;
    }

    feSpacePtr_Type pFESpacePtr ( new feSpace_Type ( meshPtr, pOrder, 1, d->comm ) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    if (verbose)
    {
        std::cout << "Building the P0 pressure FE space ... " << std::flush;
    }

    feSpacePtr_Type p0FESpacePtr ( new feSpace_Type ( meshPtr, feTetraP0, quadRuleTetra1pt,
                                                      quadRuleTria1pt, 1, d->comm ) );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }



    UInt totalVelDof   = uFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt totalPressDof = pFESpacePtr->map().map (Unique)->NumGlobalElements();
    UInt totalP0PresDof = p0FESpacePtr->map().map (Unique)->NumGlobalElements();

    if (verbose)
    {
        std::cout << "Total Velocity DOF = " << totalVelDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total Pressure DOF = " << totalPressDof << std::endl;
    }
    if (verbose)
    {
        std::cout << "Total P0Press  DOF = " << totalP0PresDof << std::endl;
    }

    if (verbose)
    {
        std::cout << "Calling the fluid constructor ... ";
    }

    OseenSolver< mesh_Type > fluid (oseenData,
                                    *uFESpacePtr,
                                    *pFESpacePtr,
                                    d->comm);
    MapEpetra fullMap (fluid.getMap() );

    if (verbose)
    {
        std::cout << "ok." << std::endl;
    }

    fluid.setUp (dataFile);

    // Initialization
    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();

    boost::shared_ptr< Exporter<mesh_Type > > exporter;
    boost::shared_ptr< Exporter<mesh_Type > > importer;

    std::string const exporterType =  dataFile ( "exporter/type", "hdf5");
    std::string const exporterName =  dataFile ( "exporter/filename", "ethiersteinman");
    std::string const exportDir    =  dataFile ( "exporter/exportDir", "./");

    std::string const importerType =  dataFile ( "importer/type", "ensight");
    std::string const importerName =  dataFile ( "importer/filename", "ethiersteinman");
    std::string const importDir    =  dataFile ( "importer/importDir", "importDir/");

#ifdef HAVE_HDF5
    if (exporterType.compare ("hdf5") == 0)
    {
        exporter.reset ( new ExporterHDF5<mesh_Type > ( dataFile, exporterName ) );
    }
    else
#endif
        exporter.reset ( new ExporterEnsight<mesh_Type > ( dataFile, exporterName ) );

    exporter->setPostDir ( exportDir ); // This is a test to see if M_post_dir is working
    exporter->setMeshProcId ( meshPtr, d->comm->MyPID() );

#ifdef HAVE_HDF5
    if (importerType.compare ("hdf5") == 0)
    {
        importer.reset ( new ExporterHDF5<mesh_Type > ( dataFile, importerName ) );
    }
    else
#endif
        importer.reset ( new ExporterEnsight<mesh_Type > ( dataFile, importerName ) );

    // todo this will not work with the ExporterEnsight filter (it uses M_importDir, a private variable)
    importer->setPostDir ( importDir ); // This is a test to see if M_post_dir is working
    importer->setMeshProcId ( meshPtr, d->comm->MyPID() );

    vectorPtr_Type velAndPressureExport ( new vector_Type (*fluid.solution(), exporter->mapType() ) );
    vectorPtr_Type velAndPressureImport ( new vector_Type (*fluid.solution(), importer->mapType() ) );

    if ( exporter->mapType() == Unique )
    {
        velAndPressureExport->setCombineMode (Zero);
    }

    importer->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                            velAndPressureImport, UInt (0) );
    importer->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                            velAndPressureImport, UInt (3 * uFESpacePtr->dof().numTotalDof() ) );
    importer->import ( t0 );

    *velAndPressureExport = *velAndPressureImport;


    vectorPtr_Type P0pres ( new vector_Type (p0FESpacePtr->map() ) );
    MPI_Barrier (MPI_COMM_WORLD);
    computeP0pressure (pFESpacePtr, p0FESpacePtr, uFESpacePtr, *velAndPressureImport, *P0pres, t0);

    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpacePtr,
                            velAndPressureExport, UInt (0) );

    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpacePtr,
                            velAndPressureExport, UInt (3 * uFESpacePtr->dof().numTotalDof() ) );

    exporter->addVariable (ExporterData<mesh_Type>::ScalarField, "P0pressure", p0FESpacePtr,
                           P0pres, UInt (0),
                           ExporterData<mesh_Type>::SteadyRegime, ExporterData<mesh_Type>::Cell );
    exporter->postProcess ( t0 );

    // Temporal loop
    LifeChrono chrono;
    int iter = 1;

    for ( Real time = t0 + dt ; time <= tFinal + dt / 2.; time += dt, iter++)
    {
        chrono.stop();
        importer->import ( time );

        *velAndPressureExport = *velAndPressureImport;
        MPI_Barrier (MPI_COMM_WORLD);
        computeP0pressure (pFESpacePtr, p0FESpacePtr, uFESpacePtr, *velAndPressureImport, *P0pres, time);

        exporter->postProcess ( time );

        chrono.stop();
        if (verbose)
        {
            std::cout << "Total iteration time " << chrono.diff() << " s." << std::endl;
        }
    }
}


void
computeP0pressure (const feSpacePtr_Type& pFESpacePtr,
                   const feSpacePtr_Type& p0FESpacePtr,
                   const feSpacePtr_Type& uFESpacePtr,
                   const vector_Type&     velAndPressure,
                   vector_Type&           P0pres,
                   Real /*time*/)
{

    int MyPID;
    MPI_Comm_rank (MPI_COMM_WORLD, &MyPID);
    UInt offset = 3 * uFESpacePtr->dof().numTotalDof();

    std::vector<int> gid0Vec (0);
    gid0Vec.reserve (p0FESpacePtr->mesh()->numVolumes() );
    std::vector<Real> val0Vec (0);
    val0Vec.reserve (p0FESpacePtr->mesh()->numVolumes() );

    for (UInt ivol = 0; ivol < pFESpacePtr->mesh()->numVolumes(); ++ivol)
    {

        pFESpacePtr->fe().update ( pFESpacePtr->mesh()->volumeList ( ivol ), UPDATE_DPHI );
        p0FESpacePtr->fe().update ( p0FESpacePtr->mesh()->volumeList ( ivol ) );

        UInt eleID = pFESpacePtr->fe().currentLocalId();

        double tmpsum = 0.;
        for (UInt iNode = 0; iNode < (UInt) pFESpacePtr->fe().nbFEDof(); iNode++)
        {
            int ig = pFESpacePtr->dof().localToGlobalMap ( eleID, iNode );
            tmpsum += velAndPressure (ig + offset);
            gid0Vec.push_back ( p0FESpacePtr->fe().currentId() );
            val0Vec.push_back ( tmpsum / (double) pFESpacePtr->fe().nbFEDof() );
        }
    }
    P0pres.setCoefficients (gid0Vec, val0Vec);
    P0pres.globalAssemble();
    MPI_Barrier (MPI_COMM_WORLD);

}
