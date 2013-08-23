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

@author Davide Forti <davide.forti@epfl.ch>
@date 02-21-2013
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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/filter/Exporter.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <sys/stat.h>

#define PI 3.14159265

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace LifeV;

int main (int argc, char** argv )
{
    typedef VectorEpetra                          vector_Type;
    typedef boost::shared_ptr<vector_Type >       vectorPtr_Type;
    typedef RegionMesh<LinearTetra >              mesh_Type;
    typedef boost::shared_ptr< mesh_Type >        meshPtr_Type;
    typedef RBFInterpolation<mesh_Type>           interpolation_Type;
    typedef boost::shared_ptr<interpolation_Type> interpolationPtr_Type;

    boost::shared_ptr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    Comm.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    Comm.reset ( new Epetra_SerialComm() );
#endif

    // DATAFILE
    GetPot command_line (argc, argv);
    GetPot dataFile ( command_line.follow ("data_rbf3d", 2, "-f", "--file" ) );

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf3d.xml" );
    
    std::string problemFolder = command_line.follow ( "Output", 2, "-o", "--output" );
    if ( problemFolder.compare ("./") )
    {
        problemFolder += "/";
        
        if ( Comm->MyPID() == 0 )
        {
            mkdir ( problemFolder.c_str(), 0777 );
        }
    }

    // LOADING MESHES
    MeshData MeshCoarseData;
    meshPtr_Type MeshCoarse ( new mesh_Type ( Comm ) );
    MeshCoarseData.setup (dataFile, "coarse/space_discretization");
    readMesh (*MeshCoarse, MeshCoarseData);

    MeshData MeshFineData;
    meshPtr_Type MeshFine ( new mesh_Type ( Comm ) );
    MeshFineData.setup (dataFile, "fine/space_discretization");
    readMesh (*MeshFine, MeshFineData);
     
    // PARTITIONING MESHES
    MeshPartitioner<mesh_Type>   MeshCoarsePart;
    boost::shared_ptr<mesh_Type> MeshCoarseLocal;
    MeshCoarsePart.setPartitionOverlap (0);
    MeshCoarsePart.doPartition (MeshCoarse, Comm);
    MeshCoarseLocal = MeshCoarsePart.meshPartition();

    MeshPartitioner<mesh_Type>   MeshFinePart;
    boost::shared_ptr<mesh_Type> MeshFineLocal;
    MeshFinePart.setPartitionOverlap (0);
    MeshFinePart.doPartition (MeshFine, Comm);
    MeshFineLocal = MeshFinePart.meshPartition();
    
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > MeshCoarseFESpace (new FESpace<mesh_Type, MapEpetra> (MeshCoarseLocal, "P1", 3, Comm) );
    vectorPtr_Type fibers (new vector_Type (MeshCoarseFESpace->map(), Unique) );
    
    //SETTING IMPORTER
    typedef ExporterData<mesh_Type> exporterData_Type;
    exporterData_Type impData (exporterData_Type::VectorField, "fibers.00000", MeshCoarseFESpace ,
                               fibers, UInt (0), exporterData_Type::UnsteadyRegime);
    typedef boost::shared_ptr< Exporter<mesh_Type> > filterPtr_Type;
    typedef ExporterHDF5<mesh_Type> hdf5Filter_Type;
    filterPtr_Type importer ( new hdf5Filter_Type() );
    importer -> setMeshProcId ( MeshCoarseLocal, Comm->MyPID() );
    importer-> setPrefix ("fibers-60+60lasca");
    importer -> readVariable (impData);
    importer -> closeFile();
    
    /*
    // EXPORTING THE DEFINED FIELD
    ExporterHDF5<mesh_Type> Solid_exporter (dataFile, Solid_localMesh, "Input field", Comm->MyPID() );
    Solid_exporter.setMeshProcId (Solid_localMesh, Comm->MyPID() );
    Solid_exporter.exportPID (Solid_localMesh, Comm, true );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "f(x,y,z)", Solid_fieldFESpace, Solid_vector, UInt (0) );
    Solid_exporter.postProcess (0);
    Solid_exporter.closeFile();
    */

    int nFlags = 1;
    std::vector<int> flags (nFlags);
    flags[0] = -1;

    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > MeshFineFESpace (new FESpace<mesh_Type, MapEpetra> (MeshFineLocal, "P1", 3, Comm) );
    
    vectorPtr_Type fibers_sol (new vector_Type (MeshFineFESpace->map(), Unique) );
    
    interpolationPtr_Type RBFinterpolant;

    RBFinterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject (dataFile("interpolation/interpolation_Type","none")));

    RBFinterpolant->setup(MeshCoarse, MeshCoarseLocal, MeshFine, MeshFineLocal, flags);

    RBFinterpolant->setupRBFData (fibers, fibers_sol, dataFile, belosList);

    RBFinterpolant->buildOperators();

    RBFinterpolant->interpolate();

    RBFinterpolant->solution (fibers_sol);

    // EXPORTING THE SOLUTION
    ExporterHDF5<mesh_Type> Fluid_exporter (dataFile, MeshFineLocal, "Results", Comm->MyPID() );
    Fluid_exporter.setMeshProcId (MeshFineLocal, Comm->MyPID() );
    Fluid_exporter.setPostDir(problemFolder);
    Fluid_exporter.exportPID (MeshFineLocal, Comm, true );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Solution", MeshFineFESpace, fibers_sol, UInt (0) );
    Fluid_exporter.postProcess (0);
    Fluid_exporter.closeFile();
    
     
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return (EXIT_SUCCESS);

}
