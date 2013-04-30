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
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFhtpVectorial.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFvectorial.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#define PI 3.14159265

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace LifeV;

double fx (double x, double y, double z)
{
    return sin (2 * PI * x) * cos (3 * PI * y) + exp (x * y) + z;
}

double fy (double x, double y, double z)
{
    return sin (4 * PI * y) * cos ( PI * x * y) + exp (x * 2) + z;
}

double fz (double x, double y, double z)
{
    return sin (0.5 * PI * x) * cos (3 * PI * z) + exp (y * y);
}

int main (int argc, char** argv )
{
    typedef VectorEpetra                          vector_Type;
    typedef boost::shared_ptr<vector_Type >       vectorPtr_Type;
    typedef RegionMesh<LinearTriangle >           mesh_Type;
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
    GetPot dataFile ( command_line.follow ("data_rbf2d", 2, "-f", "--file" ) );

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf2d.xml" );

    // LOADING MESHES
    MeshData Solid_mesh_data;
    meshPtr_Type Solid_mesh_ptr ( new mesh_Type ( Comm ) );
    Solid_mesh_data.setup (dataFile, "solid/space_discretization");
    readMesh (*Solid_mesh_ptr, Solid_mesh_data);

    MeshData Fluid_mesh_data;
    meshPtr_Type Fluid_mesh_ptr ( new mesh_Type ( Comm ) );
    Fluid_mesh_data.setup (dataFile, "fluid/space_discretization");
    readMesh (*Fluid_mesh_ptr, Fluid_mesh_data);

    // PARTITIONING MESHES
    MeshPartitioner<mesh_Type>   Solid_mesh_part;
    boost::shared_ptr<mesh_Type> Solid_localMesh;
    Solid_mesh_part.setPartitionOverlap (2);
    Solid_mesh_part.doPartition (Solid_mesh_ptr, Comm);
    Solid_localMesh = Solid_mesh_part.meshPartition();

    MeshPartitioner<mesh_Type>   Fluid_mesh_part;
    boost::shared_ptr<mesh_Type> Fluid_localMesh;
    Fluid_mesh_part.setPartitionOverlap (2);
    Fluid_mesh_part.doPartition (Fluid_mesh_ptr, Comm);
    Fluid_localMesh = Fluid_mesh_part.meshPartition();

    // CREATING A FE-SPACE FOR THE GRID ON WHICH WE ASSUME TO KNOW THE INTERFACE FIELD. DEFINING AN INTERFACE VECTOR TO BE INTERPOLATED.
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Solid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Solid_localMesh, "P1", 3, Comm) );
    vectorPtr_Type Solid_vector (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Solid_vector_one (new vector_Type (Solid_fieldFESpace->map(), Unique) );

    for ( UInt j = 0; j < 3; ++j)
        for ( UInt i = 0; i < Solid_vector->epetraVector().MyLength(); ++i )
            if( Solid_vector->blockMap().GID (i) < Solid_mesh_ptr->pointList.size())
                if (Solid_vector->blockMap().LID (Solid_vector->blockMap().GID(i) + Solid_vector->size()/3*j ) != -1)
                {
                    switch (j)
                    {
                    case(0):
                        (*Solid_vector) [Solid_vector->blockMap().GID (i) + Solid_vector->size()/3*j ] = fx ( Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).x(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).y(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).z() );
                        break;
                    case(1):
                        (*Solid_vector) [Solid_vector->blockMap().GID (i) + Solid_vector->size()/3*j ] = fy ( Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).x(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).y(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).z() );
                        break;
                    case(2):
                        (*Solid_vector) [Solid_vector->blockMap().GID (i) + Solid_vector->size()/3*j ] = fz ( Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).x(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).y(),
                                                                                                              Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).z() );
                        break;
                    }
                }

    // EXPORTING THE DEFINED FIELD
    ExporterHDF5<mesh_Type> Solid_exporter (dataFile, Solid_localMesh, "Input field", Comm->MyPID() );
    Solid_exporter.setMeshProcId (Solid_localMesh, Comm->MyPID() );
    Solid_exporter.exportPID (Solid_localMesh, Comm, true );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "f(x,y,z)", Solid_fieldFESpace, Solid_vector, UInt (0) );
    Solid_exporter.postProcess (0);
    Solid_exporter.closeFile();

    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Fluid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Fluid_localMesh, "P1", 3, Comm) );
    vectorPtr_Type Fluid_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Fluid_solution_rbf (new vector_Type (Fluid_fieldFESpace->map(), Unique) );

    int nFlags = 1;
    std::vector<int> flags (nFlags);
    flags[0] = -1;

    interpolationPtr_Type RBFinterpolant;

    RBFinterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject (dataFile("interpolation/interpolation_Type","none")));

    RBFinterpolant->setup(Solid_mesh_ptr, Solid_localMesh, Fluid_mesh_ptr, Fluid_localMesh, flags);

    if(dataFile("interpolation/interpolation_Type","none")=="RBFvectorial")
        RBFinterpolant->setBasis(dataFile("interpolation/basis","none"));

    if(dataFile("interpolation/interpolation_Type","none")!="RBFlocallyRescaledVectorial")
        RBFinterpolant->setRadius((double) MeshUtility::MeshStatistics::computeSize (*Solid_mesh_ptr).maxH);

    RBFinterpolant->setupRBFData (Solid_vector, Fluid_solution, dataFile, belosList);

    RBFinterpolant->buildOperators();

    RBFinterpolant->interpolate();

    RBFinterpolant->solution (Fluid_solution);

    // COMPUTING THE ERROR
    vectorPtr_Type Fluid_exact_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type myError (new vector_Type (Fluid_fieldFESpace->map(), Unique) );

    for ( UInt j = 0; j < 3; ++j)
        for ( UInt i = 0; i < Fluid_exact_solution->epetraVector().MyLength(); ++i )
            if ( Fluid_exact_solution->blockMap().GID (i) < Fluid_mesh_ptr->pointList.size())
                if (Fluid_exact_solution->blockMap().LID (Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ) != -1)
                {
                    switch (j)
                    {
                    case(0):
                        (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fx ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                        break;
                    case(1):
                        (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fy ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                        break;
                    case(2):
                        (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fz ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                        break;
                    }

                    (*myError) [myError->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] - (*Fluid_solution) [Fluid_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ];
		}

    Real err_Inf = myError->normInf();
    Real err_L2  = myError->norm2()/Solid_mesh_ptr->numVertices();

    if(Comm->MyPID()==0)
    {
        std::cout << "Error, norm_Inf = " <<  err_Inf  << std::endl;
        std::cout << "Error, normL2   = " <<  err_L2   << std::endl;
    }

    // EXPORTING THE SOLUTION
    ExporterHDF5<mesh_Type> Fluid_exporter (dataFile, Fluid_localMesh, "Results", Comm->MyPID() );
    Fluid_exporter.setMeshProcId (Fluid_localMesh, Comm->MyPID() );
    Fluid_exporter.exportPID (Fluid_localMesh, Comm, true );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Exact solution", Fluid_fieldFESpace, Fluid_exact_solution, UInt (0) );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Solution", Fluid_fieldFESpace, Fluid_solution, UInt (0) );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Error", Fluid_fieldFESpace, myError, UInt (0) );
    Fluid_exporter.postProcess (0);
    Fluid_exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return (EXIT_SUCCESS);

}
