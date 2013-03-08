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
#include <lifev/fsi/solver/RBFInterpolation.hpp>         // Nuovo approccio per interpolazione con RBF
#include <lifev/fsi/solver/RBFInterpolationRadius.hpp>   // Interpolazione con RBF con scelta raggio e divisione per interpolante della costante uguale a 1
#include <lifev/fsi/solver/RBFInterpolationStandard.hpp> // Interpolazione con RBF standard
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

#define PI 3.14159265

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace LifeV;

// f -> function to be interpolated
double f (double x, double y, double z)
{
    return sin (2 * PI * x) * cos (3 * PI * y) + exp (x * y);
}

typedef VectorEpetra                             vector_Type;
typedef boost::shared_ptr<vector_Type >          vectorPtr_Type;
typedef RegionMesh<LinearTriangle >              mesh_Type;
typedef boost::shared_ptr< mesh_Type >           meshPtr_Type;
typedef RBFInterpolation< mesh_Type >            interpolation_Type;
typedef boost::shared_ptr< interpolation_Type >  interpolationPtr_Type;

typedef RBFInterpolationRadius< mesh_Type >       interpolationR_Type;
typedef boost::shared_ptr< interpolationR_Type >  interpolationRPtr_Type;

typedef RBFInterpolationStandard< mesh_Type >     interpolationS_Type;
typedef boost::shared_ptr< interpolationS_Type >  interpolationSPtr_Type;

int main (int argc, char** argv )
{
    boost::shared_ptr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    Comm.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    comm.reset ( new Epetra_SerialComm() );
#endif

    // DATAFILE
    GetPot command_line (argc, argv);
    GetPot dataFile ( command_line.follow ("data", 2, "-f", "--file" ) );

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

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
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Solid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Solid_localMesh, "P1", 1, Comm) );
    vectorPtr_Type Solid_vector (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Solid_vector_one (new vector_Type (Solid_fieldFESpace->map(), Unique) );

    for ( UInt i = 0; i < Solid_vector->epetraVector().MyLength(); ++i )
        if (Solid_vector->blockMap().LID (Solid_vector->blockMap().GID (i) ) != -1)
            (*Solid_vector) [Solid_vector->blockMap().GID (i)] = f ( Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).x(),
                                                                     Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).y(),
                                                                     Solid_mesh_ptr->point (Solid_vector->blockMap().GID (i) ).z() );

    // EXPORTING THE DEFINED FIELD
    ExporterHDF5<mesh_Type> Solid_exporter (dataFile, Solid_localMesh, "Input field", Comm->MyPID() );
    Solid_exporter.setMeshProcId (Solid_localMesh, Comm->MyPID() );
    Solid_exporter.exportPID (Solid_localMesh, Comm, true );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "f(x,y,z)", Solid_fieldFESpace, Solid_vector, UInt (0) );
    Solid_exporter.postProcess (0);
    Solid_exporter.closeFile();

    // HERE, THE PIECE OF CODE RELATED TO THE RADIAL BASIS FUNCTION INTERPOLATION BEGINS. FIRSTLY WE DEFINE ON WHICH PART OF THE DOMAIN WE ARE
    // GOING TO INTERPOLATE
    int nFlags = 1;
    std::vector<int> flags (nFlags);
    flags[0] = -1;

    // DEFINING SOME STUFF FOR EVALUATING THE SOLUTION. I CREATE A FE SPACE ON THE MESH WHERE I WANT TO GET THE SOLUTION, A VECTOR THAT WILL
    // CONTAIN THE SOLUTION AND ANOTHER ONE CONTAINING THE SOLUTION GAINED BY RBF. THIS IS DUE TO THE FACT THAT WE HAVE SLIGHTLY MODIFIED
    // (POSITIVELY) THE RBF'S ORIGINAL APPROACH.
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Fluid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Fluid_localMesh, "P1", 1, Comm) );
    vectorPtr_Type Fluid_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Fluid_solution_rbf (new vector_Type (Fluid_fieldFESpace->map(), Unique) );


    // INITIALIZATION: THE FIRST TWO ARGUMENTS ARE RELATED TO THE MESHES WHERE WE KNOW THE FUNCTION. THE 3RD AND THE 4TH TO THE MESHES WHERE
    // WE ARE GOING TO EVALUATE THE FUNCTION. THE LAST IS USED TO SPECIFY THE FLAGS, SO THOSE PARTS OF THE MESHES THAT WE ARE CONSIDERING.
    interpolationPtr_Type RBFInterpolant ( new interpolation_Type ( Solid_mesh_ptr,
                                                                    Solid_localMesh,
                                                                    Fluid_mesh_ptr,
                                                                    Fluid_localMesh,
                                                                    flags) );

    /*
    // IN QUESTO CASO USI IL RAGGIO E NON LA SELEZIONE AUTOMATICA, SOLO CHE DIVIDI SEMPRE PER L'INTERPOLANTE DI 1
    double radius = 2*(double)MeshUtility::MeshStatistics::computeSize(*Solid_mesh_ptr).minH; // maxH, meanH
    interpolationRPtr_Type RBFInterpolant ( new interpolationR_Type ( Solid_mesh_ptr,
                                                                      Solid_localMesh,
                                                                      Fluid_mesh_ptr,
                                                                      Fluid_localMesh,
                                                                      flags,
                                                                      radius) );

    */
    /*
    // IN QUESTO CASO USI IL RAGGIO E NON LA SELEZIONE AUTOMATICA, SOLO CHE DIVIDI SEMPRE PER L'INTERPOLANTE DI 1
    double radius = (double)MeshUtility::MeshStatistics::computeSize (*Solid_mesh_ptr).minH; // maxH, meanH
    interpolationSPtr_Type RBFInterpolant ( new interpolationS_Type ( Solid_mesh_ptr,
                                                                      Solid_localMesh,
                                                                      Fluid_mesh_ptr,
                                                                      Fluid_localMesh,
                                                                      flags,
                                                                      radius) );
    */

    // LOADING INFORMATION ABOUT THE TWO VECTORS INVOLVED IN THE INTERPOLATION PROCESS
    RBFInterpolant->setupRBFData (Solid_vector, Fluid_solution, dataFile, belosList);

    LifeChrono buildOperatorChrono;
    buildOperatorChrono.start();
    // BUILDING THE OPERATORS
    RBFInterpolant->buildOperators();
    buildOperatorChrono.stop();
    std::cout << "Time to build operators = " << buildOperatorChrono.diff() << std::endl;

    // COMPUTING THE SOLUTION
    LifeChrono solutionChrono;
    solutionChrono.start();
    RBFInterpolant->interpolate();
    solutionChrono.stop();
    std::cout << "Time to solve the linear systems = " << solutionChrono.diff() << std::endl;


    // SAVE THE SOLUTION
    RBFInterpolant->solution (Fluid_solution);

    // SAVE THE RBF'S ORIGINAL SOLUTION
    // RBFInterpolant->solutionrbf (Fluid_solution_rbf);

    // COMPUTING THE ERROR
    vectorPtr_Type Fluid_exact_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type myError (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    // vectorPtr_Type rbfError (new vector_Type (Fluid_fieldFESpace->map(), Unique) );

    for ( UInt i = 0; i < Fluid_exact_solution->epetraVector().MyLength(); ++i )
        if (Fluid_exact_solution->blockMap().LID (Fluid_exact_solution->blockMap().GID (i) ) != -1)
        {
            (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i)] = f ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(),
                                                                                     Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(),
                                                                                     Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );

            (*myError) [myError->blockMap().GID (i)] = (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i)] - (*Fluid_solution) [Fluid_solution->blockMap().GID (i)];
            // (*rbfError) [rbfError->blockMap().GID (i)] = (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i)] - (*Fluid_solution_rbf) [Fluid_solution_rbf->blockMap().GID (i)];

        }

    // EXPORTING THE SOLUTION
    ExporterHDF5<mesh_Type> Fluid_exporter (dataFile, Fluid_localMesh, "Results", Comm->MyPID() );
    Fluid_exporter.setMeshProcId (Fluid_localMesh, Comm->MyPID() );
    Fluid_exporter.exportPID (Fluid_localMesh, Comm, true );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Exact solution", Fluid_fieldFESpace, Fluid_exact_solution, UInt (0) );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Solution", Fluid_fieldFESpace, Fluid_solution, UInt (0) );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Error", Fluid_fieldFESpace, myError, UInt (0) );
    // Fluid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "RBF's solution", Fluid_fieldFESpace, Fluid_solution_rbf, UInt (0) );
    // Fluid_exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "RBF's error", Fluid_fieldFESpace, rbfError, UInt (0) );
    Fluid_exporter.postProcess (0);
    Fluid_exporter.closeFile();

    std::cout << "MinH = " <<  (double)MeshUtility::MeshStatistics::computeSize(*Solid_mesh_ptr).minH << std::endl;
    std::cout << "MeanH = " <<  (double)MeshUtility::MeshStatistics::computeSize(*Solid_mesh_ptr).meanH << std::endl;
    std::cout << "MaxH = " <<  (double)MeshUtility::MeshStatistics::computeSize(*Solid_mesh_ptr).maxH << std::endl;
    std::cout << "Number of vertices = " << Solid_mesh_ptr->numVertices() << std::endl;
    std::cout << "Norm Inf = " << myError->normInf() << std::endl;
    std::cout << "Norm 2 = " << ( myError->norm2()/Fluid_exact_solution->norm2() ) << std::endl;//x* 1/Solid_mesh_ptr->numVertices() << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return (EXIT_SUCCESS);

}
