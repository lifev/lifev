//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
    @brief Test for the MeshPartitionTool class

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 14-11-2012
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/mesh/MeshPartitionTool.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;
typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
typedef std::shared_ptr<feSpace_Type> feSpacePtr_Type;
typedef MeshPartitionTool<mesh_Type> meshCutter_Type;

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);

    // Read first the data needed

    GetPot cl (argc, argv);
    const UInt numElements = cl.follow (9, "--num-elem");
    // partitionerType should be MeshPartitioner, MeshPartitionTool_ParMETIS or
    // MeshPartitionTool_Zoltan
    const std::string graphLib = cl.follow ("parmetis", "--graph-lib");

    if (verbose) std::cout << " ---> Number of elements : "
                               << numElements << std::endl;

    // Build and partition the mesh

    if (verbose)
    {
        std::cout << " -- Building the mesh ... " << std::flush;
    }
    std::shared_ptr< mesh_Type > fullMeshPtr (new RegionMesh<LinearTetra>);
    std::shared_ptr< mesh_Type > meshPart;
    regularMesh3D ( *fullMeshPtr, 1, numElements, numElements, numElements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Partitioning the mesh ... " << std::flush;
    }
    {
        Teuchos::ParameterList meshParameters;
        meshParameters.set ("num-parts", Comm->NumProc(), "");
        meshParameters.set ("graph-lib", graphLib, "");
        meshCutter_Type meshCutter (fullMeshPtr, Comm, meshParameters);
        if (! meshCutter.success() )
        {
            if (verbose)
            {
                std::cout << "Partitioning failed." << std::endl;
            }
            return EXIT_FAILURE;
        }
        meshPart = meshCutter.meshPart();
    }
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Freeing the global mesh ... " << std::flush;
    }
    fullMeshPtr.reset();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    // Build the FESpaces
    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }
    std::string uOrder ("P1");
    std::string bOrder ("P1");
    std::shared_ptr < FESpace < mesh_Type,
          MapEpetra > >
          uFESpace (new FESpace < mesh_Type,
                    MapEpetra > (meshPart,
                                 uOrder,
                                 1,
                                 Comm) );

    std::shared_ptr < FESpace < mesh_Type,
          MapEpetra > >
          betaFESpace (new FESpace < mesh_Type,
                       MapEpetra > (meshPart,
                                    bOrder,
                                    3,
                                    Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose) std::cout << " ---> Dofs: "
                               << uFESpace->dof().numTotalDof() << std::endl;

    // Build the assembler and the matrices

    if (verbose)
    {
        std::cout << " -- Building assembler ... " << std::flush;
    }
    ADRAssembler<mesh_Type, matrix_Type, vector_Type> adrAssembler;
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Setting up assembler ... " << std::flush;
    }
    adrAssembler.setup (uFESpace, betaFESpace);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }
    std::shared_ptr<matrix_Type>
    systemMatrix (new matrix_Type (uFESpace->map() ) );
    *systemMatrix *= 0.0;
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }

    // Perform the assembly of the matrix

    if (verbose)
    {
        std::cout << " -- Adding the diffusion ... " << std::flush;
    }
    adrAssembler.addDiffusion (systemMatrix, 1.0);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " Time needed : "
                  << adrAssembler.diffusionAssemblyChrono().diffCumul()
                  << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
    systemMatrix->globalAssemble();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    Real matrixNorm (systemMatrix->normFrobenius() );
    if (verbose)
    {
        std::cout << " ---> Norm 1 : " << matrixNorm << std::endl;
    }
    if ( std::fabs (matrixNorm - 35.908 ) > 1e-3)
    {
        std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}


