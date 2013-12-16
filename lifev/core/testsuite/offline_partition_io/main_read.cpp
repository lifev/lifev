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
    @brief Test for PartitionIO class - read and solve

    @author Radu Popescu <radu.popescu@epfl.ch>
    @maintainer Radu Popescu <radu.popescu@epfl.ch>
    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 10-05-2012

    Loads mesh parts from HDF5 and solves a Laplacian problem.
    Based on the ADRAssembler class unit test.
 */

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAS_HDF5

#include "Epetra_config.h"

#ifdef HAVE_MPI

#include <mpi.h>

#include <Epetra_MpiComm.h>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/PartitionIO.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>

using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef VectorEpetra vector_Type;

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

int
main ( int argc, char** argv )
{
#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_MpiComm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );

    const bool verbose (comm->MyPID() == 0);

    // Read first the data needed

    if (verbose)
    {
        std::cout << " -- Reading the data ... " << std::flush;
    }
//    GetPot dataFile ( "data" );
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    GetPot cl (argc, argv);
    // partitionerType should be MeshPartitioner, MeshPartitionTool_ParMETIS or
    // MeshPartitionTool_Zoltan
    const std::string partitionerType = cl.follow("MeshPartitioner",
    											  "--partitioner-type");
    std::string partsFile;
    partsFile.reserve(50);
    partsFile += "cube_";
    partsFile += partitionerType;
    partsFile += ".h5";

    boost::shared_ptr<mesh_Type> mesh;
    {
        PartitionIO<RegionMesh<LinearTetra> > partitionIO (partsFile, comm);
        partitionIO.read (mesh);
    }

    // Build the FESpaces

    if (verbose)
    {
        std::cout << " -- Building FESpaces ... " << std::flush;
    }
    std::string uOrder ("P1");
    std::string bOrder ("P1");
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> >
    uFESpace (new FESpace<mesh_Type, MapEpetra> (mesh, uOrder, 1, comm) );
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> >
    betaFESpace (new FESpace<mesh_Type, MapEpetra> (mesh, bOrder, 3, comm) );
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
    boost::shared_ptr<matrix_Type>
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
    adrAssembler.addDiffusion (systemMatrix, 1);
    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    if (verbose) std::cout << " Time needed : "
                               << adrAssembler.diffusionAssemblyChrono().diffCumul()
                               << std::endl;

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
    systemMatrix->globalAssemble();
    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    Real matrixNorm (systemMatrix->normFrobenius());
    if (verbose)
    {
        std::cout << " ---> Norm 2 : " << matrixNorm << std::endl;
    }
    if (std::fabs (matrixNorm - 35.908) > 1e-3)
    {
        std::cout << " <!> Matrix has changed !!! <!> " << std::endl;
        return EXIT_FAILURE;
    }

    if (verbose)
    {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }

    MPI_Finalize();

#else
    std::cout << "This test needs MPI to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* HAVE_MPI */
#else
    std::cout << "This test needs HDF5 to run. Aborting." << std::endl;
    return (EXIT_FAILURE);
#endif /* LIFEV_HAS_HDF5 */

    return ( EXIT_SUCCESS );
}


