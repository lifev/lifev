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
    @brief Test for building matrices with a static graph and ETA

    @author Radu Popescu <radu.popescu@epfl.ch>
    @date 2012-03-19
 */

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <cstdlib>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Epetra_FECrsGraph.h>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/WallClock.hpp>
#include <lifev/core/util/OpenMPParameters.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/expression/BuildGraph.hpp>


#include <boost/shared_ptr.hpp>


using namespace LifeV;

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);

    if (argc != 5) {
    	if (verbose) {
    		std::cout << "Please run program as " << argv[0]
					  << " " << "<num_elements> " << "<num_threads>"
					  << " " << "<scheduler> " << "<chunkSize>\n";
    		return EXIT_FAILURE;
    	}
    }

    OpenMPParameters ompParams;

    const UInt Nelements = std::atoi(argv[1]);
    ompParams.numThreads = std::atoi(argv[2]);
#ifdef _OPENMP
    ompParams.scheduler = static_cast<omp_sched_t>(std::atoi(argv[3]));
#endif
    ompParams.chunkSize = std::atoi(argv[4]);

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }


    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type);

    regularMesh3D ( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                    2.0,   2.0,   2.0,
                    -1.0,  -1.0,  -1.0);

//    MeshPartitioner< mesh_Type >  meshPart (fullMeshPtr, Comm);
    MeshPartitioner< mesh_Type >   meshPart;
    meshPart.setPartitionOverlap ( 1 );
    meshPart.doPartition ( fullMeshPtr, Comm );

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > uSpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPart, &feTetraP1, Comm) );

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << uSpace->dof().numTotalDof() << std::endl;
    }

    boost::shared_ptr<matrix_Type> closedSystemMatrix;

    if (verbose)
    {
        std::cout << " -- Precomputing matrix graph ... " << std::flush;
    }

    WallClock timer;
    timer.start();
    boost::shared_ptr<Epetra_FECrsGraph> matrixGraph;
    {
        using namespace ExpressionAssembly;

        // We first build a static graph for our problem matrix
        matrixGraph.reset(new Epetra_FECrsGraph(Copy, *(uSpace->map().map(Unique)), 0, true));

        buildGraph ( elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) )
                   ) >> matrixGraph;
    }
    matrixGraph->GlobalAssemble();
    timer.stop();

    if (verbose)
    {
        std::cout << " done in " << timer.elapsedTime() << "s." << std::endl;
    }

    if (verbose)
    {
        std::cout << " -- Assembling the Laplace matrix with a precomputed graph ... " << std::flush;
    }

    timer.start();
    {
        using namespace ExpressionAssembly;

        // We build the system matrix using the precomputed graph
        closedSystemMatrix.reset(new matrix_Type ( uSpace->map(), *matrixGraph , true) );
        *closedSystemMatrix *= 0.0;

        // Finally, we perform the FE assembly, which should be faster with
        // a closed matrix
        integrate (  elements (uSpace->mesh() ),
                     quadRuleTetra4pt,
                     uSpace,
                     uSpace,
                     dot ( grad (phi_i) , grad (phi_j) ), ompParams
                  ) >> closedSystemMatrix;

        closedSystemMatrix->globalAssemble();
    }
    timer.stop();

    if (verbose)
    {
        std::cout << " done in " << timer.elapsedTime() << "s." << std::endl;
    }

    Real closedMatrixNorm ( closedSystemMatrix->normInf() );

    if (verbose)
    {
        std::cout << " Closed matrix norm : " << closedMatrixNorm << std::endl;
    }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    Real closedMatrixNormDiff (std::abs (closedMatrixNorm - 3.2) );

    if (verbose)
    {
        std::cout << " Error (closed): " << closedMatrixNormDiff << std::endl;
    }

    Real testTolerance (1e-10);

    if ( closedMatrixNormDiff >= testTolerance )
    {
        return ( EXIT_FAILURE );
    }
    return ( EXIT_SUCCESS );

}
