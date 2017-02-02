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

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.cpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/filter/ParserGmsh.hpp>
#include <lifev/core/mesh/ConvertBareMesh.hpp>

#include <cstdlib>
#include <iostream>

typedef std::shared_ptr<Epetra_Comm> comm_t;

template <int dim>
struct partitioner
{
    template<typename M>
    static void dopartition (std::shared_ptr<M>& mesh_p, comm_t comm)
    {
        mesh_p->updateElementFacets (true, true);
        mesh_p->updateElementRidges (true, true);

        LifeV::MeshPartitioner<M> part (mesh_p, comm);
    }
};

template <>
struct partitioner<2>
{
    template<typename M>
    static void dopartition (std::shared_ptr<M>&, comm_t)
    {
        // TODO: not working yet, problem with boundary ridges (aka points).
        /*
        mesh_p->updateElementFacets(true, true);

        LifeV::MeshPartitioner<M> part(mesh_p, comm);
        */
    }
};

template <>
struct partitioner<1>
{
    template<typename M>
    static void dopartition (std::shared_ptr<M>&, comm_t)
    {}
};


template <typename S>
struct tester
{
    typedef typename LifeV::BareMesh<S>    baremesh_t;
    typedef typename LifeV::RegionMesh<S>  mesh_t;

    static bool test (std::string& filename, comm_t comm, bool part = true)
    {
        baremesh_t baremesh;

        if (!LifeV::MeshIO::ReadGmshFile (filename, baremesh, 1, true) )
        {
            return false;
        }

        // Convert the mesh
        std::shared_ptr<mesh_t> mesh (new mesh_t (comm) );
        LifeV::convertBareMesh (baremesh, *mesh, true);

        // Partitioning (if possible)
        if (part)
        {
            partitioner<S::S_nDimensions>::dopartition (mesh, comm);
        }

        return true;
    }
};

int main (int argc, char* argv[])
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> comm (new Epetra_SerialComm);
#endif

    bool ilead = (comm->MyPID() == 0);

    // Parse command-line
    if (argc < 5)
    {
        std::cout << "Usage: " << argv[0] << " [mesh_1D_P1] [mesh_2D_Q2] [mesh_3D_P1] [mesh_3D_Q2]" << std::endl;
        return EXIT_SUCCESS;
    }

    // Read the file
    std::vector<std::string> files;
    std::copy (argv + 1, argv + 5, std::back_inserter (files) );

    if (ilead)
    {
        std::cout << "\n [[ LINEAR LINE ]] \n\n";
    }
    tester<LifeV::LinearLine>::test (files[0], comm);

    if (ilead)
    {
        std::cout << "\n [[ QUADRATIC QUAD ]] \n\n";
    }
    tester<LifeV::QuadraticQuad>::test (files[1], comm);

    if (ilead)
    {
        std::cout << "\n [[ LINEAR TETRA ]] \n\n";
    }
    tester<LifeV::LinearTetra>::test (files[2], comm);

    if (ilead)
    {
        std::cout << "\n [[ QUADRATIC HEXA ]] \n\n";
    }
    // Online partitioning of Q2 Hexa not working yet.
    tester<LifeV::QuadraticHexa>::test (files[3], comm, false);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
