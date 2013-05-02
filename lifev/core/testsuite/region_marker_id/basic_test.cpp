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
    @file basic_test.cpp
    @brief Test the consistency of the mesh data structure

    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @contributor
    @maintainer

    @date 2012-09-14

    Colour a mesh with two different colours, count the number
    of elements equal of one of the two.

 */

// ===================================================
//! Includes
// ===================================================
// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#else
#include <lifev/core/filter/ExporterEmpty.hpp>
#endif

using namespace LifeV;

UInt colour_fun ( const Vector3D& barycentre )
{
    if ( barycentre[0] < 0.5 && barycentre[1] < 0.5 )
    {
        return 2;
    }
    return 3;
}

int main (int argc, char* argv[])
{
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif

    typedef RegionMesh<LinearTriangle> mesh_Type;
    typedef boost::shared_ptr< mesh_Type > meshPtr_Type;

    // Create the mesh.
    meshPtr_Type mesh ( new mesh_Type ( comm ) );

    // Fill the mesh with a structured mesh.
    regularMesh2D ( *mesh, 0, 8, 11 );

    // Colour the mesh according to a function.
    MeshUtility::assignRegionMarkerID ( *mesh, colour_fun );

    // Count the number of elements with colour 2
    const UInt colourElements = mesh->elementList().countElementsWithMarkerID ( 2, std::equal_to<markerID_Type>() );

    // Number of elements with colour 2
    const UInt exactNumber = 44;

    {
        // Needed to correctly destroy the exporterHDF5

        // Set the exporter for the mesh region.
#ifdef HAVE_HDF5
        ExporterHDF5< mesh_Type > exporter;
#else
        ExporterEmpty< mesh_Type > exporter;
#endif

        // Set the mesh.
        exporter.setMeshProcId ( mesh, comm->MyPID() );

        // Export the region marker ID.
        exporter.exportRegionMarkerID ( mesh, comm );

        // Do the export.
        exporter.postProcess ( 0 );

    } // Needed to correctly destroy the exporterHDF5

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( colourElements == exactNumber )
    {
        return ( EXIT_SUCCESS );
    }
    else
    {
        return ( EXIT_FAILURE );
    }
}
