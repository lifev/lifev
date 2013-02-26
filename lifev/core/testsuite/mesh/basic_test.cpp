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
    @brief Test the consistency of the mesh data structure

    @author
    @contributor
    @maintainer

    @date 00-00-0000

    Read a 3d mesh and perform internal consistency checks.

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

#include <lifev/core/filter/GetPot.hpp>

#include <lifev/core/mesh/MarkerDefinitions.hpp>
#include <lifev/core/filter/ImporterMesh3D.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshElementBare.hpp>
#include <lifev/core/array/MapEpetra.hpp>

// A dummy class to imitate a VectorEpetra
class dummyVect:
    public std::vector<double>
{
public:
    dummyVect() : std::vector<double>() {}
    dummyVect (unsigned const int& n) : std::vector<double> (n) {}
    dummyVect ( dummyVect v, LifeV::MapEpetraType ) : std::vector<double> (v) {}
    bool isGlobalIDPresent (int ) const
    {
        return true;
    }
    LifeV::MapEpetraType mapType() const
    {
        return LifeV::Repeated;
    }
};

template<typename meshEntity>
class ResetFlag
{
public:
    void operator() (meshEntity& m)
    {
        m.replaceFlag (LifeV::EntityFlags::DEFAULT);
        m.setMarkerID (0);
    }
};

int main (int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::cout << "MPI Initialization" << std::endl;
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> comm ( new Epetra_SerialComm );
#endif


    using namespace LifeV;
    using namespace LifeV::MeshUtility;
    using namespace std;

    GetPot datafile ( "data" );
    string dirname = datafile ( "mesh_dir", "." ); //"../data/mesh/mesh++/";
    string fname = dirname + datafile ( "mesh_file", "cube_47785.m++" ); //dirname+"cube_47785.m++";
    string outfile = "testBuilders.dat";
    ofstream ofile (outfile.c_str() );
    if (ofile.fail() )
    {
        cerr << " Error: Cannot creat output file" << endl;
        abort();
    }

    RegionMesh<LinearTetra> aMesh ( comm );
    typedef RegionMesh<LinearTetra> mesh_Type;
    //    aMesh.test3DBuilder();
    //    aMesh.readMppFile(mystream, id, m);
    ID m = 1;
    readMppFile (aMesh, fname, m);

    cout << " **********************************************************" << endl;

    cout << "  ****************** CHECKING MESH WITH INTERNAL CHECKER" << endl;

    aMesh.check (0, true, true);
    Switch sw;


    cout << " **********************************************************" << endl;

    cout << "  ****************** CLEANING edges and faces " << endl;

    aMesh.edgeList.clear();
    aMesh.faceList.clear();

    checkMesh3D (aMesh, sw,
                 true, true,
                 cerr, cout, ofile);


    aMesh.showMe();


    cout << " Now building local Edges/faces Stuff" << endl << endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    cout << " Now cleaning local Edges/faces Stuff" << endl << endl;
    aMesh.cleanElementRidges();
    aMesh.cleanElementFacets();
    aMesh.showMe();
    cout << " **********************************************************" << endl;

    cout << "  ****************** BUILDING ALL EDGES" << endl;
    UInt bedges_found, iedges_found;
    buildEdges (aMesh, ofile, cerr, bedges_found,
                iedges_found, true, true, true);
    cout << " **********************************************************" << endl;

    cout << "  ****************** BUILDING ALL FACES" << endl;
    UInt bfaces_found, ifaces_found;
    buildFaces (aMesh, ofile, cerr, bfaces_found,
                ifaces_found, true, true, true);

    aMesh.showMe();

    cout << " Now building again local Edges/faces Stuff" << endl << endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    checkMesh3D (aMesh, sw,
                 true, true,
                 cerr, cout, ofile);
    ofile.close();

    ///THIS PART IS ONLY TO VERIFY IF THESE ROUTINES COMPILE PROPERLY
    cerr << "Fixing bpoints" << endl;

    fixBoundaryPoints (aMesh, ofile, cerr, true);
    cerr << "Fixing edge markers" << endl;
    setBoundaryEdgesMarker (aMesh, ofile, cerr, true);
    cerr << "Fixing faces marker" << endl;
    setBoundaryFacesMarker (aMesh, ofile, cerr, true);
    cerr << "Fixing points marker" << endl;
    setBoundaryPointsMarker (aMesh, ofile, cerr, true);
    cerr << endl;
    dummyVect disp (3 * aMesh.numPoints() );
    MeshUtility::MeshTransformer<mesh_Type> transformer (aMesh);
    transformer.moveMesh (disp, 3);
    MeshUtility::MeshStatistics::meshSize sizes = MeshUtility::MeshStatistics::computeSize (aMesh);
    cerr << "Hmin =" << sizes.minH << " Hmax=" << sizes.maxH << std::endl;
    mesh_Type::points_Type newPointList = aMesh.pointList;
    mesh_Type::faces_Type newFaceList = aMesh.faceList;
    Utilities::fixAfterShallowCopy (newFaceList, newPointList);
    for (mesh_Type::faces_Type::iterator i = newFaceList.begin(); i < newFaceList.end(); ++i)
    {
        ID theId = i->point (0).id();
        if (& (i->point (0) ) != & (newPointList[theId]) )
        {
            cerr << "ERROR: Error after renumbering Faces" << std::endl;
            break;
        }
    }
    aMesh.faceList[2].replaceFlag (EntityFlags::CUTTED);
    aMesh.faceList.reorderAccordingToFlag (EntityFlags::CUTTED, &Flag::testOneSet);
    if (aMesh.faceList[0].flag() != EntityFlags::CUTTED)
    {
        cerr << "ERROR: Reordering is not working" << std::endl;
    }
    cout << "Number of cutted faces (should be 1) " <<
         aMesh.faceList.countElementsWithFlag (EntityFlags::CUTTED, &Flag::testOneSet) << std::endl;
    // Reset all flags to default
    aMesh.edgeList.changeAccordingToFunctor (ResetFlag<mesh_Type::edge_Type>() );
    aMesh.edge (0).setMarkerID (10);
    aMesh.edge (5).setMarkerID (10);
    aMesh.edge (10).setMarkerID (15);
    vector<ID> watermarks (2);
    watermarks[0] = 10;
    watermarks[1] = 15;
    // change flags according to marker iD
    SetFlagAccordingToWatermarks changeFlags (EntityFlags::CUTTED, watermarks);
    aMesh.edgeList.changeAccordingToFunctor (changeFlags);
    std::cout << "Number of cutted edges (should be 3) "
              << aMesh.edgeList.countElementsWithFlag (EntityFlags::CUTTED, &Flag::testOneSet) << std::endl;
    aMesh.edgeList.changeAccordingToFunctor ( ResetFlag<mesh_Type::edge_Type>() );
    aMesh.edge (0).setMarkerID (10);
    aMesh.edge (5).setMarkerID (12);
    aMesh.edge (10).setMarkerID (15);
    SetFlagAccordingToMarkerRanges changer ( &Flag::turnOn ); //I may use the default constructor
    changer.insert (std::make_pair (10, 12), EntityFlags::INTERNAL_INTERFACE);
    changer.insert (std::make_pair (15, 18), EntityFlags::CUTTED);
    aMesh.edgeList.changeAccordingToFunctor (changer);
    cout << "Number of cutted edges (should be 1) " <<
         aMesh.edgeList.countElementsWithFlag (EntityFlags::CUTTED, &Flag::testOneSet) << std::endl;
    cout << "Number of internal interface edges (should be 2) " <<
         aMesh.edgeList.countElementsWithFlag (EntityFlags::INTERNAL_INTERFACE, &Flag::testOneSet) << std::endl;

    SetFlagAccordingToWatermark<std::equal_to<markerID_Type> > changer2 (EntityFlags::INTERNAL_INTERFACE, 12000, Flag::turnOn);
    aMesh.faceList.changeAccordingToFunctor (changer2);

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return ( EXIT_SUCCESS );
}

