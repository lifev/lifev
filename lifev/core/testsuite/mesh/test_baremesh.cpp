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
    @brief Test the consistency of the baremesh data structure

    @author Antonio Cervone <ant.cervone@gmail.com>
    @contributor
    @maintainer Antonio Cervone <ant.cervone@gmail.com>

    @date 13-09-2011

    Read a 3d mesh via the baremesh structure and perform internal consistency checks.

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

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/BareMesh.hpp>

#include <lifev/core/filter/ParserINRIAMesh.hpp>
#include <lifev/core/mesh/ConvertBareMesh.hpp>
#include <lifev/core/array/MapEpetra.hpp>

// A dummy class to imitate a VectorEpetra
class dummyVect:
    public std::vector<double>
    {
    public:
        dummyVect(): std::vector<double>(){}
        dummyVect(unsigned const int & n): std::vector<double>(n){}
        dummyVect( dummyVect v, LifeV::MapEpetraType ):std::vector<double>(v){}
        bool isGlobalIDPresent(int )const {return true;}
        LifeV::MapEpetraType mapType() const { return LifeV::Repeated; }
    };

template<typename meshEntity>
class ResetFlag
{
public:
    void operator()(meshEntity &m){
        m.replaceFlag(LifeV::EntityFlags::DEFAULT);
        m.setMarkerID(0);
    }
};

using namespace LifeV;

int main(int argc, char** argv)
{
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    std::cout << "MPI Initialization" << std::endl;
    boost::shared_ptr<Epetra_Comm> comm( new Epetra_MpiComm( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> comm( new Epetra_SerialComm );
#endif

    GetPot datafile( "data_baremesh" );
    std::string dirname = datafile( "mesh_dir", "." ); // "../data/mesh/inria/"
    std::string fname = dirname + datafile( "mesh_file",
                                          "cartesian_cube8.mesh" ); // dirname + "cartesian_cube8.mesh";
    std::string outfile = "testBuilders_baremesh.dat";
    std::ofstream ofile( outfile.c_str() );
    if ( ofile.fail() ) {
      std::cerr << " Error: Cannot creat output file " << std::endl;
      abort();
    }

    typedef BareMesh<LinearTetra> bareMesh_Type;
    typedef RegionMesh<LinearTetra>     mesh_Type;

    ID m=1;

    bareMesh_Type bareMesh;
    MeshIO::ReadINRIAMeshFile( bareMesh, fname, m );

    mesh_Type aMesh( comm );
    convertBareMesh ( bareMesh, aMesh );

    std::cout << " **********************************************************" << std::endl;

    std::cout << "  ****************** CHECKING MESH WITH INTERNAL CHECKER" << std::endl;

    aMesh.check( 0, true, true );
    Switch sw;

    std::cout << " **********************************************************" << std::endl;

    std::cout << "  ****************** CLEANING edges and faces " << std::endl;

    aMesh.edgeList.clear();
    aMesh.faceList.clear();

    checkMesh3D(aMesh, sw,
                true, true,
                std::cerr, std::cout, ofile);


    aMesh.showMe();


    std::cout << " Now building local Edges/faces Stuff" << std::endl << std::endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    std::cout << " Now cleaning local Edges/faces Stuff" << std::endl << std::endl;
    aMesh.cleanElementRidges();
    aMesh.cleanElementFaces();
    aMesh.showMe();
    std::cout << " **********************************************************" << std::endl;

    std::cout << "  ****************** BUILDING ALL EDGES" << std::endl;
    UInt bedges_found, iedges_found;
    MeshUtility::buildEdges( aMesh, ofile, cerr, bedges_found,
               iedges_found, true, true, true );
    std::cout << " **********************************************************" << std::endl;

    std::cout << "  ****************** BUILDING ALL FACES" << std::endl;
    UInt bfaces_found, ifaces_found;
    MeshUtility::buildFaces( aMesh, ofile, cerr, bfaces_found,
               ifaces_found, true, true, true );

    aMesh.showMe();

    std::cout << " Now building again local Edges/faces Stuff" << std::endl << std::endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    checkMesh3D(aMesh, sw,
                true, true,
                cerr, cout, ofile);
    ofile.close();

    ///THIS PART IS ONLY TO VERIFY IF THESE ROUTINES COMPILE PROPERLY
    std::cerr << "Fixing bpoints" << std::endl;

    MeshUtility::fixBoundaryPoints( aMesh, ofile, std::cerr, true );
    std::cerr << "Fixing edge markers" << std::endl;
    MeshUtility::setBoundaryEdgesMarker( aMesh, ofile, std::cerr, true );
    std::cerr << "Fixing faces marker" << std::endl;
    MeshUtility::setBoundaryFacesMarker( aMesh, ofile, std::cerr, true );
    std::cerr << "Fixing points marker" << std::endl;
    MeshUtility::setBoundaryPointsMarker( aMesh, ofile, std::cerr, true );
    std::cerr << std::endl;
    dummyVect disp( 3*aMesh.numPoints() );
    MeshUtility::MeshTransformer<mesh_Type > transformer( aMesh );
    transformer.moveMesh( disp, 3 );
    MeshUtility::MeshStatistics::meshSize sizes= MeshUtility::MeshStatistics::computeSize( aMesh );
    std::cerr << "Hmin =" << sizes.minH << " Hmax=" << sizes.maxH << std::endl;
    mesh_Type::points_Type newPointList = aMesh.pointList;
    mesh_Type::faces_Type newFaceList = aMesh.faceList;
    Utilities::fixAfterShallowCopy( newFaceList, newPointList );
    for ( mesh_Type::faces_Type::iterator i=newFaceList.begin(); i<newFaceList.end(); ++i )
    {
        ID theId = i->point(0).id();
        if( &(i->point(0) ) != &( newPointList[ theId ] ) )
        {
            std::cerr << "ERROR: Error after renumbering Faces" << std::endl;
            break;
        }
    }
    aMesh.faceList[2].replaceFlag( EntityFlags::CUTTED );
    aMesh.faceList.reorderAccordingToFlag( EntityFlags::CUTTED, Flag::testOneSet );
    if ( aMesh.faceList[0].flag() != EntityFlags::CUTTED )
    {
        std::cerr << "ERROR: Reordering is not working" << std::endl;
    }
    std::cout << "Number of cutted faces (should be 1) " <<
                    aMesh.faceList.countElementsWithFlag( EntityFlags::CUTTED,Flag::testOneSet ) << std::endl;
    // Reset all flags to default
    aMesh.edgeList.changeAccordingToFunctor( ResetFlag<mesh_Type::edge_Type>() );
    aMesh.edge(  0 ).setMarkerID( 10 );
    aMesh.edge(  5 ).setMarkerID( 10 );
    aMesh.edge( 10 ).setMarkerID( 15 );
    std::vector<ID> watermarks( 2 );
    watermarks[ 0 ] = 10;
    watermarks[ 1 ] = 15;
    // change flags according to marker iD
    SetFlagAccordingToWatermarks  changeFlags( EntityFlags::CUTTED, watermarks );
    aMesh.edgeList.changeAccordingToFunctor( changeFlags );
    std::cout << "Number of cutted edges (should be 3) " <<
                     aMesh.edgeList.countElementsWithFlag(EntityFlags::CUTTED,Flag::testOneSet) << std::endl;
    aMesh.edgeList.changeAccordingToFunctor( ResetFlag<mesh_Type::edge_Type>() );
    aMesh.edge(  0 ).setMarkerID( 10 );
    aMesh.edge(  5 ).setMarkerID( 12 );
    aMesh.edge( 10 ).setMarkerID( 15 );
    SetFlagAccordingToMarkerRanges changer( Flag::turnOn ); //I may use the default constructor
    changer.insert( std::make_pair( 10, 12 ),EntityFlags::INTERNAL_INTERFACE );
    changer.insert( std::make_pair( 15, 18 ),EntityFlags::CUTTED );
    aMesh.edgeList.changeAccordingToFunctor( changer );
    std::cout << "Number of cutted edges (should be 1) " <<
                     aMesh.edgeList.countElementsWithFlag( EntityFlags::CUTTED, Flag::testOneSet ) << std::endl;
    std::cout << "Number of internal interface edges (should be 2) " <<
                      aMesh.edgeList.countElementsWithFlag( EntityFlags::INTERNAL_INTERFACE, Flag::testOneSet ) << std::endl;

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return( EXIT_SUCCESS );
}

