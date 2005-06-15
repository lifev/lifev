/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#ifndef __MESH_UTILITIES__
#define __MESH_UTILITIES__

#include <boost/shared_ptr.hpp>

#include <life/lifemesh/mesh_util_base.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifefem/currentFE.hpp>
#include <life/lifefem/currentBdFE.hpp>

//! \file mesh_util.h
//! \file mesh_util.h
/*!
  This file contains a set of functions to be used to test a 3D mesh and fix
  some possible problems. Some methods are general (2D and 3D meshes), some are
  specific to 3D meshes.

  They sametimes
  require in input  a Switch parementer sw and ostreams references.
  The switch valueas are
<UL>
<LI>"ABORT_CONDITION"  "SKIPPED_ORIENTATION_TEST"</LI>
<LI>"HAS_NEGATIVE_VOLUMES" "BFACE_STORED_MISMATCH"</LI>
<LI>"BELEMENT_COUNTER_UNSET" "BFACE_COUNTER_MISMATCH"</LI>
<LI>"BFACE_MISSING" "FIXED_MAX_NUM_FACES"</LI>
<LI>"FIXED_FACE_COUNTER" "NUM_FACES_MISMATCH"</LI>
<LI>"FIXED_MAX_NUM_EDGES" "FIXED_EDGES"</LI>
<LI>"NOT_HAS_POINTS" "FIXED_POINTS_ID"</LI>
<LI>"POINTS_MARKER_UNSET" "NOT_HAS_VOLUMES"</LI>
<LI>"FIXED_VOLUMES_ID" "VOLUMES_MARKER_UNSET"</LI>
<LI>"FIXED_VOLUME_COUNTER" "BUILD_BFACES"</LI>
<LI>"FIXED_BFACES_FIRST" "FIXED_FACES_ID"</LI>
<LI>"NOT_HAS_FACES" "FIXED_BFACE_COUNTER"</LI>
<LI>"FACE_MARKER_UNSET" "FACE_MARKER_FIXED"</LI>
<LI>"FIXED_FACE_COUNTER" "NOT_HAS_EDGES"</LI>
<LI>"BUILD_BEDGES" "FIXED_BEDGES_FIRST" </LI>
<LI>"FIXED_EDGES_ID" "EDGE_MARKER_UNSET"</LI>
<LI>"EDGE_MARKER_FIXED" "FIXED_BEDGES_COUNTER"</LI>
<LI>"DOMAIN_NOT_CLOSED" "FIXED_BOUNDARY_POINTS"</LI>
<LI>"POINT_MARKER_UNSET" "FIXED_POINT_MARKER"</LI>
<LI>"FIXED_BPOINTS_COUNTER" "NOT_EULER_OK"</LI>
</UL>

\pre All functions contained in this file require as precondition a mesh
with the points and volumes connectivity set. Some functions have also
other preconditions, which will be then specified in the function
documentation.
*/
namespace LifeV
{
/*
*****************************************************************************
                                GEOMETRY TESTS
*****************************************************************************
*/
//!  \brief Report 3D element orientation
/*!  It uses a linear representation of the Tetra/Hexa: it is only a
  orientation check.  The orientation is considered positive if it
  obeys the right-hand rule (right-hand orientation).

  \param mesh A region mesh of 3D elements

  \param elSign A vector of bool: true means positive orientation.

  \param sw The switch used to communicate test results. This function may
  set the conditions
  <TT>HAS_NEGATIVE_VOLUMES,SKIP_ORIENTATION_TEST</TT>

  The first reports that some mesh elements have negative volume, the
  second that the test has been skipped because has not yet beem
  implemented for the element under consideration.

  \return It returns the mesh computed volume.
*/
template <typename RegionMesh3D>
Real checkVolumes( RegionMesh3D const & mesh,
                   SimpleVect<bool> & elSign,
                   Switch & sw )
{
    Real meas = 0.0;
    Real lmeas = 0.0;
    elSign.clear();
    elSign.reserve( mesh.numVolumes() );
    typedef typename RegionMesh3D::VolumeShape GeoShape;

    switch ( GeoShape::Shape )
    {
    case TETRA:
        {
            CurrentFE fe( feTetraP1, geoLinearTetra, quadRuleTetra1pt );
            for ( ID i = 1; i <= mesh.numVolumes(); i++ )
            {
                fe.updateJac( mesh.volume( i ) );
                lmeas = fe.measure();
                meas += lmeas;
                elSign.push_back( lmeas > 0.0 );
            }
        }
        break;
    case HEXA:
        {
            CurrentFE fe( feHexaQ1, geoBilinearHexa, quadRuleHexa1pt );
            for ( ID i = 1; i <= mesh.numVolumes(); i++ )
            {
                fe.updateJac( mesh.volume( i ) );
                lmeas = fe.measure();
                meas += lmeas;
                elSign.push_back( lmeas > 0.0 );
            }
        }
        break;
    default:
        sw.create( "SKIP_ORIENTATION_TEST", true );

        return 0;
    }

    if ( std::find( elSign.begin(), elSign.end(), false ) != elSign.end() )
        sw.create( "HAS_NEGATIVE_VOLUMES", true );

    return meas;
}

/*!
  \brief Fixes  negative volume elements.

  Given a std::vector<bool> indicating negative elements, it inverts those that
  have been found negative.

  \param mesh A 3D mesh. It will be modified.

  \param elSign a vector of bools. The value false correspond to the elements that have to be swapped. It is created by
  checkVolumes().

  \post A mesh with all volumes with positive oreintation.
*/
template <typename RegionMesh3D>
void fixVolumes( RegionMesh3D & mesh,
                 const SimpleVect<bool> & elSign,
                 Switch & sw )
{
    typedef typename RegionMesh3D::VolumeShape GeoShape;

    static const ID otn_Tetra[ 10 ] =
        {
            2, 1, 3, 4, 5, 7, 6, 9, 8, 10
        };
    static const ID otn_Hexa[ 20 ] =
        {
            1, 4, 3, 2, 5, 8, 7, 6, 12, 11,
            10, 9, 13, 16, 15, 14, 20, 19, 18, 17
        };
    for ( ID i = 1; i <= mesh.numVolumes(); i++ )
    {

        if ( ! elSign( i ) )
        {
            switch ( GeoShape::Shape )
            {
            case TETRA:
                mesh.volume( i ).exchangePoints( otn_Tetra );
                break;
            case HEXA:
                mesh.volume( i ).exchangePoints( otn_Hexa );
                break;
            default:
                sw.create( "ABORT_CONDITION", true );
                return ;
            }
        }
    }
}
//!\brief Computes volume enclosed by boundary faces
/*!
  It computes, for $i=1,2,3$, the integral \f$\int_{\partial \Omega} x_i n_i
  d\gamma \f$, \f$n_i\f$ being the i-th component of the boundary normal. If
  the domain boundary is properly disretised they should all return (within
  discretisation and truncation errors) the quantity \f$\vert\Omega\vert\f$.

  \warning Not to be used for accurate computations (it always adopts
  linear or bilinear elements, with a simple integration rule) \param mesh
  A 3D mesh \param vols returns 3 Real corresponding to the 3 integrals
*/
template <typename RegionMesh3D>
void getVolumeFromFaces( RegionMesh3D const & mesh,
                         Real vols[ 3 ],
                         std::ostream & err = std::cerr )
{
    GetCoordComponent getx( 0 );
    GetCoordComponent gety( 1 );
    GetCoordComponent getz( 2 );
    vols[ 0 ] = 0.0;
    vols[ 1 ] = 0.0;
    vols[ 2 ] = 0.0;
    typedef typename RegionMesh3D::FaceShape GeoBShape;
    typedef typename RegionMesh3D::FaceType FaceType;
    typedef boost::shared_ptr<CurrentBdFE> current_fe_type;

    current_fe_type bdfe;

    switch ( GeoBShape::Shape )
    {
    case TRIANGLE:
        bdfe = current_fe_type( new CurrentBdFE( feTriaP1, geoLinearTria,
                                                 quadRuleTria1pt ) );
        for ( ID i = 1; i <= mesh.numBFaces(); i++ )
        {
            bdfe->updateMeasNormal( mesh.face( i ) );
            vols[ 0 ] += bdfe->integral_n( getx );
            vols[ 1 ] += bdfe->integral_n( gety );
            vols[ 2 ] += bdfe->integral_n( getz );
        }
        break;
    case QUAD:
        bdfe = current_fe_type( new CurrentBdFE( feQuadQ1, geoBilinearQuad,
                                                 quadRuleQuad1pt ) );
        for ( ID i = 1; i <= mesh.numBFaces(); i++ )
        {
            bdfe->updateMeasNormal( mesh.face( i ) );
            vols[ 0 ] += bdfe->integral_n( getx );
            vols[ 1 ] += bdfe->integral_n( gety );
            vols[ 2 ] += bdfe->integral_n( getz );
        }
        break;
    default:
        err << "Only tria and quad surface elements  may be checked for volume orientation at the moment" << std::endl;
        ASSERT0( false, "ABORT CONDITION OCCURRED" );
    }
}

//! Tests if the surface of the mesh is closed by computing surface integrals.
/*! It computes \f$\sum_{i=1}^3\int_{\partial \Omega} n_i d\gamma\f$.
  The value returned  should be very proximal to zero
 */
template <typename RegionMesh3D>
Real testClosedDomain( RegionMesh3D const & mesh,
                       std::ostream & err = std::cerr )
{
    typedef typename RegionMesh3D::FaceType FaceType;

    typedef boost::shared_ptr<CurrentBdFE> current_fe_type;
    current_fe_type bdfe;

    GetOnes ones;
    Real test( 0.0 );

    switch ( RegionMesh3D::FaceShape::Shape )
    {
    case TRIANGLE:
        bdfe = current_fe_type( new CurrentBdFE( feTriaP1, geoLinearTria,
                                                 quadRuleTria1pt ) );
        for ( ID i = 1; i <= mesh.numBFaces(); i++ )
        {
            bdfe->updateMeasNormal( mesh.face( i ) );
            test += bdfe->integral_n( ones );
        }
        break;
    case QUAD:
        bdfe = current_fe_type( new CurrentBdFE( feQuadQ1, geoBilinearQuad,
                                                 quadRuleQuad1pt ) );
        for ( ID i = 1; i <= mesh.numBFaces(); i++ )
        {
            bdfe->updateMeasNormal( mesh.face( i ) );
            test += bdfe->integral_n( ones );
        }

        break;

    default:
        err << "Only tria and quad surface elements  may be checked for volume orientation at the moment" << std::endl;
        ASSERT0( false, "ABORT CONDITION OCCURRED" );
    }

    return test;

}

/*
*****************************************************************************
                UTILITIES FOR GLOBAL CHECKINGS
*****************************************************************************
*/

//! This function performs  a lot of checks.

/*!
  The name is inappropriate since the tests that are performed are not just on
  the topological structure of the mesh. The output is directed to three
  output streams:
  <ul>
  <li> out-> usually standard output: important informative messages
  <li> err-> usually standard error: error messages
  <li> clog-> usually a file stream: informative messages which may be rather verbose.
  </ul>
 Furthermore, ths Switch sw (see switch.h) will return a set of
 keywords useful for other possible actions. If fix=true, this routines
 performes the steps needed to get an acceptable mesh, otherwise the input
 mesh is not modified .
*/

template <typename RegionMesh3D>
bool checkMesh3D( RegionMesh3D & mesh,
                  Switch & sw,
                  bool fix = true,
                  bool verbose = false,
                  std::ostream & out = std::cerr,
                  std::ostream & err = std::cerr,
                  std::ostream & clog = std::cout )
{

    if ( mesh.storedPoints() == 0 )
    {
        err << "FATAL: mesh does not store points: I cannot do anything"
            << std::endl;
        sw.create( "ABORT_CONDITION", true );
        sw.create( "NOT_HAS_POINTS", true );
        return false;
    }

    if ( !checkIdnumber( mesh.pointList ) )
    {
        err << "ERROR: points ids where wrongly set" << std::endl;
        err << "FIXED" << std::endl;
        if ( fix )
            sw.create( "FIXED_POINTS_ID", true );
        if ( fix )
            fixIdnumber( mesh.pointList );
    }

    if ( !checkMarkerSet( mesh.pointList ) )
    {
        err << "WARNING: Not all points have marker flag set" << std::endl;
        sw.create( "POINTS_MARKER_UNSET", true );
    }


    //-------------------------------------------------------------------------
    //                                    VOLUMES
    //-------------------------------------------------------------------------

    if ( mesh.storedVolumes() == 0 )
    {
        err << "FATAL: mesh does not store volumes: I cannot do anything"
            << std::endl;
        sw.create( "ABORT_CONDITION", true );
        sw.create( "NOT_HAS_VOLUMES", true );
        return false;
    }

    if ( !checkIdnumber( mesh.volumeList ) )
    {
        err << "ERROR: volume ids where wrongly set" << std::endl;
        err << "FIXED" << std::endl;
        if ( fix )
            sw.create( "FIXED_VOLUMES_ID", true );
        if ( fix )
            fixIdnumber( mesh.volumeList );
    }

    if ( !checkMarkerSet( mesh.volumeList ) )
    {
        err << "WARNING: Not all volumes have marker flag set" << std::endl;
        sw.create( "VOLUMES_MARKER_UNSET", true );
        if ( fix )
        {
            for ( typename RegionMesh3D::Volumes::iterator iv = mesh.volumeList.begin();
                    iv != mesh.volumeList.end();++iv )
            {
                if ( iv->isMarkerUnset() )
                    iv->setMarker( mesh.marker() );
            }
        }
    }

    if ( mesh.numElements() < mesh.storedVolumes() )
    {
        err << "WARNING: Mesh Volumes must be at least "
            << mesh.storedVolumes() << std::endl;
        if ( fix )
            mesh.numVolumes() = mesh.storedVolumes();
        if ( fix )
            sw.create( "FIXED_VOLUME_COUNTER", true );
    }

    // test now orientation

    boost::shared_ptr<SimpleVect<bool> > elSign( new SimpleVect<bool> );

    Real meshMeasure = checkVolumes( mesh, *elSign, sw );
    UInt positive;

    if ( sw.test( "SKIP_ORIENTATION_TEST" ) )
    {
        clog << "W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" << std::endl;
        err << "W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" << std::endl;
    }
    else if ( sw.test( "HAS_NEGATIVE_VOLUMES" ) )
    {
        positive = count( elSign->begin(), elSign->end(), true );
        clog << positive << "W: positive elements out of"
             << mesh.storedVolumes() << std::endl;
        if ( fix )
            clog << "Fixing negative elements" << std::endl;
        if ( fix )
            fixVolumes( mesh, *elSign, sw );
        if ( sw.test( "ABORT_CONDITION" ) )
        {
            err << "ABORT: Cannot fix volumes, this element is not supported"
                << std::endl;
            return false;
        }
        else
        {
            sw.unset( "HAS_NEGATIVE_VOLUMES" );
            meshMeasure = checkVolumes( mesh, *elSign, sw );
            if ( sw.test( "HAS_NEGATIVE_VOLUMES" ) )
            {
                if ( fix )
                    err << "ABORT: Cannot fix volumes: something wrong with this mesh" << std::endl;
                if ( fix )
                    sw.create( "ABORT_CONDITION", true );
                return false;
            }
        }
    }

    clog << "Volume enclosed by the mesh= " << meshMeasure << std::endl
    << "(Computed by integrating mesh elements measures)" << std::endl
    << "(Using 1 point Quadrature rule)" << std::endl;

    //-----------------------------------------------------
    //                                    BOUNDARY FACES
    //-----------------------------------------------------

    boost::shared_ptr<TempFaceContainer> bfaces(  new TempFaceContainer );
    UInt numInternalFaces, numFaces;

    UInt bFacesFound = findBoundaryFaces( mesh, *bfaces, numInternalFaces );

    numFaces = bFacesFound + numInternalFaces;

    EnquireBFace<RegionMesh3D> enquireBFace( mesh, *bfaces );


    if ( mesh.storedFaces() == 0 ||
         mesh.numBElements() > mesh.storedFaces() ||
         bFacesFound > mesh.storedFaces() )
    {
        err << "ERROR: Not all boundary faces stored" << std::endl;
        if ( fix )
            sw.create( "BUILD_BFACES", true );
        if ( fix )
            buildFaces( mesh, clog, err, bFacesFound, numInternalFaces,
                        true, false, verbose, bfaces.get() );
    }
    else
    {
        //
        // Make sure BFaces are stored first
        // Here I need to use a method that does not require the proper
        // setting of boundary Points!
      if ( mesh.numBFaces() !=  bFacesFound) {
	err<<" ERROR: Number of B faces does not correspond to real one"<<std::endl;
	if (fix)
	  {
	    err<<"FIXED Number of B faces has been fixed to:" << bFacesFound<<std::endl;
	    mesh.setNumBFaces( bFacesFound);
	  }
      }
      
      if ( mesh.numBFaces() < mesh.storedFaces() )
        {
	  if ( fix )
	    std::stable_partition( mesh.faceList.begin(),
				   mesh.faceList.end(), enquireBFace );
	  if ( fix )
	    fixIdnumber( mesh.faceList );
	  if ( fix )
	    sw.create( "FIXED_BFACES_FIRST" );
        }
      
      
        if ( !checkIdnumber( mesh.faceList ) )
        {
	  err << "ERROR: face ids where wrongly set" << std::endl;
	  err << "FIXED" << std::endl;
	  if ( fix )
	    sw.create( "FIXED_FACES_ID", true );
	  if ( fix )
	    fixIdnumber( mesh.faceList );
        }

        // Check Consistency with the mesh. Beware that this method changes *bfaces!

        if ( fix )
            fixBoundaryFaces( mesh, clog, err, sw, numFaces, bFacesFound,
                              false, verbose, bfaces.get() );

        if ( mesh.storedFaces() == 0 )
        {
            err << "ABORT CONDITION: cannot find boundary faces" << std::endl;
            sw.create( "NOT_HAS_FACES", true );
            sw.create( "ABORT_CONDITION", true );
        }


        if ( mesh.numBFaces() == 0 )
        {
            err << " MeshBFaces counter is unset" << std::endl;
            if ( fix )
                mesh.setNumBFaces( mesh.storedFaces() );
            if ( fix )
                sw.create( "FIXED_BFACE_COUNTER", true );
            if ( fix )
                mesh.setLinkSwitch( "HAS_BOUNDARY_FACES" );
        }


        if ( !checkMarkerSet( mesh.faceList ) )
        {
            err << "WARNING: Not all faces have marker flag set" << std::endl;
            sw.create( "FACE_MARKER_UNSET", true );
            if ( fix )
                setBFacesMarker( mesh, clog, err, verbose );
            if ( fix && checkMarkerSet( mesh.faceList ) )
            {
                sw.create( "FACE_MARKER_UNSET", false );
                sw.create( "FACE_MARKER_FIXED", true );
            }
        }
    }


    if ( mesh.numFaces() != bFacesFound + numInternalFaces )
    {
        err << "WARNING Number of faces incorrectly set" << std::endl;
        err << "        It was       " << mesh.numFaces() << std::endl;
        err << "        It should be " << bFacesFound + numInternalFaces
            << std::endl;
        if ( fix )
            err << "        Fixing" << std::endl;
        if ( fix )
            mesh.numFaces() = bFacesFound + numInternalFaces;
        if ( fix )
            sw.create( "FIXED_FACE_COUNTER", true );
    }

    if ( fix && mesh.storedFaces() == bFacesFound + numInternalFaces)
        mesh.setLinkSwitch( "HAS_ALL_FACES" );

    out << " Boundary faces found:" << bFacesFound << std::endl;
    out << " Num Faces Stored stored:" << mesh.storedFaces() << std::endl;
    out << " Num faces in mesh      :" << mesh.numFaces() << std::endl;
    out << " Boundary faces counter gives:" << mesh.numBFaces() << std::endl;

    //-----------------------------------------------------
    //                                    BOUNDARY EDGES
    //-----------------------------------------------------

    boost::shared_ptr<TempEdgeContainer> bedges( new TempEdgeContainer );

    UInt bEdgesFound = findBoundaryEdges( mesh, *bedges );
    EnquireBEdge<RegionMesh3D> enquireBEdge( mesh, *bedges );

    UInt intedge;
    UInt Ned = 0;
    TempEdgeContainer iedges;

    if ( mesh.storedEdges() == 0 ||
         mesh.numBEdges() > mesh.storedEdges() ||
         bEdgesFound > mesh.storedEdges() )
    {
        err << "WARNING: mesh does not store (all) boundary edges" << std::endl;
        sw.create( "NOT_HAS_EDGES", true );
        if ( fix )
	  buildEdges( mesh, clog, err, bEdgesFound, intedge, true, false,
		      false, bedges.get() );
        Ned = bEdgesFound + intedge;
        if ( fix )
	  sw.create( "BUILD_BEDGES", true );
    }
    else
    {

      // Make sure BEdges are first
      // Here I need to use a method that does not require the proper
      // setting of boundary Points!
      if ( mesh.numBEdges() !=  bEdgesFound) {
	err<<" ERROR: Number of BEdges does not correspond to real one"<<std::endl;
	if (fix)
	  {
	    err<<"FIXED Number of BEdges has been fixed to:" <<bEdgesFound<<std::endl;
	    mesh.setNumBEdges( bEdgesFound);
	  }
      }
      
        if ( fix )
            std::stable_partition( mesh.edgeList.begin(), mesh.edgeList.end(),
                                   enquireBEdge );
        if ( fix )
            fixIdnumber( mesh.edgeList );
        if ( fix )
            sw.create( "FIXED_BEDGES_FIRST" );

        if ( !checkIdnumber( mesh.edgeList ) )
        {
            err << "ERROR: edge ids where wrongly set" << std::endl;
            err << "FIXED" << std::endl;
            sw.create( "FIXED_EDGES_ID", true );
            fixIdnumber( mesh.edgeList );
        }


        if ( !checkMarkerSet( mesh.edgeList ) )
        {
            err << "WARNING: Not all edges have marker flag set" << std::endl;
            sw.create( "EDGE_MARKER_UNSET", true );
            if ( fix )
                setBEdgesMarker( mesh, clog, err, verbose );
            if ( fix && checkMarkerSet( mesh.edgeList ) )
            {
                sw.unset( "EDGE_MARKER_UNSET" );
                sw.create( "EDGE_MARKER_FIXED", true );
            }
        }
        out << "Computing internal edges";
        if ( fix )
            Ned = bEdgesFound + findInternalEdges( mesh, *bedges, iedges );
    }
    iedges.clear();
    TempEdgeContainer tmp;
    iedges.swap(tmp);

    if ( mesh.numBEdges() != bEdgesFound )
    {
        err << "WARNING: number of found boundary edges:" << bEdgesFound
            << std::endl
            << " does not match that declared in mesh, i.e. "
            << mesh.numBEdges() << std::endl;
        if ( mesh.numBEdges() == 0 )
        {
            if ( fix )
                err << "FIXING" << std::endl;
            if ( fix )
                sw.create( "FIXED_BEDGES_COUNTER", true );
            if ( fix )
                mesh.setNumBEdges( bEdgesFound );
        }
    }

    if ( Ned != mesh.numEdges() )
    {
        if ( fix )
            err << "WARNING: Counter of number of edges badly set: Should be (actual number)" << Ned << std::endl;
        err << "It is instead equal to " << mesh.numEdges();
        if ( fix )
        {
            err << " **FIXED" << std::endl;
            mesh.numEdges() = Ned;
        }
        std::cerr << std::endl;
    }
    UInt nbed;
    UInt counte = testClosedDomain_Top( mesh, nbed );
    if ( counte == 0 )
    {
        out << "**DOMAIN SURFACE IS (TOPOLOGICALLY) CLOSED" << std::endl;
    }
    else
    {
        sw.create( "DOMAIN_NOT_CLOSED", true );
        err << "WARNING: DOMAIN APPEARS TO HAVE AN OPEN BOUNDARY (TOPOLOGY CHECK)" << std::endl;
        err << "Number of inconsistent edges:" << counte << std::endl;
    }



    //-----------------------------------------------------
    //                                    POINTS
    //-----------------------------------------------------
    // Now that boundary faces have been correctly set we may work out
    // boundaty points
    if (fix) fixBPoints(mesh,clog,err,verbose);
    EnquireBPoint<RegionMesh3D> enquirebpoint( mesh );
    UInt foundBPoints = std::count_if( mesh.pointList.begin(),
                                       mesh.pointList.end(), enquirebpoint );

    if ( foundBPoints == 0 || foundBPoints < mesh.storedBPoints() )
    {
        err << "WARNING Bpoints indicator not correctly set" << std::endl;
        if ( fix )
            err << "FIXING by recomputing from boundary faces" << std::endl;
        fixBPoints( mesh, clog, err, verbose );
        if ( fix )
            foundBPoints = std::count_if( mesh.pointList.begin(),
                                          mesh.pointList.end(),
                                          enquirebpoint );
        if ( fix )
            sw.create( "FIXED_BOUNDARY_POINTS", true );
    }

    if ( ! checkMarkerSet( mesh.pointList ) )
    {
        err << "WARNING B. Points MARKER incorrectly set" << std::endl;
        if ( fix )
        {
            setBPointsMarker( mesh, clog, std::cerr, verbose );
            if ( ! checkMarkerSet( mesh.pointList ) )
            {
                err << "Cannot Fix Points MARKER" << std::endl;
                sw.create( "POINT_MARKER_UNSET", true );
            }
            else
            {
                err << "FIXED" << std::endl;
                sw.create( "FIXED_POINT_MARKER", true );
            }
        }
    }


    if ( mesh.storedBPoints() == 0 )
    {
        err << "WARNING B. Points COUNTER incorrectly set" << std::endl;
        if ( fix )
            setBPointsCounters( mesh ) ;
        if ( fix )
            err << " FIXED" << std::endl;
        if ( fix )
            sw.create( "FIXED_BPOINTS_COUNTER", true );
    }

    if ( mesh.numPoints() == 0 )
    {
        err << "WARNING Points Counter unset" << std::endl;
        if ( fix )
            mesh.numPoints() = mesh.storedPoints();
        if ( fix )
            sw.create( "FIXED_POINTS_COUNTER", true );
    }

    //-----------------------------------------------------
    //                                   FINAL CHECKS
    //-----------------------------------------------------
    out << " ********     COUNTERS CONTENT **********************************" << std::endl;

    out << " Num Volumes:  " << mesh.numVolumes() << std::endl;
    out << " Num Vertices: " << mesh.numVertices()
        << " Num B. Vertices: " << mesh.numBVertices() << std::endl;
    out << " Num Points:   " << mesh.numPoints()
        << " Num B. Points  : " << mesh.numBPoints() << std::endl;
    out << " Num Edges:    " << mesh.numEdges()
        << " Num B. Edges   : " << mesh.numBEdges() << std::endl;
    out << " Num Faces:    " << mesh.numFaces()
        << " Num B. Faces   : " << mesh.numBFaces() << std::endl;
    out << " ********     END COUNTERS **********************************"
        << std::endl;

    bool eulok1 = ( 2 * mesh.numFaces() -
                    mesh.numLocalFaces() * mesh.numVolumes() -
                    mesh.numBFaces() ) == 0;

    bool eulok2( true );

    if ( RegionMesh3D::ElementShape::Shape == TETRA )
    {
        out << std::endl << "*** CHECKING WITH EULER FORMULAE *****  "
            << std::endl;
        eulok2 = ( mesh.numEdges() -
                   mesh.numVolumes() -
                   mesh.numVertices() -
                   ( 3 * mesh.numBFaces() -
                     2 * mesh.numBVertices() ) / 4 ) == 0;
    }

    if ( !( eulok1 && eulok2 ) )
    {
        err << "WARNING: The following Euler formula(s) are not satisfied"
            << std::endl;
        sw.create( "NOT_EULER_OK" );
    }
    else
    {
        out << std::endl << "*** CHECKING IS OK *****  " << std::endl;
    }

    if ( !eulok1 )
    {
        err << "  2*nFaces = nFacesPerVolume*nVolumes + nBoundaryFaces"
            << std::endl;
        err << "  2*" << mesh.numFaces() << " != " << mesh.numLocalFaces()
            << " * "<< mesh.numVolumes() << " + " << mesh.numBFaces()
            << std::endl;
    }

    if ( !eulok2 )
    {
        err << "  nEdges = nVolumes + nVertices + (3*nBoundaryFaces - 2*nBoundaryVertices)/4" << std::endl;
        err << "  " << mesh.numEdges() << " != " << mesh.numVolumes() << " + "
            << mesh.numVertices() << " + (3*" << mesh.numBFaces() << " - 2*"
            << mesh.numBVertices() << ")/4" << std::endl;
    }

    mesh.setLinkSwitch( "HAS_BEEN_CHECKED" );

    return true;
}
}
#endif
