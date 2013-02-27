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
    @brief Base utilities operating on meshes

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>

    This file contains a set of base utilities used to test mesh entities or
    operate on them
 */
#ifndef __MESH_UTILITIES__
#define __MESH_UTILITIES__

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/CurrentBoundaryFE.hpp>

//! \file mesh_util.h
//! \file mesh_util.h
/*!
  This file contains a set of functions to be used to test a 3D mesh and fix
  some possible problems. Some methods are general (2D and 3D meshes), some are
  specific to 3D meshes.

  They sometimes require in input a Switch paramenter sw and ostreams references.
  The switch values are
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
template <typename RegionMesh>
Real checkVolumes ( RegionMesh const& mesh,
                    std::vector<bool>& elSign,
                    Switch& sw )
{
    Real meas = 0.0;
    Real lmeas = 0.0;
    elSign.clear();
    elSign.reserve ( mesh.numVolumes() );
    typedef typename RegionMesh::elementShape_Type GeoShape;

    switch ( GeoShape::S_shape )
    {
        case TETRA:
        {
            CurrentFE fe ( feTetraP1, geoLinearTetra, quadRuleTetra1pt );
            for ( ID i = 0; i < mesh.numVolumes(); i++ )
            {
                fe.updateJac ( mesh.volume ( i ) );
                lmeas = fe.measure();
                meas += lmeas;
                elSign.push_back ( lmeas > 0.0 );
            }
        }
        break;
        case HEXA:
        {
            CurrentFE fe ( feHexaQ1, geoBilinearHexa, quadRuleHexa1pt );
            for ( ID i = 0; i < mesh.numVolumes(); i++ )
            {
                fe.updateJac ( mesh.volume ( i ) );
                lmeas = fe.measure();
                meas += lmeas;
                elSign.push_back ( lmeas > 0.0 );
            }
        }
        break;
        default:
            sw.create ( "SKIP_ORIENTATION_TEST", true );

            return 0;
    }

    if ( std::find ( elSign.begin(), elSign.end(), false ) != elSign.end() )
    {
        sw.create ( "HAS_NEGATIVE_VOLUMES", true );
    }

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
template <typename RegionMesh>
void fixVolumes ( RegionMesh& mesh,
                  const std::vector<bool>& elSign,
                  Switch& sw )
{

    for ( ID i = 0; i < mesh.numVolumes(); i++ )
    {
        if ( ! elSign[ i ] )
        {
            mesh.volume (i).reversePoints();
        }
    }
    sw.create ("HAS_VOLUMES", true);
}
//!\brief Computes volume enclosed by boundary faces
/*!
  It computes, for $i=0,1,2$, the integral \f$\int_{\partial \Omega} x_i n_i
  d\gamma \f$, \f$n_i\f$ being the i-th component of the boundary normal. If
  the domain boundary is properly discretised they should all return (within
  discretisation and truncation errors) the quantity \f$\vert\Omega\vert\f$.

  \warning Not to be used for accurate computations (it always adopts
  linear or bilinear elements, with a simple integration rule)
  \param mesh A 3D mesh
  \param vols returns 3 Real corresponding to the 3 integrals
*/
template <typename RegionMesh>
void getVolumeFromFaces ( RegionMesh const& mesh,
                          Real vols[ 3 ],
                          std::ostream& err = std::cerr )
{
    MeshUtility::GetCoordComponent getx ( 0 );
    MeshUtility::GetCoordComponent gety ( 1 );
    MeshUtility::GetCoordComponent getz ( 2 );
    vols[ 0 ] = 0.0;
    vols[ 1 ] = 0.0;
    vols[ 2 ] = 0.0;
    typedef typename RegionMesh::facetShape_Type GeoBShape;
    typedef typename RegionMesh::facet_Type facet_Type;
    typedef boost::shared_ptr<CurrentBoundaryFE> current_fe_type;

    current_fe_type bdfe;

    switch ( GeoBShape::S_shape )
    {
        case TRIANGLE:
            bdfe = current_fe_type ( new CurrentBoundaryFE ( feTriaP1, geoLinearTria,
                                                             quadRuleTria1pt ) );
            for ( ID i = 0; i < mesh.numBFaces(); i++ )
            {
                bdfe->update ( mesh.face ( i ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );
                vols[ 0 ] += bdfe->normalIntegral ( getx );
                vols[ 1 ] += bdfe->normalIntegral ( gety );
                vols[ 2 ] += bdfe->normalIntegral ( getz );
            }
            break;
        case QUAD:
            bdfe = current_fe_type ( new CurrentBoundaryFE ( feQuadQ1, geoBilinearQuad,
                                                             quadRuleQuad1pt ) );
            for ( ID i = 0; i < mesh.numBFaces(); i++ )
            {
                bdfe->update ( mesh.face ( i ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );
                vols[ 0 ] += bdfe->normalIntegral ( getx );
                vols[ 1 ] += bdfe->normalIntegral ( gety );
                vols[ 2 ] += bdfe->normalIntegral ( getz );
            }
            break;
        default:
            err << "Only tria and quad surface elements  may be checked for volume orientation at the moment" << std::endl;
            ASSERT0 ( false, "ABORT CONDITION OCCURRED" );
    }
}

//! Tests if the surface of the mesh is closed by computing surface integrals.
/*! It computes \f$\sum_{i=0}^2\int_{\partial \Omega} n_i d\gamma\f$.
  The value returned  should be very proximal to zero
 */
template <typename RegionMesh>
Real testClosedDomain ( RegionMesh const& mesh,
                        std::ostream& err = std::cerr )
{
    typedef typename RegionMesh::facet_Type facet_Type;

    typedef boost::shared_ptr<CurrentBoundaryFE> current_fe_type;
    current_fe_type bdfe;

    MeshUtility::GetOnes ones;
    Real test ( 0.0 );

    switch ( RegionMesh::facetShape_Type::S_shape )
    {
        case TRIANGLE:
            bdfe = current_fe_type ( new CurrentBoundaryFE ( feTriaP1, geoLinearTria,
                                                             quadRuleTria1pt ) );
            for ( ID i = 0; i < mesh.numBFaces(); i++ )
            {
                bdfe->update ( mesh.face ( i ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );
                test += bdfe->normalIntegral ( ones );
            }
            break;
        case QUAD:
            bdfe = current_fe_type ( new CurrentBoundaryFE ( feQuadQ1, geoBilinearQuad,
                                                             quadRuleQuad1pt ) );
            for ( ID i = 0; i < mesh.numBFaces(); i++ )
            {
                bdfe->update ( mesh.face ( i ), UPDATE_W_ROOT_DET_METRIC | UPDATE_NORMALS | UPDATE_QUAD_NODES );
                test += bdfe->normalIntegral ( ones );
            }

            break;

        default:
            err << "Only tria and quad surface elements  may be checked for volume orientation at the moment" << std::endl;
            ASSERT0 ( false, "ABORT CONDITION OCCURRED" );
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

template <typename RegionMesh>
bool checkMesh3D ( RegionMesh& mesh,
                   Switch& sw,
                   bool fix = true,
                   bool verbose = false,
                   std::ostream& out = std::cerr,
                   std::ostream& err = std::cerr,
                   std::ostream& clog = std::clog )
{
    verbose = verbose && ( mesh.comm()->MyPID() == 0 );
    typedef typename RegionMesh::point_Type point_Type;

    if ( mesh.storedPoints() == 0 )
    {
        err << "FATAL: mesh does not store points: I cannot do anything"
            << std::endl;
        sw.create ( "ABORT_CONDITION", true );
        sw.create ( "NOT_HAS_POINTS", true );
        return false;
    }

    if (verbose)
    {
        out << " Check point marker ids" << std::endl;
    }
    if ( !MeshUtility::checkIsMarkerSet ( mesh.pointList ) )
    {
        if (verbose)
        {
            err << "WARNING: Not all points have marker id set" << std::endl;
        }

        sw.create ( "POINTS_MARKER_UNSET", true );
    }


    //-------------------------------------------------------------------------
    //                                    VOLUMES
    //-------------------------------------------------------------------------


    if ( mesh.storedVolumes() == 0 )
    {
        if (verbose) err << "FATAL: mesh does not store volumes: I cannot do anything"
                             << std::endl;
        sw.create ( "ABORT_CONDITION", true );
        sw.create ( "NOT_HAS_VOLUMES", true );
        return false;
    }


    if ( !MeshUtility::checkId ( mesh.volumeList ) )
    {
        if (verbose)
        {
            err << "ERROR: volume ids were wrongly set" << std::endl;
            err << "FIXED" << std::endl;
        }
        if ( fix )
        {
            sw.create ( "FIXED_VOLUMES_ID", true );
        }
        if ( fix )
        {
            MeshUtility::fixId ( mesh.volumeList );
        }
    }

    if (verbose)
    {
        out << " Check volum marker ids" << std::endl;
    }
    if ( !MeshUtility::checkIsMarkerSet ( mesh.volumeList ) )
    {
        if (verbose)
        {
            err << "WARNING: Not all volumes have marker flag set" << std::endl;
        }

        sw.create ( "VOLUMES_MARKER_UNSET", true );

        if ( fix )
        {
            if (verbose)
            {
                out << "Fixing volume marker ids" << std::endl;
            }
            for ( typename RegionMesh::volumes_Type::iterator iv = mesh.volumeList.begin();
                    iv != mesh.volumeList.end(); ++iv )
            {
                if ( iv->isMarkerUnset() )
                {
                    iv->setMarkerID ( mesh.markerID() );
                }
            }
        }
    }

    if ( mesh.numElements() < mesh.storedVolumes() )
    {
        if (verbose)
            err << "WARNING: Mesh Volumes must be at least "
                << mesh.storedVolumes() << std::endl;
        if ( fix )
        {
            mesh.setNumVolumes ( mesh.storedVolumes() );
        }
        if ( fix )
        {
            sw.create ( "FIXED_VOLUME_COUNTER", true );
        }
    }

    // test now orientation

    boost::shared_ptr<std::vector<bool> > elSign ( new std::vector<bool> );

    Real meshMeasure = checkVolumes ( mesh, *elSign, sw );
    UInt positive;

    if ( sw.test ( "SKIP_ORIENTATION_TEST" ) )
    {
        if (verbose)
        {
            clog << "W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" << std::endl;
            err << "W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" << std::endl;
        }
    }
    else if ( sw.test ( "HAS_NEGATIVE_VOLUMES" ) )
    {
        if (verbose)
        {
            out << "Checking volume orientation" << std::endl;
        }
        positive = count ( elSign->begin(), elSign->end(), true );
        if ( verbose ) clog << positive << "W: positive elements out of"
                                << mesh.storedVolumes() << std::endl;
        if ( fix )
        {
            if ( verbose )
            {
                clog << "Fixing negative elements" << std::endl;
            }
            fixVolumes ( mesh, *elSign, sw );
        }

        if ( sw.test ( "ABORT_CONDITION" ) )
        {
            if (verbose) err << "ABORT: Cannot fix volumes, this element is not supported"
                                 << std::endl;
            return false;
        }
        else
        {
            sw.unset ( "HAS_NEGATIVE_VOLUMES" );
            meshMeasure = checkVolumes ( mesh, *elSign, sw );
            if ( sw.test ( "HAS_NEGATIVE_VOLUMES" ) )
            {
                if ( fix )
                {
                    if ( verbose )
                    {
                        err << "ABORT: Cannot fix volumes: something wrong with this mesh" << std::endl;
                    }
                    sw.create ( "ABORT_CONDITION", true );
                }
                return false;
            }
        }
    }

    if ( verbose ) clog << "Volume enclosed by the mesh= " << meshMeasure << std::endl
                            << "(Computed by integrating mesh elements measures)" << std::endl
                            << "(Using 1 point Quadrature rule)" << std::endl;

    //-----------------------------------------------------
    //                                    BOUNDARY FACES
    //-----------------------------------------------------

    boost::shared_ptr<MeshUtility::temporaryFaceContainer_Type> bfaces (
        new MeshUtility::temporaryFaceContainer_Type );
    UInt numInternalFaces, numFaces;

    if (verbose)
    {
        out << "Finding boundary faces from mesh topology" << std::endl;
    }

    UInt bFacesFound = MeshUtility::findBoundaryFaces ( mesh, *bfaces, numInternalFaces );

    numFaces = bFacesFound + numInternalFaces;

    MeshUtility::EnquireBFace<RegionMesh> enquireBFace (*bfaces );

    UInt meshNumBoundaryFaces ( mesh.numBFaces() + mesh.faceList.countElementsWithFlag ( EntityFlags::SUBDOMAIN_INTERFACE, &Flag::testOneSet ) );

    if ( mesh.storedFaces() == 0 ||
            meshNumBoundaryFaces > mesh.storedFaces() ||
            bFacesFound > mesh.storedFaces() || bFacesFound > meshNumBoundaryFaces )
    {
        // Something strange with boundary faces
        if (verbose)
        {
            err << "ERROR: Not all boundary faces stored" << std::endl;
            err << "Found " << bFacesFound << " stored " << mesh.storedFaces() <<
                "B faces declared in mesh " << meshNumBoundaryFaces << std::endl;
        }
        if ( fix )
        {
            sw.create ( "BUILD_BFACES", true );
            if (verbose)
            {
                out << "Building boundary faces from topology data" << std::endl;
            }
            MeshUtility::buildFaces ( mesh, clog, err, bFacesFound, numInternalFaces,
                                      true, false, false, bfaces.get() );
        }
        if (verbose)
        {
            err << "After buildFaces" << std::endl;
            err << "Found " << bFacesFound << " stored " << mesh.storedFaces()
                << "B faces declared in mesh " << meshNumBoundaryFaces << std::endl;
        }
    }
    else
    {
        if ( mesh.storedFaces() == 0 )
        {
            if ( verbose )
            {
                err << "ABORT CONDITION: cannot find boundary faces" << std::endl;
            }
            sw.create ( "NOT_HAS_FACES", true );
            sw.create ( "ABORT_CONDITION", true );
        }
        // The mesh declares to have the correct boundary faces. Yet, are we sure that the flag has been properly set
        // and that they are stored first? If fix is set we don't trust anybody!
        //
        // Make sure that flags are set using the info in bfaces, which depends ONLY on the mesh topology. No messing bout with markers
        // Make sure BFaces are stored first

        if (fix)
        {
            if (verbose)
            {
                out << "Rearranging faces so that boundary faces are first" << std::endl;
            }
            MeshUtility::rearrangeFaces ( mesh, clog, err, sw, numFaces, bFacesFound,
                                          verbose, bfaces.get() );
        }
        if ( meshNumBoundaryFaces !=  bFacesFound)
        {
            if ( verbose )
            {
                err << " ERROR: Number of B faces does not correspond to real one" << std::endl;
            }
            if (fix)
            {
                if ( verbose )
                {
                    err << "FIXED Number of B faces has been fixed to:" << bFacesFound << std::endl;
                }
                mesh.setNumBFaces ( bFacesFound);
            }
        }

        if ( !MeshUtility::checkId ( mesh.faceList ) )
        {
            if ( verbose )
            {
                err << "ERROR: face ids were wrongly set" << std::endl;
                err << "FIXED" << std::endl;
            }
            if ( fix )
            {
                sw.create ( "FIXED_FACES_ID", true );
                MeshUtility::fixId ( mesh.faceList );
            }
        }

        // Check Consistency with the mesh.

        if ( fix )
        {
            if (verbose)
            {
                out << "Make sure that adjacent elements of boundary faces are correct" << std::endl;
            }
            MeshUtility::fixBoundaryFaces ( mesh, clog, err, sw, numFaces, bFacesFound,
                                            false, verbose, bfaces.get() );
        }



        if ( mesh.numBFaces() == 0 )
        {
            if ( verbose )
            {
                err << " MeshBFaces counter is unset" << std::endl;
            }
            if ( fix )
            {
                mesh.setNumBFaces ( mesh.storedFaces() );
                sw.create ( "FIXED_BFACE_COUNTER", true );
                mesh.setLinkSwitch ( "HAS_BOUNDARY_FACETS" );
            }
        }


        if ( !MeshUtility::checkIsMarkerSet ( mesh.faceList ) )
        {
            if (verbose)
            {
                out << "Fix markers id for faces" << std::endl;
            }
            if ( verbose )
            {
                err << "WARNING: Not all faces have marker flag set" << std::endl;
            }
            sw.create ( "FACE_MARKER_UNSET", true );
            if ( fix )
            {
                MeshUtility::setBoundaryFacesMarker ( mesh, clog, err, verbose );
            }
            if ( fix && MeshUtility::checkIsMarkerSet ( mesh.faceList ) )
            {
                sw.create ( "FACE_MARKER_UNSET", false );
                sw.create ( "FACE_MARKER_FIXED", true );
            }
        }
    }


    if ( mesh.numFaces() != bFacesFound + numInternalFaces )
    {
        if ( verbose )
        {
            err << "WARNING Number of faces incorrectly set" << std::endl;
            err << "        It was       " << mesh.numFaces() << std::endl;
            err << "        It should be " << bFacesFound + numInternalFaces
                << std::endl;
        }
        if ( fix )
        {
            if ( verbose )
            {
                err << "        Fixing" << std::endl;
            }
            mesh.setNumFaces ( bFacesFound + numInternalFaces );
            sw.create ( "FIXED_FACE_COUNTER", true );
        }
    }

    if ( fix && mesh.storedFaces() == bFacesFound + numInternalFaces)
    {
        mesh.setLinkSwitch ( "HAS_ALL_FACETS" );
    }

    if (verbose)
    {
        out << std::endl;
        out << " Boundary faces found        :" << bFacesFound << std::endl;
        out << " Num Faces Stored stored     :" << mesh.storedFaces() << std::endl;
        out << " Num faces in mesh           :" << mesh.numFaces() << std::endl;
        out << " Boundary faces counter gives:" << mesh.numBFaces() << std::endl;
        out << std::endl;
    }
    //-----------------------------------------------------
    //                                    BOUNDARY EDGES
    //-----------------------------------------------------

    boost::shared_ptr<MeshUtility::temporaryEdgeContainer_Type> bedges (
        new MeshUtility::temporaryEdgeContainer_Type );

    UInt bEdgesFound = MeshUtility::findBoundaryEdges ( mesh, *bedges );
    MeshUtility::EnquireBEdge<RegionMesh> enquireBEdge (*bedges );

    UInt intedge (0);
    UInt Ned (0);
    MeshUtility::temporaryEdgeContainer_Type iedges;

    if ( mesh.storedEdges() == 0 ||
            mesh.numBEdges() > mesh.storedEdges() ||
            bEdgesFound > mesh.storedEdges() )
    {
        if (verbose)
        {
            err << "WARNING: mesh does not store (all) boundary edges" << std::endl;
        }
        sw.create ( "NOT_HAS_EDGES", true );
        if ( fix )
        {
            if (verbose)
            {
                out << "Build boundary edges" << std::endl;
            }
            MeshUtility::buildEdges ( mesh, clog, err, bEdgesFound, intedge, true, false,
                                      false, bedges.get() );
        }
        Ned = bEdgesFound + intedge;
        if ( fix )
        {
            sw.create ( "BUILD_BEDGES", true );
        }
    }
    else
    {

        // Make sure BEdges are first
        // Here I need to use a method that does not require the proper
        // setting of boundary Points!
        // With the edges I am being a bit sloppy. I am trusting the given mesh
        // todo do the same it was done for faces!
        if ( mesh.numBEdges() !=  bEdgesFound)
        {
            if ( verbose )
            {
                err << " ERROR: Number of BEdges does not correspond to real one" << std::endl;
            }
            if (fix)
            {
                if ( verbose )
                {
                    err << "FIXED Number of BEdges has been fixed to:" << bEdgesFound << std::endl;
                }
                mesh.setNumBEdges ( bEdgesFound);
            }
        }

        if ( fix )
        {
            if (verbose)
            {
                out << "Reorder edges so that boundary are first" << std::endl;
            }
            mesh.edgeList.reorderAccordingToFlag (EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);

            sw.create ( "FIXED_BEDGES_FIRST" );
        }
        if ( !MeshUtility::checkId ( mesh.edgeList ) )
        {
            if ( verbose )
            {
                err << "ERROR: edge ids were wrongly set" << std::endl;
            }
            if (fix)
            {
                if ( verbose )
                {
                    err << "FIXED" << std::endl;
                }
                sw.create ( "FIXED_EDGES_ID", true );
                MeshUtility::fixId ( mesh.edgeList );
            }
        }


        if ( !MeshUtility::checkIsMarkerSet ( mesh.edgeList ) )
        {
            if ( verbose )
            {
                err << "WARNING: Not all edges have marker flag set" << std::endl;
            }
            sw.create ( "EDGE_MARKER_UNSET", true );
            if ( fix )
            {
                if (verbose)
                {
                    out << "Fix boundary edges marker" << std::endl;
                }
                MeshUtility::setBoundaryEdgesMarker ( mesh, clog, err, verbose );
            }
            if ( fix && MeshUtility::checkIsMarkerSet ( mesh.edgeList ) )
            {
                sw.unset ( "EDGE_MARKER_UNSET" );
                sw.create ( "EDGE_MARKER_FIXED", true );
            }
        }
        if (verbose)
        {
            out << "Computing number of internal edges";
        }
        if ( fix )
        {
            Ned = bEdgesFound + MeshUtility::findInternalEdges ( mesh, *bedges, iedges );
        }
    }
    iedges.clear();
    MeshUtility::temporaryEdgeContainer_Type tmp;
    iedges.swap (tmp);

    if ( mesh.numBEdges() != bEdgesFound )
    {
        if ( verbose ) err << "WARNING: number of found boundary edges:" << bEdgesFound
                               << std::endl
                               << " does not match that declared in mesh, i.e. "
                               << mesh.numBEdges() << std::endl;
        if ( mesh.numBEdges() == 0 )
        {
            if ( fix )
            {
                if ( verbose )
                {
                    err << "FIXING" << std::endl;
                }
                sw.create ( "FIXED_BEDGES_COUNTER", true );
                mesh.setNumBEdges ( bEdgesFound );
            }
        }
    }

    if ( Ned != mesh.numEdges() )
    {
        if ( fix )
        {
            if ( verbose ) err << "WARNING: Counter of number of edges badly set: Should be (actual number)" << Ned << std::endl
                                   << "It is instead equal to " << mesh.numEdges() << std::endl;
            err << " **FIXED" << std::endl;
            mesh.setNumEdges ( Ned );
        }
    }
    UInt nbed;
    UInt counte = MeshUtility::testDomainTopology ( mesh, nbed );
    if ( counte == 0 )
    {
        if ( verbose )
        {
            out << "**DOMAIN SURFACE IS (TOPOLOGICALLY) CLOSED" << std::endl;
        }
    }
    else
    {
        sw.create ( "DOMAIN_NOT_CLOSED", true );
        if ( verbose ) err << "WARNING: DOMAIN APPEARS TO HAVE AN OPEN BOUNDARY (TOPOLOGY CHECK)" << std::endl
                               << "Number of inconsistent edges:" << counte << std::endl;
    }



    //-----------------------------------------------------
    //                                    POINTS
    //-----------------------------------------------------

    if (verbose)
    {
        out << "Checking vertexes" << std::endl;
    }
    UInt numVerticesFound = mesh.pointList.countElementsWithFlag (EntityFlags::VERTEX, &Flag::testOneSet);
    if (numVerticesFound != mesh.numVertices() )
    {
        if ( verbose )
        {
            err << "warning: The number of Points with vertex flag on does not coincide with the declared one." << std::endl;
        }
        if (fix)
        {
            if ( verbose )
            {
                err << "It will be fixed now" << std::endl;
            }
            // unset the flag. It will be remade
            for (UInt i = 0; i < mesh.numPoints(); ++i)
            {
                mesh.point (i).unSetFlag (EntityFlags::VERTEX);
            }
            // Find the real vertices and set the flag
            for (UInt i = 0; i < mesh.numElements(); ++i)
                for (UInt j = 0; j < mesh.numLocalVertices(); ++j)
                {
                    ID k = mesh.element (i).point (j).localId();
                    mesh.pointList (k).setFlag (EntityFlags::VERTEX);
                }
            numVerticesFound = mesh.pointList.countElementsWithFlag (
                                   EntityFlags::VERTEX, &Flag::testOneSet);
            mesh.setNumVertices (numVerticesFound);
            UInt numBVerticesFound = mesh.pointList.countElementsWithFlag (
                                         EntityFlags::VERTEX | EntityFlags::PHYSICAL_BOUNDARY, &Flag::testAllSet);
            mesh.setNumBVertices (numBVerticesFound);
            sw.create ( "FIXED_VERTICES_COUNTER", true );
        }
    }

    // Now that boundary faces have been correctly set we may work out
    // boundaty points

    if (fix)
    {
        if (verbose)
        {
            out << "Fix boundary points using boundary faces info" << std::endl;
        }
        MeshUtility::fixBoundaryPoints (mesh, clog, err, verbose);
    }

    MeshUtility::EnquireBPoint<RegionMesh> enquirebpoint ( mesh );

    UInt foundBPoints = mesh.pointList.countElementsWithFlag (EntityFlags::PHYSICAL_BOUNDARY,
                                                              &Flag::testOneSet);
    if (verbose)
    {
        out << "B Points Found " << foundBPoints << std::endl;
    }
    if ( foundBPoints == 0 || foundBPoints < mesh.storedBPoints() )
    {
        if ( verbose )
        {
            err << "ERROR Bpoints indicator not correctly set" << std::endl;
        }
    }
    else
    {
        if ( fix )
        {
            sw.create ( "FIXED_BOUNDARY_POINTS", true );
        }
    }

    // Now that we are sure that (jus) boundary points are flagged as such we check if the marker id is set
    if (verbose)
    {
        out << "Chsck point marker Ids" << std::endl;
    }
    if ( ! MeshUtility::checkIsMarkerSet ( mesh.pointList ) )
    {
        if (verbose)
        {
            err << "Points MARKER id incorrectly set" << std::endl;
        }

        if ( fix )
        {
            if ( verbose )
            {
                err << " Fixing Points Marker ID" << std::endl << "If unset the boundary will inherit the strongest among faces" << std::endl;
            }
            if ( verbose )
            {
                err << " The internal will be set to the domain flag" << std::endl;
            }
            MeshUtility::setBoundaryPointsMarker ( mesh, clog, err, false );
            // fix marker at interior points. It takes
            if ( ! MeshUtility::checkIsMarkerSet ( mesh.pointList ) )
            {
                // Maybe boundary points marker is fine, this is enough
                bool boundaryIsOk (true);
                std::vector<point_Type const*>
                listOfPt = mesh.pointList.extractElementsWithFlag (
                               EntityFlags::PHYSICAL_BOUNDARY, &Flag::testOneSet);
                for (typename std::vector<point_Type const*>::const_iterator
                        it = listOfPt.begin();
                        it < listOfPt.end();
                        ++it)
                {
                    boundaryIsOk = boundaryIsOk | (*it)->isMarkerSet();
                }
                std::vector<point_Type const*>().swap (listOfPt); // save memory
                if ( verbose )
                {
                    clog << " Marker ID on boundary points is fine. Internal points may have marker unset" << std::endl;
                }
                if (verbose)
                {
                    err << "Cannot Fix Points MARKER" << std::endl;
                }
                if ( verbose && boundaryIsOk )
                {
                    err << "But boundary points are fine" << std::endl;
                }
                sw.create ( "POINT_MARKER_UNSET", true );
            }
            else
            {
                if (verbose)
                {
                    err << "FIXED" << std::endl;
                }
                sw.create ( "FIXED_POINT_MARKER", true );
            }
        }
    }


    if ( mesh.storedBPoints() == 0 )
    {
        if ( verbose )
        {
            err << "WARNING B. Points COUNTER incorrectly set" << std::endl;
        }
        if ( fix )
        {
            MeshUtility::setBoundaryPointsCounters ( mesh ) ;
            if ( verbose )
            {
                err << " FIXED" << std::endl;
            }
            sw.create ( "FIXED_BPOINTS_COUNTER", true );
        }
    }

    if ( mesh.numPoints() == 0 )
    {
        if ( verbose )
        {
            err << "WARNING Points Counter unset" << std::endl;
        }
        if ( fix )
        {
            mesh.numPoints() = mesh.storedPoints();
            sw.create ( "FIXED_POINTS_COUNTER", true );
        }
    }

    //-----------------------------------------------------
    //                                   FINAL CHECKS
    //-----------------------------------------------------
    if ( verbose ) out << " ********     COUNTERS CONTENT **********************************" << std::endl

                           << " Num Volumes    : " << mesh.numVolumes() << std::endl
                           << " Num Vertices   : " << mesh.numVertices() << std::endl
                           << " Num B. Vertices: " << mesh.numBVertices() << std::endl
                           << " Num Points     : " << mesh.numPoints() << std::endl
                           << " Num B. Points  : " << mesh.numBPoints() << std::endl
                           << " Num Edges      : " << mesh.numEdges() << std::endl
                           << " Num B. Edges   : " << mesh.numBEdges() << std::endl
                           << " Num Faces      : " << mesh.numFaces() << std::endl
                           << " Num B. Faces   : " << mesh.numBFaces() << std::endl
                           << " ********     END COUNTERS **********************************"
                           << std::endl;

    bool eulok1 = ( 2 * mesh.numFaces() -
                    mesh.numLocalFaces() * mesh.numVolumes() -
                    mesh.numBFaces() ) == 0;

    bool eulok2 ( true );

    if ( RegionMesh::elementShape_Type::S_shape == TETRA )
    {
        if ( verbose )
        {
            out << std::endl << "Checking Euler formulae: ";
        }
        eulok2 = ( mesh.numEdges() -
                   mesh.numVolumes() -
                   mesh.numVertices() -
                   ( 3 * mesh.numBFaces() -
                     2 * mesh.numBVertices() ) / 4 ) == 0;
    }

    if ( ! ( eulok1 && eulok2 ) )
    {
        if ( verbose ) err << "WARNING: The following Euler formula(s) are not satisfied"
                               << std::endl;
        sw.create ( "NOT_EULER_OK" );
    }
    else
    {
        if ( verbose )
        {
            out << std::endl << " ok." << std::endl;
        }
    }

    if ( !eulok1 )
    {
        if ( verbose ) err << "  2*nFaces = nFacesPerVolume*nVolumes + nBoundaryFaces"
                               << std::endl
                               << "  2*" << mesh.numFaces() << " != " << mesh.numLocalFaces()
                               << " * " << mesh.numVolumes() << " + " << mesh.numBFaces()
                               << std::endl;
    }

    if ( !eulok2 )
    {
        if ( verbose ) err << "  nEdges = nVolumes + nVertices + (3*nBoundaryFaces - 2*nBoundaryVertices)/4" << std::endl
                               << "  " << mesh.numEdges() << " != " << mesh.numVolumes() << " + "
                               << mesh.numVertices() << " + (3*" << mesh.numBFaces() << " - 2*"
                               << mesh.numBVertices() << ")/4" << std::endl;
    }

    mesh.setLinkSwitch ( "HAS_BEEN_CHECKED" );

    return true;
}
}
#endif
