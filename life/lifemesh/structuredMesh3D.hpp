//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2009-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Contains methods which generate structured meshes.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 16 April 2010

    Such methods will be usefull in order to test problems at different
    scales.
 */

#ifndef STRUCTUREDMESH3D_H
#define STRUCTUREDMESH3D_H 1

#include <life/lifecore/life.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <fstream>

namespace LifeV {

/*
namespace structuredMeshFlags {
// The boundary of the structured mesh
// is composed by 26 labels.

// Walls
const int LEFTWALL   = 4;
const int RIGHTWALL  = 2;
const int FRONTWALL  = 1;
const int BACKWALL   = 3;
const int TOPWALL    = 6;
const int BOTTOMWALL = 5;
// Edges
const int BOTTOMEDGE1 =  7;
const int BOTTOMEDGE2 =  8;
const int BOTTOMEDGE3 =  9;
const int BOTTOMEDGE4 = 10;
const int SIDEEDGE1   = 11;
const int SIDEEDGE2   = 12;
const int SIDEEDGE3   = 13;
const int SIDEEDGE4   = 14;
const int TOPEDGE1    = 15;
const int TOPEDGE2    = 16;
const int TOPEDGE3    = 17;
const int TOPEDGE4    = 18;
// Corners
const int BOTTOMCORNER1 = 19;
const int BOTTOMCORNER2 = 20;
const int BOTTOMCORNER3 = 21;
const int BOTTOMCORNER4 = 22;
const int TOPCORNER1    = 23;
const int TOPCORNER2    = 24;
const int TOPCORNER3    = 25;
const int TOPCORNER4    = 26;
}
*/

//! @name Methods
//@{

//! This method gives the flag for a cube.
/*!
  @param i_x Number of elements along the length
  @param i_y Number of elements along the width
  @param i_z Number of elements along the height
  @param l_x length of the mesh
  @param l_y width of the mesh
  @param l_z height of the mesh

  The internal points are labelled with 0.
  The labels 1-6 are reserved for the 6 faces.
  The labels 7-18 are reserved for the 12 edges.
  The labels 19-26 are reserved for the 8 corners.
*/
EntityFlag regularMeshPointPosition(const UInt& i_x,
                                    const UInt& i_y,
                                    const UInt& i_z,
                                    const UInt& n_x,
                                    const UInt& n_y,
                                    const UInt& n_z);

//! This method generate a structured mesh
/*!
  @param mesh The mesh that we want to generate
  @param regionFlag Flag of the region
  @param m_x Number of elements along the length
  @param m_y Number of elements along the width
  @param m_z Number of elements along the height
  @param l_x length of the mesh
  @param l_y width of the mesh
  @param l_z height of the mesh
  @param verbose Verbose mode enabled/disabled
*/
template <typename GeoShape, typename MC>
void regularMesh3D( RegionMesh3D<GeoShape,MC>& mesh,
                    EntityFlag regionFlag,
                    const UInt& m_x,
                    const UInt& m_y,
                    const UInt& m_z,
                    bool verbose=false,
                    const Real& l_x=1.0,
                    const Real& l_y=1.0,
                    const Real& l_z=1.0,
                    const Real& t_x=0.0,
                    const Real& t_y=0.0,
                    const Real& t_z=0.0
                    ){
    // output stream
    std::stringstream discardedLog;
    std::ostream& oStr = verbose ? std::cout : discardedLog;

    // discretization
    Real dx(l_x/m_x);
    Real dy(l_y/m_y);
    Real dz(l_z/m_z);

    // Number of nodes along the side of the unit cube
    UInt n_x(m_x+1);
    UInt n_y(m_y+1);
    UInt n_z(m_z+1);

    // Incremental values in order to get the indices of the nodes
    // Due to the structure, if we add N_i to the number of the
    // current node, we end up with the number of the next point
    // in the i axis.
    UInt N_x(1);
    UInt N_y(n_x);
    UInt N_z(n_x*n_y);

    // We calculate in advance the bound which
    // will determine the flip of the reference
    // cube.
    UInt mid_x(m_x/2);
    UInt mid_y(m_y/2);
    UInt mid_z(m_z/2);

    // Data about the mesh
    UInt verticesNumber(n_x*n_y*n_z);
    UInt boundaryVerticesNumber(verticesNumber-(n_x-2)*(n_y-2)*(n_z-2)); //Total-inside points
    UInt boundaryEdgesNumber(
                             2*(m_x*(n_y-2)+m_y*(n_x-2) + m_x*m_y) +
                             2*(m_x*(n_z-2)+m_z*(n_x-2) + m_x*m_z) +
                             2*(m_y*(n_z-2)+m_z*(n_y-2) + m_y*m_z) +
                             4*(m_x+m_y+m_z)
                             );
    UInt edgesNumber(
                     // Edges that draws the cuboids
                     m_x*n_y*n_z +
                     n_x*m_y*n_z +
                     n_x*n_y*m_z +
                     // Edges that divide into two faces the face of the cuboids
                     m_x*m_y*n_z +
                     n_x*m_y*m_z +
                     m_x*n_y*m_z +
                     // Edges that go accross the cuboids
                     m_x*m_y*m_z
                     );
    UInt boundaryFacesNumber(2*2*(m_x*m_y+m_y*m_z+m_x*m_z));
    // Faces per cuboids = 12 + 6 internal
    UInt facesNumber( (m_x*m_y*m_z*12+boundaryFacesNumber)/2 + 6*m_x*m_y*m_z);
    UInt volumesNumber(m_x*m_y*m_z*6);

    UInt pointsNumber(0);
    UInt boundaryPointsNumber(0);
    if(GeoShape::numPoints>4){
        std::cout<< "Quadratic Tetra Mesh (from Linear geometry)" << std::endl;
        // In this case there is one extra points on each edge
        pointsNumber = verticesNumber + edgesNumber;
        boundaryPointsNumber = boundaryVerticesNumber + boundaryEdgesNumber;
    }else{
        std::cout<< "Linear Tetra Mesh" << std::endl;
        pointsNumber = verticesNumber;
        boundaryPointsNumber = boundaryVerticesNumber;
    }

    // Set the data
    std::cout<<"initialization of mesh...";

    // Note: The vertices are the nodes of the mesh while the points
    //       are the nodes and some new points added for example for
    //       the quadratic tetra

    // About points:
    mesh.setNumBPoints(boundaryPointsNumber);
    mesh.setMaxNumPoints(pointsNumber,true);
    mesh.setMaxNumGlobalPoints(pointsNumber);

    // About vertices:
    mesh.setNumVertices(verticesNumber);
    mesh.setNumGlobalVertices(verticesNumber);
    mesh.setNumBVertices(boundaryVerticesNumber);


    // About edges:
    mesh.setNumEdges(edgesNumber);
    mesh.setNumBEdges(boundaryEdgesNumber);
    //mesh.setMaxNumEdges(edgesNumber); //mpp++
    mesh.setMaxNumEdges(boundaryEdgesNumber);
    mesh.setMaxNumGlobalEdges(edgesNumber);

    // About faces:
    mesh.setNumFaces(facesNumber);
    mesh.setNumBFaces(boundaryFacesNumber);
    //mesh.setMaxNumFaces(facesNumber); // i.g. the # of stored faces
    mesh.setMaxNumFaces(boundaryFacesNumber); //try
    mesh.setMaxNumGlobalFaces(boundaryFacesNumber);

    // About volumes
    mesh.setMaxNumVolumes(volumesNumber,true);
    mesh.setMaxNumGlobalVolumes(volumesNumber);

    mesh.setMarker(regionFlag);
    std::cout<<"done"<<std::endl;

    // Declaration of pointers on the different mesh entities
    typename RegionMesh3D<GeoShape,MC>::PointType*  pointPtr  = 0;
    //typename RegionMesh3D<GeoShape,MC>::EdgeType*   edgePtr   = 0;
    //typename RegionMesh3D<GeoShape,MC>::FaceType*   facePtr   = 0;
    typename RegionMesh3D<GeoShape,MC>::VolumeType* volumePtr = 0;

    // Build the points of the mesh
    std::cout<<"building the points of the mesh...";
    Real xPosition(0.0),yPosition(0.0),zPosition(0.0);
    EntityFlag nodeFlag(0);
    UInt nodeID(0);
    UInt P0(0),P1(0),P2(0),P3(0),P4(0),P5(0),P6(0),P7(0);

    for(UInt k(0);k<n_z;++k)
    {
        zPosition = dz*k;

        for(UInt j(0);j<n_y;++j)
        {
            yPosition = dy*j;

            for(UInt i(0);i<n_x;++i)
            {
                xPosition = dx*i;
                nodeFlag = regularMeshPointPosition(i,j,k,n_x,n_y,n_z);

                // We create the point
                if(nodeFlag>0){
                    pointPtr = &mesh.addPoint(true); //it is a boundary point
                }else{
                    pointPtr = &mesh.addPoint(false);
                }

                // We set the point properties
                nodeID = k*N_z+j*N_y+i+1;
                pointPtr->setId(nodeID);
                pointPtr->setLocalId(nodeID);

                mesh.localToGlobalNode().insert(std::make_pair(nodeID,nodeID));
                mesh.globalToLocalNode().insert(std::make_pair(nodeID,nodeID));

                pointPtr->setMarker(nodeFlag);
                pointPtr->x() = xPosition+t_x;
                pointPtr->y() = yPosition+t_y;
                pointPtr->z() = zPosition+t_z;
            }
        }
    }
    std::cout<<"done"<<std::endl;

    /*
    // Build the boundary edges (ONLY!)
    std::cout<<"building the boundary edges...";

     //edges l_x x l_y
    for(UInt j(0);j<m_y;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        y
            // P2--P3 ^
            // |   |  |
            // P0--P1 o-->x
            //       /
            //      z
            P0 = j*N_y+i+1;
            P1 = P0+N_x;
            P2 = P0+N_y;
            P3 = P0+N_x+N_y;

            // Diagonal bottom - marker = 5
            if((i+1<=mid_x && j+1<=mid_y)||(i+1>mid_x && j+1>mid_y)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 bottom
            if(i>0){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 2 bottom
            if(j>0){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P1));
            }

            // The ex,ey,ez unit vectors are oriented
            //        y
            // P2--P3 ^
            // |   |  |
            // P0--P1 o-->x
            //       /
            //      z
            P0 += N_z*(n_z-1);
            P1 = P0+N_x;
            P2 = P0+N_y;
            P3 = P0+N_x+N_y;

            // Diagonal top - marker = 6
            if((i+1<=mid_x && j+1<=mid_y)||(i+1>mid_x && j+1>mid_y)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 top
            if(i==0){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 2 top
            if(j==0){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P1));
            }

            // Edge 3 top
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P3).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P3));

            // Edge 4 top
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P2).marker(),mesh.point(P3).marker());
            edgePtr->setPoint(1,mesh.point(P2));
            edgePtr->setPoint(2,mesh.point(P3));
        }
    }
    //edges l_x x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^ y
            // |   |  |/
            // P0--P1 0-->x
            P0 = j*N_z+i+1;
            P1 = P0+N_x;
            P2 = P0+N_z;
            P3 = P0+N_x+N_z;

            // Diagonal front - marker = 1
            if((i+1<=mid_x && j+1<=mid_z)||(i+1>mid_x && j+1>mid_z)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 front
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 2 front
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P1));

            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^ y
            // |   |  |/
            // P0--P1 0-->x
            P0 += N_y*(n_y-1);
            P1 = P0+N_x;
            P2 = P0+N_z;
            P3 = P0+N_x+N_z;

            // Diagonal back - marker = 3
            if((i+1<=mid_x && j+1<=mid_z)||(i+1>mid_x && j+1>mid_z)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 back
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P3).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P3));

            // Edge 2 back
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P1));
        }
    }
    //edges l_y x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_y;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^
            // |   |  |
            // P0--P1 0-->y
            //       /
            //      x
            P0 = j*N_z+i*N_y+1;
            P1 = P0+N_y;
            P2 = P0+N_z;
            P3 = P0+N_y+N_z;

            // Diagonal left - marker = 4
            if((i+1<=mid_y && j+1<=mid_z)||(i+1>mid_y && j+1>mid_z)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 left
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P1));

            // Edge 2 left
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P3).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P3));

            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^
            // |   |  |
            // P0--P1 0-->y
            //       /
            //      x
            P0 += N_x*(n_x-1);
            P1 = P0+N_y;
            P2 = P0+N_z;
            P3 = P0+N_y+N_z;

            // Diagonal right - marker = 2
            if((i+1<=mid_y && j+1<=mid_z)||(i+1>mid_y && j+1>mid_z)){
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P3).marker());
                edgePtr->setPoint(1,mesh.point(P0));
                edgePtr->setPoint(2,mesh.point(P3));
            }else{
                edgePtr = &mesh.addEdge(true);
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 1 right
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P1).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P1));

            // Edge 2 right
            edgePtr = &mesh.addEdge(true);
            edgePtr->setWeakerMarker(mesh.point(P0).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P0));
            edgePtr->setPoint(2,mesh.point(P2));
        }
    }
    std::cout<<"done"<<std::endl;

    // Build the boundary faces (ONLY!)
    std::cout<<"building the boundary faces...";
    // Faces l_x x l_y
    for(UInt j(0);j<m_y;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        y
            // P2--P3 ^
            // |   |  |
            // P0--P1 o-->x
            //       /
            //      z
            P0 = j*N_y+i+1;
            P1 = P0+N_x;
            P2 = P0+N_y;
            P3 = P0+N_x+N_y;

            if((i+1<=mid_x && j+1<=mid_y)||(i+1>mid_x && j+1>mid_y)){
                // Face 1 bottom - marker = 5
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(5);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P1));

                // Face 2 bottom
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(5);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));
            }else{
                // Face 1 bottom - marker = 5
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(5);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P1));

                // Face 2 bottom
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(5);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));
            }

            // The ex,ey,ez unit vectors are oriented
            //        y
            // P2--P3 ^
            // |   |  |
            // P0--P1 o-->x
            //       /
            //      z
            P0 += N_z*(n_z-1);
            P1 = P0+N_x;
            P2 = P0+N_y;
            P3 = P0+N_x+N_y;

            if((i+1<=mid_x && j+1<=mid_y)||(i+1>mid_x && j+1>mid_y)){
                // Face 1 top - marker = 6
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(6);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P3));

                // Face 2 top
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(6);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));
            }else{
                // Face 1 top - marker = 6
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(6);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P2));

                // Face 2 top
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(6);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));
            }
        }
    }
    // Faces l_x x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^ y
            // |   |  |/
            // P0--P1 0-->x
            P0 = j*N_z+i+1;
            P1 = P0+N_x;
            P2 = P0+N_z;
            P3 = P0+N_x+N_z;

            if((i+1<=mid_x && j+1<=mid_z)||(i+1>mid_x && j+1>mid_z)){
                // Face 1 front - marker = 1
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(1);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P3));

                // Face 2 front
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(1);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));
            }else{
                // Face 1 front - marker = 1
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(1);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P2));

                // Face 2 front
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(1);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));
            }

            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^ y
            // |   |  |/
            // P0--P1 0-->x
            P0 += N_y*(n_y-1);
            P1 = P0+N_x;
            P2 = P0+N_z;
            P3 = P0+N_x+N_z;

            if((i+1<=mid_x && j+1<=mid_z)||(i+1>mid_x && j+1>mid_z)){
                // Face 1 back - marker = 3
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(3);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));

                // Face 2 back
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(3);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P1));
            }else{
                // Face 1 back - marker = 3
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(3);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P1));

                // Face 2 back
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(3);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));
            }
        }
    }
    // Faces l_y x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_y;++i)
        {
            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^
            // |   |  |
            // P0--P1 0-->y
            //       /
            //      x
            P0 = j*N_z+i*N_y+1;
            P1 = P0+N_y;
            P2 = P0+N_z;
            P3 = P0+N_y+N_z;

            if((i+1<=mid_y && j+1<=mid_z)||(i+1>mid_y && j+1>mid_z)){
                // Face 1 left - marker = 4
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(4);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));

                // Face 2 left
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(4);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P1));
            }else{
                // Face 1 left - marker = 4
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(4);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P1));

                // Face 2 left
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(4);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P2));
                facePtr->setPoint(3,mesh.point(P3));
            }

            // The ex,ey,ez unit vectors are oriented
            //        z
            // P2--P3 ^
            // |   |  |
            // P0--P1 0-->y
            //       /
            //      x
            P0 += N_x*(n_x-1);
            P1 = P0+N_y;
            P2 = P0+N_z;
            P3 = P0+N_y+N_z;

            if((i+1<=mid_y && j+1<=mid_z)||(i+1>mid_y && j+1>mid_z)){
                // Face 1 right - marker = 2
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(2);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));

                // Face 2 right
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(2);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P3));
            }else{
                // Face 1 right - marker = 2
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(2);
                facePtr->setPoint(1,mesh.point(P0));
                facePtr->setPoint(2,mesh.point(P1));
                facePtr->setPoint(3,mesh.point(P2));

                // Face 2 right
                facePtr = &mesh.addFace(true);
                facePtr->setMarker(2);
                facePtr->setPoint(1,mesh.point(P1));
                facePtr->setPoint(2,mesh.point(P3));
                facePtr->setPoint(3,mesh.point(P2));
            }
        }
    }
    std::cout<<"done"<<std::endl;
*/

    // Build the volumes
    std::cout<<"building the volumes...";
    UInt volumeID(0);
    for(UInt k(0);k<m_z;++k)
    {
        for(UInt j(0);j<m_y;++j)
        {
            for(UInt i(0);i<m_x;++i)
            {
                volumeID=(k*m_y*m_x+j*m_y+i)*6 ;

                if((i+1<=mid_x && j+1<=mid_y && k+1<=mid_z)||(i+1>mid_x && j+1>mid_y && k+1>mid_z)){
                    // Zone 0,7
                    nodeID = k*N_z+j*N_y+i;
                    P4=nodeID+1;
                    P6=nodeID+N_x+1;
                    P5=nodeID+N_y+1;
                    P7=nodeID+N_x+N_y+1;
                    P0=nodeID+N_z+1;
                    P2=nodeID+N_x+N_z+1;
                    P1=nodeID+N_y+N_z+1;
                    P3=nodeID+N_x+N_y+N_z+1;
                }else if((i+1>mid_x && j+1<=mid_y && k+1<=mid_z)||(i+1<=mid_x && j+1>mid_y && k+1>mid_z)){
                    // Zone 1,6
                    nodeID = k*N_z+j*N_y+i;
                    P1=nodeID+1;
                    P3=nodeID+N_x+1;
                    P0=nodeID+N_y+1;
                    P2=nodeID+N_x+N_y+1;
                    P5=nodeID+N_z+1;
                    P7=nodeID+N_x+N_z+1;
                    P4=nodeID+N_y+N_z+1;
                    P6=nodeID+N_x+N_y+N_z+1;
                }else if((i+1<=mid_x && j+1>mid_y && k+1<=mid_z)||(i+1>mid_x && j+1<=mid_y && k+1>mid_z)){
                    // Zone 2,5
                    nodeID = k*N_z+j*N_y+i;
                    P2=nodeID+1;
                    P0=nodeID+N_x+1;
                    P3=nodeID+N_y+1;
                    P1=nodeID+N_x+N_y+1;
                    P6=nodeID+N_z+1;
                    P4=nodeID+N_x+N_z+1;
                    P7=nodeID+N_y+N_z+1;
                    P5=nodeID+N_x+N_y+N_z+1;
                }else if((i+1<=mid_x && j+1<=mid_y && k+1>mid_z)||(i+1>mid_x && j+1>mid_y && k+1<=mid_z)){
                    // Zone 3,4
                    nodeID = k*N_z+j*N_y+i;
                    P0=nodeID+1;
                    P1=nodeID+N_x+1;
                    P2=nodeID+N_y+1;
                    P3=nodeID+N_x+N_y+1;
                    P4=nodeID+N_z+1;
                    P5=nodeID+N_x+N_z+1;
                    P6=nodeID+N_y+N_z+1;
                    P7=nodeID+N_x+N_y+N_z+1;
                }

                // Tetra 1
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+1);
                volumePtr->setLocalId(volumeID+1);
                volumePtr->setPoint(1, mesh.point(P0));
                volumePtr->setPoint(2, mesh.point(P1));
                volumePtr->setPoint(3, mesh.point(P3));
                volumePtr->setPoint(4, mesh.point(P4));
                volumePtr->setMarker(regionFlag);

                // Tetra 2
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+2);
                volumePtr->setLocalId(volumeID+2);
                volumePtr->setPoint(1, mesh.point(P1));
                volumePtr->setPoint(2, mesh.point(P3));
                volumePtr->setPoint(3, mesh.point(P4));
                volumePtr->setPoint(4, mesh.point(P5));
                volumePtr->setMarker(regionFlag);

                // Tetra 3
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+3);
                volumePtr->setLocalId(volumeID+3);
                volumePtr->setPoint(1, mesh.point(P4));
                volumePtr->setPoint(2, mesh.point(P5));
                volumePtr->setPoint(3, mesh.point(P3));
                volumePtr->setPoint(4, mesh.point(P7));
                volumePtr->setMarker(regionFlag);

                // Tetra 4
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+4);
                volumePtr->setLocalId(volumeID+4);
                volumePtr->setPoint(1, mesh.point(P0));
                volumePtr->setPoint(2, mesh.point(P3));
                volumePtr->setPoint(3, mesh.point(P2));
                volumePtr->setPoint(4, mesh.point(P4));
                volumePtr->setMarker(regionFlag);

                // Tetra 5
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+5);
                volumePtr->setLocalId(volumeID+5);
                volumePtr->setPoint(1, mesh.point(P6));
                volumePtr->setPoint(2, mesh.point(P3));
                volumePtr->setPoint(3, mesh.point(P4));
                volumePtr->setPoint(4, mesh.point(P2));
                volumePtr->setMarker(regionFlag);

                // Tetra 6
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+6);
                volumePtr->setLocalId(volumeID+6);
                volumePtr->setPoint(1, mesh.point(P7));
                volumePtr->setPoint(2, mesh.point(P6));
                volumePtr->setPoint(3, mesh.point(P3));
                volumePtr->setPoint(4, mesh.point(P4));
                volumePtr->setMarker(regionFlag);
            }
        }
    }
    std::cout<<"done"<<std::endl;

    // Build a P2 mesh from a P1 geometry
    if(GeoShape::numPoints > 4) p1top2(mesh);

    // Test mesh
    Switch sw;

    // Correction JFG
    if(!checkMesh3D(mesh, sw, true, false, oStr, std::cerr, oStr)) abort();

    Real vols[3];
    getVolumeFromFaces(mesh, vols,oStr);
    oStr << "   VOLUME ENCLOSED BY THE MESH COMPUTED BY INTEGRATION ON"<<
        " BOUNDARY FACES"<<std::endl;
    oStr << "INT(X)     INT(Y)      INT(Z) <- they should be equal and equal to"
         << std::endl
         << "                                 the volume enclosed by the mesh "
         << std::endl;
    oStr << vols[0] << " " << vols[1] << " " << vols[2] << std::endl;

    oStr << "   BOUNDARY FACES ARE DEFINING A CLOSED SURFACE IF "
         << testClosedDomain(mesh,oStr) << std::endl
         << " IS (ALMOST) ZERO" << std::endl;

    // Updates the connectivity of the mesh
    mesh.updateElementEdges( true, verbose );
    mesh.updateElementFaces( true, verbose );
}

//@}



} // Namespace LifeV

#endif /* STRUCTUREDMESH3D_H */
