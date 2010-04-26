//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
    mesh.setMaxNumEdges(boundaryEdgesNumber);
    //mesh.setMaxNumGlobalEdges(boundaryEdgesNumber);
    mesh.setMaxNumGlobalEdges(edgesNumber);

    // About faces:
    mesh.setNumFaces(facesNumber);
    mesh.setNumBFaces(boundaryFacesNumber);
    //mesh.setMaxNumFaces(boundaryFacesNumber);
    mesh.setMaxNumFaces(facesNumber);
    mesh.setMaxNumGlobalFaces(boundaryFacesNumber);

    // About volumes
    mesh.setMaxNumVolumes(volumesNumber,true);
    mesh.setMaxNumGlobalVolumes(volumesNumber);

    mesh.setMarker(regionFlag);
    std::cout<<"done"<<std::endl;

    // Declaration of pointers on the different mesh entities
    typename RegionMesh3D<GeoShape,MC>::PointType*  pointPtr  = 0;
    typename RegionMesh3D<GeoShape,MC>::EdgeType*   edgePtr   = 0;
    typename RegionMesh3D<GeoShape,MC>::FaceType*   facePtr   = 0;
    typename RegionMesh3D<GeoShape,MC>::VolumeType* volumePtr = 0;

    // Build the points of the mesh
    std::cout<<"building the points of the mesh...";
    Real xPosition,yPosition,zPosition;
    UInt nodeID;
    EntityFlag nodeFlag;

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


    // Build the boundary edges (ONLY!)
    std::cout<<"building the boundary edges...";
    UInt P1,P2;

     //edges l_x x l_y
    for(UInt j(0);j<m_y;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            nodeID = j*N_y+i+1;

            // Diagonal bottom - marker = 5
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 bottom
            if(i>0){
                edgePtr = &mesh.addEdge(true);
                P1=nodeID;
                P2=nodeID+N_y;
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 2 bottom
            if(j>0){
                edgePtr = &mesh.addEdge(true);
                P1=nodeID;
                P2=nodeID+N_x;
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            nodeID += N_z*(n_z-1);

            // Diagonal top - marker = 6
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 top
            if(i==0){
                edgePtr = &mesh.addEdge(true);
                P1=nodeID;
                P2=nodeID+N_y;
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 2 top
            if(j==0){
                edgePtr = &mesh.addEdge(true);
                P1=nodeID;
                P2=nodeID+N_x;
                edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
                edgePtr->setPoint(1,mesh.point(P1));
                edgePtr->setPoint(2,mesh.point(P2));
            }

            // Edge 3 top
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_x+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 4 top
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_y;
            P2=nodeID+N_x+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));
        }
    }
    //edges l_x x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            nodeID = j*N_z+i+1;

            // Diagonal front - marker = 1
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 front
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 2 front
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_x;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            nodeID += N_y*(n_y-1);

            // Diagonal back - marker = 3
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 back
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_x;
            P2=nodeID+N_x+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 2 back
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_x;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));
        }
    }
    //edges l_y x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_y;++i)
        {
            nodeID = j*N_z+i*N_y+1;

            // Diagonal left - marker = 4
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_y;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 left
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 2 left
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_y;
            P2=nodeID+N_y+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            nodeID += N_x*(n_x-1);

            // Diagonal right - marker = 2
            edgePtr = &mesh.addEdge(true);
            P1=nodeID+N_y;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 1 right
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_y;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
            edgePtr->setPoint(2,mesh.point(P2));

            // Edge 2 right
            edgePtr = &mesh.addEdge(true);
            P1=nodeID;
            P2=nodeID+N_z;
            edgePtr->setWeakerMarker(mesh.point(P1).marker(),mesh.point(P2).marker());
            edgePtr->setPoint(1,mesh.point(P1));
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
            nodeID = j*N_y+i+1;

            // Face 1 bottom - marker = 5
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(5);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_y));
            facePtr->setPoint(3,mesh.point(nodeID+N_x));

            // Face 2 bottom
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(5);
            facePtr->setPoint(1,mesh.point(nodeID+N_x));
            facePtr->setPoint(2,mesh.point(nodeID+N_y));
            facePtr->setPoint(3,mesh.point(nodeID+N_x+N_y));

            nodeID += N_z*(n_z-1);

            // Face 1 top - marker = 6
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(6);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_x));
            facePtr->setPoint(3,mesh.point(nodeID+N_y));

            // Face 2 top
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(6);
            facePtr->setPoint(1,mesh.point(nodeID+N_x));
            facePtr->setPoint(2,mesh.point(nodeID+N_x+N_y));
            facePtr->setPoint(3,mesh.point(nodeID+N_y));
        }
    }
    // Faces l_x x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_x;++i)
        {
            nodeID = j*N_z+i+1;

            // Face 1 front - marker = 1
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(1);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_x));
            facePtr->setPoint(3,mesh.point(nodeID+N_z));

            // Face 2 front
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(1);
            facePtr->setPoint(1,mesh.point(nodeID+N_x));
            facePtr->setPoint(2,mesh.point(nodeID+N_x+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_z));

            nodeID += N_y*(n_y-1);

            // Face 1 back - marker = 3
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(3);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_x));

            // Face 2 back
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(3);
            facePtr->setPoint(1,mesh.point(nodeID+N_x));
            facePtr->setPoint(2,mesh.point(nodeID+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_x+N_z));
        }
    }
    // Faces l_y x l_z
    for(UInt j(0);j<m_z;++j)
    {
        for(UInt i(0);i<m_y;++i)
        {
            nodeID = j*N_z+i*N_y+1;

            // Face 1 left - marker = 4
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(4);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_y));

            // Face 2 left
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(4);
            facePtr->setPoint(1,mesh.point(nodeID+N_y));
            facePtr->setPoint(2,mesh.point(nodeID+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_y+N_z));

            nodeID += N_x*(n_x-1);

            // Face 1 right - marker = 2
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(2);
            facePtr->setPoint(1,mesh.point(nodeID));
            facePtr->setPoint(2,mesh.point(nodeID+N_y));
            facePtr->setPoint(3,mesh.point(nodeID+N_z));

            // Face 2 right
            facePtr = &mesh.addFace(true);
            facePtr->setMarker(2);
            facePtr->setPoint(1,mesh.point(nodeID+N_y));
            facePtr->setPoint(2,mesh.point(nodeID+N_y+N_z));
            facePtr->setPoint(3,mesh.point(nodeID+N_z));
        }
    }
    std::cout<<"done"<<std::endl;

    // Build the volumes
    std::cout<<"building the volumes...";
    UInt volumeID(0);
    UInt nodes[8];
    for(UInt k(0);k<m_z;++k)
    {
        for(UInt j(0);j<m_y;++j)
        {
            for(UInt i(0);i<m_x;++i)
            {
                nodeID = k*N_z+j*N_y+i;
                nodes[0]=nodeID+1;
                nodes[1]=nodeID+N_x+1;
                nodes[2]=nodeID+N_y+1;
                nodes[3]=nodeID+N_x+N_y+1;
                nodes[4]=nodeID+N_z+1;
                nodes[5]=nodeID+N_x+N_z+1;
                nodes[6]=nodeID+N_y+N_z+1;
                nodes[7]=nodeID+N_x+N_y+N_z+1;
                volumeID=(k*m_y*m_x+j*m_y+i)*6 ;

                // Tetra 1
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+1);
                volumePtr->setLocalId(volumeID+1);
                volumePtr->setPoint(1, mesh.point(nodes[0]));
                volumePtr->setPoint(2, mesh.point(nodes[1]));
                volumePtr->setPoint(3, mesh.point(nodes[2]));
                volumePtr->setPoint(4, mesh.point(nodes[4]));
                volumePtr->setMarker(regionFlag);

                // Tetra 2
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+2);
                volumePtr->setLocalId(volumeID+2);
                volumePtr->setPoint(1, mesh.point(nodes[5]));
                volumePtr->setPoint(2, mesh.point(nodes[4]));
                volumePtr->setPoint(3, mesh.point(nodes[2]));
                volumePtr->setPoint(4, mesh.point(nodes[1]));
                volumePtr->setMarker(regionFlag);

                // Tetra 3
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+3);
                volumePtr->setLocalId(volumeID+3);
                volumePtr->setPoint(1, mesh.point(nodes[6]));
                volumePtr->setPoint(2, mesh.point(nodes[4]));
                volumePtr->setPoint(3, mesh.point(nodes[2]));
                volumePtr->setPoint(4, mesh.point(nodes[5]));
                volumePtr->setMarker(regionFlag);

                // Tetra 4
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+4);
                volumePtr->setLocalId(volumeID+4);
                volumePtr->setPoint(1, mesh.point(nodes[7]));
                volumePtr->setPoint(2, mesh.point(nodes[5]));
                volumePtr->setPoint(3, mesh.point(nodes[6]));
                volumePtr->setPoint(4, mesh.point(nodes[3]));
                volumePtr->setMarker(regionFlag);

                // Tetra 5
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+5);
                volumePtr->setLocalId(volumeID+5);
                volumePtr->setPoint(1, mesh.point(nodes[2]));
                volumePtr->setPoint(2, mesh.point(nodes[5]));
                volumePtr->setPoint(3, mesh.point(nodes[3]));
                volumePtr->setPoint(4, mesh.point(nodes[6]));
                volumePtr->setMarker(regionFlag);

                // Tetra 6
                volumePtr = &mesh.addVolume();
                volumePtr->setId(volumeID+6);
                volumePtr->setLocalId(volumeID+6);
                volumePtr->setPoint(1, mesh.point(nodes[1]));
                volumePtr->setPoint(2, mesh.point(nodes[3]));
                volumePtr->setPoint(3, mesh.point(nodes[2]));
                volumePtr->setPoint(4, mesh.point(nodes[5]));
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
}

//@}



} // Namespace LifeV

#endif /* STRUCTUREDMESH3D_H */
