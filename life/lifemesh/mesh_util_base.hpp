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
#ifndef __MESH_UTIL_BASE__
#define __MESH_UTIL_BASE__
#include "lifeV.hpp"
#include <vector>
#include <algorithm>
#include <set>
#include "regionMesh3D.hpp"

/*!
  \brief Base tilities operating on meshes

  
  This file contains a set of base utilities used to test mesh entities or
  operate on them
  
 */

//! A locally used structure, not meant for general use
typedef map<BareFace,pair<ID,ID >,cmpBareItem<BareFace> > TempFaceContainer;

//! A locally used structure, not meant for general use
typedef map<BareEdge,pair<ID,ID>,cmpBareItem<BareEdge> > TempEdgeContainer;

/*
*************************************************************************************
                            FUNCTORS
*************************************************************************************
*/
//! \defgroup Test_Functors Some useful functors to be used for test mesh entities

/*! \ingroup Test_Functors
  \briefFunctor to check if a Point, Face or Edge is on the boundary.
  
  \precond It assumes that boundary points in RegionMesh are correctly set.
  \precond   the RegionMesh must export the typenames
  PointType, FaceType and EdgeType. 
*/
template <typename RegionMesh>
class EnquireBEntity
{
public:
  EnquireBEntity(RegionMesh const & mesh):pmesh(&mesh){};
  typedef typename RegionMesh::FaceType  FaceType;
  typedef typename RegionMesh::EdgeType  EdgeType;
  typedef typename RegionMesh::PointType PointType;
  
  bool operator()(FaceType & face){
    bool isboundary=true;
    for (UInt k=1;k<=FaceType::numVertices;++k){
      isboundary=isboundary & face.point(k).boundary();
    }
    return isboundary;
  }

  bool operator()(EdgeType & edge){
    bool isboundary=true;
    for (UInt k=1;k<=EdgeType::numVertices;++k){
      isboundary=isboundary & edge.point(k).boundary();
    }
    return isboundary;
  }

  INLINE bool operator()(PointType & point){
    return point.boundary();
  }
  
private:
  EnquireBEntity(){}
  RegionMesh const * pmesh;
};

//!\ingroup Test_Functors
/*! Functor to check if a Face is on the boundary, by using the information
 contained in a TempFaceContainer produced by findBoundaryFaces(). It does
 not use the information contained in the mesh PointList, so it differs
 from EnquireBEntity.
 \precond bfaces have been previously set by a call to FindBoundaryFaces
*/
template <typename RegionMesh>
class EnquireBFace
{
public:
  typedef typename RegionMesh::FaceType FaceType;
  typedef typename RegionMesh::FaceShape FaceShape;
  
  EnquireBFace(RegionMesh const & mesh,TempFaceContainer const & bfaces):
    pmesh(&mesh),pbfaces(&bfaces) {}
  
  bool operator()(FaceType & f){
    ID i1,i2,i3,i4;
    BareFace bface;
    
    i1=f.point(1).id();
    i2=f.point(2).id();
    i3=f.point(3).id();
    if (FaceShape::numVertices == 4){
      i4=f.point(4).id();
      bface=(makeBareFace(i1,i2,i3,i4)).first;
    }
    else{
      bface=(makeBareFace(i1,i2,i3)).first;
    }
    return pbfaces->find(bface) != pbfaces->end();
  }
private:
  EnquireBFace(){};
  RegionMesh const * pmesh;
  TempFaceContainer const * pbfaces;
};

//!\ingroup Test_Functors
/*! Functor to check if an edge is on the boundary, by using the information
 contained in a TempFaceContainer produced by findBoundaryEdges(). It does
 not use the information contained in the mesh PointList, so it differs
 from EnquireBEntity.

 \precond bedges have been previously set by a call to FindBoundaryEdges()
 
*/
template <typename RegionMesh>
class EnquireBEdge
{
public:
  typedef typename RegionMesh::EdgeType EdgeType;
  typedef typename RegionMesh::EdgeShape EdgeShape;
  
  EnquireBEdge(RegionMesh const & mesh,TempEdgeContainer const & bedges):
    pmesh(&mesh),pbedges(&bedges) {}
  
  bool operator()(EdgeType & f){
    ID i1,i2;
    BareEdge bedge;
    
    i1=f.point(1).id();
    i2=f.point(2).id();
    bedge=(makeBareEdge(i1,i2)).first;
    return pbedges->find(bedge) != pbedges->end();
  }
  
private:
  EnquireBEdge(){};
  RegionMesh const * pmesh;
  TempEdgeContainer const * pbedges;
};

//! \ingroup Test_Functors
/*! Functor to check if a mesh entity with boundary indicator (for instance a GeoPoint)
  is on the boundary, by enquiring its boundary flag.
  \warning It assumes that boundary points are correctly set.
*/
template <typename RegionMesh>
class EnquireBPoint
{
public:
  EnquireBPoint(RegionMesh & mesh):pmesh(&mesh){};
  bool operator()(MeshEntityWithBoundary & e){
    return e.boundary();
  }
private:
  EnquireBPoint(){};
  RegionMesh * pmesh;
};


//! \ingroup Test_Functors
/*! This functor is used to do some geometry checks It returns a coordinate
 */
class GetCoordComponent
{
public:
  GetCoordComponent();
  GetCoordComponent(int i);
  void operator()(Real const x, Real const y, Real const z, Real ret[3]) const;
private:
  int comp;
};

//! \ingroup Test_Functors
/*! This functor is used to do some geometry checks It returns a vector of ones
 */
class GetOnes
{
public:
void operator()(Real const x, Real const y, Real const z, Real ret[3])const;
};
/*
*************************************************************************************
                            EDGES/FACES FINDERS
*************************************************************************************
*/

//! Finds boundary faces.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the boundary
faces connectivities.

\param mesh A 3D mesh.

\param NumInternalFaces. A reference to an integer returning the number of internal faces found.

\param bfaces This container will eventually contain a map whose key are
the BareFace corresponding to the boundary faces and the data a pair of
IDs: the ID of the adjacent element and the relative position of the face
in the element.

\param allFaces When this bool is set true the function will also construct the set of internale faces, stored in intfaces.

\param intfaces A container that will possibly contain a map whose keys are
the BareFace corresponding to an internal faces and the data a pair of IDs:
the ID of the two elements adjacent to the face.

\return Number of boundary faces found
*/
template <typename RegionMesh3D>
UInt findFaces(const RegionMesh3D  & mesh, TempFaceContainer & bfaces, UInt & numInternalFaces, TempFaceContainer & intfaces,bool allFaces=false){
  UInt i1,i2,i3,i4;
  BareFace bface;
  typename RegionMesh3D::VolumeShape ele;
  typedef typename RegionMesh3D::Volumes Volumes;
  TempFaceContainer::iterator fi;
  
  // clean first in case it has been alredy used
  
  bfaces.clear();
  if(allFaces)intfaces.clear();
  numInternalFaces=0;
  
  for (typename Volumes::const_iterator iv=mesh.volumeList.begin();
       iv != mesh.volumeList.end(); ++iv){
    for (ID j=1;j<=mesh.numLocalFaces();++j){
      i1=ele.fToP(j,1);
      i2=ele.fToP(j,2);
      i3=ele.fToP(j,3);
      // go to global
      i1=(iv->point(i1)).id();
      i2=(iv->point(i2)).id();
      i3=(iv->point(i3)).id();
      if (RegionMesh3D::FaceShape::numVertices == 4){
	i4=ele.fToP(j,4);
	i4=(iv->point(i4)).id();
	bface=(makeBareFace(i1,i2,i3,i4)).first;
      }
      else{
	bface=(makeBareFace(i1,i2,i3)).first;
      }
      
      if( (fi=bfaces.find(bface)) == bfaces.end() ){
	bfaces.insert(make_pair(bface,make_pair(iv->id(),j)));
      } else {
	if(allFaces && i1>i2)intfaces.insert((make_pair(bface,make_pair(iv->id(),j))));
	bfaces.erase(fi);// counted twice: internal face
	++numInternalFaces;
      }
    }
  }
  return bfaces.size();
}

template <typename RegionMesh3D>
UInt findBoundaryFaces(const RegionMesh3D  & mesh, TempFaceContainer & bfaces, UInt & numInternalFaces){
  TempFaceContainer dummy;
  return findFaces(mesh, bfaces, numInternalFaces, dummy, false);
}



//! Finds boundary edges.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the boundary
edges connectivities.

\param mesh A 3D mesh.

\param bedges This container contains a set with the BareEdge of the
boundary edges.

\return Number of boundary edges found.

\pre The list of boundary faces must be correctly set.
*/
template <typename RegionMesh3D>
UInt findBoundaryEdges(const RegionMesh3D & mesh, TempEdgeContainer & bedges){
  UInt i1,i2;
  BareEdge bedge;
  typedef typename RegionMesh3D::FaceShape FaceShape;
  typedef typename RegionMesh3D::Faces Faces;
  

  if (! mesh.hasFaces()) return 0;
  
  // clean first in case it has been alredy used
  bedges.clear();
  
  for (typename Faces::const_iterator ifa=mesh.faceList.begin();
       ifa != mesh.faceList.begin()+mesh.numBFaces(); ++ifa){
    for (ID j=1;j<=mesh.numLocalEdgesOfFace();++j){
      i1=FaceShape::eToP(j,1);
      i2=FaceShape::eToP(j,2);
      // go to global
      i1=(ifa->point(i1)).id();
      i2=(ifa->point(i2)).id();
      bedge=(makeBareEdge(i1,i2)).first;
      bedges.insert(make_pair(bedge,make_pair(ifa->id(),j)));
    }
  }
  return bedges.size();
}

//! Finds all  edges.
/*!  A low level routine, not meant to be called directly. It creates a
container with all the information needed to set up properly the edge connectivities.

\param mesh A 3D mesh.

\param bedges This container contains a set of  BareEdges for all mesh edges.

\return Number of edges found.

*/

template <typename RegionMesh3D>
UInt findInternalEdges(const RegionMesh3D & mesh, const TempEdgeContainer & boundary_edges, TempEdgeContainer & internal_edges){
  UInt i1,i2;
  BareEdge bedge;
  typedef typename RegionMesh3D::ElementShape  VolumeShape;
  typedef typename RegionMesh3D::Volumes Volumes;

  
  ASSERT0(mesh.numVolumes()>0,"We must have some 3D elements stored n the mesh to use this function!");
  
  internal_edges.clear();

  
  for (typename Volumes::const_iterator ifa=mesh.volumeList.begin();
       ifa != mesh.volumeList.end(); ++ifa){
    for (ID j=1;j<=mesh.numLocalEdges();++j){
      i1=VolumeShape::eToP(j,1);
      i2=VolumeShape::eToP(j,2);
      // go to global
      i1=(ifa->point(i1)).id();
      i2=(ifa->point(i2)).id();
      bedge=(makeBareEdge(i1,i2)).first;
      if(boundary_edges.find(bedge)==boundary_edges.end())
	internal_edges.insert(make_pair(bedge,make_pair(ifa->id(),j)));
    }
  }
  return internal_edges.size();
}
/*
*************************************************************************************
                            MARKERS HANDLERS
*************************************************************************************
*/
//! \defgroup marker_handlers Used to manage missing handlers

/*! \ingroup marker_handlers

//! \brief Sets the marker flag of a GeoElement of dimension greater one

 It gets the stronger marker of the GeoElement points. The marker
hierarchy is defined in the marker.h file.  It returns a bool indicating if
the flag has changed. If any of the vertices has an unset marker the result
is an unset flag for the GeoElement.

\warning It overrides the original marker flag.
*/
template <typename GeoElement>
EntityFlag inheritStrongerMarker(GeoElement & fp){
  ASSERT_PRE(GeoElement::nDim>0, "A GeoElement with ndim<1 cannot inherit marker flags");

  fp.setMarker(fp.point(1).marker());
  for (ID j=2;j <= GeoElement::numVertices;++j) fp.setStrongerMarker(fp.point(j).marker());
  return fp.marker();
  
}


/*! \ingroup marker_handlers

//! \brief Sets the marker flag of a GeoElement of dimension greater one

 It gets the weaker marker of the GeoElement points. The marker
hierarchy is defined in the marker.h file.  It returns a bool indicating if
the flag has changed. If any of the vertices has an unset marker the result
is an unset flag for the GeoElement.
  
  \warning It overrides the original marker flag.*/
template <typename GeoElement>
EntityFlag inheritWeakerMarker(GeoElement & fp){
  ASSERT_PRE(GeoElement::nDim>0, "A GeoElement with ndim<1 cannot inherit marker flags");
  fp.setMarker(fp.point(1).marker());
  for (ID j=2;j <= GeoElement::numVertices;++j) fp.setWeakerMarker(fp.point(j).marker());
  return fp.marker();
  
}

//! Fix mesh switches
/*!
  Using some heuristics it tries to fix mesh switches
 */

// template<typename RegionMesh>
// void
// fixSwitches(RegionMesh ^ mesh, ostream & clog=cout, bool verbose=false)
// {
  
//   clog<<" ************** FIXING MESH SWITCHES **********************"<<endl;
//   clog<<"            Mesh switches Status before fixing"<<endl;
//   mesh.showLinkSwitch(verbose,clog);
//   if (mesh.storedFaces()> mesh.numBFaces()){
//     mesh.setLinkSwitch("HAS_ALL_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_FACES");
//   }
//   if (mesh.numBFaces>0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_FACES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDARY_FACES");
//   }
//   if (mesh.storedEdges()> mesh.numBEdges()){
//     mesh.setLinkSwitch("HAS_ALL_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_ALL_EDGES");
//   }
//   if (mesh.numBEdges()> 0){
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }
//   if (mesh.storedFaces()>0){
//     if(mesh.face(1).ad_first(
//     mesh.setLinkSwitch("HAS_BOUNDARY_EDGES");
//   }  else{
//     mesh.unsetLinkSwitch("HAS_BOUNDAY_EDGES");
//   }
  
      

  
  
#endif
