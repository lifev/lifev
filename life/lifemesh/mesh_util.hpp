#ifndef __MESH_UTILITIES__
#define __MESH_UTILITIES__
#include "mesh_util_base.hpp"
#include "geoMap.hpp"
#include "currentFE.hpp"
#include "currentBdFE.hpp"
//! \file mesh_util.h
//! \file mesh_util.h
/*! This file contains a set of functions to be used to test a 3D mesh and
  fix some possible problems. Some methods are general (2D and 3D meshes), some are
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

using namespace std; // To be taken away some time

/*
*****************************************************************************
                                GEOMETRY TESTS
*****************************************************************************
*/
//!  \brief Report 3D element orientation
/*!  It uses a linear representation of the Tetra/Hexa: it is only a
  orientation check.  The orientation is considered positive if it obeys the right-hand rule (right-hand orientation).

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
Real checkVolumes(RegionMesh3D const & mesh, SimpleVect<bool> & elSign, Switch & sw)
{
  Real meas(0.0);
  Real lmeas;
  elSign.clear();
  elSign.reserve(mesh.numVolumes());
  typedef typename RegionMesh3D::VolumeShape GeoShape;
  CurrentFE* fe;
  switch(GeoShape::Shape){
  case TETRA:
    fe = new CurrentFE(feTetraP1,geoLinearTetra,quadRuleTetra1pt);
    for(ID i = 1; i<=mesh.numVolumes(); i++){
      fe->updateJac(mesh.volume(i));
      lmeas = fe->measure();
      meas+=lmeas;
      elSign.push_back(lmeas>0.0);
    }
    break;
  case HEXA:
    fe = new CurrentFE(feHexaQ1,geoBilinearHexa,quadRuleHexa1pt);
    for(ID i = 1; i<=mesh.numVolumes(); i++){
      fe->updateJac(mesh.volume(i));
      lmeas = fe->measure();
      meas+=lmeas;
      elSign.push_back(lmeas>0.0);
    }
    break;
  default:
    sw.create("SKIP_ORIENTATION_TEST",true);

    delete fe;
    return 0;
  }

  if(std::find(elSign.begin(),elSign.end(),false)!=elSign.end())sw.create("HAS_NEGATIVE_VOLUMES",true);

  return meas;
}

/*!
  \brief Fixes  negative volume elements.
  
  Given a vector<bool> indicating negative elements, it inverts those that
  have been found negative.

  \param mesh A 3D mesh. It will be modified.
  
  \param elSign a vector of bools. The value false correspond to the elements that have to be swapped. It is created by
  checkVolumes().

  \post A mesh with all volumes with positive oreintation.
*/
template <typename RegionMesh3D>
void fixVolumes(RegionMesh3D & mesh, const SimpleVect<bool> & elSign, Switch & sw)
{
  typedef typename RegionMesh3D::VolumeShape GeoShape;

  static const ID otn_Tetra[10]={2, 1, 3, 4, 5, 7, 6, 9, 8,10};
  static const ID otn_Hexa[20] ={1, 4, 3, 2, 5, 8, 7, 6,12,11,
				 10, 9,13,16,15,14,20,19,18,17};
  for(ID i = 1; i<=mesh.numVolumes(); i++){
 
    if(! elSign(i)){
      switch(GeoShape::Shape){
      case TETRA:
	mesh.volume(i).exchangePoints(otn_Tetra);
	break;
      case HEXA:
	mesh.volume(i).exchangePoints(otn_Hexa);
	break;
      default:
	sw.create("ABORT_CONDITION",true);
	return;
      }
    }   
  }   
}
//!\brief Computes volume enclosed by boundary faces
/*!  It computes, for $i=1,2,3$, the integral \f$\int_{\partial \Omega} x_i n_i
  d\gamma \f$, \f$n_i\f$ being the i-th component of the boundary normal. If the
  domain boundary is properly disretised they should all return (within
  discretisation and truncation errors) the quantity \f$\vert\Omega\vert\f$.
  
  \warning Not to be used for accurate computations (it always adopts
  linear or bilinear elements, with a simple integration rule) \param mesh
  A 3D mesh \param vols returns 3 Real corresponding to the 3 integrals
*/
template <typename RegionMesh3D>
void getVolumeFromFaces(RegionMesh3D const & mesh, Real vols[3],ostream & err=cerr)
{
  GetCoordComponent getx(0);
  GetCoordComponent gety(1);
  GetCoordComponent getz(2);
  vols[0]=0.0;
  vols[1]=0.0;
  vols[2]=0.0;
  typedef typename RegionMesh3D::FaceShape GeoBShape;
  typedef typename RegionMesh3D::FaceType FaceType;
  CurrentBdFE* bdfe;
  switch(GeoBShape::Shape){
  case TRIANGLE:
    bdfe=new CurrentBdFE(feTriaP1,geoLinearTria,quadRuleTria1pt);
    for(ID i = 1; i<=mesh.numBFaces(); i++){
      bdfe->updateMeasNormal(mesh.face(i));
      vols[0]+=bdfe->integral_n(getx);
      vols[1]+=bdfe->integral_n(gety);
      vols[2]+=bdfe->integral_n(getz);
    }
    break;
  case QUAD:
    bdfe=new CurrentBdFE(feQuadQ1,geoBilinearQuad,quadRuleQuad1pt);    
    for(ID i = 1; i<=mesh.numBFaces(); i++){
      bdfe->updateMeasNormal(mesh.face(i));
      vols[0]+=bdfe->integral_n(getx);
      vols[1]+=bdfe->integral_n(gety);
      vols[2]+=bdfe->integral_n(getz);
    }
    break;
   default:
     err<< "Only tria and quad surface elements  may be checked for volume orientation at the moment"<<endl;
     ASSERT0(false,"ABORT CONDITION OCCURRED");
  }
}

//! Tests if the surface of the mesh is closed by computing surface integrals.
/*! It computes \f$\sum_{i=1}^3\int_{\partial \Omega} n_i d\gamma\f$.
  The value returned  should be very proximal to zero
 */
template <typename RegionMesh3D>
Real testClosedDomain(RegionMesh3D const & mesh,ostream & err=cerr)
{
  typedef typename RegionMesh3D::FaceType FaceType;

  CurrentBdFE* bdfe;
  
  GetOnes ones;
  Real test(0.0);
  
  switch(RegionMesh3D::FaceShape::Shape){
  case TRIANGLE:
    bdfe=new CurrentBdFE(feTriaP1,geoLinearTria,quadRuleTria1pt);
    for(ID i = 1; i<=mesh.numBFaces(); i++){
      bdfe->updateMeasNormal(mesh.face(i)); 
      test+=bdfe->integral_n(ones);
    }
    break;
    case QUAD:
    bdfe=new CurrentBdFE(feQuadQ1,geoBilinearQuad,quadRuleQuad1pt);
    for(ID i = 1; i<=mesh.numBFaces(); i++){
      bdfe->updateMeasNormal(mesh.face(i)); 
      test+=bdfe->integral_n(ones);
    }
    
     break;
     
   default:
     err<< "Only tria and quad surface elements  may be checked for volume orientation at the moment"<<endl;
     ASSERT0(false,"ABORT CONDITION OCCURRED");
  }
  return test;
  
}
//! \brief Tests if the surface of the mesh is closed by computing surface integrals.
/*!
  This routine tests if the topological descrption of boundary face is sane. In particular
  all boundary edges must be adjacent to only 2 surface elements and the orientation must be correct.

  \param mesh a mesh
  \param numBedges  The function also returns the number of boundary edges in numBedges.
  \return It it returns 0 the test has been passed. If  not it returns the number of of wrong boundary edges.
  \warning numBEdges is properly set only if the test has been passed. 
*/
template <typename RegionMesh3D>
UInt testClosedDomain_Top(RegionMesh3D const & mesh, UInt & numBEdges){

  typedef set<BareEdge,cmpBareItem<BareEdge> > TempEdgeContainer2;
  TempEdgeContainer2 bedges;
  UInt i1,i2;
  BareEdge bedge;
  typename RegionMesh3D::BElementShape ele;
  typedef typename RegionMesh3D::Faces Faces;
  typedef typename RegionMesh3D::FaceType FaceType;
  TempEdgeContainer2::iterator ed;

  
  // clean first in case it has been alredy used
  
  typename Faces::const_iterator iv=mesh.faceList.begin();
  
  for (UInt k=0;k<mesh.numBFaces();++k){
    ASSERT(iv !=mesh.faceList.end()," Trying to get not existing face"<<k<<" "<<mesh.numBFaces());
    
    for (ID j=1;j<=FaceType::numEdges;++j){
      i1=ele.eToP(j,1);
      i2=ele.eToP(j,2);
      // go to global
      i1=(iv->point(i1)).id();
      i2=(iv->point(i2)).id();
      bedge=(makeBareEdge(i1,i2)).first;
      
      if( (ed=bedges.find(bedge)) == bedges.end() ){
	bedges.insert(bedge);
	++numBEdges;
      } else {
	bedges.erase(ed);
      }
    }
    ++iv;
  }
  return bedges.size();
}
/*
*****************************************************************************
                                MARKERS FIXING
*****************************************************************************
*/

//! Check wether all markers of a the goemetry entities stored in a list are set
template <typename MeshEntityList>
bool checkMarkerSet(const MeshEntityList & list)
{
  typedef typename MeshEntityList::const_iterator C_Iter;
  bool ok(true);
  for (C_Iter l=list.begin();l!=list.end();++l)ok=(ok & l->isMarkerSet());
  return ok;
}

//! Sets the marker flag for all boundary edges by inheriting them from boundary points.
/*! The paradigm is that an edge <B>WHOSE MARKER HAS NOT ALREADY BEEN
  SET</B> will get the WEAKER marker flag among its VERTICES. For instance
  is a vertex is assigned to an Essential B.C and the other to a Natural
  B.C. the edge will get the flag related to the Natural B.C.

  /param mesh A mesh
  /param clog ostream to which the logging of the map of the newly assigned marked will be output
  /param err ostream to which error messages will be sent 
  
  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/

template <typename RegionMesh>
void
setBEdgesMarker(RegionMesh  & mesh, ostream & clog=cout, ostream & err=cerr, bool verbose=true)
{
  typename RegionMesh::EdgeType *  fp=0;
  unsigned int count(0);
  
  if(verbose) clog<<"NEW EDGE MARKER MAP"<<endl<<" ID->New Marker"<<endl;
  
  for (ID k=1; k<=mesh.numBEdges(); ++k){
    fp = &(mesh.edge(k));
    if (fp->isMarkerUnset()){
      inheritWeakerMarker(*fp);
      if (verbose){
	clog<<fp->id()<<" -> ";
	fp->printFlag(clog);
	clog<<" ";
	if (++count % 3 ==0) clog<<endl;
      }
    } 
  }
  if (verbose) clog<<endl;
}


//! Sets the marker flag for all boundary faces by inheriting them from boundary points.
/*! The paradigm is that a face WHOSE MARKER HAS NOT ALREADY BEEN SET will
  get the WEAKER marker flag among its VERTICES. For instance if a vertex
  is assigned to a Natural B.C and the others to a Natural B.C. the face
  will get the flag related to the Natural B.C.
  
  /todo better handling of flags: all function handling flags should be
  wrapped into a class
*/
template <typename RegionMesh>
void
setBFacesMarker(RegionMesh  & mesh, ostream & clog=cout, ostream & err=cerr,bool verbose=true)
{
  typename RegionMesh::FaceType *  fp=0;
  unsigned int count(0);
  
  if(verbose) clog<<"NEW FACE MARKER MAP"<<endl<<" ID->New Marker"<<endl;
  
  for (UInt k=1;k<=mesh.numBFaces();++k){
    fp = &(mesh.face(k));
    if (fp->isMarkerUnset()){
      inheritWeakerMarker(*fp);
      if (verbose){
	clog<<fp->id()<<" -> ";
	fp->printFlag(clog);
	clog<<" ";
	if (++count % 3 ==0) clog<<endl;
      }
    }
  }
  if (verbose) clog<<endl;
}

//! It sets the marker flag of boundary points, by inheriting it from boundary elements.
/*! The paradigm is that a point whose marker flag is unset will inherhit
  the strongest marker flag of the surrounding Boundary Elements, with the
  convention that if the marker flag of one of the surrounding boundary
  elements is null is ignored.
*/
template <typename RegionMesh>
void
setBPointsMarker(RegionMesh & mesh, ostream & clog=cout, ostream& err=cerr, bool verbose=false)
{  
  // First looks at points whose marker has already been set
  vector<bool> markset(mesh.storedPoints(),false);
  
  typedef typename RegionMesh::Points::iterator PointIterator;
  typedef typename RegionMesh::BElementShape BElementShape;
  
  vector<bool>::iterator pm=markset.begin();
  
  for (PointIterator p=mesh.pointList.begin();p!=mesh.pointList.end();++p)*(pm++)=p->isMarkerSet();
  
  typename RegionMesh::BElementType *  fp=0;
  for (UInt k=1;k<=mesh.numBElements();++k){
    fp = &(mesh.bElement(k));
    if (fp->isMarkerSet()){
      for (UInt j=1;j<=BElementShape::numPoints;++j){
	if(!markset[(fp->point(j).id())-1])
	  fp->point(j).setStrongerMarker(fp->marker());
      }
    }
  }
  unsigned int count(0);
  if (verbose){
    clog<<"**** NEW POINTS MARKERS **************"<<endl;
    clog<<"id->marker    id->marker     id->marker"<<endl;
    pm=markset.begin();
    for (PointIterator p=mesh.pointList.begin();p!=mesh.pointList.end();++p){
      if (*pm++){
	clog<<p->id()<<" -> ";
	p->printFlag(clog);
	clog<<" ";
	if(++count % 3)clog<<endl;
      }
    }
  }
}
/*
*****************************************************************************
                                FIXING ID AND COUNTERS 
*****************************************************************************
*/
//! \brief Verifies if a list of mesh entities hav ethe ID properly set.
/* More precisely, the id() must correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.
   
   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method id().
*/
template <typename MeshEntityList>
bool checkIdnumber(const MeshEntityList & list)
{
  typedef typename MeshEntityList::const_iterator C_Iter;
  bool ok(true);
  unsigned int count(1);
  for (C_Iter l=list.begin();l!=list.end() && ok;++l, ++count)ok=(l->id()==count);
  return ok;
}

//! \brief Fixes a a list of mesh entities so that the ID is properly set.
/* \post  The id will correspond to the position of the entity
   in the list (starting from 1, since id=0 is reserved for unset entities.
   
   \pre The template argument MeshEntityList must be a stl
   compliant container and its elements must have the method UInt &id().
*/
template
<typename MeshEntityList>
void fixIdnumber(MeshEntityList & list)
{
  unsigned int count(0);
  typedef typename MeshEntityList::iterator Iter;
  for (Iter l=list.begin() ;l != list.end(); ++l)l->id()=++count;
}

/*! \brief Fixes boundary points counter
  It fix the boundary points counter by counting
  how many points have te boundary flag set.
  It also reset the Bpoints list.
  
  \pre It assumes that the points have the boundary flag corretly set
*/

template <typename RegionMesh>
void
setBPointsCounters(RegionMesh & mesh){

  unsigned int countBP(0);
  unsigned int countBV(0);

  mesh._bPoints.clear();
  
  for (UInt k=1;k<=mesh.numVertices();++k){
    if(mesh.isBoundaryPoint(k)){
      ++countBP;
      ++countBV;
    }
  }

  for (UInt k=mesh.numVertices()+1;k<=mesh.numPoints();++k){
    if(mesh.isBoundaryPoint(k)){
      ++countBP;
    }
  }

  mesh.numBVertices()=countBV;
  mesh.setNumBPoints(countBP);
  mesh._bPoints.reserve(countBP);

  for (UInt k=1;k<=mesh.numPoints();++k){
    if(mesh.isBoundaryPoint(k))mesh._bPoints.push_back(&mesh.point(k));
  }
}

/*
*****************************************************************************
                                BOUNDARY INDICATOR FIXING
*****************************************************************************
*/
//! It fixes boundary flag on points laying on boundary faces.
/*!
  \param mesh a mesh
  \param clog logging stream
  \param cerr error stream
  \param verbose If true you have a verbose output
  
  \pre mesh point list must exists and boundary face lsist  must have been set properly.
*/
template <typename RegionMesh>
void
fixBPoints(RegionMesh  & mesh, ostream & clog=cout,ostream & err=cerr, bool verbose=true){
  ASSERT_PRE(mesh.numPoints()>0,"The point list should not be empty");
  ASSERT_PRE(mesh.numBElements()>0,"The BElements list should not be empty");
  
  typedef typename RegionMesh::BElements BElements;
  typedef typename RegionMesh::BElementShape BElementShape;
  typename RegionMesh::BElementType * fp;

  if (verbose) clog<<"New BPoints Found "<<endl;
  for (UInt k=1;k<=mesh.numBElements();++k){
    fp = &(mesh.bElement(k));
    for (UInt j=1;j<=BElementShape::numPoints;++j){
      if (verbose && !fp->point(j).boundary() ) clog<<"ID: "<<fp->point(j).id()<<endl;
      fp->point(j).boundary()=true;
    }
  }
  // Fix now the number of vertices/points
  setBPointsCounters(mesh);
}

//!It makes sure that boundary edges are stored first
/*!
\pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
bool setBoundaryEdgesFirst(RegionMesh & mesh){

  typedef typename RegionMesh::Edges Edges;
  // set the functor
  EnquireBEntity<RegionMesh > enquireBEdge(mesh);

  std::partition(mesh.edgeList.begin(),mesh.edgeList.end(),enquireBEdge);
  fixIdnumber(mesh.edgeList);
}

//!It makes sure that boundary faces are stored first
/*!
\pre It assumes that boundary points are properly stored in the mesh
*/
template <typename RegionMesh>
bool setBoundaryFacesFirst(RegionMesh & mesh){
  
  typedef typename RegionMesh::Faces Faces;
  // set the functor
  EnquireBEntity<RegionMesh> enquireBFace(mesh);

  std::partition(mesh.faceList.begin(),mesh.faceList.end(),enquireBFace);
  fixIdnumber(mesh.faceList);
  
}

//! Tests if boundary faces are stored first
/*! \return true if boundary faces are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryFacesFirst(const RegionMesh & mesh){
  
  typedef typename RegionMesh::Faces Faces;
  
  // set the functor
  EnquireBEntity<RegionMesh> enquireBFace(mesh);
  typename RegionMesh::FaceType *  fp;
  bool ok(true);
  
  for (UInt k=1;k<=mesh.numBElements();++k) ok=ok && enquireBFace(mesh.boundaryFace(k));
  for (UInt k=mesh.numBElements()+1;k<=mesh.storedFaces();++k) ok=ok && ! enquireBFace(mesh.face(k));

  return ok;
}

//! Tests if boundary edges are stored first
/*! \return true if boundary edges are indeed stored first
  \pre It assumes that boundary points are set */
template <typename RegionMesh>
bool checkBoundaryEdgesFirst(const RegionMesh & mesh){
  
  typedef typename RegionMesh::Edges Edges;

  // set the functor
  EnquireBEntity<RegionMesh> enquireBEdge(mesh);
  typename RegionMesh::EdgeType *  fp;
  bool ok(true);
  
  for (UInt k=1;k<=mesh.numBEdges();++k) ok=ok && enquireBEdge(mesh.boundaryEdge(k));
  for (UInt k=mesh.numBEdges()+1;k<=mesh.storedEdges();++k) ok=ok && ! enquireBEdge(mesh.edge(k));
  return ok;
}

/*
*****************************************************************************
 UTILITIES TO VERIFY/CREATE FACES/EDGES                                 
*****************************************************************************
*/
//! It fixes boundary faces so that they are consistently numbered with volumes.

/*! An important step for building degrees of freedom on faces.  It also
  fixes other face related data.
\param mesh a mesh
\param err  ostream for error messages
\param sw A switch that will contain information on what has been done
Possible values are
<ol>
<li>NUM_FACES_MISMATCH</li>
<li>FIXED_FACE_COUNTER</li>
<li>BFACE_MISSING</li>
<li>BFACE_STORED_MISMATCH</li>
<li>BELEMENT_COUNTER_UNSET</li>
<li>BFACE_STORED_MISMATCH</li>
<li>FIXED_MAX_NUM_FACES</li>
</ol>

\param fixMarker If set to the true value all faces without a markerFlag set will inherit it from the points.

\param clog ostream that will all information regarding the markers

\param verbose if falso nothng is written to clog

\param numFaces It returns the number of faces found by the function

\param bfaces_found It returns the number of boundary faces found by the function

\param ext_container. If not NULL it is a pointer to an external map of bondary faces, already
  produced by a call to findBoundaryFaces(). This parameter may be used to save al lot of computational work, since
  findBoundaryFaces() is rather expensive.

\pre Boundary faces list must be properly set.
*/

template <class RegionMesh3D>
bool fixBoundaryFaces(RegionMesh3D & mesh, 
		      ostream & clog, ostream &err, Switch & sw, 
		      UInt & numFaces, UInt & bfaces_found,
		      bool fixMarker=false, bool verbose=false,
		      TempFaceContainer * ext_container)
{

  typedef typename RegionMesh3D::Volumes Volumes;
  typedef typename RegionMesh3D::VolumeType VolumeType;
  typedef typename RegionMesh3D::Faces Faces;
  typedef typename RegionMesh3D::FaceType FaceType;

  UInt i1,i2,i3,i4;
  BareFace bface;
  VolumeType * pv;
  typename Faces::iterator fit;
  typename RegionMesh3D::VolumeShape ele;
  TempFaceContainer * bfaces;
  TempFaceContainer::iterator fi;
  pair<ID,ID>info;
  ID j;
  ID vol;
  UInt numInternalFaces;
  bool notfound(false);
  bool extcont(false);

  if(extcont=(ext_container != 0 )){
    bfaces=ext_container;
    bfaces_found=bfaces->size();
  } else {
    bfaces=new TempFaceContainer;
    bfaces_found=findBoundaryFaces(mesh, *bfaces,numInternalFaces);
    numFaces=bfaces_found+numInternalFaces;
  }
  
  
  bool notEnough=mesh.storedFaces()<bfaces_found;
  

  
  if (notEnough){
    err << "WARNING: number of B. Faces stored smaller"<<endl;
    err << "than the number of bfaces found  and build is not set"<<endl;
    err << "POSSIBLE ERROR"<<endl;
    sw.create("BFACE_STORED_MISMATCH",true);
  }
  
  if (mesh.numBElements()==0){
    err<<"ERROR: Boundary Element counter was not set"<<endl;
    err<<"I Cannot proceed because the situation is ambiguous"<<endl;
    err<<"Please check and eventually either: (a) call buildBoundaryFaces()"<<endl;
    err<<"or (b) set the correct number of bfaces in the mesh using mesh.numBElements()"<<endl;
    err<<"ABORT";
    sw.create("BELEMENT_COUNTER_UNSET",true);
  }
  
  if(mesh.numBFaces()!=bfaces_found){
    err<<"WARNING: B Face counter in mesh is set to "<< mesh.numBFaces();
    err<<"While I have found "<< bfaces_found<<" B. Elements in mesh"<<endl;
    err<<"Plese check... I continue anyway"<<endl;
    sw.create("BFACE_COUNTER_MISMATCH",true);
  }
  
  if (verbose){
    clog<<"**** Marker Flags for Fixed Boundary Faces ***"<<endl;
    clog<<" (it only contains those that were fixed because unset !"<<endl;
    clog<<"id->marker   id->marker  id->marker"<<endl;
  }
  
  UInt count(0);
  
  fit=mesh.faceList.begin();
  for (UInt facid=0;facid<mesh.numBElements();++facid){
    i1=(fit->point(1)).id();
    i2=(fit->point(2)).id();
    i3=(fit->point(3)).id();
    if (RegionMesh3D::FaceShape::numVertices == 4){
      i4=(fit->point(4)).id();
      bface=(makeBareFace(i1,i2,i3,i4)).first;
    }
    else{
      bface=(makeBareFace(i1,i2,i3)).first;
    }
    fi=bfaces->find(bface);
    if (fi==bfaces->end()){
      notfound=true;
    } else {
      info=fi->second;
      vol=info.first; // Element ID
      pv=&mesh.volume(vol);// Element 
      j=info.second;       // The local ID of face on element
      // Reset face point definition to be consistent with face.
      for (UInt k=1;k<=FaceType::numPoints;++k){
	fit->setPoint(k,pv->point(ele.fToP(j,k)));
      }
      // Correct extra info
      fit->ad_first()=vol;
      fit->pos_first()=j;
      if (fit->markerUnset()){
	inheritWeakerMarker(*fit);
	if (verbose){
	  clog<<fit->id()<<" -> ";
	  fit->printFlag(clog);
	  clog<<" ";
	  if (++count % 3 ==0) clog<<endl;
	}
      }
      // Take out face from temporary container
      bfaces->erase(fi);
    }
    ++fit;
  }

  if(!extcont) delete bfaces;
  
  if(notfound){
    err<<"WARNING: At least one boundary face has not been found on the list stored in RegionMesh3D\n";
    sw.create("BFACE_MISSING",true);
  }
  
  if(verbose){
    clog<<endl<<"  *****  END OF LIST ****"<<endl;
  }
  
  // Here I calculate the number of faces,
  
  if (mesh.numFaces()!=numFaces){
    err<<"WARNING: faces counter in mesh  should be "<<numFaces<<endl;
    err<<"(bfaces->size()+numInternalFaces)"<<endl;
    err<<"it is instead "<<mesh.numFaces();
    sw.create("NUM_FACES_MISMATCH",true);
  }
    mesh.setLinkSwitch(string("HAS_BOUNDARY_FACES"));

    return true;
}

//! Builds faces
/*! This function may alternatively be used to build the compulsory boundary faces, all the mesh faces, or just add to an
  existing list of just boundary faces the internal ones.

  \param mesh A mesh
  
  \param clog Log file for information on the newly created markers

  \param err  Error stream

  \param buildbounary if true the function builds boundary faces

  \param buildinternal if true the function builds internal faces

  \param verbose. If true markerFrlags info is written on clog.

  \param numInternalFaces It returns the number of internal faces (only if ext_container is not provided!)
  
  \param bfaces_found It returns the number of boundary faces
  
  \param ext_container. If not NULL it is a pointer to an external map of bondary faces, already
  produced by a call to findBoundaryFaces(). This parameter may be used to save al lot of computational work, since
  findBoundaryFaces() is rather expensive.

  \pre If buildinternal=true and buildboundary=false the mesh must contain a proper list
  of boundary faces

  \note By setting buildinternal=true and buildboundary=true the function just fixes the counters
  with the number of faces in the mesh
 */
template <class RegionMesh3D>
bool buildFaces(RegionMesh3D & mesh, 
		ostream & clog, ostream &err, UInt & bfaces_found,
		UInt & numInternalFaces, 
		bool buildboundary=true,
		bool buildinternal=false,
		bool verbose=false,
		TempFaceContainer * ext_container=0)
{
  UInt i1,i2,i3,i4;
  typename RegionMesh3D::VolumeShape ele;
  typedef typename RegionMesh3D::Volumes Volumes;
  typedef typename RegionMesh3D::VolumeType VolumeType;
  typedef typename RegionMesh3D::Faces Faces;
  typedef typename RegionMesh3D::FaceType FaceType;
  VolumeType * pv;
  TempFaceContainer*  bfaces;
  TempFaceContainer::iterator fi;
  bool extcont(false);
  
  pair<ID,ID>info;
  ID j,id;
  ID vol;

  if(extcont=(ext_container!=0)){
    bfaces=ext_container;
    bfaces_found=bfaces->size();
  } else {
    bfaces= new TempFaceContainer;
    bfaces_found=findBoundaryFaces(mesh, *bfaces,numInternalFaces);
  }
  
  if(buildboundary)mesh.faceList.clear();
  mesh.setNumBFaces(bfaces_found);
  if (!buildinternal){
    mesh.setMaxNumFaces(bfaces_found,false);
    mesh.numFaces()=numInternalFaces+bfaces_found;
  }
  else{
    mesh.setMaxNumFaces(numInternalFaces+bfaces_found,true);
  }
  
  FaceType face;
  
  if(buildboundary){
    
    if (verbose){
      clog<<"**** Marker Flags for Newly Created Boundary Faces ***"<<endl;
      clog<<"id->marker   id->marker  id->marker"<<endl;
    }
    
    for(fi=bfaces->begin();fi!=bfaces->end();++fi){
      info=fi->second;
      vol=info.first; // Element ID
      pv=&mesh.volume(vol);// Element 
      j=info.second;       // The local ID of face on element
      
      for (UInt k=1;k<=FaceType::numPoints;++k)face.setPoint(k,pv->point(ele.fToP(j,k)));
      // Add extra info
      face.ad_first()=vol;
      face.pos_first()=j;
      // Get marker value
      inheritWeakerMarker(face);
      id=mesh.addFace(face,true).id();
      if(verbose){
	if (id %3 == 0) clog<<endl;
	clog<<id<<" -> ";
	face.printFlag(clog);
	clog<<" ";
      }
    }
    mesh.setLinkSwitch(string("HAS_BOUNDARY_FACES"));
    if (! buildinternal)mesh.unsetLinkSwitch(string("HAS_ALL_FACES"));
    mesh.setLinkSwitch(string("FACES_HAVE_ADIACENCY"));
  }
  
  if (!extcont) delete bfaces;
  
  if (! buildinternal ) return true;

  
  if(!buildboundary){
    if (mesh.storedFaces()<mesh.numBFaces()){
      err<< "ERROR: mesh has not boundary faces, cannot just create internal ones!!!"<<endl;
      err<<"ABORT CONDITION"<<endl;
      return false;
    }
    else if (mesh.storedFaces()>mesh.numBFaces()){
      mesh.faceList.resize(mesh.numBFaces());
    }
  }
  
  
  // I may get rid of the bfaces container. Unfortunately now I need a more complex structure, a BareItemsHandel,
  // in order to generate the internal faces id. An alternative would be to use the point data to identify boundary faces
  // as the ones with all point on the boundary. Yet in this function we do not want to use a priori infromation, so that
  // it might work even if the points boundary flag is not properly set.

  
  BareItemsHandler<BareFace> _be;
  pair<UInt,bool> e;
  pair<BareFace,bool> _face;
  
  for (UInt j=0; j<mesh.faceList.size();++j){
    i1=(mesh.faceList[j].point(1)).id();
    i2=(mesh.faceList[j].point(2)).id();
    i3=(mesh.faceList[j].point(3)).id();
    if (RegionMesh3D::FaceShape::numVertices == 4){
      i4=(mesh.faceList[j].point(4)).id();
      _face=makeBareFace(i1,i2,i3,i4);
    }
    else{
      _face=makeBareFace(i1,i2,i3);
    }
    _be.addIfNotThere(_face.first);
  }

  EntityFlag mm(mesh.marker());
  ID vid;
  
  for (typename Volumes::iterator iv=mesh.volumeList.begin();
       iv != mesh.volumeList.end(); ++iv){
    vid=iv->id();
    // REMEMBER: numbering from 1
    for (UInt j=1;j<=mesh.numLocalFaces();j++){
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
	_face=makeBareFace(i1,i2,i3,i4);
      }
      else{
	_face=makeBareFace(i1,i2,i3);
      }
      e=_be.addIfNotThere(_face.first);
      if(e.second)
	{
	  // a new face It must be internal. 
	  for (UInt k=1;k<=FaceType::numPoints;++k) face.setPoint(k,iv->point(ele.fToP(j,k)));
	  face.ad_first()=vid;
	  face.pos_first()=j;
	  // gets the marker from the RegionMesh
	  face.setMarker(mm);
	  mesh.addFace(face,false); //The id should be correct 
	}
      else
	{
	  if(e.first>bfaces_found) // internal
	    {
	      mesh.faceList(e.first).ad_second()=vid;
	      mesh.faceList(e.first).pos_second()=j;
	    }	      
	}
    }
  }
  mesh.setLinkSwitch(string("HAS_ALL_FACES"));
  return true;
}
  
//! It builds edges.

/*! This function may alternatively be used to build the boundary edges, all the mesh faces, or just add the internal edges
  to an existing list of just boundary edges.

  \param mesh A mesh

  \param clog Log file for information on the newly created markers for boundary edges

  \param err  Error stream

  \param bedges_found Returns the number of boundary edges
  
  \param iedges_found Returns the number of internal edges
  
  \param buildbounary if true the function builds boundary edges

  \param buildinternal if true the function builds internal edges

  \param verbose. If true markerFlags info is written on clog.

  \param ext_container. If not NULL it is a pointer to an external map of bondary edges, already
  produced by a call to findBoundaryEdges(). This parameter may be used to save al lot of computational work, since
  findBoundaryEdges() is rather expensive.
  
  \pre If buildinternal=true and buildboundary=false the mesh must contain a proper list
  of boundary edges
  \pre The mesh must copntain a proper list of boundary faces

  \note By setting buildinternal=true and buildboundary=true the function just fixes the counters
  with the number of edges in the mesh
 */

template <typename RegionMesh3D>
bool buildEdges(RegionMesh3D & mesh, ostream & clog, ostream &err, UInt & bedges_found,
		UInt & iedges_found, bool buildboundary=true, bool buildinternal=false, bool verbose=false,
		TempEdgeContainer * ext_container=0)
{
  typedef typename RegionMesh3D::Volumes Volumes;
  typedef typename RegionMesh3D::Faces Faces;
  typedef typename RegionMesh3D::VolumeType VolumeType;
  typedef typename RegionMesh3D::VolumeShape VolumeShape;
  typedef typename RegionMesh3D::Edges Edges;
  typedef typename RegionMesh3D::EdgeType EdgeType;
  typedef typename RegionMesh3D::FaceType FaceType;
  typedef typename RegionMesh3D::FaceShape FaceShape;
  typename RegionMesh3D::FaceType * pf;
  typename RegionMesh3D::VolumeType * pv;
  
  TempEdgeContainer * bedges;
  TempEdgeContainer iedges;
  pair<ID,ID>info;
  ID j,id;
  ID facID;


  bool extcont(false);
  

  if (extcont=(ext_container != 0)){
    bedges=ext_container;
    bedges_found=bedges->size();
  } else {
    bedges = new TempEdgeContainer;
    bedges_found=findBoundaryEdges(mesh,*bedges);
  }
  
  iedges_found= findInternalEdges(mesh,*bedges,iedges);
  // free some memory if not needed!
  if(!buildinternal)iedges.clear();
  if(!buildboundary && buildinternal){
    if (mesh.storedEdges()<bedges_found){
      err<<"ERROR in buildedges(): mesh does not contain boundary edges"<<endl;
      err<<"I need to set buildboundary=true"<<endl;
      err<<"ABORT CONDITION"<<endl;
      return false;
    }
    else if (mesh.storedEdges()>bedges_found){
      mesh.edgeList.resize(bedges_found);
    }
  }
  mesh.setNumBEdges(bedges_found);
  mesh.numEdges()=(bedges_found+iedges_found);
  if(buildboundary)mesh.edgeList.clear();  
  if(buildboundary && ! buildinternal) mesh.setMaxNumEdges(bedges_found,false);
  if(buildinternal) mesh.setMaxNumEdges(bedges_found+iedges_found,true);
  err<< "Building edges from scratch"<<endl;
  
  EdgeType edge;

  if(buildboundary){
    
    if (verbose){
      clog<<"**** Marker Flags for Newly Created Boundary Edges ***"<<endl;
      clog<<"id->marker   id->marker   id->marker"<<endl;
    }
    
    // First boundary.
    for(  TempEdgeContainer::iterator ei=bedges->begin();
	  ei!=bedges->end();++ei){
      info=ei->second;
      facID=info.first; // Face ID
      pf=&mesh.face(facID);// Face
      j=info.second;       // The local ID of edge on face
      for (UInt k=1;k<=EdgeType::numPoints;++k){
	edge.setPoint(k,pf->point(FaceShape::eToP(j,k)));
      }
      
      inheritWeakerMarker(edge);// Get marker value inheriting from points
      
      id=mesh.addEdge(edge,true).id();
      if(verbose){
	if (id %3 == 0) clog<<endl;
	clog<<id<<" -> ";
	edge.printFlag(clog);
	clog<<" ";
      }
    }
    
    if (verbose) clog<<endl<<"  *****  END OF LIST OF BOUNDARY EDGES ****"<<endl;
    
    mesh.setLinkSwitch(string("HAS_BOUNDARY_EDGES"));
  }

  if (!extcont) delete bedges;
  
  if(!buildinternal){
    mesh.unsetLinkSwitch(string("HAS_ALL_EDGES"));
    return true;
  }
  

  
  // Now internal edges
  // free some memory
  
  for(  TempEdgeContainer::iterator ei=iedges.begin();
	ei!=iedges.end();++ei){
    info=ei->second;
    facID=info.first; // Volume ID
    pv=&mesh.volume(facID);// Volume that generated the edge
    j=info.second;       // The local ID of edge on volume
    for (UInt k=1;k<=EdgeType::numPoints;++k)
      edge.setPoint(k,pv->point(VolumeShape::eToP(j,k)));
    edge.setMarker(mesh.marker()); // Get marker value: that of the mesh
    mesh.addEdge(edge,false);
  }
  
  mesh.setLinkSwitch(string("HAS_ALL_EDGES"));

  return true;
}


/*
*****************************************************************************
 UTILITIES TO TRANSFORM A MESH
*****************************************************************************

*/
//! It builds a P2 mesh from P1 data.
/*! \author L.Formaggia.
  \version Version 1.0
  \pre All compulsory structures in mesh must have been already set: volumes and boundary faces.
  \pre Points list MUST have been dimensioned correctly!!!
  \note the function takes advantage of the fact that 
*/ 
template<typename RegionMesh> void
p1top2(RegionMesh & mesh, ostream & out=std::cout){ 

  typedef typename RegionMesh::ElementShape GeoShape;
  typedef typename RegionMesh::BElementShape GeoBShape;
  ASSERT_PRE(GeoShape::numPoints > 4, "p1top2 ERROR: we need a P2 mesh"); 
  
  out << "Building P2 mesh points and connectivities from P1 data" <<endl;
  

  typename RegionMesh::PointType * pp=0;
  typename RegionMesh::EdgeType * pe=0;
  typename RegionMesh::ElementType * pv=0;
  typename RegionMesh::BElementType * pbe=0;
  typedef typename RegionMesh::Elements Elements;
  typedef typename RegionMesh::BElements BElements;
  
  BareItemsHandler<BareEdge> _be;
  pair<UInt,bool> _edgeid;
  UInt i1,i2,e_id;
  pair<BareEdge,bool> _edge;
  typename RegionMesh::ElementShape ele;
  out<<"Processing "<< mesh.storedEdges()<<" P1 Edges"<<endl;
  UInt nbe=mesh.numBEdges();
  for (UInt j=1; j<=mesh.storedEdges();++j){
    pe=& mesh.edge(j);
    i1=(pe->point(1)).id();
    i2=(pe->point(2)).id();
    pp=& mesh.addPoint(j<=nbe); // true for boundary points
    pp->x()=((pe->point(1)).x()+
	     (pe->point(2)).x())*.5;
    pp->y()=((pe->point(1)).y()+
	     (pe->point(2)).y())*.5;
    pp->z()=((pe->point(1)).z()+
	     (pe->point(2)).z())*.5;
    
    // If we have set a marker for the boundary edge, that marker is inherited by the new created point
    // Otherwise the edge (and the new created point) gets the WORST marker among the two end Vertices
    /*
      JFG 07/2002:
      if the mesh file do not contain the edges (inria files), they are built in
      fixBoundaryEdges(...), but I suspect that this function does not attribute the right
      marker to the edges (maybe a problem in setWorseMarkerOfEntity, or something like that...)
      If you do #undef JFG : no change, if you do #define JFG, we do not consider the (wrong) marker
      of the edge to define the marker of the added node (I arbitrarily take the marker of the first
      node)
     */
    //#define JFG
    //#ifndef JFG
    // original version: DOES NOT work when the edges are not give in the mesh file
    if (pe->markerUnset())inheritWeakerMarker(*pe);
    pp->setMarker(pe->marker());
    //#else
    // temporary version that works when the edges are not given in the mesh file
    // pe->setMarker(mesh.point(i1).marker());
    // pp->setMarker(pe->marker());
    //#endif
    // pp->id()=++i; // Indexing from 1
    // if I storealso non B. edges:
    // pointList[i].boundary()=
    // (edgeList[j].point(1)).boundary() &&(edgeList[j].point(2)).boundary();
    pe->setPoint(3,pp); //use overloaded version that takes a pointer
    _edge=makeBareEdge(i1,i2);
    _edgeid=_be.addIfNotThere(_edge.first,pp->id());
  }	
  // Now the other edges, of which I do NOT build the global stuff
  // (I would need to check the switch but I will do that part later on)
  if (GeoShape::nDim ==3 ){
    UInt nbf=mesh.numBFaces();
    UInt nbv=GeoBShape::numVertices;
    out<<"Processing "<<mesh.storedFaces()<<" Face Edges"<<endl;
    for (UInt k=1; k<=mesh.storedFaces(); ++k){
      pbe= &mesh.face(k);
      for (UInt j=1;j<=mesh.numLocalEdgesOfFace();j++){
	i1=GeoBShape::eToP(j,1);
	i2=GeoBShape::eToP(j,2);
	i1=(pbe->point(i1)).id();
	i2=(pbe->point(i2)).id();
	_edge=makeBareEdge(i1,i2);
	  e_id=_be.id(_edge.first);
	  if(e_id!=0){
	    pp=&mesh.point(e_id);
	  } else {
	    // new edge -> new Point
	    pp=&mesh.addPoint(k<=nbf); // true for boundary points
	    _edgeid=_be.addIfNotThere(_edge.first,pp->id());
	    pp->x()=(mesh.point(i1).x()+
		     mesh.point(i2).x())*.5;
	    pp->y()=(mesh.point(i1).y()+
		     mesh.point(i2).y())*.5;
	    pp->z()=(mesh.point(i1).z()+
		     mesh.point(i2).z())*.5;
	    // If we have set a marker for the face, that marker is inherited by the new created point
	    pp->setMarker(pbe->marker());
	  }
	  pbe->setPoint(nbv+j, pp );
	}
      }
    }

  out<<"Processing "<< mesh.numElements()<<" Mesh Elements"<<endl;
  UInt nev=GeoShape::numVertices;
  for (UInt k=1; k<= mesh.numElements(); ++k){
    pv= &mesh.element(k);
    for (UInt j=1;j<=mesh.numLocalEdges();j++){
      i1=ele.eToP(j,1);
      i2=ele.eToP(j,2);
      i1=(pv->point(i1)).id();
      i2=(pv->point(i2)).id();
      _edge=makeBareEdge(i1,i2);
      e_id=_be.id(_edge.first);
      if(e_id!=0){
	pp=&mesh.point(e_id);
      } else {
	pp=&mesh.addPoint(false);// cannot be on boundary is the mesh is proper!
	_edgeid=_be.addIfNotThere(_edge.first,pp->id());
	pp->x()=(mesh.point(i1).x()+
		 mesh.point(i2).x())*.5;
	pp->y()=(mesh.point(i1).y()+
		 mesh.point(i2).y())*.5;
	pp->z()=(mesh.point(i1).z()+
		 mesh.point(i2).z())*.5;
	pp->setMarker(pe->marker());
      }
      pv->setPoint(nev+j, pp );
    }
  }
  /*=============================*/
  out << " ******* Done Construction of P2 Mmesh *******" <<endl<<endl;
}

/*
*****************************************************************************
                UTILITIES FOR GLOBAL CHECKINGS
*****************************************************************************
*/

//! This function performs  a lot of checks.

/*!The name is inappropriate since the tests that are performed are not just on the topological structure
  of the mesh. The output is directed to three output streams:
  <ul>
  <li> out-> usually standard output: important informative messages
  <li> err-> usually standard error: error messages
  <li> clog-> usually a file stream: informative messages which may be rather verbose.
  </ul>
 Furthermore, ths Switch sw (see switch.h) will return a set of
 keywords useful for other possible actions. If fix=true, this routines performes the
 steps needed to get an acceptable mesh, otherwise the input mesh is not modified . */

template <typename RegionMesh3D>
bool checkMesh3D(RegionMesh3D & mesh, Switch & sw,
		       bool fix=true, bool verbose=false,
		       ostream & out=cerr, ostream & err=cerr, ostream & clog=cout){

  if(mesh.storedPoints()==0){
    err<<"FATAL: mesh does not store points: I cannot do anything"<<endl;
    sw.create("ABORT_CONDITION",true);
    sw.create("NOT_HAS_POINTS",true);
    return false;
  }
  
  if(!checkIdnumber(mesh.pointList)){
    err<<"ERROR: points ids where wrongly set"<<endl;
    err<<"FIXED"<<endl;
    if (fix) sw.create("FIXED_POINTS_ID",true);
    if (fix) fixIdnumber(mesh.pointList);
  }
  
  if(!checkMarkerSet(mesh.pointList)){
    err<<"WARNING: Not all points have marker flag set"<<endl;
    sw.create("POINTS_MARKER_UNSET",true);
  }

  
  //-------------------------------------------------------------------------------------------
  //                                    VOLUMES
  //-------------------------------------------------------------------------------------------
  
  if(mesh.storedVolumes()==0){
    err<<"FATAL: mesh does not store volumes: I cannot do anything"<<endl;
    sw.create("ABORT_CONDITION",true);
    sw.create("NOT_HAS_VOLUMES",true);
    return false;
  }

  if(!checkIdnumber(mesh.volumeList)){
    err<<"ERROR: volume ids where wrongly set"<<endl;
    err<<"FIXED"<<endl;
    if (fix)  sw.create("FIXED_VOLUMES_ID",true);
    if (fix) fixIdnumber(mesh.volumeList);
  }
  
  if(!checkMarkerSet(mesh.volumeList)){
    err<<"WARNING: Not all volumes have marker flag set"<<endl;
    sw.create("VOLUMES_MARKER_UNSET",true);
    if(fix){
      for (typename RegionMesh3D::Volumes::iterator iv=mesh.volumeList.begin();
	   iv !=mesh.volumeList.end();++iv){
	if (iv->markerUnset())iv->setMarker(mesh.marker());
	    }
    }
  }
  
  if(mesh.numElements()<mesh.storedVolumes()){
    err<<"WARNING: Mesh Volumes must be at least "<< mesh.storedVolumes()<<endl;
    if (fix) mesh.numVolumes()=mesh.storedVolumes();
    if (fix) sw.create("FIXED_VOLUME_COUNTER",true);
  }
  
  // test now orientation
  
  SimpleVect<bool> * elSign=new SimpleVect<bool>;
  
  Real meshMeasure=checkVolumes(mesh,*elSign,sw);
  UInt positive;

  if(sw.test("SKIP_ORIENTATION_TEST")){
    clog<<"W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" <<endl;
    err <<"W: ELEMENT ORIENTATION NOT IMPLEMENTED YET FOR THIS TYPE OF ELEMENTS, SKIP" <<endl;
  }else if(sw.test("HAS_NEGATIVE_VOLUMES")){
    positive=count(elSign->begin(),elSign->end(),true);
    clog<< positive<<"W: positive elements out of" <<mesh.storedVolumes()<<endl;
    if (fix) clog<<"Fixing negative elements"<<endl;
    if (fix) fixVolumes(mesh,*elSign,sw);
    if(sw.test("ABORT_CONDITION")){
      err<<"ABORT: Cannot fix volumes, this element is not supported"<<endl;
      return false;
    } else {
      sw.unset("HAS_NEGATIVE_VOLUMES");
      meshMeasure=checkVolumes(mesh,*elSign,sw);
      if(sw.test("HAS_NEGATIVE_VOLUMES")){
	if(fix) err<<"ABORT: Cannot fix volumes: something wrong with this mesh"<<endl;
	if(fix) sw.create("ABORT_CONDITION",true);
	return false;
      }
    }
  }
  delete elSign; // free some memory 
  
  clog<< "Volume enclosed by the mesh= "<<meshMeasure<<endl
     <<"(Computed by integrating mesh elements measures)"<<endl
     <<"(Using 1 point Quadrature rule)"<<endl;
  
  //-----------------------------------------------------
  //                                    BOUNDARY FACES
  //-----------------------------------------------------
  
  TempFaceContainer * bfaces=new TempFaceContainer;
  UInt numInternalFaces, numFaces;
  
  UInt bFacesFound=findBoundaryFaces(mesh, *bfaces,numInternalFaces);

  numFaces=bFacesFound+numInternalFaces;
  
  EnquireBFace<RegionMesh3D> enquireBFace(mesh,*bfaces);
  
  
  if(mesh.storedFaces()==0 || mesh.numBElements()> mesh.storedFaces() || bFacesFound>mesh.storedFaces()){
    err<< "ERROR: Not all boundary faces stored"<<endl;
    if (fix) sw.create("BUILD_BFACES",true);
    if (fix) buildFaces(mesh, clog, err, bFacesFound, numInternalFaces, true, false,  verbose, bfaces);
  }
  else{
    //
    // Make sure BFaces are first
    // Here I need to use a method that does not require the proper
    // setting of boundary Points!
      if(mesh.numBFaces()<mesh.storedFaces()){
	if (fix) std::stable_partition(mesh.faceList.begin(),mesh.faceList.end(),enquireBFace);
	if (fix) fixIdnumber(mesh.faceList);
	if (fix) sw.create("FIXED_BFACES_FIRST");
      }
      
      
      if(!checkIdnumber(mesh.faceList)){
	err<<"ERROR: face ids where wrongly set"<<endl;
	err<<"FIXED"<<endl;
	if (fix) sw.create("FIXED_FACES_ID",true);
	if (fix) fixIdnumber(mesh.faceList);
      }
      
      // Check Consistency with the mesh. Beware that this method changes *bfaces!
      
      if (fix) fixBoundaryFaces(mesh, clog, err, sw, numFaces,bFacesFound,true,verbose,bfaces);
      
      if(mesh.storedFaces()==0){
	err<<"ABORT CONDITION: cannot find boundary faces"<<endl;
	sw.create("NOT_HAS_FACES",true);
	sw.create("ABORT_CONDITION",true);
      }
      
      
      if(mesh.numBFaces()==0){
	err<<" MeshBFaces counter is unset"<<endl;
	if (fix) mesh.setNumBFaces(mesh.storedFaces());
	if (fix) sw.create("FIXED_BFACE_COUNTER",true);
	if (fix) mesh.setLinkSwitch("HAS_BOUNDARY_FACES");
      }
      
      
      if(!checkMarkerSet(mesh.faceList)){
	err<<"WARNING: Not all faces have marker flag set"<<endl;
	sw.create("FACE_MARKER_UNSET",true);
	if (fix) setBFacesMarker(mesh,clog,err,verbose);
	if (fix && checkMarkerSet(mesh.faceList)){
	  sw.create("FACE_MARKER_UNSET",false);
	  sw.create("FACE_MARKER_FIXED",true);
	}
      }
    }
  
  
  if(mesh.numFaces()!=bFacesFound+numInternalFaces){
    err<<"WARNING Number of faces incorrectly set"<<endl;
    err<<"It was "<<mesh.numFaces()<<endl;
    err<<"It should be" <<bFacesFound+numInternalFaces<<endl;
    if (fix) err<<"Fixing"<<endl;
    if (fix) mesh.numFaces()=bFacesFound+numInternalFaces;
    if (fix) sw.create("FIXED_FACE_COUNTER",true);
  }

  if(fix&& mesh.storedFaces()>mesh.numBFaces())mesh.setLinkSwitch("HAS_ALL_FACES");

  out<<" Boundary faces found:"<< bFacesFound<<endl;
  out<<" Num Faces Stored stored:"<< mesh.storedFaces()<<endl;
  out<<" Boundary faces counter gives:"<< mesh.numBFaces()<<endl;
  delete bfaces;

  
  //-----------------------------------------------------
  //                                    BOUNDARY EDGES
  //-----------------------------------------------------
  
  TempEdgeContainer * bedges = new TempEdgeContainer;
  
  UInt bEdgesFound=findBoundaryEdges(mesh,*bedges);
  EnquireBEdge<RegionMesh3D> enquireBEdge(mesh,*bedges);

  UInt intedge;
  UInt Ned=0;
  TempEdgeContainer iedges;
  
  if(mesh.storedEdges()==0|| mesh.numBEdges()> mesh.storedEdges() || bEdgesFound>mesh.storedEdges() ){
    err<<"WARNING: mesh does not store all boundary edges"<<endl;
    sw.create("NOT_HAS_EDGES",true);
    if (fix) buildEdges(mesh,clog, err, bEdgesFound,intedge,true,false,verbose,bedges);
    Ned=bEdgesFound+intedge;
    if (fix) sw.create("BUILD_BEDGES",true);}
  else {
    
    // Make sure BEdges are first
    // Here I need to use a method that does not require the proper
    // setting of boundary Points!
    
    if (fix) std::stable_partition(mesh.edgeList.begin(),mesh.edgeList.end(),enquireBEdge);
    if (fix) fixIdnumber(mesh.edgeList);
    if (fix) sw.create("FIXED_BEDGES_FIRST");
    
    if(!checkIdnumber(mesh.edgeList)){
      err<<"ERROR: edge ids where wrongly set"<<endl;
      err<<"FIXED"<<endl;
      sw.create("FIXED_EDGES_ID",true);
      fixIdnumber(mesh.edgeList);
    }
    
    
    if(!checkMarkerSet(mesh.edgeList)){
      err<<"WARNING: Not all edges have marker flag set"<<endl;
      sw.create("EDGE_MARKER_UNSET",true);
      if (fix) setBEdgesMarker(mesh,clog,err,verbose);
      if (fix && checkMarkerSet(mesh.edgeList)){
	sw.unset("EDGE_MARKER_UNSET");
	sw.create("EDGE_MARKER_FIXED",true);
      }
    }
    out<< "Computing internal edges";
    if(fix)Ned=bEdgesFound+findInternalEdges(mesh,*bedges,iedges);
  }
  iedges.clear();
  delete bedges;
  
  if (mesh.numBEdges()!=bEdgesFound){
    err<<"WARNING: number of found boundary edges:"<<bEdgesFound<<endl<<
      " does not match that declared in mesh, i.e. "<<mesh.numBEdges()<<endl;
    if(mesh.numBEdges()==0){
      if (fix) err<<"FIXING"<<endl;
      if (fix) sw.create("FIXED_BEDGES_COUNTER",true);
      if (fix) mesh.setNumBEdges(bEdgesFound);
    }
  }

  if (Ned !=mesh.numEdges()){
    if(fix) err<<"WARNING: Counter of number of edges badly set: Should be (actual number)"<<Ned<<endl;
    err<<"It is instead equal to "<<mesh.numEdges();
    if (fix){
      err<<" **FIXED"<<endl;
      mesh.numEdges()=Ned;
    }
    cerr<<endl;
  }
  UInt nbed;
  UInt counte=testClosedDomain_Top(mesh,nbed);
  if(counte==0){
    out<<"**DOMAIN SURFACE IS (TOPOLOGICALLY) CLOSED"<<endl;
  } else {
    sw.create("DOMAIN_NOT_CLOSED",true);
    err<<"WARNING: DOMAIN APPEARS TO HAVE AN OPEN BOUNDARY (TOPOLOGY CHECK)"<<endl;
    err<<"Number of inconsistent edges:"<<counte<<endl;
  }

  

  //-----------------------------------------------------
  //                                    POINTS
  //-----------------------------------------------------
  EnquireBPoint<RegionMesh3D> enquirebpoint(mesh);

  UInt foundBPoints=std::count_if(mesh.pointList.begin(),mesh.pointList.end(),enquirebpoint);
  
  if(foundBPoints==0 || foundBPoints<mesh.storedBPoints()){
    err<<"WARNING Bpoints indicator not correctly set"<<endl;
    if (fix) err<<"FIXING by recomputing from boundary faces"<<endl;    
    fixBPoints(mesh, clog, err,verbose);
    if (fix) foundBPoints=std::count_if(mesh.pointList.begin(),mesh.pointList.end(),enquirebpoint);
    if (fix) sw.create("FIXED_BOUNDARY_POINTS",true);
  }

  if (! checkMarkerSet(mesh.pointList)){
    err<<"WARNING B. Points MARKER incorrectly set"<<endl;
    if (fix){
      setBPointsMarker(mesh, clog, cerr, verbose);
      if (! checkMarkerSet(mesh.pointList)){
	err<<"Cannot Fix Points MARKER"<<endl;
	sw.create("POINT_MARKER_UNSET",true);
      } else {
	err<<"FIXED"<<endl;
	sw.create("FIXED_POINT_MARKER",true);
      }
    }
  }
  
  
  if(mesh.storedBPoints()==0){
    err<<"WARNING B. Points COUNTER incorrectly set"<<endl;
    if (fix) setBPointsCounters(mesh) ;
    if (fix) err<<" FIXED"<<endl;
    if (fix) sw.create("FIXED_BPOINTS_COUNTER",true);
  }
  
  if(mesh.numPoints()==0){
    err<<"WARNING Points Counter unset"<<endl;
    if (fix) mesh.numPoints()=mesh.storedPoints();
    if (fix) sw.create("FIXED_POINTS_COUNTER",true);
  }
  
  //-----------------------------------------------------
  //                                   FINAL CHECKS
  //-----------------------------------------------------
  out<<" ********     COUNTERS CONTENT **********************************"<<endl;
  
  out<<" Num Volumes:  "<<mesh.numVolumes()<<endl;
  out<<" Num Vertices: "<<mesh.numVertices()<<" Num B. Vertices: "<<mesh.numBVertices()<<endl;
  out<<" Num Points:   "<<mesh.numPoints()<<   " Num B. Points  : "<<mesh.numBPoints()<<endl;
  out<<" Num Edges:    "<<mesh.numEdges()<<   " Num B. Edges   : "<<mesh.numBEdges()<<endl;
  out<<" Num Faces:    "<<mesh.numFaces()<<   " Num B. Faces   : "<<mesh.numBFaces()<<endl;
  out<<" ********     END COUNTERS **********************************"<<endl;
  
  bool eulok(true);
  
  eulok=(2*mesh.numFaces()-mesh.numLocalFaces()*mesh.numVolumes()-mesh.numBFaces()==0);
  if (RegionMesh3D::ElementShape::Shape == TETRA){
    out<< endl<<"*** CHECKING WITH EULER FORMULAE *****  "<<endl;
    eulok=eulok &(mesh.numEdges()-mesh.numVolumes()-mesh.numVertices()-(3*mesh.numBFaces()-2*mesh.numBVertices())/4)==0;
  } 
  if(! eulok){
    err<<" WE DO NOT SATISFY EUREF FORMULAS "<<endl;
    err<<"2*nFa=nFxV*nVo+nBFa";
    if (RegionMesh3D::ElementShape::Shape == TETRA) err<< "or nEd=nVo+nVe+(3*nBFa-2*nBVe)/4 ";
    err<<endl;
    sw.create("NOT_EULER_OK");
  } else {
    out<< endl<<"*** CHECKING IS OK *****  "<<endl;
  }
  
  mesh.setLinkSwitch("HAS_BEEN_CHECKED");
  
  return true;

  
}
#endif  
