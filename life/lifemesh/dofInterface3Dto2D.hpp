/*!
  \file dofInterface3Dto2D.h
  \brief Class for connecting the dof of a mesh (3D) and an interface (2D)
         that lives on the boundary of the mesh.  
  \version 1.0
  \author V. Martin (copy-modif of M.A. Fernandez )
  \date 02/2003

  This file contains the class which may be used to update and hold 
  the connections between the dof of a mesh (3D) and an interface (2D)
  that lives on the boundary of the mesh. The interface is referenced
  by a flag.

*/
#ifndef _DOFINTERFACE3DTO2D_HH
#define _DOFINTERFACE3DTO2D_HH

#include "dofInterfaceBase.hpp"

#include "localDofPattern.hpp"
#include "dof.hpp"
#include <iostream>
#include <map>
#include <list>   //necessary to write vertices in order.
#include <vector>   //necessary to write faces and access them arbitrarily.

#include "markers.hpp"
#include <ext/slist> 
using namespace __gnu_cxx;

/*! 
  \class DofInterface3Dto2D

  Base class which holds the connections of the dof between a 3D mesh and its 2D interface
  
  The connections may be built by calling the update method.
  
*/
class DofInterface3Dto2D: 
  public DofInterfaceBase 
{
 public:

  //! Constructor for interfacing Dof of the same type (LocalDofPattern)
   /*!
    \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
    \param dof1 the Dof object of the mesh in which we want to make the computations
   */
  DofInterface3Dto2D( const LocalDofPattern& refFE, const Dof& dof1 );
  
  //! This method builds the Dof connections at the interface
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param flag1 the marker of the interface in the mesh1
   */
  template<typename Mesh> 
    void update( const Mesh& mesh1, const EntityFlag& flag1 );

  //! Creates an Inria medit type mesh for the interface (pseudo 3D)
  template<typename Mesh> 
    void Generate2DMesh( string fname, const Mesh& mesh1 ) const;
  
  //! Returns the reference of the interface
  EntityFlag InterfaceRef() const;

  //! Returns the identity of the i-th elements in the (finalised) face list 
  //! (counting from 0 ' a la C')
  ID operator[](const UInt& i) const;

  //! Assignment operator (we have a vector of DofInterface3Dto2D)
  DofInterface3Dto2D & operator=(const DofInterface3Dto2D& dofi);

  //! true if the lists have been updated.
  bool finalized() const;

  //! removes all unuseful list (all except _faceList). use it properly!
  void ClearLists();

  //! output
  std::ostream& showMe2D(bool verbose=false, std::ostream& out=std::cout ) const ;

 private:
  //! reference of the interface
  EntityFlag _interfRef;
  
  //! LocalDofPattern object used in the mesh in which we want to make the computations
  const LocalDofPattern * _refFE1;

  //! Dof object of the mesh in which we want to make the computations
  const Dof *  _dof1;

  /*! STL list which holds the connections between faces at the interface
        -> first  : global (3D) face number
        -> second : interface (2D) face number

  The name has changed : it was _elc before. (V.M.) 
  */
  vector< pair<ID,ID> > _faceList;  

  /*! Auxiliary STL list which holds the connections between vertices at the interface
  (vertices appear more than once. (Mathematically, it is a family))  
  Empty after calling update
  */
  list<ID> _vertexPerFaceList;

  /*! STL list which holds the connections between vertices at the interface 
        -> first  : global (3D) vertex number
        -> second : interface (2D) vertex number
  */
  list< pair<ID,ID> > _vertexList;  

  /*! Auxiliary STL list which holds the connections between edges at the interface
  (edges appear more than once. (Mathematically, it is a family))  
  Empty after calling update
  */
  // list<ID> _edgePerFaceList;

  /*! STL list which holds the connections between edges at the interface 
        -> first  : global (3D) edge number
        -> second : interface (2D) edge number
  */
  // list< pair<ID,ID> > _edgeList;  

  //!  Auxiliary STL list which holds the connections between Dof at the interface   
  //! Empty after calling update
  slist< pair<ID,ID> > _locDof;

  //!  STL iterator type for the lists
  typedef slist< pair<ID,ID> >::iterator Iterator;
  
  //! Transforms the 3d index of a vertex into its 2d (interface) index.
  //! This is a simple algorithm... Find out something better some day...?
  ID _Vtx3Dto2D( const ID& idpoint3D ) const;

  //! This method builds the connections between faces at the interface (_faceList container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param flag1 the marker of the interface in the mesh1
   */
  template<typename Mesh>
    void _updateFaceConnections( const Mesh& mesh1, const EntityFlag& flag1 );

  //! This method builds the list of vertices at the interface (_vertexList container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
  */
  template<typename Mesh> 
    void _updateVertices(const Mesh& mesh1);

  //! This method builds the list of edges at the interface (_edgeList container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
  */
  // template<typename Mesh> 
  // void _updateEdges(const Mesh& mesh1);
  
  //! This method builds the connections between Dof at the interface (_locDof container)
  /*!
    \param mesh1 the mesh in which we want to make the computations
    \param dof1 the Dof object of the mesh in which we want to make the computations
   */
  template<typename Mesh> 
    void _updateDofConnections( const Mesh& mesh1 );
  
  //! true if the lists have been updated.
  bool _finalized;

};


//! This method builds the connections between faces at the interface (_faceList container)
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param flag1 the marker of the interface in the mesh1
  \param tol tolerance for connecting points of both meshes at the interface 
*/
template<typename Mesh> 
void DofInterface3Dto2D::_updateFaceConnections(const Mesh& mesh1, const EntityFlag& flag1) {

  UInt bdnF1  = mesh1.numBFaces(); // Number of boundary faces in mesh1

  EntityFlag marker1;
  
  typedef  typename Mesh::FaceShape GeoBShape; // Shape of the faces

  ID fcounter = 1;  //! Face on the interface counter

  //! Loop on boundary faces on mesh1
  for (ID ibF1=1; ibF1<=bdnF1; ++ibF1) {

    //! The face marker
    marker1 = mesh1.boundaryFace(ibF1).marker();
    
    //! Is the face on the interface?
    if (marker1 == flag1) {

      pair<ID,ID> fp( ibF1, fcounter );
      _faceList.push_back( fp ); 
      fcounter ++;  //! local face number

    }
  }
  ASSERT_PRE( _faceList.size()== --fcounter,
	      "Local face counter and list size do not match (in face loop).");
}

//! This method builds the list of vertices at the interface (_vertexList container)
/*!
  \param mesh1 the mesh in which we want to make the computations
*/
template<typename Mesh> 
void DofInterface3Dto2D::_updateVertices(const Mesh& mesh1) {

  typedef  typename Mesh::FaceShape GeoBShape; // Shape of the faces

  UInt nVpF = GeoBShape::numVertices;

  // ID vcounter = 1;
  
  // Loop on faces at the interface (matching faces)
  for (vector< pair<ID,ID> >::iterator i=_faceList.begin(); i != _faceList.end(); ++i) {
    for (ID jvtx = 1 ; jvtx <= nVpF ; ++ jvtx ) {

      _vertexPerFaceList.push_back( mesh1.boundaryFace(i->first).point(jvtx).id() );      
      //vcounter ++;  //! number of vertices (multiple appearances)
    }
  }

  RemoveMultiple( _vertexPerFaceList , _vertexList );

  //save memory
  _vertexPerFaceList.clear();
}

//! This method builds the list of edges at the interface (_edgeList container)
/*!
  \param mesh1 the mesh in which we want to make the computations
*/
/*
template<typename Mesh> 
void DofInterface3Dto2D::_updateEdges(const Mesh& mesh1) {

  typedef  typename Mesh::VolumeShape GeoShape;
  typedef  typename Mesh::FaceShape GeoBShape; // Shape of the faces

  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges 

  ID ecounter = 1;

  ID iFaEl1, iEdEl1;
  
  // Loop on faces at the interface (matching faces)
  for (vector< pair<ID,ID> >::iterator i=_faceList.begin(); i != _faceList.end(); ++i) {

    iFaEl1 = mesh1.boundaryFace(i->first).pos_first(); // local id of the face in its adjacent element (mesh1)

    // loop on face edges (mesh1)
    for (ID iEdFa1=1; iEdFa1<=nFaceE; ++iEdFa1) {
      
      iEdEl1  = GeoShape::fToE(iFaEl1,iEdFa1).first; // local edge number (in element)
	
      _edgePerFaceList.push_back( mesh1.boundaryFace(i->first).edge(iEdFa1).id() );      
      ecounter ++;  //! number of vertices (multiple appearances)
    }
  }

  RemoveMultiple( _edgePerFaceList , _edgeList );

  //save memory
  _edgePerFaceList.clear();
}
*/

//! This method builds the connections between Dof at the interface (_locDof container)
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param dof1 the Dof object of the mesh in which we want to make the computations
*/
template<typename Mesh> 
void DofInterface3Dto2D::_updateDofConnections( const Mesh& mesh1 ) {

  typedef  typename Mesh::VolumeShape GeoShape;
  typedef  typename Mesh::FaceShape GeoBShape;
  
  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges 
  
  // UInt nDofpV1 = _refFE1.nbDofPerVertex; // number of Dof per vertices on mesh1
  // UInt nDofpE1 = _refFE1nbDofPerEdge;   // number of Dof per edges on mesh1
  // UInt nDofpF1 = _refFE1.nbDofPerFace;   // number of Dof per faces on mesh1
  UInt nDofpV1 = _refFE1->nbDofPerVertex; // number of Dof per vertices on mesh1
  UInt nDofpE1 = _refFE1->nbDofPerEdge;   // number of Dof per edges on mesh1
  UInt nDofpF1 = _refFE1->nbDofPerFace;   // number of Dof per faces on mesh1

  UInt nElemV = GeoShape::numVertices; // Number of element's vertices 
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges

  //  UInt nDofFV1 = nDofpV1 * nFaceV; // number of vertex's Dof on a face on mesh1
  //  UInt nDofFE1 = nDofpE1 * nFaceE; // number of edge's Dof on a face on mesh1

  UInt nDofElemV1 = nElemV*nDofpV1; // number of vertex's Dof on a Element on mesh1
  UInt nDofElemE1 = nElemE*nDofpE1; // number of edge's Dof on a Element on mesh1

  ID iElAd1, iVeEl1, iFaEl1, iEdEl1, gDof1;

  ID locDofCounter1 = 1;

  // Loop on faces at the interface (matching faces)
  for (vector< pair<ID,ID> >::iterator i=_faceList.begin(); i != _faceList.end(); ++i) {
    
    iElAd1 = mesh1.boundaryFace(i->first).ad_first();  // id of the element adjacent to the face (mesh1) 
    
    iFaEl1 = mesh1.boundaryFace(i->first).pos_first(); // local id of the face in its adjacent element (mesh1)
        
    // Vertex based Dof on mesh1
    if ( nDofpV1 ) { 
      
      // loop on face vertices (mesh1)
      for (ID iVeFa1=1; iVeFa1<=nFaceV; ++iVeFa1){
	
	iVeEl1 = GeoShape::fToP(iFaEl1,iVeFa1); // local vertex number (in element)
	
	// Loop number of Dof per vertex (mesh1)
	for (ID l=1; l<=nDofpV1; ++l) {

          gDof1 = _dof1->localToGlobal(iElAd1,(iVeEl1-1)*nDofpV1 + l); // Global Dof on mesh1
          
          pair<ID,ID> locDof( gDof1, locDofCounter1 );   //! May be : invert the 2 ??
          _locDof.push_front( locDof ); // Updating the list of dof connections

          locDofCounter1 ++;  //! local Dof (total dof on the interface)
        }        
      }    
    }       
    
    // Edge based Dof on mesh1
    if (nDofpE1) { 
	
      // loop on face edges (mesh1)
      for (ID iEdFa1=1; iEdFa1<=nFaceE; ++iEdFa1) {
	
	iEdEl1  = GeoShape::fToE(iFaEl1,iEdFa1).first; // local edge number (in element)
	
	// Loop number of Dof per edge (mesh1)
	for (ID l=1; l<=nDofpE1; ++l) {
	  
          gDof1 = _dof1->localToGlobal(iElAd1, nDofElemV1 + (iEdEl1-1)*nDofpE1 + l); // Global Dof on mesh1
          
          pair<ID,ID> locDof( gDof1, locDofCounter1 ); 
          _locDof.push_front( locDof ); // Updating the list of dof connections

           locDofCounter1 ++;  //! local Dof (total dof on the interface)
	}
      }      
    }
        
    // Face based Dof on mesh1
    for (ID l=1; l<=nDofpF1; ++l) {  

      gDof1 = _dof1->localToGlobal(iElAd1, nDofElemE1 + nDofElemV1 + (iFaEl1-1)*nDofpF1 + l); // Global Dof in mesh1 
      
      pair<ID,ID> locDof( gDof1, locDofCounter1 );
      _locDof.push_front( locDof ); // Updating the list of dof connections

      locDofCounter1 ++;  //! local Dof (total dof on the interface)
    }
  }   

  
  
  // Updating the map containter with the connections
  for (slist< pair<ID,ID> >::iterator i=_locDof.begin(); i!=_locDof.end(); ++i) {
    _locDofMap[i->first] = i->second; 
  }
  
  // Saving memory
  _locDof.clear();
}


//! This method builds the Dof connections at the interface
/*!
  \param mesh1 the mesh in which we want to make the computations
  \param flag1 the marker of the interface in the mesh1
*/
template<typename Mesh>  
void DofInterface3Dto2D::update( const Mesh& mesh1, const EntityFlag& flag1 ) {

  // Updating face connections at the interface
  _updateFaceConnections( mesh1, flag1 );
  
  // Updating vertex connections at the interface
  _updateVertices( mesh1 );

  // _updateEdges( mesh1 );

  // Update of the Dof connections without inperpolation
  _updateDofConnections( mesh1 );

  _interfRef = flag1;

  _finalized = true; //! the lists are updated
}

// ============ Generate2DMesh ================ 
/*! Write the 2D Inria mesh of the interface (pseudo 3D)
    referenced by the number _interfRef. 
    
    It uses a inria "me.hpp" format. (for medit)
    
    You should have filled the lists of vertices and faces before.
    (Call fillmyinterface once previously).
*/
template<typename Mesh> 
void DofInterface3Dto2D::Generate2DMesh(string fname, const Mesh& mesh1) const {  
  
  ASSERT_PRE(_finalized, "The lists of vertices and faces must be finalized before generating the interface mesh." );

  ofstream ofile(fname.c_str());  
  ASSERT(ofile,"Error: Output file cannot be open"); 
  
  ID idpoint3D;
  ID idpoint2D;
  ID idface3D;

  typedef  typename Mesh::FaceShape FaceShape;
  UInt nVpF = FaceShape::numVertices;

  ofile << "MeshVersionFormatted 1\n";
  ofile << "Dimension 3\n";
  ofile << endl;

  ofile << "Vertices\n";
  ofile << _vertexList.size() << "\n";

  //! Write the coordinates of the vertices
  for (list< pair<ID,ID> >::const_iterator i2D = _vertexList.begin(); i2D!=_vertexList.end(); ++i2D) {
    idpoint3D = i2D->first;
    ofile << mesh1.pointList(idpoint3D).x() << " "
	  << mesh1.pointList(idpoint3D).y() << " "
	  << mesh1.pointList(idpoint3D).z() << " " 
          << mesh1.pointList(idpoint3D).marker() << endl; 
  }
  ofile << endl;
  
  switch(FaceShape::Shape) {
  case QUAD: 
    ofile << "Quadrilaterals\n";
    break;
  case TRIANGLE: 
    ofile << "Triangles\n";
    break;
  default:
    ERROR_MSG("This shape is not implemented in myinterface!");
  }

  ofile << _faceList.size() << "\n";

  //! Write the face table
  for (vector< pair<ID,ID> >::const_iterator i2D=_faceList.begin(); i2D != _faceList.end(); ++i2D) {
    idface3D = i2D->first;

    for (ID vtx = 1 ; vtx <= nVpF ; ++ vtx) {
      idpoint3D = mesh1.boundaryFace( idface3D ).point( vtx ).id();
      idpoint2D = _Vtx3Dto2D( idpoint3D ); //Simple algorithm (of Search in the list...)
      ofile << idpoint2D << " ";
    }
    ofile << mesh1.boundaryFace( idface3D ).marker() << endl;
  }
}
  


//! useful function to sort a list and remove multiple numbers.
void RemoveMultiple(const list<ID> & list0, list< pair<ID,ID> > & listf );

 
#endif
