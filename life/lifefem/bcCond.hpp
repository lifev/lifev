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
/*!
  \file bcCond.h
  \brief classes to handle boundary conditions.
  \version 1.0
  \author M.A. Fernandez
  \date 06/2002
        11/2002 Introduction of data vectors and method bdUpdate (moved from Dof class).

  This file contains the classes which may be used to store boundary
  conditions. A boundary condition objet will have the following elements: a name identifying
  a specific BC, a flag identifying a specific part of the mesh boundary, a type
  (Essential, Natural, Mixte), a mode of implementation (Scalar, Full, Component, Normal, Tangential), a
  functor holding the data function, a bool vector describing the components involved in this boundary
  condition and finally a list of pointers to identifiers allowing the user to know to which DOF
  the boundary condition applies.
*/

#ifndef __BCCOND_HH__
#define __BCCOND_HH__

#include "lifeV.hpp"
#include "identifier.hpp"
#include "markers.hpp"
#include "bcVector.hpp"
#include "dof.hpp"
#include "currentFE.hpp"
#include "currentBdFE.hpp"
#include <set>
#include <map>

namespace LifeV
{
/*! Boundary condition basic types
    Essential, Natural or Mixte
*/
enum BCType{Essential, Natural, Mixte};

/*! Type for boundary conditions application modes

  Scalar: for scalar problems
  Full: for vector problems involving all components
  Component: for vector problems not involving all compontents
  Normal: for vector problems dealing with the normal component
  Tangential: for vector problems dealing with tangential components
*/
enum BCMode{Scalar, Full, Component, Normal, Tangential};


// ============ BCFunction_Base ================

/*!

 \class BCFunction_Base

 Base class (STL functor) that holds the function used for imposing BC.

  The data functions given by the user must have the following declaration
  Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
  We can use inheritance to hold specific boundary condition data. See, for instance,
  Mixed boundary conditions.

*/

class BCFunction_Base{
 public:
  //! Type for a generic user defined  function
  typedef Real (*Function)(const Real&, const Real&, const Real&, const Real&, const ID&);

  //! Default constructor
  /*!
    The user must supply a function by calling setFunction(..)
  */
  BCFunction_Base(){}; // JFG (26/10/2002)

  //! Constructing from a user defined function
  /*!
    \param g the user defined function
  */
  BCFunction_Base(Function g);

  //! Constructing from a user defined functor
  /*!
    \param bcf user defined functor
  */
  BCFunction_Base(const BCFunction_Base& bcf);

  //! Set the function
  /*!
    \param g the user defined function
  */
  void setFunction(Function g);

  //! Overloading function operator by calling _g
  /*!
    \param t time
    \param x coordinate
    \param y coordinate
    \param z coordinate
    \param i component of the vector function
    \return i-component of the user defined fonction evaluted in (t,x,y,z)
  */
  Real operator()(const Real& t, const Real& x, const Real& y,
		  const Real& z, const ID& i) const;

 protected:
  //! user defined function
  Function _g;
};


/*!

 \class BCFunction_Mixte

 Class (STL functor) that holds the user defined fonctions for a Mixte BC.

  The data funcitions given by the user must have the following declaration
 Real g(const Real& time, const Real& x, const Real& y, const Real& z, const ID& icomp)
*/

class BCFunction_Mixte:public BCFunction_Base {
 public:


  //! Default constructor
  /*!
    The user must supply a function by calling setFunction(..)
  */
  BCFunction_Mixte(){}; // VM (04/02/2003)


  //! Constructing from user defined functions
  /*!
    \param g user defined function
    \param coef user defined function
  */
  BCFunction_Mixte(Function g, Function coef);

  //! Constructing from a user defined functor
  /*!
    \param bcf user defined functor
  */
  BCFunction_Mixte(const BCFunction_Mixte& bcf);


  //! Set the functions in the mixte case (beware : plural!)
  /*!
    \param g : the user defined function
    \param coef : user defined function
  */
  void setFunctions_Mixte(Function g, Function coef);


  //! Method to call the auxiliary user defined function
  /*!
    \param t time
    \param x coordinate
    \param y coordinate
    \param z coordinate
    \param i component of the vector function
    \return i-component of the user defined fonction evaluted in (t,x,y,z)
  */
  Real coef(const Real& t, const Real& x, const Real& y,
	    const Real& z, const ID& i) const;
 private:
  //! user defined function
  Function _coef;
};



// ============ BC_Base ================



/*!
  \class BC_Base

  Base class which holds the boundary condition information

   For each boundary condtion the user must give a name, a mesh flag, a type, a mode, a data BCFuncion,
  and three (or two in 2D) bools describing the components involved in this boundary
  condition. Finally the list of pointers to identifiers will be updated in the Dof class
  (BC_Handler::bdUpdate method).
  The idea is to not use inheritance from dis class... if hope it will possible.
*/


//!
class BC_Base{
 public:

  friend class BC_Handler;

  //! iterator type in the identifiers list
  typedef std::set<Identifier_Base*,identifierComp>::iterator IDIterator0;
  typedef std::vector<Identifier_Base*>::iterator IDIterator;


  //! Constructor for BC
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condition applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
    \param std::vector<ID> storing the list of components involved in this boundary condition
   */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	  BCFunction_Base& bcf, const std::vector<ID>& comp);

  //! Constructor for BC without components for Scalar, Tangential or Normal  mode problems
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
  */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type,
	  const BCMode& mode, BCFunction_Base& bcf);

  //! Constructor for BC without list of components for Full mode problems
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
    \param nComp the number of componets involved in this boundary condition
  */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type,
	  const BCMode& mode, BCFunction_Base& bcf, const UInt& nComp);


  //! Constructor for BC with data vector
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condition applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcv data vector
    \param std::vector<ID> storing the list of components involved in this boundary condition
   */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	  BCVector_Base& bcv, const std::vector<ID>& comp);

  //! Constructor for BC with data vector, without components for Scalar, Tangential or Normal  mode problems
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcv data vector
  */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type,
	  const BCMode& mode, BCVector_Base& bcv);

  //! Constructor for BC with data vector, without list of components for Full mode problems
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcv data vector
    \param nComp the number of componets involved in this boundary condition
  */
  BC_Base(const std::string& name, const EntityFlag& flag, const BCType& type,
	  const BCMode& mode, BCVector_Base& bcv, const UInt& nComp);


  //! Copy constructor for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
  /*!
    \param BCb a boundary condition
  */
  BC_Base(const BC_Base& BCb);

  //! Assignment operator for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
  /*!
    \param BCb a boundary condition
  */
  BC_Base & operator=(const BC_Base&);

  //! Destructor (we have a vector of pointers to ID's and a pointer to user defined functions)
  ~BC_Base();

  //! Returns the BC name
  std::string name() const;

  //! Returns the BC associated flag
  EntityFlag flag() const;

  //! Returns the BC type
  BCType type() const;

  //! Returns the BC mode
  BCMode mode() const;

  //! Returns the number of components involved in this boundary condition
  UInt numberOfComponents() const;

  //! Returns the global i-th component involved in the boundary condition
  /*!
    \param i the specified "local" component (from 1 to numberOfComponents)
    \return true if the specified component component is involved in the BC
  */
  ID component(const ID i) const;

  //! Returns wether the list is finalised and the vector of ID's is then accessible
  bool finalised() const;

  //! Overloading function operator by calling the _bcf() user specified function
  /*!
    \param t time
    \param x coordinate
    \param y coordinate
    \param z coordinate
    \param i component of the vector function
    \return i-component of the user defined fonction evaluted in (t,x,y,z)
  */
  Real operator()(const Real& t, const Real& x, const Real& y,
		  const Real& z, const ID& i) const;


  //! Returns a pointer  to the user defined STL functor
  BCFunction_Base* pointerToFunctor() const;


  //! True is a data vector has been provided
  bool dataVector() const;


  //! Returns a pointer  to the i-th elements in the (finalised) list
  //! (counting from 1 ' a la FORTRAN')
  Identifier_Base*  operator()(const ID& i) const;


  //! Overloading function operator by calling the
  /*!
    \param iDpof global dof number
    \param iComp component number
  */
  Real operator()(const ID& iDof, const ID& iComp) const;


  //! Returns the value of the mixte coefficient (in BC Vector)
  Real MixteCoef() const;

  //! Returns a pointer to the i-th elements in the (finalised) list
  //! (counting from 0 ' a la C')
  Identifier_Base* operator[](const Index_t& i) const;

  //! Add a new indentifier in the list
  void addIdentifier(Identifier_Base*);

  //! Returns the liste size
  UInt list_size() const;

  //! Output
  std::ostream &  showMe(bool verbose=false, std::ostream & out=std::cout) const;

  //! overloaded operator allowing decreasing ordering operations
  friend bool operator<(const  BC_Base& a, const BC_Base& b) {
    return ( a.type() > b.type() );}

  //! overloaded operator allowing finding operations
  friend bool operator==(const  BC_Base& a, const EntityFlag flag) {
    return a.flag() == flag;}

 private:
  //! name identifying a specific BC
  std::string _name;

  //! flag identifying a specific part of the mesh boundary
  EntityFlag _flag;

  //! the boundary condition type
  BCType _type;

  //! the boundary condition mode of application
  BCMode _mode;

  //! Pointer to a user defined functor
  BCFunction_Base* _bcf;


  //! Pointer to a user given data vector
  BCVector_Base* _bcv;

  //! True is a data vector has been provided
  bool _dataVector;

  //! the list of involved in this BC
  std::vector<ID> _comp;

  //! set of pointers to identifiers allowing the user to get hold the DOF
  //! to which the BC applies
  std::set<Identifier_Base*,identifierComp> list0;

  //! container for id's when the list is finalised
  std::vector<Identifier_Base*> _idList;

  //! true, when idlist updated
  bool _finalised;

  //! Transfert between list and vector containers
  void finalise();
};





// ============ BC_Handler ================


/*!

  \class BC_Handler

  Container for BC_Base classes


BC_Handler is a container for the boundary condition classes
It just uses a stl vector to store BC_Base objets. The usage is simple:
the user creates the data functors from the de user defined functions
 \verbatim
  BCFunction_Base gv(g);
  BCFunction_Mixte gp(h,q);
 \endverbatim
Then he/she specifies the number of BC and uses the Handler for creating the actual
BC objects and storing them.

\verbatim
 BC_Handler bc_v(3);
 bc_v.addBC("inlet",10,Essential,Full,gv,3);
 bc_v.addBC("outflow",11,Mixte,Scalar);
 bc_v.addBC("wall",10,Essential,Full,gv,3);
\endverbatim

*/


class BC_Handler
{
 public:

  typedef std::vector<BC_Base>::iterator Iterator;

  //! Constructor doing nothing (the user must call setNumber(..))
  BC_Handler();
  //! Constructor taking the number of BC to be stored
  // JFG (25/10/2002)
  BC_Handler(const ID&);
  BC_Handler(const ID&, const bool& fullEssential);

  //! Set the number of BC to be stored
  // JFG (25/10/2002)
  void setNumber(const ID& nbc);


  //! How many BC stored?
  Index_t size() const;

  //! Is there no BC?
  bool empty() const;

  //! Iterators for the begining and end of the BC list
  Iterator begin() {return _bcList.begin();}
  Iterator end()  {return _bcList.end();}

  //! add new BC to the list (user defined function)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
    \param std::vector<ID> storing the list of components involved in this boundary condition
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCFunction_Base& bcf, const std::vector<ID>& comp);
  //! add new BC to the list  without specified components for Scalar, Tangential or Normal  mode problems  (user defined function)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCFunction_Base& bcf);


  //! add new BC to the list without list of components for Full mode problems  (user defined function)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcf the function holding the user defined function involved in this boundary condition
    \param nComp the number of componets involved in this boundary condition
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCFunction_Base& bcf, const UInt& nComp);


  //! add new BC to the list (data vector)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcv data vector
    \param std::vector<ID> storing the list of components involved in this boundary condition
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCVector_Base& bcv, const std::vector<ID>& comp);


  //! add new BC to the list  without specified components for Scalar, Tangential or Normal  mode problemst (data vector)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condtion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Component, Normal, Tangential
    \param bcv data vector
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCVector_Base& bcv);


  //! add new BC to the list without list of components for Full mode problems t (data vector)
  /*!
    \param name the name of the boundary condition
    \param flag the mesh flag identifying the part of the mesh where the boundary condiion applies
    \param type the boundary condition type: Natural, Essential, Mixte
    \param mode the boundary condition mode: Scalar, Full, Normal, Tangential
    \param bcv data vector
    \param nComp the number of componets involved in this boundary condition
  */
  void addBC(const std::string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
	     BCVector_Base& bcv, const UInt& nComp);


  //! Build the boundary stuff
  /*!
    \param mesh the mesh
    \param feBd the current finite element on the boundary
    \param dof give the local to global table
     This method udapte the BC classes checking the markers
     on the boundary. It builds the list of identifiers depending
     on the BC type
  */
  template <typename Mesh>
    void bdUpdate(Mesh& mesh, CurrentBdFE& feBd, const Dof& dof);

  //! returns true if the bdUpdate has been done before
  bool bdUpdateDone() const;

  //! returns  true if all the stored BC are of Essential type
  bool fullEssential() const;

  //! extracting a BC in the list
  BC_Base& operator[](const Index_t& );
  const BC_Base& operator[](const Index_t&) const;

  //! output
  std::ostream & showMe(bool verbose=false, std::ostream  & out=std::cout) const;

 protected:
  //! number of BC to be stored
  ID _nbc;

  //! true if the bdUpdate has been done
  bool _bdUpdateDone;

  //! true if all the stored BC are of Essential type
  bool _fullEssential;

  //! vector list holding the stored BC
  std::vector<BC_Base> _bcList;
};


/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

//! 11/2002: The update method must be a method of BCh: it updates
//! BC_Handler objects not Dof objects


//! Build the boundary staff
template <typename Mesh>
void BC_Handler::bdUpdate(Mesh& mesh, CurrentBdFE& feBd, const Dof& dof) {

  typedef  typename Mesh::VolumeShape GeoShape;

  // Some useful local variables, to save some typing
  UInt nDofpV = dof.fe.nbDofPerVertex; // number of Dof per vertices
  UInt nDofpE = dof.fe.nbDofPerEdge;   // number of Dof per edges
  UInt nDofpF = dof.fe.nbDofPerFace;   // number of Dof per faces

  UInt bdnF  = mesh.numBFaces();    // number of faces on boundary

  EntityFlag marker; //will store the marker of each geometric entity

  typedef typename GeoShape::GeoBShape  GeoBShape;

  UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
  UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

  UInt nElemV = GeoShape::numVertices; // Number of element's vertices
  UInt nElemE = GeoShape::numEdges;    // Number of element's edges

  UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
  UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

  UInt nDofF  = nDofFV+nDofFE+nDofpF; // number of total Dof on a face

  UInt nDofElemV = nElemV*nDofpV; // number of vertex's Dof on a Element
  UInt nDofElemE = nElemE*nDofpE; // number of edge's Dof on a Element

  SimpleVect<ID> bdltg(nDofF);
  typedef std::vector<BC_Base>::iterator Iterator;
  Iterator where;
  std::vector<Iterator> whereList;

  UInt iElAd, iVeEl, iFaEl, iEdEl;
  ID lDof, gDof;
  Real x,y,z;

  // ===================================================
  // Loop on boundary faces
  // ===================================================
  for (ID ibF=1 ; ibF<=bdnF; ++ibF) {

    iElAd = mesh.boundaryFace(ibF).ad_first();  // id of the element adjacent to the face
    iFaEl = mesh.boundaryFace(ibF).pos_first(); // local id of the face in its adjacent element
    feBd.updateMeas( mesh.boundaryFace(ibF) );  // updating finite element information

    // ===================================================
    // Vertex based Dof
    // ===================================================
    if ( nDofpV ) {

      // loop on face vertices
      for (ID iVeFa=1; iVeFa<=nFaceV; ++iVeFa){

	marker = mesh.boundaryFace(ibF).point(iVeFa).marker(); // vertex marker
	iVeEl = GeoShape::fToP(iFaEl,iVeFa); // local vertex number (in element)

	// Finding this marker on the BC list
	whereList.clear();
	where = _bcList.begin();
	while ( ( where = find(where, _bcList.end(), marker) ) != _bcList.end() ) {
	  whereList.push_back(where);
	  ++where;
	}

	// Loop number of Dof per vertex
	for (ID l=1; l<=nDofpV; ++l) {

	  lDof =   (iVeFa-1) * nDofpV + l ; // local Dof
	  gDof =  dof.localToGlobal( iElAd, (iVeEl-1)*nDofpV + l); // global Dof
	  bdltg( lDof ) =  gDof; // local to global on this face

	  // Adding identifier
	  for (UInt i=0 ; i< whereList.size(); ++i) {
	    where = whereList[i];
	    switch( where->type() ) {
	    case Essential:
	      // Why kind of data ?
	      if ( where->dataVector() ) { // With data vector
		where->addIdentifier( new Identifier_Base(gDof) ); // We only need the dof number
	      }
	      else { // With user defined functions
		feBd.coorMap(x, y, z, feBd.refFE.xi(lDof-1), feBd.refFE.eta(lDof-1));
		where->addIdentifier(new Identifier_Essential(gDof,x,y,z));
	      }
	      break;
            case Natural:
	      if ( where->dataVector()  ) { // With data vector
	      	where->addIdentifier( new Identifier_Natural(gDof) );
	      }
	      break;
            case Mixte:
	      // Why kind of data ?
              // vincent please check again for your Mixte-FE it doesn't work for Q1
//	      if ( where->dataVector()  ) { // With data vector
//	      	where->addIdentifier( new Identifier_Natural(gDof) );
//	      }
	      break;
	    default:
	      ERROR_MSG("This boundary condition type is not yet implemented");
	    }
	  }
	}
      }
    }


    // ===================================================
    // Edge based Dof
    // ===================================================
    if (nDofpE) {

      // loop on face edges
      for (ID iEdFa=1; iEdFa<=nFaceE; ++iEdFa) {

	iEdEl  = GeoShape::fToE(iFaEl,iEdFa).first; // local edge number (in element)
	marker = mesh.boundaryEdge( mesh.localEdgeId(iElAd, iEdEl) ).marker(); // edge marker

	// Finding this marker on the BC list
	whereList.clear();
	where = _bcList.begin();
	while ( ( where = find(where, _bcList.end(), marker) ) != _bcList.end() ) {
	  whereList.push_back(where);
	  ++where;
	}

	// Loop number of Dof per edge
	for (ID l=1; l<=nDofpE; ++l) {

        lDof =  nDofFV + (iEdFa-1) * nDofpE + l ; // local Dof
        gDof =  dof.localToGlobal( iElAd, nDofElemV + (iEdEl-1)*nDofpE + l); // global Dof
        bdltg( lDof ) =  gDof; // local to global on this face

        // Adding identifier
        for (UInt i=0 ; i< whereList.size(); ++i) {
            where = whereList[i];
            switch( where->type() ) {
                case Essential:
                    // Why kind of data ?
                    if ( where->dataVector() ) { // With data vector
                        where->addIdentifier( new Identifier_Base(gDof) );
                    }
                    else { // With user defined functions
                        feBd.coorMap(x, y, z, feBd.refFE.xi(lDof-1), feBd.refFE.eta(lDof-1));
                        where->addIdentifier( new Identifier_Essential(gDof,x,y,z) );
                    }
                    break;
                case Natural:
                    // Why kind of data ?
                    if ( where->dataVector() ) { // With data vector
                        where->addIdentifier( new Identifier_Natural(gDof) );
                    }
                    break;
                case Mixte:
                    // Why kind of data ?
                    if ( where->dataVector() ) { // With data vector
                        where->addIdentifier( new Identifier_Natural(gDof) );
                    }
                    break;
                default:
                    ERROR_MSG("This boundary condition type is not yet implemented");
            }
        }
	}
      }
    }

    // ===================================================
    // Face based Dof
    // ===================================================
    marker = mesh.boundaryFace(ibF).marker(); // edge marker

    // Finding this marker on the BC list
    whereList.clear();
    where = _bcList.begin();
    while ( ( where = find(where, _bcList.end(), marker) ) != _bcList.end() ) {
      whereList.push_back(where);
      ++where;
    }

    // Adding identifier
    for (UInt i=0 ; i< whereList.size(); ++i) {
      where = whereList[i];
      switch( where->type() ) {
      case Essential:
	// Loop on number of Dof per face
	for (ID l=1; l<=nDofpF; ++l) {
	  lDof = nDofFE + nDofFV + l; // local Dof
	  gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
	  // Why kind of data ?
	  if ( where->dataVector() ) { // With data vector
	    where->addIdentifier( new Identifier_Base(gDof) );
	  }
	  else { // With user defined functions
	    feBd.coorMap(x, y, z, feBd.refFE.xi(lDof-1), feBd.refFE.eta(lDof-1));
	    where->addIdentifier( new Identifier_Essential(gDof,x,y,z) );
	  }
	}
	break;
      case Natural:

	// Why kind of data ?
        // vincent please check again for your Mixte-FE it doesn't work for Q1
          if ( where->dataVector()  )
          { // With data vector
              for (ID l=1; l<=nDofpF; ++l) {
                  lDof = nDofFE + nDofFV + l; // local Dof
                  gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
                  where->addIdentifier( new Identifier_Natural(gDof) );
              }
          }
          else
          {
              // Loop on number of Dof per face
              for (ID l=1; l<=nDofpF; ++l) {
                  lDof = nDofFE + nDofFV + l; // local Dof
                  gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
                  bdltg( lDof ) =  gDof; // local to global on this face
              }
              where->addIdentifier( new Identifier_Natural(ibF,bdltg) );
          }
          break;
          case Mixte:
              // Why kind of data ?
        // vincent please check again for your Mixte-FE it doesn't work for Q1
//	if ( where->dataVector()  ) { // With data vector
//	  for (ID l=1; l<=nDofpF; ++l) {
//	    lDof = nDofFE + nDofFV + l; // local Dof
//	    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
//	    where->addIdentifier( new Identifier_Natural(gDof) );
//	  }
//	}
//	else {
	  // Loop on number of Dof per face
	  for (ID l=1; l<=nDofpF; ++l) {
	    lDof = nDofFE + nDofFV + l; // local Dof
	    gDof = dof.localToGlobal( iElAd, nDofElemE + nDofElemV + (iFaEl-1)*nDofpF + l); // global Dof
	    bdltg( lDof ) =  gDof; // local to global on this face
	  }
	  where->addIdentifier( new Identifier_Natural(ibF,bdltg) );
//	}
	break;
      default:
	ERROR_MSG("This boundary condition type is not yet implemented");
      }
    }
  }

  whereList.clear();
  // ============================================================================
  // There is no more identifiers to add to the boundary conditions
  // We finalise de set of identifiers by transfering it elements to a std::vector
  // ============================================================================
  for (Iterator it= _bcList.begin(); it!=_bcList.end(); ++it) {
    it->finalise();
  }

  _bdUpdateDone=1;
}
}


#endif
