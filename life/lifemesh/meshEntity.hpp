/*! file MeshEntity.h */
#ifndef _MESHENTITY_HH_
#define _MESHENTITY_HH_
#include "lifeV.hpp"

using namespace std; //Pretty useless

//! Base class of all Mesh Entities
/*! It contains the Entity ID and it defined the comparison operators */
class MeshEntity
{
public:
  MeshEntity():_id(0){};
  MeshEntity(ID i):_id(i){};
  ID id()const {return _id;}
  ID & id() {return _id;}
  bool operator==(MeshEntity & e)const {return _id==e.id();};
  bool operator<=(MeshEntity & e)const {return _id<=e.id();};
  bool operator>=(MeshEntity & e)const {return _id>=e.id();};
protected:
  ID _id;
};


//! Base class with boundary
/*! Contains info on boundary position */
class MeshEntityWithBoundary : public MeshEntity
{
public:
  MeshEntityWithBoundary():MeshEntity(),_boundary(false){};
  MeshEntityWithBoundary(ID i, bool boundary=false):MeshEntity(i),_boundary(boundary){};
  //! Tells if  item is on the boundary
  bool boundary() const {return _boundary;};
  //! Changes boundary indicator
  bool & boundary()  {return _boundary;};
protected:
  bool _boundary;
};


//! Mesh Entity with an orientation 
/*! \note Not Used so far! */
class OrientedMeshEntity: public MeshEntity
{                      
public:
  OrientedMeshEntity():MeshEntity(),_orient(true){};
  OrientedMeshEntity(ID i, bool o=true):MeshEntity(i), _orient(o){}
  //! Return entity orientation
  bool orientation()const { return _orient;}
  //! Assigns entity orientation
  bool & orientation(){ return _orient;}
protected:
  bool _orient;
}; 

#endif
