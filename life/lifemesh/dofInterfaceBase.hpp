/*!
  \file dofInterfaceBase.h
  \brief Base Class for interfacing dofs between two meshes  
  \version 1.0
  \author M.A. Fernandez and V. Martin
  \date 11/2002

  This file contains the class which may be used to update and hold the connections between the dof
  on two matching meshes.  
*/
#ifndef _DOFINTERFACEBASE_HH
#define _DOFINTERFACEBASE_HH

#include "lifeV.hpp"
#include <map>
#include <fstream>
#include "vecUnknown.hpp"

using namespace __gnu_cxx;

/*! 
  \class DofInterfaceBase

  Base class which holds the conections of the dof in two matching meshes
  The dof mapping (STL map) is set in the derived classes. 
  Each derived class depends on the type of interface you have.
  
  Method getInterfaceDof gives the connections.
 
*/
class DofInterfaceBase {
 public:

  //! Default Constructor
  DofInterfaceBase();

  //! Virtual Destructor
  virtual ~DofInterfaceBase() {};
  
  /*! read the dof and values in a file and fill the map and the vector
   (use it at your own risks!)
   USAGE: file datavec.txt : nb_couples
                             couple:(idof, value) repeated nb_couples times
  */
  void ReadVectorDataAndDofMap(const string filename, Vector& dataVec);

  //! This method returns the corrresponding dof number of the mesh2 at the interface 
  //! for a specific dof number at the interface in mesh1
  /*!
    \param i a dof number in mesh1 
  */
  ID getInterfaceDof(const ID& i) const;

  //! This method returns the number of dof that live on the interface
  ID nbInterfaceDof() const;

  //! output
  std::ostream& showMe(bool verbose=false, std::ostream& out=std::cout ) const;

 protected: 

  //!  STL map container which holds the connections between Dof at the interface
  map<ID,ID> _locDofMap;

};

#endif

