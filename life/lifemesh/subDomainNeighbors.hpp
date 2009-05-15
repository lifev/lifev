/*
This file is part of the LifeV library
Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#ifndef _SUBDOMAINNEIGHBORS_H
#define _SUBDOMAINNEIGHBORS_H
#include <fstream>
#include <vector>
#include <life/lifecore/life.hpp>



#define PREFIX 666
#define ANUMBER 1000

namespace LifeV
{

/*!
  \file subDomainNeighbors.h
  \author V. Martin
  \date 04/02/2003

  \note Classes containing the informations necessary for a subdomain:
        the number of neighbors, their identity and the names of the
 common interfaces.
*/

/*! \class SDomNeighborData contains the data that a subdomain needs to see one
    neighbor:

    - the Identity of the neighbor,
    - the names (references) of the interface between the subdomain
      and its neighbor.
*/
/* We don't need a class : we use a struct instead...
class SDomNeighborData{
 public:
  //! constructor
  SDomNeighborData();

  //! output
  std::ostream & showMe( std::ostream  & out=std::cout ) const;

 public:
*/

//! useful function to sort a list and remove multiple numbers.
void RemoveMultiple( const std::vector<id_type> & list0, std::vector<id_type> & listf );

struct SDomNeighborData
{
    //! identity of the neighbor (it is a subdomain!)
    id_type NeighborID;

    //! reference of the interface between the subdomain and the neighbor
    Int InterfaceRef;
};


/*! \class SubDomainNeighbors contains the data that a subdomain needs to see all
    its neighbors:

    - the Identity of the subdomain,
    - the number of subdomains,
    - a list of SDomNeighborData.

    \note
    The same neighbor can appear more than once in the list, but with different
    interface references.
    Reason: an interface between one subdomain and one neighbor can be made
    of more than one planar surface, it may have more than one reference.

*/
class SubDomainNeighbors
{
public:

    //! Constructor
    SubDomainNeighbors( id_type SDomID );

    //! Constructor taking the connectivity table in a file as input
    SubDomainNeighbors( id_type SDomID, std::string fname );

    //! Destructor
    ~SubDomainNeighbors();

    //! Fill the lists _faceList and _vertexList
    template <typename Mesh>
    void fillSubDomainNeighbors( const Mesh& mesh );

    //! How many neighbors stored?
    size_type sizeNeigh() const;

    //! Is there no neighbors?
    bool emptyNeigh() const;

    //! Return the reference of the interface i in the neighbors' list. (Beware: i starts from 0).
    Int NeighInterfaceRef( const size_type & i ) const;

    //! output
    std::ostream & showMe( bool verbose = false, std::ostream & out = std::cout ) const;

private:
    //! Identity of the SubDomain
    id_type _SDomID;

    //! number of interface references
    size_type _nbInterf;

    //! number of neighboring subdomains
    size_type _nbNeigh;

    //! std::vector list holding the neihgboring subdomains and the reference of the common interface.
    std::vector<SDomNeighborData> _neighList;

    //! std::vector list holding the reference of the interfaces.
    std::vector<id_type> _interfList;

    //! True if the lists are complete.
    bool _finalized;
};


/********************************************************************
                      IMPLEMENTATIONS
********************************************************************/

// ============ fillSubDomainNeighbors ================

/*! Fill the lists _interfList and _neighList.

    Use the following rule for interface numbering:
    "666007024"  is the interface reference (PREFIX=666)
    between the subdomain "7" (=007) and its neighbor "24" (=024).
    This interface is seen from the subdomain "7" (it comes first).
    Modify PREFIX and ANUMBER (10, 100, 1000...) to change the rule.

    Of course, to use this function, you have to properly define the mesh.

    Another possibility is to provide the connectivity table directly
    through a file that should be read. (not implemented yet).

 Should be called once (and only once).
*/
template <typename Mesh>
void SubDomainNeighbors::fillSubDomainNeighbors( const Mesh& mesh )
{
    ASSERT_PRE( !_finalized , "The lists of neighbors are already finalized." );

    typedef typename Mesh::FaceShape FaceShape;

    UInt nBdF = mesh.numBFaces();
    UInt icounter = 0;  //! interface counter

    std::vector<id_type> _tmpList;

    //! loop on boundary faces
    for ( id_type k = 1 ; k <= nBdF ; ++ k )
    {
        //! Is the face on the interface?
        if ( mesh.boundaryFace( k ).marker() > PREFIX * ANUMBER * ANUMBER )
        {
            //! This number (666000000) relies on my numbering of the interfaces...
            //! Fill the list of Faces on the interface :
            _tmpList.push_back( mesh.boundaryFace( k ).marker() );
            icounter ++ ;
        }
    }

    //! To avoid troubles: you must provide a non empty interface.
    ASSERT_PRE( icounter > 0 , \
                "There is no Boundary face with a relevant interface reference. (> 666000000)" );

    /*
    // Sorting list of vertices.
    sort(_tmpList.begin(), _tmpList.end());

    _interfList.push_back( _tmpList[0] );
    icounter = 0;

    //! We remove the multiple occurences :
    for (size_type i = 1 ;  i < _tmpList.size() ; i++){
      if ( _tmpList[ i ] != _interfList[ icounter ] ){
        //! Fill the list containing the indexes of the interfaces (appears once).
        _interfList.push_back( _tmpList[i] );
        icounter ++ ;
      }
    }
    _nbInterf = ++icounter ;   //!< number of vertices on the interface

    */

    //! sort the list and remove multiple occurences
    RemoveMultiple( _tmpList, _interfList );

    _nbInterf = _interfList.size(); //!< number of vertices on the interface

    id_type neigh;
    SDomNeighborData sdneighdata;

    //! loop on faces on the interface.
    for ( size_type i = 0 ; i < _interfList.size() ; i++ )
    {
        neigh = _interfList[ i ] - PREFIX * ANUMBER * ANUMBER - _SDomID * ANUMBER ;
        ASSERT_PRE( neigh < ANUMBER && neigh > 0 , \
                    "The references of the interface are not relevant." );
        sdneighdata.NeighborID = neigh;
        sdneighdata.InterfaceRef = _interfList[ i ];
        _neighList.push_back( sdneighdata );
    }

    _nbNeigh = _neighList.size(); //!< number of neighbors on the interface

    _finalized = true;    //! You can't touch the lists after this point.
}
}

#endif



