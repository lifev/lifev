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
  \file dofInterfaceHandler.h
  \brief Class for connecting the dof of a mesh (3D) and an interface (2D)
         that lives on the boundary of the mesh.
  \version 1.0
  \author V. Martin
  \date 02/2003
*/
#ifndef _DOFINTERFACEHANDLER_HH
#define _DOFINTERFACEHANDLER_HH

#include "lifeV.hpp"
#include "dofInterface3Dto2D.hpp"
#include "bcVector.hpp"

#include <vector>  //STL vector class

namespace LifeV
{
class DofInterfaceHandler
{

public:

    typedef boost::shared_ptr<DofInterface3Dto2D> dof_interface_type;

    /*! Default Constructor (call addNeigbor() after for each interface)
      \param NbNeigh  number of neighbouring subdomains
    */
    DofInterfaceHandler( const UInt & NbNeigh );

    //! Constructor
    /*! \param NbNeigh  number of neighbouring subdomains
        \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
        \param dof1 the Dof object of the mesh in which we want to make the computations
    */
    DofInterfaceHandler( const UInt & NbNeigh, const LocalDofPattern & refFE, const Dof & dof1 );

    //!  Add a DofInterface3Dto2D to the list of neighbors
    /*!
      \param refFe the part of the reference FE that contains the dof patterns (nbDofPerEdge...)
      \param dof1 the Dof object of the mesh in which we want to make the computations
     */
    void addNeighbor( const LocalDofPattern& refFE, const Dof& dof1 );

    //! creates the list of Interface BC Vectors
    //! (that store the interface values) :
    //! InIBC , OutIBC
    void initVectors();

    //! creates the list of the classes of BCVectorInterface
    void initBCVectorInterface();

    //! How many neighbors stored?
    UInt NbNeigh() const;

    //! Sum of the sizes of the Interface vectors
    //!  (it is the same for InIBC, OutIBC)
    //! It corresponds to the total size of the interface unknowns
    //! for the current subdomain.
    UInt NbInterfaceUnknowns() const;

    //! extracting a neighbor in the _neighList list (starts from 0)
    const DofInterface3Dto2D & operator[] ( const UInt & i ) const;
    DofInterface3Dto2D & operator[] ( const UInt & i ) ;

    //! extracting a Vector in the Input _InIBCList list (starts from 0)
    const Vector & InIBC( const UInt & i ) const;
    Vector & InIBC( const UInt & i ) ;

    //! extracting a Vector in the Input _InIBCList list (starts from 0)
    //! using the reference of the interface
    const Vector & InIBC_byRefInterf( const Int & refinterf ) const;
    Vector & InIBC_byRefInterf( const Int & refinterf ) ;

    //! extracting a Vector in the Output _OutIBCList list (starts from 0)
    const Vector & OutIBC( const UInt & i ) const;
    Vector & OutIBC( const UInt & i ) ;

    //! extracting a Vector in the Output _OutIBCList list (starts from 0)
    //! using the reference of the interface
    const Vector & OutIBC_byRefInterf( const Int & refinterf ) const;
    Vector & OutIBC_byRefInterf( const Int & refinterf ) ;

    //! extracting a BCVector in the _bcvList list (starts from 0)
    const BCVectorInterface & BCvec( const UInt & i ) const;
    BCVectorInterface & BCvec( const UInt & i ) ;

    /*! This method returns the corresponding index number in the vectors
      (_neighList, _InIBCList ...) living on the interfaces
      for a specific reference interface number.
      \param interfref : the reference of the interface.
    */
    UInt IndexOfInterfaceRef( const Int& interfref ) const;

    //! save some memory: destroy some unused lists in _neighList
    //! use it safely!
    //!(USE IT AFTER the construction of the interface mesh
    //! and the _locDofMap in DofInterfaceBase)
    void ClearSomeDofInterface3Dto2DList();

private:
    const UInt _nbNeigh; //!< number of neighbors

    //! list of DofInterface3Dto2D (each one lives on one interface)
    std::vector< dof_interface_type > _neighList;

    //! list of INPUT Vectors that live on the interfaces (one per interface)
    //! INPUT interface Boundary Conditions values.
    //! (This is the vector of unknows that lives on the interfaces of the
    //! current subdomain in a domain decomposition method).
    //! It is also the vector that contains the interface boundary conditions data.
    std::vector< Vector > _InIBCList;

    //! list of OUTPUT Vectors that live on the interfaces (one per interface)
    //! OUTPUT interface Boundary Conditions values.
    //! (This is the vector of unknows that live on the interfaces of the
    //! current subdomain in a domain decomposition method).
    std::vector< Vector > _OutIBCList;

    //! list of classes of BCVectorInterface (one per interface)
    std::vector< BCVectorInterface > _bcvList;

    //! map between the index number (second) in the vectors
    //!  (_neighList, _InIBCList ...)
    //!  and its reference interface number (first).
    std::map<Int, UInt> _indexInterfRefMap;

};
}
#endif
