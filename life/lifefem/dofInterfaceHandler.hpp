//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER


/*!
    @file
    @brief Class for connecting the dof of a mesh (3D) and an interface (2D)
    that lives on the boundary of the mesh.

    @author V. Martin
    @date 00-02-2003

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef _DOFINTERFACEHANDLER_HH
#define _DOFINTERFACEHANDLER_HH

#include <life/lifecore/life.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>
#include <life/lifefem/bcVector.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <vector>  //STL vector class

namespace LifeV
{

/*!
  CLASS TO BE DELETED
 */
class DofInterfaceHandler
{

public:

    //! @name Public Types
    //@{

    typedef boost::shared_ptr<DofInterface3Dto2D> dof_interface_type;

    //@}


    //! @name Constructor & Destructor
    //@{

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

    //@}


    //! @name Methods
    //@{

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

    //! Sum of the sizes of the Interface vectors
    //!  (it is the same for InIBC, OutIBC)
    //! It corresponds to the total size of the interface unknowns
    //! for the current subdomain.
    UInt NbInterfaceUnknowns() const;

    /*! This method returns the corresponding index number in the vectors
      (M_neighborList, M_inputBCList ...) living on the interfaces
      for a specific reference interface number.
      \param interfref : the reference of the interface.
    */
    UInt IndexOfInterfaceRef( const Int& interfref ) const;

    //! save some memory: destroy some unused lists in M_neighborList
    //! use it safely!
    //!(USE IT AFTER the construction of the interface mesh
    //! and the M_locDofMap in DofInterfaceBase)
    void ClearSomeDofInterface3Dto2DList();

    //@}

    //! @name Operators
    //@{

    //! extracting a neighbor in the M_neighborList list (starts from 0)
    const DofInterface3Dto2D & operator[] ( const UInt & i ) const;
    DofInterface3Dto2D & operator[] ( const UInt & i ) ;

    //@}


    //! @name Get Methods
    //@{

    //! How many neighbors stored?
    UInt NbNeigh() const;

    //! extracting a Vector in the Input M_inputBCList list (starts from 0)
    const EpetraVector & InIBC( const UInt & i ) const;
    EpetraVector & InIBC( const UInt & i ) ;

    //! extracting a Vector in the Input M_inputBCList list (starts from 0)
    //! using the reference of the interface
    const EpetraVector & InIBC_byRefInterf( const Int & refinterf ) const;
    EpetraVector & InIBC_byRefInterf( const Int & refinterf ) ;

    //! extracting a Vector in the Output M_outputBCList list (starts from 0)
    const EpetraVector & OutIBC( const UInt & i ) const;
    EpetraVector & OutIBC( const UInt & i ) ;

    //! extracting a Vector in the Output M_outputBCList list (starts from 0)
    //! using the reference of the interface
    const EpetraVector & OutIBC_byRefInterf( const Int & refinterf ) const;
    EpetraVector & OutIBC_byRefInterf( const Int & refinterf ) ;

    //! extracting a BCVector in the M_BCVectorList list (starts from 0)
    const BCVectorInterface & BCvec( const UInt & i ) const;
    BCVectorInterface & BCvec( const UInt & i ) ;


    //@}


private:

    // number of neighbors
    const UInt M_nbNeighbor;

    // list of DofInterface3Dto2D (each one lives on one interface)
    std::vector< dof_interface_type > M_neighborList;

    //! list of INPUT Vectors that live on the interfaces (one per interface)
    //! INPUT interface Boundary Conditions values.
    //! (This is the vector of unknows that lives on the interfaces of the
    //! current subdomain in a domain decomposition method).
    //! It is also the vector that contains the interface boundary conditions data.
    std::vector< EpetraVector > M_inputBCList;

    //! list of OUTPUT Vectors that live on the interfaces (one per interface)
    //! OUTPUT interface Boundary Conditions values.
    //! (This is the vector of unknows that live on the interfaces of the
    //! current subdomain in a domain decomposition method).
    std::vector< EpetraVector > M_outputBCList;

    //! list of classes of BCVectorInterface (one per interface)
    std::vector< BCVectorInterface > M_BCVectorList;

    //! map between the index number (second) in the vectors
    //!  (M_neighborList, M_inputBCList ...)
    //!  and its reference interface number (first).
    std::map<Int, UInt> M_referenceToIndexMap;

};
}
#endif
