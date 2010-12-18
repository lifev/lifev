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
    @brief Base Class for interfacing dofs between two meshes

    @author M.A. Fernandez and V. Martin
    @date 00-11-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This file contains the class which may be used to update and hold the connections between the dof
    on two matching meshes.
 */

#ifndef _DOFINTERFACEBASE_HH
#define _DOFINTERFACEBASE_HH

#include <life/lifecore/life.hpp>

#include <map>

namespace LifeV
{
/*!
  \class DofInterfaceBase

  Base class which holds the connections of the dof in two matching meshes
  The dof mapping (STL map) is set in the derived classes.
  Each derived class depends on the type of interface you have.

  Method getInterfaceDof gives the connections.

  To be deleted:
  buildInverse
*/
class DofInterfaceBase
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Default Constructor
    DofInterfaceBase();

    //! Virtual Destructor
    virtual ~DofInterfaceBase(){};

    //@}


    //! @name Methods
    //@{

    //! This method returns the corresponding dof number of the mesh2 at the interface
    //! for a specific dof number at the interface in mesh1
    /*!
      \param i a dof number in mesh1
    */
    ID getInterfaceDof( const ID& i ) const;

    //! This method says whether a specific dof number at the interface in mesh1 is on this processor
    /*!
      \param i a dof number in mesh1
    */
    bool isMyInterfaceDof( const ID& i ) const;

    //! output
    std::ostream& showMe( bool verbose = false, std::ostream& out = std::cout ) const;

    //! Makes this DofInterfaceBase to be the inverse map as the one defined by dofBase.
    void buildInverse( const DofInterfaceBase& dofBase);

    //@}

    //! @name Set Methods
    //@{

    //! Set value to be associated to key
    void set(const ID& key,const ID& value);

    //@}

    //! @name Get Methods
    //@{

    //! This method returns the number of dof that live on the interface
    const size_t& nbInterfaceDof() const;

    //! Return the correspondance map
    const std::map<ID, ID> & localDofMap() {return M_localDofMap;}

    //@}

protected:

    //!  STL map container which holds the connections between Dof at the interface
    std::map<ID, ID> M_localDofMap;

};

}
#endif

