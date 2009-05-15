/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file dofOneD.hpp
  \author Vincent Martin
  \date 07/2004
  \version 1.0

  \brief This file contains a very basic one dimensional dof class.
  The principle is to mimic the behaviour of the "real"
  Dof class of dof.hpp, without creating it.
  Thus one can use, with this DofOneD class, the assemblage
  function which already exists in assemb.hpp.
  To be used for instance with a mesh handler like BasicOneDMesh.

  In fact I can't use the "dof::update()" that requires
  a format of mesh that I did not build. So I turn around
  this difficulty, and create my own simple DofOneD class.

  \TODO
  The next improvement would be to construct in lifeV a mesh
  handler that could also take care of 1D meshes,
  and then to use the "real" Dof.

*/

#ifndef _DOFONED_H_
#define _DOFONED_H_


#include <life/lifecore/life.hpp>

namespace LifeV
{

/*!
  \class DofOneD: (very) simple class for one D dof
  (mimics the "real" Dof class of dof.hpp)

*/
class DofOneD
{
public:
    //! Constructor
    DofOneD( const UInt& totalDof ): _M_totalDof(totalDof) {}

    //! The total number of Dof
    INLINE UInt numTotalDof() const {return _M_totalDof;}

    /*!
      Returns the global numbering of a DOF, given an element and the local numbering
      \param ELId the element ID
      \param localNode the local DOF numbering (starting from 1)
      \return The numbering of the DOF (starting from 1)

      BEWARE: works only for 1D mesh with increasing numbering
      of nodes and elements. (The canonic numbering!)
      Has only been tested for P1Seg elements! (VM)
    */
    ID localToGlobal(const ID ElId, const ID localNode) const {
        return ElId + localNode;
    }

protected:
    //! number of total dof
    UInt _M_totalDof;
};

}
#endif
