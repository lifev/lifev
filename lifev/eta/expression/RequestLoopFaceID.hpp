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
    @brief File containing the classes to determine the kind of loops to be performed

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    @date 07-2011
 */


#ifndef REQUEST_LOOP_FACE_ID_HPP
#define REQUEST_LOOP_FACE_ID_HPP

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

/*!
  This class should be better named RequestLoopBoundaryFaceID

    @author Samuel Quinodoz
 */
template <typename MeshType>
class RequestLoopFaceID
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Simple constructor with a shared_ptr on the mesh
    RequestLoopFaceID (const std::shared_ptr<MeshType>& mesh, const UInt boundaryID)
        : M_mesh (mesh), M_boundaryIdentifier (boundaryID)
    {}

    //! Copy constructor
    RequestLoopFaceID (const RequestLoopFaceID& loop)
        : M_mesh (loop.M_mesh), M_boundaryIdentifier (loop.M_boundaryIdentifier)
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
    const std::shared_ptr<MeshType>& mesh() const
    {
        return M_mesh;
    }

    //! Getter for the identifier
    const UInt& id() const
    {
        return M_boundaryIdentifier;
    }

    //@}

private:


    //! @name Private Methods
    //@{

    //! No empty constructor
    RequestLoopFaceID();

    //@}

    // Pointer on the mesh
    std::shared_ptr<MeshType> M_mesh;

    const UInt M_boundaryIdentifier;
};


//! elements - A helper method to trigger the loop on the elements of a mesh
/*!
    @author Samuel Quinodoz

    This method simply returns a RequestLoopElement type, so that one does not
    need to explicitly call this type and the type of the mesh used.

    <b> Template parameters </b>

    <i>MeshType</i>: The type of the mesh.

    <b> Template requirements </b>

    <i>MeshType</i>: See in LifeV::RequestLoopElement

 */
template< typename MeshType >
RequestLoopFaceID<MeshType>
boundary (const std::shared_ptr<MeshType>& mesh, const UInt id)
{
    return RequestLoopFaceID<MeshType> (mesh, id);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
