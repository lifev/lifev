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


#ifndef LOOP_REQUEST_HPP
#define LOOP_REQUEST_HPP

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

//! RequestLoopElement - The class to request a loop over the elements of a mesh.
/*!
    @author Samuel Quinodoz

    The only role of this class is to trigger, at compile time (through its type), loops
    on the elements of the mesh stored internally.

    <b> Template parameters </b>

    <i>MeshType</i>: The type of the mesh.

    <b> Template requirements </b>

    <i>MeshType</i>: none

 */
template <typename MeshType>
class RequestLoopElement
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Simple constructor with a shared_ptr on the mesh
    RequestLoopElement (const boost::shared_ptr<MeshType>& mesh) : M_mesh (mesh) {}

    //! Copy constructor
    RequestLoopElement (const RequestLoopElement& loop) : M_mesh (loop.M_mesh) {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
    const boost::shared_ptr<MeshType>& mesh() const
    {
        return M_mesh;
    }

    //@}

private:


    //! @name Private Methods
    //@{

    //! No empty constructor
    RequestLoopElement();

    //@}

    // Pointer on the mesh
    boost::shared_ptr<MeshType> M_mesh;
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
RequestLoopElement<MeshType>
elements (const boost::shared_ptr<MeshType>& mesh)
{
    return RequestLoopElement<MeshType> (mesh);
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
