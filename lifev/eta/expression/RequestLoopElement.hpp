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
    RequestLoopElement ( const std::shared_ptr<MeshType>& mesh, const UInt regionFlag = 0,
                         const UInt numVolumes = 0, const UInt * volumeElements = nullptr,
                         const bool integrateOnSubdomains = false )
        : M_mesh ( mesh ), M_regionFlag ( regionFlag ), M_volumeElements( volumeElements ), M_numVolumes( numVolumes ),
          M_integrateOnSubdomains( integrateOnSubdomains )
        { }

    //! Copy constructor
    RequestLoopElement (const RequestLoopElement& loop)
        : M_mesh ( loop.M_mesh ), M_regionFlag ( loop.M_regionFlag ), M_numVolumes( loop.M_numVolumes ),
          M_volumeElements( loop.M_volumeElements ), M_integrateOnSubdomains( loop.M_integrateOnSubdomains )
          { }

    ~RequestLoopElement()
    {
        M_volumeElements = nullptr;
    }
    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
    const std::shared_ptr<MeshType>& mesh() const
    {
        return M_mesh;
    }

    //! Getter for the flag of the region of integration
    const UInt regionFlag() const
    {
        return M_regionFlag;
    }

    const UInt numVolumes() const
    {
        return M_numVolumes;
    }

    //! Getter for the flag of the region of integration
    const UInt * getElementsRegionFlag() const
    {
        return M_volumeElements;
    }

    const bool getIfSubDomain() const
    {
        return M_integrateOnSubdomains;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    RequestLoopElement();

    //@}

    // Pointer on the mesh
    std::shared_ptr<MeshType> M_mesh;

    // Data for integration on one subRegion, flag and elements on which iterate
    const UInt                  M_regionFlag;
    const UInt                  M_numVolumes;
    const UInt *                M_volumeElements;
    const bool                  M_integrateOnSubdomains;

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
elements (const std::shared_ptr<MeshType>& mesh, const UInt flag = 0, const UInt numVolumes = 0,
          const UInt * volumeElements = nullptr, const bool subDomain = false )
{

    return RequestLoopElement<MeshType> ( mesh, flag, numVolumes, volumeElements, subDomain );
}


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
