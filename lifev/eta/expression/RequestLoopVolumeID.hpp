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
    @brief File containing the class to loop over a certain set of volumes

    @author Paolo Tricerri<paolo.tricerri@epfl.ch>

    @date 07-2011
 */


#ifndef REQUEST_LOOP_VOLUME_ID_HPP
#define REQUEST_LOOP_VOLUME_ID_HPP

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>


namespace LifeV
{

namespace ExpressionAssembly
{

template <typename MeshType>
class RequestLoopVolumeID
{
public:

    typedef typename MeshType::element_Type             elementMesh_Type;
    typedef std::vector<elementMesh_Type const* >       vectorVolumes_Type;
    typedef boost::shared_ptr<vectorVolumes_Type>       vectorVolumesPtr_Type;

    //! @name Constructors & Destructor
    //@{

    //! Simple constructor with a shared_ptr on the mesh
	RequestLoopVolumeID(const vectorVolumes_Type& elementList)
        : M_selectedVolumes( elementList )
    {}

    //! Copy constructor
	RequestLoopVolumeID(const RequestLoopVolumeID& loop)
        : M_selectedVolumes(loop.M_selectedVolumes)
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
	const vectorVolumes_Type& selectedVolumes() const { return M_selectedVolumes; }

    //@}

private:


    //! @name Private Methods
    //@{

    //! No empty constructor
	RequestLoopVolumeID();

    //@}

    // Pointer on the list of the volumes
	vectorVolumesPtr_Type M_selectedVolumes;

};


//! elements - A helper method to trigger the loop on the elements of a mesh
/*!
    @author Paolo Tricerri

    This method simply returns a RequestLoopElement type, so that one does not
    need to explicitly call this type and the type of the mesh used.

    <b> Template parameters </b>

    <i>MeshType</i>: The type of the mesh.

    <b> Template requirements </b>

    <i>MeshType</i>: See in LifeV::RequestLoopElement

 */
template< typename MeshType >
RequestLoopVolumeID<MeshType>
integrationOverSelectedVolumes(const typename RequestLoopVolumeID<MeshType>::vectorVolumes_Type& elementList)
{
	return RequestLoopVolumeID<MeshType>( elementList );
};


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
