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

template <typename MeshType, typename FunctorType>
class RequestLoopVolumeID
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Simple constructor with a shared_ptr on the mesh
	RequestLoopVolumeID(const boost::shared_ptr<MeshType>& mesh,
                        const boost::shared_ptr<FunctorType>& selector)
        : M_mesh( mesh ), M_functorSelection( selector )
    {}

    //! Copy constructor
	RequestLoopVolumeID(const RequestLoopVolumeID& loop)
    : M_mesh(loop.M_mesh), M_functorSelection( loop.M_functorSelection )
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
	const boost::shared_ptr<MeshType>& mesh() const { return M_mesh; }

    //! Getter for the identifier
    const boost::shared_ptr<FunctorType>& functorSelection() const { return M_functorSelection; }

    //@}

private:


    //! @name Private Methods
    //@{

    //! No empty constructor
	RequestLoopVolumeID();

    //@}

    // Pointer on the mesh
	boost::shared_ptr<MeshType> M_mesh;

    // Pointer on the mesh
	boost::shared_ptr<FunctorType> M_functorSelection;


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
template< typename MeshType, typename FunctorType >
RequestLoopVolumeID<MeshType, FunctorType>
integrationOverSelectedVolumes(const boost::shared_ptr<MeshType>& mesh,
                               const boost::shared_ptr<FunctorType>& functorSelection)
{
	return RequestLoopVolumeID<MeshType, FunctorType>( mesh, functorSelection );
};


} // Namespace ExpressionAssembly

} // Namespace LifeV

#endif
