//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the ETFESpace.

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ETFESPACE_HPP
#define ETFESPACE_HPP

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/ReferenceFEScalar.hpp>
#include <lifev/core/fem/ReferenceFEHdiv.hpp>
#include <lifev/core/fem/ReferenceFEHybrid.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/DOF.hpp>

#include <lifev/eta/fem/MeshGeometricMap.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

namespace LifeV
{



//! class ETFESpace  A light, templated version of the FESpace
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class represents a data structure for everything a finite element space
    requires to be defined:

    <ul>
        <li> The mesh
        <li> The reference element
        <li> The geometric map (between the reference element and mesh elements)
        <li> The degree of freedom numbering
        <li> The repartition of the degrees of freedom across the processors (algebraic map)
    </ul>

    It does not contain any information about the quadrature (unlike LifeV::FESpace), as this is not strictly needed for the definition of the space.

    This class is supposed to be constant during a simulation, so that it can
    be shared across the different structures using it.

    <b>Template parameters</b>

    <i> MeshType </i> The type of the mesh used.
    <i> MapType </i> The type of the algebraic map (for distributed computations, e.g. MapEpetra).
    <i> SpaceDim </i> The space of the domain definition
    <i> FieldDim </i> The dimension of the field (1 for a scalar FE, more for a vectorial one)

    <b>Template requirements</b>

    <i> MeshType </i> Same as for boost::shared_ptr
    <i> MapType </i> empty constructor; copy constructor; constructor using the reference FE, the mesh and a communicator;
    concatenation operator +=

*/
template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
class ETFESpace
{
public:

    //! @name Public Types
    //@{

    //! Typedef for the mesh
    typedef MeshType mesh_Type;

    //! Typedef for the map (algebraic)
    typedef MapType map_Type;

    //! Typedef for a pointer on the mesh
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    //! Typedef for a pointer on the communicator
    typedef typename map_Type::comm_ptrtype commPtr_Type;

    //@}


    //! @name Static constants
    //@{

    //! Dimension of the space
    enum {space_dim = SpaceDim};

    //! Dimension of the field of the FE space
    enum {field_dim = FieldDim};

    //@}


    //! @name Constructors, destructor
    //@{

    //! Full constructor using the mesh, the referenceFE, the mapping and the communicator
    /*!
      @param mesh Pointer on the mesh
      @param refFE The reference element for the finite element
      @param geoMap The geometric mapping
      @param commptr Pointer on the communicator to be used
     */
    ETFESpace (const meshPtr_Type& mesh,
               const ReferenceFE* refFE,
               const GeometricMap* geoMap,
               commPtr_Type& commptr);

    //! Constructor where the geometric mapping is guessed
    /*!
      In this constructor, the geometric map is guessed using the shape of the
      elements of the mesh.

      @param mesh Pointer on the mesh
      @param refFE The reference element for the finite element
      @param commptr Pointer on the communicator to be used
     */
    ETFESpace (const meshPtr_Type& mesh,
               const ReferenceFE* refFE,
               commPtr_Type& commptr);

    //! Full constructor using the partitioner of the mesh
    /*!
      @param meshPartitioner The partition of the mesh
      @param refFE The reference element for the finite element
      @param geoMap The geometric mapping
      @param commptr Pointer on the communicator to be used
     */
    ETFESpace (const MeshPartitioner<MeshType>& meshPartitioner,
               const ReferenceFE* refFE,
               const GeometricMap* geoMap,
               commPtr_Type& commptr);

    //! Full constructor using the partitioner of the mesh and a guessed geometric map
    /*!
      In this constructor, the geometric map is guessed using the shape of the
      elements of the mesh.

      @param meshPartitioner The partition of the mesh
      @param refFE The reference element for the finite element
      @param commptr Pointer on the communicator to be used
     */
    ETFESpace (const MeshPartitioner<MeshType>& meshPartitioner,
               const ReferenceFE* refFE,
               commPtr_Type& commptr);


    //! Copy constructor
    /*!
      @param otherSpace The finite element space to be copied
     */
    ETFESpace (const ETFESpace<MeshType, MapType, SpaceDim, FieldDim>& otherSpace);

    //! Destructor
    virtual ~ETFESpace()
    {}

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the mesh pointer
    /*!
      @return The pointer (shared) on the mesh
     */
    meshPtr_Type mesh() const
    {
        return M_mesh;
    }

    //! Getter for the reference FE
    /*!
      @return The reference FE of this space
     */
    const ReferenceFE& refFE() const
    {
        return *M_referenceFE;
    }

    //! Getter for the geometric mapping
    /*!
      @return The geometric mapping used for this space
     */
    const GeometricMap& geoMap() const
    {
        return *M_geometricMap;
    }

    //! Getter for the dof manager
    /*!
      @return The structure retaining the dof numbering
     */
    const DOF& dof() const
    {
        return *M_dof;
    }

    //! Getter for the algebraic map
    /*!
      @return The algebraic map
     */
    const MapType& map() const
    {
        return *M_map;
    }

    //! Getter for the algebraic map
    /*!
      @return The algebraic map
     */
    MapType& map()
    {
        return *M_map;
    }

    //! Getter for the dimension of the space (geometric, ambiant space)
    /*!
      @return The dimension of the space in which this FE space is defined
     */
    const UInt spaceDim() const
    {
        return space_dim;
    }

    //! Getter for the dimension of the field (scalar vs vectorial FE)
    /*!
      @return The dimension of the field represented.
     */
    const UInt fieldDim() const
    {
        return field_dim;
    }
    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor
    ETFESpace();

    //! Creates the map from the input
    void createMap (const commPtr_Type& commptr);

    //@}


    // Mesh member
    meshPtr_Type M_mesh;

    // Reference FE
    const ReferenceFE* M_referenceFE;

    // Geometric mapping
    const GeometricMap* M_geometricMap;

    // DoF manager
    DOF* M_dof;

    // Algebraic map
    MapType* M_map;
};


// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
ETFESpace (const meshPtr_Type& mesh, const ReferenceFE* refFE, const GeometricMap* geoMap, commPtr_Type& commptr)

    : M_mesh (mesh),
      M_referenceFE (refFE),
      M_geometricMap (geoMap),
      M_dof ( new DOF ( *M_mesh, *M_referenceFE ) ),
      M_map (new MapType() )
{
    createMap (commptr);
}

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
ETFESpace (const meshPtr_Type& mesh, const ReferenceFE* refFE, commPtr_Type& commptr)

    : M_mesh (mesh),
      M_referenceFE (refFE),
      M_geometricMap (&geometricMapFromMesh<MeshType>() ),
      M_dof ( new DOF ( *M_mesh, *M_referenceFE ) ),
      M_map (new MapType() )
{

    createMap (commptr);
}

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
ETFESpace (const MeshPartitioner<MeshType>& meshPartitioner,
           const ReferenceFE* refFE,
           const GeometricMap* geoMap,
           commPtr_Type& commptr)
    : M_mesh (meshPartitioner.meshPartition() ),
      M_referenceFE (refFE),
      M_geometricMap (geoMap),
      M_dof ( new DOF ( *M_mesh, *M_referenceFE ) ),
      M_map (new MapType() )
{
    createMap (commptr);
}

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
ETFESpace (const MeshPartitioner<MeshType>& meshPartitioner,
           const ReferenceFE* refFE,
           commPtr_Type& commptr)
    : M_mesh (meshPartitioner.meshPartition() ),
      M_referenceFE (refFE),
      M_geometricMap ( &geometricMapFromMesh<MeshType>() ),
      M_dof ( new DOF ( *M_mesh, *M_referenceFE ) ),
      M_map (new MapType() )
{
    createMap (commptr);
}

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
ETFESpace (const ETFESpace<MeshType, MapType, SpaceDim, FieldDim>& otherSpace)

    : M_mesh (otherSpace.M_mesh),
      M_referenceFE (otherSpace.M_referenceFE),
      M_geometricMap (otherSpace.M_geometricMap),
      M_dof (otherSpace.M_dof),
      M_map (otherSpace.M_map)
{}

template<typename MeshType, typename MapType, UInt SpaceDim, UInt FieldDim>
void
ETFESpace<MeshType, MapType, SpaceDim, FieldDim>::
createMap (const commPtr_Type& commptr)
{
    // get globalElements list from DOF
    typename MapType::mapData_Type mapData = this->M_dof->createMapData ( *this->M_mesh );
    // Create the map
    MapType map ( mapData, commptr );

    for ( UInt ii (0); ii < FieldDim; ++ii )
    {
        *M_map += map;
    }
}


} //Namespace LifeV

#endif //ETFESPACE_HPP
