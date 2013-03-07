//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
   @file MapVector.hpp
   @brief The file contains the MapVector class

   @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2011-04-31
 */


#ifndef MAP_VECTOR_HPP
#define MAP_VECTOR_HPP

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! This class is used to store maps that will be used for block defined problems.
/*!
  This class consists in a collection of maps (e.g. LifeV::MapEpetra) to represent a problem made of several blocks.

  When a system is made of several blocks (like a Stokes problem, with the velocity block and pressure block), one might want to use a block matrix and block vector. This class represents the "map" to construct block structures.

  There are two operators to work with the maps: the "+" and the "|" operators. The "+" can be used to fusion two maps to create a bigger map (see the documentation of the related maps). The "|" is used to separate two maps for two different blocks. Starting from two maps (or a map and a vector of maps), it creates a vector of maps.

  Between the operators "+" and "|", the "+" has a higher priority, meaning that all the "+" are evaluated before the "|" are. However, it is strongly recommanded to use brackets (some compilers might issue warnings if brackets are not used).

  <b> Example of valid template argument </b>

  LifeV::MapEpetra

*/

template< typename MapType>
class MapVector
{

public:

    //! @name Public Types
    //@{

    //! Type of the map
    typedef MapType map_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    MapVector();

    //! Copy constructor
    /*!
      @param otherMapVector The vector of maps to copy
     */
    MapVector (const MapVector<map_Type>& otherMapVector);

    //! Constructor with one map
    /*!
      This constructor builds a MapVector with only
      the map given in arguement.
      @param map The map to be put in the vector
    */
    explicit MapVector (const map_Type& map);

    //! Constructor with two maps
    /*!
      The MapVector is filled with the two maps given in argument.
      @param map1 The first map of the vector
      @param map2 The second map of the vector
    */
    MapVector (const map_Type& map1, const map_Type& map2);

    //! Concatenation constructor
    /*!
      This constructor copies the maps in vector and adds the map to it.
      @param vector The vector to copy
      @param map The map to add to vector
    */
    MapVector (const MapVector<map_Type>& vector, const map_Type& map);

    //! Destructor
    ~MapVector();

    //@}


    //! @name Operators
    //@{

    //! Assignement operator
    /*
      @param vector The MapVector to copy
      @return This MapVector
     */
    MapVector<map_Type> operator= (const MapVector<map_Type>& vector);

    //! Access to the ith map (read only)
    /*
      @param i The index of the map
      @return The ith map stored
     */
    const map_Type& operator[] (UInt i) const;

    //! Access to the ith map (read-write)
    /*
      @param i The index of the map
      @return The ith map stored
     */
    map_Type& operator[] (UInt i);

    //! Juxtaposition operator for a vector of maps and a map
    /*
      @param map The map to be added
      @return MapVector resulting of the concatenation of this MapVector and the map
      given in argument
     */
    MapVector<map_Type>
    operator| (const map_Type& map) const;

    //@}


    //! @name Methods
    //@{

    //! Returns the number of maps stored
    UInt nbMap() const;

    //! Returns the size of the ith map stored (global number of ids)
    UInt mapSize (UInt i) const;

    //! Display internal state
    void showMe ( std::ostream& output = std::cout) const;

    //! Return the map made by concatenating all the maps of the vector
    map_Type totalMap() const;

    //! Add a map into this vector (at the end of the vector)
    void addMap (const map_Type& newMap);

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the ith map stored
    const map_Type& map (const UInt& i) const;

    //@}


private:

    // The vector containing the maps
    std::vector<map_Type> M_vector;
};


// ===================================================
// Constructors & Destructor
// ===================================================

template< typename MapType>
MapVector<MapType>::
MapVector()
    : M_vector()
{}

template< typename MapType>
MapVector<MapType>::
MapVector (const MapVector<map_Type>& otherMapVector)
    : M_vector (otherMapVector.M_vector)
{}

template< typename MapType>
MapVector<MapType>::
MapVector (const map_Type& map)
    : M_vector (1, map)
{}

template< typename MapType>
MapVector<MapType>::
MapVector (const map_Type& map1, const map_Type& map2)
    : M_vector()
{
    M_vector.push_back (map1);
    M_vector.push_back (map2);
}

template< typename MapType>
MapVector<MapType>::
MapVector (const MapVector<map_Type>& vector, const map_Type& map)
    : M_vector (vector.M_vector)
{
    M_vector.push_back (map);
}

template< typename MapType>
MapVector<MapType>::
~MapVector()
{}

// ===================================================
// Operators
// ===================================================

template< typename MapType>
MapVector<MapType>
MapVector<MapType>::
operator= (const MapVector<map_Type>& vector)
{
    M_vector = vector.M_vector;
    return *this;
}

template< typename MapType>
const MapType&
MapVector<MapType>::
operator[] (UInt i) const
{
    ASSERT ( i < M_vector.size() , "Index out of bound, no map to access (Read only)");
    return M_vector[i];
}

template< typename MapType>
MapType&
MapVector<MapType>::
operator[] (UInt i)
{
    ASSERT ( i < M_vector.size() , "Index out of bound, no map to access (Read write)");
    return M_vector[i];
}

template< typename MapType>
MapVector<MapType>
MapVector<MapType>::
operator| (const map_Type& map) const
{
    return MapVector<map_Type> (*this, map);
}

// ===================================================
// Methods
// ===================================================

template< typename MapType>
UInt
MapVector<MapType>::
nbMap() const
{
    return M_vector.size();
}

template< typename MapType>
UInt
MapVector<MapType>::
mapSize (UInt i) const
{
    ASSERT ( i < M_vector.size() , "Index out of bound, no map to return");
    return M_vector[i].mapSize();
}

template< typename MapType>
void
MapVector<MapType>::
showMe ( std::ostream& output) const
{
    output << " Number of map stored : " << M_vector.size() << std::endl;
}

template< typename MapType>
MapType
MapVector<MapType>::
totalMap() const
{
    ASSERT ( M_vector.size() != 0 , "No map to concatenate for the total map");
    map_Type total (M_vector[0]);

    for (UInt i (1); i < M_vector.size(); ++i)
    {
        total += M_vector[i];
    }
    return total;
}

template< typename MapType>
void
MapVector<MapType>::
addMap (const map_Type& newMap)
{
    M_vector.push_back (newMap);
}

// ===================================================
// Get Methods
// ===================================================

template< typename MapType>
const MapType&
MapVector<MapType>::
map (const UInt& i) const
{
    ASSERT ( i < M_vector.size() , "Index out of bound, no map to return");
    return M_vector[i];
}


} // End of the namespace

#endif
