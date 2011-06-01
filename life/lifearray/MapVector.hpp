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

#include <life/lifearray/MapEpetra.hpp>

#include <life/lifecore/LifeV.hpp>

#include <vector>

namespace LifeV
{

//! This class is used to store maps that will be used for block defined problems.
/*!
  This class consists in a collection of maps (e.g. MapEpetra) to represent a problem made of several blocks.

  When a system is made of several blocks (like a Stokes problem, with the velocity block and pressure block), one might want to use a block matrix and block vector. This class represents the "map" to construct block structures.

  There are two operators to work with the maps: the "+" and the "|" operators. The "+" can be used to fusion two maps to create a bigger map (see the documentation of the related maps). The "|" is used to separate two maps for two different blocks. Starting from two maps (or a map and a vector of maps), it creates a vector of maps.

  Between the operators "+" and "|", the "+" has a higher priority, meaning that all the "+" are evaluated before the "|" are. However, it is strongly recommanded to use brackets (some compilers might issue warnings if brackets are not used).

  The only valid template argument is for now MapEpetra, for the Epetra framework of Trilinos.
*/

template< typename MapType>
class MapVector
{

public:

    //! @name Public Types
    //@{

	typedef MapType map_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Destructor
    ~MapVector(){}


	//! Empty constructor
	MapVector()
        : M_vector(){}


	//! Copy constructor
	MapVector(const MapVector<map_Type>& otherMapVector)
        : M_vector(otherMapVector.M_vector) {}


	//! Constructor with one map
	/*!
      This constructor builds a MapVector with only
      the map given in arguement.
	*/
	MapVector(const map_Type& map)
        : M_vector(1,map) {}


	//! Constructor with two maps
	/*!
      The MapVector is filled with the two maps given in argument.
	*/
	MapVector(const map_Type& map1, const map_Type& map2)
        : M_vector() {
		M_vector.push_back(map1);
		M_vector.push_back(map2);
	}


	//! Concatenation constructor
	/*!
      This constructor copies the maps in vector and adds the map to it.
	*/
	MapVector(const MapVector<map_Type>& vector, const map_Type& map)
        : M_vector(vector.M_vector)
	{
		M_vector.push_back(map);
	}

    //@}


    //! @name Operators
    //@{

	//! Copy operator
	MapVector<map_Type> operator=(const MapVector<map_Type>& vector)
	{
		M_vector=vector.M_vector;
		return *this;
	}

    //@}


    //! @name Methods
    //@{


	//! Returns the number of maps stored
	inline UInt nbMap() const {return M_vector.size(); }


    //! Display internal state
    inline void showMe( std::ostream& output = std::cout) const
    {
        std::cout << " Number of map stored : " << M_vector.size() << std::endl;
    }

    //@}


    //! @name Get Methods
    //@{


	//! Getter for the ith map stored
	inline const map_Type& map(const UInt& i) const
	{
		ASSERT( i< M_vector.size() ,"Index out of bound, no map to return");
		return M_vector[i];
	}

    //@}



private:

    // The vector containing the maps
	std::vector<map_Type> M_vector;
};


//! Typedef for the MapVector containing MapEpetra
typedef MapVector<MapEpetra> MapEpetraVector;


//! Juxtaposition operator for the two maps
template<typename MapType>
MapVector<MapType>
operator|(const MapType& map1, const MapType& map2)
{
	return MapVector<MapType>(map1,map2);
}

//! Juxtaposition operator for a vector of maps and a map
template<typename MapType>
MapVector<MapType>
operator|(const MapVector<MapType>& vector, const MapType& map)
{
	return MapVector<MapType>(vector,map);
}


} // End of the namespace

#endif
