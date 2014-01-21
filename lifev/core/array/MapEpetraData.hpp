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
    @brief MapEpetraData

    @author Antonio Cervone <ant.cervone@gmail.com>
    @maintainer Antonio Cervone <ant.cervone@gmail.com>

    @date 2013-04-19

    This class stores information needed to build a MapEpetra object
 */

#ifndef EPETRAMAPDATA_HPP
#define EPETRAMAPDATA_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/EnumMapEpetra.hpp>

namespace LifeV
{


//! MapEpetraData - todo
/*!
  The MapEpetra class provides a general interface for the Epetra_Map class of Trilinos.
 */
class MapEpetraData
{
public:

    //! @name Public Types
    //@{

    typedef std::vector<int> idList_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    MapEpetraData ();

    //! Constructor
    /*!
      @param uniqueMapSize size of the Unique map
      @param repeatedMapSize size of the Repeated map
    */
    MapEpetraData ( const UInt& uniqueMapSize, const UInt& repeatedMapSize );

    //! Destructor
    ~MapEpetraData() {}

    //@}

    //! @name Get methods
    //@{

    idList_Type const& idList ( MapEpetraType listType ) const;

    //@}

    //! @name Methods
    //@{

    //! Set
    template <typename list_Type>
    void set ( list_Type const& list, MapEpetraType listType );

    //! Show informations about the map data
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

private:

    idList_Type M_uniqueList;
    idList_Type M_repeatedList;

}; // class MapEpetraData

inline MapEpetraData::MapEpetraData ( const UInt& uniqueMapSize, const UInt& repeatedMapSize ) :
    M_uniqueList ( uniqueMapSize ),
    M_repeatedList ( repeatedMapSize )
{}

template <typename list_Type>
inline void MapEpetraData::set ( list_Type const& list, MapEpetraType listType )
{
    switch ( listType )
    {
        case Unique:
            std::copy ( list.begin(), list.end(), M_uniqueList.begin() );
            break;
        case Repeated:
            std::copy ( list.begin(), list.end(), M_repeatedList.begin() );
            break;
        default:
            ERROR_MSG ( "wrong map type" );
    }
}

inline MapEpetraData::idList_Type const& MapEpetraData::idList ( MapEpetraType listType ) const
{
    switch ( listType )
    {
        case Unique:
            return M_uniqueList;
            break;
        case Repeated:
            return M_repeatedList;
            break;
        default:
            ERROR_MSG ( "wrong map type" );
    }
    return M_uniqueList;
}

} // namespace LifeV

#endif // EPETRAMAPDATA_HPP
