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
    @brief File containing the VectorBlockStructure class

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>
    @date 21-08-2012
 */

#ifndef _VECTORBLOCKSTRUCTURE_HPP_
#define _VECTORBLOCKSTRUCTURE_HPP_ 1

#include <boost/shared_ptr.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>

namespace LifeV
{

//! VectorBlockStructure - class representing the structure of a vector
/*!
  @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>

 */
class VectorBlockStructure
{
public:

    //! @name Public Types
    //@{

    //! Type of data stored
    typedef Real data_Type;

    //! Type of the map to be used
    typedef MapEpetra map_Type;

    //! Type of the MapVector to be used with this class
    typedef MapVector<map_Type> mapVector_Type;

    //! Type of the map (Unique/Repeated)
    typedef MapEpetraType mapType_Type;

    //! Combine mode
    typedef Epetra_CombineMode combine_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Default constructor
    VectorBlockStructure();

    //! Constructor with the monolithic map
    VectorBlockStructure ( const map_Type& map );

    //! Construction with a map
    /*!
      With this constructor, the block structure is automatically deduced from the maps in the
      vector. The monolithic map and vectors are also built by concatenating the different maps
      in the vector.
     */
    VectorBlockStructure ( const mapVector_Type& mapVector );

    //! Copy constructor
    VectorBlockStructure ( const VectorBlockStructure& blockStructure );

    //! Destructor
    ~VectorBlockStructure() {}

    //@}


    //! @name Set Methods
    //@{

    /*! Set the size of the blocks of the vector
     *  @param blockSizes Sizes of the blocks
     */
    void setBlockStructure ( const std::vector<UInt>& blockSizes );

    //! Reset the block structure using the blocks of a vector of map
    /*
      The resulting block structure is symmetric (same block structure
      in the rows and in the columns).

      This method does not involve big computation overhead. Remark that
      it is not possible to change the size of the block vector
      through this method, nor its map.
      @param The MapVector containing the maps
     */
    void setBlockStructure ( const mapVector_Type& mapVector );

    /*! Set the block structure using a block structure
     *  @param blockStructure Structure of the vector
     */
    void setBlockStructure ( const VectorBlockStructure& blockStructure );

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the size of the block index
    /*!
      @param index Index of the block
      @return size of the index-th block
     */
    UInt blockSize ( const UInt& index ) const
    {
        ASSERT ( index < M_blockSize.size(), "Invalid block index" );
        return M_blockSize[index];
    }

    /*!
      @param index Index of the block
      @return index of the first entry in the index-th block
     */
    UInt blockFirstIndex ( const UInt& index ) const
    {
        ASSERT ( index < M_blockFirstIndex.size(), "Invalid block index" );
        return M_blockFirstIndex[index];
    }

    /*!
       @return Number of blocks
     */
    UInt numBlocks() const
    {
        return M_blockSize.size();
    }

    /*!
       @return Number of blocks
     */
    UInt totalSize() const
    {
        return M_totalSize;
    }

    //@}

private:

    std::vector<UInt> M_blockSize;
    std::vector<UInt> M_blockFirstIndex;
    UInt              M_totalSize;

};

} // Namespace LifeV

#endif /* _VECTORBLOCKSTRUCTURE_HPP_ */
