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
    @brief File containing the MatrixBlockStructure class

    @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>
    @date 21-08-2012
 */

#ifndef _MATRIXBLOCKSTRUCTURE_HPP_
#define _MATRIXBLOCKSTRUCTURE_HPP_ 1

#include <boost/shared_ptr.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/VectorBlockStructure.hpp>

namespace LifeV
{

//! MatrixBlockStructure - class representing the structure of a vector
/*!
  @author Gwenol Grandperrin <gwenol.grandperrin@gmail.com>

 */
class MatrixBlockStructure
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

    //@}


    //! @name Constructor & Destructor
    //@{
    //! Default constructor
    MatrixBlockStructure();

    //! Constructor with the monolithic maps
    MatrixBlockStructure ( const map_Type& rowMap,
                           const map_Type& columnMap );

    //! Constructor with a unique monolithic map
    MatrixBlockStructure ( const map_Type& map );

    //! Construction with a map
    /*!
      With this constructor, the block structure is automatically deduced from the maps in the
      vectors. The monolithic map and vectors are also built by concatenating the different maps
      in the vector.
     */
    MatrixBlockStructure ( const mapVector_Type& rowMapVector,
                           const mapVector_Type& columnMapVector );

    //! Construction with a map
    /*!
      With this constructor, the block structure is automatically deduced from the maps in the
      vector. The monolithic map and vectors are also built by concatenating the different maps
      in the vector.
     */
    MatrixBlockStructure ( const mapVector_Type& mapVector );

    //! Construction with row and column structure
    MatrixBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                           const VectorBlockStructure& columnsBlockStructure );

    //! Construction with vector structure
    MatrixBlockStructure ( const VectorBlockStructure& vectorStructure );

    //! Copy constructor
    MatrixBlockStructure ( const MatrixBlockStructure& blockStructure );

    //! Destructor
    ~MatrixBlockStructure();

    //@}


    //! @name Set Methods
    //@{

    /*! Set the size of the blocks of the matrix
     *  @param blockNumRows Number of rows in the blocks
     *  @param blockNumColumns Number of columns in the blocks
     */
    void setBlockStructure ( const std::vector<UInt>& blockNumRows,
                             const std::vector<UInt>& blockNumColumns );

    /*! Set the size of the blocks of the matrix
     *  @param blocksSize Number of rows/columns in the blocks
     */
    void setBlockStructure ( const std::vector<UInt>& blocksSize );

    /*!
     * Set the size of the blocks of the matrix using the maps stored in the vector
     * of map. The resulting block structure in symmetric (same blocks in the rows
     * and in the columns).
     *
     * This method does not involve large computations. The global size and the map
     * of the matrix cannot be changed with this method (and the block structure has
     * to be compatible with the global size).
     */
    void setBlockStructure ( const MapVector<MapEpetra>& rowMapVector,
                             const MapVector<MapEpetra>& columnMapVector );

    /*!
     * Set the size of the blocks of the matrix using the maps stored in the vector
     * of map. The resulting block structure in symmetric (same blocks in the rows
     * and in the columns).
     *
     * This method does not involve large computations. The global size and the map
     * of the matrix cannot be changed with this method (and the block structure has
     * to be compatible with the global size).
     */
    void setBlockStructure ( const MapVector<MapEpetra>& mapVector );

    //! Set the block structure from a matrix structure object
    void setBlockStructure ( const MatrixBlockStructure& blockStructure );

    //! Set the block structure from row and column structures
    void setBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                             const VectorBlockStructure& columnsBlockStructure );

    //! Set the block structure from row and column structures
    void setBlockStructure ( const VectorBlockStructure& vectorBlockStructure );

    //@}


    //! @name Get Methods
    //@{
    //! Returns the number of rows of the block
    /*!
     * @param rowIndex Row index of the block
     */
    UInt blockNumRows ( const UInt& rowIndex ) const;

    //! Returns the number of columns of the block
    /*!
     * @param columnIndex Column index of the block
     */
    UInt blockNumColumns ( const UInt& columnIndex ) const;

    /*!
      @param index Index of the block
      @return index of the first entry in the index-th block
     */
    UInt rowBlockFirstIndex ( const UInt& index ) const;

    /*!
      @param index Index of the block
      @return index of the first entry in the index-th block
     */
    UInt columnBlockFirstIndex ( const UInt& index ) const;

    //! Number of row blocks
    UInt numRowBlocks() const;

    //! Number of column blocks
    UInt numColumnBlocks() const;

    //! Number of rows
    UInt numRows() const;

    //! Number of rows
    UInt numColumns() const;

    //! Get the rows block structure
    const VectorBlockStructure& rowsBlockStructure() const;

    //! Get the columns block structure
    const VectorBlockStructure& columnsBlockStructure() const;

    //@}

private:

    VectorBlockStructure M_rowsBlockStructure;
    VectorBlockStructure M_columnsBlockStructure;

};

} // Namespace LifeV

#endif /* _MATRIXBLOCKSTRUCTURE_HPP_ */
