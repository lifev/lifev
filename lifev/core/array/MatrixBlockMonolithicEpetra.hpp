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
   @file MatrixBlockMonolithicEpetra.hpp
   @brief The file contains the MatrixBlockMonolithicEpetra class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIX_BLOCK_MONOLITHIC_EPETRA_HPP_
#define _MATRIX_BLOCK_MONOLITHIC_EPETRA_HPP_

#include <boost/shared_ptr.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixBlockMonolithicEpetraView.hpp>

namespace LifeV
{

//! MatrixBlockMonolithicEpetra - class of block matrix
/*!
  @author Gwenol Grandperrin
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  The MatrixBlockMonolithicEpetra class contains data related
  to block matrix. It is an extension to MatrixEpetra where data about blocks have
  been set. For an introduction to the block structures in LifeV, see
  \ref BlockAlgebraPage "this page".

  There are mainly two ways to define a MatrixBlockMonolithicEpetra:
  <ul>
  <li> Construct it using the same syntax as for LifeV::MatrixEpetra and the use
  a setter for the structure.
  <li> Construct it directly with the maps of the blocks.
  </ul>
  Both ways are equivalent.

  To access the blocks, one uses then the blockView or block methods.

 */
template <typename DataType>
class MatrixBlockMonolithicEpetra : public MatrixEpetra<DataType>
{
public:

    typedef MatrixBlockMonolithicEpetraView<DataType> block_type;
    typedef std::shared_ptr<block_type> block_ptrType;

    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlockMonolithicEpetra ( const MapEpetra& map, int numEntries = 50);

    //! Block constructor
    /*!
    This is the most complete constructor, as it builds the whole block
    structure of the matrix, using the maps stored in the vector.
    */
    MatrixBlockMonolithicEpetra ( const MapVector<MapEpetra>& vector, int numEntries = 50);

    //! Casting constructor
    MatrixBlockMonolithicEpetra ( const MatrixEpetra<DataType>& matrix );

    //! Copy constructor
    MatrixBlockMonolithicEpetra ( const MatrixBlockMonolithicEpetra& matrix );

    //! Destructor
    ~MatrixBlockMonolithicEpetra();

    //@}


    //! @name  Set Methods
    //@{

    /*! Set the size of the blocks of the matrix
     *  @param blockNumRows Number of rows in the blocks
     *  @param blockNumColumns Number of columns in the blocks
     */
    void setBlockStructure ( const std::vector<UInt>& blockNumRows,
                             const std::vector<UInt>& blockNumColumns );

    /*! Set the size of the blocks of the matrix using the maps stored in the vector
      of map. The resulting block structure in symmetric (same blocks in the rows
      and in the columns).

      This method does not involve large computations. The global size and the map
      of the matrix cannot be changed with this method (and the block structure has
      to be compatible with the global size).

     */
    void setBlockStructure (const MapVector<MapEpetra>& mapVector);

    //@}


    //! @name  Get Methods
    //@{

    //! Returns the number of rows of the block
    /*!
     * @param rowIndex Row index of the block
     */
    UInt blockNumRows (const UInt& rowIndex) const;

    //! Returns the number of columns of the block
    /*!
     * @param columnIndex Column index of the block
     */
    UInt blockNumColumns (const UInt& columnIndex) const;

    //! Returns the MatrixBlockMonolithicEpetraView at the (rowIndex,colIndex) location
    /*!
     * @param rowIndex Row position of the block in the matrix
     * @param columnIndex Column position of the block in the matrix
     * @param mbv MatrixBlockMonolithicEpetraView to be filled
     */
    void blockView ( const UInt& rowIndex,
                     const UInt& columnIndex,
                     block_type& mbv);

    //! Returns the block (rowIndex,columnIndex) of the matrix
    /*!
      Remark that the returned block is a shared pointer. This is
      a limitation due to the non-copiability of the blocks.
      @param rowIndex The index of the block w.r. to the block structure
      @param columnIndex The index of the block w.r. to the block structure
     */
    block_ptrType block (const UInt& rowIndex, const UInt& columnIndex);

    //@}

private:

    std::vector<UInt> M_blockNumRows;
    std::vector<UInt> M_blockNumColumns;
    std::vector<UInt> M_blockFirstRows;
    std::vector<UInt> M_blockFirstColumns;
};

// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================
template <typename DataType>
MatrixBlockMonolithicEpetra<DataType>::MatrixBlockMonolithicEpetra (const MapEpetra& map, int numEntries) :
    MatrixEpetra<DataType> (map, numEntries),
    M_blockNumRows (std::vector<UInt> (1, map.map (Unique)->NumGlobalElements() ) ),
    M_blockNumColumns (std::vector<UInt> (1, map.map (Unique)->NumGlobalElements() ) ),
    M_blockFirstRows (std::vector<UInt> (1, 0) ),
    M_blockFirstColumns (std::vector<UInt> (1, 0) )
{
}



template <typename DataType>
MatrixBlockMonolithicEpetra<DataType>::MatrixBlockMonolithicEpetra ( const MapVector<MapEpetra>& vector, int numEntries)
    :
    MatrixEpetra<DataType> (vector.totalMap() )
    //MatrixEpetra<DataType>( typename MatrixEpetra<DataType>::matrix_ptrtype())
{
    ASSERT ( vector.nbMap() > 0 , "Map vector empty, impossible to construct a MatrixBlockMonolithicEpetra!");

    MapEpetra myMap (vector.map (0) );
    M_blockNumRows.push_back (vector.mapSize (0) );
    M_blockNumColumns.push_back (vector.mapSize (0) );
    M_blockFirstRows.push_back (0);
    M_blockFirstColumns.push_back (0);

    UInt totalRows (vector.mapSize (0) );
    UInt totalColumns (vector.mapSize (0) );

    for (UInt i (1); i < vector.nbMap(); ++i)
    {
        myMap += vector.map (i);
        M_blockNumRows.push_back (vector.mapSize (i) );
        M_blockNumColumns.push_back (vector.mapSize (i) );
        M_blockFirstRows.push_back (totalRows);
        M_blockFirstColumns.push_back (totalColumns);

        totalRows += vector.mapSize (i);
        totalColumns += vector.mapSize (i);
    }

    this->mapPtr().reset (new MapEpetra (myMap) );
    this->matrixPtr().reset ( new typename MatrixEpetra<DataType>::matrix_type ( Copy, *myMap.map ( Unique ), numEntries, false) );
}


template <typename DataType>
MatrixBlockMonolithicEpetra<DataType>::MatrixBlockMonolithicEpetra ( const MatrixEpetra<DataType>& matrix ) :
    MatrixEpetra<DataType> (matrix),
    M_blockNumRows(),
    M_blockNumColumns(),
    M_blockFirstRows(),
    M_blockFirstColumns()
{}

template <typename DataType>
MatrixBlockMonolithicEpetra<DataType>::MatrixBlockMonolithicEpetra ( const MatrixBlockMonolithicEpetra& matrix ) :
    MatrixEpetra<DataType> (matrix),
    M_blockNumRows (matrix.M_blockNumRows),
    M_blockNumColumns (matrix.M_blockNumColumns),
    M_blockFirstRows (matrix.M_blockFirstRows),
    M_blockFirstColumns (matrix.M_blockFirstColumns)
{}

template <typename DataType>
MatrixBlockMonolithicEpetra<DataType>::~MatrixBlockMonolithicEpetra()
{}

// ===================================================
// Set Methods
// ===================================================
template <typename DataType>
void
MatrixBlockMonolithicEpetra<DataType>::setBlockStructure (const std::vector<UInt>& blockNumRows,
                                                          const std::vector<UInt>& blockNumColumns)
{
    ASSERT ( blockNumRows.size() > 0, "No way to build a matrix with 0 block rows");
    ASSERT ( blockNumColumns.size() > 0, "No way to build a matrix with 0 block columns")


    M_blockNumRows    = blockNumRows;
    M_blockNumColumns = blockNumColumns;

    M_blockFirstRows.resize (M_blockNumRows.size() );
    M_blockFirstColumns.resize (M_blockNumColumns.size() );

    UInt currentSize (0);
    for (UInt i (0); i < blockNumRows.size(); ++i)
    {
        M_blockFirstRows[i] = currentSize;
        currentSize += M_blockNumRows[i];
    }

    currentSize = 0;
    for (UInt i (0); i < blockNumColumns.size(); ++i)
    {
        M_blockFirstColumns[i] = currentSize;
        currentSize += M_blockNumColumns[i];
    }
}

template <typename DataType>
void
MatrixBlockMonolithicEpetra<DataType>::setBlockStructure (const MapVector<MapEpetra>& mapVector)
{
    ASSERT ( mapVector.nbMap() > 0 , "Map vector empty, impossible to set the block structure");

    M_blockNumRows.resize (mapVector.nbMap() );
    M_blockNumColumns.resize (mapVector.nbMap() );

    M_blockFirstRows.resize (mapVector.nbMap() );
    M_blockFirstColumns.resize (mapVector.nbMap() );

    UInt totalSize (0);

    for (UInt i (0); i < mapVector.nbMap(); ++i)
    {
        M_blockNumRows[i] = mapVector.mapSize (i);
        M_blockNumColumns[i] = mapVector.mapSize (i);

        M_blockFirstRows[i] = totalSize;
        M_blockFirstColumns[i] = totalSize;

        totalSize += mapVector.mapSize (i);
    }

    ASSERT ( this->matrixPtr()->NumGlobalCols() == totalSize, " Incompatible block structure (global size does not match) ");
    ASSERT ( this->matrixPtr()->NumGlobalRows() == totalSize, " Incompatible block structure (global size does not match) ");
}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
UInt
MatrixBlockMonolithicEpetra<DataType>::blockNumRows (const UInt& rowIndex) const
{
    ASSERT (rowIndex < M_blockFirstRows.size(), "Row index out of bound. No block to return");

    return M_blockNumRows[rowIndex];
}

template <typename DataType>
UInt
MatrixBlockMonolithicEpetra<DataType>::blockNumColumns (const UInt& columnIndex) const
{
    ASSERT (columnIndex < M_blockFirstColumns.size(), "Column index out of bound. No block to return");

    return M_blockNumColumns[columnIndex];
}

template <typename DataType>
void
MatrixBlockMonolithicEpetra<DataType>::blockView (const UInt& rowIndex,
                                                  const UInt& columnIndex,
                                                  block_type& mbv)
{
    ASSERT (rowIndex < M_blockFirstRows.size(), "Row index out of bound. No block to return");
    ASSERT (columnIndex < M_blockFirstColumns.size(), "Column index out of bound. No block to return");

    mbv.setup (M_blockFirstRows[rowIndex],
               M_blockFirstColumns[columnIndex],
               M_blockNumRows[rowIndex],
               M_blockNumColumns[columnIndex],
               this);
}

template <typename DataType>
typename MatrixBlockMonolithicEpetra<DataType>::block_ptrType
MatrixBlockMonolithicEpetra<DataType>::block (const UInt& rowIndex, const UInt& columnIndex)
{
    ASSERT (rowIndex < M_blockFirstRows.size(), "Row index out of bound. No block to return");
    ASSERT (columnIndex < M_blockFirstColumns.size(), "Column index out of bound. No block to return");

    block_ptrType mbv (new block_type);

    mbv->setup (M_blockFirstRows[rowIndex],
                M_blockFirstColumns[columnIndex],
                M_blockNumRows[rowIndex],
                M_blockNumColumns[columnIndex],
                this);

    return mbv;
}


} // namespace LifeV

#endif /* MATRIXBLOCK_HPP */

