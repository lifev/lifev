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
   @file MatrixEpetraStructured.hpp
   @brief The file contains the MatrixEpetraStructured class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIXEPETRASTRUCTURED_HPP_
#define _MATRIXEPETRASTRUCTURED_HPP_

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixBlockStructure.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>

namespace LifeV
{

//! MatrixEpetraStructured - class of block matrix
/*!
  @author Gwenol Grandperrin
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  The MatrixEpetraStructured class contains data related
  to block matrix. It is an extension to MatrixEpetra where data about blocks have
  been set. For an introduction to the block structures in LifeV, see
  \ref BlockAlgebraPage "this page".

  There are mainly two ways to define a MatrixEpetraStructured:
  <ul>
  <li> Construct it using the same syntax as for LifeV::MatrixEpetra and the use
  a setter for the structure.
  <li> Construct it directly with the maps of the blocks.
  </ul>
  Both ways are equivalent.

  To access the blocks, one uses then the blockView or block methods.

 */
template <typename DataType>
class MatrixEpetraStructured : public MatrixEpetra<DataType>
{
public:

    typedef MatrixEpetraStructuredView<DataType> block_type;
    typedef boost::shared_ptr<block_type> block_ptrType;

    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixEpetraStructured ( const MapEpetra& map, int numEntries = 50, bool ignoreNonLocalValues = false );

    //! Block constructor
    /*!
    This is the most complete constructor, as it builds the whole block
    structure of the matrix, using the maps stored in the vector.
    */
    MatrixEpetraStructured ( const MapVector<MapEpetra>& vector, int numEntries = 50, bool ignoreNonLocalValues = false );

    //! Casting constructor
    MatrixEpetraStructured ( const MatrixEpetra<DataType>& matrix );

    //! Copy constructor
    MatrixEpetraStructured ( const MatrixEpetraStructured& matrix );

    //! Destructor
    ~MatrixEpetraStructured();

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
    void setBlockStructure ( const MapVector<MapEpetra>& mapVector );

    //! Set the block structure from a matrix structure object
    void setBlockStructure ( const MatrixBlockStructure& blockStructure );

    //! Set the block structure from row and column structures
    void setBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                             const VectorBlockStructure& columnsBlockStructure );

    //! Set the block structure from row and column structures
    void setBlockStructure ( const VectorBlockStructure& vectorBlockStructure );

    //@}


    //! @name  Get Methods
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

    //! Returns the MatrixEpetraStructuredView at the (rowIndex,colIndex) location
    /*!
     * @param rowIndex Row position of the block in the matrix
     * @param columnIndex Column position of the block in the matrix
     * @param mbv MatrixEpetraStructuredView to be filled
     */
    void blockView ( const UInt& rowIndex,
                     const UInt& columnIndex,
                     block_type& mbv );

    //! Returns the block (rowIndex,columnIndex) of the matrix
    /*!
      Remark that the returned block is a shared pointer. This is
      a limitation due to the non-copiability of the blocks.
      @param rowIndex The index of the block w.r. to the block structure
      @param columnIndex The index of the block w.r. to the block structure
     */
    block_ptrType block ( const UInt& rowIndex, const UInt& columnIndex );

    //! Get the rows block structure
    VectorBlockStructure rowsBlockStructure() const;

    //! Get the columns block structure
    VectorBlockStructure columnsBlockStructure() const;

    //! Get the matrix block structure
    MatrixBlockStructure blockStructure() const;

    //@}


    //! @name  Linear Algebra Methods
    //@{

    //! Returns a pointer to a new matrix which contains the transpose of the current matrix
    boost::shared_ptr<MatrixEpetraStructured<DataType> > transpose()
    {
        return boost::static_pointer_cast<MatrixEpetraStructured<DataType> > ( this->MatrixEpetra<DataType>::transpose() );
    }

    //@}

private:

    MatrixBlockStructure M_blockStructure;
};

// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================
template <typename DataType>
MatrixEpetraStructured<DataType>::MatrixEpetraStructured ( const MapEpetra& map, int numEntries, bool ignoreNonLocalValues ) :
    MatrixEpetra<DataType> ( map, numEntries, ignoreNonLocalValues ),
    M_blockStructure ( map )
{

}



template <typename DataType>
MatrixEpetraStructured<DataType>::MatrixEpetraStructured ( const MapVector<MapEpetra>& vector, int numEntries, bool ignoreNonLocalValues ) :
    MatrixEpetra<DataType> ( typename MatrixEpetra<DataType>::matrix_ptrtype() ), M_blockStructure ( vector )
{
    ASSERT ( vector.nbMap() > 0 , "Map vector empty, impossible to construct a MatrixBlockMonolithicEpetra!" );

    MapEpetra myMap ( vector.totalMap() );

    this->mapPtr().reset ( new MapEpetra ( myMap ) );
    this->matrixPtr().reset ( new typename MatrixEpetra<DataType>::matrix_type ( Copy, *myMap.map ( Unique ), numEntries, ignoreNonLocalValues ) );
}


template <typename DataType>
MatrixEpetraStructured<DataType>::MatrixEpetraStructured ( const MatrixEpetra<DataType>& matrix ) :
    MatrixEpetra<DataType> ( matrix ),
    M_blockStructure ( matrix.map() )
{}

template <typename DataType>
MatrixEpetraStructured<DataType>::MatrixEpetraStructured ( const MatrixEpetraStructured& matrix ) :
    MatrixEpetra<DataType> ( matrix ),
    M_blockStructure ( matrix.M_blockStructure )
{}

template <typename DataType>
MatrixEpetraStructured<DataType>::~MatrixEpetraStructured()
{}

// ===================================================
// Set Methods
// ===================================================
template <typename DataType>
void
MatrixEpetraStructured<DataType>::setBlockStructure ( const std::vector<UInt>& blockNumRows,
                                                      const std::vector<UInt>& blockNumColumns )
{
    ASSERT ( blockNumRows.size() > 0, "No way to build a matrix with 0 block rows" );
    ASSERT ( blockNumColumns.size() > 0, "No way to build a matrix with 0 block columns" );


    M_blockStructure.setBlockStructure ( blockNumRows, blockNumColumns );
}

template <typename DataType>
void
MatrixEpetraStructured<DataType>::setBlockStructure ( const MapVector<MapEpetra>& mapVector )
{
    ASSERT ( mapVector.nbMap() > 0 , "Map vector empty, impossible to set the block structure" );

    M_blockStructure.setBlockStructure ( mapVector );

    ASSERT ( this->matrixPtr()->NumGlobalCols() == M_blockStructure.numRows(), " Incompatible block structure (global size does not match) " );
    ASSERT ( this->matrixPtr()->NumGlobalRows() == M_blockStructure.numColumns(), " Incompatible block structure (global size does not match) " );
}

template <typename DataType>
void
MatrixEpetraStructured<DataType>::setBlockStructure ( const MatrixBlockStructure& blockStructure )
{
    M_blockStructure.setBlockStructure ( blockStructure );
}

template <typename DataType>
void
MatrixEpetraStructured<DataType>::setBlockStructure ( const VectorBlockStructure& rowsBlockStructure,
                                                      const VectorBlockStructure& columnsBlockStructure )
{
    M_blockStructure.setBlockStructure ( rowsBlockStructure, columnsBlockStructure );
}

template <typename DataType>
void
MatrixEpetraStructured<DataType>::setBlockStructure ( const VectorBlockStructure& vectorBlockStructure )
{
    M_blockStructure.setBlockStructure ( vectorBlockStructure, vectorBlockStructure );
}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
UInt
MatrixEpetraStructured<DataType>::blockNumRows ( const UInt& rowIndex ) const
{
    return M_blockStructure.blockNumRows ( rowIndex );
}

template <typename DataType>
UInt
MatrixEpetraStructured<DataType>::blockNumColumns ( const UInt& columnIndex ) const
{
    return M_blockStructure.blockNumColumns ( columnIndex );
}

template <typename DataType>
void
MatrixEpetraStructured<DataType>::blockView ( const UInt& rowIndex,
                                              const UInt& columnIndex,
                                              block_type& mbv )
{
    mbv.setup ( M_blockStructure.rowBlockFirstIndex ( rowIndex ),
                M_blockStructure.columnBlockFirstIndex ( columnIndex ),
                M_blockStructure.blockNumRows ( rowIndex ),
                M_blockStructure.blockNumColumns ( columnIndex ),
                this );
}

template <typename DataType>
typename MatrixEpetraStructured<DataType>::block_ptrType
MatrixEpetraStructured<DataType>::block ( const UInt& rowIndex, const UInt& columnIndex )
{
    block_ptrType matrixBlockView ( new block_type );

    matrixBlockView->setup ( M_blockStructure.rowBlockFirstIndex ( rowIndex ),
                             M_blockStructure.columnBlockFirstIndex ( columnIndex ),
                             M_blockStructure.blockNumRows ( rowIndex ),
                             M_blockStructure.blockNumColumns ( columnIndex ),
                             this );

    return matrixBlockView;
}

template <typename DataType>
VectorBlockStructure
MatrixEpetraStructured<DataType>::rowsBlockStructure() const
{
    return M_blockStructure.rowsBlockStructure();
}

template <typename DataType>
VectorBlockStructure
MatrixEpetraStructured<DataType>::columnsBlockStructure() const
{
    return M_blockStructure.columnsBlockStructure();
}

template <typename DataType>
MatrixBlockStructure
MatrixEpetraStructured<DataType>::blockStructure() const
{
    return M_blockStructure;
}


} // namespace LifeV

#endif /* _MATRIXEPETRASTRUCTURED_HPP_ */

