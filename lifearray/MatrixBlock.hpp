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
   @file MatrixBlock.hpp
   @brief The file contains the MatrixBlock class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIXBLOCK_HPP_
#define _MATRIXBLOCK_HPP_

#include <boost/shared_ptr.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <lifemc/lifearray/MatrixBlockView.hpp>

namespace LifeV {

//! MatrixBlock - class of block matrix
/*!
 *  @author Gwenol Grandperrin
 *
 *  The BlockMatrix class contains data related
 *  to block matrix. It is an extension to
 *  EpetraMatrix where data about blocks have
 *  been set.
 */
template <typename DataType>
class MatrixBlock : public MatrixEpetra<DataType>
{
public:

    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlock( const MapEpetra& map, int numEntries = 50, int indexBase = 1 );

    //! Casting constructor
    MatrixBlock( const MatrixEpetra<DataType>& matrix );

    //! Copy constructor
    MatrixBlock( const MatrixBlock& matrix );

    //! default virtual destructor
    ~MatrixBlock();

    //@}

    //! @name  Set Methods
    //@{
    /*! Set the size of the blocks of the matrix
     *  @param blockNumRows Number of rows in the blocks
     *  @param blockNumColumns Number of columns in the blocks
     */
    void setBlockStructure( const std::vector<UInt>& blockNumRows,
                            const std::vector<UInt>& blockNumColumns );

    //@}

    //! @name  Get Methods
    //@{
    //! Returns the number of rows of the block
    /*!
     * @param rowIndex Row index of the block
     */
    UInt getBlockNumRows(const UInt& rowIndex) const;

    //! Returns the number of columns of the block
    /*!
     * @param columnIndex Column index of the block
     */
    UInt getBlockNumColumns(const UInt& columnIndex) const;

    //! Returns the MatrixBlockView at the (rowIndex,colIndex) location
    /*!
     * @param rowIndex Row position of the block in the matrix
     * @param columnIndex Column position of the block in the matrix
     * @param mbv MatrixBlockView to be filled
     */
    void getMatrixBlockView( const UInt& rowIndex,
                             const UInt& columnIndex,
                             MatrixBlockView& mbv);

    //@}

private:

    std::vector<UInt> M_blockNumRows;
    std::vector<UInt> M_blockNumColumns;

};

// ===================================================
// IMPLEMENTATION
// ===================================================

// ===================================================
// Constructors & Destructor
// ===================================================
template <typename DataType>
MatrixBlock<DataType>::MatrixBlock(const MapEpetra& map, int numEntries, int indexBase) :
    MatrixEpetra<DataType>(map,numEntries,indexBase)
{

}

template <typename DataType>
MatrixBlock<DataType>::MatrixBlock( const MatrixEpetra<DataType>& matrix ) :
    MatrixEpetra<DataType>(matrix)
{

}

template <typename DataType>
MatrixBlock<DataType>::MatrixBlock( const MatrixBlock& matrix ) :
    MatrixEpetra<DataType>(matrix),
    M_blockNumRows(matrix.M_blockNumRows),
    M_blockNumColumns(matrix.M_blockNumColumns)
{

}

template <typename DataType>
MatrixBlock<DataType>::~MatrixBlock()
{

}

// ===================================================
// Set Methods
// ===================================================
template <typename DataType>
void
MatrixBlock<DataType>::setBlockStructure(const std::vector<UInt>& blockNumRows,
                                         const std::vector<UInt>& blockNumColumns)
{
    M_blockNumRows    = blockNumRows;
    M_blockNumColumns = blockNumColumns;
}

// ===================================================
// Get Methods
// ===================================================
template <typename DataType>
UInt
MatrixBlock<DataType>::getBlockNumRows(const UInt& rowIndex) const
{
    return M_blockNumRows[rowIndex];
}

template <typename DataType>
UInt
MatrixBlock<DataType>::getBlockNumColumns(const UInt& columnIndex) const
{
    return M_blockNumColumns[columnIndex];
}

template <typename DataType>
void
MatrixBlock<DataType>::getMatrixBlockView(const UInt& rowIndex,
                                          const UInt& columnIndex,
                                          MatrixBlockView& mbv)
{
    UInt firstRow(0);
    UInt firstCol(0);
    for(UInt i(0);i<rowIndex;++i)
    {
        firstRow += M_blockNumRows[i];
    }
    for(UInt i(0);i<columnIndex;++i)
    {
        firstCol += M_blockNumColumns[i];
    }

    mbv.setup(firstRow,
              firstCol,
              M_blockNumRows[rowIndex],
              M_blockNumColumns[columnIndex],
              *this);
}

} // namespace LifeV

#endif /* MATRIXBLOCK_HPP */

