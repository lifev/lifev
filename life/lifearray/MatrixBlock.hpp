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
   @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIXBLOCK_HPP_
#define _MATRIXBLOCK_HPP_

#include <boost/shared_ptr.hpp>
#include <life/lifearray/MapEpetra.hpp>
#include <life/lifearray/MapVector.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/MatrixBlockView.hpp>

namespace LifeV {

//! MatrixBlock - class of block matrix
/*!
 *  @author Gwenol Grandperrin
 *  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
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

    typedef MatrixBlockView<DataType> block_type;
    typedef boost::shared_ptr<block_type> block_ptrType;

    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlock( const MapEpetra& map, int numEntries = 50);

	//! Block constructor
	/*!
	This is the most complete constructor, as it builds the whole block
	structure of the matrix, using the maps stored in the vector.
	*/
	MatrixBlock( const MapEpetraVector& vector, int numEntries = 50);

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
                             block_type& mbv);

	//! Returns the block (rowIndex,columnIndex) of the matrix
	/*!
	@assert: rowIndex is a valid row number
	@assert: columnIndex is a valid column number
	*/
	block_ptrType block(const UInt& rowIndex, const UInt& columnIndex);
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
MatrixBlock<DataType>::MatrixBlock(const MapEpetra& map, int numEntries) :
    MatrixEpetra<DataType>(map,numEntries),
	M_blockNumRows(std::vector<UInt>(1,map.map(Unique)->NumMyElements())),
	M_blockNumColumns(std::vector<UInt>(1,map.map(Unique)->NumMyElements())),
	M_blockFirstRows(std::vector<UInt>(1,0)),
	M_blockFirstColumns(std::vector<UInt>(1,0))
{
}



template <typename DataType>
MatrixBlock<DataType>::MatrixBlock( const MapEpetraVector& vector, int numEntries):
    MatrixEpetra<DataType>( typename MatrixEpetra<DataType>::matrix_ptrtype())
{
	ASSERT( vector.nbMap() > 0 ,"Map vector empty, impossible to construct a MatrixBlock!");

	MapEpetra myMap(vector.map(0));
	M_blockNumRows.push_back(vector.map(0).map(Unique)->NumMyElements());
	M_blockNumColumns.push_back(vector.map(0).map(Unique)->NumMyElements());
	M_blockFirstRows.push_back(0);
	M_blockFirstColumns.push_back(0);

	UInt totalRows(vector.map(0).map(Unique)->NumMyElements());
	UInt totalColumns(vector.map(0).map(Unique)->NumMyElements());

	for (UInt i(1); i<vector.nbMap(); ++i)
	{
		myMap += vector.map(i);
		M_blockNumRows.push_back(vector.map(i).map(Unique)->NumMyElements());
		M_blockNumColumns.push_back(vector.map(i).map(Unique)->NumMyElements());
		M_blockFirstRows.push_back(totalRows);
		M_blockFirstColumns.push_back(totalColumns);

		totalRows+= vector.map(i).map(Unique)->NumMyElements();
		totalColumns+= vector.map(i).map(Unique)->NumMyElements();
	}

	this->mapPtr().reset(new MapEpetra(myMap));
	this->matrixPtr().reset( new typename MatrixEpetra<DataType>::matrix_type( Copy, *myMap.map( Unique ), numEntries, false));
}


template <typename DataType>
MatrixBlock<DataType>::MatrixBlock( const MatrixEpetra<DataType>& matrix ) :
    MatrixEpetra<DataType>(matrix),
	M_blockNumRows(),
	M_blockNumColumns(),
	M_blockFirstRows(),
	M_blockFirstColumns()
{}

template <typename DataType>
MatrixBlock<DataType>::MatrixBlock( const MatrixBlock& matrix ) :
    MatrixEpetra<DataType>(matrix),
    M_blockNumRows(matrix.M_blockNumRows),
    M_blockNumColumns(matrix.M_blockNumColumns),
	M_blockFirstRows(matrix.M_blockFirstRows),
	M_blockFirstColumns(matrix.M_blockFirstColumns)
{}

template <typename DataType>
MatrixBlock<DataType>::~MatrixBlock()
{}

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
                                          block_type& mbv)
{
    mbv.setup(M_blockFirstRows[rowIndex],
              M_blockFirstColumns[columnIndex],
              M_blockNumRows[rowIndex],
              M_blockNumColumns[columnIndex],
              *this);
}

template <typename DataType>
typename MatrixBlock<DataType>::block_ptrType
MatrixBlock<DataType>::block(const UInt& rowIndex, const UInt& columnIndex)
{
	ASSERT(rowIndex < M_blockFirstRows.size(), "Row index out of bound. No block to return");
	ASSERT(columnIndex < M_blockFirstColumns.size(), "Column index out of bound. No block to return");

	block_ptrType mbv(new block_type);

    mbv->setup(M_blockFirstRows[rowIndex],
              M_blockFirstColumns[columnIndex],
              M_blockNumRows[rowIndex],
              M_blockNumColumns[columnIndex],
              *this);

	return mbv;
}


} // namespace LifeV

#endif /* MATRIXBLOCK_HPP */

