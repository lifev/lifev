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
   @file MatrixBlockMonolithicEpetraView.hpp
   @brief The file contains the MatrixBlockMonolithicEpetraView class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIX_BLOCK_MONOLITHIC_EPETRA_VIEW_HPP_
#define _MATRIX_BLOCK_MONOLITHIC_EPETRA_VIEW_HPP_

#include <boost/shared_ptr.hpp>

#include <iostream>

#include <lifev/core/array/MatrixEpetra.hpp>


namespace LifeV
{

//! MatrixBlockMonolithicEpetraView - class representing a block in a MatrixBlockMonolithicEpetra
/*!
  @author Gwenol Grandperrin
  @author Samuel Quinodoz

  The MatrixBlockMonolithicEpetraView class contains data related
  to block of a matrix. It is useful to setup a clean and easy-to-use blocks management.

  For more information about the block structures in LifeV, see \ref BlockAlgebraPage "this page".

  <b> Remark </b>

  Using the operator "=" is not valid for this class! Indeed, copying the view
  would not copy the data stored in the matrix, which can be confusing. If the
  copy of the block is the intended use, one should use a method from the BlockUtils.

 */
template<typename DataType>
class MatrixBlockMonolithicEpetraView
{
public:

    /** @name Typedefs
     */
    //@{

    // Not a block matrix to avoid circular dependancies
    //typedef MatrixBlockMonolithicEpetra<DataType> matrix_Type;

    typedef MatrixEpetra<DataType> matrix_Type;

    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlockMonolithicEpetraView();

    //! Copy constructor
    MatrixBlockMonolithicEpetraView ( const MatrixBlockMonolithicEpetraView<DataType>& mbv );

    //! default virtual destructor
    ~MatrixBlockMonolithicEpetraView();

    //@}

    //! @name Methods
    //@{

    //! Print the informations about the MatrixBlockMonolithicEpetraView
    void showMe (std::ostream& output = std::cout) const;

    //! Tells if the viewed matrix is already filled
    bool filled() const;

    //! Function to assemble an elemental matrix in a block
    void addToCoefficients ( UInt const numRows, UInt const numColumns,
                             std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                             DataType* const* const localValues,
                             Int format = Epetra_FECrsMatrix::COLUMN_MAJOR ) const;

    //! Function to sum an elemental matrix in a block which is already closed
    void sumIntoCoefficients ( UInt const numRows, UInt const numColumns,
                               std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                               DataType* const* const localValues,
                               Int format = Epetra_FECrsMatrix::COLUMN_MAJOR ) const;

    //@}

    //! @name  Set Methods
    //@{
    /*! Set all the informations relative to the block
     *  @param firstRow First row in the block
     *  @param firstColumn First column in the block
     *  @param numRows Number of rows in the block
     *  @param numColumns Number of columns in the block
     *  @param A Matrix which contains the block
     */
    void setup ( const UInt& firstRow,
                 const UInt& firstColumn,
                 const UInt& numRows,
                 const UInt& numColumns,
                 matrix_Type* A );

    //@}

    //! @name  Get Methods
    //@{
    //! Returns the number of rows in the block
    UInt numRows() const
    {
        return M_numRows;
    }

    //! Returns the number of columns in the block
    UInt numColumns() const
    {
        return M_numColumns;
    }

    //! Returns the index of the first row in the block
    UInt firstRowIndex() const
    {
        return M_firstRowIndex;
    }

    //! Returns the index of the last row in the block
    UInt lastRowIndex() const
    {
        return M_lastRowIndex;
    }

    //! Returns the index of the first column in the block
    UInt firstColumnIndex() const
    {
        return M_firstColumnIndex;
    }

    //! Returns the index of the last column in the block
    UInt lastColumnIndex() const
    {
        return M_lastColumnIndex;
    }

    //! Return the pointer of the full matrix
    matrix_Type* matrixPtr() const
    {
        return M_matrix;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No assignement operator, it is missleading (would copy the views, not the blocks!)
    MatrixBlockMonolithicEpetraView<DataType> operator= ( const MatrixBlockMonolithicEpetraView& otherView);

    //@}


    UInt M_numRows;
    UInt M_numColumns;
    UInt M_firstRowIndex;
    UInt M_lastRowIndex;
    UInt M_firstColumnIndex;
    UInt M_lastColumnIndex;
    matrix_Type* M_matrix;
};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename DataType>
MatrixBlockMonolithicEpetraView<DataType>::MatrixBlockMonolithicEpetraView() :
    M_numRows ( 0 ),
    M_numColumns ( 0 ),
    M_firstRowIndex ( 0 ),
    M_lastRowIndex ( 0 ),
    M_firstColumnIndex ( 0 ),
    M_lastColumnIndex ( 0 ),
    M_matrix()
{

}

template<typename DataType>
MatrixBlockMonolithicEpetraView<DataType>::MatrixBlockMonolithicEpetraView ( const MatrixBlockMonolithicEpetraView<DataType>& mbv ) :
    M_numRows ( mbv.M_numRows ),
    M_numColumns ( mbv.M_numColumns ),
    M_firstRowIndex ( mbv.M_firstRowIndex ),
    M_lastRowIndex ( mbv.M_lastRowIndex ),
    M_firstColumnIndex ( mbv.M_firstColumnIndex ),
    M_lastColumnIndex ( mbv.M_lastColumnIndex ),
    M_matrix ( mbv.M_matrix )
{

}

template<typename DataType>
MatrixBlockMonolithicEpetraView<DataType>::~MatrixBlockMonolithicEpetraView()
{
    //M_matrix.reset();
}

// ===================================================
// Methods
// ===================================================

template<typename DataType>
void
MatrixBlockMonolithicEpetraView<DataType>::showMe ( std::ostream& output ) const
{
    output << "MatrixBlockMonolithicEpetraView informations:" << std::endl
           << "Size = " << M_numRows << " x " << M_numColumns << std::endl
           << "firstRow = " << M_firstRowIndex << std::endl
           << "lastRow = " << M_lastRowIndex << std::endl
           << "firstColumn = " << M_firstColumnIndex << std::endl
           << "lastColumn = " << M_lastColumnIndex << std::endl;
}


template<typename DataType>
bool
MatrixBlockMonolithicEpetraView<DataType>::
filled() const
{
    return M_matrix->filled();
}


template<typename DataType>
void
MatrixBlockMonolithicEpetraView<DataType>::
addToCoefficients ( UInt const numRows, UInt const numColumns,
                    std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                    DataType* const* const localValues,
                    Int format) const
{
    std::vector<Int> rowIndices (blockRowIndices);
    std::vector<Int> columnIndices (blockColumnIndices);

    for (UInt i (0); i < numRows; ++i)
    {
        rowIndices[i] += M_firstRowIndex;
    }
    for (UInt i (0); i < numColumns; ++i)
    {
        columnIndices[i] += M_firstColumnIndex;
    }

    M_matrix->addToCoefficients (numRows, numColumns,
                                 rowIndices, columnIndices,
                                 localValues, format);
}


template<typename DataType>
void
MatrixBlockMonolithicEpetraView<DataType>::
sumIntoCoefficients ( UInt const numRows, UInt const numColumns,
                      std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                      DataType* const* const localValues,
                      Int format) const
{
    std::vector<Int> rowIndices (blockRowIndices);
    std::vector<Int> columnIndices (blockColumnIndices);

    for (UInt i (0); i < numRows; ++i)
    {
        rowIndices[i] += M_firstRowIndex;
    }
    for (UInt i (0); i < numColumns; ++i)
    {
        columnIndices[i] += M_firstColumnIndex;
    }

    M_matrix->sumIntoCoefficients (numRows, numColumns,
                                   rowIndices, columnIndices,
                                   localValues, format);
}





// ===================================================
// Set Methods
// ===================================================

template<typename DataType>
void
MatrixBlockMonolithicEpetraView<DataType>::setup ( const UInt& firstRow,
                                                   const UInt& firstColumn,
                                                   const UInt& numRows,
                                                   const UInt& numColumns,
                                                   matrix_Type* A )
{
    M_numRows          = numRows;
    M_numColumns       = numColumns;
    M_firstRowIndex    = firstRow;
    M_lastRowIndex     = firstRow + numRows - 1;
    M_firstColumnIndex = firstColumn;
    M_lastColumnIndex  = firstColumn + numColumns - 1;
    M_matrix           = A;
}

// ===================================================
// Get Methods
// ===================================================

} // namespace LifeV

#endif /* MATRIXBLOCKVIEW_HPP */
