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
   @file MatrixEpetraStructuredView.hpp
   @brief The file contains the MatrixEpetraStructuredView class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIXEPETRASTRUCTUREDVIEW_HPP_
#define _MATRIXEPETRASTRUCTUREDVIEW_HPP_


#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>


namespace LifeV
{

//! MatrixEpetraStructuredView - class representing a block in a MatrixEpetraStructured
/*!
  @author Gwenol Grandperrin
  @author Samuel Quinodoz

  The MatrixEpetraStructuredView class contains data related
  to block of a matrix. It is useful to setup a clean and easy-to-use blocks management.

  For more information about the block structures in LifeV, see \ref BlockAlgebraPage "this page".

  <b> Remark </b>

  Using the operator "=" is not valid for this class! Indeed, copying the view
  would not copy the data stored in the matrix, which can be confusing. If the
  copy of the block is the intended use, one should use a method from the BlockUtils.

 */
template<typename DataType>
class MatrixEpetraStructuredView
{
public:

    /** @name Typedefs
     */
    //@{

    // Not a block matrix to avoid circular dependancies
    //typedef MatrixEpetraStructured<DataType> matrix_Type;

    typedef MatrixEpetra<DataType> matrix_Type;

    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixEpetraStructuredView();

    //! Copy constructor
    MatrixEpetraStructuredView ( const MatrixEpetraStructuredView<DataType>& matrixEpetraStructured );

    //! default virtual destructor
    ~MatrixEpetraStructuredView();

    //@}

    //! @name Methods
    //@{

    //! Print the informations about the MatrixEpetraStructuredView
    void showMe (std::ostream& output = std::cout) const;

    //! Function to assemble an elemental matrix in a block
    void addToCoefficients ( UInt const numRows, UInt const numColumns,
                             std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                             DataType* const* const localValues,
                             Int format = Epetra_FECrsMatrix::COLUMN_MAJOR ) const;
    //! Function to assemble an elemental matrix in a block (for closed matrices)
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

    //! Return the fill-complete status of the inner Epetra_FECrsMatrix
    bool filled() const
    {
        return M_matrix->matrixPtr()->Filled();
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
    MatrixEpetraStructuredView<DataType> operator= ( const MatrixEpetraStructuredView& otherView);

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
MatrixEpetraStructuredView<DataType>::MatrixEpetraStructuredView() :
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
MatrixEpetraStructuredView<DataType>::MatrixEpetraStructuredView ( const MatrixEpetraStructuredView<DataType>& matrixEpetraStructured ) :
    M_numRows ( matrixEpetraStructured.M_numRows ),
    M_numColumns ( matrixEpetraStructured.M_numColumns ),
    M_firstRowIndex ( matrixEpetraStructured.M_firstRowIndex ),
    M_lastRowIndex ( matrixEpetraStructured.M_lastRowIndex ),
    M_firstColumnIndex ( matrixEpetraStructured.M_firstColumnIndex ),
    M_lastColumnIndex ( matrixEpetraStructured.M_lastColumnIndex ),
    M_matrix ( matrixEpetraStructured.M_matrix )
{

}

template<typename DataType>
MatrixEpetraStructuredView<DataType>::~MatrixEpetraStructuredView()
{
    //M_matrix.reset();
}

// ===================================================
// Methods
// ===================================================

template<typename DataType>
void
MatrixEpetraStructuredView<DataType>::showMe ( std::ostream& output ) const
{
    output << "MatrixBlockMonolithicEpetraView informations:" << std::endl
           << "Size = " << M_numRows << " x " << M_numColumns << std::endl
           << "firstRow = " << M_firstRowIndex << std::endl
           << "lastRow = " << M_lastRowIndex << std::endl
           << "firstColumn = " << M_firstColumnIndex << std::endl
           << "lastColumn = " << M_lastColumnIndex << std::endl;
}

template<typename DataType>
void
MatrixEpetraStructuredView<DataType>::
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
MatrixEpetraStructuredView<DataType>::
sumIntoCoefficients ( UInt const numRows, UInt const numColumns,
                      std::vector<Int> const& blockRowIndices,
                      std::vector<Int> const& blockColumnIndices,
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
MatrixEpetraStructuredView<DataType>::setup ( const UInt& firstRow,
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

#endif /* _MATRIXEPETRASTRUCTUREDVIEW_HPP_ */

