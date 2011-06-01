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
   @file MatrixBlockView.hpp
   @brief The file contains the MatrixBlockView class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
   @date 2010-10-09
 */

#ifndef _MATRIXBLOCKVIEW_HPP_
#define _MATRIXBLOCKVIEW_HPP_

#include <boost/shared_ptr.hpp>

#include <iostream>

#include <Epetra_FECrsMatrix.h>

#include <life/lifearray/MatrixEpetra.hpp>

namespace LifeV {

//! MatrixBlockView - class to manage the block in a block matrix
/*!
 *  @author Gwenol Grandperrin
 *  @contributor Samuel Quinodoz
 *
 *  The BlockMatrixView class contains data related
 *  to block of a matrix. It is useful to setup a
 *  clean and easy-to-use blocks management
 */
template<typename DataType>
class MatrixBlockView
{
public:

    /** @name Typedefs
     */
    //@{
    typedef MatrixEpetra<DataType> matrix_type;
    typedef boost::shared_ptr<matrix_type>  matrix_ptrtype;
    typedef typename matrix_type::matrix_type rawMatrix_type;
    typedef boost::shared_ptr<rawMatrix_type>  rawMatrix_ptrtype;
	//@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlockView();

    //! Copy constructor
    MatrixBlockView( const MatrixBlockView<DataType>& mbv );

    //! default virtual destructor
    ~MatrixBlockView();

    //@}

    //! @name Methods
    //@{

    //! Print the informations about the MatrixBlockView
    void showMe(std::ostream& output = std::cout) const;

	//! Function to assemble an elemental matrix in a block
    void addToCoefficients( UInt const numRows, UInt const numColumns,
                            std::vector<Int> const& blockRowIndices, std::vector<Int> const& blockColumnIndices,
                            DataType* const* const localValues,
                            Int format = Epetra_FECrsMatrix::COLUMN_MAJOR )
	{
		std::vector<Int> rowIndices(blockRowIndices);
		std::vector<Int> columnIndices(blockColumnIndices);

		for (UInt i(0); i<numRows; ++i)
		{
			rowIndices[i]+=M_firstRowIndex;
		}
		for (UInt i(0); i<numColumns; ++i)
		{
			columnIndices[i]+=M_firstColumnIndex;
		}

		Int ierr = M_matrix->InsertGlobalValues( numRows, &rowIndices[0], numColumns,
                                                &columnIndices[0], localValues, format );

		ASSERT( ierr != -2, " \n <!> Error in block matrix insertion <!> \n Code : -2 \n Possible cause : try to insert a new element in a closed matrix")
		ASSERT( ierr >= 0 , " \n <!> Unknown error in block matrix insertion <!> ");
	}

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
    void setup( const UInt& firstRow,
                const UInt& firstColumn,
                const UInt& numRows,
                const UInt& numColumns,
                const matrix_type& A );

    //@}

    //! @name  Get Methods
    //@{
    //! Returns the number of rows in the block
    UInt numRows() const;

    //! Returns the number of columns in the block
    UInt numColumns() const;

    //! Returns the index of the first row in the block
    UInt firstRowIndex() const;

    //! Returns the index of the last row in the block
    UInt lastRowIndex() const;

    //! Returns the index of the first column in the block
    UInt firstColumnIndex() const;

    //! Returns the index of the last column in the block
    UInt lastColumnIndex() const;

    //! Return the shared_pointer of the Epetra_FECrsMatrix
    rawMatrix_ptrtype& getMatrixPtr(){return M_matrix;}

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const rawMatrix_ptrtype& getMatrixPtr() const{return M_matrix;}

    //@}

private:
    UInt M_numRows;
    UInt M_numColumns;
    UInt M_firstRowIndex;
    UInt M_lastRowIndex;
    UInt M_firstColumnIndex;
    UInt M_lastColumnIndex;
    rawMatrix_ptrtype M_matrix;
};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename DataType>
MatrixBlockView<DataType>::MatrixBlockView() :
    M_numRows( 0 ),
    M_numColumns( 0 ),
    M_firstRowIndex( 0 ),
    M_lastRowIndex( 0 ),
    M_firstColumnIndex( 0 ),
    M_lastColumnIndex( 0 ),
    M_matrix()
{

}

template<typename DataType>
MatrixBlockView<DataType>::MatrixBlockView( const MatrixBlockView<DataType>& mbv ) :
    M_numRows( mbv.M_numRows ),
    M_numColumns( mbv.M_numColumns ),
    M_firstRowIndex( mbv.M_firstRowIndex ),
    M_lastRowIndex( mbv.M_lastRowIndex ),
    M_firstColumnIndex( mbv.M_firstColumnIndex ),
    M_lastColumnIndex( mbv.M_lastColumnIndex ),
    M_matrix( mbv.M_matrix )
{

}

template<typename DataType>
MatrixBlockView<DataType>::~MatrixBlockView()
{
    M_matrix.reset();
}

// ===================================================
// Methods
// ===================================================

template<typename DataType>
void
MatrixBlockView<DataType>::showMe( std::ostream& output ) const
{
    output << "MatrixBlockView informations:" << std::endl
           << "Size = " << M_numRows << " x " << M_numColumns << std::endl
           << "firstRow = " << M_firstRowIndex << std::endl
           << "lastRow = " << M_lastRowIndex << std::endl
           << "firstColumn = " << M_firstColumnIndex << std::endl
           << "lastColumn = " << M_lastColumnIndex << std::endl;
}

// ===================================================
// Set Methods
// ===================================================

template<typename DataType>
void
MatrixBlockView<DataType>::setup( const UInt& firstRow,
                        const UInt& firstColumn,
                        const UInt& numRows,
                        const UInt& numColumns,
                        const matrix_type& A )
{
    M_numRows          = numRows;
    M_numColumns       = numColumns;
    M_firstRowIndex    = firstRow;
    M_lastRowIndex     = firstRow+numRows-1;
    M_firstColumnIndex = firstColumn;
    M_lastColumnIndex  = firstColumn+numColumns-1;
    M_matrix           = A.matrixPtr();
}

// ===================================================
// Get Methods
// ===================================================

template<typename DataType>
UInt
MatrixBlockView<DataType>::numRows() const
{
    return M_numRows;
}

template<typename DataType>
UInt
MatrixBlockView<DataType>::numColumns() const
{
    return M_numColumns;
}

template<typename DataType>
UInt
MatrixBlockView<DataType>::firstRowIndex() const
{
    return M_firstRowIndex;
}

template<typename DataType>
UInt
MatrixBlockView<DataType>::lastRowIndex() const
{
    return M_lastRowIndex;
}

template<typename DataType>
UInt
MatrixBlockView<DataType>::firstColumnIndex() const
{
    return M_firstColumnIndex;
}

template<typename DataType>
UInt
MatrixBlockView<DataType>::lastColumnIndex() const
{
    return M_lastColumnIndex;
}


} // namespace LifeV

#endif /* MATRIXBLOCKVIEW_HPP */

