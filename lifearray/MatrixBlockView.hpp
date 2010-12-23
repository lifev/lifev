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
   @date 2010-10-09
 */

#ifndef _MATRIXBLOCKVIEW_HPP_
#define _MATRIXBLOCKVIEW_HPP_

#include <boost/shared_ptr.hpp>
#include <Epetra_FECrsMatrix.h>
#include <life/lifearray/MatrixEpetra.hpp>

namespace LifeV {

//! MatrixBlockView - class to manage the block in a block matrix
/*!
 *  @author Gwenol Grandperrin
 *
 *  The BlockMatrixView class contains data related
 *  to block of a matrix. It is useful to setup a
 *  clean and easy-to-use blocks management
 */
class MatrixBlockView
{
public:

    /** @name Typedefs
     */
    //@{
    typedef MatrixEpetra<double> matrix_type;
    typedef boost::shared_ptr<Epetra_FECrsMatrix>  matrix_ptrtype;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
    MatrixBlockView();

    //! Copy constructor
    MatrixBlockView( MatrixBlockView& mbv );

    //! default virtual destructor
    ~MatrixBlockView();

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
    matrix_ptrtype& getMatrixPtr(){return M_matrix;}

    //! Return the const shared_pointer of the Epetra_FECrsMatrix
    const matrix_ptrtype& getMatrixPtr() const{return M_matrix;}

    //@}

private:
    UInt M_numRows;
    UInt M_numColumns;
    UInt M_firstRowIndex;
    UInt M_lastRowIndex;
    UInt M_firstColumnIndex;
    UInt M_lastColumnIndex;
    matrix_ptrtype M_matrix;
};

} // namespace LifeV

#endif /* MATRIXBLOCKVIEW_HPP */

