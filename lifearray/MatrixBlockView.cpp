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
   @file MatrixBlockView.cpp
   @brief The file contains the MatrixBlockView class

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-10-09
 */

#include <MatrixBlockView.hpp>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
MatrixBlockView::MatrixBlockView() :
    M_numRows( 0 ),
    M_numCols( 0 ),
    M_firstRowIndex( 0 ),
    M_lastRowIndex( 0 ),
    M_firstColumnIndex( 0 ),
    M_lastColumnIndex( 0 ),
    M_matrix()
{

}

MatrixBlockView::MatrixBlockView( const MatrixBlockView& mbv ) :
    M_numRows( mbv.M_numRows ),
    M_numCols( mbv.M_numCols ),
    M_firstRowIndex( mbv.M_firstRowIndex ),
    M_lastRowIndex( mbv.M_firstRowIndex ),
    M_firstColumnIndex( mbv.M_firstColumnIndex ),
    M_lastColumnIndex( mbv.M_firstColumnIndex ),
    M_matrix( mbv.M_matrix )
{

}

MatrixBlockView::~MatrixBlockView()
{
    M_matrix.reset();
}

// ===================================================
// Set Methods
// ===================================================
void
MatrixBlockView::setup( const UInt& nRows,
                    const UInt& nCols,
                    const UInt& firstRow,
                    const UInt& firstCol,
                    const matrix_type& A )
{
    M_numRows          = nRows;
    M_numCols          = nCols;
    M_firstRowIndex    = firstRow;
    M_lastRowIndex     = firstRow+nRows-1;
    M_firstColumnIndex = firstCol;
    M_lastColumnIndex  = firstCol+nCols-1;
    M_matrix           = A.getMatrixPtr();
}

// ===================================================
// Get Methods
// ===================================================
UInt
MatrixBlockView::numRows() const
{
    return M_numRows;
}

UInt
MatrixBlockView::numCols() const
{
    return M_numCols;
}

UInt
MatrixBlockView::firstRowIndex() const
{
    return M_firstRowIndex;
}

UInt
MatrixBlockView::lastRowIndex() const
{
    return M_lastRowIndex;
}

UInt
MatrixBlockView::firstColumnIndex() const
{
    return M_firstColumnIndex;
}

UInt
MatrixBlockView::lastColumnIndex() const
{
    return M_lastColumnIndex;
}

} // Namespace LifeV
