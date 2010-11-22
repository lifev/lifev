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
   @file MatrixBlockUtils.hpp
   @brief The file contains utility functions to manipulate BlockMatrixView objects

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-11-8
 */

#ifndef _MATRIXBLOCKUTILS_HPP_
#define _MATRIXBLOCKUTILS_HPP_

#include <lifemc/lifearray/MatrixBlockView.hpp>

namespace LifeV {

namespace MatrixBlockUtils {

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void copyBlock ( const MatrixBlockView& srcBlock,
                 MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param B Block where the data will be stored
*/
void createZeroBlock ( MatrixBlockView& B );

//! Copy the block specified in the given matrix
/*!
  @param B Block where the data will be stored
*/
void createIdentityBlock ( MatrixBlockView& B );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createDiagBlock ( const MatrixBlockView& srcBlock,
                       MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createInvDiagBlock ( const MatrixBlockView& srcBlock,
                          MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createUpperTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createLowerTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createLumpedBlock ( const MatrixBlockView& srcBlock,
                         MatrixBlockView& destBlock );

//! Copy the block specified in the given matrix
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
void createInvLumpedBlock ( const MatrixBlockView& srcBlock,
                            MatrixBlockView& destBlock );

} // namespace MatrixBlockUtils

} // namespace LifeV

#endif /* _MATRIXBLOCKUTILS_HPP_ */
