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
   @file MatrixEpetraStructuredUtility.hpp
   @brief The file contains utility functions to manipulate MatrixEpetraStructuredView objects

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-11-08

   @todo createDiagBlock() and createInvDiagBlock() can be reduced to a single function to avoid copy duplication with the aid of a bit of template meta programming. In fact, a lot of routines in this file bring back to a common block and a specialized work on the single line.
 */

#ifndef _MATRIXEPETRASTRUCTUREDUTILITY_HPP_
#define _MATRIXEPETRASTRUCTUREDUTILITY_HPP_

#include <boost/shared_ptr.hpp>
#include <lifev/core/array/MatrixBlockStructure.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredView.hpp>

namespace LifeV {

namespace MatrixEpetraStructuredUtility {

//! Copy the block specified in another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType>
void copyBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                 const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() <= destBlock.numRows(), "The destination block can not contain all the rows of the source block" );
	ASSERT( srcBlock.numColumns() <= destBlock.numColumns(), "The destination block can not contain all the columns of the source block" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    const Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    const Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    Int        srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    const Int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    const Int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    Int srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if(( static_cast<UInt>(srcRowElement)>=srcBlock.firstRowIndex())
           && ( static_cast<UInt>(srcRowElement)<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<Int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            Int numDestEntries(0);
            Int destRow(srcRowElement+rowsOffset);

            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if(( static_cast<UInt>(srcGlobalIndex)>=srcBlock.firstColumnIndex())
                   && ( static_cast<UInt>(srcGlobalIndex)<=srcBlock.lastColumnIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
            {
#ifdef NDEBUG
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
#else
                Int errorCode = destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
                ASSERT(errorCode >=0, " Error in block copy: insertion failed");
#endif
            }
            else
            {
#ifdef NDEBUG
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
#else
                Int errorCode = destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
                ASSERT(errorCode >=0, " Error in block copy: sum failed");
#endif
            }
        }
    }
}

//! Copy the block specified in another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType>
void copyBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                 boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    copyBlock( *srcBlock, *destBlock );
}

//! Create a block full of zeros
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createZeroBlock ( MatrixEpetraStructuredView<DataType>& /*destBlock*/ )
{
    // This method will maybe be replaced
    // by the method setBlockToZero
	ASSERT( false, "The method is not yet implemented");
}

//! Create a block full of zeros
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createZeroBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createZeroBlock( *destBlock );
}

//! Create a block with an identical value on the diagonal
/*!
  @param destBlock Block where the data will be stored
  @param diagonalValue Value to be inserted in the diagonal
*/
template< typename DataType >
void createScalarBlock ( const MatrixEpetraStructuredView<DataType>& destBlock, const DataType& diagonalValue )
{
    // SQUARE TEST
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK PTR TEST
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int destIndex(0);

    Int indexBase(0);

    Int firstRowIndex(destBlock.firstRowIndex()+indexBase);
    Int lastRowIndex(destBlock.lastRowIndex()+indexBase);
    Int firstColumnIndex(destBlock.firstColumnIndex()+indexBase);

    // Processor informations
    Int  numDestElements    = destBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* destGlobalElements = destBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    Int  destRowElement(0);

    for(Int i(0);i<numDestElements;++i)
    {
        destRowElement = destGlobalElements[i];

        // Test if the rows are in the block
        if((destRowElement>=firstRowIndex) && (destRowElement<=lastRowIndex))
        {
            destIndex = firstColumnIndex+destRowElement-firstRowIndex;
            destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRowElement,1,&diagonalValue,&destIndex);
        }
    }
}

//! Create a block with an identical value on the diagonal
/*!
  @param destBlock Block where the data will be stored
  @param diagonalValue Value to be inserted in the diagonal
*/
template< typename DataType >
void createScalarBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock, const DataType& diagonalValue )
{
    createScalarBlock( *destBlock, diagonalValue );
}

//! Create a block with ones on the diagonal
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createIdentityBlock ( const MatrixEpetraStructuredView<DataType>& destBlock )
{
    createScalarBlock(destBlock,1.0);
}

//! Create a block with ones on the diagonal
/*!
  @param destBlock Block where the data will be stored
*/
template< typename DataType >
void createIdentityBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createIdentityBlock( *destBlock );
}

//! Copy the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                       const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int indexBase(0);

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    Int srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            UInt diagIndex=srcRowElement-srcBlock.firstRowIndex();
            Int destRow = destBlock.firstRowIndex()+diagIndex;
            Int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType diagValue = 0.0;

            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    diagValue = srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

//! Copy the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createDiagBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                       boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createDiagBlock( *srcBlock, *destBlock );
}

//! Copy the inverse of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                          const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int indexBase(0);

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    Int srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            UInt diagIndex=srcRowElement-srcBlock.firstRowIndex();
            Int destRow = destBlock.firstRowIndex()+diagIndex;
            Int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType diagValue = 0.0;

            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    // ZERO ON DIAGONAL TEST
                	ASSERT( srcValues[j] != 0, "You cannot ask for inverse diagonal block when there are zeros on the diagonal" );

                    diagValue = 1./srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }

}

//! Copy the inverse of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvDiagBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                          boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createInvDiagBlock( *srcBlock, *destBlock );
}

//! Copy the inverse of the square root of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvSquaredDiagBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
								 const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int indexBase(0);

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt  srcRowElement(0);

    // Source informations handlers
    Int numSrcEntries;
    Real* srcValues;
    Int* srcIndices;
    UInt srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            UInt diagIndex=srcRowElement-srcBlock.firstRowIndex();
            Int destRow = destBlock.firstRowIndex()+diagIndex;
            Int destIndex = destBlock.firstColumnIndex()+diagIndex;
            Real diagValue = 0.0;

            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is on the diagonal of the source block
                if(srcGlobalIndex-srcBlock.firstColumnIndex()==diagIndex)
                {
                    // ZERO ON DIAGONAL TEST
                	ASSERT( srcValues[j] != 0, "You cannot ask for inverse squared diagonal block when there are zeros on the diagonal" );

                    diagValue = 1/sqrt(srcValues[j]);
                    j=numSrcEntries; //Exit the loop
                }
            }
            if(srcBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

//! Copy the inverse of the square root of the diagonal of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvSquaredDiagBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                                 boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createInvSquaredDiagBlock( *srcBlock, *destBlock );
}

//! Copy the upper part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createUpperTriangularBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                                  const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    Int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    Int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    UInt srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<Int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            Int numDestEntries(0);
            Int destRow(srcRowElement+rowsOffset);
            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()) &&
                   (srcGlobalIndex-srcBlock.firstColumnIndex()>=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
        }
    }
}

//! Copy the upper part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createUpperTriangularBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                                  boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createUpperTriangularBlock( *srcBlock, *destBlock );
}

//! Copy the lower part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLowerTriangularBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                                  const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    Int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    Int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    UInt srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            std::vector<Int> destIndices(numSrcEntries);
            std::vector<DataType> destValues(numSrcEntries);
            Int numDestEntries(0);
            Int destRow(srcRowElement+rowsOffset);
            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is:
                // a) in the block, and
                // b) on the lower triangular part of the block.
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()) &&
                   (srcGlobalIndex-srcBlock.firstColumnIndex()<=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,&destValues[0],&destIndices[0]);
        }
    }
}

//! Copy the lower part of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLowerTriangularBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                                  boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createLowerTriangularBlock( *srcBlock, *destBlock );
}


//! Copy the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLumpedBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                         const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int indexBase(0);

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    UInt srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            Int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            Int destRow = destBlock.firstRowIndex()+diagIndex;
            Int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType srcBlockRowSum = 0.0;
            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()+indexBase) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()+indexBase))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}

//! Copy the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createLumpedBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                         boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createLumpedBlock( *srcBlock, *destBlock );
}

//! Copy the inverse of the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvLumpedBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                            const MatrixEpetraStructuredView<DataType>& destBlock )
{
    // SQUARE TEST
	ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
	ASSERT( destBlock.numRows() == destBlock.numColumns() , "The destination block must be square" );

    // BLOCK COMPATIBILITY TEST
	ASSERT( srcBlock.numRows() == destBlock.numRows(), "The two blocks must have the same number of rows" );
	ASSERT( srcBlock.numColumns() == destBlock.numColumns(), "The two blocks must have the same number of columns" );

    // BLOCK PTR TEST
	ASSERT( srcBlock.matrixPtr() != 0 , "The source block does not have a valid pointer" );
	ASSERT( destBlock.matrixPtr() != 0 , "The destination block does not have a valid pointer" );

    Int indexBase(0);

    // Processor informations
    Int  numSrcElements    = srcBlock.matrixPtr()->matrixPtr()->RowMap().NumMyElements();
    Int* srcGlobalElements = srcBlock.matrixPtr()->matrixPtr()->RowMap().MyGlobalElements();
    UInt srcRowElement(0);

    // Source informations handlers
    Int numSrcEntries;
    DataType* srcValues;
    Int* srcIndices;
    UInt srcGlobalIndex(0);
    Int srcRow(0);

    for(Int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()+indexBase) && (srcRowElement<=srcBlock.lastRowIndex()+indexBase))
        {
            // Get the data of the row
            srcRow = srcBlock.matrixPtr()->matrixPtr()->LRID(srcRowElement);
            srcBlock.matrixPtr()->matrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            Int diagIndex=srcRowElement-srcBlock.firstRowIndex();
            Int destRow = destBlock.firstRowIndex()+diagIndex;
            Int destIndex = destBlock.firstColumnIndex()+diagIndex;
            DataType srcBlockRowSum = 0.0;
            for(Int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.matrixPtr()->matrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()+indexBase) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()+indexBase))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }

            // ZERO ON DIAGONAL TEST
        	ASSERT( srcBlockRowSum != 0, "You cannot ask for inverse lumped block when there are rows of zeros" );

            srcBlockRowSum = 1./srcBlockRowSum;
            if(destBlock.matrixPtr()->matrixPtr()->Map().MyGID(destRow))
                destBlock.matrixPtr()->matrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
            else
                destBlock.matrixPtr()->matrixPtr()->SumIntoGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}

//! Copy the inverse of the lumped version of the block specified to another block
/*!
  @param srcBlock Source block
  @param destBlock Destination block where the data will be stored
*/
template< typename DataType >
void createInvLumpedBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                            boost::shared_ptr< MatrixEpetraStructuredView<DataType> > destBlock )
{
    createInvLumpedBlock( *srcBlock, *destBlock );
}

//! Create a new matrix from the block specified
/*!
  @param srcBlock Source block
  @param dstMatrix Pointer to be initialized with a new matrix
  @param rowMap Row map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param closeMatrix If closeMatrix is equal to true, globalAssemble will be called.
  @warning This method is only intended to be used with square blocks!
*/
template< typename DataType>
void createMatrixFromBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                             boost::shared_ptr<MatrixEpetraStructured<DataType> >& destMatrix,
                             const MapEpetra& rowMap,
                             bool closeMatrix=true )
{
    // SQUARE TEST
    ASSERT( srcBlock.numRows() == srcBlock.numColumns() , "The source block must be square" );
    ASSERT( srcBlock.numRows() == rowMap.mapSize(), "The two blocks must have the same number of rows" );

    // Create destination matrix
    destMatrix.reset( new MatrixEpetraStructured<DataType>( rowMap, srcBlock.matrixPtr()->meanNumEntries() ) );

    // Copy the entries
    copyBlock( srcBlock, *( destMatrix->block( 0, 0 ) ) );

    // Close the matrix if requested
    if( closeMatrix )
    {
        destMatrix->globalAssemble();
    }
}

//! Create a new matrix from the block specified
/*!
  @param srcBlock Source block
  @param dstMatrix Pointer to be initialized with a new matrix
  @param rowMap Row map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param closeMatrix If closeMatrix is equal to true, globalAssemble will be called.
  @warning This method is only intended to be used with square blocks!
*/
template< typename DataType>
void createMatrixFromBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                             boost::shared_ptr<MatrixEpetraStructured<DataType> >& destMatrix,
                             const MapEpetra& rowMap,
                             bool closeMatrix=true )
{
    createMatrixFromBlock( *srcBlock, destMatrix, rowMap, closeMatrix );
}

//! Create a new matrix from the block specified (for rectangle matrix)
/*!
  @param srcBlock Source block
  @param dstMatrix Pointer to be initialized with a new matrix
  @param domainMap domain map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param rangeMap range map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param closeMatrix If closeMatrix is equal to true, globalAssemble will be called.
  @warning This method is only intended to be used with square blocks!
*/
template< typename DataType>
void createMatrixFromBlock ( const MatrixEpetraStructuredView<DataType>& srcBlock,
                             boost::shared_ptr<MatrixEpetraStructured<DataType> >& destMatrix,
                             boost::shared_ptr<MapEpetra> domainMap,
                             boost::shared_ptr<MapEpetra> rangeMap,
                             bool closeMatrix=true )
{


    // Create destination matrix
    if( domainMap->mapSize() > rangeMap->mapSize() )
        destMatrix.reset( new MatrixEpetraStructured<DataType>( *domainMap, srcBlock.matrixPtr()->meanNumEntries() ) );
    else
        destMatrix.reset( new MatrixEpetraStructured<DataType>( *rangeMap, srcBlock.matrixPtr()->meanNumEntries() ) );

    // Copy the entries
    copyBlock( srcBlock, *( destMatrix->block( 0, 0 ) ) );

    // Close the matrix if requested
    if( closeMatrix )
    {
        destMatrix->globalAssemble( domainMap, rangeMap );
    }
}

//! Create a new matrix from the block specified
/*!
  @param srcBlock Source block
  @param dstMatrix Pointer to be initialized with a new matrix
  @param domainMap domain map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param rangeMap range map. The column map will be defined in MatrixEpetraStructured<DataType>::GlobalAssemble(...,...)
  @param closeMatrix If closeMatrix is equal to true, globalAssemble will be called.
  @warning This method is only intended to be used with square blocks!
*/
template< typename DataType>
void createMatrixFromBlock ( boost::shared_ptr< MatrixEpetraStructuredView<DataType> > srcBlock,
                             boost::shared_ptr<MatrixEpetraStructured<DataType> >& destMatrix,
                             boost::shared_ptr<MapEpetra> domainMap,
                             boost::shared_ptr<MapEpetra> rangeMap,
                             bool closeMatrix=true )
{
    createMatrixFromBlock( *srcBlock, destMatrix, domainMap, rangeMap, closeMatrix );
}

//! Create a block view using an unstructured matrix and block structure informations
/*!
  @param matrixPtr Pointer on an unstructured matrix
  @param blockStructure Structure to be used to extract block view
  @param rowIndex Row position of the block in the matrix
  @param columnIndex Column position of the block in the matrix
*/
template <typename DataType>
boost::shared_ptr< MatrixEpetraStructuredView<DataType> >
createBlockView( boost::shared_ptr<MatrixEpetra<DataType> > matrixPtr,
                 const MatrixBlockStructure& blockStructure,
                 const UInt& rowIndex,
                 const UInt& columnIndex )
{
    ASSERT( matrixPtr->matrixPtr()->NumGlobalCols() == blockStructure.numRows(), " Incompatible block structure (global size does not match) " );
    ASSERT( matrixPtr->matrixPtr()->NumGlobalRows() == blockStructure.numColumns(), " Incompatible block structure (global size does not match) " );

    boost::shared_ptr< MatrixEpetraStructuredView<DataType> > matrixBlockView( new MatrixEpetraStructuredView<DataType> );

    matrixBlockView->setup( blockStructure.rowBlockFirstIndex( rowIndex ),
                            blockStructure.columnBlockFirstIndex( columnIndex ),
                            blockStructure.blockNumRows( rowIndex ),
                            blockStructure.blockNumColumns( columnIndex ),
                            matrixPtr->matrixPtr().get() );

    return matrixBlockView;
}

//! Fill a block view using an unstructured matrix and block structure informations
/*!
  @param matrixPtr Pointer on an unstructured matrix
  @param blockStructure Structure to be used to extract block view
  @param rowIndex Row position of the block in the matrix
  @param columnIndex Column position of the block in the matrix
  @param blockView block view to be filled with informations
*/
template <typename DataType>
void
fillBlockView( boost::shared_ptr<MatrixEpetra<DataType> > matrixPtr,
               const MatrixBlockStructure& blockStructure,
               const UInt& rowIndex,
               const UInt& columnIndex,
               MatrixEpetraStructuredView<DataType>& blockView )
{
    ASSERT( matrixPtr->matrixPtr()->NumGlobalCols() == blockStructure.numRows(), " Incompatible block structure (global size does not match) " );
    ASSERT( matrixPtr->matrixPtr()->NumGlobalRows() == blockStructure.numColumns(), " Incompatible block structure (global size does not match) " );

    blockView.setup( blockStructure.rowBlockFirstIndex( rowIndex ),
                     blockStructure.columnBlockFirstIndex( columnIndex ),
                     blockStructure.blockNumRows( rowIndex ),
                     blockStructure.blockNumColumns( columnIndex ),
                     matrixPtr.get() );
}

//! Fill a block view using an unstructured matrix and block structure informations
/*!
  @param matrixPtr Pointer on an unstructured matrix
  @param blockStructure Structure to be used to extract block view
  @param rowIndex Row position of the block in the matrix
  @param columnIndex Column position of the block in the matrix
  @param blockView block view to be filled with informations
*/
template <typename DataType>
void
fillBlockView( boost::shared_ptr<MatrixEpetra<DataType> > matrixPtr,
               const MatrixBlockStructure& blockStructure,
               const UInt& rowIndex,
               const UInt& columnIndex,
               boost::shared_ptr< MatrixEpetraStructuredView<DataType> >& blockView )
{
    if( blockView.get() == 0 )
    {
        blockView.reset( new MatrixEpetraStructuredView<DataType> );
    }
    fillBlockView( matrixPtr, blockStructure, rowIndex, columnIndex, *blockView );
}

} // namespace MatrixEpetraStructuredUtility

} // namespace LifeV

#endif /*_MATRIXEPETRASTRUCTUREDUTILITY_HPP_ */
