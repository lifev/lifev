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
   @file MatrixBlockUtils.cpp
   @brief The file contains utility functions to manipulate BlockMatrixView objects

   @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   @date 2010-11-8
 */

#include <math.h>
#include "MatrixBlockUtils.hpp"

namespace LifeV {

namespace MatrixBlockUtils {

void copyBlock ( const MatrixBlockView& srcBlock,
                 MatrixBlockView& destBlock )
{
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcGlobalIndex(0);
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            int destIndices[numSrcEntries];
            double destValues[numSrcEntries];
            int numDestEntries(0);
            int destRow(srcRowElement+rowsOffset);
            for(int j(0);j<numSrcEntries;++j)
            {
                srcGlobalIndex = srcBlock.getMatrixPtr()->GCID(srcIndices[j]);

                // Test if the coefficient is in the block
                if((srcGlobalIndex>=srcBlock.firstColumnIndex()) &&
                   (srcGlobalIndex<=srcBlock.lastColumnIndex()))
                {
                    destIndices[numDestEntries] = srcGlobalIndex+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            if(srcBlock.getMatrixPtr()->Map().MyGID(destRow))
                destBlock.getMatrixPtr()->InsertGlobalValues(destRow,numDestEntries,destValues,destIndices);
            else
                destBlock.getMatrixPtr()->SumIntoGlobalValues(destRow,numDestEntries,destValues,destIndices);
        }
    }
}

void createZeroBlock ( MatrixBlockView& destBlock )
{
    // This method will maybe be replaced
    // by the method setBlockToZero
}

void createIdentityBlock ( MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK PTR TEST

    int destIndex(0);
    double one(1.0);

    int firstRowIndex(destBlock.firstRowIndex());
    int lastRowIndex(destBlock.lastRowIndex());
    int firstColumnIndex(destBlock.firstColumnIndex());

    // Processor informations
    int  numDestElements    = destBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* destGlobalElements = destBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  destRowElement(0);

    for(int i(0);i<numDestElements;++i)
    {
        destRowElement = destGlobalElements[i];

        // Test if the rows are in the block
        if((destRowElement>=firstRowIndex) && (destRowElement<=lastRowIndex))
        {
            destIndex = firstColumnIndex+destRowElement-firstRowIndex;
            destBlock.getMatrixPtr()->InsertGlobalValues(destRowElement,1,&one,&destIndex);
        }
    }
}

void createDiagBlock ( const MatrixBlockView& srcBlock,
                       MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST
    // ZERO ON DIAGONAL TEST

    int indexBase(0);

    int destRow(0);
    int destIndex(0);
    double diagValue(0.0);
    int diagIndex(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            // index of the diagonal
            diagIndex=srcRowElement-srcBlock.firstRowIndex();
            diagValue = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is on the diagonal of the source block
                if(srcIndices[j]-srcBlock.firstColumnIndex()==diagIndex)
                {
                    diagValue = srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }

            destRow   = destBlock.firstRowIndex()+diagIndex+indexBase;
            destIndex = destBlock.firstColumnIndex()+diagIndex+indexBase;
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

void createInvDiagBlock ( const MatrixBlockView& srcBlock,
                          MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST
    // ZERO ON DIAGONAL TEST

    int indexBase(0);

    int destRow(0);
    int destIndex(0);
    double diagValue(0.0);
    int diagIndex(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            // index of the diagonal
            diagIndex=srcRowElement-srcBlock.firstRowIndex();
            diagValue = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is on the diagonal of the source block
                if(srcIndices[j]-srcBlock.firstColumnIndex()==diagIndex)
                {
                    diagValue = 1/srcValues[j];
                    j=numSrcEntries; //Exit the loop
                }
            }

            destRow   = destBlock.firstRowIndex()+diagIndex+indexBase;
            destIndex = destBlock.firstColumnIndex()+diagIndex+indexBase;
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,1,&diagValue,&destIndex);
        }
    }
}

void createUpperTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            int destIndices[numSrcEntries];
            double destValues[numSrcEntries];
            int numDestEntries(0);
            int destRow(destBlock.getMatrixPtr()->LRID(i+rowsOffset+indexBase));
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcIndices[j]>=srcBlock.firstColumnIndex()) &&
                   (srcIndices[j]<=srcBlock.lastColumnIndex()) &&
                   (srcIndices[j]-srcBlock.firstColumnIndex()>=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcBlock.getMatrixPtr()->GCID(srcIndices[j])+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,numDestEntries,destValues,destIndices);
        }
    }
}

void createLowerTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            int destIndices[numSrcEntries];
            double destValues[numSrcEntries];
            int numDestEntries(0);
            int destRow(destBlock.getMatrixPtr()->LRID(i+rowsOffset+indexBase));
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcIndices[j]>=srcBlock.firstColumnIndex()) &&
                   (srcIndices[j]<=srcBlock.lastColumnIndex()) &&
                   (srcIndices[j]-srcBlock.firstColumnIndex()<=srcRowElement-srcBlock.firstRowIndex()))
                {
                    destIndices[numDestEntries] = srcBlock.getMatrixPtr()->GCID(srcIndices[j])+columnsOffset;
                    destValues[numDestEntries] = srcValues[j];
                    numDestEntries++;
                }
            }
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,numDestEntries,destValues,destIndices);
        }
    }
}

void createLumpedBlock ( const MatrixBlockView& srcBlock,
                         MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    int destRow(0);
    int destIndex(0);
    int diagIndex(0);
    double srcBlockRowSum(0.0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            diagIndex=srcRowElement-srcBlock.firstRowIndex();
            srcBlockRowSum = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcIndices[j]>=srcBlock.firstColumnIndex()) &&
                   (srcIndices[j]<=srcBlock.lastColumnIndex()))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }
            destRow   = destBlock.firstRowIndex()+diagIndex+indexBase;
            destIndex = destBlock.firstColumnIndex()+diagIndex+indexBase;
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}

void createInvLumpedBlock ( const MatrixBlockView& srcBlock,
                            MatrixBlockView& destBlock )
{
    // SQUARE TEST
    // BLOCK COMPATIBILITY TEST
    // BLOCK PTR TEST

    int indexBase(0);

    // Processor informations
    int  numSrcElements    = srcBlock.getMatrixPtr()->RowMap().NumMyElements();
    int* srcGlobalElements = srcBlock.getMatrixPtr()->RowMap().MyGlobalElements();
    int  srcRowElement(0);

    //Offset between the first row/column of the source and destination blocks
    int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
    int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

    // Source informations handlers
    int numSrcEntries;
    double* srcValues;
    int* srcIndices;
    int srcRow(0);

    int destRow(0);
    int destIndex(0);
    int diagIndex(0);
    double srcBlockRowSum(0.0);

    for(int i(0);i<numSrcElements;++i)
    {
        // Collecting the data from the source
        srcRowElement = srcGlobalElements[i];

        // Test if the rows are in the source block
        if((srcRowElement>=srcBlock.firstRowIndex()) && (srcRowElement<=srcBlock.lastRowIndex()))
        {
            // Get the data of the row
            srcRow = srcBlock.getMatrixPtr()->LRID(srcRowElement+indexBase);
            srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

            diagIndex=srcRowElement-srcBlock.firstRowIndex();
            srcBlockRowSum = 0.0;
            for(int j(0);j<numSrcEntries;++j)
            {
                // Test if the coefficient is:
                // a) in the block, and
                // b) on the upper triangular part of the block.
                if((srcIndices[j]>=srcBlock.firstColumnIndex()) &&
                   (srcIndices[j]<=srcBlock.lastColumnIndex()))
                {
                    srcBlockRowSum += abs(srcValues[j]);
                }
            }
            srcBlockRowSum = 1/srcBlockRowSum;
            destRow   = destBlock.firstRowIndex()+diagIndex+indexBase;
            destIndex = destBlock.firstColumnIndex()+diagIndex+indexBase;
            destBlock.getMatrixPtr()->InsertGlobalValues(destRow,1,&srcBlockRowSum,&destIndex);
        }
    }
}

} // namespace MatrixBlockUtils

} // namespace LifeV
