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

#include "MatrixBlockUtils.hpp"

namespace LifeV {

namespace MatrixBlockUtils {

void copyBlock ( const MatrixBlockView& srcBlock,
                 MatrixBlockView& destBlock )
{
    // BLOCK COMPATIBILITY TEST

    int indexBase(0);

    // The matrix A should be filled
    if(srcBlock.getMatrixPtr()->Filled())
    {

        //Variables to store the informations
        int numSrcEntries;
        double* srcValues;
        int* srcIndices;
        int srcRow(0);

        //Offset between the first row/column of the source and destination blocks
        int rowsOffset(destBlock.firstRowIndex()-srcBlock.firstRowIndex());
        int columnsOffset(destBlock.firstColumnIndex()-srcBlock.firstColumnIndex());

        for(int i(0);i<srcBlock.getMatrixPtr()->NumGlobalRows();++i)
        {
            if((i>=srcBlock.firstRowIndex())&&(i<=srcBlock.lastRowIndex()))
            {
                srcRow = srcBlock.getMatrixPtr()->LRID(i+indexBase);
                srcBlock.getMatrixPtr()->ExtractMyRowView(srcRow, numSrcEntries, srcValues, srcIndices);

                int destIndices[numSrcEntries];
                double destValues[numSrcEntries];
                int numDestEntries(0);
                int destRow(destBlock.getMatrixPtr()->LRID(i+rowsOffset+indexBase));

                for(int j(0);j<numSrcEntries;++j)
                {
                    if((srcIndices[j]>=srcBlock.firstColumnIndex())&&(srcIndices[j]<=srcBlock.lastColumnIndex()))
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
}

void createZeroBlock ( MatrixBlockView& B )
{
    // This method will maybe be replaced
    // by the method setBlockToZero
}

void createIdentityBlock ( MatrixBlockView& B )
{
    // SQUARE TEST

    int from(0);
    int to(10);
    //B.getMatrixPtr()->insertOneDiagonal(from,to);
}

void createDiagBlock ( const MatrixBlockView& srcBlock,
                       MatrixBlockView& destBlock )
{
    // Usefull command
    //M->set_mat_inc(i,j,value);
}

void createInvDiagBlock ( const MatrixBlockView& srcBlock,
                          MatrixBlockView& destBlock )
{
    // Usefull command
    //M->set_mat_inc(i,j,1/value);
}

void createUpperTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock )
{

}

void createLowerTriangularBlock ( const MatrixBlockView& srcBlock,
                                  MatrixBlockView& destBlock )
{

}

void createLumpedBlock ( const MatrixBlockView& srcBlock,
                         MatrixBlockView& destBlock )
{
    // Usefull command
    //M->set_mat_inc(i,j,value);
}

void createInvLumpedBlock ( const MatrixBlockView& srcBlock,
                            MatrixBlockView& destBlock )
{
    // Usefull command
    //M->set_mat_inc(i,j,1/value);
}

} // namespace MatrixBlockUtils

} // namespace LifeV
