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
    int indexBase(0);

    // The matrix A should be filled
    if(srcBlock.getMatrixPtr()->Filled())
    {

        //Variables to store the informations
        int NumEntries;
        double* Values;
        int* Indices;
        int row(0);

        for(int i(0);i<srcBlock.getMatrixPtr()->NumGlobalRows();++i)
        {
            if((i>=srcBlock.firstRowIndex())&&(i<=srcBlock.lastRowIndex()))
            {
                row = srcBlock.getMatrixPtr()->LRID(i+indexBase);
                srcBlock.getMatrixPtr()->ExtractMyRowView(row, NumEntries, Values, Indices);

                int Indices2[NumEntries];
                double Values2[NumEntries];
                int NumEntries2(0);

                for(int j(0);j<NumEntries;++j)
                {
                    if((Indices[j]>=srcBlock.firstRowIndex())&&(Indices[j]<=srcBlock.lastRowIndex()))
                    {
                        Indices2[NumEntries2] = srcBlock.getMatrixPtr()->GCID(Indices[j]);
                        Values2[NumEntries2] = Values[j];
                        NumEntries2++;
                    }
                }
                destBlock.getMatrixPtr()->InsertGlobalValues(row,NumEntries2,Values2,Indices2);
            }
        }
	}
}

void createZeroBlock ( MatrixBlockView& B )
{

}

void createIdentityBlock ( MatrixBlockView& B )
{
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
