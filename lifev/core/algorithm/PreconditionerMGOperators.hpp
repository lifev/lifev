//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This files contains the description and the implementation of operators for Multigrid

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 18-10-2011

    @mantainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

 */

#ifndef _PRECONDITIONERMGOPERATORS_HPP_
#define _PRECONDITIONERMGOPERATORS_HPP_

#include <boost/shared_ptr.hpp>
#include <iostream>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

namespace
{
typedef MapEpetra                        map_Type;
typedef RegionMesh<LinearTetra>          mesh_Type;
typedef MatrixEpetra<Real>               matrix_Type;
typedef boost::shared_ptr<matrix_Type>   matrixPtr_Type;
typedef FESpace<mesh_Type, map_Type>      FESpace_Type;
typedef boost::shared_ptr<FESpace_Type>  FESpacePtr_Type;
typedef boost::shared_ptr<Epetra_Comm>   commPtr_Type;
}


//! @name Restriction operators
//@{

//! Build a restriction operator that pass from P2 dofs to P1 dofs
/*!
    This method interpolates the dofs at the vertcies
    @param uFESPace finite element space for the velocity
    @param pFESPace finite element space for the pressure
    @param R matrix of the retriction to be filled
    @param comm communicator
 */
void
buildRestrictionP2ToP1 ( FESpacePtr_Type uFESpace,
                         FESpacePtr_Type pFESpace,
                         matrixPtr_Type  R )
{

    int numVelocityDofs = uFESpace->dof().numTotalDof();
    int numPressureDofs = pFESpace->dof().numTotalDof();
    int diff = numVelocityDofs - numPressureDofs;

    int lowerBoundX1 = 0                , upperBoundX1 = numPressureDofs - 1;
    int lowerBoundX2 = numPressureDofs  , upperBoundX2 = 2 * numPressureDofs - 1;
    int lowerBoundX3 = 2 * numPressureDofs, upperBoundX3 = 3 * numPressureDofs - 1;

    // Creating the row and column map
    Epetra_Map colMap ( * ( ( uFESpace->map() + pFESpace->map() ).map (Unique) ) );
    Epetra_Map rowMap ( * ( ( pFESpace->map() + pFESpace->map() + pFESpace->map() + pFESpace->map() ).map (Unique) ) ); // Only P1 elements
    //Epetra_Map rowMap( *( ( uFESpace->map() + pFESpace->map() ).map(Unique) ) );
    //Epetra_Map rowMap( 4*numPressureDofs, indexBase, *Comm );

    boost::shared_ptr<Epetra_FECrsMatrix> tmpR;
    tmpR.reset ( new Epetra_FECrsMatrix ( Copy, rowMap, 1, true ) );

    int NumMyElements = rowMap.NumMyElements();
    int* MyGlobalElements = rowMap.MyGlobalElements();

    double one = 1.0;

    int numErr = 0;

    int index[1];

    for ( int i = 0; i < NumMyElements; i++ )
    {
        if ( MyGlobalElements[i] >= lowerBoundX1 && MyGlobalElements[i] <= upperBoundX1 )
        {
            index[0] = MyGlobalElements[i];
            numErr += tmpR->InsertGlobalValues ( MyGlobalElements[i], 1, &one, index );
        }
        else if ( MyGlobalElements[i] >= lowerBoundX2 && MyGlobalElements[i] <= upperBoundX2 )
        {
            index[0] = MyGlobalElements[i] + diff;
            numErr += tmpR->InsertGlobalValues ( MyGlobalElements[i], 1, &one, index );
        }
        else if ( MyGlobalElements[i] >= lowerBoundX3 && MyGlobalElements[i] <= upperBoundX3 )
        {
            index[0] = MyGlobalElements[i] + 2 * diff;
            numErr += tmpR->InsertGlobalValues ( MyGlobalElements[i], 1, &one, index );
        }
        else if ( MyGlobalElements[i] > upperBoundX3 )
        {
            index[0] = MyGlobalElements[i] + 3 * diff;
            numErr += tmpR->InsertGlobalValues ( MyGlobalElements[i], 1, &one, index );
        }
    }

    tmpR->FillComplete (colMap, rowMap);
    R->swapCrsMatrix ( tmpR );
}

//@}

//! @name Prolongation operators
//@{

//! Build a prolongator operator that pass from P1 dofs to P2 dofs
/*!
    This method interpolates the dofs at the vertcies
    @param uFESPace finite element space for the velocity
    @param pFESPace finite element space for the pressure
    @param P matrix of the prolongator to be filled
    @param comm communicator
 */
void
buildProlongationP1ToP2 ( FESpacePtr_Type uFESpace,
                          FESpacePtr_Type pFESpace,
                          matrixPtr_Type  P )
{
    int numErr = 0;

    int    index [2];
    double values[2];
    values[0] = 1.0;
    values[1] = 1.0;

    // Space dimension
    const UInt fieldDim ( uFESpace->fieldDim() );

    // Number of elements
    const UInt numElements ( uFESpace->mesh()->numElements() );

    // Total number of velocity dofs
    const UInt numUTotalDof ( uFESpace->dof().numTotalDof() );

    // Total number of pressure dofs
    const UInt numPTotalDof ( pFESpace->dof().numTotalDof() );

    const UInt diff ( numUTotalDof - numPTotalDof );

    // Total number of dofs
    const UInt numTotalDof ( numUTotalDof * fieldDim + numPTotalDof );

    // Creating the row and column map
    Epetra_Map colMap ( * ( ( pFESpace->map() + pFESpace->map() + pFESpace->map() + pFESpace->map() ).map (Unique) ) );
    //Epetra_Map colMap( *( ( uFESpace->map() + pFESpace->map() ).map(Unique) ) );
    //Epetra_Map colMap( 4*numPTotalDof, indexBase, *Comm );
    Epetra_Map rowMap ( * ( ( uFESpace->map() + pFESpace->map() ).map (Unique) ) ); // Only P1 elements

    boost::shared_ptr<Epetra_FECrsMatrix> tmpP;
    tmpP.reset ( new Epetra_FECrsMatrix ( Copy, rowMap, 10, false ) );

    // Step 1: We take care of the velocity part

    // Loop over the elements
    for (UInt iterElement ( 0 ); iterElement < numElements; ++iterElement)
    {

        // Loop over the component
        for ( Int iComponent = 0; iComponent < fieldDim; ++iComponent )
        {
            // Loop over the local dofs
            // Create and insert rows in the matrix
            // Only the first node are important for P1
            for ( UInt iNode = 0 ; iNode < uFESpace->fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = uFESpace->fe().patternFirst ( iNode ); // iLocal = iNode

                UInt iGlobalRow = uFESpace->dof().localToGlobalMap ( iterElement, iLocal ) + iComponent * numUTotalDof;

                // If the row is on this processor
                if ( rowMap.MyGID ( iGlobalRow ) )
                {
                    // If it is a P1 node, it is interpolated using itself
                    if ( iLocal < 4 )
                    {
                        index[0]  = iGlobalRow - diff * iComponent;
                        numErr += tmpP->InsertGlobalValues ( iGlobalRow, 1, values, index );
                    }
                    else
                    {
                        if ( iLocal == 4 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 0 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 1 ) + iComponent * numPTotalDof;
                        }
                        else if ( iLocal == 5 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 1 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 2 ) + iComponent * numPTotalDof;
                        }
                        else if ( iLocal == 6 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 0 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 2 ) + iComponent * numPTotalDof;
                        }
                        else if ( iLocal == 7 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 0 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 3 ) + iComponent * numPTotalDof;
                        }
                        else if ( iLocal == 8 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 1 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 3 ) + iComponent * numPTotalDof;
                        }
                        else if ( iLocal == 9 )
                        {
                            index[0] = uFESpace->dof().localToGlobalMap ( iterElement, 2 ) + iComponent * numPTotalDof;
                            index[1] = uFESpace->dof().localToGlobalMap ( iterElement, 3 ) + iComponent * numPTotalDof;
                        }
                        numErr += tmpP->InsertGlobalValues ( iGlobalRow, 2, values, index );
                    }
                }
            } // end - loop over the dofs
        } // end - loop over the components
    } // end - loop over the elements

    // Step 2: Now we take care of the pressure block
    for ( Int iRow ( numUTotalDof * fieldDim ); iRow < numTotalDof; ++iRow )
    {
        // If the row is on this processor
        if ( rowMap.MyGID ( iRow ) )
        {
            index[0]  = iRow - diff * fieldDim;
            numErr += tmpP->InsertGlobalValues ( iRow, 1, values, index );
        }
    }

    tmpP->FillComplete (colMap, rowMap);
    tmpP->PutScalar ( 0.5 );

    // Step 3: Now we take care to set one in the proper position

    for ( Int iComp ( 0 ); iComp <= fieldDim; ++iComp )
        //<= Since we also take into account the pressure
    {
        Int start ( numUTotalDof * iComp );
        Int end  ( start + numPTotalDof );
        for ( Int iRow ( start ); iRow < end; ++iRow )
        {
            // If the row is on this processor
            if ( rowMap.MyGID ( iRow ) )
            {
                index[0]  = iRow - diff * iComp;
                numErr += tmpP->ReplaceGlobalValues ( iRow, 1, values, index );
            }
        }
    }

    P->swapCrsMatrix ( tmpP );
}

//@}

} // end of the namespace
#endif // _PRECONDITIONERMGOPERATORS_HPP_
