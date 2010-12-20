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
    @brief File containing method to insert local contributions in the global system.

    @author Alessandro Veneziani
    @date 31-01-2001

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    To delete:
    extract_vec
    assemb_mat
    assemb_vec

    To create:
    namespace Assemble

    To change:
    template parameter
    DOF -> globalMatrix and localMatrix
 */


#ifndef _ASSEMBLY
#define _ASSEMBLY

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/CurrentFE.hpp>
#include <life/lifefem/localDofPattern.hpp>

#include <vector>

namespace LifeV
{



//! Assembly procedure for vectors
/*!
  This function allows the user to transfer
  a local contribution to a global vector.
 */
template <typename DofType>
void
assembleVector( EpetraVector&    globalVector,
                ElemVec&         localVector,
                const CurrentFE& currentFE,
                const DofType    dof,
                Int              block,
                Int              offset = 0 )

{
    UInt elementID = currentFE.currentLocalId();
    assembleVector(globalVector,elementID,localVector,currentFE.nbFEDof(),dof,block,offset);
}

//! Assembly procedure for vectors
/*!
  This function can transfer a local
  contribution to a global vector.
 */
template <typename DofType>
void
assembleVector( EpetraVector&    globalVector,
                const UInt&      elementID,
                ElemVec&         localVector,
                const UInt&      feNbDof,
                const DofType&   dof,
                Int              block,
                Int              offset = 0)

{
    ElemVec::vector_view localView = localVector.block( block );

    UInt iGlobalID;

    for ( UInt i(0) ; i < feNbDof ; ++i )
    {
        iGlobalID = dof.localToGlobal( elementID, i + 1 ) + offset;
        globalVector.sumIntoGlobalValues( iGlobalID, localView( i ) );
    }
}

//! Assembly procedure for the matrix
/*!
  This method allows to transfer local contributions
  to a global matrix.
 */
template <typename DofType>
void
assembleMatrix( EpetraMatrix<Real>&   globalMatrix,
                const UInt&      	  elementID,
                ElemMat&          	  localMatrix,
                const UInt&           feNbDof,
                const DofType&    	  dof,
                Int                   iblock,
                Int                   jblock,
                Int              	  iOffset,
                Int              	  jOffset)

{

    ElemMat::matrix_view localView = localMatrix.block( iblock, jblock );

    assembleMatrix( globalMatrix,
                    elementID,
                    elementID,
                    localView,
                    feNbDof,
                    feNbDof,
                    dof,
                    dof, iOffset, jOffset);
}

//! Assembly procedure for the matrix
/*!
  This method allows to transfer local contributions
  to a global matrix.
 */
template <typename DofType>
void
assembleMatrix( EpetraMatrix<Real>& globalMatrix,
                ElemMat&            localMatrix,
                const CurrentFE&    currentFE,
                const DofType&      dof,
                Int                 iblock,
                Int                 jblock,
                Int                 iOffset,
                Int                 jOffset)

{
    return assembleMatrix( globalMatrix, localMatrix, currentFE, currentFE, dof, dof,
                           iblock, jblock, iOffset, jOffset);
}

//! Assembly procedure for the matrix
/*!
  This method allows to transfer local contributions
  to a global matrix.
 */
template <typename DofType1, typename DofType2, typename LocalMatrixType>
void
assembleMatrix( EpetraMatrix<Real>& globalMatrix,
                UInt const&         elementID1,
                UInt const&         elementID2,
                LocalMatrixType&    localMatrix,
                const CurrentFE&    currentFE1,
                const CurrentFE&    currentFE2,
                const DofType1&     dof1,
                const DofType2&     dof2,
                Int                 iOffset,
                Int                 jOffset)

{
    assembleMatrix(globalMatrix, elementID1, elementID2, localMatrix, currentFE1.nbFEDof(), currentFE2.nbFEDof(),
                   dof1, dof2, iOffset, jOffset);

}

//! Assembly procedure for the matrix
/*!
  This method allows to transfer local contributions
  to a global matrix.
 */
template <typename DofType1, typename DofType2, typename LocalMatrixType>
void
assembleMatrix( EpetraMatrix<Real>&   globalMatrix,
                UInt const&           elementID1,
                UInt const&           elementID2,
                LocalMatrixType&      localMatrix,
                const UInt&           fe1NbDof,
                const UInt&           fe2NbDof,
                const DofType1&       dof1,
                const DofType2&       dof2,
                Int                   iOffset,
                Int                   jOffset )

{
    // Global ID of the dofs
    std::vector<Int> iList(fe1NbDof);
    std::vector<Int> jList(fe2NbDof);

    // Raw data to insert in the matrix
    Real* matPtr[fe2NbDof];


    for ( UInt k1 (0) ; k1 < fe1NbDof ; k1++ )
    {
        iList[k1] = dof1.localToGlobal( elementID1, k1 + 1 ) - 1 + iOffset ;
    }

    for ( UInt k2 (0) ; k2 < fe2NbDof ; k2++ )
    {
        jList[k2]  = dof2.localToGlobal( elementID2, k2 + 1 ) - 1 + jOffset ;
        matPtr[k2] = &(localMatrix(static_cast<UInt>(0),k2));
    }

    assert(localMatrix.indexij( Int (1), Int(0) ) == 1);

    globalMatrix.addToCoefficients( fe1NbDof, fe2NbDof, iList, jList, matPtr, Epetra_FECrsMatrix::COLUMN_MAJOR );
}

//! Assembly procedure for the matrix
/*!
  This method allows to transfer local contributions
  to a global matrix.
 */
template <typename DofType1, typename DofType2>
void
assembleMatrix( EpetraMatrix<Real>& globalMatrix,
                ElemMat&            localMatrix,
                const CurrentFE&    currentFE1,
                const CurrentFE&    currentFE2,
                const DofType1&     dof1,
                const DofType2&     dof2,
                Int                 iblock,
                Int                 jblock,
                Int                 iOffset,
                Int                 jOffset )

{
    ElemMat::matrix_view localView = localMatrix.block( iblock, jblock );

    UInt elementID1 = currentFE1.currentLocalId();
    UInt elementID2 = currentFE2.currentLocalId();

    assembleMatrix( globalMatrix, elementID1, elementID2, localView,
                    currentFE1, currentFE2, dof1,  dof2, iOffset, jOffset );

    return;

}

//! Assembly procedure for the transposed matrix
/*!
  This method allows to transfer local contributions
  to a transposed global matrix.

  The coefficient is used to multiply the values.
 */
template <typename DofType1, typename DofType2>
void
assembleTransposeMatrix( EpetraMatrix<Real>&   globalMatrix,
                         Real                  coefficient,
                         ElemMat&              localMatrix,
                         const CurrentFE&      currentFE1,
                         const CurrentFE&      currentFE2,
                         const DofType1&       dof1,
                         const DofType2&       dof2,
                         Int                   iblock,
                         Int                   jblock,
                         Int                   iOffset ,
                         Int                   jOffset )

{
    ElemMat::matrix_type localView(localMatrix.block( jblock, iblock ));
    localView *= coefficient;

    Int i, j;
    UInt k1, k2;

    UInt elementID1 = currentFE1.currentLocalId();
    UInt elementID2 = currentFE2.currentLocalId();

    std::vector<Int> ilist(currentFE1.nbFEDof());
    std::vector<Int> jlist(currentFE2.nbFEDof());

    Real* matPtr[currentFE1.nbFEDof()];

    for ( k1 = 0 ; k1 < currentFE1.nbFEDof() ; k1++ )
    {
        i =  k1;
        ilist[k1] = dof1.localToGlobal( elementID1, i + 1 ) - 1 + iOffset ;
        matPtr[k1] = &(localView(0,i));
    }

    for ( k2 = 0 ; k2 < currentFE2.nbFEDof() ; k2++ )
    {
        j = k2;
        jlist[k2]  = dof2.localToGlobal( elementID2, j + 1 ) - 1 + jOffset ;
    }

    assert(localView.indexij( Int (1), Int(0) ) == 1);

    globalMatrix.addToCoefficients( currentFE1.nbFEDof(), currentFE2.nbFEDof(),
                              ilist, jlist, matPtr, Epetra_FECrsMatrix::ROW_MAJOR );

}

template <typename DOF, typename ElemVec>
void
extract_vec( EpetraVector& V,
             ElemVec& elvec,
             const LocalDofPattern& fe,
             const DOF& dof,
             const UInt feId,
             Int iblock )
{
    UInt totdof (dof.numTotalDof());
    typename ElemVec::vector_view vec = elvec.block( iblock );
    UInt ig;
    for ( UInt i (0) ; i < fe.nbLocalDof() ; ++i )
    {
        ig = dof.localToGlobal( feId, i + 1 ) + iblock * totdof;
        vec( i ) = V[ ig ];
    }
}

}

#endif
