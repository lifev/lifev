/*
 * BlockOperator.cpp
 *
 *  Created on: Sep 20, 2010
 *      Author: uvilla
 */

#include <lifev/core/linear_algebra/BlockOperator.hpp>

namespace LifeV
{
namespace Operators
{
BlockOperator::BlockOperator():
        M_name("BlockOperator"),
        M_useTranspose(false),
        M_structure(NoStructure)
{

}

void BlockOperator::setUp(const boost::shared_ptr<BlockEpetra_Map> & map, const commPtr_Type & comm)
{
    M_comm = comm;

    M_nBlockRows = map->nBlocks();
    M_nBlockCols = M_nBlockRows;
    M_domainMap = map;
    M_rangeMap = M_domainMap;
    M_oper.resize(M_nBlockRows, M_nBlockCols);
}


//! SetUp for a "rectangular operator"
void BlockOperator::setUp(const boost::shared_ptr<BlockEpetra_Map> & domainMap,
                             const boost::shared_ptr<BlockEpetra_Map> & rangeMap,
                             const commPtr_Type & comm)
{
    M_comm = comm;

    M_nBlockRows = rangeMap->nBlocks();
    M_nBlockCols = domainMap->nBlocks();

    M_domainMap = domainMap;
    M_rangeMap =  rangeMap;

    M_oper.resize(M_nBlockRows, M_nBlockCols);

}

//! SetUp when the operator is given like a boost::matrix
void BlockOperator::setUp(const operatorPtrContainer_Type & blockOper, const commPtr_Type & comm)
{
    M_comm = comm;
    M_nBlockRows = blockOper.size1();
    M_nBlockCols = blockOper.size2();

    BlockEpetra_Map::mapPtrContainer_Type rangeBlockMaps(M_nBlockRows);
    BlockEpetra_Map::mapPtrContainer_Type domainBlockMaps(M_nBlockCols);

    for(UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        for(UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        {
            if(blockOper(iblock,jblock) != 0 && rangeBlockMaps[iblock]==0)
            {
                rangeBlockMaps[iblock].reset(new Epetra_Map(blockOper(iblock,jblock)->OperatorRangeMap() ));
                jblock = M_nBlockCols;
            }
        }

    for(UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        for(UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        {
            if(blockOper(iblock,jblock) != 0 && domainBlockMaps[jblock]==0)
            {
                domainBlockMaps[jblock].reset(new Epetra_Map(blockOper(iblock,jblock)->OperatorDomainMap() ));
                iblock = M_nBlockRows;
            }
        }

    M_domainMap.reset(new BlockEpetra_Map(domainBlockMaps));
    M_rangeMap.reset(new BlockEpetra_Map(rangeBlockMaps));

    M_oper = blockOper;
    fillComplete();

}

//! set a component of the block operator
void BlockOperator::setBlock(UInt iblock, UInt jblock, const operatorPtr_Type & operBlock)
{
    ASSERT_PRE(M_rangeMap->blockMap(iblock)->PointSameAs(operBlock->OperatorRangeMap()), "Wrong range map");
    ASSERT_PRE(M_domainMap->blockMap(jblock)->PointSameAs(operBlock->OperatorDomainMap()), "Wrong domain map");

    M_oper(iblock, jblock) = operBlock;
}


void BlockOperator::fillComplete()
{
    // Check for Structure with respect to the main diagonal
    bool thereAreUpperDiagonalBlocks(false);
    for(UInt iblock = 0; iblock < M_nBlockRows - 1; ++iblock)
        for(UInt jblock = iblock+1; jblock < M_nBlockCols; ++jblock)
            if(M_oper(iblock,jblock).get() != 0 && (M_oper(iblock, jblock)->HasNormInf() ? M_oper(iblock, jblock)->NormInf()!=0 : true) )
                thereAreUpperDiagonalBlocks = true;

    bool thereAreLowerDiagonalBlocks(false);
    for(UInt iblock = 1; iblock < M_nBlockRows; ++iblock)
        for(UInt jblock = 0; jblock < iblock; ++jblock)
            if(M_oper(iblock,jblock).get() != 0  && (M_oper(iblock, jblock)->HasNormInf() ? M_oper(iblock, jblock)->NormInf()!=0 : true) )
                thereAreLowerDiagonalBlocks = true;

    // Check for Structure with respect to the antidiagonal
    bool thereAreUpperAntiDiagonalBlocks(false);
    for(UInt iblock = 0; iblock < M_nBlockRows-1; ++iblock)
        for(UInt jblock = 0; jblock < iblock; ++jblock)
            if(M_oper(iblock,jblock).get() != 0 && (M_oper(iblock, jblock)->HasNormInf() ? M_oper(iblock, jblock)->NormInf()!=0 : true) )
                thereAreUpperAntiDiagonalBlocks = true;

    bool thereAreLowerAntiDiagonalBlocks(false);
    for(UInt iblock = 1; iblock < M_nBlockRows; ++iblock)
        for(UInt jblock = iblock+1; jblock < M_nBlockCols; ++jblock)
            if(M_oper(iblock,jblock).get() != 0  && (M_oper(iblock, jblock)->HasNormInf() ? M_oper(iblock, jblock)->NormInf()!=0 : true) )
                thereAreLowerAntiDiagonalBlocks = true;

    // Filling the empty blocks with null operators
    for(UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        for(UInt jblock = 0; jblock < M_nBlockCols; ++jblock)
        {
            if(M_oper(iblock,jblock).get() == 0)
            {
                NullOperator * nullOp(new NullOperator);
                nullOp->setUp(M_domainMap->blockMap(jblock), M_rangeMap->blockMap(iblock) );
                M_oper(iblock,jblock).reset(nullOp);
            }
            ASSERT(M_rangeMap->blockMap(iblock)->PointSameAs(M_oper(iblock,jblock)->OperatorRangeMap()), "Wrong range map");
            ASSERT(M_domainMap->blockMap(jblock)->PointSameAs(M_oper(iblock,jblock)->OperatorDomainMap()), "Wrong domain map");
        }

    if(M_nBlockRows != M_nBlockCols)
    {
        M_structure = Rectangular;
        return;
    }

    if(!thereAreLowerDiagonalBlocks && !thereAreUpperDiagonalBlocks)
    {
        M_structure = Diagonal;
        return;
    }
    if(thereAreLowerDiagonalBlocks && !thereAreUpperDiagonalBlocks)
    {
        M_structure = LowerTriangular;
        return;
    }
    if(!thereAreLowerDiagonalBlocks && thereAreUpperDiagonalBlocks)
    {
        M_structure = UpperTriangular;
        return;
    }
    if(!thereAreLowerAntiDiagonalBlocks && !thereAreUpperAntiDiagonalBlocks)
    {
        M_structure = AntiDiagonal;
        return;
    }
    if(thereAreLowerAntiDiagonalBlocks && !thereAreUpperAntiDiagonalBlocks)
    {
        M_structure = LowerAntiTriangular;
        return;
    }
    if(!thereAreLowerAntiDiagonalBlocks && thereAreUpperAntiDiagonalBlocks)
    {
        M_structure = UpperAntiTriangular;
        return;
    }

    M_structure = NoStructure;
}

int BlockOperator::SetUseTranspose(bool useTranspose)
{
    M_useTranspose = useTranspose;
#ifdef HAVE_LIFEV_DEBUG
    // Checking that all the blocks support the transpose option
    for (UInt i(0); i<M_nBlockRows; ++i)
        for (UInt j(0); j<M_nBlockCols; ++j)
        {
            EPETRA_CHK_ERR(M_oper(i,j)->SetUseTranspose(true));
            EPETRA_CHK_ERR(M_oper(i,j)->SetUseTranspose(false));
        }
#endif
    return 0;
}

int BlockOperator::Apply(const vector_Type & X, vector_Type & Y) const
{
    int error(-1);
    if (M_useTranspose)
        error = applyTranspose(X,Y);
    else
        error = applyNoTranspose(X,Y);
    return error;
}

int BlockOperator::ApplyInverse(const vector_Type & X, vector_Type & Y) const
{
    ASSERT_PRE(M_nBlockRows == M_nBlockCols, "The operator must be squared");
    ASSERT_PRE(M_domainMap->monolithicMap()->SameAs(*M_rangeMap->monolithicMap() ), "The operator must be squared");
    ASSERT_PRE(Y.Map().SameAs(*M_domainMap->monolithicMap()),"The map of Y is not conforming with domain map.");
    ASSERT_PRE(X.Map().SameAs(*M_rangeMap->monolithicMap()), "The map of X is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    switch(M_structure)
    {
    case Diagonal:
    case AntiDiagonal:
        return blockJacobi(X,Y);
        break;
    case LowerTriangular:
    case LowerAntiTriangular:
        return blockLowerTriangularSolve(X,Y);
        break;
    case UpperTriangular:
    case UpperAntiTriangular:
        return blockUpperTriangularSolve(X,Y);
        break;
    case NoStructure:
        Y.Scale(1.0/0.0);
        return -1;
        break;
    case Rectangular:
        Y.Scale(1.0/0.0);
        return -1;
        break;
    default:
        Y.Scale(1.0/0.0);
        return -1;
        break;
    }
    return -1;
}

const BlockOperator::operatorPtr_Type& BlockOperator::block(UInt iblock, UInt jblock) const
{
    ASSERT (iblock<M_nBlockRows,"Error! Index out of bounds.\n");
    ASSERT (jblock<M_nBlockCols,"Error! Index out of bounds.\n");
    return M_oper(iblock,jblock);
}

//===========================================================================//
//===========================================================================//
//  Private Methods                                                          //
//===========================================================================//
//===========================================================================//

int BlockOperator::applyNoTranspose(const vector_Type & X, vector_Type & Y) const
{
    ASSERT_PRE(X.Map().SameAs(*(M_domainMap->monolithicMap())),"The map of X is not conforming with domain map.");
    ASSERT_PRE(Y.Map().SameAs(*(M_rangeMap->monolithicMap())), "The map of Y is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    const std::unique_ptr<BlockEpetra_MultiVector> Xview( createBlockView(X, *M_domainMap) );
    const std::unique_ptr<BlockEpetra_MultiVector> Yview( createBlockView(Y, *M_rangeMap) );
    BlockEpetra_MultiVector tmpY(*M_rangeMap, X.NumVectors(), true);

    Yview->PutScalar(0.0);


    // Perform the mat-vec multiplications
    for(UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        for(UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        {
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->Apply(Xview->block(jblock), tmpY.block(iblock) ));
            EPETRA_CHK_ERR(Yview->block(iblock).Update(1.0, tmpY.block(iblock), 1.0));
        }

    return 0;
}

int BlockOperator::applyTranspose(const vector_Type & X, vector_Type & Y) const
{
    ASSERT_PRE(X.Map().SameAs(*(M_rangeMap->monolithicMap())),"The map of X is not conforming with domain map.");
    ASSERT_PRE(Y.Map().SameAs(*(M_domainMap->monolithicMap())), "The map of Y is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    const std::unique_ptr<BlockEpetra_MultiVector> Xview( createBlockView(X, *M_rangeMap) );
    const std::unique_ptr<BlockEpetra_MultiVector> Yview( createBlockView(Y, *M_domainMap) );
    BlockEpetra_MultiVector tmpY(*M_domainMap, X.NumVectors(), true);

    Yview->PutScalar(0.0);

    // Perform the mat-vec multiplications
    for(UInt iblock=0; iblock < M_nBlockCols; ++iblock)
        for(UInt jblock=0; jblock < M_nBlockRows; ++jblock)
        {
            EPETRA_CHK_ERR(M_oper(iblock,jblock)->SetUseTranspose(true));
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->Apply(Xview->block(jblock), tmpY.block(iblock) ));
            EPETRA_CHK_ERR(Yview->block(iblock).Update(1.0, tmpY.block(iblock), 1.0));
            EPETRA_CHK_ERR(M_oper(iblock,jblock)->SetUseTranspose(false));
        }

    return 0;
}

int BlockOperator::blockJacobi(const vector_Type & X, vector_Type & Y) const
{
    const std::unique_ptr<BlockEpetra_MultiVector> Xcopy( new BlockEpetra_MultiVector(Copy, X, *M_rangeMap) );
    const std::unique_ptr<BlockEpetra_MultiVector> Yview( new BlockEpetra_MultiVector(View, Y, *M_rangeMap) );

    Yview->PutScalar(0.0);

    if (M_structure==Diagonal)
    {
        // Diagonal case
        for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
            EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(Xcopy->block(iblock), Yview->block(iblock) ));
    }
    else
    {
        // Anti Diagonal case
        for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        {
            UInt jblock = M_nBlockCols-iblock-1;
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->ApplyInverse(Xcopy->block(jblock), Yview->block(iblock) ));
        }
    }
    return 0;
}

int BlockOperator::blockLowerTriangularSolve(const vector_Type & X, vector_Type & Y) const
{
    const std::unique_ptr<BlockEpetra_MultiVector> Xcopy( new BlockEpetra_MultiVector(Copy, X, *M_rangeMap) );
    const std::unique_ptr<BlockEpetra_MultiVector> Yview( new BlockEpetra_MultiVector(View, Y, *M_rangeMap) );
    BlockEpetra_MultiVector Z(*M_rangeMap, X.NumVectors(), true);

    Yview->PutScalar(0.0);
    Z.PutScalar(0.0);

    if (M_structure==LowerTriangular)
    {
        // Lower Triangular case
        for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        {
            EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(Xcopy->block(iblock), Yview->block(iblock) ));
            for(UInt kblock = iblock+1; kblock < M_nBlockRows; ++kblock)
            {
                EPETRA_CHK_ERR( M_oper(kblock, iblock)->Apply(Yview->block(iblock), Z.block(kblock) ) );
                EPETRA_CHK_ERR( Xcopy->block(kblock).Update(-1.0, Z.block(kblock), 1.0) );
            }
        }
    }
    else
    {
        // Lower Triangular with respect to the anti-diagonal case
        for (UInt iblock = 0; iblock<M_nBlockRows; ++iblock)
        {
            UInt jblock = M_nBlockCols-iblock-1;
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->ApplyInverse(Xcopy->block(iblock), Yview->block(jblock) ));
            for(UInt kblock = iblock+1; kblock<M_nBlockRows; ++kblock)
            {
                EPETRA_CHK_ERR( M_oper(kblock, jblock)->Apply(Yview->block(jblock), Z.block(kblock) ) );
                EPETRA_CHK_ERR( Xcopy->block(kblock).Update(-1.0, Z.block(kblock), 1.0) );
            }
        }
    }

    return 0;
}

int BlockOperator::blockUpperTriangularSolve(const vector_Type & X, vector_Type & Y) const
{
    const std::unique_ptr<BlockEpetra_MultiVector> Xcopy( new BlockEpetra_MultiVector(Copy, X, *M_rangeMap) );
    const std::unique_ptr<BlockEpetra_MultiVector> Yview( new BlockEpetra_MultiVector(View, Y, *M_rangeMap) );
    BlockEpetra_MultiVector Z(*M_rangeMap, X.NumVectors(), true);

    Yview->PutScalar(0.0);
    Z.PutScalar(0.0);

    if (M_structure==UpperTriangular)
    {
        // Upper Triangular case
        for(int iblock = M_nBlockRows - 1 ; iblock > -1 ; --iblock)
        {
            EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(Xcopy->block(iblock), Yview->block(iblock) ));
            for(int kblock = 0; kblock < iblock; ++kblock)
            {
                EPETRA_CHK_ERR( M_oper(kblock, iblock)->Apply(Yview->block(iblock), Z.block(kblock) ) );
                EPETRA_CHK_ERR( Xcopy->block(kblock).Update(-1.0, Z.block(kblock), 1.0) );
            }
        }
    }
    else
    {
        // Upper Triangular with respect to the anti-diagonal case
        for(int iblock = M_nBlockRows - 1 ; iblock > -1 ; --iblock)
        {
            UInt jblock = M_nBlockCols-iblock-1;
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->ApplyInverse(Xcopy->block(iblock), Yview->block(jblock) ));
            for(int kblock = 0; kblock < iblock; ++kblock)
            {
                EPETRA_CHK_ERR( M_oper(kblock, jblock)->Apply(Yview->block(jblock), Z.block(kblock) ) );
                EPETRA_CHK_ERR( Xcopy->block(kblock).Update(-1.0, Z.block(kblock), 1.0) );
            }
        }
    }

    return 0;
}

} /* end namespace Operators */

} /*end namespace */
