/*
 * BlockOperator.cpp
 *
 *  Created on: Sep 20, 2010
 *      Author: dforti
 */

#include <lifev/fsi_blocks/solver/FSIApplyOperator.hpp>

namespace LifeV
{
namespace Operators
{
FSIApplyOperator::FSIApplyOperator():
M_name("FSIApplyOperator"),
M_useTranspose(false)
{

}

void FSIApplyOperator::setUp(const boost::shared_ptr<BlockEpetra_Map> & map, const commPtr_Type & comm)
{
    M_comm = comm;

    M_nBlockRows = map->nBlocks();
    M_nBlockCols = M_nBlockRows;
    M_domainMap = map;
    M_rangeMap = M_domainMap;
    M_oper.resize(M_nBlockRows, M_nBlockCols);
}


//! SetUp for a "rectangular operator"
void FSIApplyOperator::setUp(const boost::shared_ptr<BlockEpetra_Map> & domainMap,
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
void FSIApplyOperator::setUp(const operatorPtrContainer_Type & blockOper, const commPtr_Type & comm)
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
void FSIApplyOperator::setBlock(UInt iblock, UInt jblock, const operatorPtr_Type & operBlock)
{
    ASSERT_PRE(M_rangeMap->blockMap(iblock)->PointSameAs(operBlock->OperatorRangeMap()), "Wrong range map");
    ASSERT_PRE(M_domainMap->blockMap(jblock)->PointSameAs(operBlock->OperatorDomainMap()), "Wrong domain map");

    M_oper(iblock, jblock) = operBlock;
}


void FSIApplyOperator::fillComplete()
{
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
}

int FSIApplyOperator::SetUseTranspose(bool useTranspose)
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

int FSIApplyOperator::Apply(const vector_Type & X, vector_Type & Y) const
{
    int error(-1);
    if (M_useTranspose)
        error = applyTranspose(X,Y);
    else
    	error = applyNoTranspose(X,Y);
    return error;
}

int FSIApplyOperator::ApplyInverse(const vector_Type & X, vector_Type & Y) const
{
    // We just have to code the Apply method.
    return -1;
}

const FSIApplyOperator::operatorPtr_Type& FSIApplyOperator::block(UInt iblock, UInt jblock) const
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

int FSIApplyOperator::applyNoTranspose(const vector_Type & X, vector_Type & Y) const
{
    ASSERT_PRE(X.Map().SameAs(*(M_domainMap->monolithicMap())),"The map of X is not conforming with domain map.");
    ASSERT_PRE(Y.Map().SameAs(*(M_rangeMap->monolithicMap())), "The map of Y is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    const std::auto_ptr<BlockEpetra_MultiVector> Xview( createBlockView(X, *M_domainMap) );
    const std::auto_ptr<BlockEpetra_MultiVector> Yview( createBlockView(Y, *M_rangeMap) );
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

int FSIApplyOperator::applyTranspose(const vector_Type & X, vector_Type & Y) const
{
    ASSERT_PRE(X.Map().SameAs(*(M_rangeMap->monolithicMap())),"The map of X is not conforming with domain map.");
    ASSERT_PRE(Y.Map().SameAs(*(M_domainMap->monolithicMap())), "The map of Y is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    const std::auto_ptr<BlockEpetra_MultiVector> Xview( createBlockView(X, *M_rangeMap) );
    const std::auto_ptr<BlockEpetra_MultiVector> Yview( createBlockView(Y, *M_domainMap) );
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

} /* end namespace Operators */

} /*end namespace */
