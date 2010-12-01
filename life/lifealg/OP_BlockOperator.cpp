/*
 * OP_BlockOperator.cpp
 *
 *  Created on: Sep 20, 2010
 *      Author: uvilla
 */

#include <life/lifealg/OP_BlockOperator.hpp>

namespace LifeV
{
namespace Operators
{
BlockOperator::BlockOperator():
        M_name("BlockOperator"),
        M_useTranspose(false),
        M_structure(NoStructure)
{};

void BlockOperator::setUp(UInt nBlocks,
                          const std::vector<boost::shared_ptr<Epetra_Map> > domainMap,
                          const boost::shared_ptr<Epetra_Comm> & comm)
{
    M_comm = comm;

    M_nBlockRows = nBlocks;
    M_nBlockCols = nBlocks;
    M_domainBlockMaps = domainMap;

    M_domainBlockMapsShift.resize(nBlocks);
    for (UInt iblock=0; iblock < nBlocks; ++iblock)
        M_domainBlockMapsShift[iblock].reset(new Epetra_Map(*M_domainBlockMaps[iblock]));

    buildImporter(nBlocks, M_domainBlockMapsShift, M_domainMap,
                  M_block2monoDomain, M_mono2blockDomain);

    //Since the operator is square I'll just copy the pointers not the content of the pointers
    M_rangeBlockMaps = M_domainBlockMaps;
    M_rangeBlockMapsShift = M_domainBlockMapsShift;
    M_rangeMap = M_domainMap;
    M_block2monoRange = M_block2monoDomain;
    M_mono2blockRange = M_mono2blockDomain;

    M_oper.resize(nBlocks, nBlocks);
}


//! SetUp for a "rectangular operator"
void BlockOperator::setUp(UInt nRowBlocks, UInt nColBlocks,
                          const std::vector<boost::shared_ptr<Epetra_Map> > & domainMap,
                          const std::vector<boost::shared_ptr<Epetra_Map> > & rangeMap,
                          const boost::shared_ptr<Epetra_Comm> & comm)
{
    M_comm = comm;

    M_nBlockRows = nRowBlocks;
    M_nBlockCols = nColBlocks;

    M_rangeBlockMaps = rangeMap;
    M_domainBlockMaps = domainMap;

    M_rangeBlockMapsShift.resize(M_nBlockRows);
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        M_rangeBlockMapsShift[iblock].reset(new Epetra_Map(*M_rangeBlockMaps[iblock]));

    M_domainBlockMapsShift.resize(M_nBlockCols);
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        M_domainBlockMapsShift[jblock].reset(new Epetra_Map(*M_domainBlockMaps[jblock]));

    buildImporter(M_nBlockRows, M_rangeBlockMapsShift, M_rangeMap,
                  M_block2monoRange, M_mono2blockRange);

    buildImporter(M_nBlockCols, M_domainBlockMapsShift, M_domainMap,
                  M_block2monoDomain, M_mono2blockDomain);

    M_oper.resize(M_nBlockRows, M_nBlockCols);

}

//! SetUp when the operator is given like a boost::matrix
void BlockOperator::setUp(const BlockOper & blockOper, const boost::shared_ptr<Epetra_Comm> & comm)
{
    M_nBlockRows = blockOper.size1();
    M_nBlockCols = blockOper.size2();

    M_rangeBlockMaps.resize(M_nBlockRows);
    M_domainBlockMaps.resize(M_nBlockCols);

    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        M_rangeBlockMaps[iblock].reset(new Epetra_Map(blockOper(iblock,0)->OperatorRangeMap() ));


    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        M_domainBlockMaps[jblock].reset(new Epetra_Map(blockOper(0,jblock)->OperatorDomainMap() ));

    setUp(M_nBlockRows, M_nBlockCols, M_domainBlockMaps, M_rangeBlockMaps, comm);
    M_oper = blockOper;
    fillComplete();

}

//! set a component of the block operator
void BlockOperator::setBlock(UInt iblock, UInt jblock, const operator_ptr & operBlock)
{
    ASSERT_PRE(M_rangeBlockMaps[iblock]->PointSameAs(operBlock->OperatorRangeMap()), "Wrong range map");
    ASSERT_PRE(M_domainBlockMaps[jblock]->PointSameAs(operBlock->OperatorDomainMap()), "Wrong domain map");

    M_oper(iblock, jblock) = operBlock;
}


void BlockOperator::fillComplete()
{
//Check for Structure
    bool thereAreUpperDiagonalBlocks(false);
    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock = iblock+1; jblock < M_nBlockCols; ++jblock)
            if (M_oper(iblock,jblock).get() != 0 && M_oper(iblock, jblock)->HasNormInf() && M_oper(iblock, jblock)->NormInf()!=0 )
                thereAreUpperDiagonalBlocks = true;

    bool thereAreLowerDiagonalBlocks(false);
    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock = 0; jblock < iblock; ++jblock)
            if (M_oper(iblock,jblock).get() != 0  && M_oper(iblock, jblock)->HasNormInf() && M_oper(iblock, jblock)->NormInf()!=0 )
                thereAreLowerDiagonalBlocks = true;

    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock = 0; jblock < M_nBlockCols; ++jblock)
        {
            if (M_oper(iblock,jblock).get() == 0)
            {
                NullOperator * nullOp(new NullOperator);
                nullOp->setUp(M_domainBlockMaps[jblock], M_rangeBlockMaps[iblock]);
                M_oper(iblock,jblock).reset(nullOp);
            }
            ASSERT(M_rangeBlockMaps[iblock]->PointSameAs(M_oper(iblock,jblock)->OperatorRangeMap()), "Wrong range map");
            ASSERT(M_domainBlockMaps[jblock]->PointSameAs(M_oper(iblock,jblock)->OperatorDomainMap()), "Wrong domain map");
        }

    if (M_nBlockRows != M_nBlockCols)
    {
        M_structure = Rectangular;
        return;
    }

    if (!thereAreLowerDiagonalBlocks && !thereAreUpperDiagonalBlocks)
    {
        M_structure = Diagonal;
        return;
    }
    if (thereAreLowerDiagonalBlocks && !thereAreUpperDiagonalBlocks)
    {
        M_structure = LowerTriangular;
        return;
    }
    if (!thereAreLowerDiagonalBlocks && thereAreUpperDiagonalBlocks)
    {
        M_structure = UpperTriangular;
        return;
    }
    if (thereAreLowerDiagonalBlocks && thereAreUpperDiagonalBlocks)
    {
        M_structure = NoStructure;
        return;
    }


}


int BlockOperator::Apply(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{

    ASSERT_PRE(X.Map().SameAs(*M_domainMap),"The map of X is not conforming with domain map.");
    ASSERT_PRE(Y.Map().SameAs(*M_rangeMap), "The map of Y is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    int nMultiVectors( X.NumVectors() );
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Xblock(M_nBlockCols);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Yblock(M_nBlockRows);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > tmpYblock(M_nBlockRows);

    // Split the vector
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
    {
        Xblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMapsShift[jblock], nMultiVectors) );
        Xblock[jblock]->PutScalar(0.0);
        EPETRA_CHK_ERR( Xblock[jblock]->Import(X, *M_mono2blockDomain[jblock], Insert) );
        EPETRA_CHK_ERR(Xblock[jblock]->ReplaceMap(*M_domainBlockMaps[jblock]) );
    }

    // Allocate Space for the Solution
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
    {
        Yblock[iblock].reset( new Epetra_MultiVector(*M_rangeBlockMaps[iblock], nMultiVectors) );
        Yblock[iblock]->PutScalar(0.0);
        tmpYblock[iblock].reset( new Epetra_MultiVector(*M_rangeBlockMaps[iblock], nMultiVectors) );
    }

    //Perform the mat-vec multiplications
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
        for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
        {
            EPETRA_CHK_ERR(M_oper(iblock, jblock)->Apply(*Xblock[jblock], *tmpYblock[iblock]));
            EPETRA_CHK_ERR(Yblock[iblock]->Update(1.0, *tmpYblock[iblock], 1.0));
        }

    //Reassemble the results
    Y.PutScalar(0.0);
    for (UInt iblock=0; iblock<M_nBlockRows; ++iblock)
    {
        EPETRA_CHK_ERR(Yblock[iblock]->ReplaceMap(*M_rangeBlockMapsShift[iblock]) );
        EPETRA_CHK_ERR(Y.Import(*Yblock[iblock], *M_block2monoRange[iblock], Insert) );
    }
    return 0;

}

int BlockOperator::ApplyInverse(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{
    switch (M_structure)
    {
    case Diagonal:
        return blockJacobi(X,Y);
        break;
    case LowerTriangular:
        return blockLowerTriangularSolve(X,Y);
        break;
    case UpperTriangular:
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
}

//Merge two vectors using the domain map
int BlockOperator::merge(const Epetra_MultiVector & vBlock, Epetra_MultiVector & vMono, UInt jblock) const
{

    ASSERT_PRE(vBlock.Map().SameAs(*M_domainBlockMaps[jblock]), "vBlock does not have the correct map" );
    ASSERT_PRE(vMono.Map().SameAs(*M_domainMap), "vMono does not have the correct map");
    ASSERT_PRE(vBlock.NumVectors() == vMono.NumVectors(), "The number of vectors in vBlock and vMono is different" );

    Epetra_MultiVector tmp(vBlock);
    EPETRA_CHK_ERR( tmp.ReplaceMap(*M_domainBlockMapsShift[jblock]) );

    return vMono.Import(tmp,*M_block2monoDomain[jblock], Insert);

}

//Extract vectors using the range map
int BlockOperator::extract(Epetra_MultiVector & vBlock, const Epetra_MultiVector & vMono, UInt iblock) const
{
    ASSERT_PRE(M_rangeBlockMaps[iblock]->SameAs(vBlock.Map()), "vBlock does not have the correct map");
    ASSERT_PRE(M_rangeMap->SameAs(vMono.Map()), "vMono does not have the correct map");
    ASSERT_PRE(vBlock.NumVectors() == vMono.NumVectors(), "The number of vectors in vBlock and vMono is different" );

    EPETRA_CHK_ERR( vBlock.ReplaceMap(*M_rangeBlockMapsShift[iblock]) );
    EPETRA_CHK_ERR( vBlock.Import(vMono, *M_mono2blockDomain[iblock], Insert) );
    EPETRA_CHK_ERR( vBlock.ReplaceMap(*M_rangeBlockMaps[iblock]) );
    return 0;
}

int BlockOperator::split(const Epetra_MultiVector & up,
                         vector_container & vi) const
{
    UInt block(0);
    for (vector_container::iterator it=vi.begin(); it != vi.end(); ++it, ++block)
        EPETRA_CHK_ERR( extract(**it, up, block) );
    return 0;
}

int BlockOperator::merge( Epetra_MultiVector & up, const vector_container & vi) const
{
    UInt block(0);
    for (vector_container::const_iterator it=vi.begin(); it != vi.end(); ++it, ++block)
        EPETRA_CHK_ERR( merge(**it, up, block) );
    return 0;

}

//===========================================================================//
//===========================================================================//
//  Private Methods                                                          //
//===========================================================================//
//===========================================================================//
// Private Functions
void BlockOperator::buildImporter(UInt nblock,
                                  std::vector<boost::shared_ptr<Epetra_Map> > & blockMaps,
                                  boost::shared_ptr<Epetra_Map> & fullMap,
                                  std::vector< boost::shared_ptr<Epetra_Import> > & block2mono,
                                  std::vector< boost::shared_ptr<Epetra_Import> > & mono2block)
{

    std::vector<int> numLocIdBlock(nblock), numGlobalElBlock(nblock);
    int numLocIdMono(0), numGlobalElMono(0);

    for (UInt iblock=0; iblock < nblock; ++iblock)
    {
        numLocIdBlock[iblock] = blockMaps[iblock]->NumMyElements();
        numGlobalElBlock[iblock] = blockMaps[iblock]->NumGlobalElements();
    }

    numLocIdMono = std::accumulate(numLocIdBlock.begin(), numLocIdBlock.end(), numLocIdMono);
    numGlobalElMono = std::accumulate(numGlobalElBlock.begin(), numGlobalElBlock.end(), numGlobalElMono);

    std::vector<int> shift(nblock+1);
    std::vector<int>::iterator itShift(shift.begin());
    *itShift = 0;
    itShift++;
    std::partial_sum(numGlobalElBlock.begin(), numGlobalElBlock.end(), itShift);

    std::vector<std::vector<int> > myGlobIdBlock(nblock);
    std::vector<int> myGlobIdMono(numGlobalElMono);

    std::vector<int>::iterator itMono;
    itMono = myGlobIdMono.begin();

    for (std::vector<std::vector<int> >::iterator myGlobIdBlockIt = myGlobIdBlock.begin();
            myGlobIdBlockIt != myGlobIdBlock.end();
            ++myGlobIdBlockIt)
    {
        UInt iblock = (UInt) (myGlobIdBlockIt - myGlobIdBlock.begin());
        myGlobIdBlockIt->resize(numLocIdBlock[iblock]);
        blockMaps[iblock]->MyGlobalElements( &(*myGlobIdBlockIt)[0] );

        for (std::vector<int>::iterator it=myGlobIdBlockIt->begin(); it!=myGlobIdBlockIt->end(); ++it)
            *it += shift[iblock];

        itMono = std::copy(myGlobIdBlockIt->begin(), myGlobIdBlockIt->end(), itMono);

        blockMaps[iblock].reset(new Epetra_Map(numGlobalElBlock[iblock], myGlobIdBlockIt->size(),
                                               &(*myGlobIdBlockIt)[0], shift[iblock]+1, *M_comm));

    }

    fullMap.reset(new Epetra_Map(numGlobalElMono, myGlobIdMono.size(), &myGlobIdMono[0], 1, *M_comm));

    block2mono.resize(nblock);
    mono2block.resize(nblock);

    for (UInt iblock=0; iblock<nblock; ++iblock)
    {
        block2mono[iblock].reset(new Epetra_Import(*fullMap, *blockMaps[iblock]));
        mono2block[iblock].reset(new Epetra_Import(*blockMaps[iblock], *fullMap));
    }

}


int BlockOperator::blockJacobi(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{

    ASSERT_PRE(M_nBlockRows == M_nBlockCols, "The operator must be squared");
    ASSERT_PRE(M_domainMap->SameAs(*M_rangeMap), "The operator must be squared");
    ASSERT_PRE(Y.Map().SameAs(*M_domainMap),"The map of Y is not conforming with domain map.");
    ASSERT_PRE(X.Map().SameAs(*M_rangeMap), "The map of X is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    int nMultiVectors( X.NumVectors() );

    std::vector< boost::shared_ptr<Epetra_MultiVector> > Xblock(M_nBlockCols);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Yblock(M_nBlockRows);

    // Split the vector
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
    {
        Xblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMapsShift[jblock], nMultiVectors) );
        Xblock[jblock]->PutScalar(0.0);
        EPETRA_CHK_ERR(Xblock[jblock]->Import(X, *M_mono2blockDomain[jblock], Insert) );
        EPETRA_CHK_ERR(Xblock[jblock]->ReplaceMap(*M_domainBlockMaps[jblock]) );
    }

    // Allocate Space for the Solution
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
    {
        Yblock[iblock].reset( new Epetra_MultiVector(*M_rangeBlockMaps[iblock], nMultiVectors) );
        Yblock[iblock]->PutScalar(0.0);
    }

    //Perform the mat-vec multiplications
    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
        EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(*Xblock[iblock], *Yblock[iblock]));

    //Reassemble the results
    Y.PutScalar(0.0);
    for (UInt iblock=0; iblock<M_nBlockRows; ++iblock)
    {
        EPETRA_CHK_ERR(Yblock[iblock]->ReplaceMap(*M_rangeBlockMapsShift[iblock]) );
        EPETRA_CHK_ERR( Y.Import(*Yblock[iblock], *M_block2monoRange[iblock], Insert) );
    }

    return 0;
}

int BlockOperator::blockLowerTriangularSolve(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{

    ASSERT_PRE(M_nBlockRows == M_nBlockCols, "The operator must be squared");
    ASSERT_PRE(M_domainMap->SameAs(*M_rangeMap), "The operator must be squared");
    ASSERT_PRE(Y.Map().SameAs(*M_domainMap),"The map of Y is not conforming with domain map.");
    ASSERT_PRE(X.Map().SameAs(*M_rangeMap), "The map of X is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    int nMultiVectors( X.NumVectors() );

    std::vector< boost::shared_ptr<Epetra_MultiVector> > Xblock(M_nBlockCols);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Yblock(M_nBlockRows);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Zblock(M_nBlockCols);

    // Split the vector
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
    {
        Xblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMapsShift[jblock], nMultiVectors) );
        Xblock[jblock]->PutScalar(0.0);
        EPETRA_CHK_ERR(Xblock[jblock]->Import(X, *M_mono2blockDomain[jblock], Insert) );
        EPETRA_CHK_ERR(Xblock[jblock]->ReplaceMap(*M_domainBlockMaps[jblock]) );

        Zblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMaps[jblock], nMultiVectors) );
    }

    // Allocate Space for the Solution
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
    {
        Yblock[iblock].reset( new Epetra_MultiVector(*M_rangeBlockMaps[iblock], nMultiVectors) );
        Yblock[iblock]->PutScalar(0.0);
    }

    //Perform the mat-vec multiplications
    for (UInt iblock = 0; iblock < M_nBlockRows; ++iblock)
    {
        EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(*Xblock[iblock], *Yblock[iblock]));
        for (UInt kblock = iblock; kblock < M_nBlockRows; ++kblock)
        {
            EPETRA_CHK_ERR( M_oper(kblock, iblock)->Apply(*Yblock[iblock], *Zblock[kblock]) );
            Xblock[kblock]->Update(-1.0, *Zblock[kblock], 1.0);
        }
    }

    //Reassemble the results
    Y.PutScalar(0.0);
    for (UInt iblock=0; iblock<M_nBlockRows; ++iblock)
    {
        EPETRA_CHK_ERR(Yblock[iblock]->ReplaceMap(*M_rangeBlockMapsShift[iblock]) );
        EPETRA_CHK_ERR( Y.Import(*Yblock[iblock], *M_block2monoRange[iblock], Insert) );
    }

    return 0;
}

int BlockOperator::blockUpperTriangularSolve(const Epetra_MultiVector & X, Epetra_MultiVector & Y) const
{

    ASSERT_PRE(M_nBlockRows == M_nBlockCols, "The operator must be squared");
    ASSERT_PRE(M_domainMap->SameAs(*M_rangeMap), "The operator must be squared");
    ASSERT_PRE(Y.Map().SameAs(*M_domainMap),"The map of Y is not conforming with domain map.");
    ASSERT_PRE(X.Map().SameAs(*M_rangeMap), "The map of X is not conforming with range  map.");
    ASSERT_PRE(X.NumVectors() == Y.NumVectors(), "The number of vectors in X and Y is different" );

    int nMultiVectors( X.NumVectors() );

    std::vector< boost::shared_ptr<Epetra_MultiVector> > Xblock(M_nBlockCols);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Yblock(M_nBlockRows);
    std::vector< boost::shared_ptr<Epetra_MultiVector> > Zblock(M_nBlockCols);

    // Split the vector
    for (UInt jblock=0; jblock < M_nBlockCols; ++jblock)
    {
        Xblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMapsShift[jblock], nMultiVectors) );
        Xblock[jblock]->PutScalar(0.0);
        EPETRA_CHK_ERR(Xblock[jblock]->Import(X, *M_mono2blockDomain[jblock], Insert) );
        EPETRA_CHK_ERR(Xblock[jblock]->ReplaceMap(*M_domainBlockMaps[jblock]) );

        Zblock[jblock].reset( new Epetra_MultiVector(*M_domainBlockMaps[jblock], nMultiVectors) );
    }

    // Allocate Space for the Solution
    for (UInt iblock=0; iblock < M_nBlockRows; ++iblock)
    {
        Yblock[iblock].reset( new Epetra_MultiVector(*M_rangeBlockMaps[iblock], nMultiVectors) );
        Yblock[iblock]->PutScalar(0.0);
    }

    //Perform the mat-vec multiplications
    for (int iblock = M_nBlockRows - 1 ; iblock > -1 ; --iblock)
    {
        EPETRA_CHK_ERR(M_oper(iblock, iblock)->ApplyInverse(*Xblock[iblock], *Yblock[iblock]));
        for (int kblock = 0; kblock < iblock; ++kblock)
        {
            EPETRA_CHK_ERR( M_oper(kblock, iblock)->Apply(*Yblock[iblock], *Zblock[kblock]) );
            Xblock[kblock]->Update(-1.0, *Zblock[kblock], 1.0);
        }
    }

    //Reassemble the results
    Y.PutScalar(0.0);
    for (UInt iblock=0; iblock<M_nBlockRows; ++iblock)
    {
        EPETRA_CHK_ERR(Yblock[iblock]->ReplaceMap(*M_rangeBlockMapsShift[iblock]) );
        EPETRA_CHK_ERR(Y.Import(*Yblock[iblock], *M_block2monoRange[iblock], Insert) );
    }

    return 0;
}


} /* end namespace Operators */
} /*end namespace */
