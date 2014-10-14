/*
 * BlockEpetra_MultiVector.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */

#include <lifev/core/linear_algebra/BlockEpetra_MultiVector.hpp>

namespace LifeV
{
BlockEpetra_MultiVector::BlockEpetra_MultiVector(const BlockEpetra_Map& map, int numVectors, bool zeroOut):
		Epetra_MultiVector(*map.monolithicMap(), numVectors, zeroOut),
		M_nBlocks(map.nBlocks()),
		M_blockMap(map),
		M_myLocalOffsets(map.M_myLocalOffsets),
		M_blocks(map.nBlocks())
{
	createBlockViews();
}

BlockEpetra_MultiVector::BlockEpetra_MultiVector(const BlockEpetra_MultiVector& source):
		Epetra_MultiVector(source),
		M_nBlocks(source.M_nBlocks),
		M_blockMap(source.M_blockMap),
		M_myLocalOffsets(source.M_myLocalOffsets),
		M_blocks(source.M_blocks)
{
	createBlockViews();
}

BlockEpetra_MultiVector::BlockEpetra_MultiVector(Epetra_DataAccess CV, const vector_Type & source, const BlockEpetra_Map& map):
		Epetra_MultiVector(CV, source, 0, source.NumVectors()),
		M_nBlocks(map.nBlocks()),
		M_blockMap(map),
		M_myLocalOffsets(map.M_myLocalOffsets),
		M_blocks(map.nBlocks())
{
	createBlockViews();
}


Epetra_MultiVector & BlockEpetra_MultiVector::block(UInt iblock)
{
	ASSERT_PRE(iblock < M_nBlocks, "Asking a block out of range!");

	double ** arrayOfPointers;
	ExtractView(&arrayOfPointers);
	double ** arrayOfShiftedPointers = new double*[NumVectors()];

	for(Int iVect(0); iVect < NumVectors(); ++iVect)
		arrayOfShiftedPointers[iVect] = arrayOfPointers[iVect] + M_myLocalOffsets[iblock];

	M_blocks[iblock]->ResetView(arrayOfShiftedPointers);

	delete[] arrayOfShiftedPointers;

	return *(M_blocks[iblock]);
}

const Epetra_MultiVector & BlockEpetra_MultiVector::block(UInt iblock) const
{
	ASSERT_PRE(iblock < M_nBlocks, "Asking a block out of range!");

	double ** arrayOfPointers;
	ExtractView(&arrayOfPointers);
	double ** arrayOfShiftedPointers = new double*[NumVectors()];

	for(Int iVect(0); iVect < NumVectors(); ++iVect)
		arrayOfShiftedPointers[iVect] = arrayOfPointers[iVect] + M_myLocalOffsets[iblock];

	M_blocks[iblock]->ResetView(arrayOfShiftedPointers);

	delete[] arrayOfShiftedPointers;

	return *(M_blocks[iblock]);
}

const BlockEpetra_Map & BlockEpetra_MultiVector::blockEpetraMap() const
{
	return M_blockMap;
}

//-----------------------------------------------------/
void BlockEpetra_MultiVector::createBlockViews()
{
	double ** arrayOfPointers;
	ExtractView(&arrayOfPointers);
	double ** arrayOfShiftedPointers = new double*[NumVectors()];
	for(UInt iblock(0); iblock < M_nBlocks; ++iblock)
	{
		for(Int iVect(0); iVect < NumVectors(); ++iVect)
			arrayOfShiftedPointers[iVect] = arrayOfPointers[iVect] + M_myLocalOffsets[iblock];
		M_blocks[iblock].reset(
				new Epetra_MultiVector(View, *(M_blockMap.blockMap(iblock)), arrayOfShiftedPointers, NumVectors())
				);
	}
	delete[] arrayOfShiftedPointers;
}

//------------------------------------------------------/
// Generate a BlockEpetra_MultiVector from a list of Epetra_MultiVector.
BlockEpetra_MultiVector * stride(std::vector<const BlockEpetra_MultiVector::vector_Type *> vectors)
{
	//Step 0: Check input
	const UInt nBlocks(vectors.size());
	const UInt numVectors(vectors[0]->NumVectors());

	for(UInt i(0); i<nBlocks; ++i)
		ASSERT_PRE( numVectors == static_cast<UInt>(vectors[i]->NumVectors()), "error in BlockEpetra_MultiVector * stride\n");

	//Step 1: Create the block map
	BlockEpetra_Map blockMap;

	BlockEpetra_Map::mapPtrContainer_Type mapPtrContainer(nBlocks);
	for(UInt i(0); i < nBlocks; ++i)
		mapPtrContainer[i].reset( blockMap2Map(&(vectors[i]->Map())) );

	blockMap.setUp(mapPtrContainer);

	//Step 2: Allocate the BlockEpetra_MultiVector
	BlockEpetra_MultiVector * stridedVector;
	stridedVector = new BlockEpetra_MultiVector(blockMap, numVectors, false);

	stridedVector->PutScalar(404.3142713);

	//Step 3: Copy the contents of the vectors
	for(UInt i(0); i < nBlocks; ++i)
		stridedVector->block(i) = *(vectors[i]);

	//Step 4: return the stridedVector
	return stridedVector;
}

// Generate a BlockEpetra_MultiVector from two Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
								 const BlockEpetra_MultiVector::vector_Type & v2)
{
	std::vector<const BlockEpetra_MultiVector::vector_Type *> vectors(2);
	vectors[0] = &v1;
	vectors[1] = &v2;
	return stride(vectors);
}

//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
		                         const BlockEpetra_MultiVector::vector_Type & v2,
		                         const BlockEpetra_MultiVector::vector_Type & v3)
{
	std::vector<const BlockEpetra_MultiVector::vector_Type *> vectors(3);
	vectors[0] = &v1;
	vectors[1] = &v2;
	vectors[2] = &v3;
	return stride(vectors);
}

//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
		                         const BlockEpetra_MultiVector::vector_Type & v2,
		                         const BlockEpetra_MultiVector::vector_Type & v3,
		                         const BlockEpetra_MultiVector::vector_Type & v4)
{
	std::vector<const BlockEpetra_MultiVector::vector_Type *> vectors(4);
	vectors[0] = &v1;
	vectors[1] = &v2;
	vectors[2] = &v3;
	vectors[3] = &v4;
	return stride(vectors);
}

//! Generate a BlockEpetra_MultiVector from a Epetra_MultiVector and a BlockEpetra_Map
BlockEpetra_MultiVector * createBlockView(const BlockEpetra_MultiVector::vector_Type & source, const BlockEpetra_Map& map)
{
	BlockEpetra_MultiVector * bMV;
	bMV = new BlockEpetra_MultiVector(View, source, map);
	return bMV;
}



} /* end namespace */
