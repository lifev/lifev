/*
 * BlockEpetra_MultiVector.hpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */

#ifndef BLOCKEPETRA_MULTIVECTOR_HPP_
#define BLOCKEPETRA_MULTIVECTOR_HPP_

#include <Epetra_MultiVector.h>
#include <lifev/operator/linear_algebra/BlockEpetra_Map.hpp>

namespace LifeV
{
class BlockEpetra_MultiVector : public Epetra_MultiVector
{

public:

	typedef Epetra_MultiVector vector_Type;
	typedef boost::shared_ptr<vector_Type> vectorPtr_Type;
	typedef std::vector<vectorPtr_Type>    vectorPtrContainer_Type;

	//! Generate a BlockEpetra_MultiVector from a BlockEpetra_Map
	BlockEpetra_MultiVector(const BlockEpetra_Map& map, int numVectors, bool zeroOut = true);
	//! Copy Constructor
	BlockEpetra_MultiVector(const BlockEpetra_MultiVector& Source);
	//! Cast
	BlockEpetra_MultiVector(Epetra_DataAccess CV, const vector_Type & source, const BlockEpetra_Map& map);

	Epetra_MultiVector & block(UInt iblock);
	const Epetra_MultiVector & block(UInt iblock) const;

	const BlockEpetra_Map & blockEpetraMap() const;

private:

	void createBlockViews();

	UInt M_nBlocks;
	BlockEpetra_Map M_blockMap;
	std::vector<UInt> M_myLocalOffsets;
	vectorPtrContainer_Type M_blocks;
};

//! Generate a BlockEpetra_MultiVector from a list of Epetra_MultiVector.
BlockEpetra_MultiVector * stride(std::vector<const BlockEpetra_MultiVector::vector_Type *> vector);
//! Generate a BlockEpetra_MultiVector from two Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
								 const BlockEpetra_MultiVector::vector_Type & v2);
//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
		                         const BlockEpetra_MultiVector::vector_Type & v2,
		                         const BlockEpetra_MultiVector::vector_Type & v3);
//! Generate a BlockEpetra_MultiVector from three Epetra_MultiVectors
BlockEpetra_MultiVector * stride(const BlockEpetra_MultiVector::vector_Type & v1,
		                         const BlockEpetra_MultiVector::vector_Type & v2,
		                         const BlockEpetra_MultiVector::vector_Type & v3,
		                         const BlockEpetra_MultiVector::vector_Type & v4);
//! Generate a BlockEpetra_MultiVector from a Epetra_MultiVector and a BlockEpetra_Map
BlockEpetra_MultiVector * createBlockView(const BlockEpetra_MultiVector::vector_Type & source, const BlockEpetra_Map& map);

} /* end namespace */


#endif /* BLOCKEPETRA_MULTIVECTOR_HPP_ */
