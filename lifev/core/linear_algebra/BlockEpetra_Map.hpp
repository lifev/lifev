/*
 * BlockEpetra_Map.hpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */

#ifndef BLOCKEPETRA_MAP_HPP_
#define BLOCKEPETRA_MAP_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>
#include <Epetra_Import.h>


namespace LifeV
{

class BlockEpetra_MultiVector;

class BlockEpetra_Map
{
public:

	typedef Operators::LinearOperator::comm_Type comm_Type;
	typedef Operators::LinearOperator::map_Type    map_Type;
	typedef Operators::LinearOperator::mapPtr_Type mapPtr_Type;
	typedef Operators::LinearOperator::constMapPtr_Type constMapPtr_Type;
	typedef std::vector<mapPtr_Type> mapPtrContainer_Type;
	typedef mapPtrContainer_Type::iterator        mapPtrIterator_Type;
	typedef mapPtrContainer_Type::const_iterator  mapPtrConstIterator_Type;

	typedef Epetra_Import import_Type;
	typedef boost::shared_ptr<import_Type> importPtr_Type;
	typedef std::vector<importPtr_Type>    importPtrContainer_Type;
	typedef importPtrContainer_Type::iterator       importPtrIterator_Type;
	typedef importPtrContainer_Type::const_iterator importPtrConstIterator_Type;


	BlockEpetra_Map();
	BlockEpetra_Map(const BlockEpetra_Map & map);
	BlockEpetra_Map(const mapPtrContainer_Type & mapPtrContainer);
	void setUp(const mapPtr_Type & map1, const mapPtr_Type & map2);
	void setUp(const mapPtr_Type & map1, const mapPtr_Type & map2, const mapPtr_Type & map3);
	void setUp(const mapPtrContainer_Type & mapPtrContainer);
	void setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2);
	void setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2, const Epetra_BlockMap & map3);


	UInt nBlocks() const;

	mapPtr_Type & monolithicMap();
	const constMapPtr_Type monolithicMap() const;

	mapPtr_Type & blockMap(UInt iblock);
	const constMapPtr_Type blockMap(UInt iblock) const;

	mapPtr_Type & blockShiftedMap(UInt iblock);
	const constMapPtr_Type blockShiftedMap(UInt iblock) const;

	const import_Type & block2monoImporter(UInt iblock) const;
	const import_Type & mono2blockImporter(UInt iblock) const;

	void showMe();

	friend class BlockEpetra_MultiVector;



private:

	void build();

	UInt M_nBlocks;
	std::vector<UInt> M_myLocalOffsets;
	mapPtr_Type M_monolithicMap;
	mapPtrContainer_Type M_blockMap;
	mapPtrContainer_Type M_blockShiftedMap;

	//! @name Importers
	//@{
	//! merge block vectors in the monolithic vector
	importPtrContainer_Type M_block2mono;
	//! split the monolithic vector in the block vectors
	importPtrContainer_Type M_mono2block;
	//@}

};

Epetra_Map * blockMap2Map(const Epetra_BlockMap * blockMap);

}
#endif /* BLOCKEPETRA_MAP_HPP_ */
