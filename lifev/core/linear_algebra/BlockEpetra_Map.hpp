/*
 * BlockEpetra_Map.hpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */

//@ HEADER

/*!
 * \file BlockEpetra_Map.hpp
 * \author Umberto Villa
 * \date 2011-08-13
 * A special map object to handle parallel block structured Vectors.
 */

#ifndef BLOCKEPETRA_MAP_HPP_
#define BLOCKEPETRA_MAP_HPP_

#include <lifev/operator/linear_algebra/LinearOperator.hpp>
#include <Epetra_Import.h>


namespace LifeV
{

class BlockEpetra_MultiVector;

//! @class
/*! @brief This class handles block access to parallel monolithic Vectors with an underling block structure.
 *
 * The goal of BlockEpetra_Map is to provide block or monolithic access to block structured parallel Vectors.
 * The BlockEpetra_Map takes as input a list of Epetra_BlockMap and construct a monolithic map by striding
 * them together.
 */
class BlockEpetra_Map
{
public:

	//@name Public Typedefs
	//@{
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
	//@}

	//@name Constructors
	//@{
	//! Default Constructor
	BlockEpetra_Map();
	//! Copy constructor
	BlockEpetra_Map(const BlockEpetra_Map & map);
	//! Construct a BlockEpetra_Map from an ordered list of Epetra_BlockMap objects
	BlockEpetra_Map(const mapPtrContainer_Type & mapPtrContainer);
	//@}

	//@name SetUp
	//@{
	//! simplified setup for the case of only 2 blocks
	void setUp(const mapPtr_Type & map1, const mapPtr_Type & map2);
	//! simplified setup for the case of only 3 blocks
	void setUp(const mapPtr_Type & map1, const mapPtr_Type & map2, const mapPtr_Type & map3);
	//! general setup for arbitrary number of blocks. The maps of the block to stride are collected in a mapPtrContainer_Type object
	void setUp(const mapPtrContainer_Type & mapPtrContainer);
	//! other set up routines to overcome a design problem in Trilinos (2 blocks)
	void setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2);
	//! other set up routines to overcome a design problem in Trilinos (3 blocks)
	void setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2, const Epetra_BlockMap & map3);
	//@}

	//@name Getters
	//@{
	//! get the number of blocks
	UInt nBlocks() const;

	//! returns the monolithicMap
	mapPtr_Type & monolithicMap();
	//! const version
	const constMapPtr_Type monolithicMap() const;

	//! returns the original map of block iblock
	mapPtr_Type & blockMap(UInt iblock);
	//! const version
	const constMapPtr_Type blockMap(UInt iblock) const;

	//! returns the map relative to block iblock in the numbering of the monolithic map
	mapPtr_Type & blockShiftedMap(UInt iblock);
	//! const version
	const constMapPtr_Type blockShiftedMap(UInt iblock) const;

	//! return an Trilinos import_Type object from block \c iblock map to the monolithic map
	const import_Type & block2monoImporter(UInt iblock) const;
	//! return an Trilinos import_Type object from the monolithic map to block \c iblock map.
	const import_Type & mono2blockImporter(UInt iblock) const;
	//@}

	//! Print debug info
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

//! @func
/*!
 * Utility function to transform a Epetra_BlockMap in a Epetra_Map.
 * Used to overcome a bad design in Trilinos.
 */
Epetra_Map * blockMap2Map(const Epetra_BlockMap * blockMap);

}
#endif /* BLOCKEPETRA_MAP_HPP_ */
