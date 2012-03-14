/*
 * BlockEpetra_Map.cpp
 *
 *  Created on: Aug 14, 2011
 *      Author: uvilla
 */

#include <lifev/operator/linear_algebra/BlockEpetra_Map.hpp>

namespace LifeV
{
BlockEpetra_Map::BlockEpetra_Map():
			M_nBlocks(0),
			M_myLocalOffsets(0),
			M_monolithicMap(),
			M_blockMap(0),
			M_blockShiftedMap(0),
			M_block2mono(0),
			M_mono2block(0)
			{ }

BlockEpetra_Map::BlockEpetra_Map(const BlockEpetra_Map & map):
		M_nBlocks(map.M_nBlocks),
		M_myLocalOffsets(map.M_myLocalOffsets),
		M_monolithicMap(new map_Type(*(map.M_monolithicMap))),
		M_blockMap(map.M_nBlocks),
		M_blockShiftedMap(map.M_nBlocks),
		M_block2mono(map.M_nBlocks),
		M_mono2block(map.M_nBlocks)
{
	for(UInt iblock(0); iblock < M_nBlocks; ++iblock)
	{
		M_blockMap[iblock].reset(new map_Type(*(map.M_blockMap[iblock])));
		M_blockShiftedMap[iblock].reset(new map_Type(*(map.M_blockShiftedMap[iblock])));
		M_mono2block[iblock].reset(new import_Type(*(map.M_mono2block[iblock])));
		M_block2mono[iblock].reset(new import_Type(*(map.M_block2mono[iblock])));
	}
}

BlockEpetra_Map::BlockEpetra_Map(const mapPtrContainer_Type & mapPtrContainer):
		M_nBlocks(0),
		M_myLocalOffsets(0),
		M_monolithicMap(),
		M_blockMap(0),
		M_blockShiftedMap(0),
		M_block2mono(0),
		M_mono2block(0)
{
	setUp(mapPtrContainer);
}

void BlockEpetra_Map::setUp(const mapPtr_Type & map1, const mapPtr_Type & map2)
{
	mapPtrContainer_Type tmpContainer(2);
	tmpContainer[0] = map1;
	tmpContainer[1] = map2;
	setUp(tmpContainer);
}

void BlockEpetra_Map::setUp(const mapPtr_Type & map1, const mapPtr_Type & map2, const mapPtr_Type & map3)
{
	mapPtrContainer_Type tmpContainer(3);
	tmpContainer[0] = map1;
	tmpContainer[1] = map2;
	tmpContainer[2] = map3;
	setUp(tmpContainer);
}

void BlockEpetra_Map::setUp(const mapPtrContainer_Type & mapPtrContainer)
{
	M_nBlocks = mapPtrContainer.size();
	M_blockMap.resize(M_nBlocks);
	M_myLocalOffsets.resize(M_nBlocks);
	M_blockShiftedMap.resize(M_nBlocks);
	M_mono2block.resize(M_nBlocks);
	M_block2mono.resize(M_nBlocks);

	for(UInt iblock(0); iblock<M_nBlocks; ++iblock)
	{
		M_blockMap[iblock].reset(new map_Type(*(mapPtrContainer[iblock])));
		M_blockShiftedMap[iblock].reset(new map_Type(*(mapPtrContainer[iblock])));
	}

	build();
}

void BlockEpetra_Map::setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2)
{
	mapPtrContainer_Type tmpContainer(2);
	tmpContainer[0].reset(blockMap2Map(&map1));
	tmpContainer[1].reset(blockMap2Map(&map2));
	setUp(tmpContainer);
}

void BlockEpetra_Map::setUp(const Epetra_BlockMap & map1, const Epetra_BlockMap & map2, const Epetra_BlockMap & map3)
{
	mapPtrContainer_Type tmpContainer(3);
	tmpContainer[0].reset(blockMap2Map(&map1));
	tmpContainer[1].reset(blockMap2Map(&map2));
	tmpContainer[2].reset(blockMap2Map(&map3));
	setUp(tmpContainer);
}


UInt BlockEpetra_Map::nBlocks() const
{
	return M_nBlocks;
}

BlockEpetra_Map::mapPtr_Type & BlockEpetra_Map::monolithicMap()
{
	return M_monolithicMap;
}
const BlockEpetra_Map::constMapPtr_Type BlockEpetra_Map::monolithicMap() const
{
	return M_monolithicMap;
}
BlockEpetra_Map::mapPtr_Type & BlockEpetra_Map::blockMap(UInt iblock)
{
	return M_blockMap[iblock];
}

const BlockEpetra_Map::constMapPtr_Type BlockEpetra_Map::blockMap(UInt iblock) const
{
	return M_blockMap[iblock];
}

BlockEpetra_Map::mapPtr_Type & BlockEpetra_Map::blockShiftedMap(UInt iblock)
{
	return M_blockShiftedMap[iblock];
}
const BlockEpetra_Map::constMapPtr_Type BlockEpetra_Map::blockShiftedMap(UInt iblock) const
{
	return M_blockShiftedMap[iblock];
}

const BlockEpetra_Map::import_Type & BlockEpetra_Map::block2monoImporter(UInt iblock) const
{
	return *M_block2mono[iblock];
}

const BlockEpetra_Map::import_Type & BlockEpetra_Map::mono2blockImporter(UInt iblock) const
{
	return *M_mono2block[iblock];
}

void BlockEpetra_Map::showMe()
{
	std::cout<<"Number of Blocks = " << M_nBlocks << "\n";
	for(UInt i(0); i < M_nBlocks; ++i)
		std::cout<<"MyLocalOffsets[" << i<< "] = " << M_myLocalOffsets[i] << "\n";
}

//-------------------------------------------------------//
void BlockEpetra_Map::build()
{

	const comm_Type & comm(M_blockMap[0]->Comm());

	// STEP 1
	// numMyElBlock(i) contains the number of element in map-i associated with my PID.
	// numGlobalElBlock(i) contains the global number of element in map-i
	// numMyElMono contains the number in the monolithic map associated with my PID.
	// numGlobalElMono contains the global number of element in the monolithic map
	std::vector<int> numMyElBlock(M_nBlocks), numGlobalElBlock(M_nBlocks);
	int numMyElMono(0), numGlobalElMono(0);

	for(UInt iblock=0; iblock < M_nBlocks; ++iblock)
	{
		numMyElBlock[iblock] = M_blockMap[iblock]->NumMyElements();
		numGlobalElBlock[iblock] = M_blockMap[iblock]->NumGlobalElements();
	}

	// I sum all over the blocks and I write the result in num*ElMono
	numMyElMono = std::accumulate(numMyElBlock.begin(), numMyElBlock.end(), numMyElMono);
	numGlobalElMono = std::accumulate(numGlobalElBlock.begin(), numGlobalElBlock.end(), numGlobalElMono);

    // STEP 2
    // Define the shift for each block
	std::vector<int> shift(M_nBlocks+1);
	std::vector<int>::iterator itShift(shift.begin());
	//the first shift is 0
	*itShift = 0; itShift++;
	// the other shifts are the sums of the global size of the previous blocks
	std::partial_sum(numGlobalElBlock.begin(), numGlobalElBlock.end(), itShift);

	// STEP 2 bis
	// Define the LocalShift for block in each processor
	std::vector<UInt>::iterator itMyShift(M_myLocalOffsets.begin());
	*itMyShift = 0; ++itMyShift;
	std::partial_sum(numMyElBlock.begin(), numMyElBlock.end()-1, itMyShift);


    // STEP 3 build the shifted and the monolithic maps
	// myGlobalElBlock(i)(j) will contain the Global Id in the shifted map-i corresponding to the local id j
	// myGlobalElMono(j) will contain the Global Id in the monilithic map corresponding to the local id j
	std::vector<std::vector<int> > myGlobalElBlock(M_nBlocks);
	std::vector<int> myGlobalElMono(numMyElMono);

	// Iterator for the vector myGlobalElMono
	std::vector<int>::iterator myGlobalElMonoIt;
	myGlobalElMonoIt = myGlobalElMono.begin();

	//I do a loop on the blocks
	for(std::vector<std::vector<int> >::iterator myGlobalElBlockIt = myGlobalElBlock.begin();
			myGlobalElBlockIt != myGlobalElBlock.end();
			++myGlobalElBlockIt)
	{
		// Get the block number from the iterator.
		UInt iblock = (UInt) (myGlobalElBlockIt - myGlobalElBlock.begin());
		// I will copy the Global Elements of the map-iblock into the vector *myGlobalElBlockIt
		myGlobalElBlockIt->resize(numMyElBlock[iblock]);
		M_blockShiftedMap[iblock]->MyGlobalElements( &(*myGlobalElBlockIt)[0] );

		// Now I shift the Global Elements
		for(std::vector<int>::iterator it=myGlobalElBlockIt->begin(); it!=myGlobalElBlockIt->end(); ++it)
			*it += shift[iblock];

		// the vector myGlobalElMonoIt is the concatenation of the *myGlobalElBlockIt vectors
		myGlobalElMonoIt = std::copy(myGlobalElBlockIt->begin(), myGlobalElBlockIt->end(), myGlobalElMonoIt);

		// Now I say that blockMaps[iblock] is the map with the Shifted Global Elements
		M_blockShiftedMap[iblock].reset(new map_Type(numGlobalElBlock[iblock], myGlobalElBlockIt->size(),
				&(*myGlobalElBlockIt)[0], shift[iblock], comm));

	}

	// I create the monolithic map
	M_monolithicMap.reset(new map_Type(numGlobalElMono, myGlobalElMono.size(), &myGlobalElMono[0], 0, comm));

	// I create the importer from the block maps to the monolithic map and viceversa
	for(UInt iblock=0; iblock<M_nBlocks; ++iblock)
	{
		M_block2mono[iblock].reset(new import_Type(*M_monolithicMap, *M_blockShiftedMap[iblock]));
		M_mono2block[iblock].reset(new import_Type(*M_blockShiftedMap[iblock], *M_monolithicMap));
	}
}


Epetra_Map * blockMap2Map(const Epetra_BlockMap * blockMap)
{
	Epetra_Map * map = new Epetra_Map
			        (blockMap->NumGlobalElements(),
					 blockMap->NumMyElements(),
					 blockMap->MyGlobalElements(),
					 blockMap->IndexBase(),
					 blockMap->Comm()
					 );
	return map;
}

} /*end name space*/
