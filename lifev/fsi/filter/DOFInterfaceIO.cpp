//@HEADER
/*
*******************************************************************************

Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
Copyright (C) 2010, 2011, 2012 EPFL, Politecnico di Milano, Emory University

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
  @brief Class that handles I/O of fluid-solid DOF interfaces in FSI

  @date 2013-04-30
  @author Radu Popescu <radu.popescu@epfl.ch>
  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#include <sstream>
#include <map>
#include <algorithm>

#include <lifev/fsi/filter/DOFInterfaceIO.hpp>

#ifdef LIFEV_HAS_HDF5
#ifdef HAVE_MPI

namespace LifeV
{

DOFInterfaceIO::DOFInterfaceIO (const std::string& fileName,
                                const commPtr_Type& comm) :
    M_comm (comm),
    M_fileName (fileName)
{
    M_myRank = M_comm->MyPID();
}

void DOFInterfaceIO::write (const interfaceVectorPtr_Type& interfaces)
{
    M_HDF5IO.openFile (M_fileName, M_comm, false);

	const UInt numParts = interfaces->size();
	std::vector<UInt> sizes(numParts, 0);

	UInt maxPoints = 0;
	for (UInt i = 0; i < numParts; ++i) {
		interface_Type& curIf = *(interfaces->at(i));
		UInt curSize = curIf.localDofMap().size();
		sizes[i] = curSize;
		if (maxPoints < curSize) {
			maxPoints = curSize;
		}
	}

	hsize_t sizeSpaceDims[2] = {numParts, 1};
	hsize_t sizeCount[2] = {1, 1};

	hsize_t currentSpaceDims[2];
	hsize_t currentCount[2];
	currentSpaceDims[0] = numParts;
	currentSpaceDims[1] = maxPoints;
	currentCount[0] = 1;
	currentCount[1] = currentSpaceDims[1];

	// Create a mini table storing the size of each interface map
	M_HDF5IO.createTable ("Sizes", H5T_STD_U32BE, sizeSpaceDims);

	M_HDF5IO.createTable ("Keys", H5T_STD_U32BE, currentSpaceDims);
	M_HDF5IO.createTable ("Values", H5T_STD_U32BE, currentSpaceDims);

	M_keyBuffer.resize (currentCount[0] * currentCount[1], 0);
	M_valueBuffer.resize (currentCount[0] * currentCount[1], 0);

	for (UInt i = 0; i < numParts; ++i) {
		const std::map<UInt, UInt>& curMap =
				interfaces->at(i)->localDofMap();
		UInt j = 0;
		for (std::map<UInt, UInt>::const_iterator it = curMap.begin();
			 it != curMap.end(); ++it) {
			M_keyBuffer[j] = it->first;
			M_valueBuffer[j] = it->second;
			++j;
		}

		hsize_t sizeOffset[2] = {i, 0};
		M_HDF5IO.write ("Sizes", H5T_NATIVE_UINT, sizeCount,
						sizeOffset, &sizes[i]);

		hsize_t currentOffset[2] = {i* currentCount[0], 0};
		M_HDF5IO.write ("Keys", H5T_NATIVE_UINT, currentCount,
						currentOffset, &M_keyBuffer[0]);
		M_HDF5IO.write ("Values", H5T_NATIVE_UINT, currentCount,
						currentOffset, &M_valueBuffer[0]);
	}

	M_HDF5IO.closeTable("Sizes");
    M_HDF5IO.closeTable("Keys");
    M_HDF5IO.closeTable("Values");

    M_HDF5IO.closeFile();

    M_keyBuffer.resize(0);
    M_valueBuffer.resize(0);
}

void DOFInterfaceIO::read (dofMapPtr_Type& interface)
{
	interface.reset(new dofMap_Type());

    M_HDF5IO.openFile (M_fileName, M_comm, true);

    hsize_t sizeSpaceDims[2];

	M_HDF5IO.openTable ("Sizes", sizeSpaceDims);

	UInt mySize;
    hsize_t sizeCount[2] = {1, 1};
    hsize_t sizeOffset[2] = {M_myRank, 0};
	M_HDF5IO.read ("Sizes", H5T_NATIVE_UINT, sizeCount, sizeOffset,
				   &mySize);

	M_HDF5IO.closeTable("Sizes");

	hsize_t currentSpaceDims[2];
	hsize_t currentCount[2];
	hsize_t currentOffset[2];

	M_HDF5IO.openTable ("Keys", currentSpaceDims);
	M_HDF5IO.openTable ("Values", currentSpaceDims);

	currentCount[0] = 1;
	currentCount[1] = currentSpaceDims[1];
	currentOffset[0] = currentCount[0] * M_myRank;
	currentOffset[1] = 0;

    M_keyBuffer.resize(currentCount[0] * currentCount[1], 0);
    M_valueBuffer.resize(currentCount[0] * currentCount[1], 0);

	M_HDF5IO.read ("Keys", H5T_NATIVE_UINT, currentCount, currentOffset,
				   &M_keyBuffer[0]);
	M_HDF5IO.read ("Values", H5T_NATIVE_UINT, currentCount, currentOffset,
				   &M_valueBuffer[0]);

	M_HDF5IO.closeTable ("Keys");
	M_HDF5IO.closeTable ("Values");

    M_HDF5IO.closeFile();

    for (UInt i = 0; i < mySize; ++i)
    {
        interface->insert (std::make_pair (M_keyBuffer[i], M_valueBuffer[i]) );
    }

    M_keyBuffer.resize(0);
    M_valueBuffer.resize(0);
}

} /* namespace LifeV */

#endif /* HAVE_MPI */
#endif /* LIFEV_HAS_HDF5 */

