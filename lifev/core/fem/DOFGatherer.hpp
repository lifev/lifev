//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010, 2011, 2012, 2013 EPFL, Politecnico di Milano, Emory University

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
    @brief Class that produces a list of dof GID from a list of element LIDs

    @date 2013-12-08
    @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef DOFGATHERER_HPP_
#define DOFGATHERER_HPP_

#include <algorithm>
#include <iostream>
#include <vector>
#include <set>

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/mesh/GraphUtil.hpp>

namespace LifeV {

using namespace GraphUtil;

//! Class that produces a list of dof GID from a list of element LIDs
/*!
    @author Radu Popescu <radu.popescu@epfl.ch>

    Objects of this class are constructed using an table of element LIDs
    (as an idTablePtr_Type object) and a finite element space (standard FESpace).
    The element LIDs are stored in a vector or vectors (see idTablePtr_Type in
    lifev/core/mesh/GraphUtil.hpp), where each vector contains the IDs of
    the elements in a second stage partition (for ShyLU_MT).

    The DOFGatherer performs two operations, relevant to the second stage
    partitioning method used by multi-threaded solvers such as ShyLU_MT:

    1. It gathers for each set of element IDs, from the FESpace, all the dof
       GIDs that are associated with the elements.
    2. It will classify the dofs into unique sets, one for each element ID set,
       and a shared set, representing the dofs on the interface between the
       second stage parts. These shared dofs are removed from the unique sets
       and are place into a separate set, in the same table.

    The classified dof sets can be accessed through the dofGIDTable method.
 */
template <typename MeshType>
class DOFGatherer {
public:
	//! Public typedefs
	typedef FESpace<MeshType, MapEpetra>       feSpace_Type;
	typedef boost::shared_ptr<feSpace_Type>    feSpacePtr_Type;

	//! Constructors and destructor
	//@{
	//! Constructor
	/*!
	 * \param elementIds - table of element LIDs from the second stage mesh partition
	 * \param feSpace - finite element space
	 */
	DOFGatherer(const idTablePtr_Type& elementIds,
				const feSpacePtr_Type& feSpace);

	//! Destructor
	virtual ~DOFGatherer();
	//@}

	//! Get methods
	//@{
	//! Return a pointer to the classified dof table
	const idSetGroupPtr_Type& dofGIDTable() const
	{
		return M_dofGIDTable;
	}
	//@}

	//! Public methods
	//@{
	//! Display the contents of the dof table
	void showMe(std::ostream& s = std::cout);
	//@}
private:
	//! Private methods
	//@{
	//! This method performs all the computations
	void run();

	//! Helper method to obtain dof GIDs associated with element LIDs
	/*!
	 * \param elementLIDs - input list of element LIDs for which the dof lookup is
	 * 				 desired
	 * \param dofGIDs - output set of dof GIDs
	 */
	void getDofGIDs(const idList_Type& elementLIDs,
					idSet_Type& dofGIDs);

	//! Private methods that computes unique and shared dofs in M_dofGIDTable
	void dofClassification();

	//! Helper method that removes a subset of IDs from an idSet_Type
	void removeValues(const idSet_Type& valuesToRemove,
					  idSet_Type& targetSet);

	//! Method that computes whether two (ordered) sets of IDs intersect
	const bool setsCouldIntersect(const idSet_Type& first,
								  const idSet_Type& second) const
	{
		Int x0 = *(first.begin());
		Int x1 = *(first.rbegin());
		Int y0 = *(second.begin());
		Int y1 = *(second.rbegin());
		return (x0 < y1) && (y0 < x1);
	}

	//! Method that checks if two ordered sets of IDs have a single shared value
	const bool singleIntersect(const idSet_Type& first,
							   const idSet_Type& second,
							   Int& returnValue) const
	{
		Int x0 = *(first.begin());
		Int x1 = *(first.rbegin());
		Int y0 = *(second.begin());
		Int y1 = *(second.rbegin());

		if (y0 == x1) {
			returnValue = y0;
			return true;
		} else if (x0 == y1) {
			returnValue = x0;
			return true;
		} else {
			returnValue = -1;
			return false;
		}
	}
	//@}

	// Private data
	const idTablePtr_Type M_elementIds;
	const boost::shared_ptr<FESpace<MeshType, MapEpetra> > M_feSpace;
	idSetGroupPtr_Type M_dofGIDTable;
};

// Method implementation

template <typename MeshType>
DOFGatherer<MeshType>::DOFGatherer(const idTablePtr_Type& elementIds,
								   const feSpacePtr_Type& feSpace)
	: M_elementIds(elementIds),
	  M_feSpace(feSpace),
	  M_dofGIDTable()
{
	run();
}

template <typename MeshType>
DOFGatherer<MeshType>::~DOFGatherer()
{
}

template <typename MeshType>
void DOFGatherer<MeshType>::run()
{
	// We need to cycle on each element LID list and obtain the dof GIDs
	// that are associated with the given element sets
	const Int numParts = M_elementIds->size();
	M_dofGIDTable.reset(new idSetGroup_Type(numParts));

	for (Int i = 0; i < numParts; ++i) {
		M_dofGIDTable->at(i).reset(new idSet_Type);
		getDofGIDs(*(M_elementIds->at(i)), *(M_dofGIDTable->at(i)));
	}

	// From the newly produces dof GID lists, we need to remove shared dof GIDs,
	// which represent the interface between the second stage parts and place
	// them into a separate GID list

	// Add a new dof set to the output table, representing the id of dofs shared
	// among the other sets
	M_dofGIDTable->push_back(idSetPtr_Type());
	M_dofGIDTable->at(numParts).reset(new idSet_Type);

	// Call the classification method
	dofClassification();
}

template <typename MeshType>
void DOFGatherer<MeshType>::getDofGIDs(const idList_Type& elementLIDs,
									   idSet_Type& dofGIDs)
{
	const UInt nbDof = M_feSpace->refFE().nbDof();
	const UInt fieldDim = M_feSpace->fieldDim();

	for (Int iElement = 0; iElement < elementLIDs.size(); ++iElement) {
		Int elemLID = elementLIDs[iElement];
		for (UInt iBlock = 0; iBlock < fieldDim; ++iBlock)
		{
			// Set the row global indices in the local matrix
			for (UInt iDof = 0; iDof < nbDof; ++iDof)
			{
				dofGIDs.insert(M_feSpace->dof().localToGlobalMap(elemLID, iDof)
					+ iBlock * M_feSpace->dof().numTotalDof());
			}
		}
	}
}

template <typename MeshType>
void DOFGatherer<MeshType>::showMe(std::ostream& s)
{
	for (Int i = 0; i < M_dofGIDTable->size(); ++i) {
		s << "DOF set " << i << " (size " << M_dofGIDTable->at(i)->size() << "): ";
		idSet_Type& currentDofs = *(M_dofGIDTable->at(i));
		for (idSet_Type::const_iterator it = currentDofs.begin();
			 it != currentDofs.end(); ++it) {
			s << *it << " ";
		}
		s << std::endl;
	}
}

template<typename MeshType>
void DOFGatherer<MeshType>::dofClassification()
{
	// TODO: try to do this with one pass over GIDs, to see if it's faster
	Int numSets = M_dofGIDTable->size() - 1;
	idSet_Type& commonSet = *(M_dofGIDTable->at(numSets));

	for (Int i = 0; i < numSets - 1; ++i) {
		idSet_Type& currentSet = *(M_dofGIDTable->at(i));
		for (Int j = i + 1; j < numSets; ++j) {
			idSet_Type& comparisonSet = *(M_dofGIDTable->at(j));

			if (setsCouldIntersect(currentSet, comparisonSet)) {
				Int temp = 0;
				if (singleIntersect(currentSet, comparisonSet, temp)) {
					currentSet.erase(temp);
					comparisonSet.erase(temp);
					commonSet.insert(temp);
				} else {
					idSet_Type intersectionSet;
					std::set_intersection(currentSet.begin(), currentSet.end(),
										  comparisonSet.begin(), comparisonSet.end(),
										  std::inserter(intersectionSet, intersectionSet.begin()));
					removeValues(intersectionSet, currentSet);
					removeValues(intersectionSet, comparisonSet);

					commonSet.insert(intersectionSet.begin(), intersectionSet.end());
				}
			}
		}
	}
}

template<typename MeshType>
void DOFGatherer<MeshType>::removeValues(const idSet_Type& valuesToRemove,
										 idSet_Type& targetSet)
{
	for (idSet_Type::const_iterator it = valuesToRemove.begin();
		 it != valuesToRemove.end(); ++it) {
		targetSet.erase(*it);
	}
}

} /* namespace LifeV */

#endif /* DOFGATHERER_HPP_ */
