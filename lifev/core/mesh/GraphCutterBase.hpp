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
 @brief Graph cutter base class (abstract)

 @date 06-02-2013
 @author Radu Popescu <radu.popescu@epfl.ch>
 */

#ifndef GRAPH_CUTTER_BASE_H
#define GRAPH_CUTTER_BASE_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_ptr.hpp>
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

namespace LifeV {

//! Graph cutter base class (abstract)
/*!
 @author Radu Popescu <radu.popescu@epfl.ch>

 TODO: detailed class decription
 */
template<typename MeshType>
class GraphCutterBase {
public:
	//! @name Public Types
	//@{
	typedef Teuchos::ParameterList pList_Type;
	typedef boost::shared_ptr<std::vector<std::vector<Int> > > graph_Type;
	//@}

	//! @name Constructor & Destructor
	//@{
	//! Default constructor
	GraphCutterBase() {
	}

	//! Destructor
	virtual ~GraphCutterBase() {
	}
	//@}

	//! @name Public methods
	//@{
	//! Performs the graph partitioning
	virtual Int run() = 0;
	//@}

	//! @name Get Methods
	//@{
	//! Get a pointer to one of the partitions
	virtual const std::vector<Int>& getPart(const UInt i) const = 0;
	virtual std::vector<Int>& getPart(const UInt i) = 0;

	//! Get the entire partitioned graph, wrapped in a smart pointer
	virtual graph_Type getGraph() = 0;

	//! Return the number of parts
	virtual const UInt numParts() const = 0;
	//@}

private:
	//! @name Private methods
	//@{
	//! Set values for all the parameters, with default values where needed
	virtual void setParameters(pList_Type& parameters) = 0;

	//@}
};

} // Namespace LifeV

#endif // GRAPH_CUTTER_BASE_H
