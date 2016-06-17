//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the MatrixGraph class.

     This function is used to perform an efficient assembly of FE matrices
     graph.

     @date 06/2016
     @author Davide Forti <davide.forti@epfl.ch>
 */

#ifndef MATRIXGRAPH_HPP
#define MATRIXGRAPH_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <Epetra_FECrsGraph.h>

namespace LifeV
{

class MatrixGraph
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
    typedef boost::shared_ptr<mesh_Type>  meshPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr< comm_Type > commPtr_Type;

    typedef Epetra_FECrsGraph graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra> fespace_Type;
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;

    //! Constructor
    /*!
     * @param mesh - input mesh
     * @param comm - communicator
     */
    MatrixGraph( const meshPtr_Type& mesh, const commPtr_Type& comm, const ReferenceFE* refFE );

    //! Destructor
	~MatrixGraph();

	//! @name Methods
	//@{

	//! Build the matrix graph
	/*!
	 * @param numElements - data file
	 * @param fe - current FE
	 * @param fespace - FE space
	 */
	void buildGraph( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace, graphPtr_Type& matrix_graph );

	//@}

private:

	meshPtr_Type M_mesh;
	commPtr_Type M_comm;

	int M_numElements;
	int M_numScalarDofs;

    double ** M_elements;
    const ReferenceFE* M_referenceFE;

	int** M_rows;
	int** M_cols;

    // methods

};

} // Namespace LifeV

#endif // MATRIXGRAPH_HPP
