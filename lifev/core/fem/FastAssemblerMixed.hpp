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
     @brief This file contains the definition of the FastAssemblerMixed class.

     This function is used to perform an efficient assembly of FE matrices
     where the test and trial functions are not the same.

     @date 06/2016
     @author Davide Forti <davide.forti@epfl.ch>
 */

#ifndef FASTASSEMBLERMIXED_HPP
#define FASTASSEMBLERMIXED_HPP

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/fem/FESpace.hpp>

namespace LifeV
{

class FastAssemblerMixed
{
public:

	typedef RegionMesh< LinearTetra > mesh_Type;
    typedef boost::shared_ptr<mesh_Type>  meshPtr_Type;

    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr< comm_Type > commPtr_Type;

    typedef QuadratureRule qr_Type;
    typedef boost::shared_ptr< qr_Type > qrPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra> fespace_Type;
    typedef boost::shared_ptr<fespace_Type> fespacePtr_Type;


    //! Constructor
    /*!
     * @param mesh - input mesh
     * @param comm - communicator
     * @param refFE_test - reference FE space of test functions
     * @param refFE_trial - reference FE space of trial functions
     * @param qr_integration - quadrature rule to be used for the integration
     */
    FastAssemblerMixed( const meshPtr_Type& mesh, const commPtr_Type& comm,
            			const ReferenceFE* refFE_test, const ReferenceFE* refFE_trial,
            			const qr_Type* qr_integration );

    //! Destructor
	~FastAssemblerMixed();

	//! @name Methods
	//@{

	//! Allocate space for members before the assembly
	/*!
	 * @param numElements - data file
	 * @param fe_test - current FE test functions
	 * @param fespace_test - FE space test functions
	 * @param fe_trial - current FE trial functions
	 * @param fespace_trial - FE space trial functions
	 */
	void allocateSpace( const int& numElements,
            			CurrentFE* fe_test, const fespacePtr_Type& fespace_test,
            			CurrentFE* fe_trial, const fespacePtr_Type& fespace_trial );

	//! FE Assembly block (0,1) of Navier-Stokes
	/*!
	 * @param matrix - global matrix, in this case the block (0,1) of Navier-Stokes
	 */
	void assemble_NS_block01 ( matrixPtr_Type& matrix );

	//! FE Assembly block (1,0) of Navier-Stokes
	/*!
	 * @param matrix - global matrix, in this case the block (1,0) of Navier-Stokes
	 */
	void assemble_NS_block10 ( matrixPtr_Type& matrix );

	//@}

private:

	meshPtr_Type M_mesh;
	commPtr_Type M_comm;

	int M_numElements;
	int M_numScalarDofs_test;
	int M_numScalarDofs_trial;

	double * M_detJacobian;
	double *** M_invJacobian;

	double ** M_phi_test;
	double *** M_dphi_test;
	double ** M_phi_trial;
	double *** M_dphi_trial;
    double ** M_elements_test;
    double ** M_elements_trial;

	const qr_Type* M_qr_integration;
	const ReferenceFE* M_referenceFE_test;
	const ReferenceFE* M_referenceFE_trial;

	double**** M_vals;
	int** M_rows;
	int** M_cols;

};

} // Namespace LifeV

#endif // FASTASSEMBLERMIXED_HPP
