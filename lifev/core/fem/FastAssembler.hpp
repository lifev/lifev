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
     @brief This file contains the definition of the FastAssembler class.

     This function is used to perform an efficient assembly of FE matrices
     where the test and trial functions are the same.

     @date 06/2016
     @author Davide Forti <davide.forti@epfl.ch>
 */

#ifndef FASTASSEMBLER_HPP
#define FASTASSEMBLER_HPP

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

class FastAssembler
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
     * @param refFE - reference FE space of test functions
     * @param qr - quadrature rule to be used for the integration
     */
    FastAssembler( const meshPtr_Type& mesh, const commPtr_Type& comm, const ReferenceFE* refFE, const qr_Type* qr );

    //! Destructor
	~FastAssembler();

	//! @name Methods
	//@{

	//! Allocate space for members before the assembly
	/*!
	 * @param numElements - data file
	 * @param fe - current FE
	 * @param fespace - FE space
	 */
	void allocateSpace( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace );

	//! Allocate space for members before the assembly
	/*!
	 * @param numElements - data file
	 * @param fe - current FE
	 * @param fespace - FE space
	 * @param meshSub_elements - list of indices if one wants to allocate space only for a portion of the elements of the mesh
	 */
	void allocateSpace( const int& numElements, CurrentFE* fe, const fespacePtr_Type& fespace, const UInt* meshSub_elements );

	//! Allocate space for supg before the assembly
	void allocateSpace_SUPG( );

	//! FE Assembly of scalar grad-grad
	/*!
	 * @param matrix - global matrix
	 */
	void assembleGradGrad_scalar( matrixPtr_Type& matrix );

	//! FE Assembly of vectorial grad-grad
	/*!
	 * @param matrix - global matrix
	 */
	void assembleGradGrad_vectorial( matrixPtr_Type& matrix );

	//! FE Assembly of vectorial mass matrix
	/*!
	 * @param matrix - global matrix
	 */
	void assembleMass_vectorial( matrixPtr_Type& matrix );

	//! FE Assembly of scalar mass matrix
	/*!
	 * @param matrix - global matrix
	 */
	void assembleMass_scalar( matrixPtr_Type& matrix );

	//! FE Assembly of NS constant terms (no scaling by coefficients like density or bdf)
	/*!
	 * @param matrix - global matrix
	 */
	void NS_constant_terms_00( matrixPtr_Type& matrix ); // Navier-Stokes constant terms belonging to block (0,0): mass + stiffness

	//! FE Assembly of NS constant terms (no scaling by coefficients like viscosity)
	/*!
	 * @param matrix - global matrix
	 * @param matrix - velocity vector
	 */
	void assembleConvective( matrix_Type& matrix, const vector_Type& u_h );

	//! FE Assembly of NS constant terms (no scaling by coefficients like viscosity)
	/*!
	 * @param matrix - global matrix
	 * @param matrix - velocity vector
	 */
	void assembleConvective( matrixPtr_Type& matrix, const vector_Type& u_h );

	//! FE Assembly of SUPG terms - work in progress
	/*!
	 * @param matrix - global matrix
	 * @param matrix - vector extrapolapolated velocity
	 */
	void assemble_SUPG_block00( matrixPtr_Type& matrix, const vector_Type& u_h );

	//@}

private:

	meshPtr_Type M_mesh;
	commPtr_Type M_comm;

	int M_numElements;
	int M_numScalarDofs;
	int M_numElementsMerked;

	double * M_detJacobian;
	double *** M_invJacobian;

	double ** M_phi;
	double *** M_dphi;
    double ** M_elements;

	const qr_Type* M_qr;
	const ReferenceFE* M_referenceFE;

	double*** M_vals;
	double***** M_vals_supg;
	int** M_rows;
	int** M_cols;

    bool M_useSUPG;

};

} // Namespace LifeV

#endif // FASTASSEMBLER_HPP
