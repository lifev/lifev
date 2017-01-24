#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H 1

//@HEADER
/*
 *******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

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
   @brief Implementation of Neohookean structural model
   @author Davide Forti <davide.forti@epfl.ch>
   @maintainer Antonello Gerbi <antonello.gerbi@epfl.ch>
   @date 25-05-2016

   This file contains an ETA implementation of a neohookean structural model. This class
   temporary as it will disappear when the new structure model will be released.

*/

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/fsi_blocks/solver/AssemblyElementalStructure.hpp>

#include <lifev/fsi_blocks/solver/ExpressionDefinitions.hpp>

namespace LifeV
{

class NeoHookean
{

	typedef Epetra_Comm comm_Type;

	typedef boost::shared_ptr< comm_Type > commPtr_Type;

    typedef VectorEpetra vector_Type;

    typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

    typedef MatrixEpetra<Real> matrix_Type;

    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef boost::shared_ptr<BCHandler> bcPtr_Type;

    typedef RegionMesh<LinearTetra> mesh_Type;

    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    typedef MapEpetra map_Type;

    typedef ETFESpace< mesh_Type, map_Type, 3, 3 > ETFESpace_displacement;

    typedef ExpressionDefinitions::deformationGradient_Type tensorF_Type;

    typedef ExpressionDefinitions::determinantTensorF_Type determinantF_Type;

    typedef ExpressionDefinitions::rightCauchyGreenTensor_Type tensorC_Type;

    typedef ExpressionDefinitions::minusTransposedTensor_Type minusT_Type;

    typedef ExpressionDefinitions::traceTensor_Type traceTensor_Type;

public:

    // Constructor
    /*!
     *  @param communicator Epetra communicator
     */
    NeoHookean( const commPtr_Type& communicator );

    // empty destructor
    ~NeoHookean();

    //! Set the physical parameters of the structure
    /*!
     *  @param density Density of the structure
     *  @param young Young modulus of the structure
     *  @param poisson Poisson ratio of the structure
     */
    void setCoefficients ( const Real density, const Real young, const Real poisson);

    //! Setup: build the FE space and the ETA fespace
    /*!
     *  @param mesh Computational mesh
     *  @param dOrder Degree of the finite element used
     */
    void setup( const meshPtr_Type& mesh, const std::string dOrder );

    //! Evaluate the residual of the problem
    /*!
     *  @param solution Solution vector at the previous Newton iterate
     *  @param coefficient Coefficient which goes in front of the mass matrix depending on the time discretization scheme used
     *  @param csi Vector containing the solution at the previous time steps which depends on the time discretization scheme used
     *  @param residual Residual vector to be assembled
     */
    void evaluate_residual( const vectorPtr_Type& solution, const Real& coefficient, const vectorPtr_Type& csi, vectorPtr_Type& residual );

    //! Updates the Jacobian matrix
    /*!
     *  @param solution Solution vector at the previous Newton iterate
     *  @param coefficient Coefficient which goes in front of the mass matrix depending on the time discretization scheme used
     *  @param jacobian Jacobian matrix to be updated
     */
    void update_jacobian(const vectorPtr_Type& solution, const Real& coefficient, matrixPtr_Type& jacobian );

    //! Getter of the standard FE space
    /*!
     *  @return M_displacementFESpace FE space used for the displacement
     *
     */
    const boost::shared_ptr<FESpace<mesh_Type, map_Type> >& fespace ( ) const { return M_displacementFESpace; };

    //! Getter of the ETA FE space
    /*!
     *  @return M_displacementFESpace_ETA ETA FE space used for the displacement
     *
     */
    const boost::shared_ptr<ETFESpace_displacement >& et_fespace ( ) const { return M_displacementFESpace_ETA; };

private:

    // communicator
    commPtr_Type M_comm;

    // Coefficients for the structure
    Real M_density;

    // Lame coefficients computed on M_young and M_poisson
    Real M_bulk;
    Real M_mu;
    Real M_offset;

    matrixSmall_Type M_identity;

    // order FE space displacement
    std::string M_dOrder;

    // FE space
    boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_displacementFESpace;

    // ET FE Space
    boost::shared_ptr<ETFESpace_displacement > M_displacementFESpace_ETA;

    // jacobian matrix
    matrixPtr_Type M_jacobian;

};

} // end namespace LifeV

#endif
