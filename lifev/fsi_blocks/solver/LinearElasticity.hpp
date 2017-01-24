#ifndef LINEARELASTICITY_H
#define LINEARELASTICITY_H 1

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
   @brief Implementation of Linear elastic structural model
   @author Davide Forti <davide.forti@epfl.ch>
   @maintainer Antonello Gerbi <antonello.gerbi@epfl.ch>
   @date 25-05-2016

   This file contains an ETA implementation of a linear elastic structural model. This class
   temporary as it will disappear when the new structure model will be released.

   This class contains also an implementation of a thin layer structural model which can be
   used together with the 3D structure to form a multilayered structure. As a reference for applications
   of such multilayered computational model in the context of fluid-structure interaction, please refer to:

   * D. Forti, M. Bukac,  A. Quaini, S. Čanić,  and S. Deparis, "A Monolithic Approach to Fluid-Composite Structure
   Interaction" available as NA & SC Preprint Series No. 47 - February, 2016.
   (http://www.uh.edu/nsm/_docs/math/NASC-preprint-series/2015_2016/Preprint_No16-47.pdf)

*/

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>

#include <lifev/core/fem/BCManage.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/eta/expression/Integrate.hpp>

namespace LifeV
{

class LinearElasticity
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

public:

    // Constructor
    /*!
     *  @param communicator Epetra communicator
     */
    LinearElasticity( const commPtr_Type& communicator );

    // empty destructor
    ~LinearElasticity();

    //! Set the physical parameters of the structure
    /*!
     *  @param density Density of the structure
     *  @param young Young modulus of the structure
     *  @param poisson Poisson ratio of the structure
     */
    void setCoefficients ( const Real density, const Real young, const Real poisson);

    //! Set the physical parameters of the thin layered structure
    /*!
     *  @param density Density of the thin layered structure
     *  @param young Young modulus of the thin layered structure
     *  @param poisson Poisson ratio of the thin layered structure
     *  @param thickness Thickness of the thin layered structure
     *  @param interface Flag of the thin layered structure
     */
    void setCoefficientsThinLayer ( const Real density, const Real young, const Real poisson, const Real thickness, const UInt interface );

    //! Setup: build the FE space and the ETA fespace
    /*!
     *  @param mesh Computational mesh
     *  @param dOrder Degree of the finite element used
     */
    void setup( const meshPtr_Type& mesh, const std::string dOrder );

    //! Assembling matrices (mass and stiffness) which for this model are constant in a time dependent simulation
    /*!
     *  @param timestep Value of the time step used
     *  @param bc Boundary conditions
     *  @param useBDF Boolean variable which specifies if BDF is used for the structure. By default it is false (therefore, Newmark is used).
     */
    void assemble_matrices ( const Real timestep, const Real beta, bcPtr_Type & bc, bool useBDF = false );

    //! Getter of the mass matrix whithout boundary conditions applied
    /*!
     *  @return M_mass_no_bc mass matrix whithout boundary conditions applied
     */
    matrixPtr_Type const& mass_matrix_no_bc ( ) const { return M_mass_no_bc; };

    //! Getter of the stiffness matrix whithout boundary conditions applied
    /*!
     *  @return M_stiffness_no_bc stiffness matrix whithout boundary conditions applied
     */
    matrixPtr_Type const& stiffness_matrix_no_bc ( ) const { return M_stiffness_no_bc; };

    //! Getter of the Jacobian matrix
    /*!
     *  @return M_jacobian stiffness matrix
     */
    matrixPtr_Type const& jacobian ( ) const { return M_jacobian; };

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
    Real M_young;
    Real M_poisson;

    // Lame coefficients computed on M_young and M_poisson
    Real M_lambda;
    Real M_mu;

    // order FE space displacement
    std::string M_dOrder;

    // FE space
    boost::shared_ptr<FESpace<mesh_Type, map_Type> > M_displacementFESpace;

    // ET FE Space
    boost::shared_ptr<ETFESpace_displacement > M_displacementFESpace_ETA;

    // Matrices without bc applied
    matrixPtr_Type M_mass_no_bc;
    matrixPtr_Type M_stiffness_no_bc;

    // Matrices without bc applied thin layer
    matrixPtr_Type M_mass_no_bc_thin;
    matrixPtr_Type M_stiffness_no_bc_thin;

    // jacobian matrix
    matrixPtr_Type M_jacobian;

    // Parameters used to deal with the thin layer
    bool				   M_thinLayer;
    Real 				   M_thinLayerThickness;
    Real 				   M_thinLayerDensity;
    Real 				   M_thinLayerLameI;
    Real 				   M_thinLayerLameII;
    UInt 				   M_interfaceFlag;

};

} // end namespace LifeV

#endif
