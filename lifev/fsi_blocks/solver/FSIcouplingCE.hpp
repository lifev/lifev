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
@brief Settings - File which handles the coupling blocks of fluid and solid variables 
 at the fluid-structure interface when conforming discretizations are used.

@author Davide Forti <davide.forti@epfl.ch>
@date 28-10-2014

@maintainer Simone Deparis <simone.deparis@epfl.ch>
*/


#ifndef FSICOUPLINGCE_H
#define FSICOUPLINGCE_H

#include <lifev/core/LifeV.hpp>

// datafile
#include <lifev/core/filter/GetPot.hpp>

// includes for matrices and vector
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>

namespace LifeV
{

//! FSIcouplingCE - File handling the coupling blocks when conforming discretizations are used.
/*!
*  @author Davide Forti
*/

class FSIcouplingCE
{
public:

	// Public typedefs

    typedef Epetra_Comm comm_Type;
	typedef boost::shared_ptr< comm_Type > commPtr_Type;

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    typedef MapEpetra map_Type;
	typedef boost::shared_ptr<map_Type> mapPtr_Type;

	typedef MatrixEpetra<Real> matrix_Type;
	typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

	typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    typedef FESpace< mesh_Type, map_Type > FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type> FESpacePtr_Type;

    //! Constructor
	FSIcouplingCE(const commPtr_Type& communicator);

    //! Destructor
    ~FSIcouplingCE();

    //! @name Methods
    //@{
    
    //! Builds the coupling blocks
    /*!
     \param locDofMap map with dofs at the fluid-structure interface
     \param lambda_num_structure if true the DOFs at the interface are numbered wrt the solid mesh
     \param useBDF if true supposes that for the structure a BDF scheme is used
     */
    void buildBlocks ( std::map<ID, ID> const& locDofMap, const bool& lambda_num_structure, bool useBDF = false );

    //! Set parameters. To be used when Newmark is used on the structure.
    /*!
     \param timeStep value of the timestep used
     \param interfaceDofs number of interface dofs
     \param beta beta coefficient Newmark scheme
     \param gamma gamma coefficient Newmark scheme
     \param interfaceMap map interface dofs
     \param fluidVelocityFESpace FE space fluid velocity
     \param structureDisplacementFESpace FE space solid displacement
     \param numerationInterface vector with global numeration of dofs at the interface
     */
    void setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& beta, const Real& gamma,
    			 const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
    			 const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface );

    //! Set parameters. To be used when Newmark is used on the structure.
    /*!
     \param timeStep value of the timestep used
     \param interfaceDofs number of interface dofs
     \param coefficientBDF coefficient BDF scheme for first derivative structure
     \param interfaceMap map interface dofs
     \param fluidVelocityFESpace FE space fluid velocity
     \param structureDisplacementFESpace FE space solid displacement
     \param numerationInterface vector with global numeration of dofs at the interface
     */
    void setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& coefficientBDF,
    		     const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
    		     const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface );

    //@}
    
    //! @name Get Methods
    //@{
    
    matrixPtr_Type lambdaToFluidMomentum() const
    {
    	return M_lambdaToFluidMomentum;
    }

    matrixPtr_Type lambdaToStructureMomentum() const
    {
    return M_lambdaToStructureMomentum;
    }

    matrixPtr_Type fluidVelocityToLambda() const
    {
    	return M_fluidVelocityToLambda;
    }

    matrixPtr_Type structureDisplacementToLambda() const
    {
    	return M_structureDisplacementToLambda;
    }

    matrixPtr_Type structureDisplacementToFluidDisplacement() const
    {
    	return M_structureDisplacementToFluidDisplacement;
    }

    //@}

private:

    //! communicator
    commPtr_Type M_comm;

    matrixPtr_Type M_lambdaToFluidMomentum;
    matrixPtr_Type M_lambdaToStructureMomentum;
    matrixPtr_Type M_fluidVelocityToLambda;
    matrixPtr_Type M_structureDisplacementToLambda;
    matrixPtr_Type M_structureDisplacementToFluidDisplacement;

    FESpacePtr_Type M_fluidVelocityFESpace;
    FESpacePtr_Type M_structureDisplacementFESpace;

    mapPtr_Type M_interfaceMap;

    vectorPtr_Type M_numerationInterface;

    Real M_interface;
    Real M_timeStep;
    Real M_beta;
    Real M_gamma;
    Real M_coefficientBDF;
};

} // end namespace LifeV

#endif // end FSICOUPLINGCE_H
