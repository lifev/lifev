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
@brief Settings - File handling the solution of the FSI problem

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

//! FSIcouplingCE - File handling the coupling when using the monolithic CE time-discretization
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

    void buildBlocks ( std::map<ID, ID> const& locDofMap, const bool& lambda_num_structure );

    void setUp ( const Real& timeStep, const Real& interfaceDofs, const Real& coefficientFirstDerivative,
    			 const mapPtr_Type& interfaceMap, const FESpacePtr_Type& fluidVelocityFESpace,
    			 const FESpacePtr_Type& structureDisplacementFESpace, const vectorPtr_Type& numerationInterface);

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
    matrixPtr_Type M_fluidVelocityToLambda; // I need to test if it faster to assemble or to transpose M_lambdaToFluidVelocity
    matrixPtr_Type M_structureDisplacementToLambda; // I need to test if it faster to assemble or to transpose M_lambdaToStructureDisplacement
    matrixPtr_Type M_structureDisplacementToFluidDisplacement;

    FESpacePtr_Type M_fluidVelocityFESpace;
    FESpacePtr_Type M_structureDisplacementFESpace;

    mapPtr_Type M_interfaceMap;

    vectorPtr_Type M_numerationInterface;

    Real M_interface;
    Real M_timeStep;
    Real M_coefficientFirstDerivative;
};

} // end namespace LifeV

#endif // end FSICOUPLINGCE_H
