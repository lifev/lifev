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

/**
   @file
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/

#ifndef __user_fun_H
#define __user_fun_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/darcy/solver/DarcySolverTransientNonLinear.hpp>

namespace dataProblem
{

using namespace LifeV;

// ===================================================
//!                   Typedef
// ===================================================

typedef LinearTriangle geoElement_Type;
typedef RegionMesh < geoElement_Type > regionMesh_Type;
typedef std::shared_ptr < regionMesh_Type > regionMeshPtr_Type;

typedef MeshPartitioner < regionMesh_Type > meshPartitioner_Type;

typedef DarcySolverTransientNonLinear < regionMesh_Type > darcySolver_Type;

typedef darcySolver_Type::darcySolverLinear_Type darcySolverLinear_Type;

typedef darcySolver_Type::data_Type darcyData_Type;
typedef darcySolver_Type::dataPtr_Type darcyDataPtr_Type;

typedef darcySolverLinear_Type::fESpace_Type fESpace_Type;
typedef darcySolverLinear_Type::fESpacePtr_Type fESpacePtr_Type;

typedef darcySolverLinear_Type::bcHandler_Type bcHandler_Type;
typedef darcySolverLinear_Type::bcHandlerPtr_Type bcHandlerPtr_Type;

typedef darcySolverLinear_Type::vectorField_Type vectorField_Type;
typedef darcySolverLinear_Type::vectorFieldPtr_Type vectorFieldPtr_Type;

typedef darcySolverLinear_Type::scalarField_Type scalarField_Type;
typedef darcySolverLinear_Type::scalarFieldPtr_Type scalarFieldPtr_Type;

typedef darcySolverLinear_Type::vectorFct_Type vectorFct_Type;
typedef darcySolverLinear_Type::vectorFctPtr_Type vectorFctPtr_Type;

typedef darcySolverLinear_Type::matrixFct_Type matrixFct_Type;
typedef darcySolverLinear_Type::matrixFctPtr_Type matrixFctPtr_Type;

typedef darcySolverLinear_Type::scalarFct_Type scalarFct_Type;
typedef darcySolverLinear_Type::scalarFctPtr_Type scalarFctPtr_Type;

typedef darcySolverLinear_Type::vector_Type vector_Type;
typedef darcySolverLinear_Type::vectorPtr_Type vectorPtr_Type;

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
class inversePermeability : public matrixFct_Type
{
public:
    virtual Matrix eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// Reaction term
class reactionTerm : public scalarFct_Type
{
public:
    virtual Real eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// Scalar source term
class scalarSource : public scalarFct_Type
{
public:
    virtual Real eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// Vector source term
class vectorSource : public vectorFct_Type
{
public:
    virtual Vector eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// Initial time primal variable for transient and non-linear transient solvers
class initialCondition : public scalarFct_Type
{
public:
    virtual Real eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// Mass function for time dependent problem
class massFunction : public scalarFct_Type
{
public:
    virtual Real eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};

// ===================================================
//!                    Boundary data
// ===================================================
namespace BCFlags
{
// Flags for structured meshes
const UInt BOTTOM = 1;
const UInt LEFT   = 2;
const UInt TOP    = 3;
const UInt RIGHT  = 4;
}

void setBoundaryConditions ( bcHandlerPtr_Type& bcDarcy) ;

// Boundary condition of Dirichlet
Real dirichlet ( const Real&, const Real&, const Real&, const Real&, const ID& );

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution ( const Real&, const Real&, const Real&, const Real&, const ID& );

// Gradient of the analytical solution
Real analyticalFlux ( const Real&, const Real&, const Real&, const Real&, const ID& );

} // namespace dataProblem

#endif /* __user_fun_H */
