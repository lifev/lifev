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
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author Anna Scotti <anna.scotti@mail.polimi.it>

    @date 2012-03-30
*/

#ifndef __user_fun_H
#define __user_fun_H 1

#include <vector>

#include <lifev/core/LifeV.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace dataProblem
{

using namespace LifeV;

typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::matrix<Real> Matrix;

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
Matrix inversePermeability ( const Real&, const Real&, const Real&, const Real&, const std::vector<Real>& );

// Source term
Real source ( const Real&, const Real&, const Real&, const Real&, const ID& );

// Vector source term
Vector vectorSource ( const Real&, const Real&, const Real&, const Real&, const ID& );

// ===================================================
//!                    Boundary data
// ===================================================

// Boundary condition of Dirichlet
Real dirichlet ( const Real&, const Real&, const Real&, const Real&, const ID& );

// Boundary condition of Neumann
Real neumann1 ( const Real&, const Real&, const Real&, const Real&, const ID& );

Real neumann2 ( const Real&, const Real&, const Real&, const Real&, const ID& );

// Boundary condition of Robin
Real robin ( const Real&, const Real&, const Real&, const Real&, const ID& );

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution ( const Real&, const Real&, const Real&, const Real&, const ID& );

// Gradient of the analytical solution
Real analyticalFlux ( const Real&, const Real&, const Real&, const Real&, const ID& );

} // namespace LifeV

#endif /* __user_fun_H */
