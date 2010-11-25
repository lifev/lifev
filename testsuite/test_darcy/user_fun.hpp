/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
       Date: 2010-07-29

  Copyright (C) 2010 EPFL, Politecnico di Milano

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/

/**
   @file user_fun.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/

#ifndef __user_fun_H
#define __user_fun_H 1

#include <vector>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>
#include <life/lifearray/tab.hpp>

namespace dataProblem
{

using namespace LifeV;

// ===================================================
//!                    Problem data
// ===================================================

// Inverse of permeability matrix
Matrix inversePermeability( const Real&, const Real&, const Real&, const Real&, const std::vector<Real>& );

// Source term
Real source_in( const Real&, const Real&, const Real&, const Real&, const ID& );

// Initial time primal variable for transient and non-linear transient solvers
Real initialPrimal( const Real&, const Real&, const Real&, const Real&, const ID& );

// Zero iteration primal variable for non-linear solver
Real zeroItarationPrimal( const Real&, const Real&, const Real&, const Real&, const ID& );

// Mass function for time dependent problem
Real mass( const Real&, const Real&, const Real&, const Real&, const ID& );

// ===================================================
//!                    Boundary data
// ===================================================

// Boundary condition of Dirichlet
Real dirichlet( const Real&, const Real&, const Real&, const Real&, const ID& );

// Boundary condition of Neumann
Real neumann1( const Real&, const Real&, const Real&, const Real&, const ID& );

Real neumann2( const Real&, const Real&, const Real&, const Real&, const ID& );

// Boundary condition of Robin
Real mixte( const Real&, const Real&, const Real&, const Real&, const ID& );

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution( const Real&, const Real&, const Real&, const Real&, const ID& );

// Gradient of the analytical solution
Real analyticalFlux( const Real&, const Real&, const Real&, const Real&, const ID& );

}

#endif /* __user_fun_H */
