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
   @file user_fun.cpp
   @author L. Oldani <luca.oldani@mail.polimi.it>
   @date 2014-09-20
*/

#include "user_fun.hpp"

namespace DataProblem
{

// ==================================================
//!                    Problem data
// ===================================================

// Initial time condition
Real initialCondition::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 0.0;
} // eval

// Mass function
Real massTerm::eval ( const UInt& /*iElem*/, const Vector3D& /*P*/, const Real& /*time*/ ) const
{
    return 1.;
} // eval

// ===================================================
//!                    Boundary data
// ===================================================

void setBoundaryConditions ( bcHandlerPtr_Type & bcHyperbolic )
{

    BCFunctionBase dirichletBDfunL;
    dirichletBDfunL.setFunction ( dirichletL );

	BCFunctionBase dirichletBDfunB;
    dirichletBDfunB.setFunction ( dirichletB );
    
	bcHyperbolic->addBC( "Left", Structured2DLabel::LEFT, Essential, Scalar, dirichletBDfunL );
    bcHyperbolic->addBC( "Bottom", Structured2DLabel::BOTTOM, Essential, Scalar, dirichletBDfunB );
} // setBoundaryConditions

// Boundary condition of Dirichlet
Real dirichletL ( const Real& t,
                 const Real& x,
                 const Real& /*y*/,
                 const Real& /*z*/,
                 const ID&   /*ic*/ )
{
    if ( t < 0.25 )
    {	
        return 1;
    }
    else
    {
        return 0.;
    }
} // dirichlet

// Boundary condition of Dirichlet
Real dirichletB ( const Real& t,
                 const Real& x,
                 const Real& /*y*/,
                 const Real& /*z*/,
                 const ID&   /*ic*/ )
{
    if ( t < 0.3 )
    {	
        return 1;
    }
    else
    {
        return 0.;
    }
} // dirichlet

// ===================================================
//!                 Analytical solution
// ===================================================

// Analytical solution
Real analyticalSolution ( const Real& t,
                          const Real& x,
                          const Real& /*y*/,
                          const Real& /* z */,
                          const ID&  /* ic */ )
{
    Real sol;
    const Real tinj=0.5;
	const Real BC=0.1;
	const Real Vm=2.;
	const Real rhom=6.;
    if (x>Vm*t)
    {
        sol = 0;
    }
    else if (t<=tinj)
    {
	if (x<(Vm*t*(1-(2*BC/rhom))))
	{
		sol = BC;
	}
	else 
	{
	        sol = -rhom/(2*Vm*t)*x+rhom/2;
	}
    }
    else if (x>(Vm*t*(1-(2*BC/rhom))))
	{
		sol = -rhom/(2*Vm*t)*x+rhom/2;
	}
    else if (x< rhom/(2*Vm)*(t-tinj))
	{
		sol = 0; 
	}
    else
	{
		sol = BC;
	}	

    return sol;
} // analyticalSolution

} // Namespace DataProblem
