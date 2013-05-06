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
    @file AssemblyPolicyNavierStokesSemiImplicit class
    @brief This class contains all the informations necessary
           to assemble a Navier-Stokes problem using a
           semi implicit scheme

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 06-12-2012
 */

#include <lifev/navier_stokes/solver/NavierStokesSolver/AssemblyPolicyNavierStokesSemiImplicit.hpp>

#include <string>
#include <lifev/core/util/LifeChrono.hpp>
// #include <lifev/navier_stokes/algorithm/PreconditionerPCD.hpp>

namespace LifeV
{

void
AssemblyPolicyNavierStokesSemiImplicit::initAssembly ( Teuchos::ParameterList& list )
{
    AssemblyPolicyStokes::initAssembly ( list );

    LifeChrono assemblyChrono;
    assemblyChrono.start();

    displayer().leaderPrint ( "Creating the mass matrix... " );
    map_Type solutionMap ( uFESpace()->map() + pFESpace()->map() );
    M_massMatrix.reset ( new matrix_Type ( solutionMap ) );
    M_massMatrix->zero();
    M_assembler->addMass ( *M_massMatrix, 1.0 );
    M_massMatrix->globalAssemble();
    displayer().leaderPrint ( "done\n" );

    assemblyChrono.stop();
    displayer().leaderPrintMax ("Matrices assembly time: ", assemblyChrono.diff(), " s.\n");
}

void
AssemblyPolicyNavierStokesSemiImplicit::assembleSystem ( matrixPtr_Type systemMatrix,
                                                         vectorPtr_Type rhs,
                                                         vectorPtr_Type solution,
                                                         preconditionerPtr_Type preconditioner )
{
    AssemblyPolicyStokes::assembleSystem ( systemMatrix, rhs, solution, preconditioner );

    bdf()->updateRHSContribution ( timestep() );
    *rhs += *M_massMatrix * bdf()->rhsContributionFirstDerivative();

    double alpha = bdf()->coefficientFirstDerivative ( 0 ) / timestep();
    *systemMatrix += *M_massMatrix * alpha;

    vector_Type beta ( systemMatrix->map(), Repeated );
    bdf()->extrapolation (beta);
    M_assembler->addConvection ( *systemMatrix, 1.0, beta );

//    if ( preconditioner->preconditionerType() == "PCD" )
//    {
//        PreconditionerPCD* pcdPtr = dynamic_cast<PreconditionerPCD*> ( preconditioner.get() );
//        pcdPtr->updateBeta ( beta );
//    }
}

} // namespace LifeV
