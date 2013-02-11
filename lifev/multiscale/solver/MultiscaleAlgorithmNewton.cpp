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
 *  @file
 *  @brief File containing the Multiscale Newton Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleAlgorithmNewton.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmNewton::MultiscaleAlgorithmNewton() :
    multiscaleAlgorithm_Type   (),
    M_solver                   (),
    M_jacobian                 ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8013 ) << "MultiscaleAlgorithmNewton::MultiscaleAlgorithmNewton() \n";
#endif

    M_type = Newton;
}

// ===================================================
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmNewton::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8013 ) << "MultiscaleAlgorithmNewton::setupData( fileName ) \n";
#endif

    // Read parameters
    multiscaleParameterListPtr_Type solverParametersList = Teuchos::rcp ( new Teuchos::ParameterList );
    solverParametersList = Teuchos::getParametersFromXmlFile ( fileName );

    // Set main parameters
    setAlgorithmName ( solverParametersList->sublist ( "Multiscale", true, "" ) );
    setAlgorithmParameters ( solverParametersList->sublist ( "Multiscale Algorithm", true, "" ) );

    // Set main parameters
    M_solver.setCommunicator ( M_comm );
    M_solver.setParameters ( solverParametersList->sublist ( "Linear Solver List", true, "" ) );
}

void
MultiscaleAlgorithmNewton::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8013 ) << "MultiscaleAlgorithmNewton::subIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( checkResidual ( 0 ) )
    {
        return;
    }

    M_multiscale->exportCouplingVariables ( *M_couplingVariables );

    multiscaleVectorPtr_Type delta ( new multiscaleVector_Type ( *M_couplingResiduals, Unique ) );
    *delta = 0.0;

    for ( UInt subIT (1); subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
        // Compute the Jacobian
        assembleJacobianMatrix();

        // Set matrix and RHS
        M_solver.setOperator ( M_jacobian );
        M_solver.setRightHandSide ( M_couplingResiduals );

        // Solve Newton (without changing the sign of the residual)
        M_solver.solve ( delta );

        // Changing the sign of the solution
        *delta *= -1;

        // Update Coupling Variables using the Newton Method
        *M_couplingVariables += *delta;

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->importCouplingVariables ( *M_couplingVariables );

        // Verify tolerance
        if ( checkResidual ( subIT ) )
        {
            return;
        }
    }

    save ( M_subiterationsMaximumNumber, M_couplingResiduals->norm2() );

    multiscaleErrorCheck ( Tolerance, "Newton algorithm residual: " + number2string ( M_couplingResiduals->norm2() ) +
                           " (required: " + number2string ( M_tolerance ) + ")\n", M_multiscale->communicator() == 0 );
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleAlgorithmNewton::assembleJacobianMatrix()
{
    // Compute the Jacobian matrix (we completely delete the previous matrix)
    M_jacobian.reset ( new multiscaleMatrix_Type ( M_couplingVariables->map(), 50 ) );
    M_multiscale->exportJacobian ( *M_jacobian );
    M_jacobian->globalAssemble();

    //M_jacobian->spy( multiscaleProblemFolder + multiscaleProblemPrefix + "_AlgorithmJacobianNewtonExported" + "_" + number2string( multiscaleProblemStep ) + "_" + number2string( M_multiscale->globalData()->dataTime()->timeStepNumber() ) );
}

} // Namespace Multiscale
} // Namespace LifeV
