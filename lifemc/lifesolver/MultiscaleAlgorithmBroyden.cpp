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
 *  @brief File containing the Multiscale Broyden Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleAlgorithmBroyden.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmBroyden::MultiscaleAlgorithmBroyden() :
        multiscaleAlgorithm_Type     (),
        M_solver                     (),
        M_jacobian                   (),
        M_initializeAsIdentityMatrix ( false ),
        M_iterationsLimitReached     ( false ),
        M_iterationsLimitForReset    ( 1 ),
        M_orthogonalization          ( false ),
        M_orthogonalizationSize      ( 1 ),
        M_orthogonalizationContainer ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8014 ) << "MultiscaleAlgorithmBroyden::MultiscaleAlgorithmBroyden() \n";
#endif

    M_type = Broyden;
}

// ===================================================
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmBroyden::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8014 ) << "MultiscaleAlgorithmBroyden::setupData( fileName ) \n";
#endif

    multiscaleAlgorithm_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_initializeAsIdentityMatrix = dataFile( "Parameters/initializeAsIdentityMatrix", false );
    M_iterationsLimitForReset = static_cast <UInt> ( M_subiterationsMaximumNumber * dataFile( "Parameters/iterationsLimitForReset", 1.0 ) );

    M_orthogonalization = dataFile( "Parameters/orthogonalization", false );
    M_orthogonalizationSize = dataFile( "Parameters/orthogonalizationSize", 1 );

    M_solver.setCommunicator( M_comm );
    M_solver.setDataFromGetPot( dataFile, "Solver/AztecOO" );
    //M_solver.setupPreconditioner( DataFile, "Solver/Preconditioner" );
}

void
MultiscaleAlgorithmBroyden::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8014 ) << "MultiscaleAlgorithmBroyden::subIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( checkResidual( 0 ) )
        return;

    M_multiscale->exportCouplingVariables( *M_couplingVariables );

    multiscaleVector_Type delta( *M_couplingResiduals );
    delta = 0.0;
    multiscaleVector_Type minusCouplingResidual( *M_couplingResiduals );
    minusCouplingResidual = 0.0;

    for ( UInt subIT(1); subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
        // Compute the Jacobian (we completery delete the previous matrix)
        if ( subIT == 1 )
        {
            if ( M_jacobian.get() == 0 || M_iterationsLimitReached || M_multiscale->topologyChange() )
            {
                assembleJacobianMatrix();
                M_iterationsLimitReached = false;
            }
        }
        else
            broydenJacobianUpdate( delta );

        // To be moved in a post-processing class
        //std::cout << " MS-  CouplingVariables:\n" << std::endl;
        //M_couplingVariables->showMe();
        //std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        //M_couplingResiduals->showMe();

        //Compute delta using -R
        minusCouplingResidual = -( *M_couplingResiduals );

        M_solver.setMatrix( *M_jacobian );
        M_solver.solve( delta, minusCouplingResidual );
        //M_solver.solveSystem( minusCouplingResidual, delta, M_jacobian, false );

        // Update Coupling Variables using the Broyden Method
        *M_couplingVariables += delta;

        //std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        //M_couplingVariables->showMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->importCouplingVariables( *M_couplingVariables );

        // solveModel
        M_multiscale->solveModel();

        if ( subIT >= M_iterationsLimitForReset )
            M_iterationsLimitReached = true;

        // Verify tolerance
        if ( checkResidual( subIT ) )
            return;
    }

    save( M_subiterationsMaximumNumber, M_couplingResiduals->norm2() );

    multiscaleErrorCheck( Tolerance, "Broyden algorithm residual: " + number2string( M_couplingResiduals->norm2() ) +
                        " (required: " + number2string( M_tolerance ) + ")\n" );
}

void
MultiscaleAlgorithmBroyden::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleAlgorithm_Type::showMe();

        std::cout << "Initialize as identity matrix        = " << M_initializeAsIdentityMatrix << std::endl;
        std::cout << "Reset limit                          = " << M_iterationsLimitForReset << std::endl;
        std::cout << "Enable orthogonalization             = " << M_orthogonalization << std::endl;
        std::cout << "Orthogonalization memory size        = " << M_orthogonalizationSize << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleAlgorithmBroyden::assembleJacobianMatrix()
{
    // Compute the Jacobian matrix
    M_jacobian.reset( new multiscaleMatrix_Type( M_couplingVariables->map(), 50 ) );

    if ( M_initializeAsIdentityMatrix )
        M_jacobian->insertValueDiagonal( 1 );
    else
        M_multiscale->exportJacobian( *M_jacobian );

    M_jacobian->globalAssemble();

    //M_jacobian->spy( "Jacobian" );

    // Reset orthogonalization
    M_orthogonalizationContainer.clear();
}

void
MultiscaleAlgorithmBroyden::broydenJacobianUpdate( const multiscaleVector_Type& delta )
{
    M_jacobian->openCrsMatrix();

    // Compute the Broyden update
    if ( M_orthogonalization )
    {
        // Orthogonalize the vector
        multiscaleVector_Type orthogonalization ( delta );
        for ( containerIterator_Type i = M_orthogonalizationContainer.begin(); i != M_orthogonalizationContainer.end() ; ++i )
            orthogonalization -= orthogonalization.dot( *i ) * *i;
        orthogonalization /= orthogonalization.norm2();

        // Update orthogonalization
        orthogonalizationUpdate( delta );

        // Update the Jacobian
        M_jacobian->addDyadicProduct( ( *M_couplingResiduals ) / orthogonalization.dot( delta ), orthogonalization );
    }
    else
    {
        // Update the Jacobian
        //M_jacobian->addDyadicProduct( ( *M_couplingResiduals + minusCouplingResidual - *M_jacobian * delta ) / delta.dot( delta ), delta );
        M_jacobian->addDyadicProduct( ( *M_couplingResiduals ) / delta.dot( delta ), delta );
    }

    M_jacobian->globalAssemble();

    //M_jacobian->spy( "Jacobian" )
}

void
MultiscaleAlgorithmBroyden::orthogonalizationUpdate( const multiscaleVector_Type& delta )
{
    if ( M_orthogonalizationContainer.size() < M_orthogonalizationSize )
        M_orthogonalizationContainer.push_back( delta );
    else
    {
        M_orthogonalizationContainer.pop_front();
        M_orthogonalizationContainer.push_back( delta );
    }
}

} // Namespace Multiscale
} // Namespace LifeV
