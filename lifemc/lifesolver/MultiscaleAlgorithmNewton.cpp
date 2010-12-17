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
 *  @brief File containing the MultiScale Newton Algorithm
 *
 *  @date 26-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleAlgorithmNewton.hpp>

namespace LifeV
{
namespace multiscale
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
    Debug( 8013 ) << "MultiscaleAlgorithmNewton::MultiscaleAlgorithmNewton() \n";
#endif

    M_type = Newton;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmNewton::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8013 ) << "MultiscaleAlgorithmNewton::SetupData( fileName ) \n";
#endif

    multiscaleAlgorithm_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_solver.setCommunicator( M_comm );
    M_solver.setDataFromGetPot( dataFile, "Solver/Algorithm/Newton_method/AztecOO" );
    //M_solver.setUpPrec( DataFile, "Solver/Algorithm/Newton_method/Preconditioner" );
}

void
MultiscaleAlgorithmNewton::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8013 ) << "MultiscaleAlgorithmNewton::SubIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( toleranceSatisfied() )
    {
        save( 0, M_couplingResiduals->norm2() );
        return;
    }

    M_multiscale->exportCouplingVariables( *M_couplingVariables );

    multiscaleVector_Type delta( *M_couplingResiduals );
    delta = 0.0;
    multiscaleVector_Type minusCouplingResidual( *M_couplingResiduals );
    minusCouplingResidual = 0.0;

    for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
        // Compute the Jacobian (we completery delete the previous matrix)
        M_jacobian.reset( new multiscaleMatrix_Type( M_couplingVariables->map(), 50, 0 ) );
        M_multiscale->exportJacobian( *M_jacobian );
        M_jacobian->globalAssemble();
        M_solver.setMatrix( *M_jacobian );

        //M_jacobian->spy( "Jacobian" );

        // To be moved in a post-processing class
        //std::cout << " MS-  CouplingVariables:\n" << std::endl;
        //M_couplingVariables->showMe();
        //std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        //M_couplingResiduals->showMe();

        //Compute delta using -R
        minusCouplingResidual = -( *M_couplingResiduals );

        M_solver.solve( delta, minusCouplingResidual );
        //M_solver.solveSystem( minusCouplingResidual, delta, M_jacobian, false );

        // Update Coupling Variables using the Newton Method
        *M_couplingVariables += delta;

        //std::cout << " MS-  New CouplingVariables:\n" << std::endl;
        //M_couplingVariables->showMe();

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->importCouplingVariables( *M_couplingVariables );

        // solveSystem
        M_multiscale->solveSystem();

        // Display subiteration information
        if ( M_displayer->isLeader() )
            std::cout << " MS-  Sub-iteration n.:                        " << subIT << std::endl;

        // Verify tolerance
        if ( toleranceSatisfied() )
        {
            save( subIT, M_couplingResiduals->norm2() );
            return;
        }
    }

    save( M_subiterationsMaximumNumber, M_couplingResiduals->norm2() );

    multiscaleErrorCheck( Tolerance, "Newton algorithm residual: " + number2string( M_couplingResiduals->norm2() ) +
                        " (required: " + number2string( M_tolerance ) + ")\n" );
}

} // Namespace multiscale
} // Namespace LifeV
