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
 *  @brief File containing the MultiScale Aitken Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleAlgorithmAitken.hpp>

namespace LifeV
{
namespace multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleAlgorithmAitken::MultiscaleAlgorithmAitken() :
        multiscaleAlgorithm_Type   (),
        M_methodMap                (),
        M_method                   (),
        M_generalizedAitken        ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MultiscaleAlgorithmAitken::MultiscaleAlgorithmAitken() \n";
#endif

    M_type = Aitken;

    //Set Map
    M_methodMap["Scalar"]  = Scalar;
    M_methodMap["Vectorial"] = Vectorial;
    //M_methodMap["VectorialBlock"] = VectorialBlock;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmAitken::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MultiscaleAlgorithmAitken::SetupData( algorithm ) \n";
#endif

    multiscaleAlgorithm_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_generalizedAitken.setDefaultOmega( dataFile( "Solver/Algorithm/Aitken_method/Omega", 1.e-3 ) );
    M_generalizedAitken.useDefaultOmega( dataFile( "Solver/Algorithm/Aitken_method/fixedOmega",   false ) );
    M_generalizedAitken.setOmegaMin( dataFile( "Solver/Algorithm/Aitken_method/range", M_generalizedAitken.defaultOmegaFluid()/1024, 0 ) );
    M_generalizedAitken.setOmegaMax( dataFile( "Solver/Algorithm/Aitken_method/range", M_generalizedAitken.defaultOmegaFluid()*1024, 1 ) );
    M_generalizedAitken.setMinimizationType( dataFile( "Solver/Algorithm/Aitken_method/inverseOmega", true ) );
    M_method = M_methodMap[ dataFile( "Solver/Algorithm/Aitken_method/method", "Vectorial" ) ];
}

void
MultiscaleAlgorithmAitken::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MultiscaleAlgorithmAitken::SubIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( toleranceSatisfied() )
    {
        save( 0, M_couplingResiduals->norm2() );
        return;
    }

    M_multiscale->exportCouplingVariables( *M_couplingVariables );

    M_generalizedAitken.restart();

    // Temporary Computation of a Block Vector - Testing purpose
    //VectorType blocksVector( M_couplingVariables ); blocksVector = 0.0;
    //for ( UInt i = 1 ; i < blocksVector.size() ; i = i+2)
    //    blocksVector[i] = 1.0;
    //std::cout << "blocksVector: " << std::endl;
    //blocksVector.showMe();

    for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
        // To be moved in a post-processing class
        //std::cout << " MS-  CouplingVariables:\n" << std::endl;
        //M_couplingVariables->showMe();
        //std::cout << " MS-  CouplingResiduals:\n" << std::endl;
        //M_couplingResiduals->showMe();

        // Update Coupling Variables
        switch ( M_method )
        {
        case Scalar:

            *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaScalar( *M_couplingVariables, *M_couplingResiduals );

            break;

        case Vectorial:

            *M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVector( *M_couplingVariables, *M_couplingResiduals, true );

            break;

        case VectorialBlock:

            //*M_couplingVariables += M_generalizedAitken.computeDeltaLambdaVectorBlock( *M_couplingVariables, *M_couplingResiduals, blocksVector, 2 );

            break;
        }

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

    multiscaleErrorCheck( Tolerance, "Aitken algorithm residual: " + number2string( M_couplingResiduals->norm2() ) +
                        " (required: " + number2string( M_tolerance ) + ")\n" );
}

void
MultiscaleAlgorithmAitken::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleAlgorithm_Type::showMe();

        std::cout << "Aitken Method       = " << Enum2String( M_method, M_methodMap ) << std::endl;
        std::cout << std::endl << std::endl;
    }
}

} // Namespace multiscale
} // Namespace LifeV
