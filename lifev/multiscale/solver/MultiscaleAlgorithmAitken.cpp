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
 *  @brief File containing the Multiscale Aitken Algorithm
 *
 *  @date 23-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleAlgorithmAitken.hpp>

namespace LifeV
{
namespace Multiscale
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
    debugStream( 8011 ) << "MultiscaleAlgorithmAitken::MultiscaleAlgorithmAitken() \n";
#endif

    M_type = Aitken;

    //Set Map
    M_methodMap["Scalar"]  = Scalar;
    M_methodMap["Vectorial"] = Vectorial;
    //M_methodMap["VectorialBlock"] = VectorialBlock;
}

// ===================================================
// Multiscale Algorithm Virtual Methods
// ===================================================
void
MultiscaleAlgorithmAitken::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8011 ) << "MultiscaleAlgorithmAitken::setupData( fileName ) \n";
#endif

    // Read parameters
    multiscaleParameterListPtr_Type solverParametersList = Teuchos::rcp( new Teuchos::ParameterList );
    solverParametersList = Teuchos::getParametersFromXmlFile( fileName );

    // Set main parameters
    setAlgorithmName( solverParametersList->sublist( "Multiscale", true, "" ) );
    setAlgorithmParameters( solverParametersList->sublist( "Multiscale Algorithm", true, "" ) );
}

void
MultiscaleAlgorithmAitken::subIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream( 8011 ) << "MultiscaleAlgorithmAitken::subIterate() \n";
#endif

    multiscaleAlgorithm_Type::subIterate();

    // Verify tolerance
    if ( checkResidual( 0 ) )
        return;

    M_multiscale->exportCouplingVariables( *M_couplingVariables );

    M_generalizedAitken.restart();

    // Temporary Computation of a Block Vector - Testing purpose
//    VectorType blocksVector( M_couplingVariables ); blocksVector = 0.0;
//    for ( UInt i = 1 ; i < blocksVector.size() ; i = i+2)
//        blocksVector[i] = 1.0;
//    std::cout << "blocksVector: " << std::endl;
//    blocksVector.showMe();

    for ( UInt subIT = 1; subIT <= M_subiterationsMaximumNumber; ++subIT )
    {
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

        // Import Coupling Variables inside the coupling blocks
        M_multiscale->importCouplingVariables( *M_couplingVariables );

        // Verify tolerance
        if ( checkResidual( subIT ) )
            return;
    }

    save( M_subiterationsMaximumNumber, M_couplingResiduals->norm2() );

    multiscaleErrorCheck( Tolerance, "Aitken algorithm residual: " + number2string( M_couplingResiduals->norm2() ) +
                        " (required: " + number2string( M_tolerance ) + ")\n", M_multiscale->communicator() == 0  );
}

void
MultiscaleAlgorithmAitken::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleAlgorithm_Type::showMe();

        std::cout << "Aitken Method                        = " << enum2String( M_method, M_methodMap ) << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
void
MultiscaleAlgorithmAitken::setAlgorithmParameters( const multiscaleParameterList_Type& parameterList )
{
    multiscaleAlgorithm_Type::setAlgorithmParameters( parameterList );

    M_generalizedAitken.setDefaultOmega( parameterList.get<Real>( "Omega" ) );
    M_generalizedAitken.useDefaultOmega( parameterList.get<bool>( "Fixed omega" ) );
    M_generalizedAitken.setOmegaMin( parameterList.get<Real>( "Range minimum" ) );
    M_generalizedAitken.setOmegaMax( parameterList.get<Real>( "Range maximum" ) );
    M_generalizedAitken.setMinimizationType( parameterList.get<bool>( "Inverse omega" ) );
    M_method = M_methodMap[ parameterList.get<std::string>( "Method" ) ];
}

} // Namespace Multiscale
} // Namespace LifeV
