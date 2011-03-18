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
 *  @brief File containing the Multiscale Windkessel 0D
 *
 *  @date 08-02-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @author Mahmoud Jafargholi <mahmoud.jafargholi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include "MultiscaleModelWindkessel0D.hpp"

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelWindkessel0D::MultiscaleModelWindkessel0D() :
        multiscaleModel_Type           (),
        M_outputFile                   (),
        M_bc                           ( new bcInterface_Type() ),
        M_pressure_tn                  (),
        M_flowRate_tn                  (),
        M_pressure                     (),
        M_flowRate                     (),
        M_nIntegration                 (),
        M_resistance1                  (),
        M_resistance2                  (),
        M_capacitance                  (),
        M_venousPressure               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::MultiscaleModelWindkessel0D() \n";
#endif

    M_type = Windkessel0D;
}

// ===================================================
// Multiscale PhysicalModel Virtual Methods
// ===================================================
void
MultiscaleModelWindkessel0D::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData( fileName );

    GetPot dataFile( fileName );

    M_nIntegration      = dataFile( "Coefficients/IntegrationSteps", 1   );

    M_resistance1       = dataFile( "Coefficients/Resistance1"     , 1.0 );
    M_resistance2       = dataFile( "Coefficients/Resistance2"     , 1.0 );
    M_capacitance       = dataFile( "Coefficients/Capacitance"     , 1.0 );
    M_venousPressure    = dataFile( "Coefficients/VenousPressure"  , 0.0 );

    setupExporterImporter();
}

void
MultiscaleModelWindkessel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupProblem() \n";
#endif

    initializeSolution();
}

void
MultiscaleModelWindkessel0D::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::buildSystem() \n";
#endif

    M_pressure_tn   = M_pressure;
    M_flowRate_tn   = M_flowRate;
}

void
MultiscaleModelWindkessel0D::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::updateSystem() \n";
#endif

    M_pressure_tn  = M_pressure;
    M_flowRate_tn  = M_flowRate;
}

void
MultiscaleModelWindkessel0D::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveSystem() \n";
#endif

    displayModelStatus( "Solve" );

    switch ( M_bc->handler()->bcType() )
    {
    case OneDimensional::Q:

        M_flowRate = M_bc->handler()->evaluate( M_globalData->dataTime()->time() );
        M_pressure = solveForPressure();

        break;

    case OneDimensional::P:

        M_pressure = M_bc->handler()->evaluate( M_globalData->dataTime()->time() );
        M_flowRate = solveForFlowRate();

        break;

    case OneDimensional::S:

        M_pressure = -M_bc->handler()->evaluate( M_globalData->dataTime()->time() );
        M_flowRate = solveForFlowRate();

        break;

    default:

        std::cout << "Warning: bcType \"" << M_bc->handler()->bcType() << "\"not available!" << std::endl;
    }
}

void
MultiscaleModelWindkessel0D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::saveSolution() \n";
#endif

    M_outputFile << "  " << M_globalData->dataTime()->time()
                 << "    " << M_flowRate
                 << "    " << M_pressure << std::endl;

    if ( M_globalData->dataTime()->isLastTimeStep() )
        M_outputFile.close();
}

void
MultiscaleModelWindkessel0D::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleModel_Type::showMe();

        std::cout << "Integration Steps   = " << M_nIntegration << std::endl << std::endl;

        std::cout << "Resistance1         = " << M_resistance1 << std::endl
                  << "Resistance2         = " << M_resistance2 << std::endl
                  << "Capacitance         = " << M_capacitance << std::endl
                  << "Venous Pressure     = " << M_venousPressure << std::endl << std::endl;
    }
}

// ===================================================
// Methods
// ===================================================
void
MultiscaleModelWindkessel0D::setupLinearModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupLinearModel( ) \n";
#endif

}

void
MultiscaleModelWindkessel0D::updateLinearModel()
{

}

void
MultiscaleModelWindkessel0D::solveLinearModel( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
        return;

    //Solve the linear problem
    displayModelStatus( "Solve linear" );


    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

// ===================================================
// Get Methods (couplings)
// ===================================================
Real
MultiscaleModelWindkessel0D::boundaryStress( const bcFlag_Type& flag, const stress_Type& stressType ) const
{
    switch ( stressType )
    {
    case Pressure:
    {
        return -boundaryPressure( flag );
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return 0;
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaPressure( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return 0;
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem, const stress_Type& stressType )
{
    switch ( stressType )
    {
    case Pressure:
    {
        return -boundaryDeltaPressure( flag, solveLinearSystem );
    }

    default:
    {
        std::cout << "ERROR: Invalid stress type [" << enum2String( stressType, multiscaleStressesMap ) << "]" << std::endl;

        return 0.0;
    }
    }
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelWindkessel0D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        ;// Importer here
    }
    else
    {
        M_pressure = 0;
        M_flowRate = 0;

        switch ( M_bc->handler()->bcType() )
        {
        case OneDimensional::Q:

            M_flowRate = M_bc->handler()->evaluate( M_globalData->dataTime()->time() );

            break;

        case OneDimensional::P:

            M_pressure = M_bc->handler()->evaluate( M_globalData->dataTime()->time() );

            break;

        case OneDimensional::S:

            M_pressure = -M_bc->handler()->evaluate( M_globalData->dataTime()->time() );

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bcType() << "\"not available!" << std::endl;
        }
    }
}

Real
MultiscaleModelWindkessel0D::solveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveForPressure() \n";
#endif

    Real dt     = M_globalData->dataTime()->timeStep();
    Real dti    = dt / M_nIntegration;


    Real RC     = 1.0 / ( M_resistance2 * M_capacitance );  //Dummy Variable
    Real mu     = exp( RC * dt);                            //Dummy Variable

    Real K      = M_pressure_tn / mu;                       //Dummy Variable
    Real dQ     = (M_flowRate - M_flowRate_tn) / dt;

    Real rhs    = 0;            //Dummy Variable (Right Hand Side)
    Real t      = 0;            //Time Integration
    Real integrationSum      = 0;            //sum of Integration

    //%------Integration in dt-------------------------------------------------------------

    rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_venousPressure, dt);
    integrationSum   = dti / 2.0 * rhs;

    //for ( Real t=dti; t<dti+dt; t+=dti )
    for ( UInt i(1) ; i<(M_nIntegration) ; ++i )
    {
        t   = t + dti;
        rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_venousPressure, dt);
        integrationSum   = integrationSum + dti * rhs ;
    }
    t   = t + dti;
    rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_venousPressure, dt);
    integrationSum   = integrationSum + dti / 2.0 * rhs;
    return (K+integrationSum);


}

Real
MultiscaleModelWindkessel0D::solveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveForFlowRate() \n";
#endif
    Real dt     = M_globalData->dataTime()->timeStep();
    Real dti    = dt / M_nIntegration;


    Real RC     = 1.0 / (M_resistance1 * M_resistance2 * M_capacitance);  //Dummy Variable
    Real mu     = exp( ( M_resistance1 + M_resistance2 ) * RC * dt);      //Dummy Variable


    Real K      = M_flowRate_tn / mu;                                     //Dummy Variable
    Real dP     = (M_pressure -M_pressure_tn)/ dt;                        //Dummy Variable


    Real rhs    = 0;        //Dummy Variable (Right Hand Side)
    Real t      = 0;        //Time Integration
    Real integrationSum      = 0;        //sum of Integration

    //%------Integration in dt-------------------------------------------------------------
    t           = 0;

    rhs         = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
    integrationSum           = dti / 2.0 * rhs;
    for ( UInt i(1) ; i<(M_nIntegration) ; ++i )
    {
        t   = t + dti;
        rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
        integrationSum   = integrationSum + dti * rhs ;
    }
    t   = t + dti;
    rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
    integrationSum   = integrationSum + dti / 2.0 * rhs;


    return (K+integrationSum);

}

void
MultiscaleModelWindkessel0D::setupExporterImporter()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupExporterImporter() \n";
#endif

    std::string file = multiscaleProblemFolder + "/Step_" + number2string( multiscaleProblemStep ) + "_Model_" + number2string( M_ID ) + ".m";
    M_outputFile.open( file.c_str(), std::ios::trunc );
    M_outputFile << std::scientific << std::setprecision( 15 )
                 << "% TIME                     FLOW RATE                PRESSURE" << std::endl;
}


} // Namespace multiscale
} // Namespace LifeV
