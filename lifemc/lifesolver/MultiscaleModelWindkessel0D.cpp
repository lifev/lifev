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
        MultiscaleInterfaceFluid       (),
        M_outputFile                   (),
        M_bc                           ( new bcInterface_Type() ),
        M_pressure_tn                  (),
        M_flowRate_tn                  (),
        M_pressure                     (),
        M_flowRate                     (),
        M_tangentPressure              (),
        M_tangentFlowRate              (),
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
// MultiscaleModel Methods
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

    // We need to create the BCHandler before using it
    M_bc->createHandler();
}

void
MultiscaleModelWindkessel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupModel() \n";
#endif

    initializeSolution();
}

void
MultiscaleModelWindkessel0D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::buildModel() \n";
#endif

    M_pressure_tn   = M_pressure;
    M_flowRate_tn   = M_flowRate;
}

void
MultiscaleModelWindkessel0D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::updateModel() \n";
#endif

    M_pressure_tn  = M_pressure;
    M_flowRate_tn  = M_flowRate;
}

void
MultiscaleModelWindkessel0D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveModel() \n";
#endif

    displayModelStatus( "Solve" );

    switch ( M_bc->handler()->bc( OneDimensional::left ).bcType() )
    {
    case OneDimensional::Q:

        M_flowRate = M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );
        M_pressure = solveForPressure();

        break;

    case OneDimensional::P:

        M_pressure = M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );
        M_flowRate = solveForFlowRate();

        break;

    case OneDimensional::S:

        M_pressure = -M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );
        M_flowRate = solveForFlowRate();

        break;

    default:

        std::cout << "Warning: bcType \"" << M_bc->handler()->bc( OneDimensional::left ).bcType() << "\"not available!" << std::endl;
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
// MultiscaleInterfaceFluid Methods
// ===================================================
void
MultiscaleModelWindkessel0D::imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function ) const
{
    M_bc->handler()->setBC( flagConverter( flag ), OneDimensional::Q, boost::bind( function, _1, _1, _1, _1, _1 ) );
}

void
MultiscaleModelWindkessel0D::imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function ) const
{
    M_bc->handler()->setBC( flagConverter( flag ), OneDimensional::S, boost::bind( function, _1, _1, _1, _1, _1 ) );
}

Real
MultiscaleModelWindkessel0D::boundaryStress( const bcFlag_Type& flag ) const
{
    return -boundaryPressure( flag );
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return M_tangentFlowRate;
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem )
{
    solveLinearModel( solveLinearSystem );

    return -M_tangentPressure;
}

// ===================================================
// Private Methods
// ===================================================
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

void
MultiscaleModelWindkessel0D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        ;// TODO: Code the importer
    }
    else
    {
        //TODO Reference pressure
        M_pressure = 0;
        M_flowRate = 0;

        switch ( M_bc->handler()->bc( OneDimensional::left ).bcType() )
        {
        case OneDimensional::Q:

            M_flowRate = M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );

            break;

        case OneDimensional::P:

            M_pressure = M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );

            break;

        case OneDimensional::S:

            M_pressure = -M_bc->handler()->bc( OneDimensional::left ).evaluate( M_globalData->dataTime()->time() );

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bc( OneDimensional::left ).bcType() << "\"not available!" << std::endl;
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
    for ( UInt i(1) ; i< M_nIntegration  ; ++i )
    {
        t   = t + dti;
        rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_venousPressure, dt);
        integrationSum   = integrationSum + dti * rhs;
    }

    t   = t + dti;
    rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_venousPressure, dt);
    integrationSum   = integrationSum + dti / 2.0 * rhs;

    return K + integrationSum;
}

Real
MultiscaleModelWindkessel0D::solveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveForFlowRate() \n";
#endif

    // TODO There is a lot of code duplication with method solveForPressure()

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
    rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
    integrationSum = dti / 2.0 * rhs;

    for ( UInt i(1) ; i< M_nIntegration  ; ++i )
    {
        t   = t + dti;
        rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
        integrationSum   = integrationSum + dti * rhs;
    }

    t   = t + dti;
    rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_venousPressure, dt);
    integrationSum   = integrationSum + dti / 2.0 * rhs;


    return K + integrationSum;
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
    switch ( M_bc->handler()->bc( OneDimensional::left ).bcType() )
    {
    case OneDimensional::Q:

        M_tangentFlowRate = 1.;
        //M_tangentPressure = tangentSolveForPressure(); //TODO Code this method

        break;

    case OneDimensional::P:
    case OneDimensional::S:

        M_tangentPressure = 1.;
        //M_tangentFlowRate = tangentSolveForFlowRate(); //TODO Code this method

        break;

    default:

        std::cout << "Warning: bcType \"" << M_bc->handler()->bc( OneDimensional::left ).bcType() << "\"not available!" << std::endl;
    }

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

} // Namespace multiscale
} // Namespace LifeV
