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

    M_resistance1       = dataFile( "Coefficients/Resistance1"     , 1.0 );
    M_resistance2       = dataFile( "Coefficients/Resistance2"     , 1.0 );
    M_capacitance       = dataFile( "Coefficients/Capacitance"     , 1.0 );
    M_venousPressure    = dataFile( "Coefficients/VenousPressure"  , 0.0 );

    if ( M_globalData.get() )
        setupGlobalData( fileName );

    // We need to create the BCHandler before using it
    M_bc->createHandler();

    // Exporter/Importer
    setupExporterImporter();
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
MultiscaleModelWindkessel0D::setupGlobalData( const std::string& fileName )
{
    GetPot dataFile( fileName );

    if ( !dataFile.checkVariable( "Coefficients/VenousPressure" ) )
        M_venousPressure = M_globalData->fluidReferencePressure();
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
        M_pressure = M_globalData->fluidReferencePressure();
        M_flowRate = 0;

/*        switch ( M_bc->handler->bc( OneDimensional::left )->bcType() )
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
        }*/
    }
}

Real
MultiscaleModelWindkessel0D::solveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveForPressure() \n";
#endif


    if (fabs(M_capacitance)<(1.0e-15))
        return (-(M_resistance1+M_resistance2)*M_flowRate+M_venousPressure);
    else
    {
        Real dt     = M_globalData->dataTime()->timeStep();
        Real dQ     = -(M_flowRate - M_flowRate_tn) / dt;
        Real K1     = 1.0 / (M_resistance2 * M_capacitance);                              //Dummy Variable
        Real K2     = (M_resistance1 + M_resistance2) / (M_resistance2 *M_capacitance);   //Dummy Variable

        return (-1.0/(K1*K1)*(K2*dQ-exp(-K1*dt)*(K2*dQ+(K1*K1)*M_pressure_tn-(K1*K1)*M_venousPressure+K1*K2*M_flowRate_tn-K1*M_resistance1*dQ))+(M_resistance1*dQ+K1*M_venousPressure-K2*M_flowRate_tn+K2*dQ*dt)/K1);
    }

}

Real
MultiscaleModelWindkessel0D::solveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveForFlowRate() \n";
#endif

    // TODO There is a lot of code duplication with method solveForPressure()

    if (fabs(M_capacitance)<(1.0e-15))
        return (-(M_pressure-M_venousPressure)/(M_resistance1+M_resistance2));
    else
    {
        Real dt     =   M_globalData->dataTime()->timeStep();
        Real dP     = ( M_pressure    - M_pressure_tn ) / dt;
        Real K1     = ( M_resistance1 + M_resistance2 ) / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K2     =   1.0 / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K3     =   1.0 /   M_resistance1;

        return (-(K3*dP+K2*M_pressure_tn-K2*M_venousPressure+K2*dP*dt)/K1+1.0/(K1*K1)*(K2*dP+exp(-K1*dt)*(-K2*dP+(K1*K1)*M_flowRate_tn+K1*K2*M_pressure_tn-K1*K2*M_venousPressure+K1*K3*dP)));
    }
}

Real
MultiscaleModelWindkessel0D::tangentSolveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::tangentSolveForFlowRate() \n";
#endif

    // TODO There is a lot of code duplication with method solveForPressure()

    if (fabs(M_capacitance)<(1.0e-15))
        return (-1.0/(M_resistance1+M_resistance2));
    else
    {
        Real dt     =   M_globalData->dataTime()->timeStep();
        Real K1     = ( M_resistance1 + M_resistance2 ) / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K2     =   1.0 / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K3     =   1.0 /   M_resistance1;
        return (-(K2+K3/dt)/K1-1.0/(K1*K1)*(exp(-K1*dt)*(K2/dt-(K1*K3)/dt)-K2/dt));
    }

}

Real
MultiscaleModelWindkessel0D::tangentSolveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::tangentSolveForPressure() \n";
#endif

    // TODO There is a lot of code duplication with method solveForPressure()
    if (fabs(M_capacitance)<(1.0e-15))
        return (-(M_resistance1+M_resistance2));
    else
    {
        Real dt     = M_globalData->dataTime()->timeStep();
        Real K1     = 1.0 / (M_resistance2 * M_capacitance);                              //Dummy Variable
        Real K2     = (M_resistance1 + M_resistance2) / (M_resistance2 *M_capacitance);   //Dummy Variable
        return (-(K2+M_resistance1/dt)/K1-1.0/(K1*K1)*(exp(-K1*dt)*(K2/dt-(K1*M_resistance1)/dt)-K2/dt));
    }

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
    case OneDimensional::Q: // dP/dQ

        M_tangentFlowRate = 1.;
        M_tangentPressure = tangentSolveForPressure();

        break;

    case OneDimensional::P: // dQ/dP

        M_tangentPressure = 1.;
        M_tangentFlowRate = tangentSolveForFlowRate();

        break;

    case OneDimensional::S: // dQ/dS

        M_tangentPressure = 1.;
        M_tangentFlowRate = -tangentSolveForFlowRate();

        break;

    default:

        std::cout << "Warning: bcType \"" << M_bc->handler()->bc( OneDimensional::left ).bcType() << "\"not available!" << std::endl;
    }

    //This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

} // Namespace multiscale
} // Namespace LifeV
