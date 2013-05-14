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
 *  @mantainer    Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/models/MultiscaleModelWindkessel0D.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleModelWindkessel0D::MultiscaleModelWindkessel0D() :
    multiscaleModel_Type           (),
    MultiscaleInterface            (),
    M_outputFile                   (),
    M_bc                           ( new bcInterface_Type() ),
    M_data                         ( new data_Type() ),
    M_pressureLeft_tn              (),
    M_flowRateLeft_tn              (),
    M_pressureLeft                 (),
    M_flowRateLeft                 (),
    M_pressureRight                (),
    M_tangentPressureLeft          (),
    M_tangentFlowRateLeft          (),
    M_resistance1                  (),
    M_resistance2                  (),
    M_capacitance                  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::MultiscaleModelWindkessel0D() \n";
#endif

    M_type = Windkessel0D;
}

// ===================================================
// MultiscaleModel Methods
// ===================================================
void
MultiscaleModelWindkessel0D::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::setupData( fileName ) \n";
#endif

    multiscaleModel_Type::setupData ( fileName );

    GetPot dataFile ( fileName );

    // R1, R2, and C can have global scaling factors applied to them, these are overridden by local scaling factors
    Real resistanceScalingFactor = dataFile ( "ScalingFactors/Resistance", M_globalData->scalingFactorResistance() );
    Real complianceScalingFactor = dataFile ( "ScalingFactors/Compliance", M_globalData->scalingFactorCompliance() );

    M_resistance1       = dataFile ( "Coefficients/Resistance1"     , 1.0 ) * resistanceScalingFactor;
    M_resistance2       = dataFile ( "Coefficients/Resistance2"     , 1.0 ) * resistanceScalingFactor;
    M_capacitance       = dataFile ( "Coefficients/Capacitance"     , 1.0 ) * complianceScalingFactor;

    if ( M_globalData.get() )
    {
        setupGlobalData ( fileName );
    }

    // We need to create the BCHandler before using it
    M_bc->createHandler();

    // Exporter/Importer
    setupExporterImporter();
}

void
MultiscaleModelWindkessel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::setupModel() \n";
#endif

    initializeSolution();

    M_bc->setPhysicalSolver ( M_data );

    // Safety check
    if ( M_bc->handler()->bc ( 1 ).bcType() != Voltage )
    {
        std::cout << "!!! Error: the Windkessel model support only stress boundary conditions on the right at the present time !!!" << std::endl;
    }
}

void
MultiscaleModelWindkessel0D::buildModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::buildModel() \n";
#endif

    updateModel();
}

void
MultiscaleModelWindkessel0D::updateModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::updateModel() \n";
#endif

    M_pressureLeft_tn = M_pressureLeft;
    M_flowRateLeft_tn = M_flowRateLeft;

    // Update BCInterface solver variables
    M_bc->updatePhysicalSolverVariables();

    M_pressureRight   = -M_bc->handler()->bc ( 1 ).evaluate ( M_globalData->dataTime()->time() );
}

void
MultiscaleModelWindkessel0D::solveModel()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::solveModel() \n";
#endif

    displayModelStatus ( "Solve" );

    switch ( M_bc->handler()->bc ( 0 ).bcType() )
    {
        case Current:

            M_flowRateLeft = M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );
            M_pressureLeft = solveForPressure();

            break;

        case Voltage:

            M_pressureLeft = -M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );
            M_flowRateLeft = solveForFlowRate();

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

            break;
    }
}

void
MultiscaleModelWindkessel0D::updateSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::updateSolution() \n";
#endif

}

void
MultiscaleModelWindkessel0D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::saveSolution() \n";
#endif

    M_outputFile << "    " << M_globalData->dataTime()->time()
                 << "    " << M_flowRateLeft
                 << "    " << M_pressureLeft << std::endl;

    if ( M_globalData->dataTime()->isLastTimeStep() )
    {
        M_outputFile.close();
    }
}

void
MultiscaleModelWindkessel0D::showMe()
{
    if ( M_comm->MyPID() == 0 )
    {
        multiscaleModel_Type::showMe();

        std::cout << "Resistance1         = " << M_resistance1 << std::endl
                  << "Resistance2         = " << M_resistance2 << std::endl
                  << "Capacitance         = " << M_capacitance << std::endl << std::endl;
    }
}

Real
MultiscaleModelWindkessel0D::checkSolution() const
{
    return M_pressureLeft + M_flowRateLeft;
}

// ===================================================
// MultiscaleInterface Methods
// ===================================================
void
MultiscaleModelWindkessel0D::imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    ZeroDimensionalFunction base;
    base.setFunction ( boost::bind ( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC ( boundaryFlag ( boundaryID ), Current, base );
}

void
MultiscaleModelWindkessel0D::imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function )
{
    ZeroDimensionalFunction base;
    base.setFunction ( boost::bind ( function, _1, _1, _1, _1, _1 ) );

    M_bc->handler()->setBC ( boundaryFlag ( boundaryID ), Voltage, base );
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    if ( boundaryFlag ( boundaryID ) == 1 )
    {
        return 0;
    }

    solveLinearModel ( solveLinearSystem );

    return M_tangentFlowRateLeft;
}

Real
MultiscaleModelWindkessel0D::boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem )
{
    if ( boundaryFlag ( boundaryID ) == 1 )
    {
        return 0;
    }

    solveLinearModel ( solveLinearSystem );

    return -M_tangentPressureLeft;
}

// ===================================================
// Private Methods
// ===================================================
void
MultiscaleModelWindkessel0D::setupGlobalData ( const std::string& fileName )
{
    GetPot dataFile ( fileName );

    //Global data time
    M_data->setTimeData ( M_globalData->dataTime() );

    if ( !dataFile.checkVariable ( "Coefficients/VenousPressure" ) )
    {
        M_pressureRight = M_globalData->fluidVenousPressure();
    }
    M_data->setVenousPressure ( M_pressureRight );
}

void
MultiscaleModelWindkessel0D::initializeSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::initializeSolution() \n";
#endif

    if ( multiscaleProblemStep > 0 )
    {
        std::string fileName = multiscaleProblemFolder + multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep - 1 ) + ".mfile";

        std::ifstream inputFile;
        inputFile.open ( fileName.c_str(), std::ios::in );

        if ( inputFile.is_open() )
        {
            // Define some variables
            std::string line;
            std::vector<std::string> stringsVector;
            Real deltaT (1e15);

            // Read the first line with comments
            std::getline ( inputFile, line, '\n' );

            // Read one-by-one all the others lines of the fileName
            while ( std::getline ( inputFile, line, '\n' ) )
            {
                // Split the three entries
                boost::split ( stringsVector, line, boost::is_any_of ( " " ), boost::token_compress_on );

                // Import values
                if ( std::abs ( string2number ( stringsVector[1] ) - M_globalData->dataTime()->initialTime() ) <= deltaT )
                {
                    deltaT = std::abs ( string2number ( stringsVector[1] ) - M_globalData->dataTime()->initialTime() );

                    M_flowRateLeft = string2number ( stringsVector[2] );
                    M_pressureLeft = string2number ( stringsVector[3] );
                }
            }

            // Close file
            inputFile.close();
        }
        else
        {
            std::cerr << " !!! Error: cannot open fileName: " << fileName.c_str() << " !!!" << std::endl;
        }
    }
    else
    {
        M_flowRateLeft = 0;
        M_pressureLeft = M_globalData->solidExternalPressure();

        switch ( M_bc->handler()->bc ( 0 ).bcType() )
        {
            case Current:

                M_flowRateLeft = M_bc->handler()->bc ( 0 ).evaluate ( M_globalData->dataTime()->time() );

                break;

            case Voltage:

                //            if ( std::abs( M_globalData->solidExternalPressure() + M_bc->handler()->bc( 0 ).evaluate( M_globalData->dataTime()->time() ) ) > 1e-14 )
                //                std::cout << "!!! Warning: external pressure should be equal to the initial pressure !!! " << std::endl;

                break;

            default:

                std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

                break;
        }
    }
}

void
MultiscaleModelWindkessel0D::setupExporterImporter()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::setupExporterImporter() \n";
#endif

    std::string file = multiscaleProblemFolder + multiscaleProblemPrefix + "_Model_" + number2string ( M_ID ) + "_" + number2string ( multiscaleProblemStep ) + ".mfile";
    M_outputFile.open ( file.c_str(), std::ios::trunc );
    M_outputFile << std::scientific << std::setprecision ( 15 )
                 << "%   MODEL: " << M_modelName << std::endl
                 << "%   TIME                     FLOW RATE                PRESSURE" << std::endl;
}

Real
MultiscaleModelWindkessel0D::solveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::solveForFlowRate() \n";
#endif

    if ( std::abs ( M_capacitance ) < 1e-14 )
    {
        return - ( M_pressureLeft - M_pressureRight ) / ( M_resistance1 + M_resistance2 );
    }
    else
    {
        Real dt     =   M_globalData->dataTime()->timeStep();
        Real dP     = - ( M_pressureLeft    - M_pressureLeft_tn ) / dt;
        Real K1     = ( M_resistance1 + M_resistance2 ) / ( M_resistance1 * M_resistance2 * M_capacitance );
        Real K2     =   1.0 / ( M_resistance1 * M_resistance2 * M_capacitance );
        Real K3     =   1.0 / M_resistance1;

        return - 1.0 / ( K1 * K1 )
               * ( K2 * dP - std::exp ( -K1 * dt ) * ( K2 * dP + (K1 * K1) * M_flowRateLeft_tn - K1 * K2 * M_pressureRight + K1 * K2 * M_pressureLeft_tn - K1 * K3 * dP ) )
               + ( K3 * dP - K2 * M_pressureLeft_tn + K2 * M_pressureRight + K2 * dP * dt ) / K1;
    }
}

Real
MultiscaleModelWindkessel0D::solveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::solveForPressure() \n";
#endif

    if ( std::abs ( M_capacitance ) < 1e-14 )
    {
        return - ( M_resistance1 + M_resistance2 ) * M_flowRateLeft + M_pressureRight;
    }
    else
    {
        Real dt     = M_globalData->dataTime()->timeStep();
        Real dQ     = - ( M_flowRateLeft - M_flowRateLeft_tn ) / dt;
        Real K1     = 1.0 / ( M_resistance2 * M_capacitance );
        Real K2     = ( M_resistance1 + M_resistance2 ) / ( M_resistance2 * M_capacitance );
        Real K3     = M_resistance1;

        return - 1.0 / ( K1 * K1 )
               * ( K2 * dQ - std::exp ( -K1 * dt ) * ( K2 * dQ + ( K1 * K1 ) * M_pressureLeft_tn - K1 * K1 * M_pressureRight + K1 * K2 * M_flowRateLeft_tn - K1 * K3 * dQ) )
               + ( K3 * dQ - K2 * M_flowRateLeft_tn + K1 * M_pressureRight + K2 * dQ * dt ) / K1;
    }
}

void
MultiscaleModelWindkessel0D::solveLinearModel ( bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::solveLinearModel() \n";
#endif

    if ( !solveLinearSystem )
    {
        return;
    }

    //Solve the linear problem
    displayModelStatus ( "Solve linear" );
    switch ( M_bc->handler()->bc ( 0 ).bcType() )
    {
        case Current: // dP/dQ

            M_tangentFlowRateLeft = 1.;
            M_tangentPressureLeft = tangentSolveForPressure();

            break;

        case Voltage: // dQ/dS

            M_tangentPressureLeft = 1.;
            M_tangentFlowRateLeft = -tangentSolveForFlowRate();

            break;

        default:

            std::cout << "Warning: bcType \"" << M_bc->handler()->bc ( 0 ).bcType() << "\"not available!" << std::endl;

            break;
    }

    // This flag avoid recomputation of the same system
    solveLinearSystem = false;
}

Real
MultiscaleModelWindkessel0D::tangentSolveForFlowRate()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::tangentSolveForFlowRate() \n";
#endif

    if ( std::abs ( M_capacitance ) < 1e-14 )
    {
        return -1.0 / ( M_resistance1 + M_resistance2 );
    }
    else
    {
        Real dt     =   M_globalData->dataTime()->timeStep();
        Real K1     = ( M_resistance1 + M_resistance2 ) / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K2     =   1.0 / ( M_resistance1 * M_resistance2 * M_capacitance);
        Real K3     =   1.0 /   M_resistance1;

        return - ( K2 + K3 / dt ) / K1
               -  1.0 / ( K1 * K1 ) * ( std::exp ( -K1 * dt ) * ( K2 / dt - ( K1 * K3 ) / dt ) - K2 / dt );
    }
}

Real
MultiscaleModelWindkessel0D::tangentSolveForPressure()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8150 ) << "MultiscaleModelWindkessel0D::tangentSolveForPressure() \n";
#endif

    if ( std::abs ( M_capacitance ) < 1e-14 )
    {
        return - ( M_resistance1 + M_resistance2 );
    }
    else
    {
        Real dt     = M_globalData->dataTime()->timeStep();
        Real K1     = 1.0 / ( M_resistance2 * M_capacitance );
        Real K2     = ( M_resistance1 + M_resistance2 ) / ( M_resistance2 * M_capacitance );
        Real K3     = M_resistance1;

        return - ( K2 + K3 / dt ) / K1
               -  1.0 / ( K1 * K1 ) * ( std::exp ( -K1 * dt ) * ( K2 / dt - ( K1 * K3 ) / dt ) - K2 / dt );
    }
}

} // Namespace multiscale
} // Namespace LifeV
