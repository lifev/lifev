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
        multiscaleModel_Type           ()
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

    M_resistance1       = dataFile("Coefficients/Resistance1"       ,1.0);
    M_resistance2       = dataFile("Coefficients/Resistance2"       ,1.0);
    M_capacitance       = dataFile("Coefficients/capacitance"       ,1.0);
    M_nIntegration      = dataFile("Coefficients/IntegrationSteps"  ,1.0);
    M_backPressure      = dataFile("Coefficients/BackPressure"      ,0.0);
    M_solveForPressure  = dataFile("Coefficients/SolveForPressure"  ,1  );


    // Read data from a GetPot file, and store them the model.
    // At this level you don't know which re the BC. You don't know if you re imposing  Q or  P
    // Here reAd only the DATA! FOR NOW, PUT THE BC HERE..
}

void
MultiscaleModelWindkessel0D::setupModel()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupProblem() \n";
#endif

    // AT THIS POINT YOU KNOW THE BC TYPE. SO YOU KNOW WHETHER YOU ARE SOLVING A PROBLEM FOR THE Q OR FOR THE P. YOU CAN INITIALIZE WHAT YOU NEED.
    // HERE YOU INITIALIZE P^n AND Q^n = 0
    //M_pressure_tn   = 0.0;
    M_pressure   = 0.0;
    M_flowRate   = 100000.0;
}

void
MultiscaleModelWindkessel0D::buildSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::buildSystem() \n";
#endif
    // THIS IS CALLED INSIDE THE TIME LOOP, BUT ONLY THE VERY FIRST TIME STEP.
    // MAYBE THIS IS EMPTY FOR THIS MODEL
    M_pressure_tn   = M_pressure;
    M_flowRate_tn   = M_flowRate;
}

void
MultiscaleModelWindkessel0D::updateSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::updateSystem() \n";
#endif
    // THIS IS CALLED INSIDE THE TIME LOOP, BUT ONLY FROM THE SECOND TIME STEP.
    // UPDATE P^n AND Q^n.
    M_pressure_tn  = M_pressure;
    M_flowRate_tn  = M_flowRate;
}

void
MultiscaleModelWindkessel0D::solveSystem()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::solveSystem() \n";
#endif
    // THIS IS CALLED INSIDE THE TIME LOOP
    // HERE YOU USE THE UPDATED VARIABLES TO SOLVE A PROBLEM FROM T^n TO T^n+1
    if (M_solveForPressure)
    {
        //M_flowRate = sin(M_globalData->dataTime()->time());
        M_flowRate = 100000.0*exp(-M_globalData->dataTime()->time());
        M_pressure = solveForPressure();
    }
    else
    {
        M_pressure = sin(M_globalData->dataTime()->time());
        //M_pressure = 100000.0*exp(-M_globalData->dataTime()->time());
        M_flowRate = solveForFlowRate();
    }

}

void
MultiscaleModelWindkessel0D::saveSolution()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::saveSolution() \n";
#endif


    static std::ofstream myfile;
    if (!(M_globalData->dataTime()->timeStepNumber()))
    {
        myfile.open ("myout.txt");
        myfile << "Time \t"<<"Pressure \t"<<"flowRate \n";
    }

    myfile << M_globalData->dataTime()->time()<<"\t"<< M_pressure <<"\t"<<M_flowRate <<"\n";

    if ( M_globalData->dataTime()->isLastTimeStep())
        myfile.close();
    // THIS IS CALLED INSIDE THE TIME LOOP
    // A txt file, with 2 columns [ t^n+1, Q, P ]

}

void
MultiscaleModelWindkessel0D::showMe()
{
    if ( M_displayer->isLeader() )
    {
        multiscaleModel_Type::showMe();
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

    //Solve here

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
    case StaticPressure:
    {
        return boundaryPressure( flag );
    }

    case TotalPressure:
    {
        return boundaryPressure( flag ) + boundaryDynamicPressure( flag ) * ( ( boundaryFlowRate( flag ) > 0.0 ) ? 1 : -1 );
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
    case StaticPressure:
    {
        return boundaryDeltaPressure( flag, solveLinearSystem );
    }

    case TotalPressure:
    {
        return boundaryDeltaPressure( flag, solveLinearSystem ) + boundaryDeltaDynamicPressure( flag, solveLinearSystem ); //Verify the sign of DynamicPressure contribute!
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
MultiscaleModelWindkessel0D::setupGlobalData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8150 ) << "MultiscaleModelWindkessel0D::setupGlobalData( fileName ) \n";
#endif

}

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
        ;// Initialization
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

    rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_backPressure, dt);
    integrationSum   = dti / 2.0 * rhs;

    //for ( Real t=dti; t<dti+dt; t+=dti )
    for ( Int i=1;i<(M_nIntegration);++i)
    {
        t   = t + dti;
        rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_backPressure, dt);
        integrationSum   = integrationSum + dti * rhs ;
    }
    t   = t + dti;
    rhs = rightHandSideP(M_resistance1, M_resistance2, RC, t, M_flowRate_tn, M_flowRate, dQ, M_backPressure, dt);
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
    int i       = 0;
    t           = 0;

    rhs         = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_backPressure, dt);
    integrationSum           = dti / 2.0 * rhs;
    for (i=1;i<(M_nIntegration);++i)
    {
        t   = t + dti;
        rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_backPressure, dt);
        integrationSum   = integrationSum + dti * rhs ;
    }
    t   = t + dti;
    rhs = rightHandSideQ(M_resistance1, M_resistance2, RC, t, M_pressure_tn, M_pressure, dP, M_backPressure, dt);
    integrationSum   = integrationSum + dti / 2.0 * rhs;


    return (K+integrationSum);

}


} // Namespace multiscale
} // Namespace LifeV
