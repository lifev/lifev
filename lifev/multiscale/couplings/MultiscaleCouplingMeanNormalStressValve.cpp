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
 *  @brief File containing the multiscale mean normal stress coupling class with simple valve
 *
 *  @date 05-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/couplings/MultiscaleCouplingMeanNormalStressValve.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingMeanNormalStressValve::MultiscaleCouplingMeanNormalStressValve() :
    multiscaleCoupling_Type     (),
    super_Type                  (),
    M_valveIsOpen               ( true ),
    M_topologyChange            ( false )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::MultiscaleCouplingMeanNormalStressValve() \n";
#endif

    M_type = MeanNormalStressValve;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanNormalStressValve::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::setupCoupling() \n";
#endif

    super_Type::setupCoupling();

    // Preliminary checks to use the valve
    if ( myModelsNumber() > 0 )
    {
        if ( modelsNumber() > 2 )
        {
            std::cout << "!!! WARNING: MultiscaleCouplingMeanNormalStressValve does not work with more than two models !!!" << std::endl;
        }
        if ( M_flowRateInterfaces != static_cast<Int>(modelsNumber()) )
        {
            std::cout << "!!! WARNING: MultiscaleCouplingMeanNormalStressValve does not work with stress boundary data !!!" << std::endl;
        }
    }
}

void
MultiscaleCouplingMeanNormalStressValve::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::initializeCouplingVariables() \n";
#endif

    super_Type::initializeCouplingVariables();

    // We determine the initial position of the valve on the leader process of model 0,
    // then we share the information with the other processes.

    Int localValvePosition ( 0 );
    Int globalValvePosition ( 0 );
    if ( myModel ( 0 ) )
        if ( isModelLeaderProcess ( 0 ) )
        {
            if ( localCouplingVariables ( 0 ) [M_flowRateInterfaces] <= 1e-10 )
            {
                std::cout << " MS-  Valve closed at coupling " << M_ID << std::endl;
                localValvePosition = 0;
            }
            else
            {
                std::cout << " MS-  Valve open at coupling " << M_ID << std::endl;
                localValvePosition = 1;
            }
        }

    // We use the SumAll() instead of the Broadcast() because this way we don't need the id of the leader process.
    M_comm->SumAll ( &localValvePosition, &globalValvePosition, 1 );

    // If the valve is closed the coupling variables are set to zero
    if ( globalValvePosition == 0 )
    {
        M_valveIsOpen = false;
        if ( myModelsNumber() > 0 )
            for ( Int i ( 0 ); i < M_flowRateInterfaces; ++i ) // Only the flow rate is set to zero
            {
                localCouplingVariables ( 0 ) [i] = 0;
            }
    }
    else
    {
        M_valveIsOpen = true;
    }
}

void
MultiscaleCouplingMeanNormalStressValve::updateCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::updateCoupling() \n";
#endif

    super_Type::updateCoupling();

    // We determine if there is a topology change on the leader process of model 0,
    // then we share the information with the other processes.
    Int localTopology ( 0 );
    Int globalTopology ( 0 );
    if ( M_valveIsOpen )
    {
        if ( myModel ( 0 ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterface > ( M_models[0] )->boundaryFlowRate ( M_boundaryIDs[0] );
            if ( isModelLeaderProcess ( 0 ) && myValue < 0 )
            {
                std::cout << " MS-  Closing the valve at coupling " << M_ID << std::endl;
                localTopology = 1;
            }
        }
    }
    else
    {
        Real globalSum ( 0 );
        Real localSum ( 0 );

        if ( myModel ( 1 ) )
        {
            Real myValue = multiscaleDynamicCast< MultiscaleInterface > ( M_models[1] )->boundaryMeanNormalStress ( M_boundaryIDs[1] );
            if ( isModelLeaderProcess ( 1 ) )
            {
                localSum = myValue;
            }
        }

        // We use the SumAll() instead of the Broadcast() because this way we don't need the id of the leader process.
        M_comm->SumAll ( &localSum, &globalSum, 1 );

        if ( myModel ( 0 ) )
        {
            Real myValue = globalSum - multiscaleDynamicCast< MultiscaleInterface > ( M_models[0] )->boundaryMeanNormalStress ( M_boundaryIDs[0] );
            if ( isModelLeaderProcess ( 0 ) && myValue > 0 )
            {
                std::cout << " MS-  Opening the valve at coupling " << M_ID << std::endl;
                localTopology = 1;
            }
        }
    }

    // We use the SumAll() instead of the Broadcast() because this way we don't need the id of the leader process.
    M_comm->SumAll ( &localTopology, &globalTopology, 1 );

    if ( globalTopology > 0 )
    {
        M_topologyChange = true;
        M_valveIsOpen = !M_valveIsOpen;

        // Reset coupling variable history
        resetCouplingHistory();
    }
    else
    {
        M_topologyChange = false;
    }
}

void
MultiscaleCouplingMeanNormalStressValve::computeCouplingResiduals()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::computeCouplingResiduals()  \n";
#endif

    if ( M_valveIsOpen )
    {
        super_Type::computeCouplingResiduals();
    }
    else
    {
        *M_localCouplingResiduals = 0.;
    }
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
void
MultiscaleCouplingMeanNormalStressValve::insertJacobianConstantCoefficients ( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    if ( M_valveIsOpen )
    {
        super_Type::insertJacobianConstantCoefficients ( jacobian );
    }
    else
    {
        // The constant coefficients are added by the leader process of model 0.
        if ( myModel ( 0 ) )
            if ( isModelLeaderProcess ( 0 ) )
            {
                UInt row    = M_couplingVariablesOffset;
                UInt column = M_couplingVariablesOffset;

                for ( Int i ( 0 ); i < M_flowRateInterfaces + 1; ++i )
                {
                    jacobian.addToCoefficient ( row + i, column + i, 1 );
                }
            }
    }
}

void
MultiscaleCouplingMeanNormalStressValve::insertJacobianDeltaCoefficients ( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8240 ) << "MultiscaleCouplingMeanNormalStressValve::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    if ( M_valveIsOpen )
    {
        super_Type::insertJacobianDeltaCoefficients ( jacobian, column, ID, solveLinearSystem );
    }
}

} // Namespace Multiscale
} // Namespace LifeV
