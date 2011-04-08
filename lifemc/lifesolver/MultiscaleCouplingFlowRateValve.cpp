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
 *  @brief File containing the Multiscale Coupling FlowRateValve
 *
 *  @date 05-04-2011
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MultiscaleCouplingFlowRateValve.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingFlowRateValve::MultiscaleCouplingFlowRateValve() :
        multiscaleCoupling_Type     (),
        super_Type                  (),
        M_valveIsOpen               ( true ),
        M_topologyChange            ( false )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::MultiscaleCouplingFlowRateValve() \n";
#endif

    M_type = FlowRateValve;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingFlowRateValve::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8240 ) << "MultiscaleCouplingFlowRateValve::setupCoupling() \n";
#endif

    super_Type::setupCoupling();

    if ( M_couplingIndex.first > 2)
        std::cout << "!!! WARNING: MultiscaleCouplingFlowRateValve does not work with more than two models !!!" << std::endl;
}

void
MultiscaleCouplingFlowRateValve::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::initializeCouplingVariables() \n";
#endif

    super_Type::initializeCouplingVariables();

    if ( localCouplingVariables( 0 )[0] <= 1e-10 )
    {
        if ( M_comm->MyPID() == 0 )
            std::cout << " MS-  Valve closed at coupling " << M_ID << std::endl;

        M_valveIsOpen = false;
        localCouplingVariables( 0 ) = 0;
    }
    else
    {
        if ( M_comm->MyPID() == 0 )
            std::cout << " MS-  Valve open at coupling " << M_ID << std::endl;

        M_valveIsOpen = true;
    }
}

void
MultiscaleCouplingFlowRateValve::updateCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::updateCoupling() \n";
#endif

    super_Type::updateCoupling();

    M_topologyChange = false;
    if ( M_valveIsOpen )
    {
        if ( multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryFlowRate( M_flags[0] ) < 0 )
        {
            if ( M_comm->MyPID() == 0 )
                std::cout << " MS-  Opening the valve at coupling " << M_ID << std::endl;

            M_valveIsOpen = false;
            M_topologyChange = true;

            // Reset coupling variable history
            resetCouplingHistory();
        }
    }
    else
    {
        if (   multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[1] )->boundaryStress( M_flags[1] )
             - multiscaleDynamicCast< MultiscaleInterfaceFluid >( M_models[0] )->boundaryStress( M_flags[0] ) > 0 )
        {
            if ( M_comm->MyPID() == 0 )
                std::cout << " MS-  Closing the valve at coupling " << M_ID << std::endl;

            M_valveIsOpen = true;
            M_topologyChange = true;

            // Reset coupling variable history
            resetCouplingHistory();
        }
    }
}

void
MultiscaleCouplingFlowRateValve::exportCouplingResiduals( multiscaleVector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::exportCouplingResiduals()  \n";
#endif

    if ( M_valveIsOpen )
        super_Type::exportCouplingResiduals( couplingResiduals );
    else
    {
        *M_localCouplingResiduals = 0.;
        exportCouplingVector( *M_localCouplingResiduals, couplingResiduals );
    }

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8250 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

void
MultiscaleCouplingFlowRateValve::showMe()
{
    if ( M_displayer->isLeader() )
    {
        super_Type::showMe();

        std::cout << "Valve position      = " << M_valveIsOpen << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private MultiscaleCoupling Implementation
// ===================================================
void
MultiscaleCouplingFlowRateValve::insertJacobianConstantCoefficients( multiscaleMatrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::insertJacobianConstantCoefficients( jacobian )  \n";
#endif

    if ( M_valveIsOpen )
        super_Type::insertJacobianConstantCoefficients( jacobian );
    else
    {
        UInt row    = M_couplingIndex.second;
        UInt column = M_couplingIndex.second;

        if ( M_comm->MyPID() == 0 )
            for ( UInt i( 0 ); i < modelsNumber(); ++i )
                jacobian.addToCoefficient( row + i, column + i, 1 );
    }
}

void
MultiscaleCouplingFlowRateValve::insertJacobianDeltaCoefficients( multiscaleMatrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8250 ) << "MultiscaleCouplingFlowRateValve::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem )  \n";
#endif

    if ( M_valveIsOpen )
        super_Type::insertJacobianDeltaCoefficients( jacobian, column, ID, solveLinearSystem );
}

} // Namespace Multiscale
} // Namespace LifeV
