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
 *  @brief File containing the Multiscale Coupling BoundaryCondition
 *
 *  @date 02-09-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/multiscale/solver/MultiscaleCouplingBoundaryCondition.hpp>

namespace LifeV
{
namespace Multiscale
{

// ===================================================
// Constructors & Destructor
// ===================================================
MultiscaleCouplingBoundaryCondition::MultiscaleCouplingBoundaryCondition() :
    multiscaleCoupling_Type       (),
    M_fileName                    (),
    M_list                        (),
    M_listSize                    ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8210 ) << "MultiscaleCouplingBoundaryCondition::MultiscaleCouplingBoundaryCondition() \n";
#endif

    M_type = BoundaryCondition;
}

// ===================================================
// Multiscale PhysicalCoupling Implementation
// ===================================================
void
MultiscaleCouplingBoundaryCondition::setupData ( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8210 ) << "MultiscaleCouplingBoundaryCondition::setupData() \n";
#endif

    multiscaleCoupling_Type::setupData ( fileName );

    if ( modelsNumber() > 0 )
    {
        M_fileName = fileName;
        GetPot dataFile ( fileName );

        //Load the list of boundary conditions
        M_listSize = dataFile.vector_variable_size ( "boundary_conditions/list" );

        M_list.reserve ( M_listSize );
        for ( UInt i ( 0 ); i < M_listSize; ++i )
        {
            M_list.push_back ( dataFile ( "boundary_conditions/list", " ", i ) );
        }
    }
}

void
MultiscaleCouplingBoundaryCondition::setupCouplingVariablesNumber()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8210 ) << "MultiscaleCouplingBoundaryCondition::setupCouplingVariablesNumber() \n";
#endif

    M_couplingVariablesNumber = 0;
}

void
MultiscaleCouplingBoundaryCondition::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 8210 ) << "MultiscaleCouplingBoundaryCondition::setupCoupling() \n";
#endif

    if ( myModelsNumber() > 0 )
    {
        // Impose boundary conditions on all the models
        for ( UInt i ( 0 ); i < modelsNumber(); ++i )
            if ( myModel ( i ) )
                switch ( M_models[i]->type() )
                {
#if defined(LIFEV_HAS_ZERODIMENSIONAL)
                    case Windkessel0D:

                        applyBoundaryConditions0D< MultiscaleModelWindkessel0D > ( i );

                        break;

                    case ZeroDimensional:

                        applyBoundaryConditions0D< MultiscaleModel0D > ( i );

                        break;
#endif
#if defined(LIFEV_HAS_ONEDFSI)
                    case FSI1D:

                        applyBoundaryConditions1D< MultiscaleModelFSI1D > ( i );

                        break;
#endif
#if defined(LIFEV_HAS_NAVIERSTOKES)
                    case Fluid3D:

                        applyBoundaryConditions3D< MultiscaleModelFluid3D > ( i );

                        break;
#endif
#if defined(LIFEV_HAS_FSI)
                    case FSI3D:

                        applyBoundaryConditions3D< MultiscaleModelFSI3D > ( i );

                        break;
#endif
                    default:

                        switchErrorMessage ( M_models[i] );

                        break;
                }
    }
}

} // Namespace Multiscale
} // Namespace LifeV
