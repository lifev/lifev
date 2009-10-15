/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-08-24

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Coupling_FluxStress.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-08-24
 */

#include <lifemc/lifesolver/MS_Coupling_FluxStress.hpp>

namespace LifeV {

std::map< std::string, MS_Coupling_FluxStress::pressureTypes > MS_Coupling_FluxStress::M_mapPressure;

// ===================================================
//! Constructors
// ===================================================
MS_Coupling_FluxStress::MS_Coupling_FluxStress() :
    super           (),
    M_baseFlux      (),
    M_baseStress    (),
    M_pressureType  ()
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress() \n";
#endif

    M_type = FluxStress;
}

MS_Coupling_FluxStress::MS_Coupling_FluxStress( const MS_Coupling_FluxStress& FluxStress ) :
    super           ( FluxStress ),
    M_baseFlux      ( FluxStress.M_baseFlux ),
    M_baseStress    ( FluxStress.M_baseStress ),
    M_pressureType  ( FluxStress.M_pressureType )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress( FluxStress ) \n";
#endif

}

// ===================================================
//! Methods
// ===================================================
MS_Coupling_FluxStress&
MS_Coupling_FluxStress::operator=( const MS_Coupling_FluxStress& FluxStress )
{
    if ( this != &FluxStress )
    {
        super::operator=( FluxStress );
        M_baseFlux     = FluxStress.M_baseFlux;
        M_baseStress   = FluxStress.M_baseStress;
        M_pressureType = FluxStress.M_pressureType;
    }
    return *this;
}

// ===================================================
//! MultiScale Physical Coupling
// ===================================================
void
MS_Coupling_FluxStress::SetupData( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::SetupData() \n";
#endif

    //Set type of pressure coupling: Static, Total
    M_mapPressure["Static"] = Static;
    M_mapPressure["Total"]  = Total;
    M_pressureType = M_mapPressure[M_dataFile( "MultiScale/pressureType", "Static" )];

    //Create coupling Vectors & Functions
    int MyGlobalIElements[2]; MyGlobalIElements[0] = 0; MyGlobalIElements[1] = 1;
    EpetraMap map( 2, 2, &MyGlobalIElements[0], 0, *M_comm );

    M_couplingVariables.reset( new EpetraVector( map, Unique ) );
    M_couplingResiduals.reset( new EpetraVector( map, Unique ) );

    //Set Functions
    M_baseFlux.setFunction  ( boost::bind( &MS_Coupling_FluxStress::functionFlux,   this, _1, _2, _3, _4, _5 ) );
    M_baseStress.setFunction( boost::bind( &MS_Coupling_FluxStress::functionStress, this, _1, _2, _3, _4, _5 ) );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxStress::SetupCoupling( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::SetupCoupling() \n";
#endif

    // Impose flux
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:

            imposeFlux< MS_Model_Fluid3D > ();

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }

    // Impose pressure
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                imposeStress< MS_Model_Fluid3D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxStress::SetupImplicitCoupling( ContainerOfVectors< EpetraVector >& couplingVariables,
                                               ContainerOfVectors< EpetraVector >& couplingResiduals )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::SetupImplicitCoupling( couplingVariables, couplingResiduals ) \n";
#endif

    couplingVariables.push_back( M_couplingVariables );
    couplingResiduals.push_back( M_couplingResiduals );
}

void
MS_Coupling_FluxStress::UpdateCouplingVariables( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::UpdateCouplingVariables() ";
#endif

    *M_couplingVariables = CouplingFunction();
}

void
MS_Coupling_FluxStress::UpdateCouplingResiduals( void )
{
#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_FluxStress::UpdateCouplingResidual()  ";
#endif

    *M_couplingResiduals  = CouplingFunction();
    *M_couplingResiduals -= *M_couplingVariables;
}

void
MS_Coupling_FluxStress::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Pressure Type       = " << Enum2String( M_pressureType, M_mapPressure ) << std::endl;
        std::cout << "Coupling Flux       = " << ( *M_couplingVariables )[0] << std::endl
                  << "Coupling Stress     = " << ( *M_couplingVariables )[1] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}

// ===================================================
//! Private Methods
// ===================================================
EpetraVector
MS_Coupling_FluxStress::CouplingFunction( void )
{

#ifdef DEBUG
    std::stringstream output;
    output << std::scientific << std::setprecision( 6 );
    output << "MS_Coupling_FluxStress::CouplingFunction() \n";
#endif

    EpetraVector couplingFunction = *M_couplingVariables;
    couplingFunction = 0.;

    // Compute the Flux
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                couplingFunction[0] -= computeFlux< MS_Model_Fluid3D > ( i );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    // Compute the Stress
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:
        {
            couplingFunction[1] = -selectPressureType< MS_Model_Fluid3D > (); //Coupling through the stress

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }

#ifdef DEBUG
    output << "                                                       Coupling ID       = " << M_ID << "\n";
    output << "                                                       Coupling Flux     = " << couplingFunction[0] << "\n";
    output << "                                                       Coupling Stress   = " << couplingFunction[1] << "\n\n";

    output << "                                                       Model | Flux            Static Pressure     Total Pressure\n";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        output << "                                                       ";
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                output << M_models[i]->GetID() << "       " << computeFlux< MS_Model_Fluid3D > ( i ) << "    " << computePressure< MS_Model_Fluid3D > ( i ) << "        " << computeTotalPressure< MS_Model_Fluid3D > ( i ) << "\n";

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    }
    Debug( 8220 ) << output.str() << "\n";
#endif

    return couplingFunction;
}

Real
MS_Coupling_FluxStress::functionFlux( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_couplingVariables )[0];
}

Real
MS_Coupling_FluxStress::functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_couplingVariables )[1];
}

} // Namespace LifeV
