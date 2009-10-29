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

std::map< std::string, MS_Coupling_FluxStress::stressTypes > MS_Coupling_FluxStress::M_mapStress;

// ===================================================
//! Constructors
// ===================================================
MS_Coupling_FluxStress::MS_Coupling_FluxStress() :
    super              (),
    M_baseFlux         (),
    M_baseStress       (),
    M_baseDeltaFlux    (),
    M_baseDeltaStress  (),
    M_stressType       ()
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress() \n";
#endif

    M_type = FluxStress;
}

MS_Coupling_FluxStress::MS_Coupling_FluxStress( const MS_Coupling_FluxStress& FluxStress ) :
    super              ( FluxStress ),
    M_baseFlux         ( FluxStress.M_baseFlux ),
    M_baseStress       ( FluxStress.M_baseStress ),
    M_baseDeltaFlux    ( FluxStress.M_baseDeltaFlux ),
    M_baseDeltaStress  ( FluxStress.M_baseDeltaStress ),
    M_stressType       ( FluxStress.M_stressType )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::MS_Coupling_FluxStress( FluxStress ) \n";
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
        M_baseFlux         = FluxStress.M_baseFlux;
        M_baseStress       = FluxStress.M_baseStress;
        M_baseDeltaFlux    = FluxStress.M_baseDeltaFlux;
        M_baseDeltaStress  = FluxStress.M_baseDeltaStress;
        M_stressType       = FluxStress.M_stressType;
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
    Debug( 8230 ) << "MS_Coupling_FluxStress::SetupData() \n";
#endif

    //Set type of pressure coupling: Static, Total
    M_mapStress["StaticPressure"] = StaticPressure;
    M_mapStress["TotalPressure"]  = TotalPressure;
    M_stressType = M_mapStress[M_dataFile( "MultiScale/pressureType", "StaticPressure" )];

    //Set number of coupling variables
    M_couplingIndex.first  = 2;

    int MyGlobalElements[M_couplingIndex.first];
    for ( int i(0) ; i < static_cast< int > ( M_couplingIndex.first ) ; ++i )
        MyGlobalElements[i] = i;
    EpetraMap map( M_couplingIndex.first, M_couplingIndex.first, &MyGlobalElements[0], 0, *M_comm );

    M_LocalCouplingVariables.reset( new EpetraVector( map, Unique ) );
    M_LocalCouplingResiduals.reset( new EpetraVector( map, Unique ) );
    M_LocalDeltaCouplingVariables.reset( new EpetraVector( map, Unique ) );
    M_LocalJacobianProduct.reset( new EpetraVector( map, Unique ) );

    //Set Functions
    M_baseFlux.setFunction  ( boost::bind( &MS_Coupling_FluxStress::functionFlux,   this, _1, _2, _3, _4, _5 ) );
    M_baseStress.setFunction( boost::bind( &MS_Coupling_FluxStress::functionStress, this, _1, _2, _3, _4, _5 ) );

    M_baseDeltaFlux.setFunction  ( boost::bind( &MS_Coupling_FluxStress::functionDeltaFlux,   this, _1, _2, _3, _4, _5 ) );
    M_baseDeltaStress.setFunction( boost::bind( &MS_Coupling_FluxStress::functionDeltaStress, this, _1, _2, _3, _4, _5 ) );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxStress::SetupCoupling( void )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::SetupCoupling() \n";
#endif

    // Impose flux
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:

            imposeFlux< MS_Model_Fluid3D > ();
            imposeDeltaFlux< MS_Model_Fluid3D > ();

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
                imposeDeltaStress< MS_Model_Fluid3D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_FluxStress::InitializeCouplingVariables( void )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::InitializeCouplingVariables() \n";
#endif

    *M_LocalCouplingVariables = CouplingFunction();
    *M_LocalDeltaCouplingVariables = 0.;

#ifdef DEBUG
    PrintQuantities();
#endif
}

void
MS_Coupling_FluxStress::ExportCouplingResiduals( VectorType& CouplingResiduals )
{
#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::ExportCouplingResiduals() \n";
#endif

    *M_LocalCouplingResiduals  = CouplingFunction();
    *M_LocalCouplingResiduals -= *M_LocalCouplingVariables;

    ExportCouplingVector( *M_LocalCouplingResiduals, CouplingResiduals );

#ifdef DEBUG
    PrintQuantities();
#endif
}

void
MS_Coupling_FluxStress::ComputeJacobianProduct( const VectorType& deltaCouplingVariables )
{

#ifdef DEBUG
    Debug( 8230 ) << "MS_Coupling_FluxStress::ComputeJacobianProduct() \n";
#endif

    // Impose delta stress and compute delta flux
    ( *M_LocalDeltaCouplingVariables )[1] = deltaCouplingVariables[1];
    ( *M_LocalJacobianProduct )[0] = 0.;
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                ( *M_LocalJacobianProduct )[0] -= computeDeltaFlux< MS_Model_Fluid3D > ( i );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    // Restore delta stress = 0
    ( *M_LocalDeltaCouplingVariables )[1] = 0;

    // Impose delta flux and compute delta stress
    ( *M_LocalDeltaCouplingVariables )[0] = deltaCouplingVariables[0];
    switch ( M_models[0]->GetType() )
    {
        case Fluid3D:
        {
            ( *M_LocalJacobianProduct )[1] = computeDeltaStress< MS_Model_Fluid3D > (); //Coupling through the stress

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }
    // Restore delta flux = 0
    ( *M_LocalDeltaCouplingVariables )[0] = 0;

    *M_LocalJacobianProduct -= deltaCouplingVariables;
}

void
MS_Coupling_FluxStress::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Pressure Type       = " << Enum2String( M_stressType, M_mapStress ) << std::endl;
        std::cout << "Coupling Flux       = " << ( *M_LocalCouplingVariables )[0] << std::endl
                  << "Coupling Stress     = " << ( *M_LocalCouplingVariables )[1] << std::endl << std::endl;
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
    Debug( 8230 ) << "MS_Coupling_FluxStress::CouplingFunction() \n";
#endif

    EpetraVector couplingFunction ( *M_LocalCouplingVariables ); couplingFunction = 0;

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
            couplingFunction[1] = computeStress< MS_Model_Fluid3D > (); //Coupling through the stress

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[0] );
    }

    return couplingFunction;
}

Real
MS_Coupling_FluxStress::functionFlux( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalCouplingVariables )[1];
}

Real
MS_Coupling_FluxStress::functionDeltaFlux( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[0];
}

Real
MS_Coupling_FluxStress::functionDeltaStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[1];
}

void
MS_Coupling_FluxStress::PrintQuantities( void )
{
    std::stringstream output;
    output << std::scientific << std::setprecision( 6 );
    output << "MS_Coupling_FluxStress::PrintQuantities() \n";

    output << "                                                       Coupling ID       = " << M_ID << "\n";
    output << "                                                       Coupling Flux     = " << ( *M_LocalCouplingVariables )[0] << "\n";
    output << "                                                       Coupling Stress   = " << ( *M_LocalCouplingVariables )[1] << "\n";

    output << "                                                       Model | Flux            Static Pressure     Dynamic Pressure\n";
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        output << "                                                       ";
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                output << M_models[i]->GetID() << "       " << computeFlux< MS_Model_Fluid3D > ( i ) << "    " << computeStaticPressure< MS_Model_Fluid3D > ( i ) << "        " << computeDynamicPressure< MS_Model_Fluid3D > ( i ) << "\n";

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    }
    Debug( 8230 ) << output.str() << "\n";
}

} // Namespace LifeV
