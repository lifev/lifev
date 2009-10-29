/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-20

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
 \file MS_Coupling_Stress.cpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-20
 */

#include <lifemc/lifesolver/MS_Coupling_Stress.hpp>

namespace LifeV {

std::map< std::string, MS_Coupling_Stress::stressTypes > MS_Coupling_Stress::M_mapStress;

// ===================================================
//! Constructors
// ===================================================
MS_Coupling_Stress::MS_Coupling_Stress() :
    super              (),
    M_baseStress       (),
    M_baseDeltaStress  (),
    M_stressType       ()
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::MS_Coupling_Stress() \n";
#endif

    M_type = Stress;
}

MS_Coupling_Stress::MS_Coupling_Stress( const MS_Coupling_Stress& Stress ) :
    super             ( Stress ),
    M_baseStress      ( Stress.M_baseStress ),
    M_baseDeltaStress ( Stress.M_baseDeltaStress ),
    M_stressType      ( Stress.M_stressType )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::MS_Coupling_Stress( Stress ) \n";
#endif

}

// ===================================================
//! Methods
// ===================================================
MS_Coupling_Stress&
MS_Coupling_Stress::operator=( const MS_Coupling_Stress& Stress )
{
    if ( this != &Stress )
    {
        super::operator=( Stress );
        M_baseStress      = Stress.M_baseStress;
        M_baseDeltaStress = Stress.M_baseDeltaStress;
        M_stressType      = Stress.M_stressType;
    }
    return *this;
}

// ===================================================
//! MultiScale Physical Coupling
// ===================================================
void
MS_Coupling_Stress::SetupData( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::SetupData() \n";
#endif

    //Set type of pressure coupling: StaticPressure, TotalPressure
    M_mapStress["StaticPressure"] = StaticPressure;
    M_mapStress["TotalPressure"]  = TotalPressure;
    M_stressType = M_mapStress[M_dataFile( "MultiScale/pressureType", "StaticPressure" )];

    //Set number of coupling variables
    M_couplingIndex.first = modelsNumber();

    int MyGlobalElements[M_couplingIndex.first];
    for ( int i(0) ; i < static_cast< int > ( M_couplingIndex.first ) ; ++i )
        MyGlobalElements[i] = i;
    EpetraMap map( M_couplingIndex.first, M_couplingIndex.first, &MyGlobalElements[0], 0, *M_comm );

    M_LocalCouplingVariables.reset( new EpetraVector( map, Unique ) );
    M_LocalCouplingResiduals.reset( new EpetraVector( map, Unique ) );
    M_LocalDeltaCouplingVariables.reset( new EpetraVector( map, Unique ) );
    M_LocalJacobianProduct.reset( new EpetraVector( map, Unique ) );

    //Set Functions
    M_baseStress.setFunction( boost::bind( &MS_Coupling_Stress::functionStress, this, _1, _2, _3, _4, _5 ) );

    M_baseDeltaStress.setFunction( boost::bind( &MS_Coupling_Stress::functionDeltaStress, this, _1, _2, _3, _4, _5 ) );

    //MPI Barrier
    M_comm->Barrier();
}

void
MS_Coupling_Stress::SetupCoupling( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::SetupCoupling() \n";
#endif

    // Impose pressure
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
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
MS_Coupling_Stress::InitializeCouplingVariables( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InitializeCouplingVariables() \n";
#endif

    *M_LocalCouplingVariables = CouplingFunction();
    *M_LocalDeltaCouplingVariables = 0.;

#ifdef DEBUG
    PrintQuantities();
#endif
}

void
MS_Coupling_Stress::ExportCouplingResiduals( VectorType& CouplingResiduals )
{
#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::ExportCouplingResiduals()  \n";
#endif

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ( *M_LocalCouplingResiduals )[i] = computeFlux< MS_Model_Fluid3D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
    {
        ( *M_LocalCouplingResiduals )[0]  += ( *M_LocalCouplingVariables )[i];
        ( *M_LocalCouplingResiduals )[i]  -= ( *M_LocalCouplingVariables )[i];
    }

    ExportCouplingVector( *M_LocalCouplingResiduals, CouplingResiduals );

#ifdef DEBUG
    PrintQuantities();
#endif
}

void
MS_Coupling_Stress::ComputeJacobianProduct( const VectorType& deltaCouplingVariables )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::ComputeJacobianProduct()  \n";
#endif


    // Impose delta stress and compute delta flux
    ( *M_LocalDeltaCouplingVariables )[0] = deltaCouplingVariables[0];
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:

                ( *M_LocalJacobianProduct )[i] = computeDeltaFlux< MS_Model_Fluid3D > ( i );

                break;

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    // Restore delta stress = 0
    ( *M_LocalDeltaCouplingVariables )[0] = 0;

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
    {
        ( *M_LocalJacobianProduct )[0]  += deltaCouplingVariables[i];
        ( *M_LocalJacobianProduct )[i]  -= deltaCouplingVariables[i];
    }
}

void
MS_Coupling_Stress::ShowMe( void )
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();

        std::cout << "Pressure Type       = " << Enum2String( M_stressType, M_mapStress ) << std::endl;
        std::cout << "Coupling Stress     = " << ( *M_LocalCouplingVariables )[0] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }

    //MPI Barrier
    M_comm->Barrier();
}

// ===================================================
//! Private Methods
// ===================================================
EpetraVector
MS_Coupling_Stress::CouplingFunction( void )
{

#ifdef DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::CouplingFunction() \n";
#endif

    EpetraVector couplingFunction( *M_LocalCouplingVariables ); couplingFunction = 0.;

    // Compute the Stress coupling variable as an average of all the stresses
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                couplingFunction[0] += computeStress< MS_Model_Fluid3D > ( i );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }
    couplingFunction[0] /=modelsNumber();

    // Compute Flux coupling variables
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->GetType() )
        {
            case Fluid3D:
            {
                couplingFunction[i] = computeFlux< MS_Model_Fluid3D > ( i );

                break;
            }

            default:

                if ( M_displayer->isLeader() )
                    switchErrorMessage( M_models[i] );
        }

    return couplingFunction;
}

Real
MS_Coupling_Stress::functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ )
{
    return ( *M_LocalCouplingVariables )[0];
}

Real
MS_Coupling_Stress::functionDeltaStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/)
{
    return ( *M_LocalDeltaCouplingVariables )[0];
}

void
MS_Coupling_Stress::PrintQuantities( void )
{
    std::stringstream output;
    output << std::scientific << std::setprecision( 6 );
    output << "MS_Coupling_Stress::PrintQuantities() \n";

    output << "                                                       Coupling ID       = " << M_ID << "\n";
    output << "                                                       Coupling Stress   = " << ( *M_LocalCouplingVariables )[0] << "\n";

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
    Debug( 8220 ) << output.str() << "\n";
}

} // Namespace LifeV
