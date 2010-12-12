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
 *  @brief File containing the MultiScale Coupling Stress
 *
 *  @date 20-10-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/MS_Coupling_Stress.hpp>

namespace LifeV
{

std::map< std::string, stress_Type > MS_stressesMap;

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Coupling_Stress::MS_Coupling_Stress() :
        MS_Coupling_Type     (),
        M_baseStress3D       (),
        M_baseStress1D       (),
        M_stressType         ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::MS_Coupling_Stress() \n";
#endif

    M_type = Stress;
}

// ===================================================
// MultiScale PhysicalCoupling Implementation
// ===================================================
void
MS_Coupling_Stress::setupData( const std::string& fileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::SetupData() \n";
#endif

    MS_Coupling_Type::setupData( fileName );

    GetPot DataFile( fileName );

    //Set type of stress coupling
    M_stressType = MS_stressesMap[DataFile( "MultiScale/stressType", "StaticPressure" )];
}

void
MS_Coupling_Stress::setupCoupling()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::setupCoupling() \n";
#endif

    //Set number of coupling variables
    M_couplingIndex.first = modelsNumber();

    //Create local vectors
    createLocalVectors();

    //Set Functions
    M_baseStress3D.setFunction( boost::bind( &MS_Coupling_Stress::functionStress, this, _1, _2, _3, _4, _5 ) );

    M_baseStress1D.setFunction( boost::bind( &MS_Coupling_Stress::functionStress, this, _1, _1, _1, _1, _1 ) );

    // Impose stress
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:

            imposeStress3D< MS_Model_Fluid3D > ( i );

            break;

        case FSI3D:

            imposeStress3D< MS_Model_FSI3D > ( i );

            break;

        case OneDimensional:

            imposeStress1D< MS_Model_1D > ( i );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }
}

void
MS_Coupling_Stress::initializeCouplingVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InitializeCouplingVariables() \n";
#endif

    *M_localCouplingVariables[0] = 0.;

    // Compute the Stress coupling variable as an average of all the stresses
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            ( *M_localCouplingVariables[0] )[0] += MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryStress( M_flags[i], M_stressType );

            break;
        }

        case FSI3D:
        {
            ( *M_localCouplingVariables[0] )[0] += MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryStress( M_flags[i], M_stressType );

            break;
        }

        case OneDimensional:
        {
            ( *M_localCouplingVariables[0] )[0] += MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryStress( M_flags[i], M_stressType );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }
    ( *M_localCouplingVariables[0] )[0] /= modelsNumber();

    // Compute Flux coupling variables
    for ( UInt i( 1 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {

            ( *M_localCouplingVariables[0] )[i] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case FSI3D:
        {

            ( *M_localCouplingVariables[0] )[i] = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        case OneDimensional:
        {

            ( *M_localCouplingVariables[0] )[i] = MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8220 ) << "C(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingVariables[0] )[i]  << "\n";
#endif

}

void
MS_Coupling_Stress::exportCouplingResiduals( MS_Vector_Type& couplingResiduals )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::ExportCouplingResiduals()  \n";
#endif

    for ( UInt i( 0 ); i < modelsNumber(); ++i )
        switch ( M_models[i]->type() )
        {
        case Fluid3D:

            ( *M_localCouplingResiduals )[i] = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;

        case FSI3D:

            ( *M_localCouplingResiduals )[i] = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;

        case OneDimensional:

            ( *M_localCouplingResiduals )[i] = MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );

            break;

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

    for ( UInt i( 1 ); i < modelsNumber(); ++i )
    {
        ( *M_localCouplingResiduals )[0]  += ( *M_localCouplingVariables[0] )[i];
        ( *M_localCouplingResiduals )[i]  -= ( *M_localCouplingVariables[0] )[i];
    }

    exportCouplingVector( *M_localCouplingResiduals, couplingResiduals );

#ifdef HAVE_LIFEV_DEBUG
    for ( UInt i( 0 ); i < M_couplingIndex.first; ++i )
        Debug( 8220 ) << "R(" << M_couplingIndex.second + i << ") = " << ( *M_localCouplingResiduals )[i]  << "\n";
#endif

}

void
MS_Coupling_Stress::showMe()
{
    if ( M_displayer->isLeader() )
    {
        MS_Coupling_Type::showMe();

        std::cout << "Stress Type         = " << Enum2String( M_stressType, MS_stressesMap ) << std::endl;
        std::cout << "Coupling Stress     = " << ( *M_localCouplingVariables[0] )[0] << std::endl << std::endl;
        std::cout << std::endl << std::endl;
    }
}

// ===================================================
// Private MultiScale PhysicalCoupling Implementation
// ===================================================
MS_ModelsVector_Type
MS_Coupling_Stress::listOfPerturbedModels( const UInt& localCouplingVariableID )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::GetListOfPerturbedModels( localCouplingVariableID ) \n";
#endif

    MS_ModelsVector_Type perturbedModelsList(0);

    if ( localCouplingVariableID == 0 )
    {
        perturbedModelsList.resize( modelsNumber() );

        for ( UInt i( 0 ); i < modelsNumber(); ++i )
            perturbedModelsList[i] = M_models[i];
    }

    return perturbedModelsList;
}

void
MS_Coupling_Stress::insertJacobianConstantCoefficients( MS_Matrix_Type& jacobian )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InsertJacobianConstantCoefficients( jacobian )  \n";
#endif

    UInt row    = M_couplingIndex.second;
    UInt column = M_couplingIndex.second;

    if ( M_comm->MyPID() == 0 )
        for ( UInt i( 1 ); i < modelsNumber(); ++i )
        {
            jacobian.set_mat_inc( row,     column + i,  1 );
            jacobian.set_mat_inc( row + i, column + i, -1 );
        }
}

void
MS_Coupling_Stress::insertJacobianDeltaCoefficients( MS_Matrix_Type& jacobian, const UInt& column, const UInt& ID, bool& solveLinearSystem )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "MS_Coupling_Stress::InsertJacobianDeltaCoefficients( jacobian, column, ID, LinearSystemSolved )  \n";
#endif

    // Definitions
    Real coefficient  = 0;
    UInt row          = 0;
    UInt modelLocalID = modelGlobalToLocalID( ID );

    switch ( M_models[modelLocalID]->type() )
    {
    case Fluid3D:

        coefficient = MS_DynamicCast< MS_Model_Fluid3D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;

    case FSI3D:

        coefficient = MS_DynamicCast< MS_Model_FSI3D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;

    case OneDimensional:

        coefficient = MS_DynamicCast< MS_Model_1D >( M_models[modelLocalID] )->boundaryDeltaFlowRate( M_flags[modelLocalID], solveLinearSystem );

        break;

    default:

        if ( M_displayer->isLeader() )
            switchErrorMessage( M_models[modelLocalID] );
    }

    // Add coefficient to the matrix
    row = M_couplingIndex.second + modelLocalID;
    if ( M_comm->MyPID() == 0 )
        jacobian.set_mat_inc( row, column, coefficient );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8220 ) << "J(" << row << "," << column << ") = " << coefficient  << "\n";
#endif

}

void
MS_Coupling_Stress::displayCouplingValues( std::ostream& output )
{
    Real flowRate(0), stress(0), pressure(0), dynamicPressure(0);
    for ( UInt i( 0 ); i < modelsNumber(); ++i )
    {
        switch ( M_models[i]->type() )
        {
        case Fluid3D:
        {
            flowRate        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[0];
            pressure        = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = MS_DynamicCast< MS_Model_Fluid3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case FSI3D:
        {
            flowRate        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[0];
            pressure        = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = MS_DynamicCast< MS_Model_FSI3D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        case OneDimensional:
        {
            flowRate        = MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryFlowRate( M_flags[i] );
            stress          = ( *M_localCouplingVariables[0] )[0];
            pressure        = MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryPressure( M_flags[i] );
            dynamicPressure = MS_DynamicCast< MS_Model_1D >( M_models[i] )->boundaryDynamicPressure( M_flags[i] );

            break;
        }

        default:

            if ( M_displayer->isLeader() )
                switchErrorMessage( M_models[i] );
        }

        if ( M_comm->MyPID() == 0 )
            output << "  " << M_globalData->dataTime()->getTime() << "    " << M_models[i]->ID()
            << "    " << M_flags[i]
            << "    " << flowRate
            << "    " << stress
            << "    " << pressure
            << "    " << dynamicPressure << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================
Real
MS_Coupling_Stress::functionStress( const Real& t, const Real&, const Real&, const Real&, const UInt& )
{
    MS_Vector_Type interpolatedCouplingVariables( *M_localCouplingVariables[0] );

    timeContainer_Type timeContainer( M_timeInterpolationOrder + 1 );
    for ( UInt i(0) ; i <= M_timeInterpolationOrder ; ++i )
        timeContainer[i] = M_globalData->dataTime()->getTime() - i * M_globalData->dataTime()->getTimeStep();

    interpolateCouplingVariables( timeContainer, t, interpolatedCouplingVariables );

    return interpolatedCouplingVariables[0];
}

} // Namespace LifeV
