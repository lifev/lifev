//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for 1D model data handling.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @date 01-07-2004
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-04-2010
 */

#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
OneDimensionalModel_Data::OneDimensionalModel_Data():
    M_PhysicsType               (),
    M_FluxType                  (),
    M_SourceType                (),
    M_Time                      (),
    M_Mesh                      ( new Mesh_Type() ),
    M_postprocessingDirectory   (),
    M_postprocessingFile        (),
    M_verbose                   (),
    M_UW                        (),
    M_inertial_wall             (),
    M_viscoelastic_wall         (),
    M_linearize_string_model    (),
    M_linearize_equations       (),
    M_longitudinal_wall         (),
    M_flux_second_der           (),
    M_dP_dt_steps               (),
    M_CFLmax                    (),
    M_initialVariable           (),
    M_initialValue              (),
    M_restValue                 (),
    M_multiplier                (),
    M_ComputeCoefficients       (),
    M_PowerlawCoefficient       (),
    M_Density                   (),
    M_Viscosity                 (),
    M_DensityWall               (),
    M_ThickVessel               (),
    M_Thickness                 (),
    M_Young                     (),
    M_Poisson                   (),
    M_ViscoelasticModulus       (),
    M_InertialModulus           (),
    M_RobertsonCorrection       (),
    M_Area0                     (),
    M_dArea0dz                  (),
    M_PressBeta0                (),
    M_dPressBeta0dz             (),
    M_PressBeta1                (),
    M_dPressBeta1dz             (),
    M_AlphaCoriolis             (),
    M_dAlphaCoriolisdz          (),
    M_FrictionKr                (),
    M_Flux11                    (),
    M_Flux12                    (),
    M_Flux21                    (),
    M_Flux22                    (),
    M_Celerity1                 (),
    M_Celerity2                 (),
    M_Celerity1LeftEigenvector1 (),
    M_Celerity1LeftEigenvector2 (),
    M_Celerity2LeftEigenvector1 (),
    M_Celerity2LeftEigenvector2 (),
    M_Source10                  (),
    M_Source20                  (),
    M_Source11                  (),
    M_Source12                  (),
    M_Source21                  (),
    M_Source22                  ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Data::setup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_PhysicsType = OneDimensionalModel_PhysicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_FluxType    = OneDimensionalModel_FluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_SourceType  = OneDimensionalModel_SourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_Time.get() )
        M_Time.reset( new Time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    M_Mesh->setup( dataFile( ( section + "/space_discretization/Length"           ).data(), 1. ),
                   dataFile( ( section + "/space_discretization/NumberOfElements" ).data(), 10 ) );

    //std::cout << " 1D- Mesh nodes:                               " << M_Mesh->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_Mesh->numElements() << std::endl;

    // Miscellaneous
    M_postprocessingDirectory= dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile     = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_verbose                = dataFile( ( section + "/miscellaneous/verbose"                        ).data(), 1 );
    M_UW                     = dataFile( ( section + "/miscellaneous/alternate_solver"               ).data(), false );
    M_inertial_wall          = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelastic_wall      = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_linearize_string_model = dataFile( ( section + "/miscellaneous/linearize_string_model"         ).data(), true );
    M_linearize_equations    = dataFile( ( section + "/miscellaneous/linearize_equations"            ).data(), false );
    M_longitudinal_wall      = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );
    M_flux_second_der        = dataFile( ( section + "/miscellaneous/compute_flux_second_derivative" ).data(), false );
    M_dP_dt_steps            = dataFile( ( section + "/miscellaneous/pressure_derivative_steps"      ).data(), 1 );
    M_CFLmax                 = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3)/3. );

    if ( M_CFLmax > std::sqrt(3)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Initialize
    std::map< std::string, OneD_Initialize > initializeMap;
    initializeMap["A"]       = OneD_InitializeArea;
    initializeMap["Q"]       = OneD_InitializeFlux;
    initializeMap["W1"]      = OneD_InitializeRiemann1;
    initializeMap["W2"]      = OneD_InitializeRiemann2;
    initializeMap["P"]       = OneD_InitializePressure;

    M_initialVariable        = initializeMap[ dataFile( ( section + "/initialize/variable"           ).data(), "P" ) ];
    M_initialValue           = dataFile( ( section + "/initialize/initialValue"                      ).data(), 0. );
    M_restValue              = dataFile( ( section + "/initialize/restValue"                         ).data(), 0. );
    M_multiplier             = dataFile( ( section + "/initialize/multiplier"                        ).data(), 1. );

    // Physical Parameters
    M_ComputeCoefficients    = dataFile( ( section + "/PhysicalParameters/ComputeCoefficients"       ).data(), false );
    M_PowerlawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_Density                = dataFile( ( section + "/PhysicalParameters/density"                   ).data(), 1. );
    M_Viscosity              = dataFile( ( section + "/PhysicalParameters/viscosity"                 ).data(), 0.035 );

    M_DensityWall            = dataFile( ( section + "/PhysicalParameters/densityWall"               ).data(), 1. );
    M_ThickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_Thickness              = dataFile( ( section + "/PhysicalParameters/thickness"                 ).data(), 0. );
    M_Young                  = dataFile( ( section + "/PhysicalParameters/young"                     ).data(), 4.0E6 );
    M_Poisson                = dataFile( ( section + "/PhysicalParameters/poisson"                   ).data(), 0.5 );

    M_ViscoelasticModulus    = dataFile( ( section + "/PhysicalParameters/gamma"                     ).data(), 0. );
    M_InertialModulus        = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_RobertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_ComputeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );

    M_Area0.resize( M_Mesh->numPoints() );
    M_dArea0dz.resize( M_Mesh->numPoints() );
    M_PressBeta0.resize( M_Mesh->numPoints() );
    M_dPressBeta0dz.resize( M_Mesh->numPoints() );
    M_PressBeta1.resize( M_Mesh->numPoints() );
    M_dPressBeta1dz.resize( M_Mesh->numPoints() );
    M_FrictionKr.resize( M_Mesh->numPoints() );

    M_AlphaCoriolis.resize( M_Mesh->numPoints() );
    M_dAlphaCoriolisdz.resize( M_Mesh->numPoints() );

    M_Flux11.resize( M_Mesh->numPoints() );
    M_Flux12.resize( M_Mesh->numPoints() );
    M_Flux21.resize( M_Mesh->numPoints() );
    M_Flux22.resize( M_Mesh->numPoints() );
    M_Celerity1.resize( M_Mesh->numPoints() );
    M_Celerity2.resize( M_Mesh->numPoints() );
    M_Celerity1LeftEigenvector1.resize( M_Mesh->numPoints() );
    M_Celerity1LeftEigenvector2.resize( M_Mesh->numPoints() );
    M_Celerity2LeftEigenvector1.resize( M_Mesh->numPoints() );
    M_Celerity2LeftEigenvector2.resize( M_Mesh->numPoints() );
    M_Source10.resize( M_Mesh->numPoints() );
    M_Source20.resize( M_Mesh->numPoints() );
    M_Source11.resize( M_Mesh->numPoints() );
    M_Source12.resize( M_Mesh->numPoints() );
    M_Source21.resize( M_Mesh->numPoints() );
    M_Source22.resize( M_Mesh->numPoints() );

    // Linear Parameters
    for ( UInt i = 0; i < M_Mesh->numPoints() ; ++i )
    {
        // Physical Parameters
        M_Area0[i]                     = dataFile( ( section + "/PhysicalParameters/Area0"           ).data(), M_PI );
        M_dArea0dz[i]                  = 0;
        M_PressBeta0[i]                = dataFile( ( section + "/PhysicalParameters/Beta0"           ).data(), 1.e6 );
        M_dPressBeta0dz[i]             = 0.;
        M_PressBeta1[i]                = dataFile( ( section + "/PhysicalParameters/Beta1"           ).data(), 0.5 );
        M_dPressBeta1dz[i]             = 0.;
        M_FrictionKr[i]                = dataFile( ( section + "/PhysicalParameters/Kr"              ).data(), 1. );

        M_AlphaCoriolis[i]             = dataFile( ( section + "/PhysicalParameters/AlphaCoriolis"   ).data(), 1. / M_RobertsonCorrection );
        M_dAlphaCoriolisdz[i]          = 0.;

        // Linear Parameters
        M_Flux11[i]                    = dataFile( ( section + "/LinearParameters/Flux11"            ).data(), 1. );
        M_Flux12[i]                    = dataFile( ( section + "/LinearParameters/Flux12"            ).data(), 0. );
        M_Flux21[i]                    = dataFile( ( section + "/LinearParameters/Flux21"            ).data(), 0. );
        M_Flux22[i]                    = dataFile( ( section + "/LinearParameters/Flux22"            ).data(), 1. );
        M_Celerity1[i]                 = dataFile( ( section + "/LinearParameters/Celerity1"         ).data(), 1. );
        M_Celerity2[i]                 = dataFile( ( section + "/LinearParameters/Celerity2"         ).data(), 1. );
        M_Celerity1LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector11" ).data(), 1. );
        M_Celerity1LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector12" ).data(), 0. );
        M_Celerity2LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector21" ).data(), 0. );
        M_Celerity2LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector22" ).data(), 1. );
        M_Source10[i]                  = dataFile( ( section + "/LinearParameters/Source10"          ).data(), 0. );
        M_Source20[i]                  = dataFile( ( section + "/LinearParameters/Source20"          ).data(), 0. );
        M_Source11[i]                  = dataFile( ( section + "/LinearParameters/Source11"          ).data(), 0. );
        M_Source12[i]                  = dataFile( ( section + "/LinearParameters/Source12"          ).data(), 0. );
        M_Source21[i]                  = dataFile( ( section + "/LinearParameters/Source21"          ).data(), 0. );
        M_Source22[i]                  = dataFile( ( section + "/LinearParameters/Source22"          ).data(), 0. );
    }

    UpdateCoefficients();
}

void
OneDimensionalModel_Data::oldStyleSetup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_PhysicsType = OneDimensionalModel_PhysicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_FluxType    = OneDimensionalModel_FluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_SourceType  = OneDimensionalModel_SourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_Time.get() )
        M_Time.reset( new Time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    Real length = dataFile( ( section + "/discretization/x_right" ).data(), 1. ) -
                  dataFile( ( section + "/discretization/x_left"  ).data(), 0. );

    M_Mesh->setup( length,
                   dataFile( ( section + "/discretization/nb_elem" ).data(), 10 ) );

    //std::cout << " 1D- Mesh nodes:                               " << M_Mesh->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_Mesh->numElements() << std::endl;

    // Miscellaneous
    M_postprocessingDirectory= dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile     = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_verbose                = dataFile( ( section + "/miscellaneous/verbose"                        ).data(), 1 );
    M_UW                     = dataFile( ( section + "/miscellaneous/alternate_solver"               ).data(), false );
    M_inertial_wall          = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelastic_wall      = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_linearize_string_model = dataFile( ( section + "/miscellaneous/linearize_string_model"         ).data(), true );
    M_linearize_equations    = dataFile( ( section + "/miscellaneous/linearize_equations"            ).data(), false );
    M_longitudinal_wall      = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );
    M_flux_second_der        = dataFile( ( section + "/miscellaneous/compute_flux_second_derivative" ).data(), false );
    M_dP_dt_steps            = dataFile( ( section + "/miscellaneous/pressure_derivative_steps"      ).data(), 1 );
    M_CFLmax                 = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3)/3. );

    if ( M_CFLmax > std::sqrt(3)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Initialize
    std::map< std::string, OneD_Initialize > initializeMap;
    initializeMap["A"]       = OneD_InitializeArea;
    initializeMap["Q"]       = OneD_InitializeFlux;
    initializeMap["W1"]      = OneD_InitializeRiemann1;
    initializeMap["W2"]      = OneD_InitializeRiemann2;
    initializeMap["P"]       = OneD_InitializePressure;

    M_initialVariable        = initializeMap[ dataFile( ( section + "/initialize/variable"           ).data(), "P" ) ];
    M_initialValue           = dataFile( ( section + "/initialize/initialValue"                      ).data(), 0. );
    M_restValue              = dataFile( ( section + "/initialize/restValue"                         ).data(), 0. );
    M_multiplier             = dataFile( ( section + "/initialize/multiplier"                        ).data(), 1. );

    // Physical Parameters
    M_ComputeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );
    M_PowerlawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_Density                = dataFile( ( section + "/1d_physics/density"                           ).data(), 1. );
    M_Viscosity              = dataFile( ( section + "/1d_physics/viscosity"                         ).data(), 0.035 );

    M_DensityWall            = dataFile( ( section + "/PhysicalParameters/densityWall"               ).data(), 1. );
    M_ThickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_Thickness              = dataFile( ( section + "/1d_physics/thickness"                         ).data(), 0. );
    M_Young                  = dataFile( ( section + "/1d_physics/young"                             ).data(), 4.0E6 );
    M_Poisson                = dataFile( ( section + "/1d_physics/ksi"                               ).data(), 0.5 );

    M_ViscoelasticModulus    = dataFile( ( section + "/PhysicalParameters/gamma"                     ).data(), 0. );
    M_InertialModulus        = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_RobertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_ComputeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );

    M_Area0.resize( M_Mesh->numPoints() );
    M_dArea0dz.resize( M_Mesh->numPoints() );
    M_PressBeta0.resize( M_Mesh->numPoints() );
    M_dPressBeta0dz.resize( M_Mesh->numPoints() );
    M_PressBeta1.resize( M_Mesh->numPoints() );
    M_dPressBeta1dz.resize( M_Mesh->numPoints() );
    M_FrictionKr.resize( M_Mesh->numPoints() );

    M_AlphaCoriolis.resize( M_Mesh->numPoints() );
    M_dAlphaCoriolisdz.resize( M_Mesh->numPoints() );

    M_Flux11.resize( M_Mesh->numPoints() );
    M_Flux12.resize( M_Mesh->numPoints() );
    M_Flux21.resize( M_Mesh->numPoints() );
    M_Flux22.resize( M_Mesh->numPoints() );
    M_Celerity1.resize( M_Mesh->numPoints() );
    M_Celerity2.resize( M_Mesh->numPoints() );
    M_Celerity1LeftEigenvector1.resize( M_Mesh->numPoints() );
    M_Celerity1LeftEigenvector2.resize( M_Mesh->numPoints() );
    M_Celerity2LeftEigenvector1.resize( M_Mesh->numPoints() );
    M_Celerity2LeftEigenvector2.resize( M_Mesh->numPoints() );
    M_Source10.resize( M_Mesh->numPoints() );
    M_Source20.resize( M_Mesh->numPoints() );
    M_Source11.resize( M_Mesh->numPoints() );
    M_Source12.resize( M_Mesh->numPoints() );
    M_Source21.resize( M_Mesh->numPoints() );
    M_Source22.resize( M_Mesh->numPoints() );

    // Linear Parameters
    for ( UInt i = 0; i < M_Mesh->numPoints() ; ++i )
    {
        // Physical Parameters
        if ( M_ComputeCoefficients )
            M_Area0[i]                 = std::pow( dataFile( ( section + "/1d_physics/radius"                  ).data(), 0.5 ), 2 ) * M_PI;
        else
            M_Area0[i]                 = dataFile( ( section + "/parameters/Area0"                   ).data(), M_PI );
        M_dArea0dz[i]                  = 0;
        M_PressBeta0[i]                = dataFile( ( section + "/parameters/beta0"                   ).data(), 1.e6 );
        M_dPressBeta0dz[i]             = 0.;
        M_PressBeta1[i]                = dataFile( ( section + "/parameters/beta1"                   ).data(), 0.5 );
        M_dPressBeta1dz[i]             = 0.;
        M_FrictionKr[i]                = dataFile( ( section + "/parameters/Kr"                      ).data(), 1. );

        M_AlphaCoriolis[i]             = dataFile( ( section + "/parameters/alphaCor"                ).data(), 1. / M_RobertsonCorrection );
        M_dAlphaCoriolisdz[i]          = 0.;

        // Linear Parameters
        M_Flux11[i]                    = dataFile( ( section + "/LinearParameters/Flux11"            ).data(), 1. );
        M_Flux12[i]                    = dataFile( ( section + "/LinearParameters/Flux12"            ).data(), 0. );
        M_Flux21[i]                    = dataFile( ( section + "/LinearParameters/Flux21"            ).data(), 0. );
        M_Flux22[i]                    = dataFile( ( section + "/LinearParameters/Flux22"            ).data(), 1. );
        M_Celerity1[i]                 = dataFile( ( section + "/LinearParameters/Celerity1"         ).data(), 1. );
        M_Celerity2[i]                 = dataFile( ( section + "/LinearParameters/Celerity2"         ).data(), 1. );
        M_Celerity1LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector11" ).data(), 1. );
        M_Celerity1LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector12" ).data(), 0. );
        M_Celerity2LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector21" ).data(), 0. );
        M_Celerity2LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector22" ).data(), 1. );
        M_Source10[i]                  = dataFile( ( section + "/LinearParameters/Source10"          ).data(), 0. );
        M_Source20[i]                  = dataFile( ( section + "/LinearParameters/Source20"          ).data(), 0. );
        M_Source11[i]                  = dataFile( ( section + "/LinearParameters/Source11"          ).data(), 0. );
        M_Source12[i]                  = dataFile( ( section + "/LinearParameters/Source12"          ).data(), 0. );
        M_Source21[i]                  = dataFile( ( section + "/LinearParameters/Source21"          ).data(), 0. );
        M_Source22[i]                  = dataFile( ( section + "/LinearParameters/Source22"          ).data(), 0. );
    }

    UpdateCoefficients();
}

void
OneDimensionalModel_Data::UpdateCoefficients()
{
    if ( M_ComputeCoefficients )
    {
        Real radius(0);
        Real profileIntegral(0);
        for ( UInt i = 0; i < M_Mesh->numPoints(); ++i )
        {
            // PowerlawProfile: s(r) = (1+2/gamma)*(1-r^gamma)
            radius = 1.; //std::sqrt( M_Area0[i] / M_PI );

            // Compute Coriolis Coefficient: Alpha = 2*pi/Area0*Int(s(r)^2)
            profileIntegral    =     std::pow(1+2./M_PowerlawCoefficient, 2) *
                                 (   std::pow(radius,2) / 2 +
                                     std::pow(radius,2*M_PowerlawCoefficient+2) / (2*M_PowerlawCoefficient+2) -
                                   2*std::pow(radius,  M_PowerlawCoefficient+2) / (  M_PowerlawCoefficient+2) );
            M_AlphaCoriolis[i] = 2 / std::pow(radius,2) * profileIntegral;

            // Compute Friction Coefficient: Kr = -2*pi*mu/rho*s'(R)
            M_FrictionKr[i]    = 2 * M_PI * M_Viscosity / M_Density * ( M_PowerlawCoefficient + 2 ) * std::pow( radius, M_PowerlawCoefficient - 1 );

            // Compute Beta0
            if( M_ThickVessel ) // see Amadori, Ferrari, Formaggia (MOX report 86)
            {
                M_PressBeta0[i] = - M_Thickness * M_Young * std::sqrt(M_PI) / ( std::sqrt( M_Area0[i] ) * ( (1 - M_Poisson * M_Poisson )
                                  + M_Poisson * ( 1 + M_Poisson ) * ( M_Thickness * std::sqrt(M_PI) / std::sqrt( M_Area0[i] ) ) ) );
                M_PressBeta1[i] = - 0.5;
            }
            else
            {
                M_PressBeta0[i] = M_Thickness * M_Young * std::sqrt( M_PI ) / ( std::sqrt( M_Area0[i] ) * (1 - M_Poisson * M_Poisson) );
                M_PressBeta1[i] = 0.5;
            }
        }
    }
}

void
OneDimensionalModel_Data::initLinearParam( const GetPot& /*dataFile*/ )  // CHECK THIS!!!
{
    /*
      The linearization of Euler model yields

      F = [ Q; A * (c_L)^2];

      B = [ 0; k_R / A0];

      c_L = sqrt( beta0 * beta1 / rho );
    */
    for( UInt indz=0; indz < M_Mesh->numPoints(); ++indz )
    {
        M_Celerity1[indz] = std::sqrt( M_PressBeta0[indz] * M_PressBeta1[indz] / M_Density );
        M_Flux21[indz]    = std::pow( M_Celerity1[indz], 2 );
        M_Source22[indz]  = M_FrictionKr[indz] / M_Area0(indz);
    }

    M_Celerity2 = - M_Celerity1;

    M_Celerity1LeftEigenvector1 = M_Celerity1;
    M_Celerity1LeftEigenvector2 = ScalarVector(M_Mesh->numPoints(), 1.);
    M_Celerity2LeftEigenvector1 = M_Celerity2;
    M_Celerity2LeftEigenvector2 = ScalarVector(M_Mesh->numPoints(), 1.);

    M_Flux11 = ZeroVector(M_Mesh->numPoints());
    M_Flux12 = ScalarVector(M_Mesh->numPoints(), 1.);
    M_Flux22 = ZeroVector(M_Mesh->numPoints());

    M_Source10 = ZeroVector(M_Mesh->numPoints());
    M_Source11 = ZeroVector(M_Mesh->numPoints());
    M_Source12 = ZeroVector(M_Mesh->numPoints());

    M_Source20 = ZeroVector(M_Mesh->numPoints());
    M_Source21 = ZeroVector(M_Mesh->numPoints());
}

void
OneDimensionalModel_Data::showMe( std::ostream& output ) const
{
    // Model
    output << "\n*** Values for data [Model]\n\n";
    output << "Physics Type           = " << Enum2String( M_PhysicsType, OneDimensionalModel_PhysicsMap ) << std::endl;
    output << "Flux Type              = " << Enum2String( M_FluxType,    OneDimensionalModel_FluxMap    ) << std::endl;
    output << "Source Type            = " << Enum2String( M_SourceType,  OneDimensionalModel_SourceMap  ) << std::endl;

    // Time
    output << "\n*** Values for data [time_discretization]\n\n";
    M_Time->showMe( output );

    // Space Discretization
    output << "\n*** Values for data [space_discretization]\n\n";
    output << "Length                 = " << Length() << std::endl;
    output << "NumberOfElements       = " << M_Mesh->numElements() << std::endl;

    // Miscellaneous
    output << "\n*** Values for data [miscellaneous]\n\n";
    output << "Postprocessing Dir.    = " << M_postprocessingDirectory << std::endl;
    output << "Postprocessing File    = " << M_postprocessingFile << std::endl;
    output << "verbose                = " << M_verbose << std::endl;
    output << "UW                     = " << M_UW << std::endl;
    output << "Use Inertial Wall      = " << M_inertial_wall << std::endl;
    output << "Use Viscoelastic Wall  = " << M_viscoelastic_wall << std::endl;
    output << "Linearize Model        = " << M_linearize_string_model << std::endl;
    output << "Linearize Equations    = " << M_linearize_equations << std::endl;
    output << "Longitudinal Wall      = " << M_longitudinal_wall << std::endl;
    output << "Flux Second Derivative = " << M_flux_second_der << std::endl;
    output << "Pressure Derivative    = " << M_dP_dt_steps << std::endl;
    output << "Maximum admissible CFL = " << M_CFLmax << std::endl;

    // Initialize
    output << "\n*** Values for data [initialize]\n\n";
    output << "Initial Variable       = " << M_initialVariable << std::endl;
    output << "Initial Value          = " << M_initialValue << std::endl;
    output << "Rest Value             = " << M_restValue << std::endl;
    output << "Multiplier             = " << M_multiplier << std::endl;

    // Physical Parameters
    output << "\n*** Values for data [PhysicalParameters]\n\n";
    output << "Compute Coefficients   = " << M_ComputeCoefficients << "\n";
    output << "Powerlaw Coefficient   = " << M_PowerlawCoefficient << "\n";
    output << "Fluid density          = " << M_Density << "\n";
    output << "Fluid dyn. viscosity   = " << M_Viscosity << "\n";
    output << "Wall density           = " << M_DensityWall << "\n";
    output << "Thick vessel           = " << M_ThickVessel << "\n";
    output << "Thickness              = " << M_Thickness << "\n";
    output << "Young modulus          = " << M_Young << "\n";
    output << "Poisson                = " << M_Poisson << "\n";
    output << "Viscoelastic modulus   = " << M_ViscoelasticModulus << "\n";
    output << "Inertial modulus       = " << M_InertialModulus << "\n";
    output << "Robertson correction   = " << M_RobertsonCorrection << "\n";

    output << "Area0                  = " << M_Area0 << "\n";
    output << "dArea0dz               = " << M_dArea0dz << "\n";
    output << "Beta0                  = " << M_PressBeta0 << "\n";
    output << "dBeta0dz               = " << M_dPressBeta0dz << "\n";
    output << "Beta1                  = " << M_PressBeta1 << "\n";
    output << "dBeta1dz               = " << M_dPressBeta1dz << "\n";

    output << "Alpha (Coriolis)       = " << M_AlphaCoriolis << "\n";
    output << "dAlpha (Coriolis)      = " << M_dAlphaCoriolisdz << "\n";

    output << "Friction               = " << M_FrictionKr << "\n";

    // Linear Parameters
    output << "\n*** Values for data [LinearParameters]\n\n";
    output << "Flux11                 = " << M_Flux11 << "\n";
    output << "Flux12                 = " << M_Flux12 << "\n";
    output << "Flux21                 = " << M_Flux21 << "\n";
    output << "Flux22                 = " << M_Flux22 << "\n";
    output << "Celerity1              = " << M_Celerity1 << "\n";
    output << "Celerity2              = " << M_Celerity2 << "\n";
    output << "Eigenvector11          = " << M_Celerity1LeftEigenvector1 << "\n";
    output << "Eigenvector12          = " << M_Celerity1LeftEigenvector2 << "\n";
    output << "Eigenvector21          = " << M_Celerity2LeftEigenvector1 << "\n";
    output << "Eigenvector22          = " << M_Celerity2LeftEigenvector2 << "\n";
    output << "Source10               = " << M_Source10 << "\n";
    output << "Source20               = " << M_Source20 << "\n";
    output << "Source11               = " << M_Source11 << "\n";
    output << "Source12               = " << M_Source12 << "\n";
    output << "Source21               = " << M_Source21 << "\n";
    output << "Source22               = " << M_Source22 << "\n";
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Data::setPostprocessingDirectory( const std::string& directory )
{
    M_postprocessingDirectory = directory;
}

void
OneDimensionalModel_Data::setPostprocessingFile( const std::string& file )
{
    M_postprocessingFile = file;
}

void
OneDimensionalModel_Data::setDataTime( const Time_PtrType DataTime )
{
    M_Time = DataTime;
}

void
OneDimensionalModel_Data::setDensity( const Real& Density )
{
    M_Density = Density;
}

void
OneDimensionalModel_Data::setViscosity( const Real& Viscosity )
{
    M_Viscosity = Viscosity;
}

void
OneDimensionalModel_Data::setDensityWall( const Real& DensityWall )
{
    M_DensityWall = DensityWall;
}

void
OneDimensionalModel_Data::setThickness( const Real& Thickness )
{
    M_Thickness = Thickness;
}

void
OneDimensionalModel_Data::setYoung( const Real& Young )
{
    M_Young = Young;
}

void
OneDimensionalModel_Data::setPoisson( const Real& Poisson )
{
    M_Poisson = Poisson;
}

void
OneDimensionalModel_Data::setArea0( const Real& Area0, const UInt& i )
{
    M_Area0[i] = Area0;
}

void
OneDimensionalModel_Data::setBeta0( const Real& Beta0, const UInt& i )
{
    M_PressBeta0[i] = Beta0;
}

void
OneDimensionalModel_Data::setdBeta0dz( const Real& dBeta0dz, const UInt& i )
{
    M_dPressBeta0dz[i] = dBeta0dz;
}

// ===================================================
// Get Methods
// ===================================================
const OneDimensionalModel_PhysicsTypes&
OneDimensionalModel_Data::PhysicsType() const
{
    return M_PhysicsType;
}

const OneDimensionalModel_FluxTypes&
OneDimensionalModel_Data::FluxType() const
{
    return M_FluxType;
}

const OneDimensionalModel_SourceTypes&
OneDimensionalModel_Data::SourceType() const
{
    return M_SourceType;
}

OneDimensionalModel_Data::Time_PtrType
OneDimensionalModel_Data::dataTime( void ) const
{
    return M_Time;
}

OneDimensionalModel_Data::Mesh_PtrType
OneDimensionalModel_Data::mesh() const
{
    return M_Mesh;
}

Real
OneDimensionalModel_Data::Length() const
{
    return M_Mesh->pointList( M_Mesh->numVertices() ).x() - M_Mesh->pointList( 1 ).x();
}

UInt
OneDimensionalModel_Data::NumberOfElements() const
{
    return M_Mesh->numElements();
}

UInt
OneDimensionalModel_Data::NumberOfNodes() const
{
    return M_Mesh->numPoints();
}

const std::string&
OneDimensionalModel_Data::postprocessingDirectory() const
{
    return M_postprocessingDirectory;
}

const std::string&
OneDimensionalModel_Data::postprocessingFile() const
{
    return M_postprocessingFile;
}

const int&
OneDimensionalModel_Data::verbose() const
{
    return M_verbose;
}

const bool&
OneDimensionalModel_Data::UW() const
{
    return M_UW;
}

const bool&
OneDimensionalModel_Data::inertialWall() const
{
    return M_inertial_wall;
}

const bool&
OneDimensionalModel_Data::viscoelasticWall() const
{
    return M_viscoelastic_wall;
}

const bool&
OneDimensionalModel_Data::linearizeStringModel() const
{
    return M_linearize_string_model;
}

const bool&
OneDimensionalModel_Data::linearizeEquations() const
{
    return M_linearize_equations;
}

const bool&
OneDimensionalModel_Data::longitudinalWall() const
{
    return M_longitudinal_wall;
}

const bool&
OneDimensionalModel_Data::fluxSecondDer() const
{
    return M_flux_second_der;
}

const int&
OneDimensionalModel_Data::DPdtSteps() const
{
    return M_dP_dt_steps;
}

const Real&
OneDimensionalModel_Data::CFLmax() const
{
    return M_CFLmax;
}

const OneD_Initialize&
OneDimensionalModel_Data::initialVariable() const
{
    return M_initialVariable;
}

const Real&
OneDimensionalModel_Data::initialValue() const
{
    return M_initialValue;
}

const Real&
OneDimensionalModel_Data::restValue() const
{
    return M_restValue;
}

const Real&
OneDimensionalModel_Data::multiplier() const
{
    return M_multiplier;
}

// ===================================================
// Get Methods - Physical Parameters
// ===================================================
const Real&
OneDimensionalModel_Data::DensityRho() const
{
    return M_Density;
}

const Real&
OneDimensionalModel_Data::Viscosity() const
{
    return M_Viscosity;
}

const Real&
OneDimensionalModel_Data::DensityWall() const
{
    return M_DensityWall;
}

const Real&
OneDimensionalModel_Data::Thickness() const
{
    return M_Thickness;
}

const Real&
OneDimensionalModel_Data::Young() const
{
    return M_Young;
}

const Real&
OneDimensionalModel_Data::Poisson() const
{
    return M_Poisson;
}

const Real&
OneDimensionalModel_Data::ViscoelasticModulus() const
{
    return M_ViscoelasticModulus;
}

const Real&
OneDimensionalModel_Data::InertialModulus() const
{
    return M_InertialModulus;
}

const Real&
OneDimensionalModel_Data::RobertsonCorrection() const
{
    return M_RobertsonCorrection;
}

const Real&
OneDimensionalModel_Data::Area0( const UInt& i ) const
{
    return M_Area0[i];
}

const Real&
OneDimensionalModel_Data::Beta0( const UInt& i ) const
{
    return M_PressBeta0[i];
}

const Real&
OneDimensionalModel_Data::Beta1( const UInt& i ) const
{
    return M_PressBeta1[i];
}

const Real&
OneDimensionalModel_Data::dArea0dz( const UInt& i ) const
{
    return M_dArea0dz[i];
}

const Real&
OneDimensionalModel_Data::dBeta0dz( const UInt& i ) const
{
    return M_dPressBeta0dz[i];
}

const Real&
OneDimensionalModel_Data::dBeta1dz( const UInt& i ) const
{
    return M_dPressBeta1dz[i];
}

const Real&
OneDimensionalModel_Data::AlphaCor( const UInt& i ) const
{
    return M_AlphaCoriolis[i];
}

const Real&
OneDimensionalModel_Data::dAlphaCordz( const UInt& i ) const
{
    return M_dAlphaCoriolisdz[i];
}

const Real&
OneDimensionalModel_Data::FrictionKr( const UInt& i ) const
{
    return M_FrictionKr[i];
}

// ===================================================
// Get Methods - Linear Parameters
// ===================================================
const Real&
OneDimensionalModel_Data::Flux11( const UInt& i ) const
{
    return M_Flux11[i];
}
const Real&
OneDimensionalModel_Data::Flux12( const UInt& i ) const
{
    return M_Flux12[i];
}
const Real&
OneDimensionalModel_Data::Flux21( const UInt& i ) const
{
    return M_Flux21[i];
}
const Real&
OneDimensionalModel_Data::Flux22( const UInt& i ) const
{
    return M_Flux22[i];
}

const Real&
OneDimensionalModel_Data::Celerity1( const UInt& i ) const
{
    return M_Celerity1[i];
}

const Real&
OneDimensionalModel_Data::Celerity2( const UInt& i ) const
{
    return M_Celerity2[i];
}

const Real&
OneDimensionalModel_Data::LeftEigenVector11( const UInt& i ) const
{
    return M_Celerity1LeftEigenvector1[i];
}

const Real&
OneDimensionalModel_Data::LeftEigenVector12( const UInt& i ) const
{
    return M_Celerity1LeftEigenvector2[i];
}

const Real&
OneDimensionalModel_Data::LeftEigenVector21( const UInt& i ) const
{
    return M_Celerity2LeftEigenvector1[i];
}

const Real&
OneDimensionalModel_Data::LeftEigenVector22( const UInt& i ) const
{
    return M_Celerity2LeftEigenvector2[i];
}

const Real&
OneDimensionalModel_Data::Source10( const UInt& i ) const
{
    return M_Source10[i];
}
const Real&
OneDimensionalModel_Data::Source20( const UInt& i ) const
{
    return M_Source20[i];
}
const Real&
OneDimensionalModel_Data::Source11( const UInt& i ) const
{
    return M_Source11[i];
}
const Real&
OneDimensionalModel_Data::Source12( const UInt& i ) const
{
    return M_Source12[i];
}
const Real&
OneDimensionalModel_Data::Source21( const UInt& i ) const
{
    return M_Source21[i];
}
const Real&
OneDimensionalModel_Data::Source22( const UInt& i ) const
{
    return M_Source22[i];
}

}
