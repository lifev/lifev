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
 *  @brief File containing a class for 1D model data handling.
 *
 *  @version 1.0
 *  @date 01-07-2004
 *  @author Vincent Martin
 *
 *  @version 2.0
 *  @date 12-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifemc/lifesolver/OneDimensionalData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
OneDimensionalData::OneDimensionalData():
    M_physicsType               (),
    M_fluxType                  (),
    M_sourceType                (),
    M_time                      (),
    M_mesh                      ( new mesh_Type() ),
    M_postprocessingDirectory   (),
    M_postprocessingFile        (),
    M_verbose                   (),
    M_inertialWall              (),
    M_viscoelasticWall          (),
    M_linearizeStringModel      (),
    M_linearizeEquations        (),
    M_longitudinalWall          (),
    M_fluxSecondDer             (),
    M_dP_dt_steps               (),
    M_CFLmax                    (),
    M_jacobianPerturbationArea  (),
    M_jacobianPerturbationFlowRate(),
    M_jacobianPerturbationPressure(),
    M_computeCoefficients       (),
    M_powerLawCoefficient       (),
    M_density                   (),
    M_viscosity                 (),
    M_densityWall               (),
    M_thickVessel               (),
    M_young                     (),
    M_poisson                   (),
    M_externalPressure          (),
    M_viscoelasticModulus       (),
    M_inertialModulus           (),
    M_robertsonCorrection       (),
    M_thickness                 (),
    M_friction                  (),
    M_area0                     (),
    M_alpha                     (),
    M_beta0                     (),
    M_beta1                     (),
    M_dArea0dz                  (),
    M_dAlphadz                  (),
    M_dBeta0dz                  (),
    M_dBeta1dz                  (),
    M_flux11                    (),
    M_flux12                    (),
    M_flux21                    (),
    M_flux22                    (),
    M_celerity1                 (),
    M_celerity2                 (),
    M_celerity1LeftEigenvector1 (),
    M_celerity1LeftEigenvector2 (),
    M_celerity2LeftEigenvector1 (),
    M_celerity2LeftEigenvector2 (),
    M_source10                  (),
    M_source20                  (),
    M_source11                  (),
    M_source12                  (),
    M_source21                  (),
    M_source22                  ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalData::setup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_physicsType = OneDimensional::physicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_fluxType    = OneDimensional::fluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_sourceType  = OneDimensional::sourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    M_mesh->setup( dataFile( ( section + "/space_discretization/Length"           ).data(), 1. ),
                   dataFile( ( section + "/space_discretization/NumberOfElements" ).data(), 10 ) );

    //std::cout << " 1D- Mesh nodes:                               " << M_mesh->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_mesh->numElements() << std::endl;

    // Miscellaneous
    M_postprocessingDirectory= dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile     = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_verbose                = dataFile( ( section + "/miscellaneous/verbose"                        ).data(), 1 );
    M_inertialWall          = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelasticWall      = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_linearizeStringModel = dataFile( ( section + "/miscellaneous/linearize_string_model"         ).data(), true );
    M_linearizeEquations    = dataFile( ( section + "/miscellaneous/linearize_equations"            ).data(), false );
    M_longitudinalWall      = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );
    M_fluxSecondDer        = dataFile( ( section + "/miscellaneous/compute_flux_second_derivative" ).data(), false );
    M_dP_dt_steps            = dataFile( ( section + "/miscellaneous/pressure_derivative_steps"      ).data(), 1 );
    M_CFLmax                 = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3)/3. );

    if ( M_CFLmax > std::sqrt(3)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Jacobian perturbation
    M_jacobianPerturbationArea     = dataFile( ( section + "/JacobianPerturbation/deltaArea"         ).data(), 0.001 );
    M_jacobianPerturbationFlowRate = dataFile( ( section + "/JacobianPerturbation/deltaFlowRate"     ).data(), 0.001 );
    M_jacobianPerturbationPressure = dataFile( ( section + "/JacobianPerturbation/deltaPressure"     ).data(), 1 );

    // Physical Parameters
    M_computeCoefficients    = dataFile( ( section + "/PhysicalParameters/ComputeCoefficients"       ).data(), false );
    M_powerLawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_density                = dataFile( ( section + "/PhysicalParameters/density"                   ).data(), 1. );
    M_viscosity              = dataFile( ( section + "/PhysicalParameters/viscosity"                 ).data(), 0.035 );

    M_densityWall            = dataFile( ( section + "/PhysicalParameters/densityWall"               ).data(), 1. );
    M_thickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_young                  = dataFile( ( section + "/PhysicalParameters/young"                     ).data(), 4.0E6 );
    M_poisson                = dataFile( ( section + "/PhysicalParameters/poisson"                   ).data(), 0.5 );

    M_externalPressure       = dataFile( ( section + "/PhysicalParameters/externalPressure"          ).data(), 0 );
    M_friction               = dataFile( ( section + "/PhysicalParameters/Kr"                        ).data(), 1. );

    M_viscoelasticModulus    = dataFile( ( section + "/PhysicalParameters/gamma"                     ).data(), 0. );
    M_inertialModulus        = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_robertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_computeCoefficients    = dataFile( ( section + "/PhysicalParameters/ComputeCoefficients"       ).data(), false );

    // Physical Parameters defined along the 1D model
    M_thickness.resize( M_mesh->numPoints() );

    M_area0.resize( M_mesh->numPoints() );
    M_beta0.resize( M_mesh->numPoints() );
    M_beta1.resize( M_mesh->numPoints() );
    M_alpha.resize( M_mesh->numPoints() );

    M_dArea0dz.resize( M_mesh->numPoints() );
    M_dBeta0dz.resize( M_mesh->numPoints() );
    M_dBeta1dz.resize( M_mesh->numPoints() );
    M_dAlphadz.resize( M_mesh->numPoints() );

    // Linear Parameters defined along the 1D model
    M_flux11.resize( M_mesh->numPoints() );
    M_flux12.resize( M_mesh->numPoints() );
    M_flux21.resize( M_mesh->numPoints() );
    M_flux22.resize( M_mesh->numPoints() );
    M_celerity1.resize( M_mesh->numPoints() );
    M_celerity2.resize( M_mesh->numPoints() );
    M_celerity1LeftEigenvector1.resize( M_mesh->numPoints() );
    M_celerity1LeftEigenvector2.resize( M_mesh->numPoints() );
    M_celerity2LeftEigenvector1.resize( M_mesh->numPoints() );
    M_celerity2LeftEigenvector2.resize( M_mesh->numPoints() );
    M_source10.resize( M_mesh->numPoints() );
    M_source20.resize( M_mesh->numPoints() );
    M_source11.resize( M_mesh->numPoints() );
    M_source12.resize( M_mesh->numPoints() );
    M_source21.resize( M_mesh->numPoints() );
    M_source22.resize( M_mesh->numPoints() );

    // Select a law for the coefficients
    std::map< string, OneD_distributionLaw > distributionLawMap;
    distributionLawMap["none"]      = none;
    distributionLawMap["linear"]    = linear;
    distributionLawMap["pointwise"] = pointwise;

    OneD_distributionLaw distributionLaw = distributionLawMap[ dataFile( ( section + "/PhysicalParameters/DistributionLaw" ).data(), "none" ) ];
    switch ( distributionLaw )
    {
    case none:

        for ( UInt i = 0; i < M_mesh->numPoints() ; ++i )
        {
            // Physical Parameters
            M_thickness[i]                 = dataFile( ( section + "/PhysicalParameters/thickness"       ).data(), 0. );

            M_area0[i]                     = dataFile( ( section + "/PhysicalParameters/Area0"           ).data(), M_PI );
            M_alpha[i]                     = dataFile( ( section + "/PhysicalParameters/AlphaCoriolis"   ).data(), 1. / M_robertsonCorrection );
            M_beta0[i]                     = dataFile( ( section + "/PhysicalParameters/Beta0"           ).data(), 1.e6 );
            M_beta1[i]                     = dataFile( ( section + "/PhysicalParameters/Beta1"           ).data(), 0.5 );

            // Linear Parameters
            M_flux11[i]                    = dataFile( ( section + "/LinearParameters/Flux11"            ).data(), 1. );
            M_flux12[i]                    = dataFile( ( section + "/LinearParameters/Flux12"            ).data(), 0. );
            M_flux21[i]                    = dataFile( ( section + "/LinearParameters/Flux21"            ).data(), 0. );
            M_flux22[i]                    = dataFile( ( section + "/LinearParameters/Flux22"            ).data(), 1. );
            M_celerity1[i]                 = dataFile( ( section + "/LinearParameters/Celerity1"         ).data(), 1. );
            M_celerity2[i]                 = dataFile( ( section + "/LinearParameters/Celerity2"         ).data(), 1. );
            M_celerity1LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector11" ).data(), 1. );
            M_celerity1LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector12" ).data(), 0. );
            M_celerity2LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector21" ).data(), 0. );
            M_celerity2LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector22" ).data(), 1. );
            M_source10[i]                  = dataFile( ( section + "/LinearParameters/Source10"          ).data(), 0. );
            M_source20[i]                  = dataFile( ( section + "/LinearParameters/Source20"          ).data(), 0. );
            M_source11[i]                  = dataFile( ( section + "/LinearParameters/Source11"          ).data(), 0. );
            M_source12[i]                  = dataFile( ( section + "/LinearParameters/Source12"          ).data(), 0. );
            M_source21[i]                  = dataFile( ( section + "/LinearParameters/Source21"          ).data(), 0. );
            M_source22[i]                  = dataFile( ( section + "/LinearParameters/Source22"          ).data(), 0. );
        }

        break;

    case linear:

        linearInterpolation( M_thickness, dataFile, section + "/PhysicalParameters/thickness", 0. );

        linearInterpolation( M_area0, dataFile, section + "/PhysicalParameters/Area0", M_PI );
        linearInterpolation( M_alpha, dataFile, section + "/PhysicalParameters/AlphaCoriolis", 1. / M_robertsonCorrection, true );
        linearInterpolation( M_beta0, dataFile, section + "/PhysicalParameters/Beta0", 1.e6 );
        linearInterpolation( M_beta1, dataFile, section + "/PhysicalParameters/Beta1", 0.5 );

        linearInterpolation( M_flux11, dataFile, section + "/PhysicalParameters/Flux11", 1. );
        linearInterpolation( M_flux12, dataFile, section + "/PhysicalParameters/Flux12", 0. );
        linearInterpolation( M_flux21, dataFile, section + "/PhysicalParameters/Flux21", 0. );
        linearInterpolation( M_flux22, dataFile, section + "/PhysicalParameters/Flux22", 1. );
        linearInterpolation( M_celerity1, dataFile, section + "/PhysicalParameters/Celerity1", 1. );
        linearInterpolation( M_celerity2, dataFile, section + "/PhysicalParameters/Celerity2", 1. );
        linearInterpolation( M_celerity1LeftEigenvector1, dataFile, section + "/PhysicalParameters/LeftEigenvector11", 1. );
        linearInterpolation( M_celerity1LeftEigenvector2, dataFile, section + "/PhysicalParameters/LeftEigenvector12", 0. );
        linearInterpolation( M_celerity2LeftEigenvector1, dataFile, section + "/PhysicalParameters/LeftEigenvector21", 0. );
        linearInterpolation( M_celerity2LeftEigenvector2, dataFile, section + "/PhysicalParameters/LeftEigenvector22", 1. );
        linearInterpolation( M_source10, dataFile, section + "/PhysicalParameters/Source10", 0. );
        linearInterpolation( M_source20, dataFile, section + "/PhysicalParameters/Source20", 0. );
        linearInterpolation( M_source11, dataFile, section + "/PhysicalParameters/Source11", 0. );
        linearInterpolation( M_source12, dataFile, section + "/PhysicalParameters/Source12", 0. );
        linearInterpolation( M_source21, dataFile, section + "/PhysicalParameters/Source21", 0. );
        linearInterpolation( M_source22, dataFile, section + "/PhysicalParameters/Source22", 0. );

        break;

    case pointwise:

        for ( UInt i = 0; i < M_mesh->numPoints() ; ++i )
        {
            // Physical Parameters
            M_thickness[i]                 = dataFile( ( section + "/PhysicalParameters/thickness"       ).data(), 0., i );

            M_area0[i]                     = dataFile( ( section + "/PhysicalParameters/Area0"           ).data(), M_PI, i );
            M_alpha[i]                     = dataFile( ( section + "/PhysicalParameters/AlphaCoriolis"   ).data(), 1. / M_robertsonCorrection, i );
            M_beta0[i]                     = dataFile( ( section + "/PhysicalParameters/Beta0"           ).data(), 1.e6, i );
            M_beta1[i]                     = dataFile( ( section + "/PhysicalParameters/Beta1"           ).data(), 0.5, i );

            // Linear Parameters
            M_flux11[i]                    = dataFile( ( section + "/LinearParameters/Flux11"            ).data(), 1., i );
            M_flux12[i]                    = dataFile( ( section + "/LinearParameters/Flux12"            ).data(), 0., i );
            M_flux21[i]                    = dataFile( ( section + "/LinearParameters/Flux21"            ).data(), 0., i );
            M_flux22[i]                    = dataFile( ( section + "/LinearParameters/Flux22"            ).data(), 1., i );
            M_celerity1[i]                 = dataFile( ( section + "/LinearParameters/Celerity1"         ).data(), 1., i );
            M_celerity2[i]                 = dataFile( ( section + "/LinearParameters/Celerity2"         ).data(), 1., i );
            M_celerity1LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector11" ).data(), 1., i );
            M_celerity1LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector12" ).data(), 0., i );
            M_celerity2LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector21" ).data(), 0., i );
            M_celerity2LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector22" ).data(), 1., i );
            M_source10[i]                  = dataFile( ( section + "/LinearParameters/Source10"          ).data(), 0., i );
            M_source20[i]                  = dataFile( ( section + "/LinearParameters/Source20"          ).data(), 0., i );
            M_source11[i]                  = dataFile( ( section + "/LinearParameters/Source11"          ).data(), 0., i );
            M_source12[i]                  = dataFile( ( section + "/LinearParameters/Source12"          ).data(), 0., i );
            M_source21[i]                  = dataFile( ( section + "/LinearParameters/Source21"          ).data(), 0., i );
            M_source22[i]                  = dataFile( ( section + "/LinearParameters/Source22"          ).data(), 0., i );
        }

        break;

    default:
        std::cout << "Warning: taperLaw \"" << distributionLaw << "\"not available!" << std::endl;
    }

    updateCoefficients();
}

void
OneDimensionalData::oldStyleSetup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_physicsType = OneDimensional::physicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_fluxType    = OneDimensional::fluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_sourceType  = OneDimensional::sourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    Real length = dataFile( ( section + "/discretization/x_right" ).data(), 1. ) -
                  dataFile( ( section + "/discretization/x_left"  ).data(), 0. );

    M_mesh->setup( length, dataFile( ( section + "/discretization/nb_elem" ).data(), 10 ) );

    //std::cout << " 1D- Mesh nodes:                               " << M_mesh->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_mesh->numElements() << std::endl;

    // Miscellaneous
    M_postprocessingDirectory= dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile     = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_verbose                = dataFile( ( section + "/miscellaneous/verbose"                        ).data(), 1 );
    M_inertialWall          = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelasticWall      = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_linearizeStringModel = dataFile( ( section + "/miscellaneous/linearize_string_model"         ).data(), true );
    M_linearizeEquations    = dataFile( ( section + "/miscellaneous/linearize_equations"            ).data(), false );
    M_longitudinalWall      = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );
    M_fluxSecondDer        = dataFile( ( section + "/miscellaneous/compute_flux_second_derivative" ).data(), false );
    M_dP_dt_steps            = dataFile( ( section + "/miscellaneous/pressure_derivative_steps"      ).data(), 1 );
    M_CFLmax                 = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3)/3. );

    if ( M_CFLmax > std::sqrt(3)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Jacobian perturbation
    M_jacobianPerturbationArea     = dataFile( ( section + "JacobianPerturbation/deltaArea"          ).data(), 0.001 );
    M_jacobianPerturbationFlowRate = dataFile( ( section + "JacobianPerturbation/deltaFlowRate"      ).data(), 0.001 );
    M_jacobianPerturbationPressure = dataFile( ( section + "JacobianPerturbation/deltaPressure"      ).data(), 1 );

    // Physical Parameters
    M_computeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );
    M_powerLawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_density                = dataFile( ( section + "/1d_physics/density"                           ).data(), 1. );
    M_viscosity              = dataFile( ( section + "/1d_physics/viscosity"                         ).data(), 0.035 );

    M_densityWall            = dataFile( ( section + "/PhysicalParameters/densityWall"               ).data(), 1. );
    M_thickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_young                  = dataFile( ( section + "/1d_physics/young"                             ).data(), 4.0E6 );
    M_poisson                = dataFile( ( section + "/1d_physics/ksi"                               ).data(), 0.5 );

    M_externalPressure       = dataFile( ( section + "/PhysicalParameters/externalPressure"          ).data(), 0 );

    M_viscoelasticModulus    = dataFile( ( section + "/PhysicalParameters/gamma"                     ).data(), 0. );
    M_inertialModulus        = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_robertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_friction               = dataFile( ( section + "/parameters/Kr"                                ).data(), 1. );
    M_computeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );

    M_thickness.resize( M_mesh->numPoints() );

    M_area0.resize( M_mesh->numPoints() );
    M_beta0.resize( M_mesh->numPoints() );
    M_beta1.resize( M_mesh->numPoints() );
    M_alpha.resize( M_mesh->numPoints() );

    M_dArea0dz.resize( M_mesh->numPoints() );
    M_dBeta0dz.resize( M_mesh->numPoints() );
    M_dBeta1dz.resize( M_mesh->numPoints() );
    M_dAlphadz.resize( M_mesh->numPoints() );

    M_flux11.resize( M_mesh->numPoints() );
    M_flux12.resize( M_mesh->numPoints() );
    M_flux21.resize( M_mesh->numPoints() );
    M_flux22.resize( M_mesh->numPoints() );
    M_celerity1.resize( M_mesh->numPoints() );
    M_celerity2.resize( M_mesh->numPoints() );
    M_celerity1LeftEigenvector1.resize( M_mesh->numPoints() );
    M_celerity1LeftEigenvector2.resize( M_mesh->numPoints() );
    M_celerity2LeftEigenvector1.resize( M_mesh->numPoints() );
    M_celerity2LeftEigenvector2.resize( M_mesh->numPoints() );
    M_source10.resize( M_mesh->numPoints() );
    M_source20.resize( M_mesh->numPoints() );
    M_source11.resize( M_mesh->numPoints() );
    M_source12.resize( M_mesh->numPoints() );
    M_source21.resize( M_mesh->numPoints() );
    M_source22.resize( M_mesh->numPoints() );

    // Linear Parameters
    for ( UInt i = 0; i < M_mesh->numPoints() ; ++i )
    {
        // Physical Parameters
        if ( M_computeCoefficients )
            M_area0[i]                 = std::pow( dataFile( ( section + "/1d_physics/radius"        ).data(), 0.5 ), 2 ) * M_PI;
        else
            M_area0[i]                 = dataFile( ( section + "/parameters/Area0"                   ).data(), M_PI );
        M_beta0[i]                     = dataFile( ( section + "/parameters/beta0"                   ).data(), 1.e6 );
        M_beta1[i]                     = dataFile( ( section + "/parameters/beta1"                   ).data(), 0.5 );
        M_thickness[i]                 = dataFile( ( section + "/1d_physics/thickness"               ).data(), 0. );

        M_alpha[i]                     = dataFile( ( section + "/parameters/alphaCor"                ).data(), 1. / M_robertsonCorrection );

        // Linear Parameters
        M_flux11[i]                    = dataFile( ( section + "/LinearParameters/Flux11"            ).data(), 1. );
        M_flux12[i]                    = dataFile( ( section + "/LinearParameters/Flux12"            ).data(), 0. );
        M_flux21[i]                    = dataFile( ( section + "/LinearParameters/Flux21"            ).data(), 0. );
        M_flux22[i]                    = dataFile( ( section + "/LinearParameters/Flux22"            ).data(), 1. );
        M_celerity1[i]                 = dataFile( ( section + "/LinearParameters/Celerity1"         ).data(), 1. );
        M_celerity2[i]                 = dataFile( ( section + "/LinearParameters/Celerity2"         ).data(), 1. );
        M_celerity1LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector11" ).data(), 1. );
        M_celerity1LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector12" ).data(), 0. );
        M_celerity2LeftEigenvector1[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector21" ).data(), 0. );
        M_celerity2LeftEigenvector2[i] = dataFile( ( section + "/LinearParameters/LeftEigenvector22" ).data(), 1. );
        M_source10[i]                  = dataFile( ( section + "/LinearParameters/Source10"          ).data(), 0. );
        M_source20[i]                  = dataFile( ( section + "/LinearParameters/Source20"          ).data(), 0. );
        M_source11[i]                  = dataFile( ( section + "/LinearParameters/Source11"          ).data(), 0. );
        M_source12[i]                  = dataFile( ( section + "/LinearParameters/Source12"          ).data(), 0. );
        M_source21[i]                  = dataFile( ( section + "/LinearParameters/Source21"          ).data(), 0. );
        M_source22[i]                  = dataFile( ( section + "/LinearParameters/Source22"          ).data(), 0. );
    }

    updateCoefficients();
}

void
OneDimensionalData::updateCoefficients()
{
    if ( M_computeCoefficients )
    {
        // PowerlawProfile: s(r) = (1+2/gamma)*(1-r^gamma)
        Real radius (1.); //std::sqrt( M_Area0[i] / M_PI );

        Real profileIntegral =   std::pow(1+2./M_powerLawCoefficient, 2) *
                                 (   std::pow(radius,2) / 2 +
                                     std::pow(radius,2*M_powerLawCoefficient+2) / (2*M_powerLawCoefficient+2) -
                                     2*std::pow(radius,  M_powerLawCoefficient+2) / (  M_powerLawCoefficient+2) );

        // Compute Friction Coefficient: Kr = -2*pi*mu/rho*s'(R)
        M_friction = 2 * M_PI * M_viscosity / M_density * ( M_powerLawCoefficient + 2 ) * std::pow( radius, M_powerLawCoefficient - 1 );

        for ( UInt i = 0; i < M_mesh->numPoints(); ++i )
        {
            // Compute Coriolis Coefficient: Alpha = 2*pi/Area0*Int(s(r)^2)
            M_alpha[i] = 2 / std::pow(radius,2) * profileIntegral;

            // Compute Beta0
            if ( M_thickVessel ) // see Amadori, Ferrari, Formaggia (MOX report 86)
            {
                M_beta0[i] = - M_thickness[i] * M_young * std::sqrt(M_PI) / ( std::sqrt( M_area0[i] ) * ( (1 - M_poisson * M_poisson )
                                                                              + M_poisson * ( 1 + M_poisson ) * ( M_thickness[i] * std::sqrt(M_PI) / std::sqrt( M_area0[i] ) ) ) );
                M_beta1[i] = - 0.5;
            }
            else
            {
                M_beta0[i] = M_thickness[i] * std::sqrt( M_PI / M_area0[i] ) *  M_young / (1 - M_poisson * M_poisson);
                M_beta1[i] = 0.5;
            }
        }
    }

    // Compute the derivatives of the coefficients
    computeDerivatives();
}

void
OneDimensionalData::initLinearParam( const GetPot& /*dataFile*/ )  // CHECK THIS!!!
{
    /*
      The linearization of Euler model yields

      F = [ Q; A * (c_L)^2];

      B = [ 0; k_R / A0];

      c_L = sqrt( beta0 * beta1 / rho );
    */
    for ( UInt indz=0; indz < M_mesh->numPoints(); ++indz )
    {
        M_celerity1[indz] = std::sqrt( M_beta0[indz] * M_beta1[indz] / M_density );
        M_flux21[indz]    = std::pow( M_celerity1[indz], 2 );
        M_source22[indz]  = M_friction / M_area0(indz);
    }

    M_celerity2 = - M_celerity1;

    M_celerity1LeftEigenvector1 = M_celerity1;
    M_celerity1LeftEigenvector2 = ScalarVector(M_mesh->numPoints(), 1.);
    M_celerity2LeftEigenvector1 = M_celerity2;
    M_celerity2LeftEigenvector2 = ScalarVector(M_mesh->numPoints(), 1.);

    M_flux11 = ZeroVector(M_mesh->numPoints());
    M_flux12 = ScalarVector(M_mesh->numPoints(), 1.);
    M_flux22 = ZeroVector(M_mesh->numPoints());

    M_source10 = ZeroVector(M_mesh->numPoints());
    M_source11 = ZeroVector(M_mesh->numPoints());
    M_source12 = ZeroVector(M_mesh->numPoints());

    M_source20 = ZeroVector(M_mesh->numPoints());
    M_source21 = ZeroVector(M_mesh->numPoints());
}

void
OneDimensionalData::showMe( std::ostream& output ) const
{
    // Model
    //output << std::scientific << std::setprecision(15);
    output << "\n*** Values for data [Model]\n\n";
    output << "Physics Type           = " << Enum2String( M_physicsType, OneDimensional::physicsMap ) << std::endl;
    output << "Flux Type              = " << Enum2String( M_fluxType,    OneDimensional::fluxMap    ) << std::endl;
    output << "Source Type            = " << Enum2String( M_sourceType,  OneDimensional::sourceMap  ) << std::endl;

    // Time
    output << "\n*** Values for data [time_discretization]\n\n";
    M_time->showMe( output );

    // Space Discretization
    output << "\n*** Values for data [space_discretization]\n\n";
    output << "Length                 = " << length() << std::endl;
    output << "NumberOfElements       = " << M_mesh->numElements() << std::endl;

    // Miscellaneous
    output << "\n*** Values for data [miscellaneous]\n\n";
    output << "Postprocessing Dir.    = " << M_postprocessingDirectory << std::endl;
    output << "Postprocessing File    = " << M_postprocessingFile << std::endl;
    output << "verbose                = " << M_verbose << std::endl;
    output << "Use Inertial Wall      = " << M_inertialWall << std::endl;
    output << "Use Viscoelastic Wall  = " << M_viscoelasticWall << std::endl;
    output << "Linearize Model        = " << M_linearizeStringModel << std::endl;
    output << "Linearize Equations    = " << M_linearizeEquations << std::endl;
    output << "Longitudinal Wall      = " << M_longitudinalWall << std::endl;
    output << "Flux Second Derivative = " << M_fluxSecondDer << std::endl;
    output << "Pressure Derivative    = " << M_dP_dt_steps << std::endl;
    output << "Maximum admissible CFL = " << M_CFLmax << std::endl;

    // Jacobian perturbation
    output << "Jacobian perturbation Area      = " << M_jacobianPerturbationArea << std::endl;
    output << "Jacobian perturbation Flow Rate = " << M_jacobianPerturbationFlowRate << std::endl;
    output << "Jacobian perturbation Pressure  = " << M_jacobianPerturbationPressure << std::endl;

    // Physical Parameters
    output << "\n*** Values for data [PhysicalParameters]\n\n";
    output << "Compute Coefficients   = " << M_computeCoefficients << "\n";
    output << "Powerlaw Coefficient   = " << M_powerLawCoefficient << "\n";
    output << "Fluid density          = " << M_density << "\n";
    output << "Fluid dyn. viscosity   = " << M_viscosity << "\n";
    output << "Wall density           = " << M_densityWall << "\n";
    output << "Thick vessel           = " << M_thickVessel << "\n";
    output << "Young modulus          = " << M_young << "\n";
    output << "Poisson                = " << M_poisson << "\n";
    output << "External pressure      = " << M_externalPressure << "\n";
    output << "Viscoelastic modulus   = " << M_viscoelasticModulus << "\n";
    output << "Inertial modulus       = " << M_inertialModulus << "\n";
    output << "Robertson correction   = " << M_robertsonCorrection << "\n";

    output << "Area0                  = " << M_area0 << "\n";
    output << "dArea0dz               = " << M_dArea0dz << "\n";
    output << "Beta0                  = " << M_beta0 << "\n";
    output << "dBeta0dz               = " << M_dBeta0dz << "\n";
    output << "Beta1                  = " << M_beta1 << "\n";
    output << "dBeta1dz               = " << M_dBeta1dz << "\n";

    output << "Alpha (Coriolis)       = " << M_alpha << "\n";
    output << "dAlpha (Coriolis)      = " << M_dAlphadz << "\n";

    output << "Thickness              = " << M_thickness << "\n";
    output << "Friction               = " << M_friction << "\n";

    // Linear Parameters
    output << "\n*** Values for data [LinearParameters]\n\n";
    output << "Flux11                 = " << M_flux11 << "\n";
    output << "Flux12                 = " << M_flux12 << "\n";
    output << "Flux21                 = " << M_flux21 << "\n";
    output << "Flux22                 = " << M_flux22 << "\n";
    output << "Celerity1              = " << M_celerity1 << "\n";
    output << "Celerity2              = " << M_celerity2 << "\n";
    output << "Eigenvector11          = " << M_celerity1LeftEigenvector1 << "\n";
    output << "Eigenvector12          = " << M_celerity1LeftEigenvector2 << "\n";
    output << "Eigenvector21          = " << M_celerity2LeftEigenvector1 << "\n";
    output << "Eigenvector22          = " << M_celerity2LeftEigenvector2 << "\n";
    output << "Source10               = " << M_source10 << "\n";
    output << "Source20               = " << M_source20 << "\n";
    output << "Source11               = " << M_source11 << "\n";
    output << "Source12               = " << M_source12 << "\n";
    output << "Source21               = " << M_source21 << "\n";
    output << "Source22               = " << M_source22 << "\n";
}

// ===================================================
// Get Methods - Physical Parameters
// ===================================================
const Real&
OneDimensionalData::robertsonCorrection() const
{
    if ( M_robertsonCorrection != 1. )
        std::cout << "!!! WARNING: Robertson corretion has not been checked in this version of the code !!!" << std::endl;
    return M_robertsonCorrection;
}
// ===================================================
// Private methods
// ===================================================
void
OneDimensionalData::linearInterpolation( scalarVector_Type& vector,
                                               const GetPot& dataFile,
                                               const std::string& quantity,
                                               const Real& defaultValue,
                                               const bool& isArea )
{
    Real a  = dataFile( quantity.data(), defaultValue, 0 );
    Real b  = dataFile( quantity.data(), a, 1 );

    // linearInterpolation disabled as tapering is not working!
    for ( UInt i(1); i <= M_mesh->numPoints() ; ++i )
        if ( isArea )
            vector[i-1] = (a + b + 2 * std::sqrt(a*b)) / 4;
        else
            vector[i-1] = (a + b) / 2;
}

void
OneDimensionalData::computeDerivatives()
{
    Real nodes = M_mesh->numPoints();

    // We use 2° order finite differences to compute the derivatives (it is coded only for homogeneous discretizations)
    for ( UInt i = 1 ; i < nodes-1 ; ++i )
    {
        M_dArea0dz[i] = ( M_area0[i+1] - M_area0[i-1] ) / 2;
        M_dBeta0dz[i] = ( M_beta0[i+1] - M_beta0[i-1] ) / 2;
        M_dBeta1dz[i] = ( M_beta1[i+1] - M_beta1[i-1] ) / 2;
        M_dAlphadz[i] = ( M_alpha[i+1] - M_alpha[i-1] ) / 2;
    }

    M_dArea0dz[0] = ( -M_area0[2] + 4*M_area0[1] - 3*M_area0[0] ) / 2;
    M_dBeta0dz[0] = ( -M_beta0[2] + 4*M_beta0[1] - 3*M_beta0[0] ) / 2;
    M_dBeta1dz[0] = ( -M_beta1[2] + 4*M_beta1[1] - 3*M_beta1[0] ) / 2;
    M_dAlphadz[0] = ( -M_alpha[2] + 4*M_alpha[1] - 3*M_alpha[0] ) / 2;

    M_dArea0dz[nodes-1] = ( 3*M_area0[nodes-1] - 4*M_area0[nodes-2] + M_area0[nodes-3] ) / 2;
    M_dBeta0dz[nodes-1] = ( 3*M_beta0[nodes-1] - 4*M_beta0[nodes-2] + M_beta0[nodes-3] ) / 2;
    M_dBeta1dz[nodes-1] = ( 3*M_beta1[nodes-1] - 4*M_beta1[nodes-2] + M_beta1[nodes-3] ) / 2;
    M_dAlphadz[nodes-1] = ( 3*M_alpha[nodes-1] - 4*M_alpha[nodes-2] + M_alpha[nodes-3] ) / 2;

    M_dArea0dz /= M_mesh->meanH();
    M_dBeta0dz /= M_mesh->meanH();
    M_dBeta1dz /= M_mesh->meanH();
    M_dAlphadz /= M_mesh->meanH();
}

} // LifeV namespace
