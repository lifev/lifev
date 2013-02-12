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
 *  @maintainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/one_d_fsi/solver/OneDFSIData.hpp>

namespace LifeV
{

// ===================================================
// Constructors
// ===================================================
OneDFSIData::OneDFSIData():
    M_physicsType               (),
    M_fluxType                  (),
    M_sourceType                (),
    M_timeDataPtr               (),
    M_meshPtr                   ( new mesh_Type ),
    M_viscoelasticWall          (),
    M_viscoelasticAngle         (),
    M_viscoelasticPeriod        (),
    M_viscoelasticCoefficient   (),
    M_inertialWall              (),
    M_densityWall               (),
    M_inertialModulus           (),
    M_longitudinalWall          (),
    M_postprocessingDirectory   (),
    M_postprocessingFile        (),
    M_CFLmax                    (),
    M_jacobianPerturbationArea  (),
    M_jacobianPerturbationFlowRate(),
    M_jacobianPerturbationStress(),
    M_computeCoefficients       (),
    M_powerLawCoefficient       (),
    M_density                   (),
    M_viscosity                 (),
    M_thickVessel               (),
    M_young                     (),
    M_poisson                   (),
    M_externalPressure          (),
    M_venousPressure            (),
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
OneDFSIData::setup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_physicsType = OneDFSI::physicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_fluxType    = OneDFSI::fluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_sourceType  = OneDFSI::sourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_timeDataPtr.get() )
        M_timeDataPtr.reset( new time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    Real length = dataFile( ( section + "/space_discretization/Length"           ).data(), 1. );
    Real numberOfElements = dataFile( ( section + "/space_discretization/NumberOfElements" ).data(), 10 );

    regularMesh1D( *M_meshPtr, 1, numberOfElements, false, length, 0 );

    //std::cout << " 1D- Mesh nodes:                               " << M_meshPtr->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_meshPtr->numElements() << std::endl;

    // Physical Wall
    M_viscoelasticWall        = dataFile( ( section + "/PhysicalWall/ViscoelasticWall"                ).data(), false );
    M_viscoelasticAngle       = dataFile( ( section + "/PhysicalWall/ViscoelasticAngle"               ).data(), 5. ) * M_PI / 180.;
    M_viscoelasticPeriod      = dataFile( ( section + "/PhysicalWall/ViscoelasticPeriod"              ).data(), 0.3 * 0.8 );

    M_inertialWall            = dataFile( ( section + "/PhysicalWall/InertialWall"                    ).data(), false );
    M_densityWall             = dataFile( ( section + "/PhysicalWall/DensityWall"                     ).data(), 1. );
    M_inertialModulus         = dataFile( ( section + "/PhysicalWall/coeffA"                          ).data(), 0. );

    M_longitudinalWall        = dataFile( ( section + "/PhysicalWall/LongitudinalWall"                ).data(), false );

    // Miscellaneous
    M_postprocessingDirectory = dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile      = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_CFLmax                  = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3.)/3. );

    if ( M_CFLmax > std::sqrt(3.)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Jacobian perturbation
    M_jacobianPerturbationArea     = dataFile( ( section + "/JacobianPerturbation/deltaArea"         ).data(), 0.001 );
    M_jacobianPerturbationFlowRate = dataFile( ( section + "/JacobianPerturbation/deltaFlowRate"     ).data(), 0.001 );
    M_jacobianPerturbationStress   = dataFile( ( section + "/JacobianPerturbation/deltaStress"       ).data(), 1.0 );

    // Physical Parameters
    M_computeCoefficients    = dataFile( ( section + "/PhysicalParameters/ComputeCoefficients"       ).data(), false );
    M_powerLawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_density                = dataFile( ( section + "/PhysicalParameters/density"                   ).data(), 1. );
    M_viscosity              = dataFile( ( section + "/PhysicalParameters/viscosity"                 ).data(), 0.035 );

    M_thickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_young                  = dataFile( ( section + "/PhysicalParameters/young"                     ).data(), 4.0E6 );
    M_poisson                = dataFile( ( section + "/PhysicalParameters/poisson"                   ).data(), 0.5 );

    M_externalPressure       = dataFile( ( section + "/PhysicalParameters/externalPressure"          ).data(), 0. );
    M_venousPressure         = dataFile( ( section + "/PhysicalParameters/venousPressure"            ).data(), 0. );
    M_friction               = dataFile( ( section + "/PhysicalParameters/Kr"                        ).data(), 1. );

    M_robertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_computeCoefficients    = dataFile( ( section + "/PhysicalParameters/ComputeCoefficients"       ).data(), false );

    // Reset all the containers
    resetContainers();

    // Select a law for the coefficients
    std::map< string, OneD_distributionLaw > distributionLawMap;
    distributionLawMap["uniform"]   = uniform;
    distributionLawMap["linear"]    = linear;
    distributionLawMap["pointwise"] = pointwise;

    OneD_distributionLaw distributionLaw = distributionLawMap[ dataFile( ( section + "/PhysicalParameters/DistributionLaw" ).data(), "uniform" ) ];
    switch ( distributionLaw )
    {
    case uniform:

        for ( UInt i = 0; i < M_meshPtr->numPoints() ; ++i )
        {
            // Physical Parameters
            M_thickness[i]                 = dataFile( ( section + "/PhysicalParameters/thickness"       ).data(), 0. );
            M_viscoelasticCoefficient[i]   = dataFile( ( section + "/PhysicalWall/ViscoelasticCoefficient" ).data(), 0. );

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
        linearInterpolation( M_viscoelasticCoefficient, dataFile, section + "/PhysicalWall/ViscoelasticCoefficient", 0. );

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

        for ( UInt i = 0; i < M_meshPtr->numPoints() ; ++i )
        {
            // Physical Parameters
            M_thickness[i]                 = dataFile( ( section + "/PhysicalParameters/thickness"       ).data(), 0., i );
            M_viscoelasticCoefficient[i]   = dataFile( ( section + "/PhysicalWall/ViscoelasticCoefficient" ).data(), 0., i );

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

        std::cout << "Warning: distributionLaw \"" << distributionLaw << "\"not available!" << std::endl;

        break;
    }

    updateCoefficients();
}

void
OneDFSIData::oldStyleSetup( const GetPot& dataFile, const std::string& section )
{
    std::cerr << std::endl << "Warning: OneDFSIData::oldStyleSetup( ... ) is deprecated!" << std::endl
              << "         You should use OneDFSIData::setup( ... ) instead!" << std::endl;

    // Model Type
    M_physicsType = OneDFSI::physicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_1DLinearPhysics" ) ];
    M_fluxType    = OneDFSI::fluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_sourceType  = OneDFSI::sourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_timeDataPtr.get() )
        M_timeDataPtr.reset( new time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Mesh setup - Space Discretization
    Real length = dataFile( ( section + "/discretization/x_right" ).data(), 1. ) -
                  dataFile( ( section + "/discretization/x_left"  ).data(), 0. );

    Real numberOfElements = dataFile( ( section + "/discretization/nb_elem" ).data(), 10 );
    regularMesh1D( *M_meshPtr, 1, numberOfElements, false, length, 0 );

    //std::cout << " 1D- Mesh nodes:                               " << M_meshPtr->numPoints() << std::endl;
    //std::cout << " 1D- Mesh elements:                            " << M_meshPtr->numElements() << std::endl;

    // Physical Wall Model
    M_inertialWall            = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelasticWall        = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_longitudinalWall        = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );

    // Miscellaneous
    M_postprocessingDirectory = dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_postprocessingFile      = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_inertialWall            = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelasticWall        = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_CFLmax                  = dataFile( ( section + "/miscellaneous/CFLmax"                         ).data(), std::sqrt(3.)/3. );

    if ( M_CFLmax > std::sqrt(3.)/3. )
        std::cout << "!!! WARNING: CFLmax greater than the theoretical value (see MOX21, eq. 1.47) - CONVERGENCE NOT GUARANTEED  !!!" << std::endl;

    // Jacobian perturbation
    M_jacobianPerturbationArea     = dataFile( ( section + "JacobianPerturbation/deltaArea"          ).data(), 0.001 );
    M_jacobianPerturbationFlowRate = dataFile( ( section + "JacobianPerturbation/deltaFlowRate"      ).data(), 0.001 );
    M_jacobianPerturbationStress   = dataFile( ( section + "JacobianPerturbation/deltaStress"        ).data(), 1.0 );

    // Physical Parameters
    M_computeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );
    M_powerLawCoefficient    = dataFile( ( section + "/PhysicalParameters/PowerlawCoefficient"       ).data(), 2 );

    M_density                = dataFile( ( section + "/1d_physics/density"                           ).data(), 1. );
    M_viscosity              = dataFile( ( section + "/1d_physics/viscosity"                         ).data(), 0.035 );

    M_densityWall            = dataFile( ( section + "/PhysicalParameters/densityWall"               ).data(), 1. );
    M_thickVessel            = dataFile( ( section + "/PhysicalParameters/thickVessel"               ).data(), false );
    M_young                  = dataFile( ( section + "/1d_physics/young"                             ).data(), 4.0E6 );
    M_poisson                = dataFile( ( section + "/1d_physics/ksi"                               ).data(), 0.5 );

    M_externalPressure       = dataFile( ( section + "/PhysicalParameters/externalPressure"          ).data(), 0. );
    M_venousPressure         = dataFile( ( section + "/PhysicalParameters/venousPressure"            ).data(), 0. );

    M_inertialModulus        = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_robertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

    M_friction               = dataFile( ( section + "/parameters/Kr"                                ).data(), 1. );
    M_computeCoefficients    = dataFile( ( section + "/parameters/use_physical_values"               ).data(), false );

    // Reset all the containers
    resetContainers();

    // Linear Parameters
    for ( UInt i = 0; i < M_meshPtr->numPoints() ; ++i )
    {
        M_viscoelasticCoefficient[i]   = dataFile( ( section + "/PhysicalParameters/gamma" ).data(), 0. );

        // Physical Parameters
        if ( M_computeCoefficients )
            M_area0[i]                 = OneDFSI::pow20( dataFile( ( section + "/1d_physics/radius"        ).data(), 0.5 ), 2 ) * M_PI;
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
OneDFSIData::updateCoefficients()
{
    if ( M_computeCoefficients )
    {
        // PowerlawProfile: s(r) = (1+2/gamma)*(1-r^gamma)
        Real radius (1.); //std::sqrt( M_Area0[i] / M_PI );

        Real profileIntegral =   ( 1+2./M_powerLawCoefficient ) * ( 1+2./M_powerLawCoefficient ) *
                                 (     radius * radius / 2 +
                                       std::pow(radius, 2*M_powerLawCoefficient+2) / (2*M_powerLawCoefficient+2 ) -
                                     2*std::pow(radius,   M_powerLawCoefficient+2) / (  M_powerLawCoefficient+2 )
                                 );

        // Compute Friction Coefficient: Kr = -2*pi*mu/rho*s'(R)
        M_friction = 2 * M_PI * M_viscosity / M_density * ( M_powerLawCoefficient + 2 ) * std::pow( radius, M_powerLawCoefficient - 1 );

        for ( UInt i = 0; i < M_meshPtr->numPoints(); ++i )
        {
            // Compute Coriolis Coefficient: Alpha = 2*pi/Area0*Int(s(r)^2)
            M_alpha[i] = 2 / ( radius * radius ) * profileIntegral;

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

            // Compute Viscoelastic Coefficient: gamma = h * E / ( 1 - nu^2 ) * T * tan(phi) / ( 4 * \sqrt(pi))
            M_viscoelasticCoefficient[i] = M_thickness[i] * M_young / (1 - M_poisson * M_poisson) * M_viscoelasticPeriod *
                                           std::tan( M_viscoelasticAngle ) / ( 4 * std::sqrt(M_PI) );
        }
    }

    // Compute the derivatives of the coefficients
    computeSpatialDerivatives();
}

//void
//OneDFSIData::initializeLinearParameters()
//{
//    for ( UInt indz=0; indz < M_meshPtr->numPoints(); ++indz )
//    {
//        M_celerity1[indz] = std::sqrt( M_beta0[indz] * M_beta1[indz] / M_density );
//        M_flux21[indz]    =  M_celerity1[indz] *  M_celerity1[indz];
//        M_source22[indz]  = M_friction / M_area0(indz);
//    }
//
//    M_celerity2 = - M_celerity1;
//
//    M_celerity1LeftEigenvector1 = M_celerity1;
//    M_celerity1LeftEigenvector2 = scalarVector_Type(M_meshPtr->numPoints(), 1.);
//    M_celerity2LeftEigenvector1 = M_celerity2;
//    M_celerity2LeftEigenvector2 = scalarVector_Type(M_meshPtr->numPoints(), 1.);
//
//    M_flux11 = ZeroVector(M_meshPtr->numPoints());
//    M_flux12 = scalarVector_Type(M_meshPtr->numPoints(), 1.);
//    M_flux22 = ZeroVector(M_meshPtr->numPoints());
//
//    M_source10 = ZeroVector(M_meshPtr->numPoints());
//    M_source11 = ZeroVector(M_meshPtr->numPoints());
//    M_source12 = ZeroVector(M_meshPtr->numPoints());
//
//    M_source20 = ZeroVector(M_meshPtr->numPoints());
//    M_source21 = ZeroVector(M_meshPtr->numPoints());
//}

// TODO Many changes to be done in these two methods:
// 1) variables name is very bad
// 2) do loop with explicit i++ is bad, use for loop instead
// 3) rename indices following other classes (iNode, iElement, etc..)
// 4) probably 3/4 of the code can be shared between the left and right method instead of duplications
//void
//OneDFSIData::stiffenVesselLeft( const Real& xl,        const Real& xr,
//                                       const Real& factor,    const Real& alpha,
//                                       const Real& delta,     const Real& n,
//                                       const Real& minDeltaX, const UInt& yesAdaptive )
//{
//    if ( yesAdaptive )
//    {
//        UInt iz=0, alpha_iz;
//        Real ratio, n_elem_delta,n_elem_r,n_elem_l;
//
//        alpha_iz = static_cast<Int>( std::floor( (alpha - delta / 2 ) / minDeltaX + 0.5 ) ) + ( ( numberOfElements() - 1 ) -
//                   static_cast<Int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) ) -
//                   static_cast<Int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) ) / 2;
//
//        n_elem_r = ( ( numberOfElements() - 1 ) - alpha_iz ) -
//                   static_cast<Int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) );
//
//        n_elem_l = alpha_iz -
//                   static_cast<Int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) );
//
//        n_elem_delta = static_cast<Real>( numberOfElements() - 1 ) / ( xr - xl ) * delta;
//
//        Real x_current,deltax,deltax_adaptive,deltax_uniform;
//        x_current = alpha;
//        do
//        {
//            ratio = ( alpha + delta / 2 - x_current ) / delta;
//
//            M_dBeta0dz[alpha_iz + iz] = M_beta0[alpha_iz + iz] * ( factor * (- n / delta) * ( std::pow(2,(n-1)) * std::pow(ratio, (n-1)) ) );
//            M_dBeta0dz[alpha_iz - iz] = M_dBeta0dz[alpha_iz + iz];
//
//            M_beta0[alpha_iz + iz] = M_beta0[alpha_iz + iz] * ( 1 + factor * ( std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//            M_beta0[alpha_iz - iz] = M_beta0[alpha_iz + iz] / ( 1 + factor * ( std::pow(2,(n-1)) * std::pow(ratio,n) ) ) * ( 1 + factor * ( 1 - ( std::pow(2,(n-1)) * std::pow(ratio,n) ) ) );
//
//            deltax_adaptive = ( -1 / n_elem_delta ) * ( 1 / ( ( -n / delta) * std::pow(2,(n-1) ) * std::pow( ratio , (n-1) ) ) );
//            deltax_uniform = ( ( alpha + delta / 2) - x_current ) / ( n_elem_r - iz );
//
//            ++iz;
//            deltax = ( ( deltax_adaptive < deltax_uniform ) && ( iz < n_elem_r ) ) ? deltax_adaptive : deltax_uniform;
//            ASSERT_PRE( deltax > 0 , "The left point is on the right..." );
//            x_current += deltax;
//        }
//        while ( ( x_current < ( alpha + delta/2 ) ) && ( ( alpha_iz - ( iz - 1 ) ) > 0) );
//
//        if ( ( alpha_iz - ( iz - 1 ) ) > 0 )
//        {
//            do
//            {
//                M_beta0[alpha_iz - iz] = M_beta0[alpha_iz - iz] * ( 1 + factor );
//                ++iz;
//            }
//            while ( (alpha_iz - ( iz - 1 ) ) > 0 );
//        }
//        else
//            std::cout << "[stiffenVesselRight] error! out of left boundary" << std::endl;
//    }
//
//    else
//    {
//        UInt iz=0;
//        Real ratio, x_current=xl, deltax;
//
//        deltax = ( xr - xl ) / static_cast<Real>(numberOfElements() - 1);
//
//        while ( ( x_current < ( alpha - delta / 2 ) ) && ( iz < numberOfElements() ) )
//        {
//            M_beta0[iz] =  M_beta0[iz] * ( 1 + factor );
//
//            ++iz;
//            x_current+=deltax;
//        }
//
//        while ( (x_current < alpha) && (iz < numberOfElements()) )
//        {
//            ratio = ( x_current - alpha - delta / 2 ) / delta;
//
//            M_dBeta0dz[iz] = M_beta0[iz] * ( factor * (- n / delta ) * ( std::pow(2,(n-1)) * std::pow(ratio,(n-1) ) ) );
//            M_beta0[iz]    = M_beta0[iz] * ( 1 + factor * ( 1 - std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//
//            ++iz;
//            x_current+=deltax;
//        }
//
//        while ( ( x_current < ( alpha + delta / 2 ) ) && (iz < numberOfElements()) )
//        {
//            ratio = ( alpha + delta / 2 - x_current ) / delta;
//
//            M_dBeta0dz[iz] = M_beta0[iz] * ( factor * ( -n / delta) * ( std::pow(2,(n-1)) * std::pow(ratio,(n-1) ) ) );
//            M_beta0[iz]    = M_beta0[iz] * ( 1 + factor * ( std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//
//            ++iz;
//            x_current += deltax;
//        }
//    }
//}

//void
//OneDFSIData::stiffenVesselRight( const Real& xl,        const Real& xr,
//                                        const Real& factor,    const Real& alpha,
//                                        const Real& delta,     const Real& n,
//                                        const Real& minDeltaX, const UInt& yesAdaptive )
//{
//    if ( yesAdaptive )
//    {
//        Real ratio, n_elem_delta,n_elem_r,n_elem_l;
//
//        UInt iz=0, alpha_iz;
//
//        alpha_iz = static_cast<Int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) + ( (numberOfElements() - 1) -
//                   static_cast<Int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) ) -
//                   static_cast<Int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) ) ) / 2;
//
//        n_elem_delta = static_cast<Real>( numberOfElements() - 1 ) / ( xr - xl ) * delta;
//
//        n_elem_r = ( ( numberOfElements() - 1 ) - alpha_iz ) -
//                   static_cast<Int>( std::floor( ( xr - ( alpha + delta / 2 ) ) / minDeltaX + 0.5 ) );
//
//        n_elem_l = alpha_iz -
//                   static_cast<Int>( std::floor( ( alpha - delta / 2 ) / minDeltaX + 0.5 ) );
//
//        Real x_current,deltax,deltax_adaptive,deltax_uniform;
//        x_current = alpha;
//        do
//        {
//            ratio = ( alpha + delta / 2 - x_current ) / delta;
//
//            M_dBeta0dz[alpha_iz + iz] = M_beta0[alpha_iz + iz] * ( factor * ( n / delta) * ( std::pow(2,(n-1)) * std::pow(ratio,(n-1)) ) );
//            M_dBeta0dz[alpha_iz - iz] = M_dBeta0dz[alpha_iz + iz];
//
//            M_beta0[alpha_iz + iz] = M_beta0[alpha_iz + iz] * ( 1 + factor * ( 1 - ( std::pow(2,(n-1)) * std::pow(ratio,n) ) ) );
//            M_beta0[alpha_iz - iz] = M_beta0[alpha_iz + iz] / ( 1 + factor * ( 1 - ( std::pow(2,(n-1)) * std::pow(ratio,n) ) ) ) * ( 1 + factor * ( std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//
//            deltax_adaptive = ( -1 / n_elem_delta ) * ( 1 / ( ( -n / delta) * std::pow(2,(n-1) ) * std::pow( ratio , (n-1) ) ) );
//            deltax_uniform = ( ( alpha + delta / 2) - x_current ) / ( n_elem_r - iz );
//
//            ++iz;
//            deltax = ( ( deltax_adaptive < deltax_uniform ) && ( iz < n_elem_r) ) ? deltax_adaptive : deltax_uniform;
//            ASSERT_PRE( deltax > 0 , "The left point is on the right..." );
//            x_current += deltax;
//
//        }
//        while ( x_current < ( alpha + delta / 2 ) && ( ( alpha_iz - ( iz - 1 ) ) > 0) );
//
//        if ( ( alpha_iz + iz ) <= ( numberOfElements() -1 ) )
//        {
//            do
//            {
//                M_beta0[alpha_iz + iz] = M_beta0[alpha_iz + iz] * ( 1 + factor );
//                ++iz;
//            }
//            while ( ( alpha_iz + iz - 1 ) < ( numberOfElements() -1 ) );
//
//        }
//        else
//            std::cout << "\n[stiffenVesselRight] error! out of right boundary" << std::endl;
//    }
//    else
//    {
//        UInt iz = numberOfElements()-1;
//        Real ratio, x_current=xr, deltax;
//
//        deltax = ( xr - xl ) / static_cast<Real>(numberOfElements() - 1 );
//
//        while ( ( x_current > ( alpha + delta / 2 ) ) && ( ( iz + 1 ) > 0 ) )
//        {
//            M_beta0[iz] =  M_beta0[iz] * ( 1 + factor );
//
//            --iz;
//            x_current -= deltax;
//        }
//
//        while ( ( x_current > alpha ) && ( ( iz + 1 ) > 0 ) )
//        {
//            ratio=( ( ( alpha + delta / 2 ) - x_current ) / delta );
//
//            M_dBeta0dz[iz] = M_beta0[iz] * ( factor * ( n / delta) *  ( std::pow(2,(n-1)) * std::pow(ratio,(n-1)) ) );
//            M_beta0[iz]    = M_beta0[iz] * ( 1 + factor * ( 1 - std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//
//            --iz;
//            x_current -= deltax;
//        }
//
//        while ( ( x_current > ( alpha - delta / 2 ) ) && ( ( iz + 1 ) > 0 ) )
//        {
//            ratio = ( x_current - alpha - delta / 2 ) / delta;
//
//            M_dBeta0dz[iz] = M_beta0[iz] * ( factor * ( n / delta) * ( std::pow(2,(n-1)) * std::pow(ratio,(n-1)) ) );
//            M_beta0[iz]    = M_beta0[iz] * ( 1 + factor * ( std::pow(2,(n-1)) * std::pow(ratio,n) ) );
//
//            --iz;
//            x_current -= deltax;
//        }
//    }
//}

void
OneDFSIData::showMe( std::ostream& output ) const
{
    //output << std::scientific << std::setprecision(15);

    // Model
    output << "\n*** Values for data [Model]" << std::endl << std::endl;
    output << "Physics Type           = " << enum2String( M_physicsType, OneDFSI::physicsMap ) << std::endl;
    output << "Flux Type              = " << enum2String( M_fluxType,    OneDFSI::fluxMap    ) << std::endl;
    output << "Source Type            = " << enum2String( M_sourceType,  OneDFSI::sourceMap  ) << std::endl;

    // Time
    output << "\n*** Values for data [time_discretization]" << std::endl << std::endl;
    M_timeDataPtr->showMe( output );

    // Space Discretization
    output << "\n*** Values for data [space_discretization]" << std::endl << std::endl;
    output << "Length                 = " << length() << std::endl;
    output << "NumberOfElements       = " << M_meshPtr->numElements() << std::endl;

    // Physical Wall Model
    output << "\n*** Values for data [PhysicalWallModel]" << std::endl << std::endl;
    output << "Viscoelastic wall      = " << M_viscoelasticWall << std::endl;
    output << "Viscoelastic angle     = " << M_viscoelasticAngle << std::endl;
    output << "Viscoelastic period    = " << M_viscoelasticPeriod << std::endl;
    output << "Viscoelastic coeff.    = " << M_viscoelasticCoefficient << std::endl;
    output << "Inertial wall          = " << M_inertialWall << std::endl;
    output << "Wall density           = " << M_densityWall << std::endl;
    output << "Inertial modulus       = " << M_inertialModulus << std::endl;
    output << "Longitudinal wall      = " << M_longitudinalWall << std::endl;

    // Miscellaneous
    output << "\n*** Values for data [miscellaneous]" << std::endl << std::endl;
    output << "Postprocessing Dir.    = " << M_postprocessingDirectory << std::endl;
    output << "Postprocessing File    = " << M_postprocessingFile << std::endl;
    output << "Maximum admissible CFL = " << M_CFLmax << std::endl;

    // Jacobian perturbation
    output << "Jacobian perturbation Area      = " << M_jacobianPerturbationArea << std::endl;
    output << "Jacobian perturbation Flow Rate = " << M_jacobianPerturbationFlowRate << std::endl;
    output << "Jacobian perturbation Stress    = " << M_jacobianPerturbationStress << std::endl;

    // Physical Parameters
    output << "\n*** Values for data [PhysicalParameters]" << std::endl << std::endl;
    output << "Compute Coefficients   = " << M_computeCoefficients << std::endl;
    output << "Powerlaw Coefficient   = " << M_powerLawCoefficient << std::endl;
    output << "Fluid density          = " << M_density << std::endl;
    output << "Fluid dyn. viscosity   = " << M_viscosity << std::endl;
    output << "Thick vessel           = " << M_thickVessel << std::endl;
    output << "Young modulus          = " << M_young << std::endl;
    output << "Poisson                = " << M_poisson << std::endl;
    output << "External pressure      = " << M_externalPressure << std::endl;
    output << "Venous pressure        = " << M_venousPressure << std::endl;
    output << "Robertson correction   = " << M_robertsonCorrection << std::endl;

    output << "Area0                  = " << M_area0 << std::endl;
    output << "dArea0dz               = " << M_dArea0dz << std::endl;
    output << "Beta0                  = " << M_beta0 << std::endl;
    output << "dBeta0dz               = " << M_dBeta0dz << std::endl;
    output << "Beta1                  = " << M_beta1 << std::endl;
    output << "dBeta1dz               = " << M_dBeta1dz << std::endl;

    output << "Alpha (Coriolis)       = " << M_alpha << std::endl;
    output << "dAlpha (Coriolis)      = " << M_dAlphadz << std::endl;

    output << "Thickness              = " << M_thickness << std::endl;
    output << "Friction               = " << M_friction << std::endl;

    // Linear Parameters
    output << "\n*** Values for data [LinearParameters]" << std::endl << std::endl;
    output << "Flux11                 = " << M_flux11 << std::endl;
    output << "Flux12                 = " << M_flux12 << std::endl;
    output << "Flux21                 = " << M_flux21 << std::endl;
    output << "Flux22                 = " << M_flux22 << std::endl;
    output << "Celerity1              = " << M_celerity1 << std::endl;
    output << "Celerity2              = " << M_celerity2 << std::endl;
    output << "Eigenvector11          = " << M_celerity1LeftEigenvector1 << std::endl;
    output << "Eigenvector12          = " << M_celerity1LeftEigenvector2 << std::endl;
    output << "Eigenvector21          = " << M_celerity2LeftEigenvector1 << std::endl;
    output << "Eigenvector22          = " << M_celerity2LeftEigenvector2 << std::endl;
    output << "Source10               = " << M_source10 << std::endl;
    output << "Source20               = " << M_source20 << std::endl;
    output << "Source11               = " << M_source11 << std::endl;
    output << "Source12               = " << M_source12 << std::endl;
    output << "Source21               = " << M_source21 << std::endl;
    output << "Source22               = " << M_source22 << std::endl;
}

// ===================================================
// Get Methods - Physical Parameters
// ===================================================
const Real&
OneDFSIData::robertsonCorrection() const
{
    if ( M_robertsonCorrection != 1. )
        std::cout << "!!! WARNING: Robertson corretion has not been checked in this version of the code !!!" << std::endl;

    return M_robertsonCorrection;
}
// ===================================================
// Private methods
// ===================================================
void
OneDFSIData::linearInterpolation( scalarVector_Type& vector,
                                  const GetPot& dataFile,
                                  const std::string& quantity,
                                  const Real& defaultValue,
                                  const bool& isArea )
{
    Real a  = dataFile( quantity.data(), defaultValue, 0 );
    Real b  = dataFile( quantity.data(), a, 1 );

    Real xa = M_meshPtr->firstPoint().x();
    Real xb = M_meshPtr->lastPoint().x();

    for ( UInt i(0) ; i < M_meshPtr->numPoints() ; ++i )
        if ( isArea )
        {
            vector[i] = std::sqrt(a / M_PI) + ( std::sqrt(b / M_PI) - std::sqrt(a / M_PI) ) /
                                              ( xb - xa ) * ( M_meshPtr->point( i ).x() - xa );
            vector[i] *= vector[i] * M_PI;
        }
        else
            vector[i] = a + (b - a) / ( xb - xa ) * ( M_meshPtr->point( i ).x() - xa );

    // linearInterpolation to disable tapering (when needed)
//    for ( UInt i(0); i < M_meshPtr->numPoints() ; ++i )
//        if ( isArea )
//            vector[i] = (a + b + 2 * std::sqrt(a*b)) / 4;
//        else
//            vector[i] = (a + b) / 2;
}

void
OneDFSIData::computeSpatialDerivatives()
{
    for ( UInt iNode( 0 ) ; iNode < M_meshPtr->numPoints() ; ++iNode )
    {
        M_dArea0dz[iNode] = computeSpatialDerivativeAtNode( M_area0, iNode );
        M_dBeta0dz[iNode] = computeSpatialDerivativeAtNode( M_beta0, iNode );
        M_dBeta1dz[iNode] = computeSpatialDerivativeAtNode( M_beta1, iNode );
        M_dAlphadz[iNode] = computeSpatialDerivativeAtNode( M_alpha, iNode );
    }
}

void
OneDFSIData::resetContainers()
{
    M_viscoelasticCoefficient.resize( M_meshPtr->numPoints() );
    M_thickness.resize( M_meshPtr->numPoints() );

    M_area0.resize( M_meshPtr->numPoints() );
    M_beta0.resize( M_meshPtr->numPoints() );
    M_beta1.resize( M_meshPtr->numPoints() );
    M_alpha.resize( M_meshPtr->numPoints() );

    M_dArea0dz.resize( M_meshPtr->numPoints() );
    M_dBeta0dz.resize( M_meshPtr->numPoints() );
    M_dBeta1dz.resize( M_meshPtr->numPoints() );
    M_dAlphadz.resize( M_meshPtr->numPoints() );

    // Linear Parameters defined along the 1D model
    M_flux11.resize( M_meshPtr->numPoints() );
    M_flux12.resize( M_meshPtr->numPoints() );
    M_flux21.resize( M_meshPtr->numPoints() );
    M_flux22.resize( M_meshPtr->numPoints() );
    M_celerity1.resize( M_meshPtr->numPoints() );
    M_celerity2.resize( M_meshPtr->numPoints() );
    M_celerity1LeftEigenvector1.resize( M_meshPtr->numPoints() );
    M_celerity1LeftEigenvector2.resize( M_meshPtr->numPoints() );
    M_celerity2LeftEigenvector1.resize( M_meshPtr->numPoints() );
    M_celerity2LeftEigenvector2.resize( M_meshPtr->numPoints() );
    M_source10.resize( M_meshPtr->numPoints() );
    M_source20.resize( M_meshPtr->numPoints() );
    M_source11.resize( M_meshPtr->numPoints() );
    M_source12.resize( M_meshPtr->numPoints() );
    M_source21.resize( M_meshPtr->numPoints() );
    M_source22.resize( M_meshPtr->numPoints() );
}

} // LifeV namespace
