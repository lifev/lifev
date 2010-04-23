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
    M_Mesh                      (),
    M_x_left                    (),
    M_x_right                   (),
    M_Area0                     (),
    M_dArea0dz                  (),
    M_PressBeta0                (),
    M_dPressBeta0dz             (),
    M_PressBeta1                (),
    M_dPressBeta1dz             (),
    M_AlphaCoriolis             (),
    M_dAlphaCoriolisdz          (),
    M_FrictionKr                (),
    M_DensityRho                (),
    M_DensityWall               (),
    M_Thickness                 (),
    M_Gamma                     (),
    M_CoeffA                    (),
    M_RobertsonCorrection       (),
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
    M_Source22                  (),
    M_post_dir                  (),
    M_post_file                 (),
    M_verbose                   (),
    M_CFL                       (),
    M_UW                        (),
    M_inertial_wall             (),
    M_viscoelastic_wall         (),
    M_linearize_string_model    (),
    M_linearize_equations       (),
    M_longitudinal_wall         (),
    M_flux_second_der           (),
    M_dP_dt_steps               (),
    M_firstNode                 (),
    M_lastNode                  (),
    M_initVar                   (),
    M_restValue                 (),
    M_initValue                 (),
    M_multiplier                (),
    M_width                     ()
{
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_Data::setup( const GetPot& dataFile, const std::string& section )
{
    // Model Type
    M_PhysicsType = OneDimensionalModel_PhysicsMap[ dataFile( ( section + "/Model/PhysicsType" ).data(), "OneD_LinearPhysics" ) ];
    M_FluxType    = OneDimensionalModel_FluxMap[    dataFile( ( section + "/Model/FluxType"    ).data(), "OneD_1DLinearFlux" ) ];
    M_SourceType  = OneDimensionalModel_SourceMap[  dataFile( ( section + "/Model/SourceType"  ).data(), "OneD_1DLinearSource" ) ];

    // If data time has not been set
    if ( !M_Time.get() )
        M_Time.reset( new Time_Type( dataFile, (section + "/time_discretization" ).data() ) );

    // Space Discretization
    M_x_left                 = dataFile( ( section + "/space_discretization/x_left"                  ).data(), 0. );
    M_x_right                = dataFile( ( section + "/space_discretization/x_right"                 ).data(), 1. );

    // Mesh setup
    std::cout << " 1D- Mesh setup ...                            " << std::flush;
    Chrono chrono;
    chrono.start();

    M_Mesh.reset( new Mesh_Type );
    M_Mesh->setUp( M_x_left, M_x_right, dataFile( ( section + "/space_discretization/NumberOfElements" ).data(), 10 ) );

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s."<< std::endl;

    std::cout << " 1D- Mesh nodes:                               " << M_Mesh->numPoints() << std::endl;
    std::cout << " 1D- Mesh elements:                            " << M_Mesh->numElements() << std::endl;

    // Physical Parameters
    M_DensityRho             = dataFile( ( section + "/PhysicalParameters/rho"                       ).data(), 1. );
    M_DensityWall            = dataFile( ( section + "/PhysicalParameters/rho_w"                     ).data(), 1. );
    M_Thickness              = dataFile( ( section + "/PhysicalParameters/thickness"                 ).data(), 0. );
    M_Gamma                  = dataFile( ( section + "/PhysicalParameters/gamma"                     ).data(), 0. );
    M_CoeffA                 = dataFile( ( section + "/PhysicalParameters/coeffA"                    ).data(), 0. );
    M_RobertsonCorrection    = dataFile( ( section + "/PhysicalParameters/RobertsonCorrection"       ).data(), 1. );

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

    // Miscellaneous
    M_post_dir               = dataFile( ( section + "/miscellaneous/post_dir"                       ).data(), "./" );
    M_post_file              = dataFile( ( section + "/miscellaneous/post_file"                      ).data(), "sol" );
    M_verbose                = dataFile( ( section + "/miscellaneous/verbose"                        ).data(), 1 );
    M_CFL                    = dataFile( ( section + "/miscellaneous/showCFL"                        ).data(), 0 );
    M_UW                     = dataFile( ( section + "/miscellaneous/alternate_solver"               ).data(), false );
    M_inertial_wall          = dataFile( ( section + "/miscellaneous/inertial_wall"                  ).data(), false );
    M_viscoelastic_wall      = dataFile( ( section + "/miscellaneous/viscoelastic_wall"              ).data(), false );
    M_linearize_string_model = dataFile( ( section + "/miscellaneous/linearize_string_model"         ).data(), true );
    M_linearize_equations    = dataFile( ( section + "/miscellaneous/linearize_equations"            ).data(), false );
    M_longitudinal_wall      = dataFile( ( section + "/miscellaneous/longitudinal_wall"              ).data(), false );
    M_flux_second_der        = dataFile( ( section + "/miscellaneous/compute_flux_second_derivative" ).data(), false );
    M_dP_dt_steps            = dataFile( ( section + "/miscellaneous/pressure_derivative_steps"      ).data(), 1 );

    // Initialize
    M_firstNode              = dataFile( ( section + "/initialize/firstnode"                         ).data(), 1 );
    M_lastNode               = dataFile( ( section + "/initialize/lastnode"                          ).data(), 2 );
    M_initVar                = dataFile( ( section + "/initialize/var"                               ).data(), "P" );
    M_initValue              = dataFile( ( section + "/initialize/init_value"                        ).data(), 0. );
    M_restValue              = dataFile( ( section + "/initialize/rest_value"                        ).data(), 0. );
    M_multiplier             = dataFile( ( section + "/initialize/multiplier"                        ).data(), 1. );
    M_width                  = dataFile( ( section + "/initialize/width"                             ).data(), 5. );

    if( dataFile( ( section + "/PhysicalParameters/use_physical_values" ).data(), false ) ) // CHECK THIS!!!
    {
        Debug( 6320 ) << "[OneDimensionalModel_Data_Operator] initializing from physical values\n";
        initParam( dataFile );
    }
}

void
OneDimensionalModel_Data::initParam( const GetPot& dataFile )  // CHECK THIS!!!
{
    Real Young_modulus    = dataFile("1d_physics/young"          , 5.e6);
    Real thickness        = dataFile("1d_physics/thickness"      , 0.05);
    Real reference_radius = dataFile("1d_physics/radius"         , 1.);
    Real viscosity        = dataFile("1d_physics/viscosity"      , 0.035);  //???
    Real ksi              = dataFile("1d_physics/ksi"            , 0.5);  //???
    Real friction_factor  = dataFile("1d_physics/friction_factor", 8.);  //???
    bool thick_vessel     = dataFile("1d_physics/thick_vessel"   , 0);

    Real _A0( M_PI*reference_radius*reference_radius );

    Real _beta0, _beta1;
    if( thick_vessel ){ // see Amadori, Ferrari, Formaggia (MOX report 86)
        //! beta0
        _beta0 = - thickness*Young_modulus*sqrt(M_PI) /
            ( sqrt(_A0) *
              ( (1 - ksi * ksi)
                + ksi * (1 + ksi) * (thickness * sqrt(M_PI)
                                     / sqrt(_A0) )
                )
              );
        //! beta1
        _beta1 = - 0.5;
    }
    else{
        //! beta0
        _beta0 = thickness * Young_modulus * sqrt(M_PI) /
            ( sqrt(_A0) * (1 - ksi * ksi) );
        //! beta1
        _beta1 = 0.5;
    }

    Real _kr( friction_factor * M_PI * viscosity );

    //-------------------------------------------
    //! Initialisation of the parameter variables
    //-------------------------------------------
    //! A0
    M_Area0.resize( M_Mesh->numPoints() );
    //! beta0
    M_PressBeta0.resize( M_Mesh->numPoints() );
    //! beta1
    M_PressBeta1.resize( M_Mesh->numPoints() );
    //! Kr
    M_FrictionKr.resize( M_Mesh->numPoints() );

    for ( UInt i = 0; i < M_Mesh->numPoints(); ++i )
    {
        M_Area0[i]         = _A0;
        M_PressBeta0[i]    = _beta0;
        M_PressBeta1[i]    = _beta1;
        M_FrictionKr[i]    = _kr;
    }

    M_Thickness  = thickness;
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
        M_Celerity1[indz] = std::sqrt( M_PressBeta0[indz] * M_PressBeta1[indz] / M_DensityRho );
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
    output << "x_left                 = " << M_x_left << std::endl;
    output << "x_right                = " << M_x_right << std::endl;
    output << "nb_elem                = " << M_Mesh->numPoints() << std::endl;

    // Physical Parameters
    output << "\n*** Values for data [PhysicalParameters]\n\n";
    output << "Area0                  = " << M_Area0 << "\n";
    output << "dArea0dz               = " << M_dArea0dz << "\n";
    output << "Beta0                  = " << M_PressBeta0 << "\n";
    output << "dBeta0dz               = " << M_dPressBeta0dz << "\n";
    output << "Beta1                  = " << M_PressBeta1 << "\n";
    output << "dBeta1dz               = " << M_dPressBeta1dz << "\n";

    output << "Alpha (Coriolis)       = " << M_AlphaCoriolis << "\n";
    output << "dAlpha (Coriolis)      = " << M_dAlphaCoriolisdz << "\n";

    output << "Fluid density          = " << M_DensityRho << "\n";
    output << "Friction               = " << M_FrictionKr << "\n";
    output << "Wall density           = " << M_DensityWall << "\n";
    output << "Viscoelastic modulus   = " << M_Gamma << "\n";
    output << "Inertial modulus       = " << M_CoeffA << "\n";
    output << "Thickness              = " << M_Thickness << "\n";

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

    // Miscellaneous
    output << "\n*** Values for data [miscellaneous]\n\n";
    output << "post_dir               = " << M_post_dir << std::endl;
    output << "post_file              = " << M_post_file << std::endl;
    output << "verbose                = " << M_verbose << std::endl;
    output << "CFL                    = " << M_CFL << std::endl;
    output << "UW                     = " << M_UW << std::endl;
    output << "Use Inertial Wall      = " << M_inertial_wall << std::endl;
    output << "Use Viscoelastic Wall  = " << M_viscoelastic_wall << std::endl;
    output << "Linearize Model        = " << M_linearize_string_model << std::endl;
    output << "Linearize Equations    = " << M_linearize_equations << std::endl;
    output << "Longitudinal Wall      = " << M_longitudinal_wall << std::endl;
    output << "Flux Second Derivative = " << M_flux_second_der << std::endl;
    output << "Pressure Derivative    = " << M_dP_dt_steps << std::endl;

    // Initialize
    output << "\n*** Values for data [initialize]\n\n";
    output << "First Node             = " << M_firstNode << std::endl;
    output << "Last Node              = " << M_lastNode << std::endl;
    output << "Init Variable          = " << M_initVar << std::endl;
    output << "Init Value             = " << M_initValue << std::endl;
    output << "Rest Value             = " << M_restValue << std::endl;
    output << "Multiplier             = " << M_multiplier << std::endl;
    output << "Width                  = " << M_width << std::endl;
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_Data::setDataTime( const Time_ptrType DataTime )
{
    M_Time = DataTime;
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

OneDimensionalModel_Data::Time_ptrType
OneDimensionalModel_Data::dataTime( void ) const
{
    return M_Time;
}

OneDimensionalModel_Data::Mesh_ptrType
OneDimensionalModel_Data::mesh() const
{
    return M_Mesh;
}

const Real&
OneDimensionalModel_Data::xLeft() const
{
    return M_x_left;
}

const Real&
OneDimensionalModel_Data::xRight() const
{
    return M_x_right;
}

Real
OneDimensionalModel_Data::nbElem() const
{
    return M_Mesh->numElements();
}

const std::string&
OneDimensionalModel_Data::PostDirectory() const
{
    return M_post_dir;
}

const std::string&
OneDimensionalModel_Data::PostFile() const
{
    return M_post_file;
}

const int&
OneDimensionalModel_Data::verbose() const
{
    return M_verbose;
}

const Real&
OneDimensionalModel_Data::CFL() const
{
    return M_CFL;
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

const int&
OneDimensionalModel_Data::firstNode() const
{
    return M_firstNode;
}

const int&
OneDimensionalModel_Data::lastNode() const
{
    return M_lastNode;
}

const std::string&
OneDimensionalModel_Data::initVar() const
{
    return M_initVar;
}

const Real&
OneDimensionalModel_Data::restValue() const
{
    return M_restValue;
}

const Real&
OneDimensionalModel_Data::initValue() const
{
    return M_initValue;
}

const Real&
OneDimensionalModel_Data::multiplier() const
{
    return M_multiplier;
}

const Real&
OneDimensionalModel_Data::width() const
{
    return M_width;
}

// ===================================================
// Get Methods - Physical Parameters
// ===================================================
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

const Real&
OneDimensionalModel_Data::DensityRho() const
{
    return M_DensityRho;
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
OneDimensionalModel_Data::Gamma() const
{
    return M_Gamma;
}

const Real&
OneDimensionalModel_Data::CoeffA() const
{
    return M_CoeffA;
}

const Real&
OneDimensionalModel_Data::RobertsonCorrection() const
{
    return M_RobertsonCorrection;
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
