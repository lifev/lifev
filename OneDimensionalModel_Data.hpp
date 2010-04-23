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

#ifndef ONEDIMENSIONALMODEL_DATA_H
#define ONEDIMENSIONALMODEL_DATA_H

// LIFEV
#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/dataTime.hpp>
#include <lifemc/lifemesh/regionMesh1D.hpp>

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>

namespace LifeV {

//! OneDimensionalModel_Data - Class which read and holds all the data for the One Dimensional Model Solver.
/*!
 *  @author Vincent Martin, Cristiano Malossi
 *
 *  NOTE: Physical Parameters
 *  =========================
 *
 *  Parameters: Area0, alpha, beta0, beta1, Kr, rho.
 *
 *  Euler equations
 *  dA/dt + dQ/dz = 0
 *  dQ/dt + d/dz(alpha * Q^2/A) + A/rho * dP/dz + Kr * Q/A = 0
 *
 *  with
 *  P - P_ext = beta0 [ ( A / Area0 )^{beta1} - 1 ]
 *
 *  BEWARE: there are at least 2 or 3 different ways of defining it!!!
 *
 *  CONVENTIONS used here:
 *  Parameter homogeneous to a pressure:
 *  P - P_ext = PressBeta0 [ ( A / Area0 )^{PressBeta1} - 1 ]
 *
 *  This PressBeta0 is homogeneous to a pressure.
 *  In most cases PressBeta1 is taken equal to 1/2.
 *
 *  PressBeta0 = ( \sqrt{\pi} h_0 E ) / ( ( 1 - \ksi^2 ) * \sqrt{Area0} )
 *  OTHER CONVENTION not used here:
 *
 *  a) from Formaggia and Veneziani (p. 1.10, MOX report No 21 - june 2003)
 *  P - P_ext = \tilde{\beta_0} ( \sqrt{A} - \sqrt{A_0} ) / A_0
 *  with
 *  \beta0 = ( \sqrt{\pi} h_0 E ) / ( 1 - \ksi^2 )
 *
 *  link with PressBeta0: \tilde{\beta_0} = PressBeta0 * \sqrt{A_0}
 *
 *  b) Auxiliary Parameter often used in the model1D code (by J-FG or D Lamponi)
 *  (ONLY when beta1=1/2 !!)
 *  P - P_ext = 2 * rho * AuxBeta ( \sqrt{A} - \sqrt{A_0} )
 *
 *  link with PressBeta0:          AuxBeta = PressBeta0 * PressBeta1 / ( rho * Area0^(PressBeta1) )
 *  or whenever PressBeta1 = 1/2 : AuxBeta = PressBeta0 / ( 2 * rho * \sqrt{A_0} )
 *
 *
 *
 *  NOTE: Linear Parameters
 *  =======================
 *
 *  Parameters: F11, F12, F21, F22, celerity1, celerity2
 *
 *  Equations:
 *  dU1/dt + F11 dU1/dz + F12 dU2/dz = 0
 *  dU2/dt + F21 dU1/dz + F22 dU2/dz = 0
 *
 *  The flux matrix F = [F11, F12 ; F21 F22] has the eigenvalues
 *  celerity1, celerity2.
 */
class OneDimensionalModel_Data
{
public:

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            Time_ptrType;

    typedef RegionMesh1D< LinearLine >                                Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >                            Mesh_ptrType;

    //! @name Constructors & Destructor
    //@{

    OneDimensionalModel_Data();

    //! Destructor
    ~OneDimensionalModel_Data() {}

    //@}


    //! @name Methods
    //@{

    void setup( const GetPot& dataFile, const std::string& section = "1D_Model" );

    void initParam( const GetPot& dataFile ); // TO BE CHECKED - DON'T USE IT!

    void initLinearParam( const GetPot& dataFile ); // TO BE CHECKED - DON'T USE IT!

    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    void setDataTime( const Time_ptrType DataTime );

    void setArea0( const Real& Area0, const UInt& i );

    void setBeta0( const Real& Beta0, const UInt& i );

    void setdBeta0dz( const Real& dBeta0dz, const UInt& i );

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics type
    /*!
     * @return Physics Type
     */
    const OneDimensionalModel_PhysicsTypes& PhysicsType() const;

    //! Get the flux type
    /*!
     * @return Flux Type
     */
    const OneDimensionalModel_FluxTypes& FluxType() const;

    //! Get the source type
    /*!
     * @return Source Type
     */
    const OneDimensionalModel_SourceTypes& SourceType() const;

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    Time_ptrType       dataTime() const;

    Mesh_ptrType       mesh() const;

    const Real&        xLeft() const;
    const Real&        xRight() const;
          Real         nbElem() const;

    const std::string& PostDirectory() const;
    const std::string& PostFile() const;

    const int&         verbose() const;

    const Real&        CFL() const;
    const bool&        UW() const;
    const bool&        inertialWall() const;
    const bool&        viscoelasticWall() const;
    const bool&        linearizeStringModel() const;
    const bool&        linearizeEquations() const;
    const bool&        longitudinalWall() const;

    const bool&        fluxSecondDer() const;

    const int&         DPdtSteps() const;

    const int&         firstNode() const;
    const int&         lastNode() const;
    const std::string& initVar() const;
    const Real&        restValue() const;
    const Real&        initValue() const;
    const Real&        multiplier() const;
    const Real&        width() const;

    // Physical Parameters
    const Real& Area0( const UInt& i ) const;
    const Real& Beta0( const UInt& i ) const;
    const Real& Beta1( const UInt& i ) const;
    const Real& dArea0dz( const UInt& i ) const;
    const Real& dBeta0dz( const UInt& i ) const;
    const Real& dBeta1dz( const UInt& i ) const;

    const Real& AlphaCor( const UInt& i ) const;
    const Real& dAlphaCordz( const UInt& i ) const;

    const Real& FrictionKr( const UInt& i ) const;
    const Real& DensityRho() const;
    const Real& DensityWall() const;
    const Real& Gamma() const;
    const Real& CoeffA() const;
    const Real& RobertsonCorrection() const;
    const Real& Thickness() const;

    // Linear Parameters
    const Real& Flux11( const UInt& i ) const;
    const Real& Flux12( const UInt& i ) const;
    const Real& Flux21( const UInt& i ) const;
    const Real& Flux22( const UInt& i ) const;

    const Real& Celerity1( const UInt& i ) const;
    const Real& Celerity2( const UInt& i ) const;

    const Real& LeftEigenVector11( const UInt& i ) const;
    const Real& LeftEigenVector12( const UInt& i ) const;
    const Real& LeftEigenVector21( const UInt& i ) const;
    const Real& LeftEigenVector22( const UInt& i ) const;

    const Real& Source10( const UInt& i ) const;
    const Real& Source20( const UInt& i ) const;
    const Real& Source11( const UInt& i ) const;
    const Real& Source12( const UInt& i ) const;
    const Real& Source21( const UInt& i ) const;
    const Real& Source22( const UInt& i ) const;

    //@}

protected:

    //! Model
    OneDimensionalModel_PhysicsTypes M_PhysicsType;
    OneDimensionalModel_FluxTypes    M_FluxType;
    OneDimensionalModel_SourceTypes  M_SourceType;

    //! Data containers for time and mesh
    Time_ptrType      M_Time;
    Mesh_ptrType      M_Mesh;

    //! Space Discretization
    Real              M_x_left;  //! left coordinate
    Real              M_x_right; //! right coordinate

    //! Physical Parameters
    ScalVec M_Area0;
    ScalVec M_dArea0dz;
    //! P - P_ext = PressBeta0 [ ( A / Area0 )^{PressBeta1} - 1 ]
    ScalVec M_PressBeta0;    // homogeneous to a pressure
    ScalVec M_dPressBeta0dz; // homogeneous to a pressure
    ScalVec M_PressBeta1;    // power coeff (>0, often=1/2)
    ScalVec M_dPressBeta1dz; // power coeff (>0, often=1/2)
    ScalVec M_AlphaCoriolis; // Coriolis coefficient (often called alpha)
    ScalVec M_dAlphaCoriolisdz;
    //! Friction parameter Kr
    ScalVec M_FrictionKr;
    Real M_DensityRho;     // Density rho (always taken constant along the vessel)
    Real M_DensityWall;
    Real M_Thickness;
    Real M_Gamma;
    Real M_CoeffA;
    Real M_RobertsonCorrection;

    //! Flux matrix
    ScalVec M_Flux11;
    ScalVec M_Flux12;
    ScalVec M_Flux21;
    ScalVec M_Flux22;

    //! Celerities of the linear problem (eigenvalues of the flux matrix)
    ScalVec M_Celerity1;
    ScalVec M_Celerity2;

    //! Eigenvector for first eigenvalue
    ScalVec M_Celerity1LeftEigenvector1;
    ScalVec M_Celerity1LeftEigenvector2;
    //! Eigenvector for second eigenvalue
    ScalVec M_Celerity2LeftEigenvector1;
    ScalVec M_Celerity2LeftEigenvector2;

    //! Source matrix
    ScalVec M_Source10;
    ScalVec M_Source20;
    ScalVec M_Source11;
    ScalVec M_Source12;
    ScalVec M_Source21;
    ScalVec M_Source22;

    //! Miscellaneous
    std::string       M_post_dir; //! full directory name (including path)
    std::string       M_post_file; //! output file name
    int               M_verbose;
    Real              M_CFL;
    bool              M_UW;
    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool              M_inertial_wall;
    bool              M_viscoelastic_wall;
    bool              M_linearize_string_model;
    bool              M_linearize_equations;
    bool              M_longitudinal_wall;
    //! boolean: compute second spatial derivative of flux?
    bool              M_flux_second_der;
    //! approximation of pressure temporal derivative: how many time steps?
    int               M_dP_dt_steps;

    //! initialize
    int               M_firstNode;
    int               M_lastNode;
    std::string       M_initVar;
    Real              M_restValue;
    Real              M_initValue;
    Real              M_multiplier;
    Real              M_width;
};

}

#endif //ONEDIMENSIONALMODEL_DATA_H
