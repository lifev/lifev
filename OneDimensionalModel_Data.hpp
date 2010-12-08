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
    @file
    @brief File containing a class for 1D model data handling.

    @version 1.0
    @author Vincent Martin
    @date 01-07-2004

    @version 2.0
    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @date 12-04-2010
 */

#ifndef ONEDIMENSIONALMODEL_DATA_H
#define ONEDIMENSIONALMODEL_DATA_H

// LIFEV
#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifemesh/regionMesh1D.hpp>

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>

namespace LifeV
{

/*
enum OneD_Initialize { OneD_InitializeArea,
                       OneD_InitializeFlux,
                       OneD_InitializeRiemann1,
                       OneD_InitializeRiemann2,
                       OneD_InitializePressure };
*/

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
 *  P - P_ext = Beta0 [ ( A / Area0 )^{Beta1} - 1 ]
 *
 *  This Beta0 is homogeneous to a pressure.
 *  In most cases Beta1 is taken equal to 1/2.
 *
 *  Beta0 = ( \sqrt{\pi} h_0 E ) / ( ( 1 - \ksi^2 ) * \sqrt{Area0} )
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

    //! @name Type definitions
    //@{

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            Time_PtrType;

    typedef RegionMesh1D< LinearLine >                                Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >                            Mesh_PtrType;

    enum OneD_distributionLaw
    {
        none,
        linear,
        pointwise
    };

    //@}


    //! @name Constructors & Destructor
    //@{

    OneDimensionalModel_Data();

    //! Destructor
    virtual ~OneDimensionalModel_Data() {}

    //@}


    //! @name Methods
    //@{

    //! Setup method.
    /*!
     * @param dataFile GetPot dataFile
     * @return section section in the dataFile
     */
    void setup( const GetPot& dataFile, const std::string& section = "1D_Model" );

    //! Deprecated setup method.
    /*!
     * This method has been implemented for compatibility reason with some old applications.
     * Please, for any new application use the setup method in place of this one. Future
     * compatibility for this method is not guaranteed.
     *
     * @param dataFile GetPot dataFile
     * @return section section in the dataFile
     */
    void oldStyleSetup( const GetPot& dataFile, const std::string& section = "1dnetwork" );

    void UpdateCoefficients();

    void initLinearParam( const GetPot& dataFile ); // TO BE CHECKED - DON'T USE IT!

    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    void setPostprocessingDirectory( const std::string& directory );

    void setPostprocessingFile( const std::string& file );

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    void setDataTime( const Time_PtrType DataTime );

    void setDensity( const Real& Density );

    void setViscosity( const Real& Viscosity );

    void setDensityWall( const Real& DensityWall );

    void setThickness( const Real& Thickness, const UInt& i );

    void setYoung( const Real& Young );

    void setPoisson( const Real& Poisson );

    void setExternalPressure( const Real& externalPressure );

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
    Time_PtrType       dataTime() const;

    Mesh_PtrType       mesh() const;
    Real               Length() const;
    UInt               NumberOfElements() const;
    UInt               NumberOfNodes() const;

    const std::string& postprocessingDirectory() const;
    const std::string& postprocessingFile() const;

    const Int&         verbose() const;

    const bool&        UW() const;
    const bool&        inertialWall() const;
    const bool&        viscoelasticWall() const;
    const bool&        linearizeStringModel() const;
    const bool&        linearizeEquations() const;
    const bool&        longitudinalWall() const;

    const bool&        fluxSecondDer() const;

    const Int&         DPdtSteps() const;
    const Real&        CFLmax() const;

//    const OneD_Initialize& initialVariable() const;
//    const Real&        initialValue() const;
//    const Real&        restValue() const;
//    const Real&        multiplier() const;

    // Physical Parameters
    const Real& DensityRho() const;
    const Real& Viscosity() const;

    const Real& DensityWall() const;
    const Real& Young() const;
    const Real& Poisson() const;

    const Real& externalPressure() const;

    const Real& ViscoelasticModulus() const;
    const Real& InertialModulus() const;
    const Real& RobertsonCorrection() const;

    const Real& Thickness( const UInt& i ) const;
    const Real& Friction() const;

    const Real& Area0( const UInt& i ) const;
    const Real& Alpha( const UInt& i ) const;
    const Real& Beta0( const UInt& i ) const;
    const Real& Beta1( const UInt& i ) const;

    const Real& dArea0dz( const UInt& i ) const;
    const Real& dAlphadz( const UInt& i ) const;
    const Real& dBeta0dz( const UInt& i ) const;
    const Real& dBeta1dz( const UInt& i ) const;

    // Jacobian perturbation parameters
    const Real& JacobianPerturbationArea() const;
    const Real& JacobianPerturbationFlowRate() const;
    const Real& JacobianPerturbationPressure() const;

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

private:

    //! @name Private Methods
    //@{

    //! Compute the linear interpolation of a quantity.
    /*!
     * Useful for tapering.
     */
    void linearInterpolation( ScalVec& vector, const GetPot& dataFile, const std::string& quantity, const Real& defaultValue, const bool& isArea = false );

    //! Compute the derivatives of alpha, area0, beta0, and beta1 using centered differences.
    /*!
     * Note: works only for homogeneous discretizations.
     */
    void computeDerivatives();

    //@}

    //! Model
    OneDimensionalModel_PhysicsTypes M_PhysicsType;
    OneDimensionalModel_FluxTypes    M_FluxType;
    OneDimensionalModel_SourceTypes  M_SourceType;

    //! Data containers for time and mesh
    Time_PtrType      M_Time;
    Mesh_PtrType      M_Mesh;

    //! Miscellaneous
    std::string       M_postprocessingDirectory; //! full directory name (including path)
    std::string       M_postprocessingFile;      //! output file name
    Int               M_verbose;
    bool              M_UW;
    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool              M_inertialWall;
    bool              M_viscoelasticWall;
    bool              M_linearizeStringModel;
    bool              M_linearizeEquations;
    bool              M_longitudinalWall;
    //! boolean: compute second spatial derivative of flux?
    bool              M_fluxSecondDer;
    //! approximation of pressure temporal derivative: how many time steps?
    Int               M_dP_dt_steps;
    Real              M_CFLmax;

    //! initialize
//    OneD_Initialize   M_initialVariable;
//    Real              M_initialValue;
//    Real              M_restValue;
//    Real              M_multiplier;

    //! Jacobian perturbation
    Real M_JacobianPerturbationArea;
    Real M_JacobianPerturbationFlowRate;
    Real M_JacobianPerturbationPressure;

    //! Physical Parameters
    bool M_ComputeCoefficients;
    Int  M_PowerlawCoefficient;

    Real M_Density;     // Density rho (always taken constant along the vessel)
    Real M_Viscosity;

    Real M_DensityWall;
    bool M_ThickVessel;
    Real M_Young;
    Real M_Poisson;

    Real M_externalPressure;

    Real M_ViscoelasticModulus;
    Real M_InertialModulus;
    Real M_RobertsonCorrection;

    ScalVec M_Thickness;
    Real M_Friction; // Friction parameter

    ScalVec M_Area0;      // area
    ScalVec M_Alpha;      // Coriolis coefficient (often called alpha)
    ScalVec M_Beta0;      // homogeneous to a pressure
    ScalVec M_Beta1;      // power coeff (>0, often=1/2)

    //! Derivatives of main coefficients
    ScalVec M_dArea0dz;
    ScalVec M_dAlphadz;
    ScalVec M_dBeta0dz; // homogeneous to a pressure
    ScalVec M_dBeta1dz; // power coeff (>0, often=1/2)

    //! Flux matrix
    ScalVec M_Flux11;
    ScalVec M_Flux12;
    ScalVec M_Flux21;
    ScalVec M_Flux22;

    //! Celerities of the linear problem (eigenvalues of the flux matrix)
    ScalVec M_Celerity1;
    ScalVec M_Celerity2;

    //! Eigenvector for first and second eigenvalue
    ScalVec M_Celerity1LeftEigenvector1;
    ScalVec M_Celerity1LeftEigenvector2;
    ScalVec M_Celerity2LeftEigenvector1;
    ScalVec M_Celerity2LeftEigenvector2;

    //! Source matrix
    ScalVec M_Source10;
    ScalVec M_Source20;
    ScalVec M_Source11;
    ScalVec M_Source12;
    ScalVec M_Source21;
    ScalVec M_Source22;
};

}

#endif //ONEDIMENSIONALMODEL_DATA_H
