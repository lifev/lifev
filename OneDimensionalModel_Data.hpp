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
    @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer  Cristiano Malossi <cristiano.malossi@epfl.ch>
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

    void updateCoefficients();

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

    void setDensity( const Real& density );

    void setViscosity( const Real& discosity );

    void setDensityWall( const Real& densityWall );

    void setThickness( const Real& thickness, const UInt& i );

    void setYoung( const Real& young );

    void setPoisson( const Real& poisson );

    void setExternalPressure( const Real& externalPressure );

    void setArea0( const Real& area0, const UInt& i );

    void setBeta0( const Real& beta0, const UInt& i );

    void setdBeta0dz( const Real& dBeta0dz, const UInt& i );

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics type
    /*!
     * @return Physics Type
     */
    const OneDimensionalModel_PhysicsTypes& physicsType() const;

    //! Get the flux type
    /*!
     * @return Flux Type
     */
    const OneDimensionalModel_FluxTypes& fluxType() const;

    //! Get the source type
    /*!
     * @return Source Type
     */
    const OneDimensionalModel_SourceTypes& sourceType() const;

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    Time_PtrType       dataTime() const;

    Mesh_PtrType       mesh() const;
    Real               length() const;
    UInt               numberOfElements() const;
    UInt               numberOfNodes() const;

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
    const Real& densityRho() const;
    const Real& viscosity() const;

    const Real& densityWall() const;
    const Real& young() const;
    const Real& poisson() const;

    const Real& externalPressure() const;

    const Real& viscoelasticModulus() const;
    const Real& inertialModulus() const;
    const Real& robertsonCorrection() const;

    const Real& thickness( const UInt& i ) const;
    const Real& friction() const;

    const Real& area0( const UInt& i ) const;
    const Real& alpha( const UInt& i ) const;
    const Real& beta0( const UInt& i ) const;
    const Real& beta1( const UInt& i ) const;

    const Real& dArea0dz( const UInt& i ) const;
    const Real& dAlphadz( const UInt& i ) const;
    const Real& dBeta0dz( const UInt& i ) const;
    const Real& dBeta1dz( const UInt& i ) const;

    // Jacobian perturbation parameters
    const Real& jacobianPerturbationArea() const;
    const Real& jacobianPerturbationFlowRate() const;
    const Real& jacobianPerturbationPressure() const;

    // Linear Parameters
    const Real& flux11( const UInt& i ) const;
    const Real& flux12( const UInt& i ) const;
    const Real& flux21( const UInt& i ) const;
    const Real& flux22( const UInt& i ) const;

    const Real& celerity1( const UInt& i ) const;
    const Real& celerity2( const UInt& i ) const;

    const Real& leftEigenVector11( const UInt& i ) const;
    const Real& leftEigenVector12( const UInt& i ) const;
    const Real& leftEigenVector21( const UInt& i ) const;
    const Real& leftEigenVector22( const UInt& i ) const;

    const Real& source10( const UInt& i ) const;
    const Real& source20( const UInt& i ) const;
    const Real& source11( const UInt& i ) const;
    const Real& source12( const UInt& i ) const;
    const Real& source21( const UInt& i ) const;
    const Real& source22( const UInt& i ) const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Compute the linear interpolation of a quantity.
    /*!
     * Useful for tapering.
     */
    void linearInterpolation( scalVec_Type& vector,
                              const GetPot& dataFile,
                              const std::string& quantity,
                              const Real& defaultValue,
                              const bool& isArea = false );

    //! Compute the derivatives of alpha, area0, beta0, and beta1 using centered differences.
    /*!
     * Note: works only for homogeneous discretizations.
     */
    void computeDerivatives();

    //@}

    //! Model
    OneDimensionalModel_PhysicsTypes M_physicsType;
    OneDimensionalModel_FluxTypes    M_fluxType;
    OneDimensionalModel_SourceTypes  M_sourceType;

    //! Data containers for time and mesh
    Time_PtrType      M_time;
    Mesh_PtrType      M_mesh;

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
    Real M_jacobianPerturbationArea;
    Real M_jacobianPerturbationFlowRate;
    Real M_jacobianPerturbationPressure;

    //! Physical Parameters
    bool M_computeCoefficients;
    Int  M_powerLawCoefficient;

    Real M_density;     // Density rho (always taken constant along the vessel)
    Real M_viscosity;

    Real M_densityWall;
    bool M_thickVessel;
    Real M_young;
    Real M_poisson;

    Real M_externalPressure;

    Real M_viscoelasticModulus;
    Real M_inertialModulus;
    Real M_robertsonCorrection;

    scalVec_Type M_thickness;
    Real M_friction; // Friction parameter

    scalVec_Type M_area0;      // area
    scalVec_Type M_alpha;      // Coriolis coefficient (often called alpha)
    scalVec_Type M_beta0;      // homogeneous to a pressure
    scalVec_Type M_beta1;      // power coeff (>0, often=1/2)

    //! Derivatives of main coefficients
    scalVec_Type M_dArea0dz;
    scalVec_Type M_dAlphadz;
    scalVec_Type M_dBeta0dz; // homogeneous to a pressure
    scalVec_Type M_dBeta1dz; // power coeff (>0, often=1/2)

    //! Flux matrix
    scalVec_Type M_flux11;
    scalVec_Type M_flux12;
    scalVec_Type M_flux21;
    scalVec_Type M_flux22;

    //! Celerities of the linear problem (eigenvalues of the flux matrix)
    scalVec_Type M_celerity1;
    scalVec_Type M_celerity2;

    //! Eigenvector for first and second eigenvalue
    scalVec_Type M_celerity1LeftEigenvector1;
    scalVec_Type M_celerity1LeftEigenvector2;
    scalVec_Type M_celerity2LeftEigenvector1;
    scalVec_Type M_celerity2LeftEigenvector2;

    //! Source matrix
    scalVec_Type M_source10;
    scalVec_Type M_source20;
    scalVec_Type M_source11;
    scalVec_Type M_source12;
    scalVec_Type M_source21;
    scalVec_Type M_source22;
};

}

#endif //ONEDIMENSIONALMODEL_DATA_H
