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

#ifndef OneDimensionalData_H
#define OneDimensionalData_H

// LIFEV
#include <life/lifecore/GetPot.hpp>
#include <life/lifefem/dataTime.hpp>
#include <life/lifemesh/regionMesh1D.hpp>

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>

namespace LifeV
{

//! OneDimensionalData - Class which read and holds all the data for the One Dimensional Model Solver.
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
class OneDimensionalData
{
public:

    //! @name Type definitions
    //@{

    typedef DataTime                                                  time_Type;
    typedef boost::shared_ptr< time_Type >                            timePtr_Type;

    typedef RegionMesh1D< LinearLine >                                mesh_Type;
    typedef boost::shared_ptr< mesh_Type >                            meshPtr_Type;

    // ScalVec SHOULD BE REPLACED EVERYWHERE BY EPETRAVECTOR FOR PARALLEL COMPUTATION
    typedef ublas::vector< Real >                                     scalarVector_Type;
    typedef boost::array< Real, 2 >                                   container2D_Type;

    enum OneD_distributionLaw
    {
        none,
        linear,
        pointwise
    };

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit OneDimensionalData();

    //! Destructor
    virtual ~OneDimensionalData() {}

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

    void setPostprocessingDirectory( const std::string& directory ) {  M_postprocessingDirectory = directory; }

    void setPostprocessingFile( const std::string& file ) { M_postprocessingFile = file; }

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    void setDataTime( const timePtr_Type dataTime ) { M_time = dataTime; }

    void setDensity( const Real& density ) { M_density = density; }

    void setViscosity( const Real& viscosity ) { M_viscosity = viscosity; }

    void setDensityWall( const Real& densityWall ) { M_densityWall = densityWall; }

    void setThickness( const Real& thickness, const UInt& i ) { M_thickness[i] = thickness; }

    void setYoung( const Real& young ) { M_young = young; }

    void setPoisson( const Real& poisson ) { M_poisson = poisson; }

    void setExternalPressure( const Real& externalPressure ) { M_externalPressure = externalPressure; }

    void setArea0( const Real& area0, const UInt& i ) { M_area0[i] = area0; }

    void setBeta0( const Real& beta0, const UInt& i ) { M_beta0[i] = beta0; }

    void setdBeta0dz( const Real& dBeta0dz, const UInt& i ) { M_dBeta0dz[i] = dBeta0dz; }

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics type
    /*!
     * @return Physics Type
     */
    const OneDimensional::physicsType_Type& physicsType() const { return M_physicsType; }

    //! Get the flux type
    /*!
     * @return Flux Type
     */
    const OneDimensional::fluxTerm_Type& fluxType() const { return M_fluxType; }

    //! Get the source type
    /*!
     * @return Source Type
     */
    const OneDimensional::sourceTerm_Type& sourceType() const { return M_sourceType; }

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    timePtr_Type dataTime() const { return M_time; }

    meshPtr_Type mesh() const { return M_mesh; }
    Real length() const { return M_mesh->pointList( M_mesh->numVertices() ).x() - M_mesh->pointList( 1 ).x(); }
    UInt numberOfElements() const { return M_mesh->numElements(); }
    UInt numberOfNodes() const { return M_mesh->numPoints(); }

    const std::string& postprocessingDirectory() const { return M_postprocessingDirectory; }
    const std::string& postprocessingFile() const { return M_postprocessingFile; }

    const Int& verbose() const { return M_verbose; }

    const bool& inertialWall() const { return M_inertialWall; }
    const bool& viscoelasticWall() const { return M_viscoelasticWall; }
    const bool& linearizeStringModel() const { return M_linearizeStringModel; }
    const bool& linearizeEquations() const { return M_linearizeEquations; }
    const bool& longitudinalWall() const { return M_longitudinalWall; }

    const bool& fluxSecondDer() const { return M_fluxSecondDer; }

    const Int& dPdtSteps() const { return M_dP_dt_steps; }
    const Real& CFLmax() const { return M_CFLmax; }

//    const OneD_Initialize& initialVariable() const;
//    const Real&        initialValue() const;
//    const Real&        restValue() const;
//    const Real&        multiplier() const;

    // Physical Parameters
    const Real& densityRho() const { return M_density; }
    const Real& viscosity() const { return M_viscosity; }

    const Real& densityWall() const { return M_densityWall; }
    const Real& young() const { return M_young; }
    const Real& poisson() const { return M_poisson; }

    const Real& externalPressure() const { return M_externalPressure; }

    const Real& viscoelasticModulus() const { return M_viscoelasticModulus; }
    const Real& inertialModulus() const { return M_inertialModulus; }
    const Real& robertsonCorrection() const;

    const Real& thickness( const UInt& i ) const {  return M_thickness[i]; }
    const Real& friction() const { return M_friction; }

    const Real& area0( const UInt& i ) const { return M_area0[i]; }
    const Real& alpha( const UInt& i ) const { return M_alpha[i]; }
    const Real& beta0( const UInt& i ) const { return M_beta0[i]; }
    const Real& beta1( const UInt& i ) const { return M_beta1[i]; }

    const Real& dArea0dz( const UInt& i ) const { return M_dArea0dz[i]; }
    const Real& dAlphadz( const UInt& i ) const { return M_dAlphadz[i]; }
    const Real& dBeta0dz( const UInt& i ) const { return M_dBeta0dz[i]; }
    const Real& dBeta1dz( const UInt& i ) const { return M_dBeta1dz[i]; }

    // Jacobian perturbation parameters
    const Real& jacobianPerturbationArea() const { return M_jacobianPerturbationArea; }
    const Real& jacobianPerturbationFlowRate() const { return M_jacobianPerturbationFlowRate; }
    const Real& jacobianPerturbationPressure() const { return M_jacobianPerturbationPressure; }

    // Linear Parameters
    const Real& flux11( const UInt& i ) const { return M_flux11[i]; }
    const Real& flux12( const UInt& i ) const { return M_flux12[i]; }
    const Real& flux21( const UInt& i ) const { return M_flux21[i]; }
    const Real& flux22( const UInt& i ) const { return M_flux22[i]; }

    const Real& celerity1( const UInt& i ) const { return M_celerity1[i]; }
    const Real& celerity2( const UInt& i ) const { return M_celerity2[i]; }

    const Real& leftEigenVector11( const UInt& i ) const { return M_celerity1LeftEigenvector1[i]; }
    const Real& leftEigenVector12( const UInt& i ) const { return M_celerity1LeftEigenvector2[i]; }
    const Real& leftEigenVector21( const UInt& i ) const { return M_celerity2LeftEigenvector1[i]; }
    const Real& leftEigenVector22( const UInt& i ) const { return M_celerity2LeftEigenvector2[i]; }

    const Real& source10( const UInt& i ) const { return M_source10[i]; }
    const Real& source20( const UInt& i ) const { return M_source20[i]; }
    const Real& source11( const UInt& i ) const { return M_source11[i]; }
    const Real& source12( const UInt& i ) const { return M_source12[i]; }
    const Real& source21( const UInt& i ) const { return M_source21[i]; }
    const Real& source22( const UInt& i ) const { return M_source22[i]; }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Compute the linear interpolation of a quantity.
    /*!
     * Useful for tapering.
     */
    void linearInterpolation( scalarVector_Type& vector,
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
    OneDimensional::physicsType_Type M_physicsType;
    OneDimensional::fluxTerm_Type    M_fluxType;
    OneDimensional::sourceTerm_Type  M_sourceType;

    //! Data containers for time and mesh
    timePtr_Type M_time;
    meshPtr_Type M_mesh;

    //! Miscellaneous
    std::string M_postprocessingDirectory; //! full directory name (including path)
    std::string M_postprocessingFile;      //! output file name
    Int         M_verbose;
    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool        M_inertialWall;
    bool        M_viscoelasticWall;
    bool        M_linearizeStringModel;
    bool        M_linearizeEquations;
    bool        M_longitudinalWall;
    //! boolean: compute second spatial derivative of flux?
    bool        M_fluxSecondDer;
    //! approximation of pressure temporal derivative: how many time steps?
    Int         M_dP_dt_steps;
    Real        M_CFLmax;

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

    scalarVector_Type M_thickness;
    Real M_friction; // Friction parameter

    scalarVector_Type M_area0;      // area
    scalarVector_Type M_alpha;      // Coriolis coefficient (often called alpha)
    scalarVector_Type M_beta0;      // homogeneous to a pressure
    scalarVector_Type M_beta1;      // power coeff (>0, often=1/2)

    //! Derivatives of main coefficients
    scalarVector_Type M_dArea0dz;
    scalarVector_Type M_dAlphadz;
    scalarVector_Type M_dBeta0dz; // homogeneous to a pressure
    scalarVector_Type M_dBeta1dz; // power coeff (>0, often=1/2)

    //! Flux matrix
    scalarVector_Type M_flux11;
    scalarVector_Type M_flux12;
    scalarVector_Type M_flux21;
    scalarVector_Type M_flux22;

    //! Celerities of the linear problem (eigenvalues of the flux matrix)
    scalarVector_Type M_celerity1;
    scalarVector_Type M_celerity2;

    //! Eigenvector for first and second eigenvalue
    scalarVector_Type M_celerity1LeftEigenvector1;
    scalarVector_Type M_celerity1LeftEigenvector2;
    scalarVector_Type M_celerity2LeftEigenvector1;
    scalarVector_Type M_celerity2LeftEigenvector2;

    //! Source matrix
    scalarVector_Type M_source10;
    scalarVector_Type M_source20;
    scalarVector_Type M_source11;
    scalarVector_Type M_source12;
    scalarVector_Type M_source21;
    scalarVector_Type M_source22;
};

} // OneDimensional namespace

#endif //OneDimensionalData_H
