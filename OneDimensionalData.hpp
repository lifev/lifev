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
#include <life/lifefilters/GetPot.hpp>
#include <life/lifefem/TimeData.hpp>
#include <life/lifemesh/RegionMesh1D.hpp>

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalDefinitions.hpp>

namespace LifeV
{

//! OneDimensionalData - Class which read and holds all the data for the One Dimensional Model Solver.
/*!
 *  @authors Vincent Martin, Cristiano Malossi
 *
 *  <b>Physical Parameters</b>
 *
 *  Main parameters: \f$A^0, \alpha, \beta_0, \beta_1, K_r, \rho\f$.
 *
 *  Euler equations:
 *
 *  \f[
 *  \left\{\begin{array}{l}
 *  \displaystyle\frac{dA}{dt} + \frac{dQ}{dz} = 0\\[2ex]
 *  \displaystyle\frac{dQ}{dt} + \frac{d}{dz}\left(\alpha \frac{Q^2}{A}\right) + \frac{A}{\rho} \frac{dP}{dz} + K_r \frac{Q}{A} = 0
 *  \end{array}\right.
 *  \f]
 *
 *  with
 *
 *  \f[
 *  P - P_\mathrm{ext} = \beta_0 \left( \left( \frac{A}{A^0} \right)^{\beta_1} - 1 \right)
 *  \f]
 *
 *  <b>Linear Parameters</b>
 *
 *  Parameters: \f$F_{11}, F_{12}, F_{21}, F_{22}, \lambda_1, \lambda_2\f$
 *
 *  Equations:
 *
 *  \f[
 *  \left\{\begin{array}{l}
 *  \displaystyle\frac{dU_1}{dt} + F_{11} \frac{dU_1}{dz} + F_{12} \frac{dU_2}{dz} = 0\\[2ex]
 *  \displaystyle\frac{dU_2}{dt} + F_{21} \frac{dU_1}{dz} + F_{22} \frac{dU_2}{dz} = 0
 *  \end{array}\right.
 *  \f]
 *
 *  The flux matrix \f$F = [F_{11}, F_{12}; F_{21}, F_{22}]\f$ has the eigenvalues \f$\lambda_1, \lambda_2\f$.
 */
class OneDimensionalData
{
public:

    //! @name Type definitions
    //@{

    typedef TimeData                                                  time_Type;
    typedef boost::shared_ptr< time_Type >                            timePtr_Type;

    typedef RegionMesh1D< LinearLine >                                mesh_Type;
    typedef boost::shared_ptr< mesh_Type >                            meshPtr_Type;

    // ScalVec SHOULD BE REPLACED EVERYWHERE BY EPETRAVECTOR FOR PARALLEL COMPUTATION
    typedef ublas::vector< Real >                                     scalarVector_Type;
    typedef boost::array< Real, 2 >                                   container2D_Type;

    enum OneD_distributionLaw
    {
        uniform,
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

    //! Compute the spatial derivative of a quantity at a node.
    /*!
     * Note: works only for homogeneous discretizations.
     * @param vector the quantity vector
     * @param iNode node
     * @return spatial derivative
     */
    template< typename VectorType >
    Real computeSpatialDerivativeAtNode( const VectorType& vector, const UInt& iNode, const UInt& bcFiniteDifferenceOrder = 2 );

    void initLinearParam( const GetPot& dataFile ); // TO BE CHECKED - DON'T USE IT!

    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    void setPostprocessingDirectory( const std::string& directory ) {  M_postprocessingDirectory = directory; }

    void setPostprocessingFile( const std::string& file ) { M_postprocessingFile = file; }

    //! Set data time container
    /*!
     * @param TimeData shared_ptr to TimeData container
     */
    void setTimeData( const timePtr_Type TimeData ) { M_time = TimeData; }

    void setDensity( const Real& density ) { M_density = density; }

    void setViscosity( const Real& viscosity ) { M_viscosity = viscosity; }

    void setDensityWall( const Real& densityWall ) { M_densityWall = densityWall; }

    void setThickness( const Real& thickness, const UInt& i ) { M_thickness[i] = thickness; }

    void setYoung( const Real& young ) { M_young = young; }

    void setPoisson( const Real& poisson ) { M_poisson = poisson; }

    void setExternalPressure( const Real& externalPressure ) { M_externalPressure = externalPressure; }

    void setVenousPressure( const Real& venousPressure ) { M_venousPressure = venousPressure; }

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
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const { return M_time; }

    meshPtr_Type mesh() const { return M_mesh; }
    Real length() const { return M_mesh->pointList( M_mesh->numVertices() - 1).x() - M_mesh->pointList( 0 ).x(); }
    UInt numberOfElements() const { return M_mesh->numElements(); }
    UInt numberOfNodes() const { return M_mesh->numPoints(); }

    const bool& viscoelasticWall() const { return M_viscoelasticWall; }
    const Real& viscoelasticCoefficient( const UInt& i ) const { return M_viscoelasticCoefficient[i]; }
    const bool& inertialWall() const { return M_inertialWall; }
    const Real& densityWall() const { return M_densityWall; }
    const Real& inertialModulus() const { return M_inertialModulus; }
    const bool& longitudinalWall() const { return M_longitudinalWall; }

    const std::string& postprocessingDirectory() const { return M_postprocessingDirectory; }
    const std::string& postprocessingFile() const { return M_postprocessingFile; }

    const Real& CFLmax() const { return M_CFLmax; }

    // Jacobian perturbation parameters
    const Real& jacobianPerturbationArea() const { return M_jacobianPerturbationArea; }
    const Real& jacobianPerturbationFlowRate() const { return M_jacobianPerturbationFlowRate; }
    const Real& jacobianPerturbationStress() const { return M_jacobianPerturbationStress; }

    // Physical Parameters
    const Real& densityRho() const { return M_density; }
    const Real& viscosity() const { return M_viscosity; }

    const Real& young() const { return M_young; }
    const Real& poisson() const { return M_poisson; }

    const Real& externalPressure() const { return M_externalPressure; }

    //! Get the venous pressure
    /*!
     * @return venous pressure.
     */
    const Real& venousPressure() const { return M_venousPressure; }

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
    void computeSpatialDerivatives();

    //! Reset all the containers.
    void resetContainers();

    //@}

    //! Model
    OneDimensional::physicsType_Type M_physicsType;
    OneDimensional::fluxTerm_Type    M_fluxType;
    OneDimensional::sourceTerm_Type  M_sourceType;

    //! Data containers for time and mesh
    timePtr_Type M_time;
    meshPtr_Type M_mesh;

    //! Physical Wall Model
    bool M_viscoelasticWall;
    Real M_viscoelasticAngle;
    Real M_viscoelasticPeriod;
    scalarVector_Type M_viscoelasticCoefficient;

    bool M_inertialWall;
    Real M_densityWall;
    Real M_inertialModulus;
    bool M_longitudinalWall;

    //! Miscellaneous
    std::string M_postprocessingDirectory; //! full directory name (including path)
    std::string M_postprocessingFile;      //! output file name
    Real        M_CFLmax;

    //! Jacobian perturbation
    Real M_jacobianPerturbationArea;
    Real M_jacobianPerturbationFlowRate;
    Real M_jacobianPerturbationStress;

    //! Physical Parameters
    bool M_computeCoefficients;
    Int  M_powerLawCoefficient;

    Real M_density;     // Density rho (always taken constant along the vessel)
    Real M_viscosity;

    bool M_thickVessel;
    Real M_young;
    Real M_poisson;

    Real M_externalPressure;
    Real M_venousPressure;
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


// ===================================================
// Template implementation
// ===================================================
template< typename VectorType >
inline Real
OneDimensionalData::computeSpatialDerivativeAtNode( const VectorType& vector, const UInt& iNode, const UInt& bcFiniteDifferenceOrder )
{
    // This method is coded only for homogeneous discretizations
    switch ( bcFiniteDifferenceOrder  )
    {
    case 1:

        // We use 1° order finite differences at the boundaries to compute the derivatives
        if ( iNode == 0 )
        {
            return ( -vector[0] + vector[1] ) / ( M_mesh->meanH() );
        }
        else if ( iNode == M_mesh->numPoints() - 1 )
        {
            return ( vector[iNode] - vector[iNode-1] ) / ( M_mesh->meanH() );
        }
        else
        {
            return ( vector[iNode+1] - vector[iNode-1] ) / ( 2.0 * M_mesh->meanH() );
        }

        break;

    case 2:

        // We use 2° order finite differences at the boundaries to compute the derivatives
        if ( iNode == 0 )
        {
            return ( -1.5 * vector[0] + 2.0*vector[1] - 0.5 * vector[2] ) / ( M_mesh->meanH() );
        }
        else if ( iNode == M_mesh->numPoints() - 1 )
        {
            return ( 1.5 * vector[iNode] - 2.0*vector[iNode-1] + 0.5 * vector[iNode-2] ) / ( M_mesh->meanH() );
        }
        else
        {
            return ( vector[iNode+1] - vector[iNode-1] ) / ( 2.0 * M_mesh->meanH() );
        }

        break;

    default:

        std::cout << "!!! Warning: finite difference order \"" << bcFiniteDifferenceOrder << "\"not available!" << std::endl;
        return 0;
    }
}

} // LifeV namespace

#endif //OneDimensionalData_H
