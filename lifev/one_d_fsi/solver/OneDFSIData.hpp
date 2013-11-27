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

#ifndef OneDFSIData_H
#define OneDFSIData_H


#include <Epetra_SerialComm.h>


#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/TimeData.hpp>
#include <lifev/core/mesh/RegionMesh1DBuilders.hpp>

#include <lifev/one_d_fsi/solver/OneDFSIDefinitions.hpp>

namespace LifeV
{

//! OneDFSIData - Class which read and holds all the data for the One Dimensional Model Solver.
/*!
 *  @authors Vincent Martin, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite BonnemainMalossi2012LVAD
 *
 *  <b>Physical Parameters</b>
 *
 *  Main parameters: \f$A^0, \alpha, \beta_0, \beta_1, \gamma, K_r, \rho\f$.
 *
 *  Euler equations:
 *
 *  \f[
 *  \left\{\begin{array}{l}
 *  \displaystyle\frac{\partial A}{\partial t} + \frac{\partial Q}{\partial z} = 0, \\[2ex]
 *  \displaystyle\frac{\partial Q}{\partial t} +
 *  \alpha \frac{\partial}{\partial z}\left(\frac{Q^2}{A}\right) +
 *  \frac{A}{\rho}\frac{\partial P}{\partial z} + K_r \frac{Q}{A} = 0
 *  \end{array}\right.
 *  \f]
 *
 *  with
 *
 *  \f[
 *  P-P_\mathrm{ext} = \psi(A,A^0,\beta_0, \beta_1, \gamma) =
 *  \underbrace{\sqrt{\frac{\pi}{A^0}}\frac{h E}{1-\nu^2}}_{\beta_0} \left(\left(\frac{A}{A^0}\right)^{\beta_1}-1\right) +
 *  \underbrace{\frac{T \tan\phi}{4 \sqrt{\pi}}\frac{h E}{1-\nu^2}}_{\displaystyle\gamma} \frac{1}{A\sqrt{A}} \frac{\partial A}{\partial t},
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
 *
 *  The flux matrix \f$\mathbf F = [F_{11}, F_{12}; F_{21}, F_{22}]\f$ has the eigenvalues \f$\lambda_1, \lambda_2\f$.
 */
class OneDFSIData
{
public:

    //! @name Type definitions
    //@{

    typedef TimeData                                                  time_Type;
    typedef boost::shared_ptr< time_Type >                            timePtr_Type;

    typedef RegionMesh< LinearLine >                                  mesh_Type;
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

    //! Empty constructor
    explicit OneDFSIData();

    //! Destructor
    virtual ~OneDFSIData() {}

    //@}


    //! @name Methods
    //@{

    //! Setup method.
    /*!
     * @param dataFile GetPot dataFile
     * @param section section in the dataFile
     */
    void setup ( const GetPot& dataFile, const std::string& section = "1D_Model" );

    //! <b>Deprecated</b> setup method. (<b>DO NOT USE THIS</b>)
    /*!
     * This method has been implemented for compatibility reason with some old applications.
     * Please, for any new application use the setup method in place of this one. Future
     * compatibility for this method is not guaranteed.
     *
     * @param dataFile GetPot dataFile
     * @return section section in the dataFile
     */
    void LIFEV_DEPRECATED ( oldStyleSetup ( const GetPot& dataFile, const std::string& section = "1dnetwork" ) );

    //! Update all the physical coefficients
    /*!
     * This method can be called after any update in the physical coefficient.
     * It recomputes the main coefficients \f$\alpha, \beta_0, \beta_1, \gamma, K_r\f$
     */
    void updateCoefficients();

    //! Compute the spatial derivative of a quantity at a node.
    /*!
     * Note: works only for homogeneous discretizations.
     * @param vector the quantity vector
     * @param iNode node
     * @return spatial derivative
     */
    template< typename VectorType >
    Real computeSpatialDerivativeAtNode ( const VectorType& vector, const UInt& iNode, const UInt& bcFiniteDifferenceOrder = 2 );

    //! Initialize linear parameters (<b>NOT WORKING</b>)
    /*!
     * The linearization of the Euler model yields
     * \f[
     * \begin{array}{l}
     * \displaystyle F = [ Q; A c^2]\\[2ex]
     * \displaystyle B = [ 0; k_R / A^0]\\[2ex]
     * \displaystyle c = \sqrt{\frac{\beta_0\beta_1}{\rho} }
     * \end{array}
     * \f]
     */
    //    void initializeLinearParameters();

    //! Make the vessel stiffer on the left side of interval [xl, xr]
    /*!
     *  \cond \TODO improve doxygen description with latex equation and other features \endcond
     *  These routines change the elastic modulus along the vessel
     *
     *  When x < alpha - delta/2, the Young modulus is E * factor
     *  When x > alpha + delta/2, the Young modulus is E
     *  When alpha - delta/2 < x < alpha + delta/2, the Young modulus changes
     *  smoothly from the larger to the smaller value, according to a
     *  polynomial law of order n.
     *
     *  The grid size can be adapted (yesadaptive=1) in the nieghborhood of alpha,
     *  where the spatial derivative of the parameter will be maximum.
     *  However, the grid size is not allowed to be smaller than min_deltax
     *
     *  \cond \TODO add doxygen description for the parameters \endcond
     */
    //    void stiffenVesselLeft( const Real& xl,          const Real& xr,
    //                            const Real& factor,      const Real& alpha,
    //                            const Real& delta,       const Real& n,
    //                            const Real& minDeltaX=1, const UInt& yesAdaptive=0 );

    //! Make the vessel stiffer on the right side of interval [xl, xr]
    /*!
     * \sa stiffenVesselLeft
     *
     *  \cond \TODO add doxygen description for the parameters \endcond
     */
    //    void stiffenVesselRight( const Real& xl,          const Real& xr,
    //                             const Real& factor,      const Real& alpha,
    //                             const Real& delta,       const Real& n,
    //                             const Real& minDeltaX=1, const UInt& yesAdaptive=0  );

    //! Display some information about the model.
    /*!
     * @param output Stream where the informations must be printed
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the post-processing directory
    /*!
     * @param directory post-processing directory
     */
    void setPostprocessingDirectory ( const std::string& directory )
    {
        M_postprocessingDirectory = directory;
    }

    //! Set the post-processing file name
    /*!
     * @param file post-processing file name
     */
    void setPostprocessingFile ( const std::string& file )
    {
        M_postprocessingFile = file;
    }

    //! Set the area perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @param jacobianPerturbationArea area perturbation parameter
     */
    void setJacobianPerturbationArea ( const Real& jacobianPerturbationArea )
    {
        M_jacobianPerturbationArea = jacobianPerturbationArea;
    }

    //! Set the flow rate perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @param jacobianPerturbationFlowRate flow rate perturbation parameter
     */
    void setJacobianPerturbationFlowRate ( const Real& jacobianPerturbationFlowRate )
    {
        M_jacobianPerturbationFlowRate = jacobianPerturbationFlowRate;
    }

    //! Set the stress perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @param jacobianPerturbationStress stress perturbation parameter
     */
    void setJacobianPerturbationStress ( const Real& jacobianPerturbationStress )
    {
        M_jacobianPerturbationStress = jacobianPerturbationStress;
    }

    //! Set data time container
    /*!
     * @param timeDataPtr shared_ptr to TimeData container
     */
    void setTimeData ( const timePtr_Type timeDataPtr )
    {
        M_timeDataPtr = timeDataPtr;
    }

    //! Set the fluid density
    /*!
     * @param density fluid density
     */
    void setDensity ( const Real& density )
    {
        M_density = density;
    }

    //! Set the fluid viscosity
    /*!
     * @param viscosity fluid viscosity
     */
    void setViscosity ( const Real& viscosity )
    {
        M_viscosity = viscosity;
    }

    //! Set the wall density
    /*!
     * @param densityWall wall density
     */
    void setDensityWall ( const Real& densityWall )
    {
        M_densityWall = densityWall;
    }

    //! Set the wall thickness
    /*!
     * @param i node id
     * @param thickness wall thickness
     */
    void setThickness ( const Real& thickness, const UInt& i )
    {
        M_thickness[i] = thickness;
    }

    //! Set the wall Young modulus
    /*!
     * @param young wall Young modulus
     */
    void setYoung ( const Real& young )
    {
        M_young = young;
    }

    //! Set the wall Poisson number
    /*!
     * @param poisson wall Poisson number
     */
    void setPoisson ( const Real& poisson )
    {
        M_poisson = poisson;
    }

    //! Set the wall external pressure
    /*!
     * @param externalPressure wall external pressure
     */
    void setExternalPressure ( const Real& externalPressure )
    {
        M_externalPressure = externalPressure;
    }

    //! Set the venous pressure at the terminals
    /*!
     * @param venousPressure venous pressure at the terminals
     */
    void setVenousPressure ( const Real& venousPressure )
    {
        M_venousPressure = venousPressure;
    }

    //! Set the reference area
    /*!
     * @param i node id
     * @param area0 reference area
     */
    void setArea0 ( const Real& area0, const UInt& i )
    {
        M_area0[i] = area0;
    }

    //! Set the wall \f$\beta_0\f$
    /*!
     * @param i node id
     * @param beta0 wall \f$\beta_0\f$
     */
    void setBeta0 ( const Real& beta0, const UInt& i )
    {
        M_beta0[i] = beta0;
    }

    //! Set the wall \f$\frac{d\beta_0}{dz}\f$
    /*!
     * @param i node id
     * @param thickness wall \f$\frac{d\beta_0}{dz}\f$
     */
    void setdBeta0dz ( const Real& dBeta0dz, const UInt& i )
    {
        M_dBeta0dz[i] = dBeta0dz;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics type
    /*!
     * @return Physics type
     */
    const OneDFSI::physicsType_Type& physicsType() const
    {
        return M_physicsType;
    }

    //! Get the flux type
    /*!
     * @return Flux type
     */
    const OneDFSI::fluxTerm_Type& fluxType() const
    {
        return M_fluxType;
    }

    //! Get the source type
    /*!
     * @return Source type
     */
    const OneDFSI::sourceTerm_Type& sourceType() const
    {
        return M_sourceType;
    }

    //! Get data time container
    /*!
     * @return shared_ptr to TimeData container
     */
    timePtr_Type dataTime() const
    {
        return M_timeDataPtr;
    }

    //! Get the mesh container
    /*!
     * @return shared_ptr to the mesh
     */
    meshPtr_Type mesh() const
    {
        return M_meshPtr;
    }

    //! Get the length of the 1D segment
    /*!
     * @return length of the 1D segment
     */
    Real length() const
    {
        return M_meshPtr->pointList ( M_meshPtr->numVertices() - 1).x() - M_meshPtr->pointList ( 0 ).x();
    }

    //! Get the number of elements in the 1D segment
    /*!
     * @return number of elements in the 1D segment
     */
    UInt numberOfElements() const
    {
        return M_meshPtr->numElements();
    }

    //! Get the number of nodes in the 1D segment
    /*!
     * @return number of nodes in the 1D segment
     */
    UInt numberOfNodes() const
    {
        return M_meshPtr->numPoints();
    }

    //! Get the flag identifying if the wall is viscoelastic
    /*!
     * @return true if the wall is viscoelastic, false otherwise
     */
    const bool& viscoelasticWall() const
    {
        return M_viscoelasticWall;
    }

    //! Get the viscoelastic coefficient \f$\gamma\f$
    /*!
     * @return viscoelastic coefficient \f$\gamma\f$
     */
    const Real& viscoelasticCoefficient ( const UInt& i ) const
    {
        return M_viscoelasticCoefficient[i];
    }

    //! Get the flag identifying if the wall has inertia
    /*!
     * @return true if the wall has intertia, false otherwise
     */
    const bool& inertialWall() const
    {
        return M_inertialWall;
    }

    //! Get the density of the wall
    /*!
     * @return density of the wall
     */
    const Real& densityWall() const
    {
        return M_densityWall;
    }

    //! Get the inertial coefficient (to be defined)
    /*!
     * @return inertial coefficient (to be defined)
     */
    const Real& inertialModulus() const
    {
        return M_inertialModulus;
    }

    //! Get the flag identifying if the wall has a longitudinal pre-stress
    /*!
     * @return true if the wall has a longitudinal pre-stress, false otherwise
     */
    const bool& longitudinalWall() const
    {
        return M_longitudinalWall;
    }

    //! Get the post-processing directory
    /*!
     * @return post-processing directory
     */
    const std::string& postprocessingDirectory() const
    {
        return M_postprocessingDirectory;
    }

    //! Get the post-processing file
    /*!
     * @return post-processing file
     */
    const std::string& postprocessingFile() const
    {
        return M_postprocessingFile;
    }

    //! Get the imposed CFL condition
    /*!
     * @return imposed CFL condition
     */
    const Real& CFLmax() const
    {
        return M_CFLmax;
    }

    //! Get the area perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @return area perturbation parameter
     */
    const Real& jacobianPerturbationArea() const
    {
        return M_jacobianPerturbationArea;
    }

    //! Get the flow rate perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @return flow rate perturbation parameter
     */
    const Real& jacobianPerturbationFlowRate() const
    {
        return M_jacobianPerturbationFlowRate;
    }

    //! Get the stress perturbation parameter to compute the Jacobian matrix (in the Multiscale framework)
    /*!
     * @return stress perturbation parameter
     */
    const Real& jacobianPerturbationStress() const
    {
        return M_jacobianPerturbationStress;
    }

    //! Get the fluid density
    /*!
     * @return fluid density
     */
    const Real& densityRho() const
    {
        return M_density;
    }

    //! Get the fluid viscosity
    /*!
     * @return fluid viscosity
     */
    const Real& viscosity() const
    {
        return M_viscosity;
    }

    //! Get the wall Young modulus
    /*!
     * @return wall Young modulus
     */
    const Real& young() const
    {
        return M_young;
    }

    //! Get the wall Poisson number
    /*!
     * @return wall Poisson number
     */
    const Real& poisson() const
    {
        return M_poisson;
    }

    //! Get the wall external pressure
    /*!
     * @return wall external pressure
     */
    const Real& externalPressure() const
    {
        return M_externalPressure;
    }

    //! Get the venous pressure
    /*!
     * @return venous pressure.
     */
    const Real& venousPressure() const
    {
        return M_venousPressure;
    }

    //! Get the Robertson correction coefficient (<b>Not tested: maybe wrong in the code</b>)
    /*!
     * @return Robertson correction coefficient
     */
    const Real& robertsonCorrection() const;

    //! Get the wall thickness
    /*!
     * @return wall thickness
     */
    const Real& thickness ( const UInt& i ) const
    {
        return M_thickness[i];
    }

    //! Get the wall friction coefficient \f$K_r\f$
    /*!
     * @return wall friction coefficient \f$K_r\f$
     */
    const Real& friction() const
    {
        return M_friction;
    }

    //! Get the reference area \f$A^0\f$
    /*!
     * @param i node id
     * @return reference area \f$A^0\f$
     */
    const Real& area0 ( const UInt& i ) const
    {
        return M_area0[i];
    }

    //! Get the Coriolis coefficient \f$\alpha\f$
    /*!
     * @param i node id
     * @return Coriolis coefficient \f$\alpha\f$
     */
    const Real& alpha ( const UInt& i ) const
    {
        return M_alpha[i];
    }

    //! Get the \f$\beta_0\f$ coefficient
    /*!
     * @param i node id
     * @return \f$\beta_0\f$ coefficient
     */
    const Real& beta0 ( const UInt& i ) const
    {
        return M_beta0[i];
    }

    //! Get the \f$\beta_1\f$ coefficient
    /*!
     * @param i node id
     * @return \f$\beta_1\f$ coefficient
     */
    const Real& beta1 ( const UInt& i ) const
    {
        return M_beta1[i];
    }

    //! Get \f$\frac{dA^0}{dz}\f$
    /*!
     * @param i node id
     * @return \f$\frac{dA^0}{dz}\f$
     */
    const Real& dArea0dz ( const UInt& i ) const
    {
        return M_dArea0dz[i];
    }

    //! Get \f$\frac{d\alpha}{dz}\f$
    /*!
     * @param i node id
     * @return \f$\frac{d\alpha}{dz}\f$
     */
    const Real& dAlphadz ( const UInt& i ) const
    {
        return M_dAlphadz[i];
    }

    //! Get \f$\frac{d\beta_0}{dz}\f$
    /*!
     * @param i node id
     * @return \f$\frac{d\beta_0}{dz}\f$
     */
    const Real& dBeta0dz ( const UInt& i ) const
    {
        return M_dBeta0dz[i];
    }

    //! Get \f$\frac{d\beta_1}{dz}\f$
    /*!
     * @param i node id
     * @return \f$\frac{d\beta_1}{dz}\f$
     */
    const Real& dBeta1dz ( const UInt& i ) const
    {
        return M_dBeta1dz[i];
    }

    //! Get the flux coefficient \f$F_{11}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$F_{11}\f$
     */
    const Real& flux11 ( const UInt& i ) const
    {
        return M_flux11[i];
    }

    //! Get the flux coefficient \f$F_{12}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$F_{12}\f$
     */
    const Real& flux12 ( const UInt& i ) const
    {
        return M_flux12[i];
    }

    //! Get the flux coefficient \f$F_{21}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$F_{21}\f$
     */
    const Real& flux21 ( const UInt& i ) const
    {
        return M_flux21[i];
    }

    //! Get the flux coefficient \f$F_{22}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$F_{22}\f$
     */
    const Real& flux22 ( const UInt& i ) const
    {
        return M_flux22[i];
    }

    //! Get the first eigenvector \f$\lambda_{1}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$\lambda_{1}\f$
     */
    const Real& celerity1 ( const UInt& i ) const
    {
        return M_celerity1[i];
    }

    //! Get the second eigenvector \f$\lambda_{2}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$\lambda_{2}\f$
     */
    const Real& celerity2 ( const UInt& i ) const
    {
        return M_celerity2[i];
    }

    //! Get the left eigenvector coefficient \f$L_{11}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$L_{11}\f$
     */
    const Real& leftEigenVector11 ( const UInt& i ) const
    {
        return M_celerity1LeftEigenvector1[i];
    }

    //! Get the left eigenvector coefficient \f$L_{12}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$L_{12}\f$
     */
    const Real& leftEigenVector12 ( const UInt& i ) const
    {
        return M_celerity1LeftEigenvector2[i];
    }

    //! Get the left eigenvector coefficient \f$L_{21}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$L_{21}\f$
     */
    const Real& leftEigenVector21 ( const UInt& i ) const
    {
        return M_celerity2LeftEigenvector1[i];
    }

    //! Get the left eigenvector coefficient \f$L_{22}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$L_{22}\f$
     */
    const Real& leftEigenVector22 ( const UInt& i ) const
    {
        return M_celerity2LeftEigenvector2[i];
    }

    //! Get the source coefficient \f$S_{10}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{10}\f$
     */
    const Real& source10 ( const UInt& i ) const
    {
        return M_source10[i];
    }

    //! Get the source coefficient \f$S_{20}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{20}\f$
     */
    const Real& source20 ( const UInt& i ) const
    {
        return M_source20[i];
    }

    //! Get the source coefficient \f$S_{11}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{11}\f$
     */
    const Real& source11 ( const UInt& i ) const
    {
        return M_source11[i];
    }

    //! Get the source coefficient \f$S_{12}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{12}\f$
     */
    const Real& source12 ( const UInt& i ) const
    {
        return M_source12[i];
    }

    //! Get the source coefficient \f$S_{21}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{21}\f$
     */
    const Real& source21 ( const UInt& i ) const
    {
        return M_source21[i];
    }

    //! Get the source coefficient \f$S_{22}\f$ (used only in the linear problem)
    /*!
     * @param i node id
     * @return \f$S_{22}\f$
     */
    const Real& source22 ( const UInt& i ) const
    {
        return M_source22[i];
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Compute the linear interpolation of a general quantity.
    /*!
     *  Very useful for tapering.
     *  @param vector interpolated vector
     *  @param dataFile data file
     *  @param quantity quantity
     *  @param defaultValue default value
     *  @param isArea flag identifying if the vector is an area (the linear interpolation is done on the radius)
     */
    void linearInterpolation ( scalarVector_Type& vector,
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
    OneDFSI::physicsType_Type M_physicsType;
    OneDFSI::fluxTerm_Type    M_fluxType;
    OneDFSI::sourceTerm_Type  M_sourceType;

    //! Data containers for time and mesh
    timePtr_Type M_timeDataPtr;
    meshPtr_Type M_meshPtr;

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
OneDFSIData::computeSpatialDerivativeAtNode ( const VectorType& vector, const UInt& iNode, const UInt& bcFiniteDifferenceOrder )
{
    // This method is coded only for homogeneous discretizations

    Real meanH = MeshUtility::MeshStatistics::computeSize (*M_meshPtr).meanH;
    switch ( bcFiniteDifferenceOrder  )
    {
        case 1:

            // We use 1° order finite differences at the boundaries to compute the derivatives
            if ( iNode == 0 )
            {
                return ( -vector[0] + vector[1] ) / ( meanH );
            }
            else if ( iNode == M_meshPtr->numPoints() - 1 )
            {
                return ( vector[iNode] - vector[iNode - 1] ) / ( meanH );
            }
            else
            {
                return ( vector[iNode + 1] - vector[iNode - 1] ) / ( 2.0 * meanH );
            }

            break;

        case 2:

            // We use 2° order finite differences at the boundaries to compute the derivatives
            if ( iNode == 0 )
            {
                return ( -1.5 * vector[0] + 2.0 * vector[1] - 0.5 * vector[2] ) / ( meanH );
            }
            else if ( iNode == M_meshPtr->numPoints() - 1 )
            {
                return ( 1.5 * vector[iNode] - 2.0 * vector[iNode - 1] + 0.5 * vector[iNode - 2] ) / ( meanH );
            }
            else
            {
                return ( vector[iNode + 1] - vector[iNode - 1] ) / ( 2.0 * meanH );
            }

            break;

        default:

            std::cout << "!!! Warning: finite difference order \"" << bcFiniteDifferenceOrder << "\"not available!" << std::endl;
            return 0;
    }
}

} // LifeV namespace

#endif //OneDFSIData_H
