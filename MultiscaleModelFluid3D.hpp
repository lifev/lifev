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
 *  @brief File containing the MultiScale Model Fluid3D
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */


#ifndef MultiscaleModelFluid3D_H
#define MultiscaleModelFluid3D_H 1

// LifeV includes
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/hdf5exporter.hpp>
#endif
#include <life/lifesolver/OseenShapeDerivative.hpp>

// Mathcard includes
#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>

namespace LifeV
{
namespace multiscale
{

//! MultiscaleModelFluid3D - MultiScale model for 3D Fluid simulations
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleModelFluid3D class is an implementation of the multiscaleModel_Type
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 */
class MultiscaleModelFluid3D: public virtual multiscaleModel_Type
{
public:

    //! @name Public Types
    //@{

    typedef RegionMesh3D< LinearTetra >           mesh_Type;
    typedef partitionMesh< mesh_Type >            partitionMesh_Type;

    typedef OseenShapeDerivative< mesh_Type >     fluid_Type;
    typedef fluid_Type::vector_type               fluidVector_Type;
    typedef boost::shared_ptr< fluidVector_Type > fluidVectorPtr_Type;

    typedef Exporter< mesh_Type >                 IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >      IOFilePtr_Type;

    typedef Ensight< mesh_Type >                  ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef Hdf5exporter< mesh_Type >             hdf5IOFile_Type;
#endif

    typedef BCHandler                             bc_Type;
    typedef boost::shared_ptr< bc_Type >          bcPtr_Type;
    typedef BCInterface< fluid_Type >             bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type > bcInterfacePtr_Type;

    typedef BdfTNS< fluidVector_Type >            bdf_Type;
    typedef DataNavierStokes		              data_Type;

    typedef FESpace< mesh_Type, EpetraMap >       FESpace_Type;
    typedef boost::shared_ptr< FESpace_Type >     FESpacePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFluid3D();

    //! Destructor
    virtual ~MultiscaleModelFluid3D() {}

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial system (matrix and vectors).
    void buildSystem();

    //! Update the system for (matrix and vectors).
    void updateSystem();

    //! Solve the problem.
    void solveSystem();

    //! save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //@}


    //! @name Methods
    //@{

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

    //@}


    //! @name Set Methods
    //@{

    //! Set the solution vector
    /*!
     * @param solution solution vector
     */
    void setSolution( const fluidVectorPtr_Type& solution );

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() { return *M_bc; }

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real boundaryDensity( const BCFlag& /*flag*/) const { return M_data->density(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const BCFlag& /*flag*/) const { return M_data->viscosity(); }

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real boundaryArea( const BCFlag& flag ) const { return M_fluid->area( flag ); }

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    Real boundaryFlowRate( const BCFlag& flag ) const { return M_fluid->flux( flag ); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const BCFlag& flag ) const { return M_fluid->pressure( flag ); }

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real boundaryDynamicPressure( const BCFlag& flag ) const { return 0.5 * boundaryDensity( flag ) * ( boundaryFlowRate( flag ) * boundaryFlowRate( flag ) ) / ( boundaryArea( flag ) * boundaryArea( flag ) ); }

    //! Get the value of the Lagrange multiplier associated to a specific boundary face
    /*!
     * @param flag flag of the boundary face
     * @return Lagrange multiplier value
     */
    Real boundaryLagrangeMultiplier( const BCFlag& flag ) const { return M_fluid->LagrangeMultiplier(flag, *M_bc->handler() ); }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const BCFlag& flag, const stress_Type& stressType = StaticPressure ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real boundaryDeltaPressure( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the total pressure (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the dynamic pressure
     */
    Real boundaryDeltaDynamicPressure( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the Lagrange multiplier associated to a specific boundary face, using the linear model
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return Lagrange multiplier value
     */
    Real boundaryDeltaLagrangeMultiplier( const BCFlag& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const BCFlag& flag, bool& solveLinearSystem, const stress_Type& stressType = StaticPressure );

    //@}


    //! @name Get Methods
    //@{

    //! Get the data container
    /*!
     * @return Data container
     */
    const data_Type& data() const { return *M_data; }

    //! Get the solution vector
    /*!
     * @return Solution vector
     */
    const fluidVector_Type& solution() const { return *M_solution; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelFluid3D( const MultiscaleModelFluid3D& model );

    MultiscaleModelFluid3D& operator=( const MultiscaleModelFluid3D& model );

    //@}


    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param fileName File name of the specific model.
     */
    void setupGlobalData( const std::string& fileName );

    //! Setup the exporter and the importer
    /*!
     * @param fileName File name of the specific model.
     */
    void setupExporterImporter( const std::string& fileName );

    //! Setup the mesh for the fluid problem
    void setupMesh();

    //! Setup the FE space for pressure and velocity
    void setupFEspace();

    //! Setup the DOF of the model
    void setupDOF();

    //! Setup the offset for fluxes boundary conditions
    void setupBCOffset( const bcPtr_Type& BC );

    //! Initialize the solution.
    void initializeSolution();

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void imposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void resetPerturbation();

    Real bcFunctionDeltaZero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/) { return 0.; }
    Real bcFunctionDeltaOne ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/) { return 1.; }

    //@}

    IOFilePtr_Type                          M_exporter;
    IOFilePtr_Type                          M_importer;

    std::string                             M_fileName; //TODO Temporary - To be removed

    // Fluid problem
    boost::shared_ptr< fluid_Type >         M_fluid;
    bcInterfacePtr_Type                     M_bc;
    boost::shared_ptr< bdf_Type >           M_bdf;
    boost::shared_ptr< data_Type >          M_data;
    boost::shared_ptr< DataMesh >           M_dataMesh;
    boost::shared_ptr< partitionMesh_Type > M_mesh;
    boost::shared_ptr< EpetraMap >          M_map;
    fluidVectorPtr_Type                     M_solution;

    // Linear Fluid problem
    bcPtr_Type                              M_linearBC;
    bool                                    M_updateLinearModel;

    // FE spaces
    FESpacePtr_Type                         M_uFESpace;
    FESpacePtr_Type                         M_pFESpace;

    // Lagrange multipliers
    UInt                                    M_lmDOF;

    // Problem coefficients
    Real                                    M_alpha;
    fluidVectorPtr_Type                     M_beta;
    fluidVectorPtr_Type                     M_rhs;

    // NS parameters
    UInt                                    M_subiterationsMaximumNumber;
    Real                                    M_tolerance;
    generalizedAitken< fluidVector_Type >   M_generalizedAitken;

    // BC Functions for tangent problem
    BCFunctionBase                          M_bcBaseDeltaZero;
    BCFunctionBase                          M_bcBaseDeltaOne;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFluid3D()
{
    return new MultiscaleModelFluid3D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelFluid3D_H */
