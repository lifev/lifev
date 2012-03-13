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
 *  @brief File containing the Multiscale Model Fluid3D
 *
 *  @date 12-03-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */


#ifndef MultiscaleModelFluid3D_H
#define MultiscaleModelFluid3D_H 1

// LifeV includes
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
<<<<<<< HEAD:lifev/multiscale/solver/MultiscaleModelFluid3D.hpp
#include <lifev/navier_stokes/solver/OseenData.hpp>
=======
#include <life/lifesolver/OseenData.hpp>
>>>>>>> Updating headers for core module:lifemc/lifesolver/MultiscaleModelFluid3D.hpp
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/bc_interface/fem/BCInterface3D.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>

// Mathcard includes

#include <lifev/multiscale/solver/MultiscaleModel.hpp>
#include <lifev/multiscale/solver/MultiscaleInterfaceFluid.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelFluid3D - Multiscale model for 3D Fluid simulations
/*!
 *  @author Cristiano Malossi
 *
 *  The MultiscaleModelFluid3D class is an implementation of the multiscaleModel_Type
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 */
class MultiscaleModelFluid3D: public virtual multiscaleModel_Type,
                              public virtual MultiscaleInterfaceFluid
{
public:

    //! @name Public Types
    //@{

    typedef RegionMesh< LinearTetra >                         mesh_Type;
    typedef MeshPartitioner< mesh_Type >                      MeshPartitioner_Type;


    typedef OseenSolverShapeDerivative< mesh_Type >           fluid_Type;
    typedef fluid_Type::vector_Type                           fluidVector_Type;
    typedef boost::shared_ptr< fluidVector_Type >             fluidVectorPtr_Type;

    typedef Exporter< mesh_Type >                             IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >                  IOFilePtr_Type;
    typedef ExporterData< mesh_Type >                         IOData_Type;

    typedef ExporterEnsight< mesh_Type >                      ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                         hdf5IOFile_Type;
#endif

    typedef BCHandler                                         bc_Type;
    typedef boost::shared_ptr< bc_Type >                      bcPtr_Type;
    typedef BCInterface3D< bc_Type, fluid_Type >              bcInterface_Type;
    typedef boost::shared_ptr< bcInterface_Type >             bcInterfacePtr_Type;

    typedef TimeAdvanceBDFNavierStokes< fluidVector_Type >    bdf_Type;
    typedef OseenData                                         data_Type;

    typedef FESpace< mesh_Type, MapEpetra >                   FESpace_Type;
    typedef boost::shared_ptr< FESpace_Type >                 FESpacePtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFluid3D();

    //! Destructor
    virtual ~MultiscaleModelFluid3D() {}

    //@}


    //! @name MultiscaleModel Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param fileName Name of data file.
     */
    void setupData( const std::string& fileName );

    //! Setup the model.
    void setupModel();

    //! Build the initial model.
    void buildModel();

    //! Update the model.
    void updateModel();

    //! Solve the model.
    void solveModel();

    //! Update the solution.
    void updateSolution();

    //! Save the solution
    void saveSolution();

    //! Display some information about the model.
    void showMe();

    //! Return a specific scalar quantity to be used for a comparison with a reference value.
    /*!
     * This method is meant to be used for night checks.
     * @return reference quantity.
     */
    Real checkSolution() const;

    //@}


    //! @name MultiscaleInterfaceFluid Methods
    //@{

    //! Impose the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate( const bcFlag_Type& flag, const function_Type& function );

    //! Impose the integral of the normal stress on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @param function boundary condition function
     */
    void imposeBoundaryStress( const bcFlag_Type& flag, const function_Type& function );

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    Real boundaryFlowRate( const bcFlag_Type& flag ) const { return M_fluid->flux( flag ); }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param stressType Type of approximation for the stress
     * @return stress value
     */
    Real boundaryStress( const bcFlag_Type& flag ) const { return -boundaryPressure( flag ); }

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return stress value
     */
    Real boundaryTotalStress( const bcFlag_Type& flag ) const { return -boundaryTotalPressure( flag ); }

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param flag flag of the boundary face on which quantity should be computed
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real boundaryDeltaStress( const bcFlag_Type& flag, bool& solveLinearSystem );

    //! Get the variation of the integral of the total normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @param stressType Type of approximation for the stress
     * @return variation of the total normal stress
     */
    Real boundaryDeltaTotalStress( const bcFlag_Type& flag, bool& solveLinearSystem );

    //@}


    //! @name Set Methods
    //@{

    //! Set the solution vector
    /*!
     * @param solution solution vector
     */
    void setSolution( const fluidVectorPtr_Type& solution );

    //@}


    //! @name Get Methods
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
    Real boundaryDensity( const bcFlag_Type& /*flag*/) const { return M_data->density(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const bcFlag_Type& /*flag*/) const { return M_data->viscosity(); }

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real boundaryArea( const bcFlag_Type& flag ) const { return M_fluid->area( flag ); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const bcFlag_Type& flag ) const;

    //! Get the integral of the total pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return total pressure value
     */
    Real boundaryTotalPressure( const bcFlag_Type& flag ) const;

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

    //! Initialize the solution.
    void initializeSolution();

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

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void imposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void resetPerturbation();

    Real bcFunctionDeltaZero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ ) { return 0.; }
    Real bcFunctionDeltaOne ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ ) { return 1.; }

    //@}

    IOFilePtr_Type                          M_exporter;
    IOFilePtr_Type                          M_importer;

    std::string                             M_fileName; //TODO Temporary - To be removed

    // Fluid problem
    boost::shared_ptr< fluid_Type >         M_fluid;
    bcInterfacePtr_Type                     M_bc;
    boost::shared_ptr< bdf_Type >           M_bdf;
    boost::shared_ptr< data_Type >          M_data;
    boost::shared_ptr< MeshData >           M_meshData;
    boost::shared_ptr< MeshPartitioner_Type > M_mesh;
    boost::shared_ptr< MapEpetra >          M_map;
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
    NonLinearAitken< fluidVector_Type >     M_generalizedAitken;
};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFluid3D()
{
    return new MultiscaleModelFluid3D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelFluid3D_H */
