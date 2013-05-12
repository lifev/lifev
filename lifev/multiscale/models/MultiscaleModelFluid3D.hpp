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

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/bc_interface/3D/bc/BCInterface3D.hpp>
#include <lifev/core/algorithm/NonLinearAitken.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>

#include <lifev/multiscale/models/MultiscaleModel.hpp>
#include <lifev/multiscale/framework/MultiscaleInterface.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelFluid3D - Multiscale model for 3D Fluid simulations
/*!
 *  @author Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleModelFluid3D class is an implementation of the multiscaleModel_Type
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 */
class MultiscaleModelFluid3D: public virtual multiscaleModel_Type,
    public virtual MultiscaleInterface
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
    void setupData ( const std::string& fileName );

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


    //! @name MultiscaleInterface Methods
    //@{

    //! Impose the flow rate on a specific interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryFlowRate ( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanNormalStress ( const multiscaleID_Type& boundaryID, const function_Type& function );

    //! Impose the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * Note: mean total normal stress cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryMeanTotalNormalStress ( const multiscaleID_Type& /*boundaryID*/, const function_Type& /*function*/ )
    {
        multiscaleErrorCheck ( ModelInterface, "Invalid interface [MeanTotalNormalStress] for model type [" + enum2String ( M_type, multiscaleModelsMap ) + "]", M_comm->MyPID() == 0 );
    }

    //! Impose the area on a specific boundary interface of the model
    /*!
     * Note: area cannot be imposed at the interfaces of this model.
     *
     * @param boundaryID ID of the boundary interface
     * @param function boundary condition function
     */
    void imposeBoundaryArea ( const multiscaleID_Type& /*boundaryID*/, const function_Type& /*function*/ )
    {
        multiscaleErrorCheck ( ModelInterface, "Invalid interface [Area] for model type [" + enum2String ( M_type, multiscaleModelsMap ) + "]", M_comm->MyPID() == 0 );
    }

    //! Get the flow rate on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return flow rate value
     */
    Real boundaryFlowRate ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->flux ( boundaryFlag ( boundaryID ) );
    }

    //! Get the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean normal stress value
     */
    Real boundaryMeanNormalStress ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->meanNormalStress ( boundaryFlag ( boundaryID ), *M_bc->handler() );
    }

    //! Get the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean total normal stress value
     */
    Real boundaryMeanTotalNormalStress ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->meanTotalNormalStress ( boundaryFlag ( boundaryID ), *M_bc->handler() );
    }

    //! Get the area on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return area value
     */
    Real boundaryArea ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->area ( boundaryFlag ( boundaryID ) );
    }

    //! Get the variation of the flow rate (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real boundaryDeltaFlowRate ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean normal stress
     */
    Real boundaryDeltaMeanNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the mean total normal stress (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean total normal stress
     */
    Real boundaryDeltaMeanTotalNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the area (on a specific boundary interface) using the linear model
    /*!
     *  Note: returns always a NaN, since this method is not used by the current interface equations.
     *
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the area
     */
    Real boundaryDeltaArea ( const multiscaleID_Type& /*boundaryID*/, bool& /*solveLinearSystem*/ )
    {
        return NaN;
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set the solution vector
    /*!
     * @param solution solution vector
     */
    void setSolution ( const fluidVectorPtr_Type& solution );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface()
    {
        return *M_bc;
    }

    //! Get the density on a specific boundary face of the model
    /*!
     * @return density value
     */
    Real boundaryDensity() const
    {
        return M_data->density();
    }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @return viscosity value
     */
    Real boundaryViscosity() const
    {
        return M_data->viscosity();
    }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @return pressure value
     */
    Real boundaryPressure ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->pressure ( boundaryFlag ( boundaryID ) );
    }

    //! Get the integral of the total pressure (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @return total pressure value
     */
    Real boundaryTotalPressure ( const multiscaleID_Type& boundaryID ) const
    {
        return M_fluid->pressure ( boundaryFlag ( boundaryID ) ) + M_fluid->kineticNormalStress ( boundaryFlag ( boundaryID ) );
    }

    //! Get the data container
    /*!
     * @return Data container
     */
    const data_Type& data() const
    {
        return *M_data;
    }

    //! Get the solution vector
    /*!
     * @return Solution vector
     */
    const fluidVector_Type& solution() const
    {
        return *M_solution;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelFluid3D ( const MultiscaleModelFluid3D& model );

    MultiscaleModelFluid3D& operator= ( const MultiscaleModelFluid3D& model );

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
    void setupGlobalData ( const std::string& fileName );

    //! Initialize the solution.
    void initializeSolution();

    //! Setup the exporter and the importer
    /*!
     * @param fileName File name of the specific model.
     */
    void setupExporterImporter ( const std::string& fileName );

    //! Setup the mesh for the fluid problem
    void setupMesh();

    //! Setup the FE space for pressure and velocity
    void setupFEspace();

    //! Setup the DOF of the model
    void setupDOF();

    //! Setup the offset for fluxes boundary conditions
    void setupBCOffset ( const bcPtr_Type& BC );

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel ( bool& solveLinearSystem );

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void imposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void resetPerturbation();

    Real bcFunctionDeltaZero ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ )
    {
        return 0.;
    }
    Real bcFunctionDeltaOne ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*id*/ )
    {
        return 1.;
    }

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
