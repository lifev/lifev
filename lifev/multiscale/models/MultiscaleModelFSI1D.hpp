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
 *  @brief File containing the Multiscale Model 1D
 *
 *  @version 1.1
 *  @date 26-02-2010
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *
 *  @version 1.2 and subsequents
 *  @date 23-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef MultiscaleModelFSI1D_H
#define MultiscaleModelFSI1D_H 1

// Jacobian coefficient approximation
#define JACOBIAN_WITH_FINITEDIFFERENCE

// Matlab post-processing
#define HAVE_MATLAB_POSTPROCESSING 1

#include <lifev/one_d_fsi/fem/OneDFSIBCHandler.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIPhysicsLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIPhysicsNonLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIFluxLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIFluxNonLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSISourceLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSISourceNonLinear.hpp>
#include <lifev/one_d_fsi/solver/OneDFSISolver.hpp>

#include <lifev/bc_interface/1D/bc/BCInterface1D.hpp>

#include <lifev/core/fem/FESpace.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif

#include <lifev/multiscale/models/MultiscaleModel.hpp>
#include <lifev/multiscale/framework/MultiscaleInterface.hpp>

namespace LifeV
{
namespace Multiscale
{

//! MultiscaleModelFSI1D - Multiscale model for 1D Fluid simulations
/*!
 *  @author Gilles Fourestey, Cristiano Malossi
 *
 *  @see Full description of the Geometrical Multiscale Framework: \cite Malossi-Thesis
 *  @see Methodology: \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite Malossi2011Algorithms3D1DFSI \cite BlancoMalossi2012
 *  @see Applications: \cite Malossi2011Algorithms3D1DFSIAortaIliac \cite LassilaMalossi2012IdealLeftVentricle \cite BonnemainMalossi2012LVAD
 *
 *  The MultiscaleModelFSI1D class is an implementation of the multiscaleModel_Type
 *  for 1D Fluid problem.
 */
class MultiscaleModelFSI1D: public virtual multiscaleModel_Type,
    public virtual MultiscaleInterface
{
public:

    //! @name Type definitions
    //@{

    typedef OneDFSIPhysics                                         physics_Type;
    typedef std::shared_ptr< physics_Type >                      physicsPtr_Type;

    typedef OneDFSIFlux                                            flux_Type;
    typedef std::shared_ptr< flux_Type >                         fluxPtr_Type;

    typedef OneDFSISource                                          source_Type;
    typedef std::shared_ptr< source_Type >                       sourcePtr_Type;

    typedef OneDFSISolver                                          solver_Type;
    typedef std::shared_ptr< solver_Type >                       solverPtr_Type;

    typedef solver_Type::data_Type                                 data_Type;
    typedef solver_Type::mesh_Type                                 mesh_Type;
    typedef solver_Type::vector_Type                               vector_Type;
    typedef solver_Type::vectorPtr_Type                            vectorPtr_Type;
    typedef solver_Type::solution_Type                             solution_Type;
    typedef solver_Type::solutionPtr_Type                          solutionPtr_Type;
    typedef solver_Type::solutionConstIterator_Type                solutionConstIterator_Type;
    typedef solver_Type::linearSolver_Type                         linearSolver_Type;

    typedef solver_Type::feSpace_Type                              feSpace_Type;
    typedef solver_Type::feSpacePtr_Type                           feSpacePtr_Type;

    typedef OneDFSIBCHandler                                       bc_Type;
    typedef std::shared_ptr< bc_Type >                           bcPtr_Type;
    typedef BCInterface1D< bc_Type, solver_Type >                  bcInterface_Type;
    typedef std::shared_ptr< bcInterface_Type >                  bcInterfacePtr_Type;

    typedef OneDFSIFunction                                        bcFunction_Type;

    typedef OneDFSI::bcType_Type                                   bcType_Type;
    typedef OneDFSI::bcSide_Type                                   bcSide_Type;
    typedef OneDFSI::bcLine_Type                                   bcLine_Type;

#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                              IOFile_Type;
    typedef ExporterData< mesh_Type >                              IOData_Type;
#endif

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModelFSI1D();

    //! Destructor
    virtual ~MultiscaleModelFSI1D() {}

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
     * TODO The area can be imposed to the 1-D model: need to be coded
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
        return M_solver->boundaryValue ( *M_solution, OneDFSI::Q, flagConverter ( boundaryID ) );
    }

    //! Get the integral of the mean normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean normal stress value
     */
    Real boundaryMeanNormalStress ( const multiscaleID_Type& boundaryID ) const
    {
        return M_solver->boundaryValue ( *M_solution, OneDFSI::S, flagConverter ( boundaryID ) );
    }

    //! Get the integral of the mean total normal stress on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return mean total normal stress value
     */
    Real boundaryMeanTotalNormalStress ( const multiscaleID_Type& boundaryID ) const
    {
        return M_solver->boundaryValue ( *M_solution, OneDFSI::T, flagConverter ( boundaryID ) );
    }

    //! Get the area on a specific boundary interface of the model
    /*!
     * @param boundaryID ID of the boundary interface
     * @return area value
     */
    Real boundaryArea ( const multiscaleID_Type& boundaryID ) const
    {
        return M_solver->boundaryValue ( *M_solution, OneDFSI::A, flagConverter ( boundaryID ) );
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

    //! Get the variation of the integral of the total normal stress (on a specific boundary face)
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the mean total normal stress
     */
    Real boundaryDeltaMeanTotalNormalStress ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //! Get the variation of the integral of the area (on a specific boundary interface) using the linear model
    /*!
     * @param boundaryID ID of the boundary interface
     * @param solveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the area
     */
    Real boundaryDeltaArea ( const multiscaleID_Type& boundaryID, bool& solveLinearSystem );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BC handler container of the boundary conditions of the model
    /*!
     * @return BC handler
     */
    bc_Type& bc() const
    {
        return * ( M_bc->handler() );
    }

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    bcInterface_Type& bcInterface() const
    {
        return *M_bc;
    }

    //! Get the density on a specific boundary face of the model
    /*!
     * @return density value
     */
    Real boundaryDensity() const
    {
        return M_data->densityRho();
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
        return M_solver->boundaryValue ( *M_solution, OneDFSI::P, flagConverter ( boundaryID ) );
    }

    //! Get the data container of the 1D model.
    /*!
     * @return 1D Model data container.
     */
    data_Type& data() const
    {
        return *M_data;
    }

    //! Get the Physics of the 1D model.
    /*!
     * @return 1D Model physics.
     */
    physicsPtr_Type physics() const
    {
        return M_physics;
    }

    //! Get the Flux of the 1D model.
    /*!
     * @return 1D Model Flux.
     */
    fluxPtr_Type flux() const
    {
        return M_flux;
    }

    //! Get the Source of the 1D model.
    /*!
     * @return 1D Model Source.
     */
    sourcePtr_Type source() const
    {
        return M_source;
    }

    //! Get the FESpace of the 1D model.
    /*!
     * @return 1D model FESpace
     */
    feSpacePtr_Type feSpace() const
    {
        return M_feSpace;
    }

    //! Get the Solver of the 1D model.
    /*!
     * @return 1D model solver.
     */
    solverPtr_Type solver() const
    {
        return M_solver;
    }

    //! Get the solution container of the 1D model.
    /*!
     * @return 1D model solution.
     */
    const solutionPtr_Type& solution() const
    {
        return M_solution;
    }

    //! Get a specific quantity of the solution container of the 1D model.
    /*!
     * @param quantity solution quantity.
     * @return 1D model solution.
     */
    const vectorPtr_Type& solution ( const std::string& quantity) const
    {
        return (*M_solution) [quantity];
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModelFSI1D ( const MultiscaleModelFSI1D& model );

    MultiscaleModelFSI1D& operator= ( const MultiscaleModelFSI1D& model );

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

    //! Setup the FE space for pressure and velocity
    void setupFESpace();

    //! Copy the solution (solution2 = solution1)
    /*!
     * NOTE: if the size of the two vector is different due to the presence of ghost nodes
     * the method automatically add/remove them from the resulting vector.
     * @param solution1 solution to be copied.
     * @param solution2 copy of solution1.
     */
    void copySolution ( const solution_Type& solution1, solution_Type& solution2 );

    //! Update BCInterface physical solver variables
    void updateBCPhysicalSolverVariables();

    //! Solve the 1D hyperbolic problem
    /*!
     * @param bc BCInterface container.
     * @param solution solution container.
     * @param solverType string containing the prefix ID to display when solving the system.
     */
    void solve ( bc_Type& bc, solution_Type& solution, const std::string& solverType = " 1D-" );

    //! Convert the boundaryID to a bcSide type
    /*!
     * @param boundaryID ID of the boundary interface
     * @return boundary condition side.
     */
    bcSide_Type flagConverter ( const multiscaleID_Type& boundaryID ) const
    {
        return ( boundaryFlag ( boundaryID ) == 0) ? OneDFSI::left : OneDFSI::right;
    }

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

    //! Update linear BC
    void createLinearBC();

    //! Update linear BC
    void updateLinearBC ( const solution_Type& solution );

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

    Real bcFunctionDelta ( const Real& t );

#else

    //! Compute Jacobian coefficients using tangent problem formulation
    /*!
     * @param bcOutputSide side of the quantity to be computed.
     * @param bcOutputType type of the quantity to be computed.
     * @return Jacobian coefficient.
     */
    Real tangentProblem ( const bcSide_Type& bcOutputSide, const bcType_Type& bcOutputType );

    //! Solve the tangent problem formulation
    /*!
     * @param flowRate rhs for the tangent problem.
     * @param bcNode node of the quantity to be computed.
     * @return solution of the tangent problem at specific node.
     */
    Real solveTangentProblem ( solver_Type::vector_Type& rhs, const UInt& bcNode );

#endif
    //@}

#ifdef HAVE_HDF5
    std::shared_ptr< IOFile_Type >       M_exporter;
    std::shared_ptr< IOFile_Type >       M_importer;

    std::shared_ptr< mesh_Type >         M_exporterMesh;
    solutionPtr_Type                       M_exporterSolution;
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    // Linear BC
    bcPtr_Type                             M_linearBC;

    solutionPtr_Type                       M_linearSolution; // Solution of the perturbed problem

    // BC perturbation
    std::vector< std::map< bcSide_Type, std::map< bcType_Type, Real > > > M_bcPreviousTimeSteps;

    // BC Functions for tangent problem
    bcFunction_Type                        M_bcBaseDelta;

    Real                                   M_bcDelta;
    bcType_Type                            M_bcDeltaType;
    bcSide_Type                            M_bcDeltaSide;
#endif

    // 1D problem
    std::shared_ptr< data_Type >         M_data;
    bcInterfacePtr_Type                    M_bc;
    physicsPtr_Type                        M_physics;
    fluxPtr_Type                           M_flux;
    sourcePtr_Type                         M_source;
    solverPtr_Type                         M_solver;

    // Linear solver
    std::shared_ptr< linearSolver_Type > M_linearSolver;
    std::shared_ptr< linearSolver_Type > M_linearViscoelasticSolver;

    // FE spaces
    std::shared_ptr< feSpace_Type >      M_feSpace;

    solutionPtr_Type                       M_solution_tn;    // Solution at time t_n
    solutionPtr_Type                       M_solution;       // Solution at time t_n+1

};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelFSI1D()
{
    return new MultiscaleModelFSI1D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModelFSI1D_H */
