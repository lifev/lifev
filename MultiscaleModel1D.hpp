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
 *  @brief File containing the MultiScale Model 1D
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

#ifndef MultiscaleModel1D_H
#define MultiscaleModel1D_H 1

// Jacobian coefficient approximation
#define JACOBIAN_WITH_FINITEDIFFERENCE
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
//#define JACOBIAN_WITH_FINITEDIFFERENCE_AREA
#endif

// Matlab post-processing
#define HAVE_MATLAB_POSTPROCESSING 1

// Mathcard includes
#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

#include <lifemc/lifesolver/BCInterface1D.hpp>
#include <lifemc/lifesolver/MultiscaleModel.hpp>

// LifeV includes
#include <life/lifefem/FESpace.hpp>
#ifdef HAVE_HDF5
#include <life/lifefilters/hdf5exporter.hpp>
#endif

namespace LifeV
{
namespace multiscale
{

//! MultiscaleModel1D - MultiScale model for 1D Fluid simulations
/*!
 *  @author Gilles Fourestey, Cristiano Malossi
 *
 *  The MultiscaleModel1D class is an implementation of the multiscaleModel_Type
 *  for 1D Fluid problem.
 */
class MultiscaleModel1D: public virtual multiscaleModel_Type
{
public:

    typedef OneDimensionalModel_Physics                            Physics_Type;
    typedef boost::shared_ptr< Physics_Type >                      Physics_PtrType;

    typedef OneDimensionalModel_Flux                               Flux_Type;
    typedef boost::shared_ptr< Flux_Type >                         Flux_PtrType;

    typedef OneDimensionalModel_Source                             Source_Type;
    typedef boost::shared_ptr< Source_Type >                       Source_PtrType;

    typedef OneDimensionalModel_Solver                             Solver_Type;
    typedef boost::shared_ptr< Solver_Type >                       Solver_PtrType;

    typedef Solver_Type::data_Type                                 Data_Type;
    typedef Solver_Type::mesh_Type                                 Mesh_Type;
    typedef Solver_Type::vector_Type                               Vector_Type;
    typedef Solver_Type::vectorPtr_Type                            Vector_PtrType;
    typedef Solver_Type::solution_Type                             Solution_Type;
    typedef Solver_Type::solutionPtr_Type                          Solution_PtrType;
    typedef Solver_Type::solutionConstIterator_Type                Solution_ConstIterator;
    typedef Solver_Type::linearSolver_Type                         LinearSolver_Type;

    typedef Solver_Type::FESpace_Type                              FESpace_Type;
    typedef Solver_Type::FESpacePtr_Type                           FESpace_PtrType;

    typedef BCInterface1D< Solver_Type >                           BCInterface_Type;
    typedef boost::shared_ptr< BCInterface_Type >                  BCInterface_PtrType;
    typedef OneDimensionalModel_BCHandler                          BC_Type;
    typedef boost::shared_ptr< BC_Type >                           BC_PtrType;

#ifdef HAVE_HDF5
    typedef Hdf5exporter< Mesh_Type >                              IOFile_Type;
#endif

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit MultiscaleModel1D();

    //! Destructor
    virtual ~MultiscaleModel1D() {}

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

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

    //! Setup the linear model
    void setupLinearModel();

    //! Update the linear system matrix and vectors
    void updateLinearModel();

    //! Solve the linear problem
    void solveLinearModel( bool& solveLinearSystem );

#endif

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BCInterface_Type& bcInterface() const { return *M_bc; }

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real boundaryDensity( const BCFlag& /*flag*/) const { return M_data->DensityRho(); }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real boundaryViscosity( const BCFlag& /*flag*/) const { return M_data->Viscosity(); }

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real boundaryArea( const BCFlag& flag ) const { return M_solver->BoundaryValue( *M_solution, OneD_A, flagConverter( flag ) ); }

    //! Get the flux on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flux value
     */
    Real boundaryFlowRate( const BCFlag& flag ) const { return M_solver->BoundaryValue( *M_solution, OneD_Q, flagConverter( flag ) ); }

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real boundaryPressure( const BCFlag& flag ) const { return M_solver->BoundaryValue( *M_solution, OneD_P, flagConverter( flag ) ); }

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real boundaryDynamicPressure( const BCFlag& flag ) const { return 0.5 * boundaryDensity( flag ) * std::pow( boundaryFlowRate( flag ) / boundaryArea( flag ), 2 ); }

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

    //! Get the BC handler container of the boundary conditions of the model
    /*!
     * @return BC handler
     */
    BC_Type& bc() const { return *(M_bc->handler()); }

    //! Get the data container of the 1D model.
    /*!
     * @return 1D Model data container.
     */
    Data_Type& data() const { return *M_data; }

    //! Get the Physics of the 1D model.
    /*!
     * @return 1D Model physics.
     */
    Physics_PtrType physics() const { return M_physics; }

    //! Get the Flux of the 1D model.
    /*!
     * @return 1D Model Flux.
     */
    Flux_PtrType flux() const { return M_flux; }

    //! Get the Source of the 1D model.
    /*!
     * @return 1D Model Source.
     */
    Source_PtrType source() const { return M_source; }

    //! Get the FESpace of the 1D model.
    /*!
     * @return 1D model FESpace
     */
    FESpace_PtrType FESpace() const { return M_FESpace; }

    //! Get the Solver of the 1D model.
    /*!
     * @return 1D model solver.
     */
    Solver_PtrType solver() const { return M_solver; }

    //! Get the solution container of the 1D model.
    /*!
     * @return 1D model solution.
     */
    const Solution_PtrType& solution() const { return M_solution; }

    //! Get a specific quantity of the solution container of the 1D model.
    /*!
     * @param quantity solution quantity.
     * @return 1D model solution.
     */
    const Vector_PtrType& solution( const std::string& quantity) const { return (*M_solution)[quantity]; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    MultiscaleModel1D( const MultiscaleModel1D& model );

    MultiscaleModel1D& operator=( const MultiscaleModel1D& model );

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

    //! Setup the FE space for pressure and velocity
    void setupFESpace();

    //! Initialize the solution.
    void initializeSolution();

    //! Update the solution (solution2 = solution1)
    /*!
     * @param solution1 solution to be copied.
     * @param solution2 copy of solution1.
     */
    void updateSolution( const Solution_Type& solution1, Solution_Type& solution2 );

    //! Solve the 1D hyperbolic problem
    /*!
     * @param bc BCInterface container.
     * @param solution solution container.
     * @param solverType string containing the prefix ID to display when solving the system.
     */
    void solve( BC_Type& bc, Solution_Type& solution, const std::string& solverType = " 1D-" );

    OneD_BCSide flagConverter( const BCFlag& flag ) const { return (flag == 0) ? OneD_left : OneD_right; }

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE


    //! Update linear BC
    void createLinearBC();

    //! Update linear BC
    void updateLinearBC( const Solution_Type& solution );

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void imposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void resetPerturbation();

    Real bcFunctionDelta( const Real& t );

#else

    //! Compute Jacobian coefficients using tangent problem formulation
    /*!
     * @param bcOutputSide side of the quantity to be computed.
     * @param bcOutputType type of the quantity to be computed.
     * @return Jacobian coefficient.
     */
    Real tangentProblem( const OneD_BCSide& bcOutputSide, const OneD_BC& bcOutputType );

#endif
    //@}

#ifdef HAVE_HDF5
    boost::shared_ptr< IOFile_Type >       M_exporter;
    boost::shared_ptr< IOFile_Type >       M_importer;

    boost::shared_ptr< Mesh_Type >         M_exporterMesh;
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    // Linear BC
    BC_PtrType                             M_linearBC;

    Solution_PtrType                       M_linearSolution; // Solution of the perturbed problem

    // BC perturbation
    std::vector< std::map< OneD_BCSide, std::map< OneD_BC, Real > > > M_bcPreviousTimeSteps;

    // BC Functions for tangent problem
    OneDimensionalModel_BCFunction         M_bcBaseDelta;

    Real                                   M_bcDelta;
    OneD_BC                                M_bcDeltaType;
    OneD_BCSide                            M_bcDeltaSide;
#endif

    // 1D problem
    boost::shared_ptr< Data_Type >         M_data;
    BCInterface_PtrType                    M_bc;
    Physics_PtrType                        M_physics;
    Flux_PtrType                           M_flux;
    Source_PtrType                         M_source;
    Solver_PtrType                         M_solver;

    // Linear solver
    boost::shared_ptr< LinearSolver_Type > M_linearSolver;

    // FE spaces
    boost::shared_ptr< FESpace_Type >      M_FESpace;

    Solution_PtrType                       M_solution_tn;    // Solution at time t_n
    Solution_PtrType                       M_solution;       // Solution at time t_n+1

};

//! Factory create function
inline multiscaleModel_Type* createMultiscaleModelOneDimensional()
{
    return new MultiscaleModel1D();
}

} // Namespace multiscale
} // Namespace LifeV

#endif /* MultiscaleModel1D_H */
