//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief MultiScale Model 1D
 *
 *  @version 1.1
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 26-02-2010
 *
 *  @version 1.2 and subsequents
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-04-2010
 */

#ifndef MS_Model_1D_H
#define MS_Model_1D_H 1

// Jacobian coefficient approximation
//#define JACOBIAN_WITH_FINITEDIFFERENCE
#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    //#define JACOBIAN_WITH_FINITEDIFFERENCE_AREA
#endif

// Matlab post-processing
#define HAVE_MATLAB_POSTPROCESSING 1

// Mathcard includes
#include <lifemc/lifesolver/BCInterface1D.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Physics_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Flux_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_Linear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Source_NonLinear.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

// LifeV includes
#include <life/lifefem/FESpace.hpp>
#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
#endif

namespace LifeV {

//! MS_Model_1D - MultiScale model for 1D Fluid simulations
/*!
 *  @author Gilles Fourestey, Cristiano Malossi
 *
 *  The MS_Model_1D class is an implementation of the MS_PhysicalModel
 *  for 1D Fluid problem.
 */
class MS_Model_1D: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                                       super;

    typedef OneDimensionalModel_Physics                            Physics_Type;
    typedef boost::shared_ptr< Physics_Type >                      Physics_PtrType;

    typedef OneDimensionalModel_Flux                               Flux_Type;
    typedef boost::shared_ptr< Flux_Type >                         Flux_PtrType;

    typedef OneDimensionalModel_Source                             Source_Type;
    typedef boost::shared_ptr< Source_Type >                       Source_PtrType;

    typedef OneDimensionalModel_Solver                             Solver_Type;
    typedef boost::shared_ptr< Solver_Type >                       Solver_PtrType;

    typedef Solver_Type::Data_Type                                 Data_Type;
    typedef Solver_Type::Mesh_Type                                 Mesh_Type;
    typedef Solver_Type::Vector_Type                               Vector_Type;
    typedef Solver_Type::Vector_PtrType                            Vector_PtrType;
    typedef Solver_Type::Solution_Type                             Solution_Type;
    typedef Solver_Type::Solution_PtrType                          Solution_PtrType;
    typedef Solver_Type::Solution_ConstIterator                    Solution_ConstIterator;
    typedef Solver_Type::LinearSolver_Type                         LinearSolver_Type;

    typedef Solver_Type::FESpace_Type                              FESpace_Type;
    typedef Solver_Type::FESpace_PtrType                           FESpace_PtrType;

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
    MS_Model_1D();

    //! Copy constructor
    /*!
     * @param OneDimensionalModel MS_Model_1D
     */
    MS_Model_1D( const MS_Model_1D& OneDimensionalModel );

    //! Destructor
    ~MS_Model_1D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param OneDimensionalModel MS_Model_1D
     * @return reference to a copy of the class
     */
    MS_Model_1D& operator = ( const MS_Model_1D& OneDimensionalModel );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the model.
    void SetupModel();

    //! Build the initial system (matrix and vectors).
    void BuildSystem();

    //! Update the system for (matrix and vectors).
    void UpdateSystem();

    //! Solve the problem.
    void SolveSystem();

    //! Save the solution
    void SaveSolution();

    //! Display some information about the model.
    void ShowMe();

    //@}


    //! @name Methods
    //@{

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE

    //! Setup the linear model
    void SetupLinearModel();

    //! Update the linear system matrix and vectors
    void UpdateLinearModel();

    //! Solve the linear problem
    void SolveLinearModel( bool& SolveLinearSystem );

#endif

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BCInterface_Type& GetBCInterface() const;

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real GetBoundaryDensity( const BCFlag& /*Flag*/) const;

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real GetBoundaryViscosity( const BCFlag& /*Flag*/) const;

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real GetBoundaryArea( const BCFlag& Flag ) const;

    //! Get the flux on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flux value
     */
    Real GetBoundaryFlowRate( const BCFlag& Flag ) const;

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real GetBoundaryPressure( const BCFlag& Flag ) const;

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real GetBoundaryDynamicPressure( const BCFlag& Flag ) const;

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param StressType Type of approximation for the stress
     * @return stress value
     */
    Real GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType = StaticPressure ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the total pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the dynamic pressure
     */
    Real GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @param StressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType = StaticPressure );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BC handler container of the boundary conditions of the model
    /*!
     * @return BC handler
     */
    BC_Type& GetBC() const;

    //! Get the data container of the 1D model.
    /*!
     * @return 1D Model data container.
     */
    Data_Type&    GetData() const;

    //! Get the Physics of the 1D model.
    /*!
     * @return 1D Model physics.
     */
    Physics_PtrType GetPhysics() const;

    //! Get the Flux of the 1D model.
    /*!
     * @return 1D Model Flux.
     */
    Flux_PtrType GetFlux() const;

    //! Get the Source of the 1D model.
    /*!
     * @return 1D Model Source.
     */
    Source_PtrType GetSource() const;

    //! Get the FESpace of the 1D model.
    /*!
     * @return 1D model FESpace
     */
    FESpace_PtrType GetFESpace() const;

    //! Get the Solver of the 1D model.
    /*!
     * @return 1D model solver.
     */
    Solver_PtrType GetSolver() const;

    //! Get the solution container of the 1D model.
    /*!
     * @return 1D model solution.
     */
    const Solution_PtrType& GetSolution() const;

    //! Get a specific quantity of the solution container of the 1D model.
    /*!
     * @param quantity solution quantity.
     * @return 1D model solution.
     */
    const Vector_PtrType& GetSolution( const std::string& quantity) const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param FileName File name of the specific model.
     */
    void SetupGlobalData( const std::string& FileName );

    //! Setup the FE space for pressure and velocity
    void SetupFESpace();

    //! Update the solution (solution2 = solution1)
    /*!
     * @param solution1 solution to be copied.
     * @param solution2 copy of solution1.
     */
    void UpdateSolution( const Solution_Type& solution1, Solution_Type& solution2 );

    //! Solve the 1D hyperbolic problem
    /*!
     * @param bc BCInterface container.
     * @param solution solution container.
     * @param solverType string containing the prefix ID to display when solving the system.
     */
    void Solve( BC_Type& bc, Solution_Type& solution, const std::string& solverType = " 1D-" );

    OneD_BCSide FlagConverter( const BCFlag& flag ) const;

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE


    //! Update linear BC
    void CreateLinearBC();

    //! Update linear BC
    void UpdateLinearBC( const Solution_Type& solution );

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void ImposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void ResetPerturbation();

    Real BCFunctionDelta( const Real& t );

#else

    //! Compute Jacobian coefficients using tangent problem formulation
    /*!
     * @param bcOutputSide side of the quantity to be computed.
     * @param bcOutputType type of the quantity to be computed.
     * @return Jacobian coefficient.
     */
    Real TangentProblem( const OneD_BCSide& bcOutputSide, const OneD_BC& bcOutputType );

#endif
    //@}

#ifdef HAVE_HDF5
    boost::shared_ptr< IOFile_Type >       M_Exporter;
    boost::shared_ptr< Mesh_Type >         M_ExporterMesh;
#endif

#ifdef JACOBIAN_WITH_FINITEDIFFERENCE
    // Linear BC
    BC_PtrType                             M_LinearBC;

    Solution_PtrType                       M_LinearSolution; // Solution of the perturbed problem

    // BC perturbation
    std::vector< std::map< OneD_BCSide, std::map< OneD_BC, Real > > > M_BCPreviousTimeSteps;

    // BC Functions for tangent problem
    OneDimensionalModel_BCFunction         M_BCBaseDelta;

    Real                                   M_BCDelta;
    OneD_BC                                M_BCDeltaType;
    OneD_BCSide                            M_BCDeltaSide;
#endif

    // 1D problem
    boost::shared_ptr< Data_Type >         M_Data;
    BCInterface_PtrType                    M_BC;
    Physics_PtrType                        M_Physics;
    Flux_PtrType                           M_Flux;
    Source_PtrType                         M_Source;
    Solver_PtrType                         M_Solver;

    // Linear solver
    boost::shared_ptr< LinearSolver_Type > M_LinearSolver;

    // FE spaces
    boost::shared_ptr< FESpace_Type >      M_FESpace;

    Solution_PtrType                       M_Solution_tn;    // Solution at time t_n
    Solution_PtrType                       M_Solution;       // Solution at time t_n+1

};

//! Factory create function
inline MS_PhysicalModel* MS_createOneDimensional()
{
    return new MS_Model_1D();
}

} // Namespace LifeV

#endif /* MS_Model_1D_H */
