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
    @brief This file contains an Oseen equation solver class.

    @author Gilles Fourestey <gilles.fourestey@imag.fr>
            Simone Deparis   <simone.deparis@epfl.ch>
    @contributor Zhen Wang <zhen.wang@emory.edu>

    @date 01-06-2007

    This file contains an Oseen equation solver class.
    The resulting linear systems are solved by GMRES on the full
    matrix ( u and p coupled ).

 */


#ifndef OSEENSOLVER_H
#define OSEENSOLVER_H 1

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/util/LifeChrono.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/core/fem/PostProcessingBoundary.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/navier_stokes/solver/StabilizationIP.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>

#include <boost/shared_ptr.hpp>

#include <list>

namespace LifeV
{
//! @class Oseen
/*!
    @brief This class contains an Oseen equation solver.

    @author Gilles Fourestey <gilles.fourestey@imag.fr>
            Simone Deparis   <simone.deparis@epfl.ch>
    @contributor Zhen Wang <zhen.wang@emory.edu>

 */

template< typename MeshType, typename SolverType = LifeV::SolverAztecOO >
class OseenSolver
{

public:

    //! @name Public Types
    //@{

    typedef MeshType                                    mesh_Type;
    typedef SolverType                                  linearSolver_Type;
    typedef boost::shared_ptr<linearSolver_Type>        linearSolverPtr_Type;
    typedef OseenData                                   data_Type;
    typedef boost::shared_ptr< data_Type >              dataPtr_Type;

    typedef boost::function < Real ( const Real& t, const Real& x, const Real& y,
                                     const Real& z, const ID& i ) > function_Type;

    typedef boost::function < Real ( const Real& t, const Real& x, const Real& y,
                                     const Real& z, const ID& i ) > source_Type;

    typedef BCHandler                                   bcHandler_Type;
    typedef boost::shared_ptr<bcHandler_Type>           bcHandlerPtr_Type;

#ifdef HAVE_NS_PREC
    typedef MatrixEpetraStructured<Real>                matrix_Type;
#else
    typedef typename linearSolver_Type::matrix_type     matrix_Type;
#endif
    typedef boost::shared_ptr<matrix_Type>              matrixPtr_Type;
    typedef typename linearSolver_Type::vector_type     vector_Type;
    typedef boost::shared_ptr<vector_Type>              vectorPtr_Type;

    typedef vector_Type                                 solution_Type;
    typedef boost::shared_ptr<solution_Type>            solutionPtr_Type;

    typedef typename linearSolver_Type::prec_raw_type   preconditioner_Type;
    typedef typename linearSolver_Type::prec_type       preconditionerPtr_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    OseenSolver();

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */

    OseenSolver ( boost::shared_ptr<data_Type>    dataType,
                  FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                  FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                  boost::shared_ptr<Epetra_Comm>& communicator,
                  const Int                       lagrangeMultiplier = 0 );

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param monolithicMap MapEpetra class
        @param offset
     */
    OseenSolver ( boost::shared_ptr<data_Type>    dataType,
                  FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                  FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                  boost::shared_ptr<Epetra_Comm>& communicator,
                  const MapEpetra                 monolithicMap,
                  const UInt                      offset = 0 );

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param lagrangeMultipliers (lagrange multipliers for the flux problem with rufaec flag)
        @param communicator MPI communicator
     */
    OseenSolver ( boost::shared_ptr<data_Type>    dataType,
                  FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                  FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                  const std::vector<Int>&         lagrangeMultipliers,
                  boost::shared_ptr<Epetra_Comm>& communicator );

    //! virtual destructor
    virtual ~OseenSolver();

    //@}

    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    virtual void setUp ( const GetPot& dataFile );

    //! Initialize with velocityFunction and pressureFunction
    /*!
        @param velocityFunction
        @param pressureFunction
     */
    void initialize ( const function_Type& velocityFunction, const function_Type& pressureFunction );

    //! Initialize with velocityInitialGuess and pressureInitialGuess
    /*!
        @param velocityInitialGuess
        @param pressureInitialGuess
     */
    void initialize ( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess );

    //! Initialize with velocityAndPressure
    /*!
        @param velocityAndPressure
     */
    void initialize ( const vector_Type& velocityAndPressure );

    //! Build linear system.
    virtual void buildSystem();

    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
     */
    virtual void updateSystem ( const Real         alpha,
                                const vector_Type& betaVector,
                                const vector_Type& sourceVector );

    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
        @param matrix
        @param un
     */
    virtual void updateSystem ( const Real         alpha,
                                const vector_Type& betaVector,
                                const vector_Type& sourceVector,
                                matrixPtr_Type     matrix,
                                const vector_Type&     un );

    //! Update stabilization term
    /*!
        @param matrixFull
     */
    void updateStabilization ( matrix_Type& matrixFull );

    //! Update the right hand side
    /*!
        @param rightHandSide right hand side
     */
    virtual void updateRightHandSide ( const vector_Type& rightHandSide )
    {
        M_rightHandSideNoBC.reset (new vector_Type (rightHandSide) );
        M_rightHandSideNoBC->globalAssemble();
    }

    //! Update the source term
    /*!
     @param sourceTerm
    */
    void updateSourceTerm (const source_Type& source );

    //! Update convective term, boundary condition and solve the linearized ns system
    /*!
        @param bcHandler BC handler
     */
    virtual void iterate ( bcHandler_Type& bcHandler );

    //! Reduce the local solution in global vectors
    /*!
        @param velocity
        @param pressure
     */
    void reduceSolution ( Vector& velocity, Vector& pressure );

    //! Reduce the residual
    /*!
        @param residual
     */
    void reduceResidual ( Vector& residual );

    //! Set a block preconditioner
    /*!
        @blockPrecconditioner Block preconditioner
     */
    void setBlockPreconditioner ( matrixPtr_Type blockPreconditioner );

    //! Update and return the coefficient matrix
    /*!
        @param matrixFull The coefficient matrix
     */
    void getFluidMatrix ( matrix_Type& matrixFull );

    //! Set up post processing
    void setupPostProc( );

    //! Compute area on a boundary face with given flag
    /*!
        @param  flag
        @return area
     */
    Real area ( const markerID_Type& flag );

    //! Compute the outgoing normal of a boundary face with given flag
    /*!
        @param  flag
        @return boundary normal
     */
    Vector normal ( const markerID_Type& flag );

    //! Compute the geometric center of a boundary face with given flag
    /*!
        @param  flag
        @return geometric center
     */
    Vector geometricCenter ( const markerID_Type& flag );

    //! Compute flux on a boundary face with given flag and a given solution
    /*!
        @param  flag
        @param  solution
        @return flux
     */
    Real flux ( const markerID_Type& flag, const vector_Type& solution );

    //! Compute flux on a boundary face with given flag
    /*!
        @param flag
        @return flux
     */
    Real flux ( const markerID_Type& flag );

    //! Compute the kinetic normal stress (i.e., the normal stress due to the kinetic energy) on a boundary face with a given flag and a given solution
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  This method computes the following quantity:
     *
     *  \f[
     *  \mathcal{K} = \frac{1}{2}\rho_\textrm{F}\frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}\left({\mathbf{u}}_\textrm{F} \mathbf{\cdot} {\mathbf{n}}_\textrm{F}\right)^2  \textrm{d} \Gamma
     *  \f]
     *
     *  @param flag boundary flag
     *  @param solution problem solution
     *  @return kinetic normal stress
     */
    Real kineticNormalStress ( const markerID_Type& flag, const vector_Type& solution );

    //! Compute the kinetic normal stress (i.e., the normal stress due to the kinetic energy) on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  This method computes the following quantity:
     *
     *  \f[
     *  \mathcal{K} = \frac{1}{2}\rho_\textrm{F}\frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}\left({\mathbf{u}}_\textrm{F} \mathbf{\cdot} {\mathbf{n}}_\textrm{F}\right)^2  \textrm{d} \Gamma
     *  \f]
     *
     *  @param flag boundary flag
     *  @return kinetic normal stress
     */
    Real kineticNormalStress ( const markerID_Type& flag );

    //! Compute average pressure on a boundary face with given flag and a given solution
    /*!
        @param  flag
        @param  solution
        @return average pressure
     */
    Real pressure ( const markerID_Type& flag, const vector_Type& solution );

    //! Compute average pressure on a boundary face with given flag
    /*!
        @param flag
        @return average pressure
     */
    Real pressure ( const markerID_Type& flag );

    //! Compute the mean normal stress on a boundary face with a given flag and a given solution
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  The mean normal stress is defined as the average of the normal component of the traction vector.
     *
     *  \f[
     *  \mathcal{S}^{3-D}_j = \frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}
     *  \left(\sigma_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right) \mathbf{\cdot} \mathbf{n}_\textrm{F} \: \textrm{d} \Gamma
     *  \f]
     *
     *  @param flag Flag of the boundary face
     *  @param bcHandler BChandler containing the boundary conditions of the problem.
     *  @param solution Vector containing the solution of the problem
     *                   (and also the Lagrange multipliers at the end).
     *  @return mean normal stress
     */
    Real meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution );

    //! Compute the mean normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  The mean normal stress is defined as the average of the normal component of the traction vector.
     *
     *  \f[
     *  \mathcal{S}^{3-D}_j = \frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}
     *  \left(\sigma_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right) \mathbf{\cdot} \mathbf{n}_\textrm{F} \: \textrm{d} \Gamma
     *  \f]
     *
     *  TODO The current version returns the exact mean normal stress if a flow rate boundary condition is imposed on the chosen boundary face.
     *  On the contrary, if other boundary conditions are applied, the mean normal stress is approximated with the mean pressure, which is a
     *  reasonable approximation for several applications.
     *
     *  @param flag Flag of the boundary face
     *  @param bcHandler BChandler containing the boundary conditions of the problem.
     *  @return mean normal stress
     */
    Real meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler );

    //! Compute the mean total normal stress on a boundary face with a given flag and a given solution
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  The mean total normal stress is defined as the average of the normal component of the traction vector minus the kinetic contribution.
     *
     *  \f[
     *  \mathcal{T}^{3-D}_j = \frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}
     *  \left(\sigma_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right) \mathbf{\cdot} \mathbf{n}_\textrm{F} \: \textrm{d} \Gamma
     *  - \frac{1}{2}\rho_\textrm{F}\frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}\left(\mathbf{u}_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right)^2 \: \textrm{d} \Gamma
     *  \f]
     *
     *  TODO The current version returns the exact mean normal stress if a flow rate boundary condition is imposed on the chosen boundary face.
     *  On the contrary, if other boundary conditions are applied, the mean normal stress is approximated with the mean pressure, which is a
     *  reasonable approximation for several applications.
     *
     *  @param flag Flag of the boundary face
     *  @param bcHandler BChandler containing the boundary conditions of the problem.
     *  @param solution Vector containing the solution of the problem
     *                   (and also the Lagrange multipliers at the end).
     *  @return mean total normal stress
     */
    Real meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution );

    //! Compute the mean total normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  The mean total normal stress is defined as the average of the normal component of the traction vector minus the kinetic contribution.
     *
     *  \f[
     *  \mathcal{T}^{3-D}_j = \frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}
     *  \left(\sigma_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right) \mathbf{\cdot} \mathbf{n}_\textrm{F} \: \textrm{d} \Gamma
     *  - \frac{1}{2}\rho_\textrm{F}\frac{1}{\left|\Gamma^t_{\textrm{F},j}\right|}\displaystyle\int_{\Gamma^t_{\textrm{F},j}}\left(\mathbf{u}_\textrm{F} \mathbf{\cdot} \mathbf{n}_\textrm{F}\right)^2 \: \textrm{d} \Gamma
     *  \f]
     *
     *  @param flag Flag of the boundary face
     *  @param bcHandler BChandler containing the boundary conditions of the problem.
     *  @return mean total normal stress
     */
    Real meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary
    /*!
        @param flag      Flag of the boundary face associated with the flux
                         and the Lagrange multiplier we want.
        @param bcHandler BChandler containing the boundary conditions of the problem.
        @return          Lagrange multiplier
     */
    Real lagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary
    /*!
        @param flag      Flag of the boundary face associated
                         with the flux and the Lagrange multiplier we want.
        @param bcHandler BChandler containing the boundary conditions of the problem.
        @param solution  Vector containing the solution of the problem
                         (and also the Lagrange multipliers at the end).
        @return          Lagrange multiplier
     */
    Real lagrangeMultiplier ( const markerID_Type&  flag,
                              bcHandler_Type& bcHandler,
                              const vector_Type& solution );

    //! Reset the preconditioner.
    /*!
        @param reset Reset preconditioner.
     */
    void resetPreconditioner ( bool reset = true )
    {
        if ( reset )
        {
            M_linearSolver->resetPreconditioner();
        }
    }

    //! Reset stabilization matrix at the same time as the preconditioner
    void resetStabilization()
    {
        M_resetStabilization = true;
    }

    //! Update
    void updateUn()
    {
        *M_un = *M_solution;
    }

    //! Update for the monolithic
    void updateUn ( const vector_Type& solution )
    {
        *M_un = solution;
    }

    //! Display general information about the content of the class
    /*!
        @param output specify the output format (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Set
    /*!
        @param recomputeMatrix
     */
    void setRecomputeMatrix ( const bool& recomputeMatrix )
    {
        M_recomputeMatrix = recomputeMatrix;
    }

    //! set the source term functor
    /*!
        @param source
     */
    void setSourceTerm ( source_Type source )
    {
        M_source = source;
    }

    //! Set the tolerance and the maximum number of iterations of the linear solver
    /*!
        @param tolerance Tolerance
        @param maxIteration maximum number of iterations
     */
    void setTolMaxIteration ( const Real& tolerance, const Int& maxIteration = -1 );

    //@}

    //! @name Get Methods
    //@{

    //! Return the data container of the fluid
    /*!
        @return data container of the fluid
     */
    const dataPtr_Type& data() const
    {
        return M_oseenData;
    }

    //! Return the density of the fluid
    /*!
        @return Density of the fluid
     */
    const Real& density() const
    {
        return M_oseenData->density();
    }

    //! Return the viscosity of the fluid
    /*!
        @return Viscosity of the fluid
     */
    const Real& viscosity() const
    {
        return M_oseenData->viscosity();
    }

    //! Return the local solution vector
    /*!
        @return vectorPtr_Type Solution vector
     */
    const vectorPtr_Type& solution() const
    {
        return M_solution;
    }

    //! Return the local residual vector
    /*!
        @return Residual vector
     */
    const vector_Type& residual() const
    {
        return *M_residual;
    }

    //! Return velocity FE space
    /*!
        @return velocity FE space
     */
    FESpace<mesh_Type, MapEpetra>& velocityFESpace()
    {
        return M_velocityFESpace;
    }

    const FESpace<mesh_Type, MapEpetra>& velocityFESpace() const
    {
        return M_velocityFESpace;
    }

    //! Return pressure FE space
    /*!
        @return pressure FE space
     */
    FESpace<mesh_Type, MapEpetra>& pressureFESpace()
    {
        return M_pressureFESpace;
    }

    const FESpace<mesh_Type, MapEpetra>& pressureFESpace() const
    {
        return M_pressureFESpace;
    }

    //! Get the source term
    /*!
        @return Source term
     */
    const source_Type& sourceTerm() const
    {
        return M_source;
    }

    //! Returns the post processing structure
    /*!
        @return Post processing
    */
    PostProcessingBoundary<mesh_Type>& postProcessing()
    {
        return *M_postProcessing;
    }

    const PostProcessingBoundary<mesh_Type>& postProcessing() const
    {
        return *M_postProcessing;
    }

    //! Return MapEpetra.
    /*!
        @return MapEpetra
     */
    const MapEpetra& getMap() const
    {
        return M_localMap;
    }

    //! Return Epetra communicator
    /*!
        @return Epetra communicator
     */
    const boost::shared_ptr<Epetra_Comm>& comm() const
    {
        return M_Displayer.comm();
    }

    //! Return displayer
    /*!
        @return
     */
    const Displayer& getDisplayer() const
    {
        return M_Displayer;
    }

    //! Return
    /*!
        @return recomputeMatrix
     */
    const bool& recomputeMatrix() const
    {
        return M_recomputeMatrix;
    }

    //! Return matrix without boundary conditions
    /*!
        @return Matrix without boundary conditions
     */
    matrix_Type& matrixNoBC()
    {
        return *M_matrixNoBC;
    }

    const matrix_Type& matrixNoBC() const
    {
        return *M_matrixNoBC;
    }

    //! Return mass matrix
    /*!
        @return Mass matrix
     */
    matrix_Type& matrixMass()
    {
        return *M_velocityMatrixMass;
    }

    const matrix_Type& matrixMass() const
    {
        return *M_velocityMatrixMass;
    }

    const matrixPtr_Type matrixMassPtr() const
    {
        return M_velocityMatrixMass;
    }

    //@}

    //@{ unused methods

    //! Set up post processing structures
    void postProcessingSetArea();

    //! Set up post processing
    void postProcessingSetNormal();

    //! Set up post processing
    void postProcessingSetPhi();

    //! Return a bool value if using diagonal block preconditioner
    bool getIsDiagonalBlockPreconditioner()
    {
        return M_isDiagonalBlockPreconditioner;
    }

    const bool& getIsDiagonalBlockPreconditioner() const
    {
        return M_isDiagonalBlockPreconditioner;
    }

    //@}

    //! Return a shared pointer to the preconditioner (of type derived from EpetraPreconditioner)
    preconditionerPtr_Type& preconditioner()
    {
        return M_linearSolver.preconditioner();
    }

protected:

    //! @name Constructor
    //@{

    //! Empty copy constructor
    OseenSolver ( const OseenSolver& oseen);

    //@}

    //! @name Private Methods
    //@{

    //! Removes mean of component of vector x
    /*!
        @param x
        @return
     */
    Real removeMean ( vector_Type& x );

    //! Apply boundary conditions.
    /*!
        @param matrix
        @param rightHandSide
        @param bcHandler
     */
    void applyBoundaryConditions ( matrix_Type&        matrix,
                                   vector_Type&        rightHandSide,
                                   bcHandler_Type& bcHandler );

    //! Echo message.
    /*!
        @param message
     */
    void echo ( std::string message );

    //! Return the dim of velocity FE space
    const UInt& dimVelocity() const
    {
        return M_velocityFESpace.dim();
    }

    //! Return the dim of pressure FE space
    const UInt& dimPressure() const
    {
        return M_pressureFESpace.dim();
    }

    //@}

    //private members

    //! data for Navier-Stokes solvers
    dataPtr_Type                   M_oseenData;

    // FE spaces
    FESpace<mesh_Type, MapEpetra>& M_velocityFESpace;
    FESpace<mesh_Type, MapEpetra>& M_pressureFESpace;

    //! MPI communicator
    Displayer                      M_Displayer;

    MapEpetra                      M_localMap;

    //! mass matrix
    matrixPtr_Type                 M_velocityMatrixMass;

    //! mass matrix
    matrixPtr_Type                 M_pressureMatrixMass;

    //! Stokes matrix: nu*stiff
    matrixPtr_Type                 M_matrixStokes;

    //! matrix to be solved
    //    matrixPtr_Type               M_matrixFull;

    //! matrix without boundary conditions
    matrixPtr_Type                 M_matrixNoBC;

    //! stabilization matrix
    matrixPtr_Type                 M_matrixStabilization;

    //! source term for Navier-Stokes equations
    source_Type                    M_source;

    //! Right hand side for the velocity component
    vectorPtr_Type                 M_rightHandSideNoBC;

    //! Global solution
    vectorPtr_Type                 M_solution;

    //! residual
    vectorPtr_Type                    M_residual;

    linearSolverPtr_Type              M_linearSolver;

    bool                           M_steady;

    //! Postprocessing class
    boost::shared_ptr<PostProcessingBoundary<mesh_Type> > M_postProcessing;

    //! Stabilization
    bool                           M_stabilization;
    bool                           M_reuseStabilization;
    bool                           M_resetStabilization;
    Int                            M_iterReuseStabilization;

    details::StabilizationIP<mesh_Type, DOF> M_ipStabilization;
    Real                           M_gammaBeta;
    Real                           M_gammaDiv;
    Real                           M_gammaPress;

    const function_Type*                M_betaFunction;

    bool                           M_divBetaUv;

    bool                           M_stiffStrain;

    //
    Real                           M_diagonalize;

    UInt                           M_count;

    bool                           M_recomputeMatrix;

    bool                           M_isDiagonalBlockPreconditioner;

    //! Elementary matrices and vectors
    MatrixElemental                        M_elementMatrixStiff;      // velocity Stokes
    MatrixElemental                        M_elementMatrixMass;       // velocity mass
    MatrixElemental                        M_elementMatrixPreconditioner;          // (p,q) bloc for preconditioners
    MatrixElemental                        M_elementMatrixDivergence;
    MatrixElemental                        M_elementMatrixGradient;
    VectorElemental                        M_elementRightHandSide;           // Elementary right hand side
    matrixPtr_Type                         M_blockPreconditioner;
    VectorElemental                        M_wLoc;
    VectorElemental                        M_uLoc;
    boost::shared_ptr<vector_Type> M_un;

}; // class OseenSolver



// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              const Int                       lagrangeMultiplier ) :
    M_oseenData       ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
    M_velocityMatrixMass     ( ),
    M_pressureMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( new vector_Type ( M_localMap ) ),
    M_residual               ( new vector_Type (M_localMap ) ),
    M_linearSolver           ( new linearSolver_Type (communicator) ),
    M_steady                 ( ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
                                                                       &M_velocityFESpace.feBd(),
                                                                       &M_velocityFESpace.dof(),
                                                                       &M_pressureFESpace.feBd(),
                                                                       &M_pressureFESpace.dof(),
                                                                       M_localMap ) ),
    M_stabilization          ( false ),
    M_reuseStabilization     ( false ),
    M_resetStabilization     ( false ),
    M_iterReuseStabilization ( -1 ),
    //        M_ipStabilization        ( M_velocityFESpace.mesh(),
    //                                   M_velocityFESpace.dof(),
    //                                   M_velocityFESpace.refFE(),
    //                                   M_velocityFESpace.feBd(),
    //                                   M_velocityFESpace.qr(),
    //                                   0., 0., 0.,
    //                                   M_oseenData->viscosity() ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_un                     ( new vector_Type (M_localMap) )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              MapEpetra                       monolithicMap,
              UInt                            /*offset*/ ) :
    M_oseenData              ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( monolithicMap ),
    M_velocityMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( ),
    M_residual               (  ),
    M_linearSolver           ( ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
                                                                      &M_velocityFESpace.feBd(),
                                                                      &M_velocityFESpace.dof(),
                                                                      &M_pressureFESpace.feBd(),
                                                                      &M_pressureFESpace.dof(),
                                                                      M_localMap ) ),
    M_stabilization          ( false ),
    M_reuseStabilization     ( false ),
    M_resetStabilization     ( false ),
    M_iterReuseStabilization ( -1 ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner                 ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
    M_un                     ( /*new vector_Type(M_localMap)*/ )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              const std::vector<Int>&         lagrangeMultipliers,
              boost::shared_ptr<Epetra_Comm>& communicator ) :
    M_oseenData       ( dataType ),
    M_velocityFESpace        ( velocityFESpace ),
    M_pressureFESpace        ( pressureFESpace ),
    M_Displayer              ( communicator ),
    M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultipliers ),
    M_velocityMatrixMass     ( ),
    M_matrixStokes           ( ),
    M_matrixNoBC             ( ),
    M_matrixStabilization    ( ),
    M_rightHandSideNoBC      ( ),
    M_solution               ( new vector_Type ( M_localMap ) ),
    M_residual               (  ),
    M_linearSolver           ( new linearSolver_Type (communicator) ),
    M_postProcessing         ( new PostProcessingBoundary<mesh_Type> (M_velocityFESpace.mesh(),
                                                                      &M_velocityFESpace.feBd(),
                                                                      &M_velocityFESpace.dof(),
                                                                      &M_pressureFESpace.feBd(),
                                                                      &M_pressureFESpace.dof(),
                                                                      M_localMap ) ),
    M_stabilization          ( false ),
    M_reuseStabilization     ( false ),
    M_resetStabilization     ( false ),
    M_iterReuseStabilization ( -1 ),
    M_betaFunction           ( 0 ),
    M_divBetaUv              ( false ),
    M_stiffStrain            ( false ),
    M_diagonalize            ( false ),
    M_count                  ( 0 ),
    M_recomputeMatrix        ( false ),
    M_elementMatrixStiff     ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixMass      ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), M_velocityFESpace.fieldDim() ),
    M_elementMatrixPreconditioner ( M_pressureFESpace.fe().nbFEDof(), 1, 1 ),
    M_elementMatrixDivergence ( M_pressureFESpace.fe().nbFEDof(), 1, 0,
                                M_velocityFESpace.fe().nbFEDof(), 0, M_velocityFESpace.fieldDim() ),
    M_elementMatrixGradient  ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim(), 0,
                               M_pressureFESpace.fe().nbFEDof(), 0, 1 ),
    M_elementRightHandSide   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_blockPreconditioner    ( ),
    M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
    M_un                     ( new vector_Type (M_localMap) )
{
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType>
OseenSolver<MeshType, SolverType>::
~OseenSolver()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setUp ( const GetPot& dataFile )
{
    if (M_linearSolver.get() )
    {
        M_linearSolver->setupPreconditioner ( dataFile, "fluid/prec" );
        M_linearSolver->setDataFromGetPot ( dataFile, "fluid/solver" );
    }

    M_stabilization = dataFile ( "fluid/ipstab/use",  (&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ) );
    M_steady        = dataFile ( "fluid/miscellaneous/steady", 0 );
    if (M_stabilization)
    {
        M_gammaBeta     = dataFile ( "fluid/ipstab/gammaBeta",  0. );
        M_gammaDiv      = dataFile ( "fluid/ipstab/gammaDiv",   0. );
        M_gammaPress    = dataFile ( "fluid/ipstab/gammaPress", 0. );
        M_reuseStabilization     = dataFile ( "fluid/ipstab/reuse", false );
        if (M_linearSolver.get() )
            M_iterReuseStabilization = dataFile ( "fluid/ipstab/max_iter_reuse",
                                                  static_cast<Int> ( M_linearSolver->maxNumIterations() * 8. / 10. ) );
    }
    // Energetic stabilization term
    M_divBetaUv   = dataFile ( "fluid/space_discretization/div_beta_u_v", false);
    // Enable grad( u )^T in stress tensor
    M_stiffStrain = dataFile ( "fluid/space_discretization/stiff_strain", false);
    M_diagonalize = dataFile ( "fluid/space_discretization/diagonalize", 1. );
    M_isDiagonalBlockPreconditioner = dataFile ( "fluid/diagonalBlockPrec", false );

    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

    M_ipStabilization.setGammaBeta ( M_gammaBeta );
    M_ipStabilization.setGammaDiv  ( M_gammaDiv );
    M_ipStabilization.setGammaPress ( M_gammaPress );
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const function_Type& velocityFunction, const function_Type& pressureFunction )
{
    vector_Type velocityInitialGuess ( M_velocityFESpace.map() );
    M_velocityFESpace.interpolate ( velocityFunction,
                                    velocityInitialGuess,
                                    M_oseenData->dataTime()->time() );

    vector_Type pressureInitialGuess ( M_pressureFESpace.map() );
    M_pressureFESpace.interpolate ( pressureFunction,
                                    pressureInitialGuess,
                                    M_oseenData->dataTime()->time() );

    initialize ( velocityInitialGuess, pressureInitialGuess );
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess )
{

    *M_solution = velocityInitialGuess;
    *M_un = velocityInitialGuess;
    M_solution->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );
    M_un->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
initialize ( const vector_Type& velocityAndPressure )
{

    *M_un = velocityAndPressure;
    if ( M_solution.get() )
    {
        *M_solution = velocityAndPressure;
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::buildSystem()
{
    M_velocityMatrixMass.reset  ( new matrix_Type ( M_localMap ) );
    M_matrixStokes.reset ( new matrix_Type ( M_localMap ) );

    M_Displayer.leaderPrint ( "  F-  Computing constant matrices ...          " );

    LifeChrono chrono;

    LifeChrono chronoDer;
    LifeChrono chronoStiff;
    LifeChrono chronoMass;
    LifeChrono chronoGrad;

    LifeChrono chronoStiffAssemble;
    LifeChrono chronoMassAssemble;
    LifeChrono chronoGradAssemble;
    LifeChrono chronoDivAssemble;
    LifeChrono chronoStab;
    LifeChrono chronoZero;

    // Number of velocity components
    UInt numVelocityComponent = M_velocityFESpace.fieldDim();

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
    //    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner.reset ( new matrix_Type ( M_localMap ) );
    }
    chrono.start();

    for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); iElement++ )
    {
        chronoDer.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update ( M_velocityFESpace.mesh()->element ( iElement ) );
        // just to provide the id number in the assem_mat_mixed
        // M_pressureFESpace.fe().updateFirstDeriv( M_velocityFESpace.mesh()->element( iElement ) );
        M_velocityFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );

        chronoDer.stop();

        chronoZero.start();
        M_elementMatrixStiff.zero();
        M_elementMatrixMass.zero();
        M_elementMatrixPreconditioner.zero();
        M_elementMatrixDivergence.zero();
        M_elementMatrixGradient.zero();
        chronoZero.stop();

        // stiffness matrix
        chronoStiff.start();
        if ( M_stiffStrain )
            stiff_strain ( 2.0 * M_oseenData->viscosity(),
                           M_elementMatrixStiff,
                           M_velocityFESpace.fe() );
        else
            stiff ( M_oseenData->viscosity(),
                    M_elementMatrixStiff,
                    M_velocityFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
        //stiff_div( 0.5*M_velocityFESpace.fe().diameter(), M_elementMatrixStiff, M_velocityFESpace.fe() );
        chronoStiff.stop();

        // mass matrix
        if ( !M_steady )
        {
            chronoMass.start();
            mass ( M_oseenData->density(),
                   M_elementMatrixMass,
                   M_velocityFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
            chronoMass.stop();
        }

        for ( UInt iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
        {
            // stiffness matrix
            chronoStiffAssemble.start();
            if ( M_isDiagonalBlockPreconditioner == true )
            {
                assembleMatrix ( *M_blockPreconditioner,
                                 M_elementMatrixStiff,
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.dof(),
                                 M_velocityFESpace.dof(),
                                 iComponent, iComponent,
                                 iComponent * velocityTotalDof, iComponent * velocityTotalDof);
            }
            else
            {
                if ( M_stiffStrain ) // sigma = 0.5 * mu (grad( u ) + grad ( u )^T)
                {
                    for ( UInt jComp = 0; jComp < numVelocityComponent; jComp++ )
                    {
                        assembleMatrix ( *M_matrixStokes,
                                         M_elementMatrixStiff,
                                         M_velocityFESpace.fe(),
                                         M_velocityFESpace.fe(),
                                         M_velocityFESpace.dof(),
                                         M_velocityFESpace.dof(),
                                         iComponent, jComp,
                                         iComponent * velocityTotalDof, jComp * velocityTotalDof);

                    }
                }
                else // sigma = mu grad( u )
                {
                    assembleMatrix ( *M_matrixStokes,
                                     M_elementMatrixStiff,
                                     M_velocityFESpace.fe(),
                                     M_velocityFESpace.fe(),
                                     M_velocityFESpace.dof(),
                                     M_velocityFESpace.dof(),
                                     iComponent, iComponent,
                                     iComponent * velocityTotalDof, iComponent * velocityTotalDof);
                }
            }
            chronoStiffAssemble.stop();

            // mass matrix
            if ( !M_steady )
            {
                chronoMassAssemble.start();
                assembleMatrix ( *M_velocityMatrixMass,
                                 M_elementMatrixMass,
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.dof(),
                                 M_velocityFESpace.dof(),
                                 iComponent, iComponent,
                                 iComponent * velocityTotalDof, iComponent * velocityTotalDof);
                chronoMassAssemble.stop();
            }

            // divergence
            chronoGrad.start();
            grad ( iComponent, 1.0,
                   M_elementMatrixGradient,
                   M_velocityFESpace.fe(),
                   M_pressureFESpace.fe(),
                   iComponent, 0 );
            chronoGrad.stop();

            chronoGradAssemble.start();
            assembleMatrix ( *M_matrixStokes,
                             M_elementMatrixGradient,
                             M_velocityFESpace.fe(),
                             M_pressureFESpace.fe(),
                             M_velocityFESpace.dof(),
                             M_pressureFESpace.dof(),
                             iComponent, 0,
                             iComponent * velocityTotalDof, numVelocityComponent * velocityTotalDof );
            chronoGradAssemble.stop();

            chronoDivAssemble.start();
            assembleTransposeMatrix ( *M_matrixStokes,
                                      -1.,
                                      M_elementMatrixGradient,
                                      M_pressureFESpace.fe(),
                                      M_velocityFESpace.fe(),
                                      M_pressureFESpace.dof(),
                                      M_velocityFESpace.dof(),
                                      0 , iComponent,
                                      numVelocityComponent * velocityTotalDof, iComponent * velocityTotalDof );
            chronoDivAssemble.stop();
        }
    }

    //    for (UInt ii = M_velocityFESpace.fieldDim()*dimVelocity(); ii < M_velocityFESpace.fieldDim()*dimVelocity() + dimPressure(); ++ii)
    //  M_matrixStokes->set_mat_inc( ii ,ii, 0. ); not scalable!!!

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        M_blockPreconditioner->globalAssemble();
        *M_matrixStokes += *M_blockPreconditioner;
    }
    comm()->Barrier();

    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    M_Displayer.leaderPrint ( "  F-  Finalizing the matrices ...              " );

    chrono.start();

    M_matrixStokes->globalAssemble();
    M_velocityMatrixMass->globalAssemble();

    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    if ( false )
        std::cout << " partial times:  \n"
                  << " Der            " << chronoDer.diffCumul() << " s.\n"
                  << " Zero           " << chronoZero.diffCumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
                  << " Mass           " << chronoMass.diffCumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
                  << " Grad           " << chronoGrad.diffCumul() << " s.\n"
                  << " Grad Assemble  " << chronoGradAssemble.diffCumul() << " s.\n"
                  << " Div Assemble   " << chronoDivAssemble.diffCumul() << " s.\n"
                  << std::endl;

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
updateSystem ( const Real         alpha,
               const vector_Type& betaVector,
               const vector_Type& sourceVector )
{
    if ( M_matrixNoBC.get() )
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    }
    else
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap ) );
    }

    updateSystem ( alpha, betaVector, sourceVector, M_matrixNoBC, *M_un );
    if ( alpha != 0. )
    {
        M_matrixNoBC->globalAssemble();
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::
updateSystem ( const Real         alpha,
               const vector_Type& betaVector,
               const vector_Type& sourceVector,
               matrixPtr_Type     matrixNoBC,
               const vector_Type&     un )
{
    LifeChrono chrono;

    // clearing pressure mass matrix in case we need it in removeMean;
    M_pressureMatrixMass.reset( );


    M_Displayer.leaderPrint ( "  F-  Updating mass term on right hand side... " );

    chrono.start();

    UInt velocityTotalDof   = M_velocityFESpace.dof().numTotalDof();
    //    UInt pressureTotalDof = M_pressureFESpace.dof().numTotalDof();

    // Right hand side for the velocity at time

    updateRightHandSide ( sourceVector );

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );


    //    M_updated = false;

    if ( M_recomputeMatrix )
    {
        buildSystem();
    }

    M_Displayer.leaderPrint ( "  F-  Copying the matrices ...                 " );

    chrono.start();

    if ( M_isDiagonalBlockPreconditioner == true )
    {
        matrixPtr_Type tempMatrix ( M_blockPreconditioner );
        M_blockPreconditioner.reset ( new matrix_Type ( M_localMap,
                                                        M_blockPreconditioner->meanNumEntries() ) );
        *M_blockPreconditioner += *tempMatrix;
    }


    chrono.stop();
    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );


    UInt numVelocityComponent = M_velocityFESpace.fieldDim();

    //! managing the convective term

    Real normInf;
    betaVector.normInf ( &normInf );

    if ( normInf != 0. )
    {
        M_Displayer.leaderPrint ( "  F-  Sharing convective term ...              " );
        chrono.start();

        // vector with repeated nodes over the processors

        vector_Type betaVectorRepeated ( betaVector, Repeated );
        vector_Type unRepeated ( un, Repeated, Add );

        chrono.stop();

        M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
        M_Displayer.leaderPrint ( "  F-  Updating the convective terms ...        " );
        chrono.start();

        for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); ++iElement )
        {
            // just to provide the id number in the assem_mat_mixed
            M_pressureFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );
            //as updateFirstDer
            M_velocityFESpace.fe().updateFirstDeriv ( M_velocityFESpace.mesh()->element ( iElement ) );

            M_elementMatrixStiff.zero();

            UInt elementID = M_velocityFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elementRightHandSide contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = M_velocityFESpace.fe().patternFirst ( iNode );
                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = M_velocityFESpace.dof().localToGlobalMap ( elementID, iLocal )
                                   + iComponent * dimVelocity();
                    M_elementRightHandSide.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = betaVectorRepeated[iGlobal];

                    M_uLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = unRepeated (iGlobal);
                    M_wLoc.vec() [ iLocal + iComponent * M_velocityFESpace.fe().nbFEDof() ]
                        = unRepeated (iGlobal) - betaVectorRepeated (iGlobal);
                }
            }

            if (M_oseenData->conservativeFormulation() )
            {
                // ALE term: - rho div(w) u v
                mass_divw ( - M_oseenData->density(),
                            M_wLoc,
                            M_elementMatrixStiff,
                            M_velocityFESpace.fe(), 0, 0, numVelocityComponent );
            }

            // ALE stab implicit: 0.5 rho div u w v
            mass_divw ( 0.5 * M_oseenData->density(),
                        M_uLoc,
                        M_elementMatrixStiff,
                        M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // Stabilising term: -rho div(u^n) u v
            if ( M_divBetaUv )
                mass_divw ( -0.5 * M_oseenData->density(),
                            M_uLoc,
                            M_elementMatrixStiff,
                            M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // compute local convective terms
            advection ( M_oseenData->density(),
                        M_elementRightHandSide,
                        M_elementMatrixStiff,
                        M_velocityFESpace.fe(), 0, 0, numVelocityComponent );

            // loop on components
            for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
            {
                // compute local convective term and assembling
                // grad( 0, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );
                // grad( 1, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );
                // grad( 2, M_elementRightHandSide, M_elementMatrixStiff, M_velocityFESpace.fe(),
                //       M_velocityFESpace.fe(), iComponent, iComponent );

                assembleMatrix ( *matrixNoBC,
                                 M_elementMatrixStiff,
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.fe(),
                                 M_velocityFESpace.dof(),
                                 M_velocityFESpace.dof(),
                                 iComponent, iComponent,
                                 iComponent * velocityTotalDof, iComponent * velocityTotalDof );
            }
        }

        chrono.stop();
        M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

        if ( M_stabilization &&
                ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) ) )
        {
            M_Displayer.leaderPrint ( "  F-  Updating the stabilization terms ...     " );
            chrono.start();
            M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
            M_ipStabilization.apply ( *M_matrixStabilization, betaVectorRepeated, false );
            M_matrixStabilization->globalAssemble();
            M_resetStabilization = false;
            chrono.stop();
            M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
        }

    }
    else
    {
        if ( M_stabilization )
        {
            M_Displayer.leaderPrint ( "  F-  Updating the stabilization terms ...     " );
            chrono.start();

            if ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )
            {
                M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
                M_ipStabilization.apply ( *M_matrixStabilization, betaVector, false );
                M_matrixStabilization->globalAssemble();
                M_resetStabilization = false;
                chrono.stop();
                M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
            }
            else
            {
                M_Displayer.leaderPrint ( "reusing\n" );
            }
        }
    }

    if ( alpha != 0. )
    {
        *matrixNoBC += (*M_velocityMatrixMass) * alpha;
        if ( M_isDiagonalBlockPreconditioner == true )
        {
            matrixNoBC->globalAssemble();
            *M_blockPreconditioner += *matrixNoBC;
            matrix_Type tempMatrix ( *matrixNoBC );
            matrixNoBC.reset ( new matrix_Type ( M_localMap, tempMatrix.meanNumEntries() ) );
            *matrixNoBC += tempMatrix;
            M_blockPreconditioner->globalAssemble();
        }
    }
    *matrixNoBC += *M_matrixStokes;
}


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::updateStabilization ( matrix_Type& matrixFull )
{

    if ( M_stabilization )
    {
        matrixFull += *M_matrixStabilization;
    }

}

template <typename Mesh, typename SolverType>
void OseenSolver<Mesh, SolverType>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs ( M_localMap );

    VectorElemental M_elvec (M_velocityFESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= M_velocityFESpace->mesh()->numVolumes(); ++i )
    {

        M_velocityFESpace->fe().updateFirstDerivQuadPt ( M_velocityFESpace->mesh()->volumeList ( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec ( source, M_elvec, M_velocityFESpace->fe(),  M_oseenData->dataTime()->time(), ic ); // compute local vector
            assembleVector ( *rhs, M_elvec, M_velocityFESpace->fe(), M_velocityFESpace->dof(), ic, ic * M_velocityFESpace->getDim() ); // assemble local vector into global one
        }
    }
    M_rightHandSideNoBC += rhs;
}




template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::iterate ( bcHandler_Type& bcHandler )
{

    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint ( "  F-  Updating the boundary conditions ...     " );

    chrono.start();

    M_matrixNoBC->globalAssemble();

    matrixPtr_Type matrixFull ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );

    updateStabilization ( *matrixFull );
    getFluidMatrix ( *matrixFull );

    vector_Type rightHandSideFull ( *M_rightHandSideNoBC );

    //     matrixFull.reset( new matrix_Type( *M_matrixNoBC ) );
    //     M_rightHandSideFull = M_rightHandSideNoBC;

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint ("  F-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions ( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();
    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    // solving the system
    M_linearSolver->setMatrix ( *matrixFull );

    boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (matrixFull);

    Int numIter = M_linearSolver->solveSystem ( rightHandSideFull, *M_solution, staticCast );

    // if the preconditioner has been rese the stab terms are to be updated
    if ( numIter < 0 || numIter > M_iterReuseStabilization )
    {
        resetStabilization();
    }

    *M_residual  = *M_rightHandSideNoBC;
    *M_residual -= (*M_matrixNoBC) * (*M_solution);

    //M_residual.spy("residual");
} // iterate()


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::reduceSolution ( Vector& velocityVector, Vector& pressureVector )
{
    vector_Type solution ( *M_solution, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            velocityVector[ iDof ] = solution[ iDof ];
        }

        for ( UInt iDof = 0; iDof < dimPressure(); ++iDof )
        {
            pressureVector[ iDof ] = solution[ iDof + M_velocityFESpace.fieldDim() * dimVelocity() ];
        }
    }

}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::reduceResidual ( Vector& residualVector )
{
    vector_Type residual ( *M_residual, 0 );

    if ( false /*S_verbose*/ )
    {
        for ( UInt iDof = 0; iDof < M_velocityFESpace.fieldDim() * dimVelocity(); ++iDof )
        {
            residualVector[ iDof ] = residual[ iDof ];
        }

    }
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setBlockPreconditioner ( matrixPtr_Type blockPreconditioner )
{
    // blockPreconditioner.reset(new matrix_Type(M_monolithicMap, M_solid->getMatrixPtr()->getMeanNumEntries()));
    *blockPreconditioner += *M_blockPreconditioner;
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::getFluidMatrix ( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetArea()
{
    M_postProcessing->set_area();
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetNormal()
{
    M_postProcessing->set_normal();
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::postProcessingSetPhi()
{
    M_postProcessing->set_phi();
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::flux ( const markerID_Type& flag )
{
    return flux ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::flux ( const markerID_Type& flag,
                                          const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated, Add );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->flux ( velocity, flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::kineticNormalStress ( const markerID_Type& flag )
{
    return kineticNormalStress ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::kineticNormalStress ( const markerID_Type& flag,
                                                         const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated, Add );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->kineticNormalStress ( velocity, M_oseenData->density(), flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::area ( const markerID_Type& flag )
{
    return M_postProcessing->measure ( flag );
}

template<typename MeshType, typename SolverType>
Vector
OseenSolver<MeshType, SolverType>::normal ( const markerID_Type& flag )
{
    return M_postProcessing->normal ( flag );
}

template<typename MeshType, typename SolverType>
Vector
OseenSolver<MeshType, SolverType>::geometricCenter ( const markerID_Type& flag )
{
    return M_postProcessing->geometricCenter ( flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::pressure ( const markerID_Type& flag )
{
    return pressure ( flag, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::pressure (const markerID_Type& flag,
                                             const vector_Type& solution)
{
    vector_Type velocityAndPressure ( solution, Repeated, Add );
    vector_Type pressure ( this->M_pressureFESpace.map(), Repeated );
    pressure.subset ( velocityAndPressure,
                      this->M_velocityFESpace.dim() *this->M_velocityFESpace.fieldDim() );

    // third argument is 1, to use the pressure finite element space (see PostProcessingBoundary docs)
    return M_postProcessing->average ( pressure, flag, 1 ) [0];
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    if ( bcHandler.findBCWithFlag ( flag ).type() == Flux )
    {
        return -lagrangeMultiplier ( flag, bcHandler, solution );
    }
    else
    {
#ifdef HAVE_LIFEV_DEBUG
        M_Displayer.leaderPrint ( " !!! WARNING - OseenSolver::meanNormalStress( flag, bcHandler, solution) is returning an approximation \n" );
#endif
        return -pressure ( flag, solution ); // TODO: This is an approximation of the mean normal stress as the pressure.
        // A proper method should be coded in the PostprocessingBoundary class
        // to compute the exact mean normal stress.
    }
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanTotalNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::meanTotalNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    return meanNormalStress ( flag, bcHandler, solution ) - kineticNormalStress ( flag, solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::lagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::lagrangeMultiplier ( const markerID_Type&  flag,
                                                        bcHandler_Type& bcHandler,
                                                        const vector_Type& solution )
{
    // Create a list of Flux bcName_Type ??
    std::vector< bcName_Type > fluxBCVector = bcHandler.findAllBCWithType ( Flux );
    bcName_Type fluxbcName_Type = bcHandler.findBCWithFlag ( flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_Type velocityPressureLambda ( solution, Repeated, Add );

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex = 0; lmIndex < static_cast <UInt> ( fluxBCVector.size() ); ++lmIndex )
        if ( fluxbcName_Type.compare ( fluxBCVector[ lmIndex ] ) == 0 )
            return velocityPressureLambda[M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof()
                                          + M_pressureFESpace.dof().numTotalDof() + lmIndex];

    // If lmIndex has not been found a warning message is printed
    M_Displayer.leaderPrint (  "!!! Warning - Lagrange multiplier for Flux BC not found!\n" );
    return 0;
}

template<typename MeshType, typename SolverType>
Real
OseenSolver<MeshType, SolverType>::removeMean ( vector_Type& x )
{

    LifeChrono chrono;
    chrono.start();

    const UInt numVelocityComponent ( M_velocityFESpace.fieldDim() );
    const UInt velocityTotalDof ( M_velocityFESpace.dof().numTotalDof() );


    if ( M_pressureMatrixMass.get() == 0 )
    {
        M_pressureMatrixMass.reset ( new matrix_Type ( M_localMap ) );
    }

    for ( UInt iElement = 0; iElement < M_velocityFESpace.mesh()->numElements(); iElement++ )
    {
        chrono.start();
        // just to provide the id number in the assem_mat_mixed
        M_pressureFESpace.fe().update ( M_pressureFESpace.mesh()->element ( iElement ) );

        M_elementMatrixPreconditioner.zero();
        // mass
        chrono.start();
        mass ( 1, M_elementMatrixPreconditioner, M_pressureFESpace.fe(), 0, 0, M_velocityFESpace.fieldDim() );
        chrono.stop();

        chrono.start();
        assembleMatrix ( *M_pressureMatrixMass,
                         M_elementMatrixPreconditioner,
                         M_pressureFESpace.fe(),
                         M_pressureFESpace.fe(),
                         M_pressureFESpace.dof(),
                         M_pressureFESpace.dof(),
                         numVelocityComponent,
                         numVelocityComponent,
                         numVelocityComponent * velocityTotalDof,
                         numVelocityComponent * velocityTotalDof );
        chrono.stop();
    }

    M_pressureMatrixMass->globalAssemble();

    vector_Type ones ( *M_solution );
    ones = 1.0;

    Real mean;
    mean = ones * ( M_pressureMatrixMass * x );
    x += ( -mean );

    ASSERT ( std::fabs ( ones * ( M_pressureMatrixMass * x ) ) < 1e-9 , "after removeMean the mean pressure should be zero!");

    return mean;

} // removeMean()

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::applyBoundaryConditions ( matrix_Type&       matrix,
                                                             vector_Type&       rightHandSide,
                                                             bcHandler_Type& bcHandler )
{
    // M_rightHandSideFull = M_rightHandSideNoBC;

    // BC manage for the velocity
    if ( !bcHandler.bcUpdateDone() || M_recomputeMatrix )
    {
        bcHandler.bcUpdate ( *M_velocityFESpace.mesh(),
                             M_velocityFESpace.feBd(),
                             M_velocityFESpace.dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    //vector_Type rightHandSideFull( rightHandSide, Repeated, Zero );
    // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rightHandSideFull ( rightHandSide, Unique );

    bcManage ( matrix, rightHandSideFull,
               *M_velocityFESpace.mesh(),
               M_velocityFESpace.dof(),
               bcHandler,
               M_velocityFESpace.feBd(),
               1.,
               M_oseenData->dataTime()->time() );

    rightHandSide = rightHandSideFull;

    if ( bcHandler.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( M_velocityFESpace.fieldDim() *dimVelocity(),
                             M_diagonalize,
                             rightHandSide,
                             0. );
    }

} // applyBoundaryCondition

// ===================================================
// Set Methods
// ===================================================


template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setupPostProc( )
{
    M_postProcessing.reset ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
                                                                     &M_velocityFESpace.feBd(),
                                                                     &M_velocityFESpace.dof(),
                                                                     &M_pressureFESpace.feBd(),
                                                                     &M_pressureFESpace.dof(),
                                                                     M_localMap ) );
}

template<typename MeshType, typename SolverType>
void
OseenSolver<MeshType, SolverType>::setTolMaxIteration ( const Real& tolerance, const Int& maxIteration )
{
    M_linearSolver->setTolerance ( tolerance );
    M_linearSolver->setMaxNumIterations ( maxIteration );
}


} // namespace LifeV

#endif // OSEENSOLVER_H
