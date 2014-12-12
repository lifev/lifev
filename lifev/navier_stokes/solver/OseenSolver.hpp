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

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/eta/expression/Integrate.hpp>

#include <lifev/navier_stokes/testsuite/basic_test/uExactFunctor.hpp>
#include <lifev/navier_stokes/testsuite/basic_test/gradUExactFunctor.hpp>
#include <lifev/navier_stokes/testsuite/basic_test/pExactFunctor.hpp>

//#include <lifev/navier_stokes/solver/StabilizationSUPG.hpp>
//#include <lifev/navier_stokes/solver/StabilizationSUPGVMS.hpp>
#include <lifev/navier_stokes/solver/StabilizationVMSLES.hpp>

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

template< typename MeshType, typename SolverType = LifeV::SolverAztecOO, typename MapType = MapEpetra, UInt SpaceDim = 3, UInt FieldDim = 3>
class OseenSolver
{

public:

    //! @name Public Types
    //@{

	typedef MapType 									map_Type;
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

    // typedefs added to use block matrices

    typedef ETFESpace<mesh_Type, map_Type, SpaceDim, SpaceDim > ETFESpace_velocity;
    typedef ETFESpace<mesh_Type, map_Type, SpaceDim, 1 >        ETFESpace_pressure;

    typedef MatrixEpetraStructured<Real>           matrix_block_Type;
    typedef boost::shared_ptr<matrix_block_Type>   matrixPtr_block_Type;
    
    typedef VectorEpetraStructured                 vector_block_Type;
    typedef boost::shared_ptr<vector_block_Type>   vectorPtr_block_Type;

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
        @param u_star
        @param sourceVector
     */
    virtual void updateSystem ( const Real         alpha,
                                const vector_Type& u_star,
                                const vector_Type& rightHandSide );

    //! Update system
    /*!
        @param alpha
        @param u_star
        @param rightHandSide
        @param matrix
        @param un
     */
    virtual void updateSystem ( const Real         alpha,
                                const vector_Type& betaVector,
                                const vector_Type& rightHandSide,
                                matrixPtr_Type     matrix,
                                const vector_Type& un );

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

    //! Update and return the coefficient matrix
    /*!
        @param matrixFull The coefficient matrix
     */
    void getFluidMatrix ( matrix_Type& matrixFull );

    //! Set up post processing
    void setupPostProc( );

    //! Get the fluid matrix and the fluid right hand side. Called by the FSIMonolithic class when assembling the fluid block
    void getMatrixAndRhs( matrixPtr_Type& matrixFull, vectorPtr_Type& rhsFull)
    {
    	// Matrix
    	matrixFull.reset( new matrix_Type ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap ) );
    	updateStabilization ( *matrixFull );
    	getFluidMatrix ( *matrixFull );
    	matrixFull->globalAssemble();

    	// Right hand side
    	rhsFull.reset( new vector_Type (*M_rightHandSideNoBC) );
    	if(M_stabilization)
    	{
    		if(M_oseenData->stabilizationType() == "SUPG" || M_oseenData->stabilizationType() == "SUPGVMS" || M_oseenData->stabilizationType() == "VMSLES" )
    		{
    			*rhsFull += *M_rhsStabilization;
    		}
    	}
    }

    //! Set the initial velocity when restarting
    void initializeVelocitySolution(const vector_Type v0)
    {
    	M_solution->subset(v0, M_velocityFESpace.map(), 0, 0);
    }

    void initializePressureSolution(const vector_Type p0)
    {
    	M_solution->subset(p0, M_pressureFESpace.map(), 0, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof());
    }


    //! Set the initial pressure when restarting


    //! Compute the energy of the fluid
    /*!
        @return energy
     */
    Real energy ( );

    void preprocessing( BCHandler& bcloads, const GetPot& dataFile, BCHandler& bcDrag, BCHandler& bcHLift, vector_Type& solDrag, vector_Type& solLift);

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

    
    //! Compute the forces
    /*!
         @param flag      Flag of the boundary face associated with the boundary of the body.
         @param bcHandler BChandler containing the boundary conditions of the problem.
         @return          Forces
     */
    VectorSmall<2> computeForces(	const markerID_Type&  flag,
    								bcHandler_Type& bcHandlerDrag,
    								bcHandler_Type& bcHandlerLift);

    VectorSmall<2> computeForcesNewTest(	const vector_Type& drag_vector,
    										const vector_Type& lift_vector);


    //! Compute the drag
    /*!
     @param flag      Flag of the boundary face associated
     with the boundary of the body.
     @param bcHandler BChandler containing the boundary conditions of the problem.

     @return          Drag
     */
    VectorSmall<2> computeDrag(const markerID_Type&  flag,
                     bcHandler_Type& bcHandlerDrag,
                     bcHandler_Type& bcHandlerLift,
                     const Real& velocityInfty,
                     const Real& Area);
    
    //! Compute the drag
    /*!
        @param bcHandler bcHandlerDrag containing the boundary conditions for the drag.
		@param bcHandler bcHandlerLift containing the boundary conditions for the lift.
        @return          Forces
    */
    VectorSmall<2> computeForces ( bcHandler_Type& bcHandlerDrag, bcHandler_Type& bcHandlerLift );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary
    /*!
     *  @param uh1error value of the h1 norm of the velocity error
     */
    void h1normVelocity(Real& uh1error );

    //! compute the stabilization contributions that will be added to the system through the method updateStabilization()
    /*!
     *  @param betaVector
     *  @param alpha
     */
    void computeStabilization ( const vector_Type& betaVector, const Real& alpha );

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
    
    //! Set the velocity to be used by the stabilizations
    /*!
         @param velocityRHS
     */
    void setVelocityRhs( const vector_Type& velocityRHS )
    {
    	M_velocityRhs.reset(new vector_Type(velocityRHS.map(), Repeated));
    	*M_velocityRhs *= 0;
    	*M_velocityRhs += velocityRHS;
    }

    //! Set the velocity to be used by the stabilizations
    /*!
       	@param velocityRHS
     */
    void setPressureExtrapolated( const vector_Type& pressureExtrapolated )
    {
    	M_pressureExtrapolated.reset(new vector_Type(pressureExtrapolated.map(), Repeated));
    	*M_pressureExtrapolated *= 0;
    	*M_pressureExtrapolated += pressureExtrapolated;
    }

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
        return *M_matrixMass;
    }

    const matrix_Type& matrixMass() const
    {
        return *M_matrixMass;
    }

    const matrixPtr_Type matrixMassPtr() const
    {
        return M_matrixMass;
    }

    //@}

    //@{ unused methods

    //! Set up post processing structures
    void postProcessingSetArea();

    //! Set up post processing
    void postProcessingSetNormal();

    //! Set up post processing
    void postProcessingSetPhi();

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

    // ET FE Spaces
    boost::shared_ptr<ETFESpace_velocity > M_fespaceUETA;
    boost::shared_ptr<ETFESpace_pressure > M_fespacePETA;
    
    //! Displayer to print in parallel (only PID 0 will print)
    Displayer                      M_Displayer;

    //! Map of the whole system (velocity+pressure+fluxes)
    MapEpetra                      M_localMap;
    
    //! Map used for fluxes (when used)
    MapEpetra                      M_fluxMap;

    //! Mass matrix
    matrixPtr_block_Type          M_velocityMatrixMass;

    //! Mass matrix built in the pressure space
    matrixPtr_Type                M_pressureMatrixMass;
    
    //! Mass matrix built in the pressure space
    matrixPtr_Type                M_matrixMass;
    
    //! Convective term
    matrixPtr_block_Type          M_convectiveMatrix;
    
    //! Convective term
    matrixPtr_block_Type          M_matrixNoBC_block;

    //! Full matrix without boundary conditions
    matrixPtr_Type                M_matrixNoBC;

    //! stabilization matrix
    matrixPtr_Type                M_matrixStabilization;

    //! stabilization matrix
    matrixPtr_block_Type          M_matrixStabilizationET;

    //! Constant terms
    matrixPtr_block_Type          M_matrixStokes;
    
    //! source term for Navier-Stokes equations
    source_Type                    M_source;

    //! Right hand side for the velocity component
    vectorPtr_Type                 M_rightHandSideNoBC;

    //! Global solution
    vectorPtr_Type                 M_solution;
    
    //! Global solution
    vectorPtr_Type                 M_conservativeTerm;
    
    //! Global solution
    vectorPtr_Type                 M_betaVector;

    //! residual
    vectorPtr_Type                 M_residual;


    vectorPtr_Type                 M_velocityRhs;

    vectorPtr_Type                 M_velocityPreviousTimestep;

    vectorPtr_Type                 M_pressurePreviousTimestep;
    
    vectorPtr_Type 				   M_pressureExtrapolated;

    //! Linear solver
    linearSolverPtr_Type           M_linearSolver;

    vectorPtr_block_Type           M_rhsStabilization;

    //! True for steady simulations
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

    const function_Type*           M_betaFunction;

    bool                           M_divBetaUv;

    //! Use stiff-strain formulation
    bool                           M_stiffStrain;

    //
    Real                           M_diagonalize;

    UInt                           M_count;

    bool                           M_recomputeMatrix;

    // SUPG stabilization
    // boost::shared_ptr<StabilizationSUPG<mesh_Type, map_Type, SpaceDim > > M_supgStabilization;
    
    // SUPGVMS stabilization
    // boost::shared_ptr<StabilizationSUPGVMS<mesh_Type, map_Type, SpaceDim > > M_supgVmsStabilization;

    // VMSLES stabilization
    boost::shared_ptr<StabilizationVMSLES<mesh_Type, map_Type, SpaceDim > > M_VMSLESStabilization;

    //! Elementary matrices and vectors
    MatrixElemental                        M_elementMatrixStiff;            // velocity Stokes
    MatrixElemental                        M_elementMatrixMass;             // velocity mass
    MatrixElemental                        M_elementMatrixPreconditioner;   // (p,q) bloc for preconditioners
    MatrixElemental                        M_elementMatrixDivergence;
    MatrixElemental                        M_elementMatrixGradient;
    VectorElemental                        M_elementRightHandSide;          // Elementary right hand side
    VectorElemental                        M_wLoc;
    VectorElemental                        M_uLoc;
    boost::shared_ptr<vector_Type> M_un;

}; // class OseenSolver


class buildVector
{
public:
    typedef VectorSmall<3> return_Type;

    inline return_Type operator() (Real a, Real b, Real c)
    {
        VectorSmall<3> o;
        o[0] = a;
        o[1] = b;
        o[2] = c;
        return o;
    }

    buildVector() {}
    buildVector (const buildVector&) {}
    ~buildVector() {}
};

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              const Int                       lagrangeMultiplier ) :
	  M_oseenData       ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
	  M_fluxMap                ( lagrangeMultiplier, communicator),
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
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_un                     ( new vector_Type (M_localMap) ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_pressureExtrapolated   ( new vector_Type (M_pressureFESpace.map() ) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      //M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      //M_supgVmsStabilization   (new StabilizationSUPGVMS<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    *M_solution *= 0;
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              boost::shared_ptr<Epetra_Comm>& communicator,
              MapEpetra                       monolithicMap,
              UInt                            offset ) :
	  M_oseenData              ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
      M_fluxMap                ( offset, communicator),
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
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), M_velocityFESpace.fieldDim() ),
	  M_un                     ( /*new vector_Type(M_localMap)*/ ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_pressureExtrapolated   ( new vector_Type (M_pressureFESpace.map() ) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      //M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      //M_supgVmsStabilization   (new StabilizationSUPGVMS<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    ASSERT(0!=0,"ENTERING IN THE MONOLITHIC CONSTRUCTOR");
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
OseenSolver ( boost::shared_ptr<data_Type>    dataType,
              FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
              FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
              const std::vector<Int>&         lagrangeMultipliers,
              boost::shared_ptr<Epetra_Comm>& communicator ) :
	  M_oseenData              ( dataType ),
	  M_velocityFESpace        ( velocityFESpace ),
	  M_pressureFESpace        ( pressureFESpace ),
	  M_Displayer              ( communicator ),
	  M_fluxMap                ( lagrangeMultipliers, communicator),
	  M_localMap               ( M_velocityFESpace.map() + M_pressureFESpace.map() + lagrangeMultiplier),
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
	  M_wLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_uLoc                   ( M_velocityFESpace.fe().nbFEDof(), velocityFESpace.fieldDim() ),
	  M_un                     ( new vector_Type (M_localMap) ),
      M_velocityPreviousTimestep (new vector_Type (M_velocityFESpace.map()) ),
      M_pressurePreviousTimestep (new vector_Type (M_pressureFESpace.map()) ),
      M_pressureExtrapolated   ( new vector_Type (M_pressureFESpace.map() ) ),
      M_fespaceUETA            ( new ETFESpace_velocity(M_velocityFESpace.mesh(), &(M_velocityFESpace.refFE()), communicator)),
      M_fespacePETA            ( new ETFESpace_pressure(M_pressureFESpace.mesh(), &(M_pressureFESpace.refFE()), communicator)),
      //M_supgStabilization      (new StabilizationSUPG<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      //M_supgVmsStabilization   (new StabilizationSUPGVMS<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace)),
      M_VMSLESStabilization    (new StabilizationVMSLES<mesh_Type, MapEpetra, SpaceDim>(velocityFESpace, pressureFESpace))
{
    *M_solution *= 0;
    // if(M_stabilization = ( &M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ))
    {
        M_ipStabilization.setFeSpaceVelocity (M_velocityFESpace);
        M_ipStabilization.setViscosity (M_oseenData->viscosity() );
    }
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
~OseenSolver()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setUp ( const GetPot& dataFile )
{
    if (M_linearSolver.get() )
    {
        M_linearSolver->setupPreconditioner ( dataFile, "fluid/prec" );
        M_linearSolver->setDataFromGetPot ( dataFile, "fluid/solver" );
    }

    M_stabilization = dataFile ( "fluid/stabilization/use",  /*(&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE() ) ,*/ false);

    // If using P1-P1 the use of the stabilization is necessary
    //if(&M_velocityFESpace.refFE() == &M_pressureFESpace.refFE())
    //	M_stabilization = true;

    M_steady        = dataFile ( "fluid/miscellaneous/steady", 0 );

    Real bdfOrder = dataFile ( "fluid/time_discretization/BDF_order", 1.0 );

    if (M_stabilization)
    {
    	if(M_oseenData->stabilizationType() == "IP")
    	{
    		M_gammaBeta     = dataFile ( "fluid/ipstab/gammaBeta",  0. );
    		M_gammaDiv      = dataFile ( "fluid/ipstab/gammaDiv",   0. );
    		M_gammaPress    = dataFile ( "fluid/ipstab/gammaPress", 0. );
    		M_reuseStabilization     = dataFile ( "fluid/ipstab/reuse", false );

    		M_ipStabilization.setGammaBeta ( M_gammaBeta );
    		M_ipStabilization.setGammaDiv  ( M_gammaDiv );
    		M_ipStabilization.setGammaPress ( M_gammaPress );

    		if (M_linearSolver.get() )
    			M_iterReuseStabilization = dataFile ( "fluid/ipstab/max_iter_reuse", static_cast<Int> ( M_linearSolver->maxNumIterations() * 8. / 10. ) );
    	}
    	else if (M_oseenData->stabilizationType() == "SUPG")
    	{
    		/*
    		int vel_order = dataFile ( "fluid/space_discretization/vel_order", 1 );
    		M_supgStabilization->setConstant ( vel_order );
    		M_supgStabilization->setETvelocitySpace(M_fespaceUETA);
    		M_supgStabilization->setETpressureSpace(M_fespacePETA);
    		M_supgStabilization->setCommunicator(M_velocityFESpace.map().commPtr());
    		M_supgStabilization->setDensity(M_oseenData->density());
    		M_supgStabilization->setBDForder(bdfOrder);
    		M_supgStabilization->setViscosity(M_oseenData->viscosity());
    		M_supgStabilization->setTimeStep(M_oseenData->dataTime()->timeStep());
    		*/

    	}
    	else if (M_oseenData->stabilizationType() == "SUPGVMS")
    	{
    		/*
    		int vel_order = dataFile ( "fluid/space_discretization/vel_order", 1 );
    		M_supgVmsStabilization->setConstant ( vel_order );
    		M_supgVmsStabilization->setETvelocitySpace(M_fespaceUETA);
    		M_supgVmsStabilization->setETpressureSpace(M_fespacePETA);
    		M_supgVmsStabilization->setCommunicator(M_velocityFESpace.map().commPtr());
    		M_supgVmsStabilization->setDensity(M_oseenData->density());
    		M_supgVmsStabilization->setBDForder(bdfOrder);
    		M_supgVmsStabilization->setViscosity(M_oseenData->viscosity());
    		M_supgVmsStabilization->setTimeStep(M_oseenData->dataTime()->timeStep());
    		*/

    	}
        else if (M_oseenData->stabilizationType() == "VMSLES")
    	{
            *M_velocityPreviousTimestep *= 0;
            *M_pressurePreviousTimestep *= 0;
    		int vel_order = dataFile ( "fluid/space_discretization/vel_order", 1 );
    		M_VMSLESStabilization->setConstant ( vel_order );
    		M_VMSLESStabilization->setETvelocitySpace(M_fespaceUETA);
    		M_VMSLESStabilization->setETpressureSpace(M_fespacePETA);
    		M_VMSLESStabilization->setCommunicator(M_velocityFESpace.map().commPtr());
    		M_VMSLESStabilization->setDensity(M_oseenData->density());
    		M_VMSLESStabilization->setBDForder(bdfOrder);
    		M_VMSLESStabilization->setViscosity(M_oseenData->viscosity());
    		M_VMSLESStabilization->setTimeStep(M_oseenData->dataTime()->timeStep());
    	}
	}
    // Energetic stabilization term
    M_divBetaUv   = dataFile ( "fluid/space_discretization/div_beta_u_v", false);
    // Enable grad( u )^T in stress tensor
    M_stiffStrain = dataFile ( "fluid/space_discretization/stiff_strain", false);
    M_diagonalize = dataFile ( "fluid/space_discretization/diagonalize", 1. );

    //    M_linearSolver.setAztecooPreconditioner( dataFile, "fluid/solver" );

}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
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


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
initialize ( const vector_Type& velocityInitialGuess, const vector_Type& pressureInitialGuess )
{

    *M_solution = velocityInitialGuess;
    *M_un = velocityInitialGuess;
    M_solution->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );
    M_un->add ( pressureInitialGuess, M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof() );

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
initialize ( const vector_Type& velocityAndPressure )
{

    *M_un = velocityAndPressure;
    if ( M_solution.get() )
    {
        *M_solution = velocityAndPressure;
    }

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::preprocessing( BCHandler& bcloads,
																				const GetPot& dataFile, BCHandler& bcDrag,
																				BCHandler& bcHLift,
																				vector_Type& solDrag,
																				vector_Type& solLift )
{
	matrixPtr_Type laplacian;
	laplacian.reset( new matrix_Type ( M_velocityFESpace.map() ) );

	{
		using namespace ExpressionAssembly;

		integrate(
				elements(M_fespaceUETA->mesh()),
				M_velocityFESpace.qr(),
				M_fespaceUETA,
				M_fespaceUETA,
				dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
		) >> laplacian;
	}

	vector_Type rhs_drag( M_velocityFESpace.map(), Unique );
	rhs_drag *= 0;

	vector_Type rhs_lift( M_velocityFESpace.map(), Unique );
	rhs_lift *= 0;

	vectorPtr_Type solution_drag;
	solution_drag.reset( new vector_Type( M_velocityFESpace.map(), Unique ) );

	vectorPtr_Type solution_lift;
	solution_lift.reset( new vector_Type( M_velocityFESpace.map(), Unique ) );

	bcloads.bcUpdate ( *M_velocityFESpace.mesh(), M_velocityFESpace.feBd(), M_velocityFESpace.dof() );

	bcManage ( *laplacian, rhs_drag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(), bcloads, M_velocityFESpace.feBd(), 1.0, 0.0 );

	bcManageRhs ( rhs_drag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcDrag, M_velocityFESpace.feBd(), 1., 0.);
	bcManageRhs ( rhs_lift, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHLift, M_velocityFESpace.feBd(), 1., 0.);

	laplacian->globalAssemble();
	rhs_drag.globalAssemble();
	rhs_lift.globalAssemble();

	// drag
	linearSolverPtr_Type linearSolver_drag;
	linearSolver_drag.reset( new linearSolver_Type (M_velocityFESpace.map().commPtr()) );

	linearSolver_drag->setupPreconditioner ( dataFile, "preprocessing/prec" );
	linearSolver_drag->setDataFromGetPot ( dataFile, "preprocessing/solver" );

	linearSolver_drag->setMatrix ( *laplacian );

	boost::shared_ptr<MatrixEpetra<Real> > staticCast_drag = boost::static_pointer_cast<MatrixEpetra<Real> > (laplacian);
	Int numIter_drag = linearSolver_drag->solveSystem ( rhs_drag, *solution_drag, staticCast_drag );

	solDrag = *solution_drag;

	// lift
	linearSolverPtr_Type linearSolver_lift;
	linearSolver_lift.reset( new linearSolver_Type (M_velocityFESpace.map().commPtr()) );

	linearSolver_lift->setupPreconditioner ( dataFile, "preprocessing/prec" );
	linearSolver_lift->setDataFromGetPot ( dataFile, "preprocessing/solver" );

	linearSolver_lift->setMatrix ( *laplacian );

	boost::shared_ptr<MatrixEpetra<Real> > staticCast_lift = boost::static_pointer_cast<MatrixEpetra<Real> > (laplacian);
	Int numIter_lift = linearSolver_lift->solveSystem ( rhs_lift, *solution_lift, staticCast_lift );

	solLift = *solution_lift;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::buildSystem()
{
	// New ETA PART
	/*
	 *  TESTING THE ASSEMBLY USING ETA
	 *
	 */

	M_Displayer.leaderPrint ( "  F-  Computing constant matrices ...          " );
	LifeChrono chrono;
	chrono.start();

	M_velocityMatrixMass.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
	M_matrixStokes.reset(new matrix_block_Type( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ));

    *M_velocityMatrixMass *= 0;
	*M_matrixStokes *= 0;

	{
		using namespace ExpressionAssembly;

        if ( !M_steady )
		{
			integrate(
					elements(M_fespaceUETA->mesh()),
					M_velocityFESpace.qr(),
					M_fespaceUETA,
					M_fespaceUETA,
					M_oseenData->density() * dot(phi_i, phi_j)
			) >> M_velocityMatrixMass->block(0,0);

            M_velocityMatrixMass->globalAssemble();

            M_matrixMass.reset(new matrix_Type(M_localMap));
            *M_matrixMass += *M_velocityMatrixMass;
            M_matrixMass->globalAssemble();
		}

        if ( M_stiffStrain )
        {
            integrate(
                    elements(M_fespaceUETA->mesh()),
                    M_velocityFESpace.qr(),
                    M_fespaceUETA,
                    M_fespaceUETA,
                    value( 0.5 * M_oseenData->viscosity() ) * dot( grad(phi_i) + transpose(grad(phi_i)) , grad(phi_j) + transpose(grad(phi_j)) )
            ) >> M_matrixStokes->block(0,0);
        }
        else
        {
            integrate(
                    elements(M_fespaceUETA->mesh()),
                    M_velocityFESpace.qr(),
                    M_fespaceUETA,
                    M_fespaceUETA,
                    M_oseenData->viscosity() * dot( grad(phi_i) , grad(phi_j) + transpose(grad(phi_j)) )
            ) >> M_matrixStokes->block(0,0);
        }


		integrate(
				elements(M_fespaceUETA->mesh()),
				M_velocityFESpace.qr(),
				M_fespaceUETA,
				M_fespacePETA,
				value(-1.0) * phi_j * div(phi_i)
		) >> M_matrixStokes->block(0,1);

		integrate(
				elements(M_fespaceUETA->mesh()),
				M_pressureFESpace.qr(),
				M_fespacePETA,
				M_fespaceUETA,
				phi_i * div(phi_j)
		) >> M_matrixStokes->block(1,0);
	}

    M_matrixStokes->globalAssemble();

	chrono.stop();
	M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
updateSystem ( const Real         alphaOverTimestep,
               const vector_Type& u_star,
               const vector_Type& rightHandSide )
{
    if ( M_matrixNoBC.get() )
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap, M_matrixNoBC->meanNumEntries() ) );
    }
    else
    {
        M_matrixNoBC.reset ( new matrix_Type ( M_localMap ) );
    }

    updateSystem ( alphaOverTimestep, u_star, rightHandSide, M_matrixNoBC, *M_un );
    if ( alphaOverTimestep != 0. )
    {
        M_matrixNoBC->globalAssemble();
    }

}

class SignFunction
{
public:
    typedef Real return_Type;

    inline return_Type operator() (const Real& a)
    {
    	if ( a < 0 )
    		return a;
    	else
    		return 0.0;
    }

    SignFunction() {}
    SignFunction (const SignFunction&) {}
    ~SignFunction() {}
};

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::
updateSystem ( const Real          alphaOverTimestep,
		       const vector_Type&  u_star,
               const vector_Type&  rightHandSide,
               matrixPtr_Type      matrixNoBC,
               const vector_Type&  un )
{
    M_Displayer.leaderPrint ( "  F-  Setting the right hand side ..." );
    LifeChrono chrono;
    chrono.start();

    updateRightHandSide ( rightHandSide );

    chrono.stop();
    M_Displayer.leaderPrintMax ( "          done in ", chrono.diff() );

    M_matrixNoBC_block.reset(new matrix_block_Type(M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap));

    if ( M_recomputeMatrix )
        buildSystem();

    //! managing the convective term : semi-implicit approximation of the convective term

    Real normInf;
    u_star.normInf ( &normInf );

    MatrixSmall<3, 3> Eye;
    Eye *= 0.0;
    Eye[0][0] = 1;
    Eye[1][1] = 1;
    Eye[2][2] = 1;

    // u_star: extrapolation of the velocity
    // NOTE:   for ALE formulation it has to be already: extrapolation of fluid velocity - extrapolation of fluid mesh velocity
    vector_Type u_starRepeated ( u_star, Repeated );

    M_convectiveMatrix.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
    *M_convectiveMatrix *= 0;

    boost::shared_ptr<SignFunction> signEvaluation(new SignFunction());

    // only if the extrapolation of the velocity differs from zero we update the convective term
    if ( normInf != 0. )
    {
        M_convectiveMatrix.reset  ( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
        *M_convectiveMatrix *= 0;

        {
            using namespace ExpressionAssembly;
            integrate(
                        elements(M_fespaceUETA->mesh()), // Mesh
                        M_velocityFESpace.qr(),               // Quadrature rule
                        M_fespaceUETA,
                        M_fespaceUETA,
                        value(M_oseenData->density())*dot( value(M_fespaceUETA, u_starRepeated)*grad(phi_j), phi_i) // semi-implicit treatment of the convective term
                        //value( -1.0*M_oseenData->density() )* dot( grad(phi_i), outerProduct( value(M_fespaceUETA, u_starRepeated), phi_j ) )
                     )
            >> M_convectiveMatrix->block(0,0);

            VectorSmall<3> normal;
            normal[0] = 1.0;
            normal[1] = 0.0;
            normal[2] = 0.0;

            // Reverse flow at the ouflow
            QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria6pt) );
            integrate (
            			boundary (M_fespaceUETA->mesh(), 3), // 3 is for the flag of the outflow
                        myBDQR,
                        M_fespaceUETA,
                        M_fespaceUETA,
                        value(-1.0*M_oseenData->density())*eval( signEvaluation, dot(value(M_fespaceUETA, u_starRepeated),normal))*dot(phi_i, phi_j)
                      )
            >> M_convectiveMatrix->block(0,0);
        }

        // Stabilising term: -rho div(u^n) u v
        if ( M_divBetaUv )
        {
            ASSERT(0!=0,"If really needed, please code it using ET");
            // OLD-STYLE assembly:
            // mass_divw ( -0.5 * M_oseenData->density(),
            //           M_uLoc,
            //           M_elementMatrixStiff,
            //           M_velocityFESpace.fe(), 0, 0, numVelocityComponent );
        }
    }

    if(M_stabilization)
        computeStabilization(u_starRepeated, alphaOverTimestep);

    if ( alphaOverTimestep != 0. )
    {
        *M_matrixNoBC_block += (*M_velocityMatrixMass) * alphaOverTimestep;
    }

    *M_convectiveMatrix += *M_matrixStokes;
    M_convectiveMatrix->globalAssemble();

    *M_matrixNoBC_block += *M_convectiveMatrix;
    M_matrixNoBC_block->globalAssemble();
    MapEpetra fullMap ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap );
    M_matrixNoBC.reset(new matrix_Type(fullMap));
    *M_matrixNoBC += *M_matrixNoBC_block;
    M_matrixNoBC->globalAssemble();

}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeStabilization ( const vector_Type& u_star, const Real& alpha )
{
	Real normInf;
	u_star.normInf ( &normInf );

	LifeChrono chrono;

	if ( normInf != 0. )
	{
		if( M_oseenData->stabilizationType() == "IP" /*&& ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )*/ )
		{
			vector_Type u_starRepeated ( u_star, Repeated );
			M_Displayer.leaderPrint ( "  F-  Updating the IP stabilization terms ... " );
			chrono.start();

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			M_ipStabilization.apply ( *M_matrixStabilization, u_starRepeated, false );
			M_matrixStabilization->globalAssemble();
			M_resetStabilization = false;
			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
 		}
		else if(M_oseenData->stabilizationType() == "SUPG")
		{
			/*
			M_Displayer.leaderPrint ( "  F-  Updating the SUPG stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgStabilization->applySUPG_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);
			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgStabilization->applySUPG_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			*/
		}
		else if(M_oseenData->stabilizationType() == "SUPGVMS")
		{
			/*
			M_Displayer.leaderPrint ( "  F-  Updating the SUPG-VMS stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgVmsStabilization->applySUPGVMS_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);
			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgVmsStabilization->applySUPGVMS_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			*/
		}
        else if(M_oseenData->stabilizationType() == "VMSLES")
		{

			M_Displayer.leaderPrint ( "  F-  Updating the VMSLES stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_VMSLESStabilization->applyVMSLES_Matrix_semi_implicit(M_matrixStabilizationET,
                                                                    u_star,
                                                                    alpha,
                                                                    *M_pressureExtrapolated,
                                                                    *M_velocityRhs);

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_VMSLESStabilization->applyVMSLES_RHS_semi_implicit(M_rhsStabilization,
                                                                 u_star,
                                                                 *M_velocityRhs,
                                                                 alpha,
                                                                 *M_pressureExtrapolated);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

		}
	}
	else
	{
		if (M_oseenData->stabilizationType() == "IP")
		{
			M_Displayer.leaderPrint ( "  F-  Updating the IP stabilization terms ... " );
			chrono.start();

			if ( M_resetStabilization || !M_reuseStabilization || ( M_matrixStabilization.get() == 0 ) )
			{
				vector_Type u_starUnique ( u_star, Unique );
				M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
				M_ipStabilization.apply ( *M_matrixStabilization, u_starUnique, false );
				M_resetStabilization = false;
				chrono.stop();
				M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			}
			else
			{
				M_Displayer.leaderPrint ( "reusing\n" );
			}
		}
		else if(M_oseenData->stabilizationType() == "SUPG")
		{
			/*
			M_Displayer.leaderPrint ( "  F-  Updating the SUPG stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgStabilization->applySUPG_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			// comment because if steady simulation this is not needed
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgStabilization->applySUPG_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			*/
		}
		else if(M_oseenData->stabilizationType() == "SUPGVMS")
		{
			/*
			M_Displayer.leaderPrint ( "  F-  Updating the SUPG-VMS stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_supgVmsStabilization->applySUPGVMS_Matrix_semi_implicit(M_matrixStabilizationET, u_star, alpha);

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			// comment because if steady simulation this is not needed
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_supgVmsStabilization->applySUPGVMS_RHS_semi_implicit(M_rhsStabilization, u_star, *M_velocityRhs);

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );
			*/
		}
        else if(M_oseenData->stabilizationType() == "VMSLES")
		{

			M_Displayer.leaderPrint ( "  F-  Updating the VMSLES stabilization terms ... " );
			chrono.start();

			// TIpically here alpha is already divided by the timestep, but I want to use the actual alfa, so I multiply
			Real alfa = alpha*M_oseenData->dataTime()->timeStep();
			M_matrixStabilizationET.reset( new matrix_block_Type ( M_fespaceUETA->map() | M_fespacePETA->map() | M_fluxMap ) );
			*M_matrixStabilizationET *= 0;
			M_VMSLESStabilization->applyVMSLES_Matrix_semi_implicit(M_matrixStabilizationET,
                    												u_star,
                    												alpha,
                    												*M_pressureExtrapolated,
                    												*M_velocityRhs);

			M_matrixStabilization.reset ( new matrix_Type ( M_localMap ) );
			*M_matrixStabilization += *M_matrixStabilizationET;
			M_matrixStabilization->globalAssemble();

			// comment because if steady simulation this is not needed
			M_rhsStabilization.reset(new vector_block_Type( M_velocityFESpace.map() | M_pressureFESpace.map() | M_fluxMap ));
			*M_rhsStabilization *= 0;
			M_VMSLESStabilization->applyVMSLES_RHS_semi_implicit( M_rhsStabilization,
                    											  u_star,
                    											  *M_velocityRhs,
                    											  alpha,
                    											  *M_pressureExtrapolated );

			chrono.stop();
			M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

		}
	}
}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::updateStabilization ( matrix_Type& matrixFull )
{
    if ( M_stabilization )
    {
    	M_matrixStabilization->globalAssemble();
        matrixFull += *M_matrixStabilization;
    }

}

template <typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap );

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
    M_rightHandSideNoBC *= 0;
    M_rightHandSideNoBC += rhs;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::iterate ( bcHandler_Type& bcHandler )
{

    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer.leaderPrint ( "  F-  Updating the boundary conditions ...     " );

    // HERE I SHOULD APPLY THE STABILIZATION ON THE RHS

    chrono.start();

    M_matrixNoBC->globalAssemble();

    matrixPtr_Type matrixFull ( new matrix_Type ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap ) );

    updateStabilization ( *matrixFull );

    getFluidMatrix ( *matrixFull );

    vector_Type rightHandSideFull ( *M_rightHandSideNoBC );

    if(M_stabilization)
    {
        if(M_oseenData->stabilizationType() == "SUPG" || M_oseenData->stabilizationType() == "SUPGVMS" || M_oseenData->stabilizationType() == "VMSLES" )
        {
            rightHandSideFull += *M_rhsStabilization;
        }
    }

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );

    // boundary conditions update
    M_Displayer.leaderPrint ("  F-  Applying boundary conditions ...         ");

    chrono.start();

    applyBoundaryConditions ( *matrixFull, rightHandSideFull, bcHandler );

    matrixFull->globalAssemble();

    chrono.stop();

    M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    // Metto il tempo attuale in una stinga
//    std::string time = static_cast<std::ostringstream*>( &(ostringstream() << M_oseenData->dataTime()->time()) )->str();

    // Spy della matrice usata per risolvere il sistema lineare
//    std::string name = "matrixLinearSystem_" + time;
//    matrixFull->spy(name);

    // solving the system
    M_linearSolver->setMatrix ( *matrixFull );

    boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (matrixFull);

    // Spy del rhs usato nel sistema lineare
//    std::string nameVec = "rhsLinearSystem_" + time;
//    rightHandSideFull.spy(nameVec);

    Int numIter = M_linearSolver->solveSystem ( rightHandSideFull, *M_solution, staticCast );

    // Spy della souzione del sistema lineare
//    std::string nameSol = "solutionLinearSystem_" + time;
//    M_solution->spy(nameSol);

    // if the preconditioner has been rese the stab terms are to be updated
    if ( numIter < 0 || numIter > M_iterReuseStabilization )
    {
        resetStabilization();
    }

    // Computation of the residual for the coefficients of lift and drag
    matrixPtr_Type matrixNoBC ( new matrix_Type ( M_velocityFESpace.map() + M_pressureFESpace.map() + M_fluxMap ) );
    *matrixNoBC *= 0;
    updateStabilization ( *matrixNoBC );
    getFluidMatrix ( *matrixNoBC );
    matrixNoBC->globalAssemble();

    vector_Type rightHandSideNoBC ( *M_rightHandSideNoBC );
    if(M_stabilization)
    {
    	if(M_oseenData->stabilizationType() == "SUPG" || M_oseenData->stabilizationType() == "SUPGVMS" || M_oseenData->stabilizationType() == "VMSLES" )
    	{
    		rightHandSideNoBC += *M_rhsStabilization;
    	}
    }
    rightHandSideNoBC.globalAssemble();

    *M_residual *= 0;

    // Spy della matrice usata per risolvere il sistema lineare
//    std::string matrixRes = "matrixResidual_" + time;
//    matrixNoBC->spy(matrixRes);

    // Spy del rhs usato per il calcolo del residuo
//    std::string nameRhsResidual = "rhsResidual_" + time;
//    rightHandSideNoBC.spy(nameRhsResidual);

    *M_residual  = rightHandSideNoBC;
    *M_residual -= (*matrixNoBC) * (*M_solution);

    // Spy della soluzioneusata nel sistema lineare
//    std::string nameResidual = "residual_" + time;
//    rightHandSideFull.spy(nameResidual);

    //M_residual.spy("residual");
} // iterate()

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::h1normVelocity(Real& uh1error )
{
	boost::shared_ptr<uExactFunctor> uExactFct( new uExactFunctor );
	boost::shared_ptr<vector_Type> uExactVec(new vector_Type(M_velocityFESpace.map(),Unique));
	M_velocityFESpace.interpolate(RossEthierSteinmanUnsteadyDec::uexact, *uExactVec, 0.0);

	vector_Type uComputed ( M_fespaceUETA->map() , Unique );
	uComputed.subset( *M_solution );

	Real errorH1SquaredLocal( 0.0 );
	Real errorH1Squared( 0.0 );

	boost::shared_ptr<gradUExactFunctor> gradUExactFct( new gradUExactFunctor );

	{
		using namespace ExpressionAssembly;
		integrate (
				elements (M_fespaceUETA->mesh() ), // Mesh
				M_velocityFESpace.qr(), // QR
				dot ( ( eval ( gradUExactFct, value(M_oseenData->dataTime()->time()), X) + (-1) * grad( M_fespaceUETA , uComputed ) ) ,
					  ( eval ( gradUExactFct, value(M_oseenData->dataTime()->time()), X) + (-1) * grad( M_fespaceUETA , uComputed ) ) )
		)
		>> errorH1SquaredLocal;

	}

	M_fespaceUETA->mesh()->comm()->Barrier();
	M_fespaceUETA->mesh()->comm()->SumAll (&errorH1SquaredLocal, &errorH1Squared, 1);
	uh1error= std::sqrt(errorH1Squared);
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::reduceSolution ( Vector& velocityVector, Vector& pressureVector )
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

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::reduceResidual ( Vector& residualVector )
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

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::getFluidMatrix ( matrix_Type& matrixFull )
{
    M_matrixNoBC->globalAssemble();
    matrixFull += *M_matrixNoBC;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetArea()
{
    M_postProcessing->set_area();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetNormal()
{
    M_postProcessing->set_normal();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::postProcessingSetPhi()
{
    M_postProcessing->set_phi();
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::flux ( const markerID_Type& flag )
{
    return flux ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::flux ( const markerID_Type& flag,
                                          const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->flux ( velocity, flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::kineticNormalStress ( const markerID_Type& flag )
{
    return kineticNormalStress ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::kineticNormalStress ( const markerID_Type& flag,
                                                         const vector_Type& solution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    return M_postProcessing->kineticNormalStress ( velocity, M_oseenData->density(), flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::area ( const markerID_Type& flag )
{
    return M_postProcessing->measure ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Vector
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::normal ( const markerID_Type& flag )
{
    return M_postProcessing->normal ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Vector
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::geometricCenter ( const markerID_Type& flag )
{
    return M_postProcessing->geometricCenter ( flag );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::pressure ( const markerID_Type& flag )
{
    return pressure ( flag, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::pressure (const markerID_Type& flag,
                                             const vector_Type& solution)
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type pressure ( this->M_pressureFESpace.map(), Repeated );
    pressure.subset ( velocityAndPressure,
                      this->M_velocityFESpace.dim() *this->M_velocityFESpace.fieldDim() );

    // third argument is 1, to use the pressure finite element space (see PostProcessingBoundary docs)
    return M_postProcessing->average ( pressure, flag, 1 ) [0];
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
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

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return meanTotalNormalStress ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::meanTotalNormalStress (const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution )
{
    return meanNormalStress ( flag, bcHandler, solution ) - kineticNormalStress ( flag, solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::lagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier ( flag, bcHandler, *M_solution );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::lagrangeMultiplier ( const markerID_Type&  flag,
                                                        bcHandler_Type& bcHandler,
                                                        const vector_Type& solution )
{
    // Create a list of Flux bcName_Type ??
    std::vector< bcName_Type > fluxBCVector = bcHandler.findAllBCWithType ( Flux );
    bcName_Type fluxbcName_Type = bcHandler.findBCWithFlag ( flag ).name();

    // Create a Repeated vector for looking to the lambda
    vector_Type velocityPressureLambda ( solution, Repeated );

    // Find the index associated to the correct Lagrange multiplier
    for ( UInt lmIndex = 0; lmIndex < static_cast <UInt> ( fluxBCVector.size() ); ++lmIndex )
        if ( fluxbcName_Type.compare ( fluxBCVector[ lmIndex ] ) == 0 )
            return velocityPressureLambda[M_velocityFESpace.fieldDim() * M_velocityFESpace.dof().numTotalDof()
                                          + M_pressureFESpace.dof().numTotalDof() + lmIndex];

    // If lmIndex has not been found a warning message is printed
    M_Displayer.leaderPrint (  "!!! Warning - Lagrange multiplier for Flux BC not found!\n" );
    return 0;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::energy ( )
{
	Real energy = 0.0;

	vector_Type tmp(M_localMap);
	tmp *= 0;
	tmp = (*M_matrixMass) * (*M_solution);
	energy = M_solution->dot(tmp);
	energy *= 0.5;

	return energy;
}


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
VectorSmall<2>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeForces ( bcHandler_Type& bcHandlerDrag,
                                                                               	 bcHandler_Type& bcHandlerLift )
{
    vector_Type onesOnBodyDrag(M_localMap, Unique);
    onesOnBodyDrag *= 0;

    vector_Type onesOnBodyLift(M_localMap, Unique);
    onesOnBodyLift *= 0;

    bcManageRhs ( onesOnBodyDrag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerDrag, M_velocityFESpace.feBd(), 1., 0.);
    bcManageRhs ( onesOnBodyLift, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerLift, M_velocityFESpace.feBd(), 1., 0.);

    Real drag (0.0);
    Real lift (0.0);

    drag = M_residual->dot(onesOnBodyDrag);
    lift = M_residual->dot(onesOnBodyLift);

    M_Displayer.leaderPrint ( "  F-  Value of the drag:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift:          ", lift );

    VectorSmall<2> Forces;
    Forces[0] = drag;
    Forces[1] = lift;

    return Forces;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
VectorSmall<2>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeForcesNewTest ( const vector_Type& drag_vector,
																						const vector_Type& lift_vector )
{
    Real drag (0.0);
    Real lift (0.0);

    drag = M_residual->dot(drag_vector);
    lift = M_residual->dot(lift_vector);

    M_Displayer.leaderPrint ( "  F-  Value of the drag:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift:          ", lift );

    VectorSmall<2> Forces;
    Forces[0] = drag;
    Forces[1] = lift;

    return Forces;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
VectorSmall<2>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeForces ( const markerID_Type&  flag,
                                                                               bcHandler_Type& bcHandlerDrag,
                                                                               bcHandler_Type& bcHandlerLift)
{
    // 1) flag che dice in quale componente calcolare la resistenza (0 in x, 1 in y, 2 in z)
    // 2) creare vettore pieno di zeri lungo quanto il vettore soluzione -> vettore soluzione u organizzato come ( ux | uy | uz | p)
    // 3) mettere degli 1 nei nodi corrispondenti al bordo nel blocco del vettore, ovvero:
    //      se la resistenza é in direzione x: vettore da creare ( 1suBody | 0 | 0 | 0 )
    // 4) prodotto scalare del vettore del residuo con questo vettore definito al punto 3)


    vector_Type onesOnBodyDrag(M_localMap, Unique);
    onesOnBodyDrag *= 0;

    vector_Type onesOnBodyLift(M_localMap, Unique);
    onesOnBodyLift *= 0;

    bcManageRhs ( onesOnBodyDrag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerDrag, M_velocityFESpace.feBd(), 1., 0.);
    bcManageRhs ( onesOnBodyLift, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerLift, M_velocityFESpace.feBd(), 1., 0.);

    Real drag (0.0);
    Real lift (0.0);

    drag = M_residual->dot(onesOnBodyDrag);
    lift = M_residual->dot(onesOnBodyLift);

    M_Displayer.leaderPrint ( "  F-  Value of the drag:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift:          ", lift );

    VectorSmall<2> Forces;
    Forces[0] = drag;
    Forces[1] = lift;

    return Forces;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
VectorSmall<2>
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::computeDrag ( const markerID_Type&  flag,
                                                                               bcHandler_Type& bcHandlerDrag,
                                                                               bcHandler_Type& bcHandlerLift,
                                                                               const Real& velocityInfty,
                                                                               const Real& Area)
{
    // 1) flag che dice in quale componente calcolare la resistenza (0 in x, 1 in y, 2 in z)
    // 2) creare vettore pieno di zeri lungo quanto il vettore soluzione -> vettore soluzione u organizzato come ( ux | uy | uz | p)
    // 3) mettere degli 1 nei nodi corrispondenti al bordo nel blocco del vettore, ovvero:
    //      se la resistenza é in direzione x: vettore da creare ( 1suBody | 0 | 0 | 0 )
    // 4) prodotto scalare del vettore del residuo con questo vettore definito al punto 3)


    vector_Type onesOnBodyDrag(M_localMap, Unique);
    onesOnBodyDrag *= 0;

    vector_Type onesOnBodyLift(M_localMap, Unique);
    onesOnBodyLift *= 0;

    bcManageRhs ( onesOnBodyDrag, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerDrag, M_velocityFESpace.feBd(), 1., 0.);
    bcManageRhs ( onesOnBodyLift, *M_velocityFESpace.mesh(), M_velocityFESpace.dof(),  bcHandlerLift, M_velocityFESpace.feBd(), 1., 0.);

    Real drag (0.0);
    Real lift (0.0);

    drag = M_residual->dot(onesOnBodyDrag);
    lift = M_residual->dot(onesOnBodyLift);

    M_Displayer.leaderPrint ( "  F-  Value of the drag:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift:          ", lift );

    M_Displayer.leaderPrint ( "\n\n" );

    drag /= (0.5*M_oseenData->density()*velocityInfty*velocityInfty*Area);
    lift /= (0.5*M_oseenData->density()*velocityInfty*velocityInfty*Area);

    M_Displayer.leaderPrint ( "  F-  Value of the drag coefficient:          ", drag );
    M_Displayer.leaderPrint ( "  F-  Value of the lift coefficient:          ", lift );

    VectorSmall<2> Coefficients;
    Coefficients[0] = drag;
    Coefficients[1] = lift;

    return Coefficients;
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
Real
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::removeMean ( vector_Type& x )
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

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::applyBoundaryConditions ( matrix_Type&       matrix,
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


template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setupPostProc( )
{
    M_postProcessing.reset ( new PostProcessingBoundary<mesh_Type> ( M_velocityFESpace.mesh(),
                                                                     &M_velocityFESpace.feBd(),
                                                                     &M_velocityFESpace.dof(),
                                                                     &M_pressureFESpace.feBd(),
                                                                     &M_pressureFESpace.dof(),
                                                                     M_localMap ) );
}

template<typename MeshType, typename SolverType, typename  MapType , UInt SpaceDim, UInt FieldDim>
void
OseenSolver<MeshType, SolverType, MapType , SpaceDim, FieldDim>::setTolMaxIteration ( const Real& tolerance, const Int& maxIteration )
{
    M_linearSolver->setTolerance ( tolerance );
    M_linearSolver->setMaxNumIterations ( maxIteration );
}

} // namespace LifeV



// #include <lifev/navier_stokes/solver/OseenSolverImplementation.cpp>

#endif // OSEENSOLVER_H
