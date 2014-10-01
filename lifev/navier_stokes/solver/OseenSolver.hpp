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

#include <lifev/navier_stokes/solver/StabilizationSUPG.hpp>
#include <lifev/navier_stokes/solver/StabilizationSUPGVMS.hpp>
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

    //! Compute the energy of the fluid
    /*!
        @return energy
     */
    Real energy ( );

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
    boost::shared_ptr<StabilizationSUPG<mesh_Type, map_Type, SpaceDim > > M_supgStabilization;
    
    // SUPGVMS stabilization
    boost::shared_ptr<StabilizationSUPGVMS<mesh_Type, map_Type, SpaceDim > > M_supgVmsStabilization;

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

} // namespace LifeV

#include <lifev/navier_stokes/solver/OseenSolverImplementation.cpp>

#endif // OSEENSOLVER_H
