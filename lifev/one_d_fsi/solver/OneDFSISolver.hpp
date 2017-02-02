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
 *  @brief File containing a solver class for the 1D model.
 *
 *  @version 1.0
 *  @date 01-10-2006
 *  @author Vincent Martin
 *  @author Tiziano Passerini
 *  @author Lucia Mirabella
 *
 *  @version 2.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 01-08-2009
 *
 *  @version 2.1
 *  @date 21-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @contributors Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */



#ifndef OneDFSISolver_H
#define OneDFSISolver_H

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/SolverAmesos.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/one_d_fsi/fem/OneDFSIBCHandler.hpp>
#include <lifev/one_d_fsi/solver/OneDFSIDefinitions.hpp>


namespace LifeV
{

//! OneDFSISolver - Solver class for the 1D model.
/*!
 *  @author Vincent Martin, Tiziano Passerini, Lucia Mirabella, Gilles Fourestey, Cristiano Malossi
 *  @see Equations and networks of 1-D models \cite FormaggiaLamponi2003
 *  @see Geometrical multiscale coupling of 1-D models \cite Malossi2011Algorithms \cite Malossi2011Algorithms1D \cite BonnemainMalossi2012LVAD
 *
 *  <b>EQUATIONS:</b> <BR>
 *  The conservative form of the generic hyperbolic problem is
 *
 *  \f[
 *  \frac{\partial \mathbf U}{\partial t} + \frac{\partial \mathbf F(\mathbf U)}{\partial z} + \mathbf S(\mathbf U) = 0,
 *  \f]
 *
 *  where \f$\mathbf U\f$ are the conservative variables, \f$\mathbf F\f$ the corresponding fluxes,
 *  and \f$\mathbf S\f$ represents the source terms.
 *
 *  <b>NUMERICAL DISCRETIZATION:</b> <BR>
 *  We discretize the problem by a second-order Taylor--Galerkin scheme.
 *  Let us consider the time interval \f$[t^n, t^{n+1}]\f$, for \f$n=0,1,2,\dots\f$, with \f$t^n=n \Delta t$, $\Delta t\f$ being the time step.
 *  Given \f$ \mathbf U_h^n\f$ we compute \f$\mathbf U_h^{n+1}\f$ by using the following scheme
 *
 *  \f[
 *  (\mathbf U^{n+1}_h,\varphi_h) = (\mathbf U^n_h,\varphi_h) + \Delta t \left( \mathbf F(\mathbf U^n_h)-
 *  \displaystyle\frac{\Delta t}{2} \displaystyle\frac{\partial \mathbf F(\mathbf U^n_h)}{\partial \mathbf U}\left(\mathbf S(\mathbf U^n_h)+
 *  \displaystyle\frac{\partial \mathbf F(\mathbf U^n_h)}{\partial z}\right), \displaystyle\frac{\partial \varphi_h}{\partial z}\right) +
 *  \Delta t \left( -\mathbf S(\mathbf U^n_h)+\displaystyle\frac{\Delta t}{2}
 *  \displaystyle\frac{\partial \mathbf S(\mathbf U^n_h)}{\partial \mathbf U}\left( \mathbf S(\mathbf U^n_h)+
 *  \displaystyle\frac{\partial \mathbf F(\mathbf U^n_h)}{\partial z}\right), \varphi_h\right).
 *  \f]
 *
 *  <b>IMPLEMENTATION:</b> <BR>
 *  The implementation of the Taylor-Galerkin scheme is the following:
 *
 *  <CODE>
 *  (Un+1, phi) =                               //! massFactor^{-1} * Un+1  <BR>
 *  (Un, phi)                                   //!            mass * U     <BR>
 *  +dt     * (       Fh(Un), dphi/dz )         //!            grad * F(U)  <BR>
 *  -dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     //! gradDiffFlux(U) * S(U)  <BR>
 *  +dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //!   divDiffSrc(U) * F(U)  <BR>
 *  -dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //!stiffDiffFlux(U) * F(U)  <BR>
 *  -dt     * (       Sh(Un), phi )             //!            mass * S(U)  <BR>
 *  +dt^2/2 * (diffSh(Un) Sh(Un), phi )         //!  massDiffSrc(U) * S(U)  <BR>
 *  </CODE>
 *
 *  Let's define:
 *
 *  \cond \TODO improve doxygen description with latex equation and other features \endcond
 *  <ol>
 *      <li> (\phi_i)_{i in nodes} is the basis of P1 (the "hat" functions)
 *      <li> (1_{i+1/2})_{i+1/2 in elements} is the basis of P0 (constant per element).
 *           The vertices of the element "i+1/2" are the nodes "i" and "i+1".
 *  </ol>
 *
 *  Then:
 *
 *  \cond \TODO improve doxygen description with latex equation and other features \endcond
 *  <ol>
 *      <li> Uh    is in P1 : U = sum_{i in nodes} U_i phi_i
 *      <li> Fh(U) is in P1 : F(U) = sum_{i in nodes} F(U_i) phi_i
 *      <li> diffFh(U) is in P0 : diffFlux(U) = sum_{i+1/2 in elements} 1/2 { dF/dU(U_i) + dF/dU(U_i+1) } 1_{i+1/2}
 *           (means of the two extremal values of the cell)
 *      <li> dF/dz(U) = sum_{i in nodes} F(U_i) d(phi_i)/dz
 *      <li> Sh(U) is in P1 : S(U) = sum_{i in nodes} S(U_i) phi_i
 *      <li> diffSh(U) is in P0 : diffSrc(U) = sum_{i+1/2 in elements} 1/2 { dS/dU(U_i) + dS/dU(U_i+1) } 1_{i+1/2}
 *           (means of the two extremal values of the cell)
 *  </ol>
 *
 *  <b>DEVELOPMENT NOTES:</b> <BR>
 *  The option taken here is to define the different tridiagonal matrix
 *  operators (div, grad, mass, stiff) and reconstruct them at each time
 *  step (as they depend on diffFlux and diffSrc). They are thus rebuilt
 *  at the element level and reassembled.
 *  Afterwards, there remains to do only some tridiagonal matrix vector
 *  products to obtain the right hand side. This procedure might appear a bit memory consuming (there are 18
 *  tridiagonal matrices stored), but it has the advantage of being
 *  very clear. If it is too costly, it should be quite easy to improve
 *  it.
 */
class OneDFSISolver
{
public:

    //! @name Typedef & Enumerator
    //@{

    typedef OneDFSIPhysics                          physics_Type;
    typedef std::shared_ptr< physics_Type >       physicsPtr_Type;

    typedef OneDFSIFlux                             flux_Type;
    typedef std::shared_ptr< flux_Type >          fluxPtr_Type;

    typedef OneDFSISource                           source_Type;
    typedef std::shared_ptr< source_Type >        sourcePtr_Type;

    typedef OneDFSIData                             data_Type;
    typedef data_Type::mesh_Type                    mesh_Type;

    typedef data_Type::container2D_Type             container2D_Type;
    typedef data_Type::scalarVector_Type            scalarVector_Type;
    typedef std::array< scalarVector_Type, 4 >    scalarVectorContainer_Type;

    typedef FESpace< mesh_Type, MapEpetra >         feSpace_Type;
    typedef std::shared_ptr< feSpace_Type >       feSpacePtr_Type;

    typedef Epetra_Comm                             comm_Type;
    typedef std::shared_ptr< comm_Type >          commPtr_Type;

    typedef SolverAmesos                            linearSolver_Type;
    typedef std::shared_ptr< linearSolver_Type >  linearSolverPtr_Type;

    typedef linearSolver_Type::vector_type          vector_Type;
    typedef std::shared_ptr< vector_Type >        vectorPtr_Type;
    typedef std::array< vectorPtr_Type, 2 >       vectorPtrContainer_Type;

    typedef linearSolver_Type::matrix_type          matrix_Type;
    typedef std::shared_ptr<matrix_Type>          matrixPtr_Type;
    typedef std::array<matrixPtr_Type, 4 >        matrixPtrContainer_Type;

    typedef std::map< std::string, vectorPtr_Type > solution_Type;
    typedef std::shared_ptr< solution_Type >      solutionPtr_Type;
    typedef solution_Type::const_iterator           solutionConstIterator_Type;

    typedef OneDFSI::bcLine_Type                    bcLine_Type;
    typedef OneDFSI::bcSide_Type                    bcSide_Type;
    typedef OneDFSI::bcType_Type                    bcType_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
     * Need a call to: \c setCommunicator(), \c setProblem(), \c setFESpace()
     */
    explicit OneDFSISolver();

    //! Destructor
    virtual ~OneDFSISolver() {}

    //@}


    //! @name Methods
    //@{

    //! Build constant matrices (mass and grad)
    void buildConstantMatrices();

    //! Setup the solution using the default FESpace map.
    /*!
     * @param solution solution container
     */
    void setupSolution ( solution_Type& solution )
    {
        setupSolution ( solution, M_feSpacePtr->map() );
    }

    //! Setup the solution using user defined FESpace map.
    /*!
     * @param solution solution container
     * @param map map for initializing the solution vectors
     * @param onlyMainQuantities if true setup only \f$Q\f$, \f$P\f$, and \f$\displaystyle\frac{A}{A^0}-1\f$
     */
    void setupSolution ( solution_Type& solution, const MapEpetra& map, const bool& onlyMainQuantities = false );

    //! Initialize all the variables of the solution to a reference condition with \f$Q=0\f$, \f$A=A^0\f$, and \f$P=P_\mathrm{ext}\f$
    /*!
     * @param solution the solution container
     */
    void initialize ( solution_Type& solution );

    //! Update the Riemann variables.
    /*!
     *  @param solution the solution container is passed with \f$A^n\f$, \f$Q^n\f$ and it is updated with \f$W_1^n\f$, \f$W_2^n\f$
     */
    void computeW1W2 ( solution_Type& solution );

    //! Update the pressure.
    /*!
     *  This method compute the value of the pressure (elastic and if necessary also viscoelastic)
     *  adding it to the solution.
     *  @param solution the solution container is passed with \f$A^n\f$, \f$Q^n\f$, \f$W_1^n\f$, \f$W_2^n\f$ and is updated with \f$P^n\f$
     *  @param timeStep time step
     */
    void computePressure ( solution_Type& solution, const Real& timeStep );

    //! Update the ratio between \f$A\f$ and \f$A^0\f$.
    /*!
     *  @param solution the solution container is passed with \f$A^n\f$, is updated with \f$\displaystyle\frac{A}{A^0}-1\f$
     */
    void computeAreaRatio ( solution_Type& solution );

    //! Compute A from the area ratio: \f$\displaystyle\frac{A}{A^0}-1\f$.
    /*!
     *  @param solution the solution container is passed with \f$\displaystyle\frac{A}{A^0}-1\f$ and is updated with \f$A^n\f$
     */
    void computeArea ( solution_Type& solution );

    //! Compute the right hand side
    /*!
     *  @param solution the solution container
     *  @param timeStep the time step.
     */
    void updateRHS ( const solution_Type& solution, const Real& timeStep );

    //! Update convective term and BC. Then solve the linearized system
    /*!
     * @param bcH the BC handler
     * @param time the time
     * @param timeStep the time step
     */
    void iterate ( OneDFSIBCHandler& bcH, solution_Type& solution, const Real& time, const Real& timeStep );

    //! Apply the viscoelastic flow rate correction.
    /*!
     * To introduce the viscoelastic component of the wall within the formulation, we use an operator-splitting technique,
     * where the flow rate is splitted into two components such that \f$Q = \hat{Q} + \tilde{Q}\f$, where \f$\hat{Q}\f$ is the solution
     * of the pure elastic problem and \f$\tilde{Q}\f$ is the viscoelastic correction.
     * On each time interval \f$[t^n, t^{n+1}]\f$ with \f$n \ge 0\f$, firstly we solve the elastic part for \f$\hat{Q}^{n+1}\f$, using \f$Q^n\f$
     * to compute the contributions in such equation to the right hand side, and then we correct the flow rate by solving the following equation
     * By using the mass conservation equation, we remove the time dependence from the viscoelastic wall term. The resulting problem is
     *
     * \f[
     * \displaystyle\frac{1}{A}\displaystyle\frac{\partial \tilde{Q}}{\partial t} -
     * \displaystyle\frac{\partial}{\partial z}\left(\displaystyle\frac{\gamma}{\rho A^{3/2}}\displaystyle\frac{\partial Q}{\partial z}\right) = 0,
     * \f]
     *
     *  which is closed by a proper set of homogeneous boundary conditions for \f$\tilde{Q}\f$.
     *  The corresponding finite element formulation reads: given \f$(A_h^{n+1},\hat{Q}_h^{n+1})\f$,
     *  find \f$\tilde{Q}_h^{n+1}\f$ such that
     *
     *
     * \f[
     * \left(\displaystyle\frac{\tilde{Q}^{n+1}_h}{A^{n+1}_h},\varphi_h\right) +
     * \Delta t \left(\displaystyle\frac{\gamma}{\rho \left(A^{n+1}_h\right)^{3/2}}\displaystyle\frac{\partial \tilde{Q}^{n+1}_h}{\partial z},
     * \displaystyle\frac{\partial \varphi_h}{\partial z}\right) = \left(\displaystyle\frac{\tilde{Q}^{n}_h}{A^{n+1}_h},\varphi_h\right) - \Delta t \left(\displaystyle\frac{\gamma}{\rho \left(A^{n+1}_h\right)^{3/2}}
     * \displaystyle\frac{\partial \hat{Q}^{n+1}_h}{\partial z},\displaystyle\frac{\partial \varphi_h}{\partial z}\right)+
     * \Delta t \left[\displaystyle\frac{\gamma}{\rho \left(A^{n+1}_h\right)^{3/2}}\displaystyle\frac{\partial\hat{Q}^{n+1}_h}{\partial z}\,\varphi_h\right]^L_0.
     * \f]
     *
     * @param area area
     * @param flowRate flow rate
     * @param timeStep the time step
     * @param bcH the BC handler
     * @param updateSystemMatrix flag for the recomputation of the system matrix
     * @return the viscoelastic flow rate correction \f$\tilde{Q}\f$
     */
    vector_Type viscoelasticFlowRateCorrection ( const vector_Type& newArea, const vector_Type& newElasticFlowRate,
                                                 const vector_Type& oldViscoelasticFlowRate, const Real& timeStep,
                                                 OneDFSIBCHandler& bcHandler, const bool& updateSystemMatrix = true );

    //! CFL computation (correct for constant mesh)
    /*!
     *  @param solution the solution container
     *  @param timeStep the time step
     *  @return CFL
     */
    Real computeCFL ( const solution_Type& solution, const Real& timeStep ) const;

    //! Reset the output files
    /*!
     *  @param solution the solution container
     */
    void resetOutput ( const solution_Type& solution );

    //! Save results on output files
    /*!
     * @param solution solution container
     * @param time solution time
     */
    void postProcess ( const solution_Type& solution, const Real& time );

    //@}


    //! @name Set Methods
    //@{

    //! Set problem classes
    /*!
     * @param physicsPtr pointer to the physics class.
     * @param fluxPtr pointer to the flux class.
     * @param sourcePtr pointer to the source class.
     */
    void setProblem ( const physicsPtr_Type& physicsPtr,
                      const fluxPtr_Type&    fluxPtr,
                      const sourcePtr_Type&  sourcePtr );

    //! Set the communicator
    /*!
     * @param commPtr pointer to the Epetra MPI communicator
     */
    void setCommunicator ( const commPtr_Type& commPtr );

    //! Set the FEspace
    /*!
     * @param feSpacePtr pointer to the FE space
     */
    void setFESpace ( const feSpacePtr_Type& feSpacePtr );

    //! Set the linear solver
    /*!
     * @param linearSolverPtr pointer to the linear solver for the hyperbolic problem
     */
    void setLinearSolver ( const linearSolverPtr_Type& linearSolverPtr );

    //! Set the viscoelastic linear solver
    /*!
     * @param linearViscoelasticSolverPtr pointer to the linear solver for the viscoelastic problem
     */
    void setLinearViscoelasticSolver ( const linearSolverPtr_Type& linearViscoelasticSolverPtr );

    //@}


    //! @name Get Methods
    //@{

    //! Get the physics class
    /*!
     *  @return shared pointer to the physics class.
     */
    const physicsPtr_Type& physics() const
    {
        return M_physicsPtr;
    }

    //! Get the flux class
    /*!
     *  @return shared pointer to the flux class.
     */
    const fluxPtr_Type& flux() const
    {
        return M_fluxPtr;
    }

    //! Get the source class
    /*!
     *  @return shared pointer to the source class.
     */
    const sourcePtr_Type& source() const
    {
        return M_sourcePtr;
    }

    //! Return the ID of the boundary node given a side.
    /*!
     *  @param bcSide Side of the boundary.
     *  @return ID of the boundary node.
     */
    UInt boundaryDOF ( const bcSide_Type& bcSide ) const;

    //! Return the value of a quantity (\f$P\f$, \f$A\f$, \f$Q\f$, \f$W_1\f$, \f$W_2\f$) on a specified boundary.
    /*!
     *  Given a bcType and a bcSide it return the value of the quantity.
     *  @param bcType Type of the asked boundary value.
     *  @param bcSide Side of the boundary.
     *  @return value of the quantity on the specified side.
     */
    Real boundaryValue ( const solution_Type& solution, const bcType_Type& bcType, const bcSide_Type& bcSide ) const;

    //! Return the value of the eigenvalues and eigenvectors on a specified boundary.
    /*!
     *  @param bcSide Side of the boundary.
     *  @param solution solution container.
     *  @param eigenvalues output eigenvalues.
     *  @param leftEigenvector1 output left eigenvector associated to the first eigenvalue.
     *  @param leftEigenvector1 output left eigenvector associated to the second eigenvalue.
     */
    void boundaryEigenValuesEigenVectors ( const bcSide_Type& bcSide, const solution_Type& solution,
                                           container2D_Type& eigenvalues,
                                           container2D_Type& leftEigenvector1,
                                           container2D_Type& leftEigenvector2 );

    //! Get the residual container
    /*!
     * @return System residual container
     */
    const vectorPtrContainer_Type& residual() const
    {
        return M_residual;
    }

    //! Get the system matrix without BC
    /*!
     * @return shared pointer to the system matrix without BC
     */
    const matrixPtr_Type& massMatrix() const
    {
        return M_homogeneousMassMatrixPtr;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Update the P1 flux vector from U: M_fluxi = F_h(Un) i=1,2 (works only for P1Seg elements)
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     */
    void updateFlux ( const solution_Type& solution );

    //! Call _updateFlux and update the P0 derivative of flux vector from U:
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *  M_diffFluxij = dF_h/dU(Un) i,j=1,2
     *  M_diffFluxij(elem) = 1/2 [ dF/dU(U(node1(elem))) + dF/dU(U(node2(elem))) ]
     *
     *  (mean value of the two extremal values of dF/dU)
     *  BEWARE: works only for P1Seg elements
     */
    void updatedFdU ( const solution_Type& solution );

    //! Update the P1 source vector from U: M_sourcei = S_h(Un) i=1,2 (works only for P1Seg elements)
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     */
    void updateSource ( const solution_Type& solution );

    //! Call _updateSource and update the P0 derivative of source vector from U:
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *  M_diffSrcij = dS_h/dU(Un) i,j=1,2
     *  M_diffSrcij(elem) = 1/2 [ dS/dU(U(node1(elem))) + dS/dU(U(node2(elem))) ]
     *
     *  (mean value of the two extremal values of dS/dU)
     *  BEWARE: works only for P1Seg elements
     */
    void updatedSdU ( const solution_Type& solution );

    //! Update the matrices
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *
     *  M_massMatrixDiffSrcij, M_stiffMatrixDiffFluxij
     *  M_gradMatrixDiffFluxij, and M_divMatrixDiffSrcij (i,j=1,2)
     *
     *  from the values of diffFlux(Un) and diffSrc(Un)
     *  that are computed with _updateMatrixCoefficients.
     *
     *  call of  _updateMatrixCoefficients,
     *  _updateMatrixElementalrices and _assemble_matrices.
     */
    void updateMatrices();

    //! Update the element matrices with the current element
    /*!
     * \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     */
    void updateElementalMatrices ( const Real& dFdU, const Real& dSdU );

    //! Assemble the matrices
    /*!
     * \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     */
    void matrixAssemble ( const UInt& ii, const UInt& jj );

    //! Update the matrices to take into account Dirichlet BC.
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *  Modify the matrix to take into account
     *  the Dirichlet boundary conditions
     *  (works for P1Seg and canonic numbering!)
     */
    void applyDirichletBCToMatrix ( matrix_Type& matrix );

    //! Apply the inertial Flux correction:
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *
     *  We use a finite element scheme for the correction term:
     *  given the solution of Taylor-Galerkin scheme, solve
     *  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  ( m / rho ) * ( dQtildeh(n+1)/dz, dphi/dz )  //! stiff * Qtilde(U)
     *  = ( m / rho ) *       ( dQhath(n+1)/dz, dphi/dz )  //! stiff * Qhat(U)
     *
     *  m = rho_w h0 / ( 2 sqrt(pi) sqrt(A0) )
     */
    vector_Type inertialFlowRateCorrection ( const vector_Type& );

    //! Apply the longitudinal Flux correction:
    /*!
     *  \cond \TODO improve doxygen description with latex equation, input/output parameter, etc... \endcond
     *  We use a finite element scheme for the correction term:
     *  given the solution of Taylor-Galerkin scheme, solve
     *  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  = ( 1/Ah(n+1) Qtildeh(n), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  + ( a / rho ) *       ( d3Ahath(n+1)/dz3, phi )  //! mass * d3Ahat(U)/dz
     */
    vector_Type longitudinalFlowRateCorrection();

    //! L2 Projection of the second derivative of Q over P1 space.
    //scalarVector_Type                       _compute_d2Q_dx2( const scalarVector_Type& );

    //@}

    physicsPtr_Type                    M_physicsPtr;
    fluxPtr_Type                       M_fluxPtr;
    sourcePtr_Type                     M_sourcePtr;
    feSpacePtr_Type                    M_feSpacePtr;
    commPtr_Type                       M_commPtr;
    Displayer                          M_displayer;

    std::shared_ptr< MatrixElemental > M_elementalMassMatrixPtr;       //!< element mass matrix
    std::shared_ptr< MatrixElemental > M_elementalStiffnessMatrixPtr;  //!< element stiffness matrix
    std::shared_ptr< MatrixElemental > M_elementalGradientMatrixPtr;   //!< element gradient matrix
    std::shared_ptr< MatrixElemental > M_elementalDivergenceMatrixPtr; //!< element divergence matrix

    //! Right hand sides of the linear system i: "mass * M_Ui = M_rhsi"
    vectorPtrContainer_Type            M_rhs;

    //! Residual of the linear system
    vectorPtrContainer_Type            M_residual;

    //! Flux F(U) (in P1)
    vectorPtrContainer_Type            M_fluxVector;

    //! Source term S (in P1)
    vectorPtrContainer_Type            M_sourceVector;

    //! diffFlux = dF(U)/dU (in P0)
    scalarVectorContainer_Type         M_dFdUVector;

    //! diffSrc = dSource(U)/dU (in P0)
    scalarVectorContainer_Type         M_dSdUVector;

    //! tridiagonal mass matrix
    matrixPtr_Type                     M_homogeneousMassMatrixPtr;

    //! tridiagonal gradient matrix
    matrixPtr_Type                     M_homogeneousGradientMatrixPtr;

    //! tridiagonal mass matrices multiplied by diffSrcij
    matrixPtrContainer_Type            M_dSdUMassMatrixPtr;

    //! tridiagonal stiffness matrices multiplied by diffFluxij
    matrixPtrContainer_Type            M_dFdUStiffnessMatrixPtr;

    //! tridiagonal gradient matrices multiplied by diffFluxij
    matrixPtrContainer_Type            M_dFdUGradientMatrixPtr;

    //! tridiagonal divergence matrices multiplied by diffSrcij
    matrixPtrContainer_Type            M_dSdUDivergenceMatrixPtr;

    //! The linear solver
    linearSolverPtr_Type               M_linearSolverPtr;
    linearSolverPtr_Type               M_linearViscoelasticSolverPtr;

private:

    //! @name Unimplemented Methods
    //@{

    explicit OneDFSISolver ( const OneDFSISolver& solver );

    OneDFSISolver& operator= ( const OneDFSISolver& solver );

    //@}
};

}

#endif // OneDFSISolver_H
