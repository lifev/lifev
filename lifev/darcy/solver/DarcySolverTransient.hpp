//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains a time dependent Darcy equation solver class

     @date 05/2010
     @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

     @contributor M. Kern <michel.kern@inria.fr>
     @maintainer M. Kern <michel.kern@inria.fr>
 */


#ifndef _DARCYSOLVERTRANSIENT_H_
#define _DARCYSOLVERTRANSIENT_H_ 1

#include <lifev/core/fem/TimeAdvanceBDF.hpp>

#include <lifev/darcy/solver/DarcySolverLinear.hpp>

// LifeV namespace.
namespace LifeV
{
//!  @class DarcySolverTransient This class implements a mixed-hybrid FE Darcy solver for transient problems

/*!
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    This class implements a transient Darcy solver in the finite time interval \f$ [0, T] \f$.
    <br>
    The classical time dependent Darcy formulation is a couple of differential equations of first order with
    the unknowns \f$ p \in C^1(0, T; C^1 (\Omega )) \f$, being the pressure or the primal unknown,
    and \f$ \sigma \in C^0(0, T; (C^1( \Omega ) )^n) \f$, being the Darcy velocity or the flux or the dual unknown,
    such that
    \f[
    \left\{
    \begin{array}{l l l}
    \Lambda^{-1}(t) \sigma + \nabla p = f_v(t)                      & \mathrm{in} & \Omega   \times [0, T]\,, \vspace{0.2cm} \\
    \displaystyle \phi \frac{\partial p}{\partial t} + \nabla \cdot \sigma + \pi(t) p - f(t) = 0  & \mathrm{in} & \Omega   \times [0, T]\,, \vspace{0.2cm} \\
    p = g_D(t)                                                      & \mathrm{on} & \Gamma_D \times [0, T]\,, \vspace{0.2cm} \\
    \sigma \cdot n + h(t) p = g_R(t)                                & \mathrm{on} & \Gamma_R \times [0, T]\,. \vspace{0.2cm} \\
    p(0) = p_0                                                      & \mathrm{in} & \Omega
    \end{array}
    \right.
    \f]
    The data in the system are:
    <ul>
        <li> \f$ p_0 \f$ the primal variable at initial time; </li>
        <li> \f$ \phi \f$ the mass term, which does not depend on time; </li>
        <li> \f$ \Lambda(t) \f$ the permeability tensor; </li>
        <li> \f$ f(t) \f$ the scalar source term; </li>
        <li> \f$ f_v(t) \f$ the vector source term; </li>
        <li> \f$ \pi(t) \f$ the reaction coefficient; </li>
        <li> \f$ \Gamma_D \f$ the subset of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions; </li>
        <li> \f$ g_D(t) \f$ Dirichlet boundary condition data; </li>
        <li> \f$ \Gamma_R \f$ that is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary conditions; </li>
        <li> \f$ h(t) \f$ and \f$ g_R(t) \f$ Robin boundary condition data; </li>
    </ul>
    We suppose that \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and \f$ \Gamma_D \cap \Gamma_R = \emptyset \f$.
    <br>
    Using the hybridization procedure, and introducing a new variable, we may split the problem from
    the entire domain to problems in each elements of the triangulation \f$ \mathcal{T}_h \f$, then we may write
    the weak formulation for the dual problem.
    The new variable is the hybrid unknown \f$ \lambda \f$, e.g. "the trace" of pressure of the primal unknown,
    it is the Lagrange multipliers forcing the continuity of the flux across each faces or edges in the triangulation.
    We introduce the following functional spaces, the first is the space of the primal variable, the third for
    the dual variable, the fourth for the hybrid variable.
    \f[
    \begin{array}{l}
    V = L^2 (\Omega ) \,,\vspace{0.2cm} \\
    H(div, K) = \displaystyle \left\{ \tau \in ( L^2(K))^n : \nabla \cdot \tau \in L^2(K)\right\}\,, \vspace{0.2cm}\\
    Z = \displaystyle \left\{ \tau \in L^2(\Omega) : \tau\vert_K \in H(div, K) \, \forall K \in  \mathcal{T}_h \right\}\,, \vspace{0.2cm}\\
    \Lambda = \displaystyle \left\{ \lambda \in \prod_{K \in \mathcal{T}_h} H^{1/2} (\partial K): \lambda_K = \lambda_{K'} \,\, \mathrm{on} \,\, e_{K-K'} \, \forall K \in \mathcal{T}_h,\, \lambda = g_D \,\, \mathrm{on} \,\, \Gamma_D \right\}\,.
    \end{array}
    \f]
    Introducing the following bilinear forms and functionals at fixed time \f$ t \in [0, T] \f$
    \f[
    \begin{array}{l l}
    a(\sigma(t), \tau) = \displaystyle \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1}(t) \sigma(t) \cdot \tau \,,  &
    b(p(t), \tau) = \displaystyle -\sum_{K \in \mathcal{T}_h} \int_K p(t) \nabla \cdot \tau\,, \vspace{0.2cm}\\
    c(\lambda(t), \tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda(t) \tau \cdot n\,,&
    d(p(t), v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \pi(t) p(t) v \,,\vspace{0.2cm}\\
    h(\lambda(t), \mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h(t) \mu \lambda(t) \,,&
    m(p(t), v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K \phi p(t) v \,,\vspace{0.2cm}\\
    F(v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f(t) v\,,&
    F_v(\tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f_v(t) \tau \,, \vspace{0.2cm}\\
    G(\mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g(t) \mu\,,
    \end{array}
    \f]
    we obtain the Darcy problem in the weak form: find \f$ (\sigma(t), \, p(t), \, \lambda(t)) \in Z \times V \times \Lambda \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma(t), \tau) + b(p(t), \tau) + c(\lambda(t), \tau) = F_v(\tau) \,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\partial }{\partial t} m(p(t), v) + b(v, \sigma(t)) - d(p(t), v) = -F(v)\,, & \forall v \in V \,,\vspace{0.2cm}\\
    c(\mu, \sigma(t)) + h(\lambda(t), \mu) = G(\mu) \,,         & \forall \mu \in \Lambda\,.
    \end{array}
    \right.
    \f]
    At semi-discrete level. i.e. only space discretization is performed,
    we introduce the polynomial space, of degree \f$ r \f$, that approximate the finite dimensional
    spaces introduced above \f$ V_h \subset V \f$, \f$ Z_h \subset Z \f$ and \f$\Lambda_h \subset \Lambda \f$
    \f[
    \begin{array}{l}
    V_h = \displaystyle \left\{ v_h \in V: v_h|_K \in P_r (K)\, \forall K \in \mathcal{T}_h \right\}\,, \vspace{0.2cm}\\
    Z_h = \displaystyle \left\{ \tau_h \in Z: \tau_h|_K \in RT_r(K) \, \forall K \in \mathcal{T}_h \right\} \,, \vspace{0.2cm} \\
    \Lambda_h = \displaystyle \left\{ \lambda_h \in \Lambda: \lambda_h|_{\partial K} \in R_r( \partial K ) \, \forall K \in \mathcal{T}_h\right\}\,,
    \end{array}
    \f]
    where \f$ P_r(K) \f$ is the space of polynomial of degree \f$ r \f$ in the element \f$ K \f$, \f$ RT_r(K) \f$ is the space of
    polynomial of Raviart-Thomas of degrees \f$ r \f$ in the element \f$ K \f$ and \f$ R_r(\partial K) \f$ is the space of polynomial of
    degree \f$ r \f$ definite on each of the boundary of the element \f$ K \f$ and discontinuous from one edge to the other.
    <br>
    The finite dimensional problem is: find \f$ (\sigma_h(t),\, p_h(t), \, \lambda_h(t)) \in Z_h \times V_h \times \Lambda_h \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h(t), \tau_h) + b(p_h(t), \tau_h) + c(\lambda_h(t), \tau_h) = F_v(\tau_h) \,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\partial }{\partial t} m(p_h(t),v_h) + b(v_h, \sigma_h(t)) - d(p_h(t), v_h) = -F(v_h)\,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h(t)) + h(\lambda_h(t), \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To obatin a fully discrete problem we discretize the time derivative via BDF (backward differentiation formulae) schemes of order \f$ m \f$,
    with \f$ \Delta t \f$ as time step. Choosing
    \f[
    p_h^0 = \prod_{V_h} p_0
    \f]
    we obtain the following system for each \f$ n = 0, \ldots, N \f$, with \f$ N = \frac{T}{\Delta t} \f$
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h^{n+1}, \tau_h) + b(p_h^{n+1}, \tau_h) + c(\lambda_h^{n+1}, \tau_h) = F_v(\tau_h)\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\alpha_0}{\Delta t} m(p_h^{n+1},v_h) + b(v_h, \sigma_h^{n+1}) - d(p_h^{n+1}, v_h)= -F(v_h)
                  - \frac{1}{\Delta t} m\left( \sum_{i=1}^m \alpha_i p_h^{n+1-i}, v_h \right) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h^{n+1}) + h(\lambda_h^{n+1}, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To solve the problem we use the static condensation procedure, i.e. the unknowns in the discrete
    weak system are not independent and \f$ p_K \f$, \f$\sigma_K \f$ may be written in function of
    \f$ \lambda_K \f$ alone. We introduce the following local matrices
    \f[
    \begin{array}{l l l}
    \left[ A \right]_{ij} = \displaystyle   \int_K \Lambda^{-1}( (n+1) \Delta t) \psi_j \cdot \psi_i \,, &
    \left[ B \right]_{ij} = \displaystyle - \int_K \phi_j \nabla \cdot \psi_i \,, &
    \left[ C \right]_{ij} = \displaystyle   \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
    \left[ D \right]_{ij} = \displaystyle   \int_K \pi(t) \phi_j \phi_i\,, &
    \left[ H \right]_{ij} = \displaystyle   \int_{\partial K \cap \Gamma_R} h( (n+1) \Delta t) \xi_i \xi_j\,, &
    \left[ M \right]_{ij} = \displaystyle   \int_{K} \phi \phi_j \phi_i \,, \vspace{0.2cm} \\
    \left[ F \right]_{j}  = \displaystyle   \int_K f( (n+1) \Delta t) \phi_j\,,&
    \left[ F_v \right]_{j}= \displaystyle   \int_K f_v( (n+1) \Delta t) \cdot \psi_j\,,  &
    \left[ G \right]_{j}  = \displaystyle   \int_{\partial K \cap \Gamma_R } g( (n+1) \Delta t) \xi_j\,,
    \end{array}
    \f]
    where we avoid to write the dependence on the triangle \f$ K \f$ and on the current time step \f$ n+1 \f$ in all the matrices and vectors.
    <br>
    bre local matrix formulation of the finite dimensional problem is
    \f[
    \left\{
    \begin{array}{l}
    A \sigma_K^{n+1} + B p_K^{n+1} + C \lambda_K^{n+1} = F_v\,, \vspace{0.2cm} \\
    \displaystyle -\frac{\alpha_0}{\Delta t} M p_K^{n+1} + B^T \sigma_K^{n+1} - D p_K^{n+1} = -F
    - \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i} \,,\vspace{0.2cm}\\
    C^T \sigma_K^{n+1} + H \lambda_K^{n+1} = G\,.
    \end{array}
    \right.
    \f]
    Or alternatively
    \f[
    \begin{array}{l l l}
    \left[
    \begin{array}{c c c}
    A   & B                    &  C \vspace{0.2cm} \\
    B^T & \displaystyle -\frac{\alpha_0}{\Delta t} M - D &  0 \vspace{0.2cm} \\
    C^T & 0                    & H
    \end{array}
    \right] \, \cdot &
    \left[
    \begin{array}{c}
    \sigma_K^{n+1}  \vspace{0.2cm}\\
    p_K^{n+1}       \vspace{0.2cm}\\
    \lambda_K^{n+1}
    \end{array}
    \right] \, &
    =
    \left[
    \begin{array}{c}
    F_v \vspace{0.2cm}\\
    \displaystyle - F - \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i} \vspace{0.2cm}\\
    G
    \end{array}
    \right]\,.
    \end{array}
    \f]
    Introducing the local hybrid matrix and local hybrid right hand side
    \f[
    \begin{array}{l}
    \displaystyle L_K = -C^T A^{-1} C + C^T A^{-1} B \left( \frac{\alpha_0}{\Delta t} M + B^T A^{-1} B  + D \right)^{-1} B^T A^{-1} C + H \,, \vspace{0.2cm} \\
    \displaystyle r_K = G + C^T A^{-1} B \left( \frac{\alpha_0}{\Delta t} M + B^T A^{-1} B + D \right)^{-1}
     \left( F + \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i} \right)
     + C^T A^{-1} B \left( \frac{\alpha_0}{\Delta t}M + B^T A^{-1} B + D \right)^{-1} B^T A^{-1} F_v - C^T A^{-1} F_v \,,
    \end{array}
    \f]
    Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
    \f[
    L \lambda^{n+1} = r \,.
    \f]
    We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
    \f[
    \begin{array}{l}
    \displaystyle p_K^{n+1} = \left( \frac{\alpha_0}{\Delta t} M + B^T A^{-1} B + D\right)^{-1} \left( F + \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i}
    - B^T A^{-1} C \lambda_K^{n+1} \right)\,, \vspace{0.2cm} \\
    \sigma_K^{n+1} = -A^{-1} \left( B p_K^{n+1} + C \lambda_K^{n+1} \right) \,.
    \end{array}
    \f]
    @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$,
    because all the boundary conditions are imposed via BCHandler class.
    @note The initial time is not fixed at zero.
    @note Example of usage can be found in darcy_nonlinear and darcy_linear.
    Coupled with an hyperbolic solver in impes.
    @todo Insert any scientific publications that use this solver.
*/
template < typename MeshType >
class DarcySolverTransient :
    virtual public DarcySolverLinear < MeshType >
{

public:

    //! @name Public Types
    //@{

    //! Typedef for mesh template.
    typedef MeshType mesh_Type;

    //! Self typedef.
    typedef DarcySolverTransient < mesh_Type > darcySolverTransient_Type;

    //! Darcy solver class.
    typedef DarcySolverLinear < mesh_Type > darcySolverLinear_Type;

    //! Typedef for the data type.
    typedef typename darcySolverLinear_Type::data_Type data_Type;

    //! Shared pointer for the data type.
    typedef typename darcySolverLinear_Type::dataPtr_Type dataPtr_Type;

    //! Shared pointer to a scalar value function.
    typedef typename darcySolverLinear_Type::scalarFctPtr_Type scalarFctPtr_Type;

    //! Shared pointer to a distributed vector.
    typedef typename darcySolverLinear_Type::vectorPtr_Type vectorPtr_Type;

    //! Distributed vector.
    typedef typename darcySolverLinear_Type::vector_Type vector_Type;

    //! Shared pointer to a scalar field.
    typedef typename darcySolverLinear_Type::scalarFieldPtr_Type scalarFieldPtr_Type;

    //! Scalar field.
    typedef typename darcySolverLinear_Type::scalarField_Type scalarField_Type;

    //! Container of matrix elemental.
    typedef std::vector< MatrixElemental > matrixElementalContainer_Type;

    //! Time advance scheme.
    typedef TimeAdvanceBDF < vector_Type > timeAdvance_Type;

    //! Shared pointer to a time advance scheme.
    typedef boost::shared_ptr < timeAdvance_Type > timeAdvancePtr_Type;

    //@}

    //! @name Constructors & destructor
    //@{

    //! Constructor for the class.
    DarcySolverTransient ();

    //! Virtual destructor.
    virtual ~DarcySolverTransient () {};

    //@}

    // Methods
    //! @name methods
    //@{

    //! Set up the linear solver and the preconditioner for the linear system.
    virtual void setup ();

    //! Solve the problem.
    virtual void solve ();

    //@}

    // Update and set Methods.
    //! @name Update and set methods
    //@{

    //! Initialize primal solution
    /*!
      Set the initial value function for the primal variable and compute the primal
      variable at step zero.
      @param primalInitial The primal initial function.
    */
    void setInitialPrimal ( const scalarFctPtr_Type& primalInitialFct );

    //! Set mass matrix
    /*!
      Set the mass term, the default source term is the function one.
      By defaul it does not depend on time.
      @param mass Mass term for the problem.
    */
    void setMass ( const scalarFctPtr_Type& massFct );

    //@}

protected:

    // Methods
    //! @name Protected methods
    //@{

    //! Compute element matrices
    /*!
      Call the Darcy solver localMatrixComputation method and
      compute the mass matrix for the time dependent term.
      @param iElem Id of the current geometrical element.
      @param elmatMix The local matrix in mixed form.
      @param elmatReactionTerm The local matrix for the reaction term.
    */
    virtual void localMatrixComputation ( const UInt& iElem,
                                          MatrixElemental& elmatMix,
                                          MatrixElemental& elmatReactionTerm );

    //! Computes local vectors
    /*!
      Call the Darc_y solver localVectorComputation method and
      compute the additional scalar vector for the time dependent term.
      @param iElem Id of the current geometrical element.
      @param elvecMix The local vector in mixed form.
    */
    virtual void localVectorComputation ( const UInt& iElem,
                                          VectorElemental& elvecMix );

    //! Do some computation after the calculation of the primal and dual variable.
    /*!
      Save into the time advance scheme the computed primal solution.
    */
    virtual void postComputePrimalAndDual ()
    {
        // Update the solution in the time advance.
        M_timeAdvance->shiftRight ( this->M_primalField->getVector() );
    }

    //! Setup the time data.
    void setupTime ();

    //@}

    // Time advance stuff
    //! @name Time advance stuff
    //@{

    //! Time advance.
    timeAdvancePtr_Type M_timeAdvance;

    //! Right hand side coming from the time advance scheme.
    vectorPtr_Type M_rhsTimeAdvance;

    //@}

private:

    // Private Constructors
    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    DarcySolverTransient ( const darcySolverTransient_Type& );

    //@}

    // Private Operators
    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    darcySolverTransient_Type& operator= ( const darcySolverTransient_Type& );

    //@}

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //! Initial time primal variable.
    scalarFctPtr_Type M_primalFieldInitialFct;

    //! Mass function, it does not depend on time.
    scalarFctPtr_Type M_massFct;

    //@}

    // Time advance stuff
    //! @name Time advance stuff
    //@{

    //! Local mass matrices.
    matrixElementalContainer_Type M_localMassMatrix;

    //@}

    // Algebraic stuff
    //! @name Algebraic stuff
    //@{

    //! Boolean that indicates if the preconditioner is re-used or not.
    bool M_reusePrec;

    //! Boolean that indicates if the matrix is updated for the current iteration.
    bool M_updated;

    //! Interger storing the max number of solver iteration with preconditioner recomputing.
    UInt M_maxIterSolver;

    //! Boolean that indicates if the matrix is recomputed for the current iteration.
    bool M_recomputeMatrix;

    //@}

}; // class DarcySolverTransient

// IMPLEMENTATION

// Complete constructor.
template < typename MeshType >
DarcySolverTransient < MeshType >::
DarcySolverTransient () :
    // Standard Darcy solver constructor.
    darcySolverLinear_Type::DarcySolverLinear (),
    // Time advance data.
    M_timeAdvance                 ( new timeAdvance_Type ),
    // Linear solver.
    M_reusePrec                   ( false ),
    M_updated                     ( false ),
    M_maxIterSolver               ( static_cast<UInt> (0) ),
    M_recomputeMatrix             ( false )
{
} // Constructor

// ===========================================================================================
// Public methods
// ==========================================================================================

// Set up the linear solver and the preconditioner.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
setup ()
{

    // Call the DarcySolverLinear setup method for setting up the linear solver.
    darcySolverLinear_Type::setup ();

    // Setup the time data
    setupTime ();

} // setup

// Set up the linear solver and the preconditioner.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
setupTime ()
{

    const typename darcySolverLinear_Type::data_Type::data_Type& dataFile = * ( this->M_data->dataFilePtr() );

    // Set if the preconditioner is re-used.
    M_reusePrec = dataFile ( ( this->M_data->section() + "/solver/reuse" ).data(), false );

    // Set if the preconditioner is reused or not.
    this->M_linearSolver.setReusePreconditioner ( M_reusePrec );

    // Set the max number of iteration to mantein the same preconditioner.
    M_maxIterSolver = dataFile ( ( this->M_data->section() + "/solver/max_iter_reuse" ).data(), static_cast<Int> (0) );

    // Set up the time advance.
    M_timeAdvance->setup ( this->M_data->dataTimeAdvancePtr()->orderBDF(), 1 );

} // setupTime

// Set the inital value
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
setInitialPrimal ( const scalarFctPtr_Type& primalInitialFct )
{
    // Set the initial value function.
    M_primalFieldInitialFct = primalInitialFct;

    // Create the interpolated initial condition for the primal variable.
    scalarField_Type primalInitialField ( this->M_primalField->getFESpacePtr(),
                                          this->M_primalField->getVector().mapType() );

    // Interpolate the primal initial value.
    M_primalFieldInitialFct->interpolate ( primalInitialField,
                                           this->M_data->dataTimePtr()->initialTime() );

    // Set the initial condition for the time advance method.
    M_timeAdvance->setInitialCondition ( primalInitialField.getVector() );

} // setInitialPrimal

// ===========================================================================================
// Protected methods
// ==========================================================================================

// Perform all the operations before doing the loop on volume elements.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
solve ()
{
    // Reset the right hand side coming from the time advance scheme.
    M_rhsTimeAdvance.reset ( new vector_Type ( this->M_primalField->getFESpace().map() ) );

    // Update the RHS
    M_timeAdvance->updateRHSFirstDerivative ();

    // Put in M_rhsTimeAdvance the contribution for the right hand side coming
    // from the time scheme, without the time step.
    *M_rhsTimeAdvance = M_timeAdvance->rhsContributionFirstDerivative ();

    // Solve the problem
    darcySolverLinear_Type::solve ();

} // solve

// Set and compute the mass matrix.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
setMass ( const scalarFctPtr_Type& massFct )
{
    // Save the mass function.
    M_massFct = massFct;

    // The total number of elements in the mesh.
    const UInt meshNumberOfElements = this->M_primalField->getFESpace().mesh()->numElements();

    // Element of current id.
    typename mesh_Type::element_Type element;

    // Reseve the memory for the mass matrices.
    M_localMassMatrix.reserve ( meshNumberOfElements );

    // Elemental matrix to store the current local mass matrix.
    MatrixElemental localMassMatrix ( this->M_primalField->getFESpace().refFE().nbDof(), 1, 1 );

    // The coordinate of the barycenter of local element.
    Vector3D barycenter;

    //! Loop on all the volume elements.
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {
        // Element of current ID.
        element = this->M_primalField->getFESpace().mesh()->element ( iElem );

        // Update the current element of ID iElem for the primal variable.
        this->M_primalField->getFESpace().fe().update ( element, UPDATE_QUAD_NODES | UPDATE_WDET );

        // Clean the local mass matrix before the computation.
        localMassMatrix.zero ();

        // Get the coordinates of the barycenter of the current element of ID iElem.
        this->M_primalField->getFESpace().fe().barycenter ( barycenter[0], barycenter[1], barycenter[2] );

        // Computes the value for the mass term.
        const Real massValue = M_massFct->eval ( iElem, barycenter );

        // Compute the mass matrix for the primal variable.
        AssemblyElemental::mass ( massValue, localMassMatrix, this->M_primalField->getFESpace().fe(), 0, 0);

        // Save the computed mass matrix.
        M_localMassMatrix.push_back ( localMassMatrix );
    }

} // setMass

// Call the Darcy solver localMatrixComputation method and compute the mass matrix for the time dependent term.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
localMatrixComputation ( const UInt& iElem, MatrixElemental& elmatMix,
                         MatrixElemental& elmatReactionTerm )
{
    // Call the Darcy solver local matrix computation.
    darcySolverLinear_Type::localMatrixComputation ( iElem, elmatMix, elmatReactionTerm );

    // Check if the mass function is set or not.
    ASSERT ( M_massFct.get(), "DarcySolverTransient : mass function not set." );

    // Elemental matrix to store the current local mass matrix.
    MatrixElemental localMassMatrix = M_localMassMatrix [ iElem ];

    // Multiply the mass matrix by the time step and the time scheme coefficient.
    localMassMatrix *= M_timeAdvance->coefficientFirstDerivative ( 0 ) /
                       this->M_data->dataTimePtr()->timeStep();

    /* Store in the reaction matrix the mass matrix,
       divided by the time step and multiplied by the time scheme coefficient. */
    elmatReactionTerm.mat() += localMassMatrix.mat();

} // localMatrixComputation

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template < typename MeshType >
void
DarcySolverTransient < MeshType >::
localVectorComputation ( const UInt& iElem, VectorElemental& elvecMix )
{
    /* Call the Darcy solver localVectorComputation to update the finite elements
       spaces and to compute the scalar and vector source term. */
    darcySolverLinear_Type::localVectorComputation ( iElem, elvecMix );

    // Create the local contribution of the time advance right hand side.
    VectorElemental localRhsTimeAdvance ( this->M_primalField->getFESpace().refFE().nbDof(), 1 );

    // Clean the local contribution before the computations.
    localRhsTimeAdvance.zero ();

    // Extract from the time advance right hand side the current contribution.
    extract_vec ( *M_rhsTimeAdvance,
                  localRhsTimeAdvance,
                  this->M_primalField->getFESpace().refFE(),
                  this->M_primalField->getFESpace().dof(),
                  this->M_primalField->getFESpace().fe().currentLocalId(), 0 );

    // Check if the mass function is setted or not.
    ASSERT ( M_massFct.get(), "DarcySolverTransient : mass function not setted." );

    // Extract the current mass matrix.
    MatrixElemental::matrix_type localMassMatrix ( M_localMassMatrix [ iElem ].mat() );

    // Divide the mass matrix by the time step.
    localMassMatrix /= this->M_data->dataTimePtr()->timeStep();

    // Add to the local scalar source term the time advance term.
    elvecMix.block ( 1 ) += localMassMatrix * localRhsTimeAdvance.vec();

} // localVectorComputation

} // namespace LifeV

#endif // _DARCYSOLVERTRANSIENT_H_

// -*- mode: c++ -*-
