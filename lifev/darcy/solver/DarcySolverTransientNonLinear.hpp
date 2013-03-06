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
     @brief This file contains a non-linear and (possibly) degenerate permeability
     term Darcy equation using expanded formulation

     @date 05/2010
     @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

     @contributor M. Kern <michel.kern@inria.fr>
     @maintainer M. Kern <michel.kern@inria.fr>
 */

#ifndef _DARCYSOLVERTRANSIENTNONLINEAR_H_
#define _DARCYSOLVERTRANSIENTNONLINEAR_H_ 1

#include <lifev/darcy/solver/DarcySolverNonLinear.hpp>
#include <lifev/darcy/solver/DarcySolverTransient.hpp>

// LifeV namespace.
namespace LifeV
{
//! @class DarcySolverTransientNonLinear This class implements a non-linear transient mixed-hybrid FE Darcy solver.
/*!
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    This class implements a non-linear and transient Darcy solver.
    <br>
    <br>
    The classical time dependant, non-linear, Darcy formulation is a couple of differential equations of first order with
    the unknowns \f$ p \in C^1 (\Omega ) \f$, being the pressure or the primal unknown,
    and \f$ \sigma \in (C^1( \Omega ) )^n \f$, being the Darcy velocity or the flux or the dual unknown,
    such that
    \f[
    \left\{
    \begin{array}{l l l}
    \Lambda^{-1}(t, p) \sigma + \nabla p = f_v(t) & \mathrm{in} & \Omega \times [0, T]\,,  \vspace{0.2cm} \\
    \displaystyle \phi \frac{\partial p}{\partial t} + \nabla \cdot \sigma + \pi(t) p - f(t) = 0 & \mathrm{in} & \Omega \times [0, T]\,,  \vspace{0.2cm} \\
    p = g_D(t)                            & \mathrm{on} & \Gamma_D \times [0, T]\,,\vspace{0.2cm} \\
    \sigma \cdot n + h(t) p = g_R(t)         & \mathrm{on} & \Gamma_R \times [0, T]\,. \vspace{0.2cm} \\
    p(0) = p_0    & \mathrm{in} & \Omega
    \end{array}
    \right.
    \f]
    The data in the system are:
    <ul>
        <li> \f$ p_0 \f$ the primal variable at initial time; </li>
        <li> \f$ \phi \f$ the mass term, which does not depend on time; </li>
        <li> \f$ \Lambda(t, p) \f$ the non linear in \f$ p \f$ permeability tensor; </li>
        <li> \f$ f(t) \f$ the scalar source term; </li>
        <li> \f$ f_v(t) \f$ the vector source term; </li>
        <li> \f$ \pi(t) \f$ the reaction coefficient; </li>
        <li> \f$ \Gamma_D \f$ the subset of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions; </li>
        <li> \f$ g_D(t) \f$ Dirichlet boundary condition data; </li>
        <li> \f$ \Gamma_R \f$ that is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary conditions; </li>
        <li> \f$ h(t) \f$ and \f$ g_R(t) \f$ Robin boundary condition data; </li>
    </ul>
    We suppose that \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and \f$ \Gamma_D \cap \Gamma_R = \emptyset \f$.
    <BR>
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
    Introducing the following bilinear forms, operator and functionals
    \f[
    \begin{array}{l l}
    a(\sigma(t), \tau, p) = \displaystyle \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1}(t, p) \sigma \cdot \tau \,,  &
    b(p, \tau) =  \displaystyle -\sum_{K \in \mathcal{T}_h} \int_K p \nabla \cdot \tau\,, \vspace{0.2cm}\\
    c(\lambda(t), \tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda(t) \tau \cdot n\,,&
    d(p(t), q) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \pi(t) p(t) q \,,\vspace{0.2cm}\\
    h(\lambda(t), \mu) =  \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda(t) \,,&
    m(p(t), v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K \phi p(t) v \,,\vspace{0.2cm}\\
    F(v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f v\,,&
    F_v(\tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f_v(t) \tau \,, \vspace{0.2cm}\\
    G(\mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
    \end{array}
    \f]
    we obtain the Darcy problem in the weak form: find \f$ (\sigma(t), \, p(t), \, \lambda(t)) \in Z \times V \times \Lambda \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma(t), \tau, p(t)) + b(p(t), \tau) + c(\lambda(t), \tau) = F_v(\tau) \,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\partial}{\partial t} m(p(t), v) + b(v, \sigma) - d(p(t), v)= -F(v)\,,                                  & \forall v \in V \,,\vspace{0.2cm}\\
    c(\mu, \sigma(t)) + h(\lambda(t), \mu) = G(\mu) \,,           & \forall \mu \in \Lambda\,.
    \end{array}
    \right.
    \f]
    At the semi-discrete level only space discretization is performed, we introduce the polynomial space, of degree \f$ r \f$, that approximate the finite dimensional
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
    The finite dimensional problem is: find \f$ (\sigma_h,\,, p_h, \, \lambda_h) \in Z_h \times V_h \times \Lambda_h \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h(t), \tau_h, p_h(t)) + b(p_h(t), \tau_h) + c(\lambda_h(t), \tau_h) = F_v(\tau_h) \,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\partial}{\partial t} m(p_h(t),v_h) + b(v_h, \sigma_h(t)) - d(p_h(t), v_h) = -F(v_h)\,,                                  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h(t)) + h(\lambda_h(t(, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To obatin a fully discrete problem we discretize the time derivative via  BDF (backward differentiation formulae) schemes of order \f$ m \f$,
    with \f$ \Delta t \f$ as time step. Choosing
    \f[
    p_h^0 = \prod_{V_h} p_0
    \f]
    we obtain the following system for each \f$ n = 0, \ldots, N \f$, with \f$ N = \frac{T}{\Delta t} \f$
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h^{n+1}, \tau_h, p_h^{n+1}) + b(p_h^{n+1}, \tau_h) + c(\lambda_h^{n+1}, \tau_h) = F_v(\tau_h)\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\alpha_0}{\Delta t} m(p_h^{n+1},v_h) + b(v_h, \sigma_h^{n+1}) - d(p_h^{n+1}, v_h)= -F(v_h) - \frac{1}{\Delta t} m \left( \sum_{i=1}^m \alpha_i p_h^{n+1-i}, v_h \right) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h^{n+1}) + h(\lambda_h^{n+1}, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To solve the non-linearity we use a fixed point scheme based on the relative difference between two consecutive iterations
    of the primal variable. We start from the, user defined, function \f$ p^{n+1,0} \f$ and solve the linearized problem for \f$ k \geq 1 \f$:
    find \f$ (\sigma_h^{n+1,k},\, p_h^{n+1,k}, \, \lambda_h^{n+1,k}) \in Z_h \times V_h \times \Lambda_h \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h^{n+1,k}, \tau_h, p_h^{n+1,k-1}) + b(p_h^{n+1,k}, \tau_h) + c(\lambda_h^{n+1,k}, \tau_h) = F_v(\tau_h)\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    \displaystyle -\frac{\alpha_0}{\Delta t} m(p_h^{n+1,k},v_h) + b(v_h, \sigma_h^{n+1,k}) - d(p_h^{n+1,k}, v_h)= -F(v_h)
    - \frac{1}{\Delta t} m\left( \sum_{i=1}^m \alpha_i p_h^{n+1-i}, v_h \right) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
      c(\mu_h, \sigma_h^{n+1,k}) + h(\lambda_h^{n+1,k}, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    At each iteration,  to solve the linearized problem we use the static condensation procedure, i.e. the unknowns in the discrete
    weak system are not independent and \f$ p_K \f$, \f$\sigma_K \f$ may be written in function of
    \f$ \lambda_K \f$ alone. We introduce the following local matrices
    \f[
    \begin{array}{l l l}
    \left[ A(p_K^{n+1,k-1}) \right]_{ij} =    \displaystyle \int_K \Lambda^{-1}(t^{n+1}, p^{n+1,k-1}) \psi_j \cdot \psi_i \,, &
    \left[ B \right]_{ij} = -  \displaystyle \int_K \phi_j \nabla \cdot \psi_i \,, &
    \left[ C \right]_{ij} =    \displaystyle \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
    \left[ D \right]_{ij} =    \displaystyle \int_K \pi(t) \phi_j \phi_i\,, &
    \left[ H \right]_{ij} =    \displaystyle \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
    \left[ M \right]_{ij} =    \displaystyle \int_{K} \phi \phi_j \phi_i \,, \vspace{0.2cm} \\
    \left[ F \right]_{j}  =    \displaystyle \int_K f \phi_j\,, &
    \left[ F_v \right]_{j}= \displaystyle   \int_K f_v( (n+1) \Delta t) \cdot \psi_j\,,  &
    \left[ G \right]_{j}  =    \displaystyle \int_{\partial K \cap \Gamma_R } g \xi_j\,,
    \end{array}
    \f]
    where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors. <BR>
    The local matrix formulation of the finite dimensional problem is
    \f[
    \left\{
    \begin{array}{l}
    A(p^{n+1,k-1}) \sigma_K^{n+1,k} + B p_K^{n+1,k} + C \lambda_K^{n+1,k} = F_v\,, \vspace{0.2cm} \\
    \displaystyle -\frac{\alpha_0}{\Delta t} M p_K^{n+1,k} +B^T \sigma_K^{n+1,k}  - D p_K^{n+1,k}= -F
    - \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i}  \,,                    \vspace{0.2cm}\\
    C^T \sigma_K^{n+1,k} + H \lambda_K^{n+1,k} = G\,.
    \end{array}
    \right.
    \f]
    Or alternatively
    \f[
    \begin{array}{l l l}
    \left[
    \begin{array}{c c c}
    A(p_K^{n+1,k-1})    & B &  C \vspace{0.2cm} \\
    B^T & \displaystyle -\frac{\alpha_0}{\Delta t} M -D&  0 \vspace{0.2cm} \\
    C^T & 0 & H
    \end{array}
    \right] \, \cdot &
    \left[
    \begin{array}{c}
    \sigma_K^{n+1,k}  \vspace{0.2cm}\\
    p_K^{n+1,k}       \vspace{0.2cm}\\
    \lambda_K^{n+1,K}
    \end{array}
    \right] \, &
    =
    \left[
    \begin{array}{c}
    F_v \vspace{0.2cm}\\
    \displaystyle -F - \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i} \vspace{0.2cm}\\
    G
    \end{array}
    \right]\,.
    \end{array}
    \f]
    Introducing the local hybrid matrix and local hybrid right hand side
    \f[
    \begin{array}{l}
    \displaystyle L_K =  -C^T A^{-1}( p_k^{n+1,k-1}) C + C^T A^{-1}( p_k^{n+1,k-1}) B \left( \frac{\alpha_0}{\Delta t} M + B^T A^{-1}( p_k^{n+1,k-1}) B +
            D \right)^{-1} B^T A^{-1}( p_k^{n+1,k-1}) C + H \,, \vspace{0.2cm} \\
    \displaystyle r_K = G + C^T A^{-1}( p_k^{n+1,k-1}) B \left( \frac{\alpha_0}{\Delta t} + B^T A^{-1}( p_k^{n+1,k-1}) B + D\right)^{-1} F
        + C^T A^{-1}( p_k^{n+1,k-1}) B \left( \frac{\alpha_0}{\Delta t}M + B^T A^{-1}( p_k^{n+1,k-1}) B + D \right)^{-1} B^T A^{-1}( p_k^{n+1,k-1}) F_v
        - C^T A^{-1}( p_k^{n+1,k-1}) F_v \,,
    \end{array}
    \f]
    Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
    \f[
    L \lambda^{n+1,k} = r \,.
    \f]
    We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
    \f[
    \begin{array}{l}
    \displaystyle p_K^{n+1,k} = \left( \frac{\alpha_0}{\Delta t} M + B^T A^{-1}( p_k^{n+1,k-1}) B + D \right)^{-1}
    \left( F + \frac{1}{\Delta t} M \sum_{i=1}^m \alpha_i p_K^{n+1-i} - B^T A^{-1}( p_k^{n+1,k-1}) C \lambda_K^{n+1,k} \right)\,, \vspace{0.2cm} \\
    \displaystyle \sigma_K^{n+1,k} = -A^{-1} ( p_k^{n+1,k-1})( B p_K^{n+1,k} + C \lambda_K^{n+1,k} )  \,.
    \end{array}
    \f]
    @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
    conditions are imposed via BCHandler class.
    @note Example of usage can be found in darcy_nonlinear and darcy_linear.
    Coupled with an hyperbolic solver in impes.
    @todo Insert any scientific publications that use this solver.
*/
template < typename MeshType >
class DarcySolverTransientNonLinear
        :
        public DarcySolverNonLinear < MeshType >,
        public DarcySolverTransient < MeshType >
{

public:

    //! @name Public Types
    //@{

    //! Typedef for mesh template.
    typedef MeshType mesh_Type;

    //! Self typedef.
    typedef DarcySolverTransientNonLinear < mesh_Type > darcySolverTransientNonLinear_Type;

    //! Darcy solver class.
    typedef DarcySolverLinear < mesh_Type > darcySolverLinear_Type;

    //! Darcy non linear solver class.
    typedef DarcySolverNonLinear < mesh_Type > darcySolverNonLinear_Type;

    //! Darcy transient solver class.
    typedef DarcySolverTransient < mesh_Type > darcySolverTransient_Type;

    //! Typedef for the data type.
    typedef typename darcySolverLinear_Type::data_Type data_Type;

    //! Shared pointer for the data type.
    typedef typename darcySolverLinear_Type::dataPtr_Type dataPtr_Type;

    //! Shared pointer to a matrix value function.
    typedef typename darcySolverLinear_Type::matrixFctPtr_Type matrixFctPtr_Type;

    //! Distributed vector.
    typedef typename darcySolverLinear_Type::vector_Type vector_Type;

    //@}

    //! @name Constructors & destructor
    //@{

    //! Constructor for the class.
    DarcySolverTransientNonLinear ();

    //! Virtual destructor.
    virtual ~DarcySolverTransientNonLinear () {};

    //@}

    //! @name Methods
    //@{

    //!  Set up the linear solver, the preconditioner for the linear system and the exporter to save the solution.
    virtual void setup ();

    //! Solve the problem calling the non linear solver.
    virtual void solve ();

    //@}

    //! @name Set methods
    //@{

    //! Set the inverse of diffusion tensor
    /*!
      Set the inverse of diffusion tensor.
      @param invPerm Inverse of the permeability tensor for the problem.
    */
    void setInversePermeability ( const matrixFctPtr_Type& invPerm )
    {
        // Call the set inverse permeability of the non-linear Darcy solver
        darcySolverNonLinear_Type::setInversePermeability ( invPerm );
    }

    //@}

protected:

    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    DarcySolverTransientNonLinear ( const darcySolverTransientNonLinear_Type& );

    //@}

    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    darcySolverTransientNonLinear_Type& operator= ( const darcySolverTransientNonLinear_Type& );

    //@}

    //! @name Protected Methods
    //@{

    //! Update all problem variables
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class.
    */
    virtual void resetVariables ()
    {
        darcySolverNonLinear_Type::resetVariables ();
    }

    //! Compute elementary matrices
    /*!
      Locally update the current finite element for the dual finite element
      space, then compute the Hdiv mass matrix.
      @param iElem Id of the current geometrical element.
      @param elmatMix The local matrix in mixed form.
      @param elmatReactionTerm The local matrix for the reaction term.
    */
    virtual void localMatrixComputation ( const UInt & iElem,
                                          MatrixElemental& elmatMix,
                                          MatrixElemental& elmatReactionTerm )
    {
        darcySolverTransient_Type::localMatrixComputation ( iElem, elmatMix, elmatReactionTerm );
    }

    //! Computes local vectors
    /*!
      Call the Darc_y solver localVectorComputation method and
      compute the additional scalar vector for the time dependent term.
      @param iElem Id of the current geometrical element.
      @param elvecMix The local vector in mixed form.
    */
    virtual void localVectorComputation ( const UInt & iElem,
                                          VectorElemental& elvecMix )
    {
        darcySolverTransient_Type::localVectorComputation ( iElem, elvecMix );
    }

    //@}

}; // class DarcySolverTransientNonLinear

//
// IMPLEMENTATION
//

// ===================================================
// Constructors & Destructor
// ===================================================

// Complete constructor.
template < typename MeshType >
DarcySolverTransientNonLinear < MeshType >::
DarcySolverTransientNonLinear ():
        // Standard Darcy solver constructor.
        darcySolverLinear_Type::DarcySolverLinear (),
        // Non-linear Darcy solver constructor.
        darcySolverNonLinear_Type::DarcySolverNonLinear (),
        // Transient Darcy solver contructor.
        darcySolverTransient_Type::DarcySolverTransient ()
{} // Constructor

// ===================================================
// Public methods
// ===================================================

// Set up the linear solver and the preconditioner.
template < typename MeshType >
void
DarcySolverTransientNonLinear < MeshType >::
setup ()
{
    // Call the DarcySolverLinear setup method for setting up the linear solver and preconditioner.
    darcySolverLinear_Type::setup ();

    // Call the DarcySolverTransient setup method for setting up the time data.
    darcySolverTransient_Type::setupTime ();

    // Call the DarcySolverNonLinear setup method for setting up the non-linear solver.
    darcySolverNonLinear_Type::setupNonLinear ();
} // setup

// Solve the problem.
template < typename MeshType >
void
DarcySolverTransientNonLinear < MeshType >::
solve ()
{
    // Reset the right hand side coming from the time advance scheme.
    this->M_rhsTimeAdvance.reset ( new vector_Type ( this->M_primalField->getFESpace().map() ) );

    // Update the RHS
    this->M_timeAdvance->updateRHSFirstDerivative();

    // Put in M_rhsTimeAdvance the contribution for the right hand side coming
    // from the time scheme, without the time step.
    *(this->M_rhsTimeAdvance) = this->M_timeAdvance->rhsContributionFirstDerivative ();

    // Solve the problem with the fixed point scheme.
    this->fixedPoint ();

} // solve

} // namespace LifeV

#endif //_DARCYSOLVERTRANSIENTNONLINEAR_H_

// -*- mode: c++ -*-
