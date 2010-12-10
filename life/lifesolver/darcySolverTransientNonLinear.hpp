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

#include <life/lifesolver/darcySolverNonLinear.hpp>
#include <life/lifesolver/darcySolverTransient.hpp>

// LifeV namespace.
namespace LifeV
{
//  @class DarcySolverTransientNonLinear This class implements a non-linear transient mixed-hybrid FE Darcy solver

/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

  This class implements a nolinear (possibly degenerate) Darcy solver. <br>
  The classical time dependant, non-linmar, Darcy formulation is a couple of differential equations of first order with
  the unknowns \f$ p \in C^1 (\Omega ) \f$, being the pressure or the primal unknown,
  and \f$ \sigma \in (C^1( \Omega ) )^n \f$, being the Darcy velocity or the flux or the dual unknown,
  such that
  \f[
  \left\{
  \begin{array}{l l l}
  \Lambda^{-1}(t, p) \sigma + \nabla p = 0 & \mathrm{in} & \Omega \times [0, T]\,,  \vspace{0.2cm} \\
  \displaystyle \frac{\partial p}{\partial t} + \nabla \cdot \sigma - f(t) = 0        & \mathrm{in} & \Omega \times [0, T]\,,  \vspace{0.2cm} \\
  p = g_D(t)                            & \mathrm{on} & \Gamma_D \times [0, T]\,,\vspace{0.2cm} \\
  \sigma \cdot n + h(t) p = g_R(t)         & \mathrm{on} & \Gamma_R \times [0, T]\,. \vspace{0.2cm} \\
  p(0) = p_0    & \mathrm{in} & \Omega
  \end{array}
  \right.
  \f]
  Where \f$ \Lambda(t, p) \f$ is the permeability tensor that depends on \f$ p \f$ and \f$ t \f$, \f$ f \f$ is the source term, \f$ \Gamma_D \f$ is the subset
  of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions with datum \f$ g_D \f$ and \f$ \Gamma_R \f$
  is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary conditions with data \f$ h \f$ and
  \f$ g_R \f$. We suppose that \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and
  \f$ \Gamma_D \cap \Gamma_R = \emptyset \f$.
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
  H(div, K) = \left\{ \tau \in ( L^2(K))^n : \nabla \cdot \tau \in L^2(K)\right\}\,, \vspace{0.2cm}\\
  Z = \left\{ \tau \in L^2(\Omega) : \tau\vert_K \in H(div, K) \, \forall K \in  \mathcal{T}_h \right\}\,, \vspace{0.2cm}\\
  \Lambda = \left\{ \lambda \in \prod_{K \in \mathcal{T}_h} H^{1/2} (\partial K): \lambda_K = \lambda_{K'} \,\, \mathrm{on} \,\, e_{K-K'} \, \forall K \in \mathcal{T}_h,\, \lambda = g_D \,\, \mathrm{on} \,\, \Gamma_D \right\}\,.
  \end{array}
  \f]
  Introducing the following bilinear forms and functionals
  \f[
  \begin{array}{l l}
  a(\sigma(t), \tau, p) = \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1}(t, p) \sigma \cdot \tau \,,  &
  b(p, \tau) = -\sum_{K \in \mathcal{T}_h} \int_K p \nabla \cdot \tau\,, \vspace{0.2cm}\\
  c(\lambda(t), \tau) = \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda(t) \tau \cdot n\,,&
  h(\lambda(t), \mu) = \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda(t) \,,\vspace{0.2cm}\\
  m(p(t), v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K p(t) v \,,\vspace{0.2cm}\\
  F(v) = \sum_{K \in \mathcal{T}_h} \int_K f v\,,&
  G(\mu) =\sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
  \end{array}
  \f]
  we obtain the Darcy problem in the weak form: find \f$ (\sigma(t), \, p(t), \, \lambda(t)) \in Z \times V \times \Lambda \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma(t), \tau, p(t)) + b(p(t), \tau) + c(\lambda(t), \tau) = 0\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
  \displaystyle -\frac{\partial }{\partial t} m(p(t), v) + b(v, \sigma) = -F(v)\,,                                  & \forall v \in V \,,\vspace{0.2cm}\\
  -c(\mu, \sigma(t)) - h(\lambda(t), \mu) = - G(\mu) \,,           & \forall \mu \in \Lambda\,.
  \end{array}
  \right.
  \f]
  At the semi-discrete level ie only space discretization is performes, we introduce the polynomial space, of degree \f$ r \f$, that approximate the finite dimensional
  spaces introduced above \f$ V_h \subset V \f$, \f$ Z_h \subset Z \f$ and \f$\Lambda_h \subset \Lambda \f$
  \f[
  \begin{array}{l}
  V_h = \left\{ v_h \in V: v_h|_K \in P_r (K)\, \forall K \in \mathcal{T}_h \right\}\,, \vspace{0.2cm}\\
  Z_h = \left\{ \tau_h \in Z: \tau_h|_K \in RT_r(K) \, \forall K \in \mathcal{T}_h \right\} \,, \vspace{0.2cm} \\
  \Lambda_h = \left\{ \lambda_h \in \Lambda: \lambda_h|_{\partial K} \in R_r( \partial K ) \, \forall K \in \mathcal{T}_h\right\}\,,
  \end{array}
  \f]
  where \f$ P_r(K) \f$ is the space of polynomial of degree \f$ r \f$ in the element \f$ K \f$, \f$ RT_r(K) \f$ is the space of
  polynomial of Raviart-Thomas of degrees \f$ r \f$ in the element \f$ K \f$ and \f$ R_r(\partial K) \f$ is the space of polynomial of
  degree \f$ r \f$ definite on each of the boundary of the element \f$ K \f$ and discontinuous from one edge to the other.
  <BR>
  The finite dimensional problem is: find \f$ (\sigma_h,\,, p_h, \, \lambda_h) \in Z_h \times V_h \times \Lambda_h \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma_h(t), \tau_h, p_h(t)) + b(p_h(t), \tau_h) + c(\lambda_h(t), \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  \displaystyle -\frac{\partial }{\partial t} m(p_h(t),v_h) + b(v_h, \sigma_h) = -F(v_h)\,,                                  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
  -c(\mu_h, \sigma_h(t)) - h(\lambda_h(t(, \mu_h) = - G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
  \end{array}
  \right.
  \f]
  To obatin a fully discrete problem we discretize the time derivative via an implicit Euler scheme, with \f$ \Delta t \f$ as time step. Choosing
  \f[
  p_h^0 = \prod_{V_h} p_0
  \f]
  we obtain the following system for each \f$ n = 0, \ldots, N \f$, with \f$ N = \frac{T}{\Delta t} \f$
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma_h^{n+1}, \tau_h, p_h^{n+1}) + b(p_h^{n+1}, \tau_h) + c(\lambda_h^{n+1}, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  \displaystyle -\frac{1}{\Delta t} m(p_h^{n+1},v_h) + b(v_h, \sigma_h^{n+1}) = -F(v_h) - \frac{1}{\Delta t} m(p_h^n, v_h) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
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
  a(\sigma_h^{n+1,k}, \tau_h, p_h^{n+1,k-1}) + b(p_h^{n+1,k}, \tau_h) + c(\lambda_h^{n+1,k}, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
\displaystyle -\frac{1}{\Delta t} m(p_h^{n+1,k},v_h) + b(v_h, \sigma_h^{n+1,k}) = -F(v_h) - \frac{1}{\Delta t} m(p_h^{n}, v_h) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
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
  \left[ H \right]_{ij} =    \displaystyle \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
  \left[ M \right]_{ij} =    \displaystyle \int_{K} \phi_j \phi_i \,, &
  \left[ F \right]_{j}  =    \displaystyle \int_K f \phi_j\,,  \vspace{0.2cm} \\
  \left[ G \right]_{j}  =    \displaystyle \int_{\partial K \cap \Gamma_R } g \xi_j\,,
  \end{array}
  \f]
  where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors. <BR>
  The local matrix formulation of the finite dimensional problem is
  \f[
  \left\{
  \begin{array}{l}
  A(p^{n+1,k-1}) \sigma_K^{n+1,k} + B p_K^{n+1,k} + C \lambda_K^{n+1,k} = 0\,, \vspace{0.2cm} \\
  \displaystyle -\frac{1}{\Delta t} M p_K^{n+1,k} +B^T \sigma_K^{n+1,k}  = -F - \frac{1}{\Delta t} M p_K^n  \,,                    \vspace{0.2cm}\\
  -C^T \sigma_K^{n+1,k} - H \lambda_K^{n+1,k} = - G\,.
  \end{array}
  \right.
  \f]
  Or alternatively
  \f[
  \begin{array}{l l l}
  \left[
  \begin{array}{c c c}
  A(p_K^{n+1,k-1})    & B &  C \vspace{0.2cm} \\
  -B^T & \displaystyle +\frac{1}{\Delta t} M &  0 \vspace{0.2cm} \\
  -C^T & 0 & -H
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
  0 \vspace{0.2cm}\\
  F + \frac{1}{\Delta t} M p_K^n \vspace{0.2cm}\\
  -G
  \end{array}
  \right]\,.
  \end{array}
  \f]
  Introducing the local hybrid matrix and local hybrid right hand side
  \f[
  \begin{array}{l}
  L_K =  C^T A^{-1}( p_k^{n+1,k-1}) C - C^T A^{-1}( p_k^{n+1,k-1}) B (  \frac{1}{\Delta t} M + B^T A^{-1}( p_k^{n+1,k-1}) B )^{-1} B^T A^{-1}( p_k^{n+1,k-1}) C - H \,, \vspace{0.2cm} \\
  r_K = - G - C^T A^{-1}( p_k^{n+1,k-1}) B ( B^T A^{-1}( p_k^{n+1,k-1}) B )^{-1} F\,,
  \end{array}
  \f]
  Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
  \f[
  L \lambda^{n+1,k} = r \,.
  \f]
  We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
  \f[
  \begin{array}{l}
  p_K^{n+1,k} = ( \frac{1}{\Delta t} M + B^T A^{-1}( p_k^{n+1,k-1}) B )^{-1} ( F + \frac{1}{\Delta t} M p_K^n - B^T A^{-1}( p_k^{n+1,k-1}) C \lambda_K^{n+1,k} )\,, \vspace{0.2cm} \\
  \sigma_K^{n+1,k} = -A^{-1} ( p_k^{n+1,k-1})( B p_K^{n+1,k} + C \lambda_K^{n+1,k} )  \,.
  \end{array}
  \f]
  @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
  conditions are imposed via BCHandler class.
  @todo Insert any scientific publications that use this solver.
*/
template< typename Mesh, typename SolverType = LifeV::SolverTrilinos >
class DarcySolverTransientNonLinear
        :
        public DarcySolverNonLinear<Mesh, SolverType>,
        public DarcySolverTransient<Mesh, SolverType>
{

public:

    //! @name Public Types
    //@{

    typedef SolverType                    solver_Type;
    typedef DarcySolver<Mesh, solver_Type> DarcySolverPolicies;

    typedef typename DarcySolverPolicies::permeability_Type  permeability_Type;

    typedef typename DarcySolverPolicies::data_Type          data_Type;

    typedef typename DarcySolverPolicies::bchandler_raw_Type bchandler_raw_Type;
    typedef typename DarcySolverPolicies::bchandler_Type     bchandler_Type;

    typedef typename DarcySolverPolicies::matrix_Type        matrix_Type;
    typedef typename DarcySolverPolicies::matrixPtr_Type     matrixPtr_Type;

    typedef typename DarcySolverPolicies::vector_Type        vector_Type;
    typedef typename DarcySolverPolicies::vectorPtr_Type     vectorPtr_Type;

    typedef typename DarcySolverPolicies::comm_Type          comm_Type;
    typedef typename DarcySolverPolicies::commPtr_Type       commPtr_Type;

    //@}

    //! @name Constructors & destructor
    //@{

    /*!
      Full constructor for the class.
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer for the Epetra communicator.
    */
    DarcySolverTransientNonLinear ( const data_Type&          dataFile,
                                    FESpace<Mesh, EpetraMap>& primal_FESpace,
                                    FESpace<Mesh, EpetraMap>& dual_FESpace,
                                    FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                                    FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                                    bchandler_raw_Type&       bcHandler,
                                    commPtr_Type&             comm );

    /*!
      Constructor for the class without the definition of the boundary handler.
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual finite element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer for the Epetra communicator.
    */
    DarcySolverTransientNonLinear ( const data_Type&          dataFile,
                                    FESpace<Mesh, EpetraMap>& primal_FESpace,
                                    FESpace<Mesh, EpetraMap>& dual_FESpace,
                                    FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                                    FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                                    commPtr_Type&             comm );

    //! Virtual destructor.
    virtual ~DarcySolverTransientNonLinear ();

    //@}
 
    //! @name Methods
    //@{

    //!  Set up the linear solver, the preconditioner for the linear system and the exporter to save the solution.
    virtual void setup ();

    //! Update the primalOld solution needed, as starting solution, in the fixed point scheme.
    inline void updatePrimalOldSolution()
    {

        // Reset the primal old vector.
        this->M_primalOld.reset( new vector_Type( this->M_primal_FESpace.map() ) );

        // Update the primal solution.
        *this->M_primalOld = *(this->M_primal);

    }

    //@}

    //! @name Set methods
    //@{

    //! Set the inverse of diffusion tensor
    /*!
      Set the inverse of diffusion tensor, the default (inverse of) permeability is the identity matrix.
      @param invPerm Inverse of the permeability tensor for the problem.
    */
    inline void setInversePermeability ( const permeability_Type& invPerm )
    {

        // Call the set inverse permeability of the non-linear Darcy solver
        DarcySolverNonLinear<Mesh, solver_Type>::setInversePermeability( invPerm );

    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Update all problem variables
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class.
    */
    virtual void updateVariables ();

    //! Compute elementary matrices
    /*!
      Locally update the current finite element for the primal
      and dual finite element space, then compute the Hdiv mass
      matrix.
      @param iElem Id of the current geometrical element.
    */
    virtual void localElementComputation ( const UInt & iElem );

    //@}

}; // class DarcySolverTransientNonLinear

//
// IMPLEMENTATION
//

// ===================================================
// Constructors & Destructor
// ===================================================

// Complete constructor.
template<typename Mesh, typename SolverType>
DarcySolverTransientNonLinear<Mesh, SolverType>::
DarcySolverTransientNonLinear ( const data_Type&           dataFile,
                                FESpace<Mesh, EpetraMap>&  primal_FESpace,
                                FESpace<Mesh, EpetraMap>&  dual_FESpace,
                                FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                                FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                                bchandler_raw_Type&        bcHandler,
                                commPtr_Type&              comm ):
        // Standard Darcy solver constructor.
        DarcySolver<Mesh, solver_Type>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, bcHandler, comm),
        // Non-linear Darcy solver constructor.
        DarcySolverNonLinear<Mesh, solver_Type>::DarcySolverNonLinear( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, bcHandler, comm),
        // Transient Darcy solver contructor.
        DarcySolverTransient<Mesh, solver_Type>::DarcySolverTransient( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, bcHandler, comm)
{

} // Constructor


// Constructor without boundary condition handler.
template<typename Mesh, typename SolverType>
DarcySolverTransientNonLinear<Mesh, SolverType>::
DarcySolverTransientNonLinear ( const data_Type&           dataFile,
                                FESpace<Mesh, EpetraMap>&  primal_FESpace,
                                FESpace<Mesh, EpetraMap>&  dual_FESpace,
                                FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                                FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                                commPtr_Type&              comm ):
        // Standard Darcy solver constructor.
        DarcySolver<Mesh, solver_Type>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, comm),
        // Non-linear Darcy solver constructor.
        DarcySolverNonLinear<Mesh, solver_Type>::DarcySolverNonLinear( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, comm),
        // Transient Darcy solver contructor.
        DarcySolverTransient<Mesh, solver_Type>::DarcySolverTransient( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, comm)
{

} // Constructor

// Virtual destructor.
template<typename Mesh, typename SolverType>
DarcySolverTransientNonLinear<Mesh, SolverType>::
~DarcySolverTransientNonLinear ( void )
{

} // Destructor

// ===================================================
// Public methods
// ===================================================

// Set up the linear solver and the preconditioner.
template<typename Mesh, typename SolverType>
void
DarcySolverTransientNonLinear<Mesh, SolverType>::
setup ()
{

    GetPot dataFile( *(this->M_data.dataFile()) );

    // Call the DarcySolverTransient setup method for setting up the linear solver.
    DarcySolverTransient<Mesh, solver_Type>::setup();

    // Set the maximum number of iteration for the fixed point iteration scheme.
    this->M_maxIterationFixedPoint = static_cast<UInt>( dataFile( ( this->M_data.section() + "/non-linear/fixed_point_iteration" ).data(), 10 ) );

    // Set the tollerance for the fixed point iteration scheme.
    this->M_tollFixedPoint = dataFile( ( this->M_data.section() + "/non-linear/fixed_point_toll" ).data(), 1.e-8 );

} // setup

// ===================================================
// Protected methods
// ===================================================

// Update all the variables of the problem.
template<typename Mesh, typename SolverType>
void
DarcySolverTransientNonLinear<Mesh, SolverType>::
updateVariables ()
{

    // Reset the primal vector at the previous iteration step.
    this->M_primalPreviousIteration.reset( new vector_Type( this->M_primal_FESpace.map() ) );

    // Update the primal vector at the previous iteration step.
    *(this->M_primalPreviousIteration) = *(this->M_primal);

    // Call the method of the DarcySolver to update all the variables defined in it.
    DarcySolver<Mesh, solver_Type>::updateVariables();

} // updateVariables

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template<typename Mesh, typename SolverType>
void
DarcySolverTransientNonLinear<Mesh, SolverType>::
localElementComputation ( const UInt & iElem )
{

    DarcySolverTransient<Mesh, solver_Type>::localElementComputation( iElem );

} // localElementComputation

} // namespace LifeV


#endif //_DARCYSOLVERTRANSIENTNONLINEAR_H_

// -*- mode: c++ -*-
