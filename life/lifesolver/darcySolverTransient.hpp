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

#include <life/lifesolver/darcySolver.hpp>

namespace
{

using namespace LifeV;

Real _One_ ( const Real&, const Real&, const Real&, const Real&, const ID& )
{
    return static_cast<Real>( 1. );
}

}

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
  \Lambda^{-1}(t) \sigma + \nabla p = 0                           & \mathrm{in} & \Omega   \times [0, T]\,, \vspace{0.2cm} \\
  \displaystyle \frac{\partial p}{\partial t} + \nabla \cdot \sigma - f(t) = 0  & \mathrm{in} & \Omega   \times [0, T]\,, \vspace{0.2cm} \\
  p = g_D(t)                                                      & \mathrm{on} & \Gamma_D \times [0, T]\,, \vspace{0.2cm} \\
  \sigma \cdot n + h(t) p = g_R(t)                                & \mathrm{on} & \Gamma_R \times [0, T]\,. \vspace{0.2cm} \\
  p(0) = p_0                                                      & \mathrm{in} & \Omega
  \end{array}
  \right.
  \f]
  Where \f$ p_0 \f$ is the primal initial value, \f$ \Lambda \f$ is the permeability tensor,
  \f$ f \f$ is the source term, \f$ \Gamma_D \f$ is the subset of the boundary of \f$ \Omega \f$ with Dirichlet
  boundary conditions with datum \f$ g_D \f$ and \f$ \Gamma_R \f$ is the part of the boundary of \f$ \Omega \f$
  with Robin, or Neumann, boundary conditions with data \f$ h \f$ and \f$ g_R \f$. We suppose that
  \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and \f$ \Gamma_D \cap \Gamma_R = \emptyset \f$.
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
  h(\lambda(t), \mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h(t) \mu \lambda(t) \,,\vspace{0.2cm}\\
  m(p(t), v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K p(t) v \,,\vspace{0.2cm}\\
  F(v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f(t) v\,,&
  G(\mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g(t) \mu\,,
  \end{array}
  \f]
  we obtain the Darcy problem in the weak form: find \f$ (\sigma(t), \, p(t), \, \lambda(t)) \in Z \times V \times \Lambda \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma(t), \tau) + b(p(t), \tau) + c(\lambda(t), \tau) = 0\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
  \displaystyle -\frac{\partial }{\partial t} m(p(t), v) + b(v, \sigma(t)) = -F(v)\,, & \forall v \in V \,,\vspace{0.2cm}\\
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
  <BR>
  The finite dimensional problem is: find \f$ (\sigma_h(t),\, p_h(t), \, \lambda_h(t)) \in Z_h \times V_h \times \Lambda_h \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma_h(t), \tau_h) + b(p_h(t), \tau_h) + c(\lambda_h(t), \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  \displaystyle -\frac{\partial }{\partial t} m(p_h(t),v_h) + b(v_h, \sigma_h(t)) = -F(v_h)\,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
  c(\mu_h, \sigma_h(t)) + h(\lambda_h(t), \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
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
  a(\sigma_h^{n+1}, \tau_h) + b(p_h^{n+1}, \tau_h) + c(\lambda_h^{n+1}, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  \displaystyle -\frac{1}{\Delta t} m(p_h^{n+1},v_h) + b(v_h, \sigma_h^{n+1}) = -F(v_h) - \frac{1}{\Delta t} m(p_h^n, v_h) \,,  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
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
  \left[ H \right]_{ij} = \displaystyle   \int_{\partial K \cap \Gamma_R} h( (n+1) \Delta t) \xi_i \xi_j\,, &
  \left[ M \right]_{ij} = \displaystyle   \int_{K} \phi_j \phi_i \,, &
  \left[ F \right]_{j}  = \displaystyle   \int_K f( (n+1) \Delta t) \phi_j\,, \vspace{0.2cm} \\
  \left[ G \right]_{j}  = \displaystyle   \int_{\partial K \cap \Gamma_R } g( (n+1) \Delta t) \xi_j\,,
  \end{array}
  \f]
  where we avoid to write the dependence on the triangle \f$ K \f$ and on the current time step \f$ n+1 \f$ in all the matrices and vectors. <BR>
  The local matrix formulation of the finite dimensional problem is
  \f[
  \left\{
  \begin{array}{l}
  A \sigma_K^{n+1} + B p_K^{n+1} + C \lambda_K^{n+1} = 0\,, \vspace{0.2cm} \\
  \displaystyle -\frac{1}{\Delta t} M p_K^{n+1} + B^T \sigma_K^{n+1} = -F - \frac{1}{\Delta t} M p_K^n \,,  \vspace{0.2cm}\\
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
  B^T & \displaystyle -\frac{1}{\Delta t} M &  0 \vspace{0.2cm} \\
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
  0 \vspace{0.2cm}\\
  \displaystyle - F - \frac{1}{\Delta t} M p_K^n \vspace{0.2cm}\\
  G
  \end{array}
  \right]\,.
  \end{array}
  \f]
  Introducing the local hybrid matrix and local hybrid right hand side
  \f[
  \begin{array}{l}
  \displaystyle L_K = -C^T A^{-1} C + C^T A^{-1} B \left( \frac{1}{\Delta t} M + B^T A^{-1} B \right)^{-1} B^T A^{-1} C + H \,, \vspace{0.2cm} \\
  \displaystyle r_K = G + C^T A^{-1} B \left( \frac{1}{\Delta t} M + B^T A^{-1} B \right)^{-1} \left( F + \frac{1}{\Delta t} M p_K^n \right) \,,
  \end{array}
  \f]
  Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
  \f[
  L \lambda^{n+1} = r \,.
  \f]
  We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
  \f[
  \begin{array}{l}
  \displaystyle p_K^{n+1} = \left( \frac{1}{\Delta t} M + B^T A^{-1} B \right)^{-1} \left( F + \frac{1}{\Delta t} M p_K^n - B^T A^{-1} C \lambda_K^{n+1} \right)\,, \vspace{0.2cm} \\
  \sigma_K^{n+1} = -A^{-1} \left( B p_K^{n+1} + C \lambda_K^{n+1} \right) \,.
  \end{array}
  \f]
  @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
  @note The initial time is not fix at zero.
  conditions are imposed via BCHandler class.
  @todo Insert any scientific publications that use this solver.
  @todo Post process for the dual variable.
  @todo Use a better mass assembler
*/
template< typename Mesh, typename SolverType = LifeV::SolverTrilinos >
class DarcySolverTransient :
        virtual public DarcySolver<Mesh, SolverType>
{

public:

    //! @name Pulbic Types
    //@{

    typedef DarcySolver<Mesh, SolverType> DarcySolverPolicies;

    typedef typename DarcySolverPolicies::Function           Function;

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
      @param comm Shared pointer of the Epetra communicator.
    */
    DarcySolverTransient ( const data_Type&          dataFile,
                           FESpace<Mesh, EpetraMap>& primal_FESpace,
                           FESpace<Mesh, EpetraMap>& dual_FESpace,
                           FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                           FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                           bchandler_Type&           bcHandler,
                           const commPtr_Type&             comm );
    /*!
      Constructor for the class without boundary condition handler.
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param comm Shared pointer of the Epetra communicator.
    */
    DarcySolverTransient ( const data_Type&          dataFile,
                           FESpace<Mesh, EpetraMap>& primal_FESpace,
                           FESpace<Mesh, EpetraMap>& dual_FESpace,
                           FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                           FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                           const commPtr_Type&             comm );

    //! Virtual destructor.
    virtual ~DarcySolverTransient ();

    //@}

    // Methods
    //! @name methods
    //@{

    //! Set up the linear solver and the preconditioner for the linear system.
    virtual void setup ();

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
    void setInitialPrimal ( const Function& primalInitial );

    //! Set mass matrix
    /*!
      Set the mass term, the default source term is the function one.
      By defaul it does not depend on time.
      @param mass Mass term for the problem.
    */

    void setMass ( const Function& mass )
    {
        M_mass = mass;
    }

    void __attribute__ ((__deprecated__)) setMassTerm ( const Function& mass )
    {
	return setMass(mass);
    }
    //@}

protected:

    // Methods
    //! @name methods
    //@{

    //! Compute element matrices
    /*!
      Locally update the current finite element for the primal
      and dual finite element space, then compute the Hdiv mass
      matrix.
      @param iElem Id of the current geometrical element.
    */
    virtual void localElementComputation ( const UInt & iElem );

    //! Perform static condensation
    /*!
      Create the local hybrid matrix and local hybrid right hand side.
    */
    virtual void staticCondensation ();

    //! Compute locally, as a post process, the primal and dual variable.
    virtual void localComputePrimalAndDual ();

    //! Update all variable
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
    */
    virtual void updateVariables ();

    //@}

    //! @name Protected data
    //@{

    //! Primal solution at previous time.
    vectorPtr_Type  M_primalOld;

    //@}


private:

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //! Initial time primal variable.
    Function    M_primalInitial;

    //! Mass function, it does not depend on time.
    Function    M_mass;

    //@}

    // Algebraic stuff.
    //! @name Algebraic stuff
    //@{

    //! Boolean that indicates if the preconditioner is re-used or not.
    bool            M_reusePrec;

    //! Boolean that indicates if the matrix is updated for the current iteration.
    bool            M_updated;

    //! Interger storing the max number of solver iteration with preconditioner recomputing.
    UInt             M_maxIterSolver;

    //! Boolean that indicates if the matrix is recomputed for the current iteration.
    bool            M_recomputeMatrix;

    //@}

    // Elementary matrices and vectors used for the static condensation.
    //! @name Elementary matrices and vectors used for the static condensation.
    //@{

    /*! Temporary array of size the square of the number of degrees of freedom of the
        primal variable. */
    ElemMat         M_elmatMassPrimal;

    //@}

}; // class DarcySolverTransient

//
// IMPLEMENTATION
//

// Complete constructor.
template<typename Mesh, typename SolverType>
DarcySolverTransient<Mesh, SolverType>::
DarcySolverTransient ( const data_Type&           dataFile,
                       FESpace<Mesh, EpetraMap>&  primal_FESpace,
                       FESpace<Mesh, EpetraMap>&  dual_FESpace,
                       FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                       FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                       bchandler_Type&            bcHandler,
                       const commPtr_Type&              comm ):
        // Standard Darcy solver constructor.
        DarcySolver<Mesh, SolverType>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, bcHandler, comm),
        // Data of the problem
        M_primalInitial          ( NULL ),
        M_mass                   ( _One_ ),
        // Linear solver.
        M_primalOld              ( new vector_Type ( this->M_primal_FESpace.map() ) ),
        M_reusePrec              ( false ),
        M_updated                ( false ),
        M_maxIterSolver          ( static_cast<UInt>(0) ),
        M_recomputeMatrix        ( false ),
        // Local matrices and vectors
        M_elmatMassPrimal        ( this->M_primal_FESpace.refFE().nbDof(), 1, 1 )
{

} // Constructor

// Constructor without boundary condition handler.
template<typename Mesh, typename SolverType>
DarcySolverTransient<Mesh, SolverType>::
DarcySolverTransient ( const data_Type&           dataFile,
                       FESpace<Mesh, EpetraMap>&  primal_FESpace,
                       FESpace<Mesh, EpetraMap>&  dual_FESpace,
                       FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                       FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                       const commPtr_Type&              comm ):
        // Standard Darcy solver constructor.
        DarcySolver<Mesh, SolverType>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, comm),
        // Data of the problem
        M_primalInitial          ( NULL ),
        M_mass                   ( _One_ ),
        // Linear solver.
        M_primalOld              ( new vector_Type ( this->M_primal_FESpace.map() ) ),
        M_reusePrec              ( false ),
        M_updated                ( false ),
        M_maxIterSolver          ( static_cast<UInt>(0) ),
        M_recomputeMatrix        ( false ),
        // Local matrices and vectors
        M_elmatMassPrimal        ( this->M_primal_FESpace.refFE().nbDof(), 1, 1 )
{

} // Constructor

// Virtual destructor.
template<typename Mesh, typename SolverType>
DarcySolverTransient<Mesh, SolverType>::
~DarcySolverTransient ( void )
{

} // Destructor

// ===========================================================================================
// Public methods
// ==========================================================================================

// Set up the linear solver and the preconditioner.
template<typename Mesh, typename SolverType>
void
DarcySolverTransient<Mesh, SolverType>::
setup ()
{

    GetPot dataFile( *(this->M_data.dataFile()) );

    // Call the DarcySolver setup method for setting up the linear solver.
    DarcySolver<Mesh, SolverType>::setup();

    // Set if the preconditioner is re-used.
    M_reusePrec = dataFile( ( this->M_data.section() + "/solver/reuse" ).data(), false );

    // Set if the preconditioner is reused or not.
    this->M_linearSolver.setReusePreconditioner( M_reusePrec );

    // Set the max number of iteration to mantein the same preconditioner.
    M_maxIterSolver = static_cast<UInt>( dataFile( ( this->M_data.section() + "/solver/max_iter_reuse" ).data(),
                                                   static_cast<Int>(0) ) );

} // setup

// Set the inital value
template<typename Mesh, typename SolverType>
inline
void
DarcySolverTransient<Mesh, SolverType>::
setInitialPrimal ( const Function& primalInitial )
{
    // Set the initial value function.
    M_primalInitial = primalInitial;

    // Interpolate the primal initial value.
    this->M_primal_FESpace.interpolate( M_primalInitial,
                                        *(this->M_primal),
                                        this->M_data.dataTime()->getInitialTime() );

} // setInitialPrimal

// ===========================================================================================
// Protected methods
// ==========================================================================================

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template <typename Mesh, typename SolverType>
void
DarcySolverTransient<Mesh, SolverType>::
localElementComputation ( const UInt & iElem )
{

    /* Call the DarcySolver localElementComputation for update the local finite elements
       and then compute the massHdiv matrix. */
    DarcySolver<Mesh, SolverType>::localElementComputation( iElem );

    // Clear the mass matrix for the primal variable.
    M_elmatMassPrimal.zero();

    // Get the coordinate of the barycenter of the current element of ID iElem.
    Real xg(0), yg(0), zg(0);
    this->M_primal_FESpace.fe().barycenter( xg, yg, zg );

    // Compute the mass matrix for the primal variable.
    mass( M_mass( 0., xg, yg, zg, 0),
          M_elmatMassPrimal,
          this->M_primal_FESpace.fe(), 0, 0);

    // Store in the mass matrix for the primal variable also the time step.
    M_elmatMassPrimal *= static_cast<Real>(1.) / this->M_data.dataTime()->getTimeStep();

} // localElementComputation

// Perform the static condensation for the local hybrid matrix.
template <typename Mesh, typename SolverType>
void
DarcySolverTransient<Mesh, SolverType>::
staticCondensation ()
{

    // Flags for the BLAS and LAPACK routine.

    int INFO[1]         = {0};
    // Number of columns of the right hand side := 1.
    int NBRHS[1]        = {1};
    // Primal variable degrees of freedom.
    int NBP[1]          = { this->M_primal_FESpace.refFE().nbDof() };
    // Dual variable degrees of freedom.
    int NBU[1]          = { this->M_dual_FESpace.refFE().nbDof() };
    // Hybrid variable degree of freedom.
    int NBL[1]          = { this->M_hybrid_FESpace.refFE().nbDof() };
    // Increment := 1.
    int INC[1]         = {1};

    double ONE[1]      = {1.0};
    double MINUSONE[1] = {-1.0};
    double ZERO[1]     = {0.0};

    // Parameter that indicate the Lower storage of matrices.
    char UPLO[1]     = {'L'};
    // Parameter that indicate the Transpose of matrices.
    char   TRANS[1]    = {'T'};
    char NOTRANS[1]    = {'N'};
    // Parameter that indicates whether the matrix has diagonal unit ('N' means no)
    char NODIAG[1] = {'N'};

    // Create and assign the local matrices A, B and C.
    ElemMat::matrix_type A = this->M_elmatMix.block( 0, 0 );
    ElemMat::matrix_type B = this->M_elmatMix.block( 0, 1 );
    ElemMat::matrix_type C = this->M_elmatMix.block( 0, 2 );

    //.................................
    //        MATRIX OPERATIONS
    //.................................

    /* Put in A the matrix L and L^T, where L and L^T is the Cholesky factorization of A.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    dpotrf_ (UPLO, NBU, A, NBU, INFO);
    ASSERT_PRE( !INFO[0], "Lapack factorization of A is not achieved." );

    /* Put in B the matrix L^{-1} * B, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBU, NBP, A, NBU, B, NBU, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation B = L^{-1} B  is not achieved." );

    /* Put in C the matrix L^{-1} * C, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBU, NBL, A, NBU, C, NBL, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation C = L^{-1} C  is not achieved." );

    /* Put in M_BtB the matrix  B^T * L^{-T} * L^{-1} * B = B^T * A^{-1} * B
       M_BtB stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f */
    dsyrk_ (UPLO, TRANS, NBP, NBU, ONE, B, NBU, ZERO, this->M_BtB, NBP);

    /* Put in M_BtB the matrix
       M_BtB + M_elmatMassPrimal / \Delta t = B^T * A^{-1} * B + M_elmatMassPrimal / \Delta t
       M_BtB stored only on lower part. */
    this->M_BtB += M_elmatMassPrimal.mat();

    /* Put in M_CtC the matrix C^T * L^{-T} * L^{-1} * C = C^T * A^{-1} * C
       M_CtC stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f */
    dsyrk_ (UPLO, TRANS, NBL, NBU, ONE, C, NBU, ZERO, this->M_CtC, NBL);

    /* Put in M_BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       M_BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBP, NBL, NBU, ONE, B, NBU, C, NBU, ZERO, this->M_BtC, NBP);

    /* Put in M_BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B + M_elmatMassPrimal / \Delta t.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f  */
    dpotrf_ (UPLO, NBP, this->M_BtB, NBP, INFO);
    ASSERT_PRE( !INFO[0],"Lapack factorization of BtB is not achieved." );

    /* Put in M_BtC the matrix LB^{-1} * M_BtC = LB^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBL, this->M_BtB, NBP, this->M_BtC, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    /* Put in M_CtC the matrix -M_CtC + M_BtC^T * M_BtC
       Result stored only on lower part, the matrix M_CtC stores
       M_CtC = -C^T * A^{-1} * C + C^T * A^{-t} * B * ( M_elmatMassPrimal / \Delta t + B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f  */
    dsyrk_ (UPLO, TRANS, NBL, NBP, ONE, this->M_BtC, NBP, MINUSONE, this->M_CtC, NBL);

    //...................................
    //      END OF MATRIX OPERATIONS
    //...................................

    /* Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       B^T stores LB and LB^T where LB and LB^T is the factorization of B^T * A^{-1} * B
       M_BtC stores LB^{-1} * B^T * A^{-1} * C
       M_CtC stores C^T * A^{-1} * C - C^T * A^{-t} * B * (B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C */

    //..........................
    //     VECTOR OPERATIONS
    //..........................

    // Clear some vectors.
    this->M_elvecSource.zero();
    this->M_elvecHyb.zero();

    // Compute the right hand side.
    source( this->M_source,
            this->M_elvecSource,
            this->M_primal_FESpace.fe(),
            this->M_data.dataTime()->getTime(), 0 );

    // Take the pressure at the previews time step in the local element.
    ElemVec elvecPrimalOldLocal( this->M_primal_FESpace.refFE().nbDof(), 1 );

    // Extract the old primal variable for the current finite element and put it into elvecPrimalOldLocal.
    extract_vec( *M_primalOld,
                 elvecPrimalOldLocal,
                 this->M_primal_FESpace.refFE(),
                 this->M_primal_FESpace.dof(),
                 this->M_primal_FESpace.fe().currentLocalId(), 0 );

    /* Put in M_elvecSource the vector M_elmatMassPrimal / \Delta t * elvecPrimalOldLocal + M_elvecSource
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBP, NBP, ONE, M_elmatMassPrimal.mat(), NBP, elvecPrimalOldLocal, INC, ZERO, this->M_elvecSource, INC);

    /* Put in M_elvecSource the vector LB^{-1} * M_elvecSource = LB^{-1} *( M_elmatMassPrimal / \Delta t + F)
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBRHS, this->M_BtB, NBP, this->M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

    /* Put in M_elvecHyb the vector
       M_BtC^T * M_elvecSource =
       C^T * A^{-1} * B^T * ( M_elmatMassPrimal / \Delta t + B^T * A^{-1} * B)^{-1} * ( F + M_elmatMassPrimal / \Delta t * M_primalOld )
       M_elvecHyb is fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBL, NBRHS, NBP, ONE, this->M_BtC, NBP, this->M_elvecSource, NBP, ZERO, this->M_elvecHyb, NBL);

    //........................
    // END OF VECTOR OPERATIONS.
    //........................

    /* Previously the matrix M_CtC is stored only in the lower part, but at the moment there is not
       a function assembleMatrix that store a lower triangular sparse matrix.
       Remind to correct these line in the future. */
    this->symmetrizeMatrix( UPLO, NBL, this->M_CtC );

    // Update the hybrid element matrix.
    this->M_elmatHyb.block(0,0) = this->M_CtC;

} // staticCondensation

// Locally compute the primal and dual variable.
template<typename Mesh, typename SolverType>
void
DarcySolverTransient<Mesh, SolverType>::
localComputePrimalAndDual ()
{

    // Flags for the BLAS and LAPACK routine.

    int INFO[1]         = {0};
    // Number of columns of the right hand side := 1.
    int NBRHS[1]        = {1};
    // Primal variable degrees of freedom.
    int NBP[1]          = { this->M_primal_FESpace.refFE().nbDof() };
    // Dual variable degrees of freedom.
    int NBU[1]          = { this->M_dual_FESpace.refFE().nbDof() };
    // Hybrid variable degree of freedom.
    int NBL[1]          = { this->M_hybrid_FESpace.refFE().nbDof() };
    // Increment := 1.
    int INC[1]         = {1};

    double ONE[1]      = {1.0};
    double MINUSONE[1] = {-1.0};
    double ZERO[1]     = {0.0};

    // Parameter that indicate the Lower storage of matrices.
    char UPLO[1]     = {'L'};
    // Parameter that indicate the Transpose of matrices.
    char   TRANS[1]    = {'T'};
    char NOTRANS[1]    = {'N'};
    // Parameter that indicates whether the matrix has diagonal unit ('N' means no)
    char NODIAG[1] = {'N'};

    // No need for CtC in this part, and last dsyrk, the only differences.
    ElemMat::matrix_type A = this->M_elmatMix.block( 0, 0 );
    ElemMat::matrix_type B = this->M_elmatMix.block( 0, 1 );
    ElemMat::matrix_type C = this->M_elmatMix.block( 0, 2 );

    //........................
    //    MATRIX OPERATIONS
    //........................

    /* Put in A the matrix L and L^T, where L and L^T is the Cholesky factorization of A.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    dpotrf_ (UPLO, NBU, A, NBU, INFO);
    ASSERT_PRE( !INFO[0], "Lapack factorization of A is not achieved." );

    /* Put in B the matrix L^{-1} * B, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBU, NBP, A, NBU, B, NBU, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation B = L^{-1} B  is not achieved." );

    /* Put in C the matrix L^{-1} * C, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBU, NBL, A, NBU, C, NBU, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation C = L^{-1} C  is not achieved." );

    /* Put in M_BtB the matrix  B^T * L^{-T} * L^{-1} * B = B^T * A^{-1} * B
       M_BtB stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f */
    dsyrk_ (UPLO, TRANS, NBP, NBU, ONE, B, NBU, ZERO, this->M_BtB, NBP);

    /* Put in M_BtB the matrix
       M_BtB - M_elmatMassPrimal / \Delta t = B^T * A^{-1} * B - M_elmatMassPrimal / \Delta t
       M_BtB stored only on lower part. */
    this->M_BtB += M_elmatMassPrimal.mat();

    /* Put in M_BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       M_BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBP, NBL, NBU, ONE, B, NBU, C, NBU, ZERO, this->M_BtC, NBP);

    /* Put in M_BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B + M_elmatMassPrimal / \Delta t.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    dpotrf_ (UPLO, NBP, this->M_BtB, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack factorization of BtB is not achieved." );

    /* Put in M_BtC the matrix LB^{-1} * M_BtC = LB^{-1} * B^T * A^{-1} * C.
           For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBL, this->M_BtB, NBP, this->M_BtC, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    //..............................
    //   END OF MATRIX OPERATIONS
    //..............................

    /* Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       B^T stores LB and LB^T where LB and LB^T is the factorization of B^T * A^{-1} * B
       M_BtC stores LB^{-1} * B^T * A^{-1} * C */


    //......................
    //    VECTOR OPERATIONS (Computation of Pressure and Velocities)
    //......................

    //...................................
    //  1) Computation of the PRESSURE
    //...................................

    // Clear the source vector.
    this->M_elvecSource.zero();

    // The source term is computed with a test function in the primal variable space.
    source( this->M_source,
            this->M_elvecSource,
            this->M_primal_FESpace.fe(),
            this->M_data.dataTime()->getTime(), 0 );

    // Take the pressure at the previews time step in the local element.
    ElemVec elvecPrimalOldLocal( this->M_primal_FESpace.refFE().nbDof(), 1 );

    // Extract the old primal variable for the current finite element and put it into elvecPrimalOldLocal.
    extract_vec( *M_primalOld,
                 elvecPrimalOldLocal,
                 this->M_primal_FESpace.refFE(),
                 this->M_primal_FESpace.dof(),
                 this->M_primal_FESpace.fe().currentLocalId(), 0 );

    /* Put in M_elvecSource the vector M_elmatMassPrimal / \Delta t * elvecPrimalOldLocal + M_elvecSource
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBP, NBP, ONE, M_elmatMassPrimal.mat(), NBP, elvecPrimalOldLocal, INC, ONE, this->M_elvecSource, INC);

    /* Put in M_elvecSource the vector LB^{-1} * M_elvecSource = LB^{-1} * ( F + M_elmatMassPrimal / \Delta t * elvecPrimalOldLocal )
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBRHS, this->M_BtB, NBP, this->M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} M_elvecSource is not achieved." );

    /* Put in M_elvecSource the vector
       -M_BtC * M_elvecHyb + M_elvecSource = -LB^{-1} * B^T * A^{-1} * C * lambda_K + LB^{-1} * ( F + M_elmatMassPrimal / \Delta t)
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (NOTRANS, NOTRANS, NBP, NBRHS, NBL, MINUSONE, this->M_BtC, NBP, this->M_elvecHyb, NBL, ONE, this->M_elvecSource, NBP);

    /* Put in M_elvecSource the vector LB^{-T} * M_elvecSource, where
       M_elvecSource stores
       - (B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C * lambda_K + (B^T * A^{-1} * B)^{-1} * (F + M_elmatMassPrimal / \Delta t)
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, TRANS, NODIAG, NBP, NBRHS, this->M_BtB, NBP, this->M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-T} M_elvecSource is not achieved." );

    // Now rhs contains the primal variable for the current element, we must put it in the global vector.

    //.....................................
    //  2) Computation of the VELOCITIES
    //.....................................

    // Clear the element dual vector.
    this->M_elvecFlux.zero();

    /* Put in M_elvecFlux the vector B * M_elvecSource = L^{-1} * B * primal_K
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBU, NBP, ONE, B, NBU, this->M_elvecSource, INC, ZERO, this->M_elvecFlux, INC);

    /* Put in M_elvecFlux the vector
       - C * M_elvecHyb - M_elvecFlux = - L^{-1} * C * lambda_K - L^{-1} * B * primal_K
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBU, NBL, MINUSONE, C, NBL, this->M_elvecHyb, INC, MINUSONE, this->M_elvecFlux, INC);

    /* Put in flux the vector
       L^{-T} * M_elvecFlux = - A^{-1} * C^T * lambda_K - A^{-1} * B^T * primal_K
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, TRANS, NODIAG, NBU, NBRHS, A, NBU, this->M_elvecFlux, NBU, INFO);
    ASSERT_PRE(!INFO[0], "Lapack Computation M_elvecFlux = L^{-T} M_elvecFlux is not achieved.");

} // localComputePrimalAndDual

// Update all the variables of the problem.
template<typename Mesh, typename SolverType>
inline
void
DarcySolverTransient<Mesh, SolverType>::
updateVariables ()
{
    // Reset the primal old vector
    M_primalOld.reset( new vector_Type ( this->M_primal_FESpace.map() ) );

    // Update the primal solution.
    *M_primalOld = *( this->M_primal );

    // Call the method of the DarcySolver to update all the variables defined in it.
    DarcySolver<Mesh, SolverType>::updateVariables();

} // updateVariables

} // namespace LifeV

#endif // _DARCYSOLVERTRANSIENT_H_

// -*- mode: c++ -*-
