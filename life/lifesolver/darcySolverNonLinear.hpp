/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
      Date: 2010-05-09

 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
 Copyright (C) 2006-2010 Politecnico di Milano

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
  @file darcySolverNonLinear.hpp
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
  @date 05/2010

  @brief This file contains a non-linear permeability term Darcy equation solver class
*/

#ifndef _darcySolverNonLinear_H_
#define _darcySolverNonLinear_H_

#include <life/lifesolver/darcySolver.hpp>

// Local namespace to store the default function for start up the fixed point scheme and the default non-linear permeability tensor.
namespace
{

/*! @class Default Darcy start up function for the fixed point scheme. Needed in DarcySolverNonLinear class.
  This class implement the default start up function for the fixed point scheme in the the non-linear Darcy problem,
  it is the zero function.
*/
struct DarcyDefaultStartUpFunction
{
    const LifeV::Real operator()( const LifeV::Real&, const LifeV::Real&,
                                  const LifeV::Real&, const LifeV::Real&,
                                  const LifeV::UInt&) const
    {
    	return static_cast<LifeV::Real>( 0 );
    }
};

}

// LifeV namespace.
namespace LifeV
{
/*!
  @class DarcySolverNonLinear

  This class implements a non-linear Darcy solver with fixed point scheme to handle the non-linearity.
  <br>
  The classical non-linear Darcy formulation is a couple of differential equations of first order with
  the unknowns \f$ p \in C^1 (\Omega ) \f$, being the pressure or the primal unknown,
  and \f$ \sigma \in (C^1( \Omega ) )^n \f$, being the Darcy velocity or the flux or the dual unknown,
  such that
  \f[
  \left\{
  \begin{array}{l l l}
  \Lambda^{-1}(p) \sigma + \nabla p = 0 & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  \nabla \cdot \sigma - f = 0        & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  p = g_D                            & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
  \sigma \cdot n + h p = g_R         & \mathrm{on} & \Gamma_R\,.
  \end{array}
  \right.
  \f]
  Where \f$ \Lambda(p) \f$ is the permeability tensor that depends on \f$ p \f$, \f$ f \f$ is the source term,
  \f$ \Gamma_D \f$ is the subset  of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions with datum
  \f$ g_D \f$ and \f$ \Gamma_R \f$  is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary
  conditions with data \f$ h \f$ and \f$ g_R \f$. We suppose that \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and
  \f$ \Gamma_D \cap \Gamma_R = \emptyset \f$.
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
  H(div, K) = \left\{ \tau \in ( L^2(K))^n : \nabla \cdot \tau \in L^2(K)\right\}\,, \vspace{0.2cm}\\
  Z = \left\{ \tau \in L^2(\Omega) : \tau\vert_K \in H(div, K) \, \forall K \in  \mathcal{T}_h \right\}\,, \vspace{0.2cm}\\
  \Lambda = \left\{ \lambda \in \prod_{K \in \mathcal{T}_h} H^{1/2} (\partial K): \lambda_K = \lambda_{K'} \,\, \mathrm{on} \,\, e_{K-K'} \, \forall K \in \mathcal{T}_h,\, \lambda = g_D \,\, \mathrm{on} \,\, \Gamma_D \right\}\,.
  \end{array}
  \f]
  Introducing the following bilinear forms, functionals and operator
  \f[
  \begin{array}{l l}
  a(\sigma, \tau, p) = \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1}(p) \sigma \cdot \tau \,,  &
  b(p, \tau) = -\sum_{K \in \mathcal{T}_h} \int_K p \nabla \cdot \tau\,, \vspace{0.2cm}\\
  c(\lambda, \tau) = \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda \tau \cdot n\,,&
  h(\lambda, \mu) = \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda \,,\vspace{0.2cm}\\
  F(v) = \sum_{K \in \mathcal{T}_h} \int_K f v\,,&
  G(\mu) =\sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
  \end{array}
  \f]
  we obtain the Darcy problem in the weak form: find \f$ (\sigma, \, p, \, \lambda) \in Z \times V \times \Lambda \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma, \tau, p) + b(p, \tau) + c(\lambda, \tau) = 0\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
  b(v, \sigma) = - F(v)\,,                                  & \forall v \in V \,,\vspace{0.2cm}\\
  c(\mu, \sigma) + h(\lambda, \mu) = G(\mu) \,,           & \forall \mu \in \Lambda\,.
  \end{array}
  \right.
  \f]
  At discrete level we introduce the polynomial space, of degree \f$ r \f$, that approximate the finite dimensional
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
  a(\sigma_h, \tau_h, p_h) + b(p_h, \tau_h) + c(\lambda_h, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  b(v_h, \sigma_h) = - F(v_h)\,,                                  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
  c(\mu_h, \sigma_h) + h(\lambda_h, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
  \end{array}
  \right.
  \f]
  To solve the non-linearity we use a fixed point scheme based on the relative difference between two consecutive iterations
  of the primal variable. We start from the, user defined, function \f$ p^0 \f$ and solve the linearized problem for \f$ n \geq 1 \f$:
  find \f$ (\sigma_h^n,\, p_h^n, \, \lambda_h^n) \in Z_h \times V_h \times \Lambda_h \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma_h^n, \tau_h, p_h^{n-1}) + b(p_h^n, \tau_h) + c(\lambda_h^n, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  b(v_h, \sigma_h^n) = - F(v_h)\,,                                  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
  c(\mu_h, \sigma_h^n) + h(\lambda_h^n, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
  \end{array}
  \right.
  \f]
  At each iteration, to solve the problem, we use the static condensation procedure, i.e. the unknowns in the discrete
  weak system are not independent and \f$ p_K^n \f$, \f$\sigma_K^n \f$ may be written in function of
  \f$ \lambda_K^n \f$ alone. We introduce the following local matrices
  \f[
  \begin{array}{l l l}
  \left[ A \left( p_K^{n-1}\right) \right]_{ij} =   \int_K \Lambda^{-1}\left( p^{n-1} \right) \psi_j \cdot \psi_i \,, &
  \left[ B \right]_{ij} = - \int_K \phi_j \nabla \cdot \psi_i \,, &
  \left[ C \right]_{ij} =   \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
  \left[ H \right]_{ij} =   \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
  \left[ F \right]_{j}  =   \int_K f \phi_j\,, &
  \left[ G \right]_{j}  =   \int_{\partial K \cap \Gamma_R } g \xi_j\,,
  \end{array}
  \f]
  where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors.
  <br>
  The local matrix formulation of the finite dimensional problem, at each iteration, is
  \f[
  \left\{
  \begin{array}{l}
  A \left( p_K^{n-1} \right) \sigma_K^n + B p_K^n + C \lambda_K^n = 0\,, \vspace{0.2cm} \\
  B^T \sigma_K^n = -F \,,                    \vspace{0.2cm}\\
  C^T \sigma_K^n + H \lambda_K^n = G\,.
  \end{array}
  \right.
  \f]
  Or alternatively
  \f[
  \begin{array}{l l l}
  \left[
  \begin{array}{c c c}
  A \left( p_K^{n-1} \right) & B & C \vspace{0.2cm} \\
  B^T                        & 0 & 0 \vspace{0.2cm} \\
  C^T                        & 0 & H
  \end{array}
  \right] \, \cdot &
  \left[
  \begin{array}{c}
  \sigma_K^n  \vspace{0.2cm}\\
  p_K^n       \vspace{0.2cm}\\
  \lambda_K^n
  \end{array}
  \right] \, &
  =
  \left[
  \begin{array}{c}
  0 \vspace{0.2cm}\\
  -F \vspace{0.2cm}\\
  G
  \end{array}
  \right]\,.
  \end{array}
  \f]
  Introducing the local hybrid matrix and local hybrid right hand side
  \f[
  \begin{array}{l}
  L_K =  -C^T A^{-1}\left( p_k^{n-1}\right) C + C^T A^{-1} \left( p_K^{n-1} \right) B \left[ B^T A^{-1}\left( p_K^{n-1}\right) B \right]^{-1} B^T A^{-1}\left( p_K^{n-1} \right) C + H \,, \vspace{0.2cm} \\
  r_K = G + C^T A^{-1}\left( p_K^{n-1} \right) B \left[ B^T A^{-1}\left( p_K^{n-1} \right) B \right]^{-1} F\,,
  \end{array}
  \f]
  Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
  \f[
  L^n \lambda^n = r \,.
  \f]
  We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
  \f[
  \begin{array}{l}
  p_K^n = \left[ B^T A^{-1}\left( p_K^{n-1} \right) B \right]^{-1} \left[ F -  B^T A^{-1}\left( p_K^{n-1}\right) C \lambda_K^n \right]\,, \vspace{0.2cm} \\
  \sigma_K^n = -A^{-1}\left(p_K^{n-1} \right) \left( B p_K^n + C \lambda_K^n \right) \,.
  \end{array}
  \f]
  @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
  conditions are imposed via BCHandler class.
  @todo Insert any scientific publications that use this solver.
  @todo Post process for the dual variable.
  @todo Attention! We have the hypothesis that we use P0 elements for the primal unknown. Change this in a future!
  @todo Add criteria to ensure convergence of the fixed point method.
*/
template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class DarcySolverNonLinear
    :
        virtual public DarcySolver<Mesh, SolverType>
{

public:

    // Policies.
    //! @name Policies
    //@{

    typedef DarcySolver<Mesh, SolverType> DarcySolverPolicies;

    typedef typename DarcySolverPolicies::Function           Function;

    typedef typename DarcySolverPolicies::permeability_type  permeability_type;

    typedef typename DarcySolverPolicies::data_type          data_type;

    typedef typename DarcySolverPolicies::bchandler_raw_type bchandler_raw_type;
    typedef typename DarcySolverPolicies::bchandler_type     bchandler_type;

    typedef typename DarcySolverPolicies::matrix_type        matrix_type;
    typedef typename DarcySolverPolicies::matrix_ptrtype     matrix_ptrtype;

    typedef typename DarcySolverPolicies::vector_type        vector_type;
    typedef typename DarcySolverPolicies::vector_ptrtype     vector_ptrtype;

    typedef typename DarcySolverPolicies::comm_type          comm_type;
    typedef typename DarcySolverPolicies::comm_ptrtype       comm_ptrtype;

    //@}

    // Constructors & destructor.
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
    DarcySolverNonLinear ( const data_type&          dataFile,
                           FESpace<Mesh, EpetraMap>& primal_FESpace,
                           FESpace<Mesh, EpetraMap>& dual_FESpace,
                           FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                           FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                           bchandler_raw_type&       bcHandler,
                           comm_ptrtype&             comm );

    /*!
      Constructor for the class without the definition of the boundary handler.
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual finite element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer of the Epetra communicator.
    */
    DarcySolverNonLinear ( const data_type&          dataFile,
                           FESpace<Mesh, EpetraMap>& primal_FESpace,
                           FESpace<Mesh, EpetraMap>& dual_FESpace,
                           FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                           FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                           comm_ptrtype&             comm );

    //! Virtual destructor.
    virtual ~DarcySolverNonLinear ();

    //@}

    // Update & set methos.
    //! @name Update & set methos
    //@{

    /*!
      Set up the linear solver and the preconditioner for the linear system.
    */
    virtual void setup ();

    /*!
      Set the function for the first iteration for the fixed point method. The default is the zero function.
      @param primalZeroIteration The function for the first iteration.
    */
    void setZeroItarationPrimal ( const Function& primalZeroIteration );
   
    /*!
      Set the inverse of diffusion tensor, the default setted inverse of permeability is the identity matrix.
      @param invPerm Inverse of the permeability tensor for the problem.
    */
    inline void setInversePermeability ( const permeability_type& invPerm )
    {

        // Call the standard set inverse permeability
        DarcySolver<Mesh, SolverType>::setInversePermeability( invPerm );

        // Set the dependence of the previous solution in the permeability
        this->M_inversePermeability->setField( M_primalPreviousIteration );
    }

    //@}

    // Get methods.
    //! @name Get methods
    //@{

    /*!
      Returns the number of iteration of the fixed point scheme.
      @return Final number of iterations of the fixed point method as a constant UInt.
    */
    const UInt IterationFixedPoint () const
    {
        return M_iterationFixedPoint;
    }

    /*!
      Returns the residual between two iterations of the fixed point scheme.
      @return Final residual of the fixed point method as a constant Real.
    */
    const Real residualFixedPoint () const
    {
        return M_residualFixedPoint;
    }

    /*!
      Returns the pointer of the primal solution vector at previous step.
      @return Constant vector_type reference of the primal solution at previous step.
    */
    inline const vector_ptrtype& primalPreviousIteration () const
    {
        return M_primalPreviousIteration;
    }

    //@}

    // Solve functions.
    //! @name Solve functions
    //@{

    /*!
      Perform the fixed point scheme.
    */

    void fixedPointScheme ();

    //@}

protected:

    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class.
    */
    virtual void updateVariables ();

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //@}

    // Non-linear stuff.
    //! @name Non-linear stuff
    //@{

    //! The maximum number of iterations for the fixed point method.
    UInt           M_maxIterationFixedPoint;

    //! The current iterations for the fixed point method.
    UInt           M_iterationFixedPoint;

    //! Tollerance for the stopping criteria for the fixed point method.
    Real           M_tollFixedPoint;

    //! The residual between two iteration of the fixed point scheme.
    Real           M_residualFixedPoint;

    //! Primal solution at previous iteration step.
    vector_ptrtype M_primalPreviousIteration;

    //! Primal solution at zero time step.
    Function       M_primalZeroIteration;

    //@}

}; // class DarcySolverNonLinear

//
// IMPLEMENTATION
//

// Complete constructor.
template<typename Mesh, typename SolverType>
DarcySolverNonLinear<Mesh, SolverType>::
DarcySolverNonLinear ( const data_type&           dataFile,
                       FESpace<Mesh, EpetraMap>&  primal_FESpace,
                       FESpace<Mesh, EpetraMap>&  dual_FESpace,
                       FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                       FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                       bchandler_raw_type&        bcHandler,
                       comm_ptrtype&              comm ):
    // Standard Darcy solver constructor.
    DarcySolver<Mesh, SolverType>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, bcHandler, comm),
    // Non-linear stuff.
    M_maxIterationFixedPoint        ( static_cast<UInt>(10) ),
    M_iterationFixedPoint           ( static_cast<UInt>(0) ),
    M_tollFixedPoint                ( static_cast<Real>(1.e-8) ),
    M_residualFixedPoint            ( M_tollFixedPoint + static_cast<Real>(1) ),
    M_primalPreviousIteration       ( new vector_type ( this->M_primal_FESpace.map() ) ),
    M_primalZeroIteration           ( DarcyDefaultStartUpFunction() )

{

    // Interpolate the primal initial value on the default initial value function, for the fixed point scheme.
    this->M_primal_FESpace.interpolate( M_primalZeroIteration,
                                        *(this->M_primal),
                                        static_cast<Real>(0) );

    CONSTRUCTOR( "DarcySolverNonLinear" );

} // Constructor

// Constructor without boundary condition handler.
template<typename Mesh, typename SolverType>
DarcySolverNonLinear<Mesh, SolverType>::
DarcySolverNonLinear ( const data_type&           dataFile,
                       FESpace<Mesh, EpetraMap>&  primal_FESpace,
                       FESpace<Mesh, EpetraMap>&  dual_FESpace,
                       FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
                       FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
                       comm_ptrtype&              comm ):
    // Standard Darcy solver constructor.
    DarcySolver<Mesh, SolverType>::DarcySolver( dataFile, primal_FESpace, dual_FESpace, hybrid_FESpace, VdotN_FESpace, comm),
    // Non-linear stuff.
    M_maxIterationFixedPoint        ( static_cast<UInt>(10) ),
    M_iterationFixedPoint           ( static_cast<UInt>(0) ),
    M_tollFixedPoint                ( static_cast<Real>(1.e-8) ),
    M_residualFixedPoint            ( M_tollFixedPoint + static_cast<Real>(1) ),
    M_primalPreviousIteration       ( new vector_type ( this->M_primal_FESpace.map() ) ),
    M_primalZeroIteration           ( DarcyDefaultStartUpFunction() )
{

    // Interpolate the primal initial value on the default initial value function, for the fixed point scheme.
    this->M_primal_FESpace.interpolate( M_primalZeroIteration,
                                        *(this->M_primal),
                                        static_cast<Real>(0) );

    CONSTRUCTOR( "DarcySolverNonLinear" );

} // Constructor

// Virtual destructor.
template<typename Mesh, typename SolverType>
DarcySolverNonLinear<Mesh, SolverType>::
~DarcySolverNonLinear ()
{

	DESTRUCTOR( "DarcySolverNonLinear" );

} // Destructor

// Set up the linear solver and the preconditioner.
template<typename Mesh, typename SolverType>
void
DarcySolverNonLinear<Mesh, SolverType>::
setup ()
{

    GetPot dataFile( *(this->M_data.dataFile()) );

    // Call the DarcySolver setup method for setting up the linear solver.
    DarcySolver<Mesh, SolverType>::setup();

    // Set the maximum number of iteration for the fixed point iteration scheme.
    M_maxIterationFixedPoint = static_cast<UInt>( dataFile( ( this->M_data.section() + "/non-linear/fixed_point_iteration" ).data(), 10 ) );

    // Set the tollerance for the fixed point iteration scheme.
    M_tollFixedPoint = dataFile( ( this->M_data.section() + "/non-linear/fixed_point_toll" ).data(), 1.e-8 );

} // setup

// Set the first value for the fixed point method. 
template<typename Mesh, typename SolverType>
void 
DarcySolverNonLinear<Mesh, SolverType>::
setZeroItarationPrimal ( const Function& primalZeroIteration )
{
    // Set the function for the first iteration.
    M_primalZeroIteration = primalZeroIteration;

    // Interpolate the primal variable for the first iteration.
    this->M_primal_FESpace.interpolate( M_primalZeroIteration,
                                        *(this->M_primal),
                                        this->M_data.dataTime()->getInitialTime() );
    
} // SetZeroIterationPrimal

// Update all the variables of the problem.
template<typename Mesh, typename SolverType>
void
DarcySolverNonLinear<Mesh, SolverType>::
updateVariables ()
{
    // Reset the primal vector at the previous iteration step.
    M_primalPreviousIteration.reset( new vector_type( this->M_primal_FESpace.map() ) );

    // Update the primal vector at the previous iteration step.
    *M_primalPreviousIteration = *(this->M_primal);

    // Call the method of the DarcySolver to update all the variables defined in it.
    DarcySolver<Mesh, SolverType>::updateVariables();

} // updateVariables

// Fixed point scheme.
template <typename Mesh, typename SolverType>
void
DarcySolverNonLinear<Mesh, SolverType>::
fixedPointScheme ()
{

    // Current iteration.
    M_iterationFixedPoint = static_cast<UInt>(0);

    // Error between two iterations, it is the relative error between two step of the primal vector
    M_residualFixedPoint =  M_tollFixedPoint + 1;

    /* A loop for the fixed point scheme, with exit condition based on stagnate of the
       primal variable and the maximum iteration. */
    while( M_residualFixedPoint > M_tollFixedPoint && M_iterationFixedPoint < M_maxIterationFixedPoint )
    {
        // Increment the iteration number.
        ++M_iterationFixedPoint;

        // Build the linear system and the right hand side.
        this->buildSystem();

        // Solve the linear system.
        this->solve();

        // Post process of the primal and dual variables.
        this->computePrimalAndDual();

        // Compute the error.
        M_residualFixedPoint = ( *(this->M_primal) - *M_primalPreviousIteration ).Norm2() / this->M_primal->Norm2();

        // The leader process prints the iteration data.
        this->M_displayer.leaderPrint( "Fixed point scheme           \n" );

        // Print the maximum number of iterations
        this->M_displayer.leaderPrint( "Maximum number of iterations ",
                                       M_maxIterationFixedPoint, "\n" );

        // Print the actual iterations
        this->M_displayer.leaderPrint( "Iteration number             ",
                                       M_iterationFixedPoint, "\n" );

        // Print the tollerance
        this->M_displayer.leaderPrint( "Tollerance                   ",
                                       M_tollFixedPoint, "\n" );

        // Print the error reached
        this->M_displayer.leaderPrint( "Error                        ",
                                       M_residualFixedPoint, "\n" );

    }

    // Check if the fixed point method reach the tollerance.
    ASSERT( M_maxIterationFixedPoint > M_iterationFixedPoint, "Attention the fixed point scheme does not reach the convergence." );

} // fixedPointScheme

} // namespace LifeV


#endif //_darcySolverNonLinear_H_
