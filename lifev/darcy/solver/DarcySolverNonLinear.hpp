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

 This file is part of the LifeV library

 Author(s): A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
      Date: 2010-05-09

 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
 Copyright (C) 2006-2010 Politecnico di Milano

*******************************************************************************

*/
//@HEADER

/*!
    @file
    @brief This file contains a non-linear permeability term Darcy equation solver class

    @date 05/2010
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    @contributor M. Kern <michel.kern@inria.fr>
    @maintainer M. Kern <michel.kern@inria.fr>
*/

#ifndef _DARCYSOLVERNONLINEAR_H_
#define _DARCYSOLVERNONLINEAR_H_ 1

#include <lifev/darcy/solver/DarcySolverLinear.hpp>

// LifeV namespace.
namespace LifeV
{
//! @class DarcySolverNonLinear This class implements a non-linear Darcy solver
/*!
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    This class implements a non-linear Darcy solver with fixed point scheme to handle the non-linearity.
    <br>
    <br>
    The classical non-linear Darcy formulation is a couple of differential equations of first order with
    the unknowns \f$ p \in C^1 (\Omega ) \f$, being the pressure or the primal unknown,
    and \f$ \sigma \in (C^1( \Omega ) )^n \f$, being the Darcy velocity or the flux or the dual unknown,
    such that
    \f[
    \left\{
    \begin{array}{l l l}
    \Lambda^{-1}(p) \sigma + \nabla p = f_v & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
    \nabla \cdot \sigma + \pi p - f = 0     & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
    p = g_D                                 & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
    \sigma \cdot n + h p = g_R              & \mathrm{on} & \Gamma_R\,.
    \end{array}
    \right.
    \f]
    The data in the system are:
    <ul>
        <li> \f$ \Lambda \f$ the permeability tensor that depends on \f$ p \f$; </li>
        <li> \f$ f \f$ the scalar source term; </li>
        <li> \f$ f_v \f$ the vector source term; </li>
        <li> \f$ \pi \f$ the reaction coefficient; </li>
        <li> \f$ \Gamma_D \f$ the subset of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions; </li>
        <li> \f$ g_D \f$ Dirichlet boundary condition data; </li>
        <li> \f$ \Gamma_R \f$ that is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary conditions; </li>
        <li> \f$ h \f$ and \f$ g_R \f$ Robin boundary condition data; </li>
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
    Introducing the following operator, bilinear forms and the functionals
    \f[
    \begin{array}{l l}
    a(\sigma, \tau, p) = \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1}(p) \sigma \cdot \tau \,,  &
    b(p, \tau) = -\sum_{K \in \mathcal{T}_h} \int_K p \nabla \cdot \tau\,, \vspace{0.2cm}\\
    c(\lambda, \tau) = \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda \tau \cdot n\,,&
    d(p, q) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \pi p q \,,\vspace{0.2cm}\\
    h(\lambda, \mu) = \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda \,,&
    F(q) = \sum_{K \in \mathcal{T}_h} \int_K f q\,,\vspace{0.2cm}\\
    F_v(\tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f_v \tau \,, \vspace{0.2cm}\,,&
    G(\mu) =\sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
    \end{array}
    \f]
    we obtain the Darcy problem in the weak form: find \f$ (\sigma, \, p, \, \lambda) \in Z \times V \times \Lambda \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma, \tau, p) + b(p, \tau) + c(\lambda, \tau) = F_v(\tau)\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
    b(q, \sigma) + d (p, q)= - F(q)\,,                                 & \forall q \in V \,,\vspace{0.2cm}\\
    c(\mu, \sigma) + h(\lambda, \mu) = G(\mu) \,,                      & \forall \mu \in \Lambda\,.
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
    a(\sigma_h, \tau_h, p_h) + b(p_h, \tau_h) + c(\lambda_h, \tau_h) = F_v(\tau)\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    b(q_h, \sigma_h) - d(p_h, q_h) = - F(q_h)\,,                             & \forall q_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h) + h(\lambda_h, \mu_h) = G(\mu_h) \,,                  & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To solve the non-linearity we use a fixed point scheme based on the relative difference between two consecutive iterations
    of the primal variable. We start from the, user defined, function \f$ p^0 \f$ and solve the linearized problem for \f$ n \geq 1 \f$:
    find \f$ (\sigma_h^n,\, p_h^n, \, \lambda_h^n) \in Z_h \times V_h \times \Lambda_h \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h^n, \tau_h; p_h^{n-1}) + b(p_h^n, \tau_h) + c(\lambda_h^n, \tau_h) = F_v(\tau) \,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    b(q_h, \sigma_h^n) - d(p_h, q_h)= - F(q_h)\,,                         & \forall q_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h^n) + h(\lambda_h^n, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    At each iteration, to solve the problem, we use the static condensation procedure, i.e. the unknowns in the discrete
    weak system are not independent and \f$ p_K^n \f$, \f$\sigma_K^n \f$ may be written in function of
    \f$ \lambda_K^n \f$ alone. We introduce the following local matrices
    \f[
    \begin{array}{l l l}
    \left[ A \left( p_K^{n-1}\right) \right]_{ij} = \displaystyle \int_K \Lambda^{-1}\left( p^{n-1} \right) \psi_j \cdot \psi_i \,, &
    \left[ B \right]_{ij} = - \displaystyle \int_K \phi_j \nabla \cdot \psi_i \,, &
    \left[ C \right]_{ij} =   \displaystyle \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
    \left[ D \right]_{ij} =   \displaystyle   \int_K \pi \phi_j \phi_i\,, &
    \left[ H \right]_{ij} =   \displaystyle \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
    \left[ F \right]_{j}  =   \displaystyle \int_K f \phi_j\,,\vspace{0.2cm} \\
    \left[ F_v \right]_{j} = \displaystyle   \int_K f_v \psi_j\,,&
    \left[ G \right]_{j}  =   \displaystyle \int_{\partial K \cap \Gamma_R } g \xi_j\,,
    \end{array}
    \f]
    where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors.
    <br>
    The local matrix formulation of the finite dimensional problem, at each iteration, is
    \f[
    \left\{
    \begin{array}{l}
    A \left( p_K^{n-1} \right) \sigma_K^n + B p_K^n + C \lambda_K^n = F_v\,, \vspace{0.2cm} \\
    B^T \sigma_K^n - D p_K^n= -F \,,                    \vspace{0.2cm}\\
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
    B^T                        & -D & 0 \vspace{0.2cm} \\
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
    F_v \vspace{0.2cm}\\
    -F \vspace{0.2cm}\\
    G
    \end{array}
    \right]\,.
    \end{array}
    \f]
    Introducing the local hybrid matrix and local hybrid right hand side
    \f[
    \begin{array}{l}
    L_K =  -C^T A^{-1}\left( p_k^{n-1}\right) C + C^T A^{-1} \left( p_K^{n-1} \right) B
           \left[ B^T A^{-1}\left( p_K^{n-1}\right) B + D \right]^{-1} B^T A^{-1}\left( p_K^{n-1} \right) C + H \,, \vspace{0.2cm} \\
    r_K = G + C^T A^{-1}\left( p_K^{n-1} \right) B \left[ B^T A^{-1}\left( p_K^{n-1} \right) B + D \right]^{-1} F
          + C^T A^{-1}\left( p_K^{n-1} \right) B \left[ B^T A^{-1}\left( p_K^{n-1} \right) B +
          D \right]^{-1} B^T A^{-1}\left( p_K^{n-1} \right) F_v - C^T A^{-1}\left( p_K^{n-1} \right) F_v\,,
    \end{array}
    \f]
    Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
    \f[
    L^n \lambda^n = r \,.
    \f]
    We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
    \f[
    \begin{array}{l}
    p_K^n = \left[ B^T A^{-1}\left( p_K^{n-1} \right) B + D\right]^{-1} \left[ F + B^T A^{-1}\left( p_K^{n-1}\right) F_v
           - B^T A^{-1}\left( p_K^{n-1}\right) C \lambda_K^n \right]\,, \vspace{0.2cm} \\
    \sigma_K^n = -A^{-1}\left(p_K^{n-1} \right) \left( B p_K^n - F_v + C \lambda_K^n \right) \,.
    \end{array}
    \f]
    @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
    conditions are imposed via BCHandler class.
    @note Example of usage can be found in darcy_nonlinear and darcy_linear.
    Coupled with an hyperbolic solver in impes.
    @todo Insert any scientific publications that use this solver.
    @todo Attention! We have the hypothesis that we use P0 elements for the primal unknown. Change this in a future!
    @todo Add criteria to ensure convergence of the fixed point method.
*/
template < typename MeshType >
class DarcySolverNonLinear :
        virtual public DarcySolverLinear < MeshType >
{

public:

    //! @name Public Types
    //@{

    //! Typedef for mesh template.
    typedef MeshType mesh_Type;

    //! Self typedef.
    typedef DarcySolverNonLinear < mesh_Type > darcySolverNonLinear_Type;

    //! Darcy solver class.
    typedef DarcySolverLinear < mesh_Type > darcySolverLinear_Type;

    //! Typedef for the data type.
    typedef typename darcySolverLinear_Type::data_Type data_Type;

    //! Shared pointer for the data type.
    typedef typename darcySolverLinear_Type::dataPtr_Type dataPtr_Type;

    //! Shared pointer to a scalar value function.
    typedef typename darcySolverLinear_Type::scalarFctPtr_Type scalarFctPtr_Type;

    //! Shared pointer to a matrix value function.
    typedef typename darcySolverLinear_Type::matrixFctPtr_Type matrixFctPtr_Type;

    //! Scalar field.
    typedef typename darcySolverLinear_Type::scalarField_Type scalarField_Type;

    //! Shared pointer to a scalar field.
    typedef typename darcySolverLinear_Type::scalarFieldPtr_Type scalarFieldPtr_Type;

    //@}

    //! @name Constructors & destructor
    //@{

    //! Constructor for the class.
    DarcySolverNonLinear ();

    //! Virtual destructor.
    virtual ~DarcySolverNonLinear () {};

    //@}

    //! @name Methods
    //@{

    //!  Set up the linear solver and the preconditioner for the linear system.
    virtual void setup ();

    //! Solve the problem performing the fixed point scheme.
    virtual void solve ()
    {
        fixedPoint ();
    }

    //@}

    //! @name Set methos
    //@{

    //! Initial guess for fixed point iteration
    /*!
      Set the function for the first iteration for the fixed point method.
      @param primalZeroIteration The function for the first iteration.
    */
    void setPrimalZeroIteration ( const scalarFctPtr_Type& primalZeroIteration );

    //! Set the inverse of diffusion tensor,
    /*!
      The non linearity of the solver is in the diffusion or permeability tensor.
      To handle this type of non linearity we use an iterative scheme and we
      write the tensor depending on the solution at previuos step.
      The non linear tensor now has as an additional scalar field, the last one,
      the solution of the problem at previuos step.
      This scalar field is automatically updated each iteration of the non
      linear solver.
      <br>
      The user should take into account that the value of the primal variable is
      the last field set, by the user, plus one. So the inverse of permeability
      can access to the primal variable, for example, with
      \code
// Inverse of permeability matrix
typedef RegionMesh < LinearTriangle > regionMesh_Type;
class inversePermeability : public FEFunction < regionMesh_Type, MapEpetra, Matrix >

{
public:
    virtual Matrix eval ( const UInt& iElem, const Vector3D& P, const Real& time = 0. ) const;
};
      \endcode
      Its implementation is
      \code
Matrix inversePermeability::eval ( const UInt& iElem, const Vector3D& P, const Real& time ) const
{
    Matrix invK ( 2, 2 );

    Real unkown_n = scalarField(0).eval( iElem, P, time );

    // Fill in of the inversePermeabilityMatrix
    invK ( 0, 0 ) = 1;
    invK ( 0, 1 ) = 0;
    invK ( 1, 0 ) = 0;
    invK ( 1, 1 ) = 1. / ( unkown_n * unkown_n + 1. );

    return invK;
}
      \endcode
      obtaining a tensor with the non linearity in the position \f$ (1,1) \f$
      the function \f$ \left[K^{-1}\right]_{1,1} = u^2 + 1 \f$.
      <br>
      In the previous
      example we have supposed that the permeability has not others scalar fields.
      If the user has set \f$ m \f$ scalar fields prior to calling setInversePermeability, then the code
      should be
      \code
const Real unkown_n = scalarField(m).eval( iElem, P, time );
      \endcode
      to access the primal at previous step.
      @note For an example of the usage see darcy_nonlinear and
      twophase_impes
      @param invPerm Inverse of the permeability tensor for the problem.
    */
    virtual void setInversePermeability ( const matrixFctPtr_Type& invPerm )
    {

        // Call the standard set inverse permeability
        darcySolverLinear_Type::setInversePermeability ( invPerm );

        // Create the primal variable for the first iteration.
        M_primalFieldPreviousIteration.reset ( new scalarField_Type ( this->M_primalField->getFESpacePtr() ) );

        // Set the dependence of the previous solution in the permeability
        this->M_inversePermeabilityFct->addScalarField ( M_primalFieldPreviousIteration );
    }

    /*!
      Set fixed point tolerance
      @param tol requested tolerance
    */
    void setFixedPointTolerance ( const Real& tol )
    {
        M_fixedPointTolerance = tol;
    }

    /*!
      Set maximum number of fixed point iterations permitted
      @param maxit requested maximum
     */
    void setFixedPointMaxIteration ( const UInt& maxit )
    {
        M_fixedPointMaxIteration = maxit;
    }

    //@}

    //! @name Get methods
    //@{

    /*!
      Returns the number of iteration of the fixed point scheme.
      @return Final number of iterations of the fixed point method as a constant UInt.
    */
    UInt fixedPointNumIteration ()
    {
        return M_fixedPointNumIteration;
    }

    //!  Returns the residual between two iterations of the fixed point scheme.
    /*!
      @return Final residual of the fixed point method as a constant Real.
    */
    Real fixedPointResidual ()
    {
        return M_fixedPointResidual;
    }

    //! Returns fixed point tolerance
    /*!
      @return M_fixedPointTolerance
     */
    const Real& fixedPointTolerance () const
    {
        return M_fixedPointTolerance;
    }

    //! Returns fixed point tolerance
    /*!
      @return M_fixedPointTolerance
    */
    Real& fixedPointTolerance ()
    {
        return M_fixedPointTolerance;
    }

    //! Returns maximum number of fixed point iterations allowed
    /*!
      @return max possible number of fixed point iterations
    */
    const UInt& fixedPointMaxIteration () const
    {
        return M_fixedPointMaxIteration;
    }

    //! Returns maximum number of fixed point iterations allowed
    /*!
      @return max possible number of fixed point iterations
    */
    UInt& fixedPointMaxIteration ()
    {
        return M_fixedPointMaxIteration;
    }

    //!  Returns the pointer of the primal solution field at previous step.
    /*!
      @return Constant scalarFieldPtr_Type reference of the primal solution field at previous step.
    */
    const scalarFieldPtr_Type& primalPreviousIterationPtr () const
    {
        return M_primalFieldPreviousIteration;
    }

    //!  Returns the pointer of the primal solution field at previous step.
    /*!
      @return Constant scalarFieldPtr_Type reference of the primal solution field at previous step.
    */
    scalarFieldPtr_Type& primalPreviousIterationPtr ()
    {
        return M_primalFieldPreviousIteration;
    }

    //@}

protected:

    // Methods
    //! @name Protected Methods
    //@{

    //! Update all problem variables
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class.
    */
    virtual void resetVariables ();

    //! Perform the fixed point loop to solve the non-linear problem.
    void fixedPoint ();

    //! Set up the data for the non-linear solver.
    void setupNonLinear ();

    //@}

    //! @name Protected Members
    //@{

    //! Primal solution at previous iteration step.
    scalarFieldPtr_Type M_primalFieldPreviousIteration;

    //@}

private:

    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    DarcySolverNonLinear ( const darcySolverNonLinear_Type& );

    //@}

    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    darcySolverNonLinear_Type& operator= ( const darcySolverNonLinear_Type& );

    //@}

    // Non-linear stuff.
    //! @name Non-linear stuff
    //@{

    //! The maximum number of iterations for the fixed point method.
    UInt M_fixedPointMaxIteration;

    //! The current iterations for the fixed point method.
    UInt M_fixedPointNumIteration;

    //! Tollerance for the stopping criteria for the fixed point method.
    Real M_fixedPointTolerance;

    //! The residual between two iteration of the fixed point scheme.
    Real M_fixedPointResidual;

    //! Primal solution at zero time step.
    scalarFctPtr_Type  M_primalFieldZeroIterationFct;

    //@}

}; // class DarcySolverNonLinear

//
// IMPLEMENTATION
//

// ===================================================
// Constructors & Destructor
// ===================================================

// Complete constructor.
template < typename MeshType >
DarcySolverNonLinear < MeshType >::
DarcySolverNonLinear ( ):
        // Standard Darcy solver constructor.
        darcySolverLinear_Type::DarcySolverLinear  ( ),
        // Non-linear stuff.
        M_fixedPointMaxIteration       ( static_cast<UInt>(10) ),
        M_fixedPointNumIteration       ( static_cast<UInt>(0) ),
        M_fixedPointTolerance          ( static_cast<Real>(1.e-8) ),
        M_fixedPointResidual           ( M_fixedPointTolerance + static_cast<Real>(1.) )
{
} // Constructor

// ===================================================
// Public methods
// ===================================================

// Set up the linear solver and the preconditioner.
template < typename MeshType >
void
DarcySolverNonLinear < MeshType >::
setup ()
{

    // Call the DarcySolverLinear setup method for setting up the linear solver.
    darcySolverLinear_Type::setup();

    // Call the setup for the non-linear data.
    setupNonLinear ();

} // setup

// Set up the non-linear data.
template < typename MeshType >
void
DarcySolverNonLinear < MeshType >::
setupNonLinear ()
{

    const typename darcySolverLinear_Type::data_Type::data_Type& dataFile = *( this->M_data->dataFilePtr() );

    // Path for the non linear stuff in the data file.
    const std::string dataPath = this->M_data->section() + "/non-linear";

    // Set the maximum number of iteration for the fixed point iteration scheme.
    const UInt maxIter = dataFile( ( dataPath + "/fixed_point_iteration" ).data(), 10 );
    setFixedPointMaxIteration ( maxIter );

    // Set the tollerance for the fixed point iteration scheme.
    const Real tol = dataFile( ( dataPath + "/fixed_point_toll" ).data(), 1.e-8 );
    setFixedPointTolerance ( tol );

} // setupNonLinear

// Fixed point scheme.
template < typename MeshType >
void
DarcySolverNonLinear < MeshType >::
fixedPoint ()
{

    // Current iteration.
    M_fixedPointNumIteration = static_cast<UInt>(0);

    // Error between two iterations, it is the relative error between two step of the primal vector
    M_fixedPointResidual = fixedPointTolerance() + 1.;

    /* A loop for the fixed point scheme, with exit condition based on stagnate of the
       primal variable and the maximum iteration. */
    while ( M_fixedPointResidual > M_fixedPointTolerance
            && M_fixedPointNumIteration < M_fixedPointMaxIteration )
    {
        // Increment the iteration number.
        ++M_fixedPointNumIteration;

        // Solve the problem.
        darcySolverLinear_Type::solve ();

        // Compute the error.
        M_fixedPointResidual = ( this->M_primalField->getVector() -
                                 M_primalFieldPreviousIteration->getVector() ).norm2()
                               / this->M_primalField->getVector().norm2();

        // The leader process prints the iteration data.
        this->M_displayer->leaderPrint ( "Fixed point scheme\n" );

        // Print the actual iterations
        this->M_displayer->leaderPrint ( "Number of iterations ", M_fixedPointNumIteration );

        // Print the maximum number of iterations
        this->M_displayer->leaderPrint ( " of ", M_fixedPointMaxIteration, "\n" );

        // Print the error reached
        this->M_displayer->leaderPrint ( "Error ", M_fixedPointResidual );

        // Print the tollerance
        this->M_displayer->leaderPrint ( " over ", M_fixedPointTolerance, "\n" );

    }

    // Check if the fixed point method reach the tolerance.
    if ( M_fixedPointMaxIteration < M_fixedPointNumIteration )
    {
        std::cerr << "Attention the fixed point scheme did not reach convergence."
                  << std::endl << "Max of iterations " << M_fixedPointMaxIteration
                  << std::endl << "Number of iterations " << M_fixedPointNumIteration
                  << std::endl;
        exit(1);
    }

} // fixedPoint

// Set the first value for the fixed point method.
template < typename MeshType >
void
DarcySolverNonLinear < MeshType >::
setPrimalZeroIteration ( const scalarFctPtr_Type& primalZeroIterationFct )
{
    // Set the function for the first iteration.
    M_primalFieldZeroIterationFct = primalZeroIterationFct;

    // Interpolate the primal variable for the first iteration.
    M_primalFieldZeroIterationFct->interpolate ( *(this->M_primalField),
                                                 this->M_data->dataTimePtr()->initialTime() );

} // SetZeroIterationPrimal

// ===================================================
// Protected Methods
// ===================================================

// Update all the variables of the problem.
template < typename MeshType >
void
DarcySolverNonLinear < MeshType >::
resetVariables ()
{

    // Update the primal vector at the previous iteration step.
    M_primalFieldPreviousIteration->setVector( this->M_primalField->getVector() );

    // Call the method of the DarcySolverLinear to update all the variables defined in it.
    darcySolverLinear_Type::resetVariables();

} // resetVariables

} // namespace LifeV

#endif //_DARCYSOLVERNONLINEAR_H_
