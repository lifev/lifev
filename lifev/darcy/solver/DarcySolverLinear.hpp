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
     @brief This file contains a Darcy solver with mixed-hybrid finite elements

     @date 05/2010
     @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

     @contributor M. Kern <michel.kern@inria.fr>
     @maintainer M. Kern <michel.kern@inria.fr>
*/

#ifndef _DARCYSOLVERLINEAR_HPP_
#define _DARCYSOLVERLINEAR_HPP_ 1


#include <Epetra_LAPACK.h>
#include <Epetra_BLAS.h>


#include <lifev/core/algorithm/LinearSolver.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>

#include <lifev/core/util/Displayer.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>

#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FEFunction.hpp>
#include <lifev/core/fem/TimeAdvance.hpp>

#include <lifev/darcy/solver/DarcyData.hpp>

// LifeV namespace.
namespace LifeV
{
//! @class DarcySolverLinear This class implements a mixed-hybrid FE Darcy solver.
/*!
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @see For applications related to two-phase flow see \cite Fumagalli2011a

    This class implements a Darcy solver.
    <br>
    <br>
    The classical formulation of this problem is a couple of differential equations of first order with two
    unknowns: \f$ p \in C^1 (\Omega ) \f$, being the pressure or the primal unknown,
    and \f$ \sigma \in (C^1( \Omega ) )^n \f$, being the Darcy velocity or the total flux or the dual unknown,
    such that
    \f[
    \left\{
    \begin{array}{l l l}
    \Lambda^{-1} \sigma + \nabla p = f_v & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
    \nabla \cdot \sigma + \pi p - f = 0  & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
    p = g_D                              & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
    \sigma \cdot n + h p = g_R           & \mathrm{on} & \Gamma_R\,.
    \end{array}
    \right.
    \f]
    The first equation is the Darcy equation and the second equation is the conservation law.
    <br>
    The data in the system are:
    <ul>
        <li> \f$ \Lambda \f$ the permeability tensor; </li>
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
    \Lambda = \displaystyle \left\{ \lambda \in \prod_{K \in \mathcal{T}_h} H^{1/2} (\partial K): \lambda_K = \lambda_{K'} \
    ,\, \mathrm{on} \,\, e_{K\cap K'} \, \forall K \in \mathcal{T}_h,\, \lambda = g_D \,\, \mathrm{on} \,\, \Gamma_D \right\}\,.
    \end{array}
    \f]
    Introducing the following bilinear forms and functionals
    \f[
    \begin{array}{l l}
    a(\sigma, \tau) = \displaystyle \sum_{K \in \mathcal{T}_h}\int_K \Lambda^{-1} \sigma \cdot \tau \,,  &
    b(p, \tau) = \displaystyle -\sum_{K \in \mathcal{T}_h} \int_K p \nabla \cdot \tau\,, \vspace{0.2cm}\\
    c(\lambda, \tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \lambda \tau \cdot n\,,&
    d(p, q) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K} \pi p q \,,\vspace{0.2cm}\\
    h(\lambda, \mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda \,,&
    F(q) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f q\,,\vspace{0.2cm}\\
    F_v(\tau) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f_v \tau \,, \vspace{0.2cm}\,,&
    G(\mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
    \end{array}
    \f]
    we obtain the Darcy problem in the weak form: find \f$ (\sigma, \, p, \, \lambda) \in Z \times V \times \Lambda \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma, \tau) + b(p, \tau) + c(\lambda, \tau) = F_v(\tau)\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
    b(q, \sigma) - d(p,q) = - F(q)\,,                               & \forall q \in V \,,\vspace{0.2cm}\\
    c(\mu, \sigma) + h(\lambda, \mu) = G(\mu) \,,                   & \forall \mu \in \Lambda\,.
    \end{array}
    \right.
    \f]
    At discrete level we introduce the polynomial space, of degree \f$ r \f$, that approximate the finite dimensional
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
    The finite dimensional problem is: find \f$ (\sigma_h,\, p_h, \, \lambda_h) \in Z_h \times V_h \times \Lambda_h \f$ such that
    \f[
    \left\{
    \begin{array}{l l}
    a(\sigma_h, \tau_h) + b(p_h, \tau_h) + c(\lambda_h, \tau_h) = F_v(\tau_h) \,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
    b(q_h, \sigma_h) - d(p_h, q_h)= -F(q_h)\,,                                     & \forall q_h \in V_h \,,\vspace{0.2cm}\\
    c(\mu_h, \sigma_h) + h(\lambda_h, \mu_h) = G(\mu_h) \,,                        & \forall \mu_h \in \Lambda_h\,.
    \end{array}
    \right.
    \f]
    To solve the problem we use the static condensation procedure, i.e. the unknowns in the discrete
    weak system are not independent and \f$ p_K \f$, \f$\sigma_K \f$ may be written in function of
    \f$ \lambda_K \f$ alone. We introduce the following local matrices
    \f[
    \begin{array}{l l l}
    \left[ A \right]_{ij}  = \displaystyle   \int_K \Lambda^{-1} \psi_j \cdot \psi_i \,, &
    \left[ B \right]_{ij}  = \displaystyle - \int_K \phi_j \nabla \cdot \psi_i \,, &
    \left[ C \right]_{ij}  = \displaystyle   \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
    \left[ D \right]_{ij}  = \displaystyle   \int_K \pi \phi_j \phi_i\,, &
    \left[ H \right]_{ij}  = \displaystyle   \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
    \left[ F \right]_{j}   = \displaystyle   \int_K f \phi_j\,, \vspace{0.2cm} \\
    \left[ F_v \right]_{j} = \displaystyle   \int_K f_v \cdot \psi_j\,,&
    \left[ G \right]_{j}   = \displaystyle   \int_{\partial K \cap \Gamma_R } g \xi_j\,,
    \end{array}
    \f]
    where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors. <BR>
    The local matrix formulation of the finite dimensional problem is
    \f[
    \left\{
    \begin{array}{l}
    A \sigma_K + B p_K + C \lambda_K = F_v\,, \vspace{0.2cm} \\
    B^T \sigma_K - D p_K = -F \,, \vspace{0.2cm}\\
    C^T \sigma_K + H \lambda_K = G\,.
    \end{array}
    \right.
    \f]
    Or alternatively
    \f[
    \begin{array}{l l l}
    \left[
    \begin{array}{c c c}
    A   & B &  C \vspace{0.2cm} \\
    B^T & -D &  0 \vspace{0.2cm} \\
    C^T & 0 & H
    \end{array}
    \right] \, \cdot &
    \left[
    \begin{array}{c}
    \sigma_K  \vspace{0.2cm}\\
    p_K       \vspace{0.2cm}\\
    \lambda_K
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
    L_K = - C^T A^{-1} C + C^T A^{-1} B \left( B^T A^{-1} B + D \right)^{-1} B^T A^{-1} C + H \,, \vspace{0.2cm} \\
    r_K = G + C^T A^{-1} B \left( B^T A^{-1} B + D \right)^{-1} F + C^T A^{-1} B \left( B^T A^{-1} B + D \right)^{-1} B^T A^{-1} F_v - C^T A^{-1} F_v \,,
    \end{array}
    \f]
    Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
    \f[
    L \lambda = r \,.
    \f]
    We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
    \f[
    \begin{array}{l}
    p_K = \left( B^T A^{-1} B + D \right)^{-1} \left( F + B^T A^{-1} F_v - B^T A^{-1} C \lambda_K \right)\,, \vspace{0.2cm} \\
    \sigma_K = -A^{-1} \left( B p_K - F_v + C \lambda_K \right) \,.
    \end{array}
    \f]
    @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
    conditions are imposed via BCHandler class.
    @note Example of usage can be found in darcy_nonlinear and darcy_linear.
    Coupled with an hyperbolic solver in impes.
    @todo Insert any scientific publications that use this solver.
    @bug If the save flag for the exporter is setted to 0 the program fails.
*/

template < typename MeshType >
class DarcySolverLinear
{

public:

    //! @name Public Types
    //@{

    //! Typedef for mesh template.
    typedef MeshType mesh_Type;

    //! Typedef for solver template.
    typedef LinearSolver solver_Type;

    //! Self typedef
    typedef DarcySolverLinear < mesh_Type > darcySolver_Type;

    //! Typedef for the data type.
    typedef DarcyData < mesh_Type > data_Type;

    //! Shared pointer for the data type.
    typedef std::shared_ptr < data_Type > dataPtr_Type;

    //! Boundary condition handler.
    typedef BCHandler bcHandler_Type;

    //! Shared pointer to a boundary condition handler.
    typedef std::shared_ptr < bcHandler_Type > bcHandlerPtr_Type;

    //! Shared pointer to a MPI communicator.
    typedef std::shared_ptr < Epetra_Comm > commPtr_Type;

    //! Map type.
    typedef MapEpetra map_Type;

    //! Shared pointer to a displayer.
    typedef std::shared_ptr < Displayer > displayerPtr_Type;

    //! Finite element space.
    typedef FESpace < mesh_Type, map_Type > fESpace_Type;

    //! Shared pointer to a finite element space.
    typedef std::shared_ptr < fESpace_Type > fESpacePtr_Type;

    //! Scalar field.
    typedef FEScalarField < mesh_Type, map_Type > scalarField_Type;

    //! Shared pointer to a scalar field.
    typedef std::shared_ptr < scalarField_Type > scalarFieldPtr_Type;

    //! Vector field.
    typedef FEVectorField < mesh_Type, map_Type > vectorField_Type;

    //! Shared pointer to a scalar field.
    typedef std::shared_ptr < vectorField_Type > vectorFieldPtr_Type;

    //! Scalar value function.
    typedef FEFunction < mesh_Type, map_Type, Real > scalarFct_Type;

    //! Shared pointer to a scalar value function.
    typedef std::shared_ptr < scalarFct_Type > scalarFctPtr_Type;

    //! Vector value function.
    typedef FEFunction < mesh_Type, map_Type, Vector > vectorFct_Type;

    //! Shared pointer to a vector value function.
    typedef std::shared_ptr < vectorFct_Type > vectorFctPtr_Type;

    //! Matrix value funcion.
    typedef FEFunction < mesh_Type, map_Type, Matrix > matrixFct_Type;

    //! Shared pointer to a matrix value function.
    typedef std::shared_ptr < matrixFct_Type > matrixFctPtr_Type;

    //! Sparse and distributed matrix.
    typedef typename solver_Type::matrix_Type matrix_Type;

    //! Shared pointer to a sparse and distributed matrix.
    typedef typename solver_Type::matrixPtr_Type matrixPtr_Type;

    //! Distributed vector.
    typedef typename solver_Type::vector_Type vector_Type;

    //! Shared pointer to a distributed vector.
    typedef typename solver_Type::vectorPtr_Type vectorPtr_Type;

    //! Shared pointer to the preconditioner.
    typedef typename solver_Type::preconditionerPtr_Type preconditionerPtr_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Constructor for the class.
    DarcySolverLinear () {};

    //! Virtual destructor.
    virtual ~DarcySolverLinear () {};

    //@}

    //! @name Methods
    //@{

    //! Build the global hybrid system, the right hand and apply the boundary conditions.
    void buildSystem ();

    //! Set up the linear solver and the preconditioner for the linear system.
    virtual void setup ();

    //! Solve the global hybrid system.
    void solveLinearSystem ();

    //! Compute primal and dual variables from the hybrid variable as a post process.
    void computePrimalAndDual ();

    //! Solve the Darcy problem grouping other public methods.
    virtual void solve ();

    //@}

    //! @name Set Methods
    //@{

    //! Set the data.
    /*!
      @param data Data for the problem.
    */
    void setData ( dataPtr_Type& data )
    {
        M_data = data;
    }

    //! Set the boundary conditions.
    /*!
      @param bcHandler Boundary condition handler for the problem.
    */
    void setBoundaryConditions ( bcHandlerPtr_Type& bcHandler )
    {
        M_boundaryConditionHandler = bcHandler;
    }

    //! Set scalar source term
    /*!
      @param scalarSourceFct Vector source term for the problem.
    */
    void setScalarSource ( const scalarFctPtr_Type& scalarSourceFct )
    {
        M_scalarSourceFct = scalarSourceFct;
    }

    //! Set vector source term
    /*!
      @param vectorSourceFct Vector source term for the problem.
    */
    void setVectorSource ( const vectorFctPtr_Type& vectorSourceFct )
    {
        M_vectorSourceFct = vectorSourceFct;
    }

    //! Set inverse diffusion tensor
    /*!
      Set the inverse of diffusion tensor, the default inverse of permeability is the identity matrix.
      @param invPermFct Inverse of the permeability tensor for the problem.
    */
    virtual void setInversePermeability ( const matrixFctPtr_Type& invPermFct )
    {
        M_inversePermeabilityFct = invPermFct;
    }

    //! Set the coefficient for the reaction term.
    /*!
      @param reactionTermFct Reaction term for the problem.
      @note Useful also for time discretization term.
    */
    void setReactionTerm ( const scalarFctPtr_Type& reactionTermFct )
    {
        M_reactionTermFct = reactionTermFct;
    }

    //! Set the hybrid field vector.
    /*!
      @param hybrid Constant scalarFieldPtr_Type reference of the hybrid vector.
     */
    void setHybridField ( const scalarFieldPtr_Type& hybridField )
    {
        M_hybridField = hybridField;
    }

    //! Set the primal field vector.
    /*!
      @param primalField Constant scalarFieldPtr_Type reference of the primal solution.
    */
    void setPrimalField ( const scalarFieldPtr_Type& primalField )
    {
        M_primalField = primalField;
    }

    //! Set the dual field vector.
    /*!
      @param dualField Constant scalarFieldPtr_Type reference of the dual solution.
    */
    void setDualField ( const scalarFieldPtr_Type& dualField )
    {
        M_dualField = dualField;
    }

    //! Set the fields.
    /*!
      @param dualField Constant scalarFieldPtr_Type reference of the dual solution.
      @param primalField Constant scalarFieldPtr_Type reference of the primal solution.
      @param hybrid Constant scalarFieldPtr_Type reference of the hybrid vector.
    */
    void setFields ( const scalarFieldPtr_Type& dualField,
                     const scalarFieldPtr_Type& primalField,
                     const scalarFieldPtr_Type& hybridField )
    {
        // Set the dual field.
        setDualField ( dualField );

        // Set the primal field.
        setPrimalField ( primalField );

        // Set the hybrid field.
        setHybridField ( hybridField );
    }

    //! Set the displayer and, implicitly, the communicator.
    /*!
      @param displayer Constant displayerPtr_Type reference of the displayer.
    */
    void setDisplayer ( const displayerPtr_Type& displayer )
    {
        M_displayer = displayer;
    }

    //! Set the communicator and, implicitly, the displayer.
    /*!
      @param comm Constant commPtr_Type reference of the communicator.
    */
    void setCommunicator ( const commPtr_Type& comm )
    {
        M_displayer.reset ( new displayerPtr_Type::element_type ( comm ) );
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns pointer to the hybrid solution field.
    /*!
      @return Constant scalarFieldPtr_Type reference of the hybrid field.
    */
    const scalarFieldPtr_Type& hybridFieldPtr () const
    {
        return M_hybridField;
    }

    //! Returns pointer to the hybrid solution field.
    /*!
      @return scalarFieldPtr_Type reference of the hybrid solution field.
    */
    scalarFieldPtr_Type& hybridFieldPtr ()
    {
        return M_hybridField;
    }

    //! Returns pointer to the primal solution field.
    /*!
      @return Constant scalarFieldPtr_Type reference of the primal solution field.
    */
    const scalarFieldPtr_Type& primalFieldPtr () const
    {
        return M_primalField;
    }

    //! Returns pointer to the primal solution field.
    /*!
      @return scalarFieldPtr_Type reference of the primal solution field.
    */
    scalarFieldPtr_Type& primalFieldPtr ()
    {
        return M_primalField;
    }

    //! Returns pointer to the dual solution field.
    /*!
      @return Constant scalarFieldPtr_Type reference of the dual solution field.
    */
    const scalarFieldPtr_Type& dualFieldPtr () const
    {
        return M_dualField;
    }

    //! Returns pointer to the dual solution field.
    /*!
      @return scalarFieldPtr_Type reference of the dual solution field.
    */
    scalarFieldPtr_Type& dualFieldPtr ()
    {
        return M_dualField;
    }

    /*!
      Returns the pointer of the residual vector.
      @return Constant vectorPtr_Type reference of the residual.
    */
    const vectorPtr_Type& residualPtr () const
    {
        return M_residual;
    }

    /*!
      Returns the pointer of the residual vector.
      @return vectorPtr_Type reference of the residual.
    */
    vectorPtr_Type& residualPtr ()
    {
        return M_residual;
    }

    //! Returns boundary conditions handler.
    /*!
      @return Constant reference of boundary conditions handler.
    */
    const bcHandlerPtr_Type& boundaryConditionHandlerPtr () const
    {
        return M_boundaryConditionHandler;
    }

    //! Returns boundary conditions handler.
    /*!
      @return Reference of boundary conditions handler.
    */
    bcHandlerPtr_Type& boundaryConditionHandlerPtr ()
    {
        return M_boundaryConditionHandler;
    }

    //! Returns Epetra local map.
    /*!
      @return Constant MapEpetra reference of the problem.
    */
    const map_Type& getMap () const
    {
        return M_hybridField->getFESpace().map();
    }

    //! Returns Epetra communicator.
    /*!
      @return Reference of the shared pointer of the communicator of the problem.
    */
    const commPtr_Type& getCommPtr () const
    {
        return M_displayer->comm();
    }

    //@}

protected:

    //! @name Private Constructors
    //@{

    //! Inhibited copy constructor.
    DarcySolverLinear ( const darcySolver_Type& );

    //@}

    //! @name Private Operators
    //@{

    //! Inhibited assign operator.
    darcySolver_Type& operator= ( const darcySolver_Type& );

    //@}

    //! @name Protected Methods
    //@{

    //! Perform all the operations before doing the loop on volume elements.
    /*!
      Computes the element independent matrices, clean all the fields and
      all the algebraich stuff. Useful for extentions in child classes.
    */
    virtual void preLoopElementsComputation ()
    {
        // Update all the variables, e.g. reset the hybrid matrix
        resetVariables ();
    }

    //! Pre-computes local (element independant) matrices
    /*!
      Compute all the local matrices that are independant
      from the geometrical element.
      @param elmatMix The local matrix in mixed form.
    */
    void computeConstantMatrices ( MatrixElemental& elmatMix );

    //! Computes local (element dependent) matrices
    /*!
      Locally update the current finite element for the dual finite element space,
      then compute the Hdiv mass matrix.
      @param iElem Id of the current geometrical element.
      @param elmatMix The local matrix in mixed form.
      @param elmatReactionTerm The local matrix for the reaction term.
    */
    virtual void localMatrixComputation ( const UInt& iElem,
                                          MatrixElemental& elmatMix,
                                          MatrixElemental& elmatReactionTerm );

    //! Computes local (element dependent) vectors
    /*!
      Locally update the current finite element for the primal
      and dual finite element space, then compute the Hdiv and L2 source
      vectors.
      @param iElem Id of the current geometrical element.
      @param elvecMix The local vector in mixed form.
    */
    virtual void localVectorComputation ( const UInt& iElem,
                                          VectorElemental& elvecMix );

    //! Performs static condensation
    /*!
      Locally eliminate pressure and velocity DOFs, create the local
      hybrid matrix and local hybrid right hand side.
      @param localMatrixHybrid The matrix which will store the hybrid local matrix.
      @param localVectorHybrid The vector which will store the hybrid local vector.
      @param elmatMix The local matrix in mixed form.
      @param elmatReactionTerm The local matrix for the reaction term.
      @param elvecMix The local vector in mixed form.
    */
    void staticCondensation ( MatrixElemental& localMatrixHybrid,
                              VectorElemental& localVectorHybrid,
                              MatrixElemental& elmatMix,
                              MatrixElemental& elmatReactionTerm,
                              VectorElemental& elvecMix );

    //! Compute locally, as a post process, the primal and dual variable given the hybrid.
    /*!
      @param localSolution A vector which stores the dual, primal and hybrid local solution.
      @param elmatMix The local matrix in mixed form.
      @param elmatReactionTerm The local matrix for the reaction term.
      @param elvecMix The local vector in mixed form.
    */
    void localComputePrimalAndDual ( VectorElemental& localSolution,
                                     MatrixElemental& elmatMix,
                                     MatrixElemental& elmatReactionTerm,
                                     VectorElemental& elvecMix );

    //! Do some computation after the calculation of the primal and dual variable.
    /*!
      The function is empty while it is useful in derived classes, i.e. DarcySolverTransient.
    */
    virtual void postComputePrimalAndDual () {};

    //! Update all problem variables
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class, in fact we
      want to clear each time step all the variables.
    */
    virtual void resetVariables ();

    //! Apply the boundary condition
    /* Apply BC to the hybrid global matrix and to the hybrid global right hand side.
    */
    void applyBoundaryConditions ();

    //! Make matrix symmetric
    /*!
      Transform a symmetric matrix that is stored only in the lower
      triangular part in a full symmetric matrix.
      @param N The size of the matrix A.
      @tparam A The matrix to be reordered.
    */
    template < typename MatrixType >
    void symmetrizeMatrix ( Int N, MatrixType& A );

    //@}

    // Parallel stuff
    //! @name Parallel stuff
    //@{

    //! Displayer for parallel cout
    displayerPtr_Type M_displayer;

    //@}

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //! Data for Darcy solvers.
    dataPtr_Type M_data;

    //! Source function.
    scalarFctPtr_Type M_scalarSourceFct;

    //! Vector source function.
    vectorFctPtr_Type M_vectorSourceFct;

    //! Reaction term function.
    scalarFctPtr_Type M_reactionTermFct;

    //! Inverse of the permeability tensor.
    matrixFctPtr_Type M_inversePermeabilityFct;

    //! Bondary conditions handler.
    bcHandlerPtr_Type M_boundaryConditionHandler;

    //@}

    // Solution fields stuff.
    //! @name Solution fields stuff
    //@{

    //! Primal solution.
    scalarFieldPtr_Type M_primalField;

    //! Dual solution.
    scalarFieldPtr_Type M_dualField;

    //! Hybrid solution.
    scalarFieldPtr_Type M_hybridField;

    //@}

    // Algebraic stuff.
    //! @name Algebraic stuff
    //@{

    //! Hybrid matrix.
    matrixPtr_Type M_matrHybrid;

    //! Right hand side.
    vectorPtr_Type M_rhs;

    //! Residual.
    vectorPtr_Type M_residual;

    //! Linear solver.
    solver_Type M_linearSolver;

    //! Epetra preconditioner for the linear system.
    preconditionerPtr_Type M_prec;

    //@}

}; // class DarcySolverLinear

//
// IMPLEMENTATION
//

// ===========================================================================================
// Public methods
// ==========================================================================================

// Build the global hybrid matrix and global hybrid right hand side, with boundary conditions.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
buildSystem ()
{

    // Check if the primal field is set or not.
    ASSERT ( M_primalField.get(), "DarcySolverLinear : primal field not set." );

    // Check if the dual field is set or not.
    ASSERT ( M_dualField.get(), "DarcySolverLinear : dual field not set." );

    // Check if the hybrid field is set or not.
    ASSERT ( M_hybridField.get(), "DarcySolverLinear : hybrid field not set." );

    // LifeChronos.
    LifeChrono chronoStaticCondensation;
    LifeChrono chronoAssemble;

    // The total number of elements in the mesh.
    const UInt meshNumberOfElements = M_primalField->getFESpace().mesh()->numElements();

    // Check if the displayer is set or not.
    ASSERT ( M_displayer.get(), "DarcySolverLinear : displayer not set." );

    M_displayer->leaderPrint ( "Perform Static Condensation..." );
    chronoStaticCondensation.start();

    /* Elemental matrix for mixed hybrid matrix, maps [A | B | C] in
      | A    B  C |
      | B^T  0  0 |
      | C^T  0  0 | */
    const UInt primalNbDof = M_primalField->getFESpace().refFE().nbDof();
    const UInt dualNbDof   = M_dualField->getFESpace().refFE().nbDof();
    const UInt hybridNbDof = M_hybridField->getFESpace().refFE().nbDof();

    MatrixElemental elmatMix ( dualNbDof, 1, 1,
                               primalNbDof, 0, 1,
                               hybridNbDof, 0, 1 );

    // Elemental vector for the source terms of the problem, lives in RTk(fe) + Qk(fe).
    VectorElemental elvecMix ( dualNbDof, 1,
                               primalNbDof, 1 );

    // Elemental matrix for the reaction term.
    MatrixElemental elmatReactionTerm ( primalNbDof, 1, 1 );

    // Elemental matrix stores the local hybrid matrix.
    MatrixElemental localMatrixHybrid ( hybridNbDof, 1, 1 );

    // Elemental vector stores the local hybrid right hand side.
    VectorElemental localVectorHybrid ( hybridNbDof, 1 );

    // Compute all the constant matrices, e.g. the matrix B and C
    computeConstantMatrices ( elmatMix );

    // Prepare all the stuff before the loop on all the volume elements.
    preLoopElementsComputation ();

    //! Loop on all the volume elements.
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {

        // Clear the local hybrid matrix and the local hybrid right hand side.
        localMatrixHybrid.zero();
        localVectorHybrid.zero();

        // Compute the Hdiv mass matrix as a local matrix depending on the current element.
        localMatrixComputation ( iElem, elmatMix, elmatReactionTerm );

        // Compute the source vectors as a local vectors depending on the current element.
        localVectorComputation ( iElem, elvecMix );

        // Perform the static condensation to compute the local hybrid matrix and the local hybrid right hand side.
        staticCondensation ( localMatrixHybrid, localVectorHybrid,
                             elmatMix, elmatReactionTerm, elvecMix );

        /* Assemble the global hybrid matrix.
           M_primal_FESpace is used instead of M_hybridField_FESpace for currentLocalId,
           because currentFE cannot store a ReferenceFEHybrid. */
        assembleMatrix ( *M_matrHybrid,
                         M_primalField->getFESpace().fe().currentLocalId(),
                         localMatrixHybrid,
                         hybridNbDof,
                         M_hybridField->getFESpace().dof(),
                         0, 0, 0, 0 );

        /* Assemble the global hybrid right hand side.
           M_primal_FESpace is used instead of M_hybridField_FESpace for currentLocalId,
           because currentFE cannot store a ReferenceFEHybrid. */
        assembleVector ( *M_rhs,
                         M_primalField->getFESpace().fe().currentLocalId(),
                         localVectorHybrid,
                         hybridNbDof,
                         M_hybridField->getFESpace().dof(), 0 );
    }
    //! End of loop volume operation.

    chronoStaticCondensation.stop();
    M_displayer->leaderPrintMax ( " done in " , chronoStaticCondensation.diff() );

    M_displayer->leaderPrint ( "Apply boundary conditions and assemble global matrix and vector..." );

    chronoAssemble.start();

    // Apply boundary conditions to the hybrid global matrix and to the hybrid global right hand side.
    applyBoundaryConditions();

    // Assemble the global hybrid matrix.
    M_matrHybrid->globalAssemble();

    // Assemble the global hybrid vector.
    M_rhs->globalAssemble();

    chronoAssemble.stop();

    M_displayer->leaderPrintMax ( " done in " , chronoAssemble.diff() );

} // buildSystem

// Set up the linear solver and the preconditioner.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
setup ()
{
    // Setup the linear solver.
    M_linearSolver.setParameters ( M_data->linearSolverList() );
    M_linearSolver.setCommunicator ( M_displayer->comm() );

    // Choose the preconditioner type.
    const std::string precType = M_data->preconditionerList().template get<std::string> ( "prectype" );

    // Create a preconditioner object.
    M_prec.reset ( PRECFactory::instance().createObject ( precType ) );
    ASSERT ( M_prec.get() != 0, "DarcySolverLinear : Preconditioner not set" );

    // Set the data for the preconditioner.
    M_prec->setParametersList ( M_data->preconditionerList().sublist ( precType ) );
} // setup

// Solve the linear system.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
solveLinearSystem ()
{

    // Set the matrix.
    M_linearSolver.setOperator ( M_matrHybrid );

    // Set the righthand side.
    M_linearSolver.setRightHandSide ( M_rhs );

    // Set and build the preconditioner.
    M_linearSolver.setPreconditioner ( M_prec );
    M_linearSolver.buildPreconditioner ();

    // Create the solution vector, it has to be of Unique type.
    vectorPtr_Type solution ( new vector_Type ( M_hybridField->getFESpace().map(), Unique ) );

    // Solve the linear system.
    M_linearSolver.solve ( solution );

    // Save the solution into the hybrid variable.
    M_hybridField->setVector ( *solution );

} // solveLinearSystem

// Compute primal and dual unknown as a post process.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
computePrimalAndDual ()
{

    // LifeChrono.
    LifeChrono chronoComputePrimalAndDual;

    M_displayer->leaderPrint ( "Compute pressure and flux..." );
    chronoComputePrimalAndDual.start();

    // The total number of elements in the mesh.
    const UInt meshNumberOfElements = M_primalField->getFESpace().mesh()->numElements();

    // Number of local faces for the geometrical element
    const UInt elementNumberFacets = mesh_Type::element_Type::S_numFacets;

    /*==========================================
      POST PROCESSING
      compute the pressure (Qk or Pk / element)
      and the velocity (RTk / element) => 2 (opposite) velocities / face
      ==========================================*/

    /* The vector M_hybridField can be scalar field with an unique epetra vector, so it does not share one layer elements.
       The reconstruction of the primal and dual variabile need this share, so we create a repeated
       epetra vector. This operation maximize the efficiency of send/receive time cost instead to
       send and receive each time a single datum. */
    vector_Type hybrid_Repeated ( M_hybridField->getVector(), Repeated );

    // Clean the vector for the primal and the dual variable
    M_primalField->cleanField();
    M_dualField->cleanField();

    /* Elemental matrix for mixed hybrid matrix, maps [A | B | C] in
     | A    B  C |
     | B^T  0  0 |
     | C^T  0  0 | */
    const UInt primalNbDof = M_primalField->getFESpace().refFE().nbDof();
    const UInt dualNbDof   = M_dualField->getFESpace().refFE().nbDof();
    const UInt hybridNbDof = M_hybridField->getFESpace().refFE().nbDof();

    MatrixElemental elmatMix ( dualNbDof, 1, 1,
                               primalNbDof, 0, 1,
                               hybridNbDof, 0, 1 );

    // Elemental vector for the source terms of the problem, lives in RTk(fe) + Qk(fe).
    VectorElemental elvecMix ( dualNbDof, 1,
                               primalNbDof, 1 );

    // Elemental matrix for the reaction term.
    MatrixElemental elmatReactionTerm ( primalNbDof, 1, 1 );

    // Compute all the constant matrices, e.g. the matrix B and C
    computeConstantMatrices ( elmatMix );

    // Element vector stores the local solution: (dual, primal, hybrid).
    VectorElemental localSolution ( dualNbDof, 1,
                                    primalNbDof, 1,
                                    hybridNbDof, 1 );

    //! Loop on all the volume elements.
    for ( UInt iElem (0); iElem < meshNumberOfElements; ++iElem )
    {
        // Clear the local solution vector.
        localSolution.zero();

        // Compute the Hdiv mass matrix as a local matrix depending on the current element.
        localMatrixComputation ( iElem,  elmatMix, elmatReactionTerm );

        // Compute the source vectors as local vector depending on the current element.
        localVectorComputation ( iElem, elvecMix );

        // Extract the computed hybrid variable for the current finite element and put it into localHybrid.
        extract_vec ( hybrid_Repeated,
                      localSolution,
                      M_hybridField->getFESpace().refFE (),
                      M_hybridField->getFESpace().dof (),
                      M_primalField->getFESpace().fe().currentLocalId (), 2 );

        // Given the local hybrid variable, computes locally the primal and dual variable.
        localComputePrimalAndDual ( localSolution, elmatMix, elmatReactionTerm, elvecMix );

        // Put the primal variable of the current finite element in the global vector M_primalField.
        assembleVector ( M_primalField->getVector (),
                         M_primalField->getFESpace().fe().currentLocalId (),
                         localSolution,
                         primalNbDof,
                         M_primalField->getFESpace().dof (), 1 );

        for ( UInt iLocalFacet (0); iLocalFacet < elementNumberFacets; ++iLocalFacet )
        {
            const UInt iGlobalFacet ( M_dualField->getFESpace().mesh()->localFacetId ( iElem, iLocalFacet ) );
            if ( M_dualField->getFESpace().mesh()->facet ( iGlobalFacet ).firstAdjacentElementIdentity() != iElem )
            {
                localSolution.block ( 0 ) [ iLocalFacet ] = 0.;
            }
        }

        // Put the dual variable of the current finite element in the global vector M_dualField.
        assembleVector ( M_dualField->getVector (),
                         M_dualField->getFESpace().fe().currentLocalId (),
                         localSolution,
                         dualNbDof,
                         M_dualField->getFESpace().dof (), 0 );

    }
    //! End of loop on the volume elements.

    // Assemble the primal variable.
    M_primalField->getVector().globalAssemble ();

    // It is wrong to assemble the dual variable, no communications required.

    postComputePrimalAndDual ();

    chronoComputePrimalAndDual.stop ();
    M_displayer->leaderPrintMax ( " done in " , chronoComputePrimalAndDual.diff() );

} // computePrimalAndDual

// Solve the Darcy problem.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
solve ()
{
    // Create the hybrid system using the static condensation.
    buildSystem ();

    // Solve the linear system.
    solveLinearSystem ();

    // Given the hybrid solution computes the primal and dual solutions.
    computePrimalAndDual ();

} // solve

// ==============================================================================
// Protected methods
// ==============================================================================

// Compute all the constant matrices.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
computeConstantMatrices ( MatrixElemental& elmatMix )
{

    // Clean the local matrix which will store the matrix B.
    elmatMix.block ( 0, 1 ) = static_cast<Real> (0.);

    // Clean the local matrix which will store the matrix C.
    elmatMix.block ( 0, 2 ) = static_cast<Real> (0.);

    /* Update the divergence matrix, it is independent of the current element
       thanks to the Piola transform. */
    AssemblyElemental::grad_Hdiv ( static_cast<Real> (1.), elmatMix, M_dualField->getFESpace().fe(),
                                   M_primalField->getFESpace().fe(), 0, 1 );

    // Select the correct element which represent ( RT0 \cdot N ) * Hybrid.
    const ReferenceFEHybrid* feRT0VdotNHyb = 0;

    // Selection, in the opt mode the switch is suppressed.
    switch ( mesh_Type::elementShape_Type::S_shape )
    {
        case TRIANGLE:
            feRT0VdotNHyb = &feTriaRT0VdotNHyb;
            break;
        case TETRA:
            feRT0VdotNHyb = &feTetraRT0VdotNHyb;
            break;
        case HEXA:
            feRT0VdotNHyb = &feHexaRT0VdotNHyb;
            break;
        default:
            ASSERT ( false, "DarcySolverLinear : Hybrid finite element not found." );
    }

    /* Update the boundary matrix, it is independent of the current element
       thanks to the Piola transform.
       Here we use the fact that a RefHybridFE "IS A" ReferenceFE and the method refFE of
       a FESpace object. In fact the method refFE return a const ReferenceFE&, but the function
       TP_VdotN_Hdiv takes two const RefHybridFE& so we must cast a const ReferenceFE&
       to a const RefHybridFE&. The cast of type is static and uses pointers. */
    AssemblyElemental::TP_VdotN_Hdiv ( 1., elmatMix,
                                       *static_cast < const ReferenceFEHybrid* > (& ( M_hybridField->getFESpace().refFE() ) ),
                                       *feRT0VdotNHyb, 0, 2 );

} // computeConstantMatrices

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
localMatrixComputation ( const UInt& iElem, MatrixElemental& elmatMix,
                         MatrixElemental& elmatReactionTerm )
{
    // Element of current id.
    const typename mesh_Type::element_Type& element = M_primalField->getFESpace().mesh()->element ( iElem );

    /* Modify the (0,0) block (A) of the matrix elmatMix. The blocks (0,1) (B)
       and (0,2) (C) are independent of the element and have already been computed. */

    // Only one block is set to zero since the other ones are not recomputed.
    elmatMix.block ( 0, 0 ) = 0.;

    // Update the current element of ID iElem for the dual variable.
    M_dualField->getFESpace().fe().update ( element, UPDATE_PHI_VECT | UPDATE_WDET );

    // Get the coordinate of the barycenter of the current element of ID iElem.
    Vector3D barycenter;
    M_dualField->getFESpace().fe().barycenter ( barycenter[0], barycenter[1], barycenter[2] );

    // Check if the inverse of permeability is set or not.
    ASSERT ( M_inversePermeabilityFct.get(), "DarcySolverLinear : inverse of the permeability tensor not set." );

    // Computes the value of the permeability tensor.
    const Matrix permeabilityValue = M_inversePermeabilityFct->eval ( iElem, barycenter, M_data->dataTimePtr()->time() );

    /* Compute the Hdiv mass matrix. We pass the time at the inverse of the permeability
       because the DarcySolverTransient needs the permeability time dependent. In this case
       we do not have a time evolution. */
    AssemblyElemental::mass_Hdiv ( permeabilityValue, elmatMix, M_dualField->getFESpace().fe(), 0, 0 );

    // Check if the reaction term is set or not.
    ASSERT ( M_reactionTermFct.get(), "DarcySolverLinear : reaction term not set." );

    // Clean the reaction matrix.
    elmatReactionTerm.zero();

    // Computes the value for the reaction term.
    const Real reactionValue = M_reactionTermFct->eval ( iElem, barycenter, M_data->dataTimePtr()->time() );

    // Update the current element of ID iElem for the primal variable.
    M_primalField->getFESpace().fe().update ( element, UPDATE_WDET );

    // Compute the reaction matrix.
    AssemblyElemental::mass ( reactionValue, elmatReactionTerm, M_primalField->getFESpace().fe(), 0, 0 );

} // localMatrixComputation

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
localVectorComputation ( const UInt& iElem, VectorElemental& elvecMix )
{
    // Element of current id.
    const typename mesh_Type::element_Type& element = M_primalField->getFESpace().mesh()->element ( iElem );

    // Update the current element of ID iElem only for the dual variable.
    M_dualField->getFESpace().fe().update ( element, UPDATE_PHI_VECT | UPDATE_WDET );

    // Get the coordinate of the barycenter of the current element of ID iElem.
    Vector3D barycenter;
    M_dualField->getFESpace().fe().barycenter ( barycenter[0], barycenter[1], barycenter[2] );

    // Clear the source vector.
    elvecMix.zero();

    // Check if the vector source term is set or not.
    ASSERT ( M_vectorSourceFct.get(), "DarcySolverLinear : vector source term not set." );

    // Computes the value of the vector source.
    const Vector vectorSourceValue = M_vectorSourceFct->eval ( iElem, barycenter, M_data->dataTimePtr()->time() );

    // Compute the vector source term.
    AssemblyElemental::source_Hdiv ( vectorSourceValue, elvecMix, M_dualField->getFESpace().fe(), 0 );

    // Check if the scalar source term is set or not.
    ASSERT ( M_scalarSourceFct.get(), "DarcySolverLinear : scalar source term not set." );

    // Computes the scalar source.
    const Real scalarSourceValue = M_scalarSourceFct->eval ( iElem, barycenter, M_data->dataTimePtr()->time() );

    /* Update the current element of ID iElem only for the primal variable,
       it is used for computing the source term. */
    M_primalField->getFESpace().fe().update ( element, UPDATE_WDET );

    // Compute the scalar source term.
    AssemblyElemental::source ( scalarSourceValue, elvecMix, M_primalField->getFESpace().fe(), 1 );

} // localVectorComputation

// Perform the static condensation for the local hybrid matrix.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
staticCondensation ( MatrixElemental& localMatrixHybrid,
                     VectorElemental& localVectorHybrid,
                     MatrixElemental& elmatMix,
                     MatrixElemental& elmatReactionTerm,
                     VectorElemental& elvecMix  )
{

    // LAPACK wrapper of Epetra.
    Epetra_LAPACK lapack;

    // BLAS wrapper of Epetra.
    Epetra_BLAS blas;

    // Flags for the BLAS and LAPACK routine.
    Int INFO[1] = {0};

    // Number of columns of the right hand side := 1.
    const Int NBRHS = 1;
    // Primal variable degrees of freedom.
    const Int primalNbDof = M_primalField->getFESpace().refFE().nbDof();
    // Dual variable degrees of freedom.
    const Int dualNbDof = M_dualField->getFESpace().refFE().nbDof();
    // Hybrid variable degree of freedom.
    const Int hybridNbDof = M_hybridField->getFESpace().refFE().nbDof();

    const Real ONE = 1.0;
    const Real MINUSONE = -1.0;
    const Real ZERO = 0.0;

    // Parameter that indicate the Lower storage of matrices.
    const char UPLO = 'L';

    // Paramater that indicate the Transpose of matrices.
    const char TRANS = 'T';
    const char NOTRANS = 'N';

    // Parameter that indicates whether the matrix has diagonal unit ('N' means no).
    const char NODIAG = 'N';

    // Create and assign the local matrices A, B and C.
    MatrixElemental::matrix_type A = elmatMix.block ( 0, 0 );
    MatrixElemental::matrix_type B = elmatMix.block ( 0, 1 );
    MatrixElemental::matrix_type C = elmatMix.block ( 0, 2 );

    // Create and assign the local vectors fv and fp.
    VectorElemental::super fv = elvecMix.block ( 0 );
    VectorElemental::super fp = elvecMix.block ( 1 );

    // Create local matrices.
    MatrixElemental::matrix_type BtB ( primalNbDof, primalNbDof );
    MatrixElemental::matrix_type CtC ( hybridNbDof, hybridNbDof );
    MatrixElemental::matrix_type BtC ( primalNbDof, hybridNbDof );

    // Clean the local matrices.
    BtB *= 0.;
    CtC *= 0.;
    BtC *= 0.;

    //! Matrix operations.
    /* Put in A the matrix L and L^T, where L and L^T is the Cholesky factorization of A.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    lapack.POTRF ( UPLO, dualNbDof, A, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack factorization of A is not achieved." );

    /* Put in B the matrix L^{-1} * B, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, primalNbDof, A, dualNbDof, B, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation B = L^{-1} B  is not achieved." );

    /* Put in C the matrix L^{-1} * C, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, hybridNbDof, A, dualNbDof, C, hybridNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation C = L^{-1} C  is not achieved." );

    /* Put in BtB the matrix  B^T * L^{-T} * L^{-1} * B = B^T * A^{-1} * B
       BtB stored only on lower part.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/BLAS/SRC/dsyrk.f */
    blas.SYRK ( UPLO, TRANS, primalNbDof, dualNbDof, ONE, B, dualNbDof, ZERO, BtB, primalNbDof );

    /* Put in BtB the matrix
       BtB + M_elmatReactionTerm = B^T * A^{-1} * B + elmatReactionTerm
       BtB stored only on lower part. */
    BtB += elmatReactionTerm.mat();

    /* Put in CtC the matrix C^T * L^{-T} * L^{-1} * C = C^T * A^{-1} * C
       CtC stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f  */
    blas.SYRK ( UPLO, TRANS, hybridNbDof, dualNbDof, ONE, C, dualNbDof, ZERO, CtC, hybridNbDof );

    /* Put in BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, NOTRANS, primalNbDof, hybridNbDof, dualNbDof, ONE, B, dualNbDof, C, dualNbDof,
                ZERO, BtC, primalNbDof );

    /* Put in BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B + elmatReactionTerm.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f  */
    lapack.POTRF ( UPLO, primalNbDof, BtB, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack factorization of BtB is not achieved." );

    /* Put in BtC the matrix LB^{-1} * BtC = LB^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, primalNbDof, hybridNbDof, BtB, primalNbDof, BtC, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    /* Put in CtC the matrix -CtC + BtC^T * BtC
       Result stored only on lower part, the matrix CtC stores
       -C^T * A^{-1} * C + C^T * A^{-t} * B * ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f  */
    blas.SYRK ( UPLO, TRANS, hybridNbDof, primalNbDof, ONE, BtC, primalNbDof, MINUSONE, CtC, hybridNbDof );

    //! End of matrix operations.

    /*
       Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       B^T stores LB and LB^T where LB and LB^T is the factorization of B^T * A^{-1} * B + elmatReactionTerm
       BtC stores LB^{-1} * B^T * A^{-1} * C
       CtC stores -C^T * A^{-1} * C + C^T * A^{-t} * B * ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * B^T * A^{-1} * C
    */

    //! Vector operations.

    /* Put in fp the vector LB^{-1} * fp = LB^{-1} Fp
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, primalNbDof, NBRHS, BtB, primalNbDof, fp, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fp = LB^{-1} fp is not achieved." );

    /* Put in localVectorHybrid the vector BtC^T * fp =
       C^T * A^{-1} * B^T * ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * Fp
       localVectorHybrid is fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, NOTRANS, hybridNbDof, NBRHS, primalNbDof, ONE, BtC, primalNbDof, fp, primalNbDof,
                ZERO, localVectorHybrid, hybridNbDof );

    /* Put in fv the vector L^{-1} * fv = L^{-1} * Fv, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, NBRHS, A, dualNbDof, fv, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fv = L^{-1} fv is not achieved." );

    /* Put in localVectorHybrid the vector - C^T * fv + localVectorHybrid =
       = C^T * A^{-1} * ( B^T * ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * Fp - Fv )
       localVectorHybrid is fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, NOTRANS, hybridNbDof, NBRHS, hybridNbDof, MINUSONE, C, hybridNbDof, fv, dualNbDof,
                ONE, localVectorHybrid, hybridNbDof );

    /* Put in fp the vector B^T * L^{-T} * fv =  B^T * A^{-1} * Fv
       fp fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, TRANS, primalNbDof, NBRHS, dualNbDof, ONE, B, dualNbDof, fv, NBRHS, ZERO, fp, NBRHS );

    /* Put in fp the vector LB^{-1} * fp = LB^{-1} * B^T * A^{-1} * Fv
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, primalNbDof, NBRHS, BtB, primalNbDof, fp, NBRHS, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fp = LB^{-1} rhs is not achieved." );

    /* Put in M_elvecHyb the vector BtC^T * fp + localVectorHybrid =
       C^T * A^{-1} * [ B^T * ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * ( B^T * A^{-1} + Fp ) - Fv ]
       localVectorHybrid is fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, NOTRANS, hybridNbDof, NBRHS, primalNbDof, ONE, BtC, primalNbDof, fp, NBRHS, ONE,
                localVectorHybrid, hybridNbDof );

    //! End of vector operations.

    /* Previously the matrix CtC is stored only in the lower part, but at the moment there is not
       a function assembleMatrix that store a lower triangular sparse matrix.
       Remind to correct these line in the future. */
    symmetrizeMatrix ( hybridNbDof, CtC );

    // Update the hybrid element matrix.
    localMatrixHybrid.block (0, 0) = CtC;

} // staticCondensation

// Locally compute the primal and dual variable.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
localComputePrimalAndDual ( VectorElemental& localSolution,
                            MatrixElemental& elmatMix,
                            MatrixElemental& elmatReactionTerm,
                            VectorElemental& elvecMix )
{

    // LAPACK wrapper of Epetra
    Epetra_LAPACK lapack;

    // BLAS wrapper of Epetra
    Epetra_BLAS blas;

    // Flags for the BLAS and LAPACK routine.

    Int INFO[1] = {0};

    // Number of columns of the right hand side := 1.
    const Int NBRHS = 1;
    // Primal variable degrees of freedom.
    const Int primalNbDof = M_primalField->getFESpace().refFE().nbDof();
    // Dual variable degrees of freedom.
    const Int dualNbDof = M_dualField->getFESpace().refFE().nbDof();
    // Hybrid variable degree of freedom.
    const Int hybridNbDof = M_hybridField->getFESpace().refFE().nbDof();

    const Real ONE = 1.0;
    const Real MINUSONE = -1.0;
    const Real ZERO = 0.0;

    // Parameter that indicate the Lower storage of matrices.
    const char UPLO = 'L';
    // Parameter that indicate the Transpose of matrices.
    const char TRANS = 'T';
    const char NOTRANS = 'N';
    // Parameter that indicates whether the matrix has diagonal unit ('N' means no)
    const char NODIAG = 'N';

    // No need for CtC in this part, and last dsyrk, the only differences.
    MatrixElemental::matrix_type A = elmatMix.block ( 0, 0 );
    MatrixElemental::matrix_type B = elmatMix.block ( 0, 1 );
    MatrixElemental::matrix_type C = elmatMix.block ( 0, 2 );

    VectorElemental::super fv = elvecMix.block ( 0 );
    VectorElemental::super fp = elvecMix.block ( 1 );

    // Matrix Operations

    // Create local matrices.
    MatrixElemental::matrix_type BtB ( primalNbDof, primalNbDof );
    MatrixElemental::matrix_type BtC ( primalNbDof, hybridNbDof );

    // Clean the local matrices.
    BtB *= 0.;
    BtC *= 0.;

    // Matrix Operations

    /* Put in A the matrix L and L^T, where L and L^T is the Cholesky factorization of A.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    lapack.POTRF ( UPLO, dualNbDof, A, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack factorization of A is not achieved." );

    /* Put in B the matrix L^{-1} * B, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, primalNbDof, A, dualNbDof, B, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation B = L^{-1} B  is not achieved." );

    /* Put in C the matrix L^{-1} * C, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, hybridNbDof, A, dualNbDof, C, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation C = L^{-1} C  is not achieved." );

    /* Put in BtB the matrix  B^T * L^{-T} * L^{-1} * B = B^T * A^{-1} * B
       BtB stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f */
    blas.SYRK ( UPLO, TRANS, primalNbDof, dualNbDof, ONE, B, dualNbDof, ZERO, BtB, primalNbDof );

    /* Put in BtB the matrix
       BtB + elmatReactionTerm = B^T * A^{-1} * B + elmatReactionTerm
       BtB stored only on lower part. */
    BtB += elmatReactionTerm.mat();

    /* Put in BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, NOTRANS, primalNbDof, hybridNbDof, dualNbDof, ONE, B, dualNbDof, C,
                dualNbDof, ZERO, BtC, primalNbDof );

    /* Put in BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B + elmatReactionTerm.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    lapack.POTRF ( UPLO, primalNbDof, BtB, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack factorization of BtB is not achieved." );

    /* Put in BtC the matrix LB^{-1} * BtC = LB^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, primalNbDof, hybridNbDof, BtB, primalNbDof, BtC, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    //! End of matrix operations.

    /* Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       BtB stores LB and LB^T where LB and LB^T is the factorization of
             B^T * A^{-1} * B + elmatReactionTerm
       BtC stores LB^{-1} * B^T * A^{-1} * C
     */

    //! Vector operations, computation of primal and dual variable.

    //! Computation of the primal variable

    // Save in localDual the values of fv useful for the dual variable.
    localSolution.block ( 0 ) = fv;

    /* Put in fv the vector L^{-1} * fv, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, NBRHS, A, dualNbDof, fv, dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fv = L^{-1} fv is not achieved." );

    /* Put in fp the vector B^T * L^{-T} * fv + fp =  B^T * A^{-1} * fv + fp
       fp fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( TRANS, TRANS, primalNbDof, NBRHS, dualNbDof, ONE, B, dualNbDof, fv, NBRHS, ONE, fp, NBRHS );

    /* Put in fp the vector LB^{-1} * fp = LB^{-1} * ( B^T * A^{-1} * fv + fp )
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, primalNbDof, NBRHS, BtB, primalNbDof, fp, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fp = LB^{-1} fp is not achieved." );

    /* Put in fp the vector
       BtC * localHybrid + fp = LB^{-1} * ( B^T * A^{-1} * C * lambda_K - B^T * A^{-1} * fv - fp )
       For more details see http://www.netlib.org/blas/dgemm.f */
    blas.GEMM ( NOTRANS, NOTRANS, primalNbDof, NBRHS, hybridNbDof, MINUSONE, BtC, primalNbDof,
                localSolution.block ( 2 ), hybridNbDof, MINUSONE, fp, primalNbDof );

    /* Put in fp the vector LB^{-T} * fp, where fp stores
       - ( B^T * A^{-1} * B + elmatReactionTerm )^{-1} * [ B^T * A^{-1} * ( C * lambda_K - fv ) - fp ]
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, TRANS, NODIAG, primalNbDof, NBRHS, BtB, primalNbDof, fp, primalNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation fp = LB^{-T} fp is not achieved." );

    // Copy the local primal stored in fp into the localPrimal.
    localSolution.block ( 1 ) = fp;

    //! Computation of the dual variable.

    /* Put in localDual the vector L^{-1} * fv, solving a triangular system.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, NOTRANS, NODIAG, dualNbDof, NBRHS, A, dualNbDof, localSolution.block ( 0 ), dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation localDual = L^{-1} localDual is not achieved." );

    /* Put in localDual the vector B * localPrimal - localDual =  L^{-1} * ( B * primal_K - fv )
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    blas.GEMV ( NOTRANS, dualNbDof, primalNbDof, ONE, B, dualNbDof, localSolution.block ( 1 ),
                MINUSONE, localSolution.block ( 0 ) );

    /* Put in localDual the vector - C * localHybrid - localDual =
       = - L^{-1} * ( C * lambda_K + B^T * primal_K - fv )
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    blas.GEMV ( NOTRANS, dualNbDof, hybridNbDof, MINUSONE, C, hybridNbDof, localSolution.block ( 2 ),
                MINUSONE, localSolution.block ( 0 ) );

    /* Put in localDual the vector L^{-T} * localDual =
       = - A^{-1} ( C^T * lambda_K - B^T * primal_K + fv )
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    lapack.TRTRS ( UPLO, TRANS, NODIAG, dualNbDof, NBRHS, A, dualNbDof, localSolution.block ( 0 ), dualNbDof, INFO );
    ASSERT_PRE ( !INFO[0], "Lapack Computation localDual = L^{-T} localDual is not achieved.");

    //! End of vector computation.

} // localComputePrimalAndDual

// Update all the variables of the problem.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
resetVariables ()
{

    const map_Type& problemMap = M_hybridField->getFESpace().map();

    // Reset the global hybrid matrix.
    M_matrHybrid.reset ( new matrix_Type ( problemMap ) );

    // Reset the global right hand side vector
    M_rhs.reset ( new vector_Type ( problemMap ) );

    // Reset the global residual vector
    M_residual.reset ( new vector_Type ( problemMap ) );

    // Reset the global primal vector
    M_primalField->cleanField();

    // Reset the global dual vector
    M_dualField->cleanField();

    // Reset the global hybrid vector
    M_hybridField->cleanField();

} // resetVariables

// Apply the boundary conditions to the global hybrid matrix and the global hybrid right hand side.
template < typename MeshType >
void
DarcySolverLinear < MeshType >::
applyBoundaryConditions ()
{

    // Check if the boundary conditions are set or not.
    ASSERT ( M_boundaryConditionHandler.get(), "DarcySolverLinear : boundary conditions not set." );

    // Check if the boundary conditions were updated.
    if ( !M_boundaryConditionHandler->bcUpdateDone() )
    {
        // Update the boundary conditions handler. We use the finite element of the boundary of the dual variable.
        M_boundaryConditionHandler->bcUpdate ( *M_dualField->getFESpace().mesh(),
                                               M_dualField->getFESpace().feBd(),
                                               M_dualField->getFESpace().dof() );
    }

    // Ignoring non-local entries, otherwise they are summed up lately.
    vector_Type rhsFull ( *M_rhs, Unique );

    /* Update the global hybrid matrix and the global hybrid right hand side with the boundary conditions.
       It takes care of the current time for DarcySolverTransient derived class. */
    bcManage ( *M_matrHybrid, rhsFull, *M_dualField->getFESpace().mesh(),
               M_dualField->getFESpace().dof(), *M_boundaryConditionHandler, M_dualField->getFESpace().feBd(), 1.,
               M_data->dataTimePtr()->time() );

    // Save the global hybrid right hand side.
    *M_rhs = rhsFull;

} // applyBoundaryConditions

// Reorder a non full stored symmetric matrix.
template < typename MeshType >
template < typename MatrixType >
void
DarcySolverLinear < MeshType >::
symmetrizeMatrix ( Int N, MatrixType& A  )
{

    for ( UInt i ( 0 ); i < static_cast<UInt> (N); ++i )
    {
        for ( UInt j ( i + 1 ); j < static_cast<UInt> (N); ++j )
        {
            A ( i, j ) = A ( j, i );
        }
    }

} //symmetrizeMatrix


} // namespace LifeV


#endif //_DARCYSOLVERLINEAR_HPP_

// -*- mode: c++ -*-
