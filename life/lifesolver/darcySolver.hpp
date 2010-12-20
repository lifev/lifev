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

#ifndef _DARCYSOLVER_H_
#define _DARCYSOLVER_H_ 1

#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/bcManage.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/clapack.hpp>
#include <life/lifealg/cblas.hpp>
//
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataDarcy.hpp>
//
#include <life/lifecore/displayer.hpp>


// Local namespace to store the default source term and the default permeability tensor.
namespace
{

/*! @class Default Darcy source term. Needed in DarcySolver class.
  This class implement the default source term for the Darcy problem, it is the zero function.
*/
struct DarcyDefaultSource
{
    LifeV::Real operator()( const LifeV::Real&, const LifeV::Real&,
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
//! class DarcySolver   This class implements a mixed-hybrid FE Darcy solver.

/*!
  @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

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
  \Lambda^{-1} \sigma + \nabla p = 0 & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  \nabla \cdot \sigma - f = 0        & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  p = g_D                            & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
  \sigma \cdot n + h p = g_R         & \mathrm{on} & \Gamma_R\,.
  \end{array}
  \right.
  \f]
  The first equation is the Darcy equation and the second equation is the conservation law. The data in the system
  are \f$ \Lambda \f$ that is the permeability tensor, \f$ f \f$ that is the source term, \f$ \Gamma_D \f$ that
  is the subset of the boundary of \f$ \Omega \f$ with Dirichlet boundary conditions with datum \f$ g_D \f$ and
  \f$ \Gamma_R \f$ that is the part of the boundary of \f$ \Omega \f$ with Robin, or Neumann, boundary conditions
  with data \f$ h \f$ and \f$ g_R \f$. We suppose that \f$ \partial \Omega = \Gamma_D \cup \Gamma_R \f$ and
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
  V = L^2 (\Omega ) \,,\vspace{0.2cm} \					\
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
  h(\lambda, \mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} h \mu \lambda \,,\vspace{0.2cm}\\
  F(v) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_K f v\,,&
  G(\mu) = \displaystyle \sum_{K \in \mathcal{T}_h} \int_{\partial K \cap \Gamma_R} g \mu\,,
  \end{array}
  \f]
  we obtain the Darcy problem in the weak form: find \f$ (\sigma, \, p, \, \lambda) \in Z \times V \times \Lambda \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma, \tau) + b(p, \tau) + c(\lambda, \tau) = 0\,,  & \forall \tau \in Z \,,\vspace{0.2cm}\\
  b(v, \sigma) = - F(v)\,,                                  & \forall v \in V \,,\vspace{0.2cm}\\
  c(\mu, \sigma) + h(\lambda, \mu) = G(\mu) \,,           & \forall \mu \in \Lambda\,.
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
  <BR>
  The finite dimensional problem is: find \f$ (\sigma_h,\, p_h, \, \lambda_h) \in Z_h \times V_h \times \Lambda_h \f$ such that
  \f[
  \left\{
  \begin{array}{l l}
  a(\sigma_h, \tau_h) + b(p_h, \tau_h) + c(\lambda_h, \tau_h) = 0\,,  & \forall \tau_h \in Z_h \,,\vspace{0.2cm}\\
  b(v_h, \sigma_h) = -F(v_h)\,,                                  & \forall v_h \in V_h \,,\vspace{0.2cm}\\
  c(\mu_h, \sigma_h) + h(\lambda_h, \mu_h) = G(\mu_h) \,,           & \forall \mu_h \in \Lambda_h\,.
  \end{array}
  \right.
  \f]
  To solve the problem we use the static condensation procedure, i.e. the unknowns in the discrete
  weak system are not independent and \f$ p_K \f$, \f$\sigma_K \f$ may be written in function of
  \f$ \lambda_K \f$ alone. We introduce the following local matrices
  \f[
  \begin{array}{l l l}
  \left[ A \right]_{ij} = \displaystyle   \int_K \Lambda^{-1} \psi_j \cdot \psi_i \,, &
  \left[ B \right]_{ij} = \displaystyle - \int_K \phi_j \nabla \cdot \psi_i \,, &
  \left[ C \right]_{ij} = \displaystyle   \int_{\partial K} \xi_i \psi_j \cdot n \,, \vspace{0.2cm} \\
  \left[ H \right]_{ij} = \displaystyle   \int_{\partial K \cap \Gamma_R} h \xi_i \xi_j\,, &
  \left[ F \right]_{j}  = \displaystyle   \int_K f \phi_j\,, &
  \left[ G \right]_{j}  = \displaystyle   \int_{\partial K \cap \Gamma_R } g \xi_j\,,
  \end{array}
  \f]
  where we avoid to write the dependence on the triangle \f$ K \f$ in all the matrices and vectors. <BR>
  The local matrix formulation of the finite dimensional problem is
  \f[
  \left\{
  \begin{array}{l}
  A \sigma_K + B p_K + C \lambda_K = 0\,, \vspace{0.2cm} \\
  B^T \sigma_K = -F \,,                    \vspace{0.2cm}\\
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
  B^T & 0 &  0 \vspace{0.2cm} \\
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
  L_K = - C^T A^{-1} C + C^T A^{-1} B \left( B^T A^{-1} B \right)^{-1} B^T A^{-1} C + H \,, \vspace{0.2cm} \\
  r_K = G + C^T A^{-1} B \left( B^T A^{-1} B \right)^{-1} F\,,
  \end{array}
  \f]
  Imposing that at each edge or face the hybrid unknown is single value we obtain a linear system for the hybrid unknown
  \f[
  L \lambda = r \,.
  \f]
  We recover the primal and dual variable as a post-process from the hybrid variable at element level, so we have
  \f[
  \begin{array}{l}
  p_K = \left( B^T A^{-1} B \right)^{-1} \left( F - B^T A^{-1} C \lambda_K \right)\,, \vspace{0.2cm} \\
  \sigma_K = -A^{-1} \left( B p_K + C \lambda_K \right) \,.
  \end{array}
  \f]
  @note In the code we do not use the matrix \f$ H \f$ and the vector \f$ G \f$, because all the boundary
  conditions are imposed via BCHandler class.
  @todo Insert any scientific publications that use this solver.
  @bug If the save flag for the exporter is setted to 0 the program fails.
*/

template< typename Mesh, typename SolverType = LifeV::SolverTrilinos >
class DarcySolver
{

public:

    //! @name Public Types
    //@{

    typedef SolverType                             solver_Type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&,
                                   const Real&, const UInt& )> Function;

    typedef inversePermeability<Mesh, SolverType>  permeability_Type;
    typedef boost::shared_ptr<permeability_Type>   permeabilityPtr_Type;

    typedef DataDarcy<Mesh>                        data_Type;

    typedef Mesh                                   mesh_Type;

    typedef BCHandler                              bchandler_raw_Type;
    typedef boost::shared_ptr<bchandler_raw_Type>  bchandler_Type;

    typedef typename solver_Type::matrix_type      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;

    typedef typename solver_Type::vector_type      vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef typename solver_Type::prec_raw_type    prec_raw_Type;
    typedef typename solver_Type::prec_type        prec_Type;

    typedef Epetra_Comm                            comm_Type;
    typedef boost::shared_ptr< comm_Type >         commPtr_Type;

    //@}

    //! @name Constructors and destructor
    //@{

    //! Full constructor for the class.

    /*!
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param bcHandler Boundary conditions for the problem.
      @param comm Shared pointer of the Epetra communicator.
    */
    DarcySolver ( const data_Type&                dataFile,
                        FESpace<Mesh, EpetraMap>& primal_FESpace,
                        FESpace<Mesh, EpetraMap>& dual_FESpace,
                        FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                        FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                  const bchandler_raw_Type&       bcHandler,
                  const commPtr_Type&             comm );

    //! Constructor for the class without the definition of the boundary handler.

    /*!
      @param dataFile Data for the problem.
      @param primal_FESpace Primal finite element space.
      @param dual_FESpace Dual finite element space.
      @param hybrid_FESpace Hybrid finite element space.
      @param VdotN_FESpace Dual basis function dot outward unit normal at each face (3D) or edge (2D) finite element space.
      @param comm Shared pointer of the Epetra communicator.
    */
    DarcySolver ( const data_Type&                dataFile,
                        FESpace<Mesh, EpetraMap>& primal_FESpace,
                        FESpace<Mesh, EpetraMap>& dual_FESpace,
                        FESpace<Mesh, EpetraMap>& hybrid_FESpace,
                        FESpace<Mesh, EpetraMap>& VdotN_FESpace,
                  const commPtr_Type&             comm );

    //! Virtual destructor.
    virtual ~DarcySolver ();

    //@}

    //! @name Methods
    //@{

    //! Build the global hybrid system, the right hand and apply the boundary conditions.
    void buildSystem ();

    //! Set up the linear solver and the preconditioner for the linear system.
    virtual void setup ();

    //! Solve the global hybrid system.
    void solve ();

    //! Compute primal and dual variables from the hybrid variable as a post process.
    void computePrimalAndDual ();

    //@}

    //! @name Set Methods
    //@{

    //! Set the boundary conditions.
    /*!
      @param bcHandler Boundary condition handler for the problem.
    */
    void setBC ( bchandler_raw_Type& bcHandler )
    {
        M_BCh   = &bcHandler;
        M_setBC = true;
    }

    //! Set source term
    /*!
      Set the source term, the default source term is the zero function.
      @param source Source term for the problem.
    */
    void setSourceTerm ( const Function& source )
    {
        M_source = source;
    }

    //! Set inverse diffusion tensor
    /*!
      Set the inverse of diffusion tensor, the default inverse of permeability is the identity matrix.
      @param invPerm Inverse of the permeability tensor for the problem.
    */
    void setInversePermeability ( const permeability_Type& invPerm )
    {
        M_inversePermeability.reset ( new permeability_Type( invPerm ) );
    }

    //! Set the hybrid solution vector.
    /*!
      @param hybrid Constant vectorPtr_Type reference of the hybrid vector.
     */
    void setHybridSolution ( const vectorPtr_Type& hybrid  )
    {
        M_hybrid = hybrid;
    }

    //! Set the primal solution vector.
    /*!
      @param primal Constant vector_Type reference of the primal solution.
    */
    void setPrimalSolution ( const vectorPtr_Type& primal )
    {
        M_primal = primal;
    }

    //! Set the dual solution vector.
    /*!
      @param dual Constant vectorPtr_Type reference of the dual solution.
    */
    void setDualSolution ( const vectorPtr_Type& dual )
    {
        M_dual = dual;
    }

    //@}

    //! @name Get Methods
    //@{

    //! Returns pointer to the hybrid solution vector.
    /*!
      @return Constant vectorPtr_Type reference of the hybrid vector.
     */
    const vectorPtr_Type& hybridSolution () const
    {
        return M_hybrid;
    }

          vectorPtr_Type& hybridSolution ()
    {
        return M_hybrid;
    }

    //! Returns pointer to the primal solution vector.
    /*!
      @return Constant vector_Type reference of the primal solution.
    */
    const vectorPtr_Type& primalSolution () const
    {
        return M_primal;
    }

          vectorPtr_Type& primalSolution ()
    {
        return M_primal;
    }

    //! Returns pointer to the dual solution vector.
    /*!
      @return Constant vectorPtr_Type reference of the dual solution.
    */
    const vectorPtr_Type& dualSolution () const
    {
        return M_dual;
    }

          vectorPtr_Type& dualSolution ()
    {
        return M_dual;
    }

    /*!
      Returns the pointer of the residual vector.
      @return Constant vectorPtr_Type reference of the residual.
    */
    const vectorPtr_Type& residual () const
    {
        return M_residual;
    }

          vectorPtr_Type& residual ()
    {
        return M_residual;
    }

    //! Returns whether bounday conditions is set
    /*!
      @return Constant boolean with value true if the boundary condition is setted,
      false otherwise
    */
    bool isBCset ()
    {
        return M_setBC;
    }

    //! Returns boundary conditions handler.
    /*!
      @return Reference of boundary conditions handler.
    */
    const bchandler_Type& bcHandler () const
    {
        return M_BCh;
    }

          bchandler_Type& bcHandler ()
    {
        return M_BCh;
    }

    //! Returns Epetra local map.
    /*!
      @return Constant EpetraMap reference of the problem.
    */
    const EpetraMap& getMap () const
    {
        return M_localMap;
    }

          EpetraMap& getMap ()
    {
        return M_localMap;
    }

    //! Returns Epetra communicator.
    /*!
      @return Constant reference of the shared pointer of the communicator of the problem.
    */
    const commPtr_Type & getComm () const
    {
        return M_displayer.comm();
    }

          commPtr_Type & getComm ()
    {
        return M_displayer.comm();
    }

    //! Returns displayer.
    /*!
      @return Constant reference of the displayer of the problem.
    */
    const Displayer & getDisplayer() const
    {
        return M_displayer;
    }

          Displayer & getDisplayer()
    {
        return M_displayer;
    }

    //@}

protected:

    //! @name Protected Methods
    //@{

    //! Pre-computes local (element independant) matrices
    /*!
      Compute all the local matrices that are independant
      from the geometrical element.
    */
    virtual void computeConstantMatrices ();

    //! Computes local (element dependent) matrices
    /*!
      Locally update the current finite element for the primal
      and dual finite element space, then compute the Hdiv mass
      matrix.
      @param iElem Id of the current geometrical element.
    */
    virtual void localElementComputation ( const UInt & iElem );

    //! Performs static condensation
    /*!
      Locally eliminate pressure and velocity DOFs, create the local
      hybrid matrix and local hybrid right hand side.
    */
    virtual void staticCondensation ();

    //! Compute locally, as a post process, the primal and dual variable.
    virtual void localComputePrimalAndDual ();

    //! Update all problem variables
    /*!
      Update all the variables of the problem before the construction of
      the global hybrid matrix, e.g. reset the global hybrid matrix.
      It is principally used for a time dependent derived class, in fact we
      want to clear each time step all the variables.
    */
    virtual void updateVariables ();

    //! Apply the boundary condition
    /* Apply BC to the hybrid global matrix and to the hybrid global right hand side.
    */
    void applyBoundaryConditions ();

    //! Make matrix symmetric
    /*!
      Transform a symmetric matrix that is stored only in the lower or upper
      triangular part in a full symmetric matrix.
      @param UPLO Character that indicate if A is in the upper triangular part
      using flag 'U' or in the lower triangular part using flag 'L'.
      @param N The size of the matrix A.
      @param A The matrix to be reordered.
    */
    template<typename matrix>
    void symmetrizeMatrix ( char* UPLO, int* N, matrix& A  );

    //@}

    //! @name Parallel stuff
    //@{

    //! MPI process identifier.
    UInt      M_me;

    //! Local map.
    EpetraMap M_localMap;

    //! Displayer.
    Displayer M_displayer;

    //@}

    // Data of the problem.
    //! @name Data of the problem
    //@{

    //! Data for Darcy solvers.
    const data_Type&     M_data;

    //! Source function.
    Function             M_source;

    //! Permeability tensor.
    permeabilityPtr_Type M_inversePermeability;

    //! Bondary conditions handler.
    bchandler_raw_Type*  M_BCh;

    //! Flag if the boundary conditions are setted or not.
    bool                 M_setBC;

    //@}

    // Finite element spaces.
    //! @name Finite element spaces
    //@{

    //! Primal finite element space.
    FESpace<Mesh, EpetraMap>&  M_primal_FESpace;

    //! Dual finite element space.
    FESpace<Mesh, EpetraMap>&  M_dual_FESpace;

    //! Hybrid finite element space.
    FESpace<Mesh, EpetraMap>&  M_hybrid_FESpace;

    //! Dual function dot outward unit vector finite element space.
    FESpace<Mesh, EpetraMap>&  M_VdotN_FESpace;

    //@}

    // Algebraic stuff.
    //! @name Algebraic stuff
    //@{

    //! Hybrid matrix.
    matrixPtr_Type                          M_matrHybrid;

    //! Right hand side.
    vectorPtr_Type                          M_rhs;

    //! Primal solution.
    vectorPtr_Type                          M_primal;

    //! Dual solution.
    vectorPtr_Type                          M_dual;

    //! Hybrid solution.
    vectorPtr_Type                          M_hybrid;

    //! Residual.
    vectorPtr_Type                          M_residual;

    //! Linear solver.
    solver_Type                             M_linearSolver;

    //! Epetra preconditioner for the linear system.
    boost::shared_ptr<EpetraPreconditioner> M_prec;

    //@}

    // Elementary matrices and vectors used for the static condensation.
    //! @name Elementary matrices and vectors used for the static condensation.
    //@{

    //! Element vector for the right hand side, lives in RTkHyb(fe).
    ElemVec M_elvecHyb;

    //! Element vector for the source term of the problem, lives in Qk(fe).
    ElemVec M_elvecSource;

    //! Element vector for the dual variabl, lives in RTk(fe).
    ElemVec M_elvecFlux;

    /*! Mixed hybrid matrix, maps
      \f[
      [A\,| B\,| C] \qquad \mathrm{in} \qquad
      \left[
      \begin{array}{c c c}
      A   & B & C \\
      B^T & 0 & 0 \\
      C^T & 0 & 0 \\
      \end{array}
      \right]
      \f]
    */
    ElemMat M_elmatMix;

    //! Element hybrid matrix.
    ElemMat M_elmatHyb;

    //! Temporary array of size the square of the number of degrees of freedom of the primal variable.
    KNM<Real> M_BtB;

    /*! Temporary array of size the square of the number of degrees of freedom of the hybrid variable.
        It stores of the final hybrid array. */
    KNM<Real> M_CtC;

    /*! Temporary array of size the prodocut between the number of degrees of freedom of the primal
        and the degrees of freedom of the hybrid variable. */
    KNM<Real> M_BtC;

    //@}

}; // class DarcySolver

//
// IMPLEMENTATION
//

// ====================================================================================================
// Constructors and destructors
// ====================================================================================================

// Complete constructor.
template<typename Mesh, typename SolverType>
DarcySolver<Mesh, SolverType>::
DarcySolver ( const data_Type&                 dataFile,
              FESpace<Mesh, EpetraMap>&  primal_FESpace,
              FESpace<Mesh, EpetraMap>&  dual_FESpace,
              FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
              FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
              const bchandler_raw_Type&        bcHandler,
              const commPtr_Type&              comm ):
        // Parallel stuff.
        M_me                     ( comm->MyPID() ),
        M_localMap               ( hybrid_FESpace.map() ),
        M_displayer              ( comm ),
        // Data of the problem.
        M_data                   ( dataFile ),
        M_source                 ( DarcyDefaultSource() ),
        M_BCh                    ( &bcHandler ),
        M_setBC                  ( true ),
        // Finite element spaces.
        M_primal_FESpace         ( primal_FESpace ),
        M_dual_FESpace           ( dual_FESpace ),
        M_hybrid_FESpace         ( hybrid_FESpace ),
        M_VdotN_FESpace          ( VdotN_FESpace ),
        // Algebraic stuff.
        M_matrHybrid             ( new matrix_Type ( M_localMap ) ),
        M_rhs                    ( new vector_Type ( M_localMap ) ),
        M_primal    		 ( new vector_Type ( M_primal_FESpace.map() ) ),
        M_dual			 ( new vector_Type ( M_dual_FESpace.map(), Repeated ) ),
        M_hybrid                 ( new vector_Type ( M_hybrid_FESpace.map() ) ),
        M_residual               ( new vector_Type ( M_localMap ) ),
        M_linearSolver           ( ),
        M_prec                   ( ),
        // Elementary matrices and vectors used for the static condensation.
        M_elvecHyb               ( M_hybrid_FESpace.refFE().nbDof(), 1 ),
        M_elvecSource            ( M_primal_FESpace.refFE().nbDof(), 1 ),
        M_elvecFlux              ( M_dual_FESpace.refFE().nbDof(), 1 ),
        M_elmatMix               ( M_dual_FESpace.refFE().nbDof(), 1, 1, M_primal_FESpace.refFE().nbDof(), 0, 1,
                                   M_hybrid_FESpace.refFE().nbDof(), 0, 1 ),
        M_elmatHyb               ( M_hybrid_FESpace.refFE().nbDof(), 1, 1 ),
        M_BtB                    ( M_primal_FESpace.refFE().nbDof(), M_primal_FESpace.refFE().nbDof() ),
        M_CtC                    ( M_hybrid_FESpace.refFE().nbDof(), M_hybrid_FESpace.refFE().nbDof() ),
        M_BtC                    ( M_primal_FESpace.refFE().nbDof(), M_hybrid_FESpace.refFE().nbDof() )
{

} // Constructor


// Constructor without boundary condition handler.
template<typename Mesh, typename SolverType>
DarcySolver<Mesh, SolverType>::
DarcySolver ( const data_Type&                 dataFile,
              FESpace<Mesh, EpetraMap>&  primal_FESpace,
              FESpace<Mesh, EpetraMap>&  dual_FESpace,
              FESpace<Mesh, EpetraMap>&  hybrid_FESpace,
              FESpace<Mesh, EpetraMap>&  VdotN_FESpace,
              const commPtr_Type&              comm ):
        // Parallel stuff.
        M_me                     ( comm->MyPID() ),
        M_localMap               ( hybrid_FESpace.map() ),
        M_displayer              ( comm ),
        // Data of the problem.
        M_data                   ( dataFile ),
        M_source                 ( DarcyDefaultSource() ),
        M_setBC                  ( false ),
        // Finite element spaces.
        M_primal_FESpace         ( primal_FESpace ),
        M_dual_FESpace           ( dual_FESpace ),
        M_hybrid_FESpace         ( hybrid_FESpace ),
        M_VdotN_FESpace          ( VdotN_FESpace ),
        // Algebraic stuff.
        M_matrHybrid             ( new matrix_Type ( M_localMap ) ),
        M_rhs                    ( new vector_Type ( M_localMap ) ),
        M_primal    	    	 ( new vector_Type ( M_primal_FESpace.map() ) ),
        M_dual		        	 ( new vector_Type ( M_dual_FESpace.map(), Repeated ) ),
        M_hybrid                 ( new vector_Type ( M_hybrid_FESpace.map() ) ),
        M_residual               ( new vector_Type ( M_localMap ) ),
        M_linearSolver           ( ),
        M_prec                   ( ),
        // Local matrices and vectors.
        M_elvecHyb               ( M_hybrid_FESpace.refFE().nbDof(), 1 ),
        M_elvecSource            ( M_primal_FESpace.refFE().nbDof(), 1 ),
        M_elvecFlux              ( M_dual_FESpace.refFE().nbDof(), 1 ),
        M_elmatMix               ( M_dual_FESpace.refFE().nbDof(), 1, 1, M_primal_FESpace.refFE().nbDof(), 0, 1,
                                   M_hybrid_FESpace.refFE().nbDof(), 0, 1 ),
        M_elmatHyb               ( M_hybrid_FESpace.refFE().nbDof(), 1, 1 ),
        M_BtB                    ( M_primal_FESpace.refFE().nbDof(), M_primal_FESpace.refFE().nbDof() ),
        M_CtC                    ( M_hybrid_FESpace.refFE().nbDof(), M_hybrid_FESpace.refFE().nbDof() ),
        M_BtC                    ( M_primal_FESpace.refFE().nbDof(), M_hybrid_FESpace.refFE().nbDof() )
{

} // Constructor

// Virtual destructor.
template<typename Mesh, typename SolverType>
DarcySolver<Mesh, SolverType>::
~DarcySolver ( void )
{

} // Destructor

// ===========================================================================================
// Public methods
// ==========================================================================================

// Build the global hybrid matrix and global hybrid right hand side, with boundary conditions.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
buildSystem ()
{

    // Chronos.
    Chrono chronoStaticCondensation;
    Chrono chronoConstantLocalMatrix;
    Chrono chronoGlobalAssembleMatrix;
    Chrono chronoGlobalAssembleVector;

    // The total number of elements in the mesh.
    UInt meshNumberOfElements = M_primal_FESpace.mesh()->numElements();

    M_displayer.leaderPrint("  Darcy solver - Computing constant matrices   ...          ");

    // Clear all the local matrices and vectors.
    M_elvecHyb.zero();
    M_elvecSource.zero();
    M_elvecFlux.zero();
    M_elmatMix.zero();
    M_elmatHyb.zero();

    chronoConstantLocalMatrix.start();

    // Compute all the constant matrices, e.g. the matrix B and C
    computeConstantMatrices();

    chronoConstantLocalMatrix.stop();

    M_displayer.leaderPrintMax( "done in " , chronoConstantLocalMatrix.diff() );

    // If setted print the constant matrices computed.
    if ( M_displayer.isLeader() &&  M_data.verbose() > static_cast<UInt>(1) )
    {
        M_displayer.leaderPrint( "elmatHyb :\n\n" );
        M_elmatHyb.showMe();

        M_displayer.leaderPrint( "elmatMix :\n\n" );
        M_elmatMix.showMe();
    }

    // Update all the variables, e.g. reset the hybrid matrix
    updateVariables();

    M_displayer.leaderPrint("                 Perform Static Condensation   ...          ");
    chronoStaticCondensation.start();

    //---------------------------------------
    //    LOOP ON ALL THE VOLUME ELEMENTS
    //---------------------------------------

    for ( UInt iElem(1); iElem <= meshNumberOfElements; ++iElem )
    {
        // Compute the Hdiv mass matrix as a local matrix depending on the current element.
        localElementComputation( iElem );

        staticCondensation();

        /* Assemble the global hybrid matrix.
           M_primal_FESpace is used instead of M_hybrid_FESpace for currentLocalId,
           because currentFE cannot store a RefFEHybrid. */
        assembleMatrix( *M_matrHybrid,
                        M_primal_FESpace.fe().currentLocalId(),
                        M_elmatHyb,
                        M_hybrid_FESpace.refFE().nbDof(),
                        M_hybrid_FESpace.dof(),
                        0,0,0,0);

        /* Assemble the global hybrid right hand side.
           M_primal_FESpace is used instead of M_hybrid_FESpace for currentLocalId,
           because currentFE cannot store a RefFEHybrid. */
        assembleVector( *M_rhs,
                        M_primal_FESpace.fe().currentLocalId(),
                        M_elvecHyb,
                        M_hybrid_FESpace.refFE().nbDof(),
                        M_hybrid_FESpace.dof(), 0 );
    }

    //-----------------------------------
    // END OF LOOP ON THE VOLUME ELEMENTS
    //-----------------------------------

    chronoStaticCondensation.stop();
    M_displayer.leaderPrintMax( "done in " , chronoStaticCondensation.diff() );

    // Apply boundary conditions to the hybrid global matrix and to the hybrid global right hand side.
    applyBoundaryConditions();

    M_displayer.leaderPrint("                 Assemble Global Hybrid Matrix ...          ");

    chronoGlobalAssembleMatrix.start();

    // Assemble the global hybrid matrix.
    M_matrHybrid->globalAssemble();

    chronoGlobalAssembleMatrix.stop();

    M_displayer.leaderPrintMax( "done in " , chronoGlobalAssembleMatrix.diff() );

    M_displayer.leaderPrint("                 Assemble Global Hybrid Vector ...          ");

    chronoGlobalAssembleVector.start();

    // Assemble the global hybrid vector.
    M_rhs->globalAssemble();

    chronoGlobalAssembleVector.stop();

    M_displayer.leaderPrintMax( "done in " , chronoGlobalAssembleVector.diff() );

} // buildSystem

// Set up the linear solver and the preconditioner.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
setup ()
{

    GetPot dataFile( *(M_data.dataFile()) );

    // Set up data for the linear solver and the preconditioner.
    M_linearSolver.setDataFromGetPot( dataFile, "darcy/solver" );
    M_linearSolver.setUpPrec( dataFile, "darcy/prec" );
    M_linearSolver.setCommunicator( M_displayer.comm() );

    // Choose the preconditioner type.
    std::string precType = dataFile( "darcy/prec/prectype", "Ifpack");

    // Create a preconditioner object.
    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT( M_prec.get() != 0, "DarcySolver : Preconditioner not set" );

} // setup

// Solve the linear system.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
solve ()
{

    // Set the matrix.
    M_linearSolver.setMatrix( *M_matrHybrid );

    // Solve the linear system, if used a iterative scheme numIter stores the number of iterations.
    UInt numIter = M_linearSolver.solveSystem( *M_rhs, *M_hybrid, M_matrHybrid );

    M_displayer.leaderPrint( "                 Number of iterations  ", numIter, "\n" );

} // solve

// Compute primal and dual unknown as a post process.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
computePrimalAndDual ()
{

    // Chrono.
    Chrono chronoComputePrimalAndDual;

    // The total number of elements in the mesh.
    UInt meshNumberOfElements = M_primal_FESpace.mesh()->numElements();

    M_displayer.leaderPrint("                 Compute pressure and flux     ...          ");
    chronoComputePrimalAndDual.start();

    /*==========================================
      POST PROCESSING
      compute the pressure (Qk or Pk / element)
      and the velocity (RTk / element) => 2 (opposite) velocities / face
      ==========================================*/

    /* The vector H_hybrid is an unique epetra vector, so it does not share one layer elements.
       The reconstruction of the primal and dual variabile need this share, so we create a repeated
       epetra vector. This operation maximize the efficiency of send/receive time cost instead to
       send and receive each time a single datum. */
    vector_Type M_hybrid_Repeated( *M_hybrid, Repeated );

    //---------------------------------------
    //    LOOP ON ALL THE VOLUME ELEMENTS
    //---------------------------------------

    for ( UInt iElem(1); iElem <= meshNumberOfElements; ++iElem )
    {
        // Clear the hybrid right hand side, it will store the extranctions from the global vecotor.
        M_elvecHyb.zero();

        // Compute the matrix A as a local matrix depending on the current element.
        localElementComputation( iElem );

        // Extract the computed hybrid variable for the current finite element and put it into elvecHyb.
        extract_vec( M_hybrid_Repeated,
                     M_elvecHyb,
                     M_hybrid_FESpace.refFE(),
                     M_hybrid_FESpace.dof(),
                     M_primal_FESpace.fe().currentLocalId(), 0 );

        // Compute locally the primal and dual variable.
        localComputePrimalAndDual();

        /* Put the primal variable of the current finite element, stored in elvecSource,
           in the global vector M_primal. */
        assembleVector( *M_primal,
                        M_primal_FESpace.fe().currentLocalId(),
                        M_elvecSource,
                        M_primal_FESpace.refFE().nbDof(),
                        M_primal_FESpace.dof(), 0 );


        for ( UInt iLocalFace(1); iLocalFace <=  M_dual_FESpace.mesh()->element( iElem ).S_numLocalFaces; ++iLocalFace )
        {
            UInt iGlobalFace( M_dual_FESpace.mesh()->localFaceId( iElem, iLocalFace ) );
            if ( M_dual_FESpace.mesh()->faceElement( iGlobalFace, 1 ) != iElem )
            {
                M_elvecFlux[ iLocalFace - 1 ] = 0;
            }
        }

        /* Put the dual variable of the current finite element, stored in elvecFlux,
           in the global vector M_dual. */
        assembleVector( *M_dual,
                        M_dual_FESpace.fe().currentLocalId(),
                        M_elvecFlux,
                        M_dual_FESpace.refFE().nbDof(),
                        M_dual_FESpace.dof(), 0);

    }


    //---------------------------------------
    // END OF THE LOOP ON THE VOLUME ELEMENTS
    //---------------------------------------

    // Assemble the primal variable.
    M_primal->globalAssemble();

    // It is wrong to assemble the dual variable, no communications required.

    chronoComputePrimalAndDual.stop();
    M_displayer.leaderPrintMax( "done in " , chronoComputePrimalAndDual.diff() );

} // computePrimalAndDual

// ==============================================================================
// Protected methods
// ==============================================================================

// Compute all the constant matrices.
template <typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
computeConstantMatrices ()
{

    /* Update the divergence matrix, it is independant of the current element
       thanks to the Piola transform. */
    grad_Hdiv( static_cast<Real>(1.),
               M_elmatMix,
               M_dual_FESpace.fe(),
               M_primal_FESpace.fe(), 0, 1 );

    /* Update the boundary matrix, it is independant of the current element
       thanks to the Piola transform.
       Here we use the fact that a RefHybridFE "IS A" RefFE and the method refFE of
       a FESpace object. In fact the method refFE return a const RefFE&, but the function
       TP_VdotN_Hdiv takes two const RefHybridFE& so we must cast a const RefFE&
       to a const RefHybridFE&. The cast of type is static and uses pointers. */
    TP_VdotN_Hdiv( static_cast<Real>(1.),
                   M_elmatMix,
                   *static_cast<const RefFEHybrid*>(&M_hybrid_FESpace.refFE()),
                   *static_cast<const RefFEHybrid*>(&M_VdotN_FESpace.refFE()),
                   0, 2 );

} // computeConstantMatrices

// Update the primal and dual variable at the current element and compute the element Hdiv mass matrix.
template <typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
localElementComputation ( const UInt & iElem )
{

    // Update the current element of ID iElem only for the dual variable.
    M_dual_FESpace.fe().update( M_primal_FESpace.mesh()->element( iElem ),
                                UPDATE_PHI_VECT | UPDATE_WDET );

    /* Update the current element of ID iElem only for the primal variable,
       it is used for computing the source term. */
    M_primal_FESpace.fe().update( M_primal_FESpace.mesh()->element( iElem ),
                                  UPDATE_QUAD_NODES | UPDATE_WDET );

    /* Modify the (0,0) block (A) of the matrix M_elmatMix. The blocks (0,1) (B)
       and (0,2) (C) are independent of the element and have already been computed. */

    // Only one block is set to zero since the other ones are not recomputed.
    M_elmatMix.block(0, 0) = static_cast<Real>( 0 );

    // Get the coordinate of the barycenter of the current element of ID iElem.
    Real xg(0), yg(0), zg(0);
    M_primal_FESpace.fe().barycenter( xg, yg, zg );

    /* Compute the Hdiv mass matrix. We pass the time at the inverse of the permeability
       because the DarcySolverTransient needs the pemeability time dependent. In this case
       we do not have a time evolution. */
    mass_Hdiv( (*M_inversePermeability)( M_data.dataTime()->time(), xg, yg, zg, iElem ),
               M_elmatMix,
               M_dual_FESpace.fe(), 0, 0 );

} // localElementComputation

// Perform the static condensation for the local hybrid matrix.
template <typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
staticCondensation ()
{

    // Flags for the BLAS and LAPACK routine.

    int INFO[1]         = {0};
    // Number of columns of the right hand side := 1.
    int NBRHS[1]        = {1};
    // Primal variable degrees of freedom.
    int NBP[1]          = { M_primal_FESpace.refFE().nbDof() };
    // Dual variable degrees of freedom.
    int NBU[1]          = { M_dual_FESpace.refFE().nbDof() };
    // Hybrid variable degree of freedom.
    int NBL[1]          = { M_hybrid_FESpace.refFE().nbDof() };

    double ONE[1]      = {1.0};
    double MINUSONE[1] = {-1.0};
    double ZERO[1]     = {0.0};

    // Parameter that indicate the Lower storage of matrices.
    char UPLO[1]     = {'L'};

    // Paramater that indicate the Transpose of matrices.
    char   TRANS[1]    = {'T'};
    char NOTRANS[1]    = {'N'};

    // Parameter that indicates whether the matrix has diagonal unit ('N' means no)
    char NODIAG[1] = {'N'};

    // Create and assign the local matrices A, B and C.
    ElemMat::matrix_type A = M_elmatMix.block( 0, 0 );
    ElemMat::matrix_type B = M_elmatMix.block( 0, 1 );
    ElemMat::matrix_type C = M_elmatMix.block( 0, 2 );

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
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/BLAS/SRC/dsyrk.f */
    dsyrk_ (UPLO, TRANS, NBP, NBU, ONE, B, NBU, ZERO, M_BtB, NBP);

    /* Put in M_CtC the matrix C^T * L^{-T} * L^{-1} * C = C^T * A^{-1} * C
       M_CtC stored only on lower part.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f  */
    dsyrk_ (UPLO, TRANS, NBL, NBU, ONE, C, NBU, ZERO, M_CtC, NBL);

    /* Put in M_BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       M_BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBP, NBL, NBU, ONE, B, NBU, C, NBU, ZERO, M_BtC, NBP);

    /* Put in M_BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f  */
    dpotrf_ (UPLO, NBP, M_BtB, NBP, INFO);
    ASSERT_PRE( !INFO[0],"Lapack factorization of BtB is not achieved." );

    /* Put in M_BtC the matrix LB^{-1} * M_BtC = LB^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBL, M_BtB, NBP, M_BtC, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    /* Put in M_CtC the matrix -M_CtC + M_BtC^T * M_BtC
       Result stored only on lower part, the matrix M_CtC stores
       M_CtC = -C^T * A^{-1} * C + C^T * A^{-t} * B * ( B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/slatec/lin/dsyrk.f  */
    dsyrk_ (UPLO, TRANS, NBL, NBP, ONE, M_BtC, NBP, MINUSONE, M_CtC, NBL);

    //...................................
    //      END OF MATRIX OPERATIONS
    //...................................

    /* Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       B^T stores LB and LB^T where LB and LB^T is the factorization of B^T * A^{-1} * B
       M_BtC stores LB^{-1} * B^T * A^{-1} * C
       M_CtC stores -C^T * A^{-1} * C + C^T * A^{-t} * B * (B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C */

    //..........................
    //     VECTOR OPERATIONS
    //..........................

    // Clear some vectors.
    M_elvecSource.zero();
    M_elvecHyb.zero();

    // Compute the right hand side.
    source( M_source,
            M_elvecSource,
            M_primal_FESpace.fe(),
            static_cast<Real>(0), 0 );

    /* Put in M_elvecSource the vector LB^{-1} * M_elvecSource = LB^{-1} F
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBRHS, M_BtB, NBP, M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} rhs is not achieved." );

    /* Put in M_elvecHyb the vector M_BtC^T * M_elvecSource = C^T * A^{-1} * B^T * (B^T * A^{-1} * B)^{-1} * F
       M_elvecHyb is fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBL, NBRHS, NBP, ONE, M_BtC, NBP, M_elvecSource, NBP, ZERO, M_elvecHyb, NBL);

    //........................
    // END OF VECTOR OPERATIONS.
    //........................

    /* Previously the matrix M_CtC is stored only in the lower part, but at the moment there is not
       a function assembleMatrix that store a lower triangular sparse matrix.
       Remind to correct these line in the future. */
    symmetrizeMatrix( UPLO, NBL, M_CtC );

    // Update the hybrid element matrix.
    M_elmatHyb.block(0,0) = M_CtC;

} // staticCondensation

// Locally compute the primal and dual variable.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
localComputePrimalAndDual ()
{

    // Flags for the BLAS and LAPACK routine.

    int INFO[1]         = {0};
    // Number of columns of the right hand side := 1.
    int NBRHS[1]        = {1};
    // Primal variable degrees of freedom.
    int NBP[1]          = { M_primal_FESpace.refFE().nbDof() };
    // Dual variable degrees of freedom.
    int NBU[1]          = { M_dual_FESpace.refFE().nbDof() };
    // Hybrid variable degree of freedom.
    int NBL[1]          = { M_hybrid_FESpace.refFE().nbDof() };
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
    ElemMat::matrix_type A = M_elmatMix.block( 0, 0 );
    ElemMat::matrix_type B = M_elmatMix.block( 0, 1 );
    ElemMat::matrix_type C = M_elmatMix.block( 0, 2 );

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
    dsyrk_ (UPLO, TRANS, NBP, NBU, ONE, B, NBU, ZERO, M_BtB, NBP);

    /* Put in M_BtC the matrix B^T * L^{-T} * L^{-1} * C = B^T * A^{-1} * C
       M_BtC fully stored.
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (TRANS, NOTRANS, NBP, NBL, NBU, ONE, B, NBU, C, NBU, ZERO, M_BtC, NBP);

    /* Put in M_BtB the matrix LB and LB^T where LB and LB^T is the cholesky
       factorization of B^T * A^{-1} * B.
       For more details see http://www.netlib.org/lapack/double/dpotrf.f */
    dpotrf_ (UPLO, NBP, M_BtB, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack factorization of BtB is not achieved." );

    /* Put in M_BtC the matrix LB^{-1} * M_BtC = LB^{-1} * B^T * A^{-1} * C.
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBL, M_BtB, NBP, M_BtC, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation BtC = LB^{-1} BtC is not achieved." );

    //..............................
    //   END OF MATRIX OPERATIONS
    //..............................

    /* Sum up of the previews steps
       A stores L and L^T where L and L^T is the Cholesky factorization of A
       B stores L^{-1} * B
       C stores L^{-1} * C
       M_BtB stores LB and LB^T where LB and LB^T is the factorization of B^T * A^{-1} * B
       M_BtC stores LB^{-1} * B^T * A^{-1} * C */


    //......................
    //    VECTOR OPERATIONS (Computation of Pressure and Velocities)
    //......................

    //...................................
    //  1) Computation of the PRESSURE
    //...................................

    // Clear the source vector.
    M_elvecSource.zero();

    // The source term is computed with a test function in the primal variable space.
    source( M_source,
            M_elvecSource,
            M_primal_FESpace.fe(),
            static_cast<Real>(0), 0 );

    /* Put in M_elvecSource the vector LB^{-1} * M_elvecSource = LB^{-1} * F
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, NOTRANS, NODIAG, NBP, NBRHS, M_BtB, NBP, M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-1} M_elvecSource is not achieved." );

    /* Put in M_elvecSource the vector
       M_BtC * M_elvecHyb + M_elvecSource = LB^{-1} * B^T * A^{-1} * C * lambda_K + LB^{-1} * F
       For more details see http://www.netlib.org/blas/dgemm.f */
    dgemm_ (NOTRANS, NOTRANS, NBP, NBRHS, NBL, MINUSONE, M_BtC, NBP, M_elvecHyb, NBL, MINUSONE, M_elvecSource, NBP);

    /* Put in M_elvecSource the vector LB^{-T} * M_elvecSource, where
       M_elvecSource stores - (B^T * A^{-1} * B)^{-1} * B^T * A^{-1} * C * lambda_K - (B^T * A^{-1} * B)^{-1} * F
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, TRANS, NODIAG, NBP, NBRHS, M_BtB, NBP, M_elvecSource, NBP, INFO);
    ASSERT_PRE( !INFO[0], "Lapack Computation M_elvecSource = LB^{-T} M_elvecSource is not achieved." );

    // Now rhs contains the primal variable for the current element, we must put it in the global vector.

    //.....................................
    //  2) Computation of the VELOCITIES
    //.....................................

    // Clear the element dual vector.
    M_elvecFlux.zero();

    /* Put in M_elvecFlux the vector B * M_elvecSource = L^{-1} * B * primal_K
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBU, NBP, ONE, B, NBU, M_elvecSource, INC, ZERO, M_elvecFlux, INC);
    /* Put in M_elvecFlux the vector
       - C * M_elvecHyb - M_elvecFlux = - L^{-1} * C * lambda_K - L^{-1} * B * primal_K
       For more details see http://www.netlib.org/slatec/lin/dgemv.f */
    dgemv_ (NOTRANS, NBU, NBL, MINUSONE, C, NBL, M_elvecHyb, INC, MINUSONE, M_elvecFlux, INC);

    /* Put in flux the vector
       L^{-T} * M_elvecFlux = - A^{-1} * C^T * lambda_K - A^{-1} * B^T * primal_K
       For more details see http://www.netlib.org/lapack/lapack-3.1.1/SRC/dtrtrs.f */
    dtrtrs_ (UPLO, TRANS, NODIAG, NBU, NBRHS, A, NBU, M_elvecFlux, NBU, INFO);
    ASSERT_PRE(!INFO[0], "Lapack Computation M_elvecFlux = L^{-T} M_elvecFlux is not achieved.");

} // localComputePrimalAndDual

// Update all the variables of the problem.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
updateVariables ()
{

    // Reset the global hybrid matrix.
    M_matrHybrid.reset( new matrix_Type( M_localMap ) );

    // Reset the global right hand side vector
    M_rhs.reset( new vector_Type( M_localMap ) );

    // Reset the global primal vector
    M_primal.reset( new vector_Type( M_primal_FESpace.map() ) );

    // Reset the global dual vector
    M_dual.reset( new vector_Type( M_dual_FESpace.map(), Repeated ) );

    // Reset the global hybrid vector
    M_hybrid.reset( new vector_Type( M_hybrid_FESpace.map() ) );

    // Reset the global residual vector
    M_residual.reset( new vector_Type( M_localMap ) );

} // updateVariables

// Apply the boundary conditions to the global hybrid matrix and the global hybrid right hand side.
template<typename Mesh, typename SolverType>
void
DarcySolver<Mesh, SolverType>::
applyBoundaryConditions ()
{
    // Chrono.
    Chrono chronoBC;

    // Check if the boundary conditions were updated.
    if ( !M_BCh->bcUpdateDone() )
    {

        M_displayer.leaderPrint( "                 Updating the BC               ...          ");
        chronoBC.start();

        // Update the boundary conditions handler. We use the finite element of the boundary of the dual variable.
        M_BCh->bcUpdate( *M_dual_FESpace.mesh(),
                         M_dual_FESpace.feBd(),
                         M_dual_FESpace.dof() );

        chronoBC.stop();
        M_displayer.leaderPrintMax( "done in " , chronoBC.diff() );
    }

    // Ignoring non-local entries, otherwise they are summed up lately.
    vector_Type rhsFull( *M_rhs, Unique );

    M_displayer.leaderPrint( "                 Managing the BC               ...          ");
    chronoBC.start();

    /* Update the global hybrid matrix and the global hybrid right hand side with the boundary conditions.
       It takes care of the current time for DarcySolverTransient derived class. */
    bcManage( *M_matrHybrid,
              rhsFull,
              *M_dual_FESpace.mesh(),
              M_dual_FESpace.dof(),
              *M_BCh,
              M_dual_FESpace.feBd(), 1.,
              M_data.dataTime()->time() );

    chronoBC.stop();
    M_displayer.leaderPrintMax( "done in " , chronoBC.diff() );

    // Save the global hybrid right hand side.
    *M_rhs = rhsFull;

} // applyBoundaryConditions

// Reorder a non full stored symmetric matrix.
template<typename Mesh, typename SolverType>
template<typename matrix>
void
DarcySolver<Mesh, SolverType>::
symmetrizeMatrix ( char* UPLO, int* N, matrix& A  )
{

    // If the matrix is stored in the lower part
    if ( UPLO[0] == 'L')
    {
        for (UInt i(0); i < static_cast<UInt>(N[0]); ++i)
        {
            for (UInt j(i + 1); j < static_cast<UInt>(N[0]); ++j)
            {
                A(i, j) = A(j, i);
            }
        }
    }
    // If the matrix is stored in the upper part
    else
    {
        for (UInt i(0); i < static_cast<UInt>(N[0]); ++i)
        {
            for (UInt j(i + 1); j < static_cast<UInt>(N[0]); ++j)
            {
                A(j, i) = A(i, j);
            }
        }
    }

} //symmetrizeMatrix


} // namespace LifeV


#endif //_DARCYSOLVER_H_

// -*- mode: c++ -*-
