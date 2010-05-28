//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a solver class for the 1D model.
 *
 *  @version 1.0
 *  @author Vincent Martin
 *  @author Tiziano Passerini
 *  @author Lucia Mirabella
 *  @date 01-10-2006
 *
 *  @version 2.0
 *  @author Gilles Fourestey <gilles.fourestey@epfl.ch>
 *  @date 01-08-2009
 *
 *  @version 2.1
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 21-04-2010
 */

#ifndef ONEDMODELSOLVER_H
#define ONEDMODELSOLVER_H

// LIFEV
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/SolverAmesos.hpp>

#include <life/lifefem/FESpace.hpp>

// LIFEV - MATHCARD
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>

namespace LifeV {

//! OneDimensionalModel_Solver - Solver class for the 1D model.
/*!
 *  @author Vincent Martin, Tiziano Passerini, Lucia Mirabella, Gilles Fourestey
 *
 *  ---------------------------------------------------
 *  1 dimensional hyperbolic equation:
 *  ---------------------------------------------------
 *  dU/dt  +  dF(U)/dz  +  S(U) = 0
 *  ---------------------------------------------------
 *
 *  with U = [U1,U2]^T in R^2.
 *
 *  The non linear flux function F(U) and source function S(U)
 *  are quite independant of this solver : they are taken into
 *  account only via two classes that define a vectorial function
 *  and its derivatives.
 *
 *  More precisely:
 *  two functions (M_fluxFun and M_sourceFun) have to be defined
 *  separately to allow the update (_updateFluxDer, _updateSourceDer)
 *  of the corresponding vectors (M_Fluxi, M_diffFluxij,
 *  M_Sourcei, M_diffSrcij). I also separated the treatment of the
 *  parameters that exist in the functions.
 *  Normally, one should be able to create a template to allow
 *  the user to select between different problems (linear, non-linear,
 *  etc.).
 *
 *  -----------------------------------------------------
 *  Solver based on the 2nd order Taylor-Galerkin scheme.
 *  -----------------------------------------------------
 *  (see for instance Formaggia-Veneziani, Mox report no 21 june 2003)
 *
 *  ---------------------------------------------------
 *  Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T )
 *  ---------------------------------------------------
 *
 *  (Un+1, phi) =                             //! massFactor^{-1} * Un+1
 *  (Un, phi)                         //!            mass * U
 *  + dt     * (       Fh(Un), dphi/dz )         //!            grad * F(U)
 *  - dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     //! gradDiffFlux(U) * S(U)
 *  + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //!   divDiffSrc(U) * F(U)
 *  - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //!stiffDiffFlux(U) * F(U)
 *  - dt     * (       Sh(Un), phi )             //!            mass * S(U)
 *  + dt^2/2 * (diffSh(Un) Sh(Un), phi )         //!  massDiffSrc(U) * S(U)
 *  ---------------------------------------------------
 *
 *  Approximation of the unknowns and non-linearities:
 *
 *  Let's define:
 *  a) (phi_i)_{i in nodes} is the basis of P1 (the "hat" functions)
 *  b) (1_{i+1/2})_{i+1/2 in elements} is the basis of P0 (constant per
 *  element). The vertices of the element "i+1/2" are the nodes "i" and "i+1".
 *
 *  Then:
 *
 *  c) Uh    is in P1 : U = sum_{i in nodes} U_i phi_i
 *
 *  d) Fh(U) is in P1 : F(U) = sum_{i in nodes} F(U_i) phi_i
 *  e) diffFh(U) is in P0 :
 *  diffFlux(U) = sum_{i+1/2 in elements} 1/2 { dF/dU(U_i) + dF/dU(U_i+1) } 1_{i+1/2}
 *  (means of the two extremal values of the cell)
 *
 *  We note that d) allows us to define easily
 *  f) dF/dz(U) = sum_{i in nodes} F(U_i) d(phi_i)/dz
 *
 *  g) Sh(U) is in P1 : S(U) = sum_{i in nodes} S(U_i) phi_i
 *  h) diffSh(U) is in P0 :
 *  diffSrc(U) = sum_{i+1/2 in elements} 1/2 { dS/dU(U_i) + dS/dU(U_i+1) } 1_{i+1/2}
 *  (means of the two extremal values of the cell)
 *
 *  -----------------------------------------------------
 *  The option taken here is to define the different tridiagonal matrix
 *  operators (div, grad, mass, stiff) and reconstruct them at each time
 *  step (as they depend on diffFlux and diffSrc). They are thus rebuilt
 *  at the element level and reassembled.
 *  Afterwards, there remains to do only some tridiagonal matrix vector
 *  products to obtain the right hand side.
 *
 *  This procedure might appear a bit memory consuming (there are 18
 *  tridiagonal matrices stored), but it has the advantage of being
 *  very clear. If it is too costly, it should be quite easy to improve
 *  it.
 *  -----------------------------------------------------
 */
class OneDimensionalModel_Solver
{
public:

    //! @name Typedef & Enumerator
    //@{

    typedef OneDimensionalModel_Physics             Physics_Type;
    typedef boost::shared_ptr< Physics_Type >       Physics_PtrType;

    typedef OneDimensionalModel_Flux                Flux_Type;
    typedef boost::shared_ptr< Flux_Type >          Flux_PtrType;

    typedef OneDimensionalModel_Source              Source_Type;
    typedef boost::shared_ptr< Source_Type >        Source_PtrType;

    typedef OneDimensionalModel_Data                Data_Type;
    typedef Data_Type::Mesh_Type                    Mesh_Type;

    typedef FESpace< Mesh_Type, EpetraMap >         FESpace_Type;
    typedef boost::shared_ptr< FESpace_Type >       FESpace_PtrType;

    typedef Epetra_Comm                             Comm_Type;
    typedef boost::shared_ptr< Comm_Type >          Comm_PtrType;

    typedef SolverAmesos                            LinearSolver_Type;
    typedef boost::shared_ptr< LinearSolver_Type >  LinearSolver_PtrType;
    typedef LinearSolver_Type::vector_type          Vector_Type;
    typedef boost::shared_ptr< Vector_Type >        Vector_PtrType;

    typedef LinearSolver_Type::matrix_type          Matrix_Type;
    typedef boost::shared_ptr<Matrix_Type>          Matrix_PtrType;

    typedef std::vector<Vector_Type>                Solution_Type;
    typedef boost::shared_ptr<Solution_Type>        Solution_PtrType;

    typedef boost::function<Real ( const Real&, const Real&, const Real&,
                                   const Real&, const ID& )> Function;

    //! member function pointer: to "automagically" manage postprocess
    typedef void (LifeV::OneDimensionalModel_Solver::*postproc_funptr)( std::string, Real, const Vector_Type& U, std::string );

    typedef std::map< std::string, UInt>::iterator              M_variable_index_iter;
    typedef std::map< std::string, postproc_funptr >::iterator  M_variable_filter_iter;
    typedef std::map< std::string, std::string>::iterator       M_variable_string_iter;

    enum OneDInitializeVar { OneDInitPressure,
                             OneDInitArea,
                             OneDInitFlux,
                             OneDInitRiemann1,
                             OneDInitRiemann2 };

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    /*!
     * Need a call to: setCommunicator, setProblem, setFESpace
     */
    OneDimensionalModel_Solver();

    //! Destructor
    ~OneDimensionalModel_Solver(){}

    //@}


    //! @name Methods
    //@{

    //! setup
    void setup();

    //! Sets initial condition for the unknowns
    /*!
     *  Initialize variables with non constant values in space.
     *  A step profile is approximated by smoothing the edges with a gaussian "bell".
     *  Read options from data file: section [initialize]
     *  var = { A, Q, P, W1, W2} the variable to be initialized
     *  rest_value = rest value
     *  firstnode =,
     *  lastnode = domain of the step
     *  multiplier = step value (= multiplier * value)
     *  amplitude = width of the gaussian bell
     *  X = rest_value ( 1 + multiplier * exp( - ( z - node ) / ( 2 * width^2 ) ) )
     *  where node is alternatively firstnode or lastnode
     */
    void initialize();

    //! Sets initial condition for the unknowns
    /*!
     * @param var if var == "physical" then A  = u10 and Q  = u20
     *            if var == "Riemann"  then W1 = u10 and W2 = u20
     */
    void initialize( const Real& u10, const Real& u20, const std::string& var = "physical" );

    //! Sets initial condition for the unknowns
    /*!
     *  Initialize with vectors containing all the nodal values
     */
    void initialize( const Vector_Type& u10, const Vector_Type& u20 );

    //! Sets initial condition for the unknowns
    /*!
     *  Initialize only Flux ( Area read from OneDNonLinParam )
     */
    void initialize( const Real& u20 );

    //! Update convective term and BC. Then solve the linearized NS system
    /*!
     * @param bcH The BC handler
     * @param Time the time
     * @param TimeStep the time step
     */
    void iterate( OneDimensionalModel_BCHandler& bcH, const Real& Time, const Real& TimeStep );

    //! CFL computation (correct for constant mesh)
    /*!
     * @param TimeStep the time step
     * @return CFL
     */
    Real ComputeCFL( const Real& timeStep ) const;

    //! Save the solution for the next timestep
    void savesol();

    //! Recover the solution at previous timestep (keep unaltered the boundary values)
    void loadsol();

    //! Save results on file
    void postProcess( const Real& time );

    //! Store solutions on matlab readable files
    void output_to_matlab( std::string fname,
                           Real time_val,
                           const Vector_Type& U,
                           std::string vname );

    //! Create matlab file for postprocessing
    void create_movie_file();

    //! Prepare ostringstream buffer to receive the solution before sending it to file stream
    void openFileBuffers();

    //! Write the solution on ostringstream buffers
    void output2FileBuffers( const Real& time_val );

    //! Empty the buffers
    void resetFileBuffers();

    //! Move the pointer to the stored position in the ostringstream
    void seekpFileBuffers();

    //! Store in private variable a desired stream position
    void tellpFileBuffers();

    //! Close the buffers
    void closeFileBuffers();

    //! Print to screen informations on the solver class
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set problem elements
    void setProblem( const Physics_PtrType Physics,
                     const Flux_PtrType    Flux,
                     const Source_PtrType  Source );

    //! Set the communicator
    void setCommunicator( const Comm_PtrType Comm );

    //! Set the FEspace
    void setFESpace( const FESpace_PtrType FESpace );

    //! Set the linear solver
    void setLinearSolver( const LinearSolver_PtrType linearSolver );

    //! Set the Dirichlet boundary conditions (left)
    void setBCValuesLeft( const Real& bcL1, const Real& bcL2 );

    //! Set the Dirichlet boundary conditions (right)
    void setBCValuesRight( const Real& bcR1, const Real& bcR2 );

    //@}


    //! @name Get Methods
    //@{

    //! Return the solution at current time step (U)
    const Solution_PtrType U_thistime() const;

    //! Return the solution at current time step (Area)
    const Vector_Type& U1_thistime() const;

    //! Return the solution at current time step (Flux)
    const Vector_Type& U2_thistime() const;

    //! Return the Riemann invariant W1 at current time step
    const Vector_Type& W1_thistime() const;

    //! Return the Riemann invarant W2 at current time step
    const Vector_Type& W2_thistime() const;

    //! Return the solution at current time step (Pressure)
    const Vector_Type& P_thistime() const;

    //! Get the Physics function
    const Physics_PtrType& Physics() const;

    //! Get the flux function
    const Flux_PtrType& Flux() const;

    //! Get the source function
    const Source_PtrType& Source() const;

    //! Get the left node identifier
    const UInt& LeftNodeId() const;

    //! Get the left internal node (neighboring node)
    const UInt& LeftInternalNodeId() const;

    //! Get the right node identifier
    const UInt& RightNodeId() const;

    //! Get the right internal node (neighboring node)
    const UInt& RightInternalNodeId() const;

    //! Get the Dirichlet boundary conditions (left)
    Container2D_Type BCValuesLeft() const;

    //! Get the value at neighboring node (left)
    Container2D_Type BCValuesInternalLeft() const;

    //! Get the Dirichlet boundary conditions (right)
    Container2D_Type BCValuesRight() const;

    //! Get the value at neighboring node (right)
    Container2D_Type BCValuesInternalRight() const;

    //! Return the selection solution (P, A, Q, W1, W2)
    Real value(std::string var, UInt pos) const;

    //! Return the value of a quantity (P, A, Q, W1, W2) on a specified boundary.
    /*!
     *  Given a bcType and a bcSide it return the value of the quantity.
     *  @param bcType Type of the asked boundary value.
     *  @param bcSide Side of the boundary.
     *  @return value of the quantity on the specified side.
     */
    Real BoundaryValue( const OneD_BC& bcType, const OneD_BCSide& bcSide ) const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Update the pressure
    void _updatePressure( const Real& TimeStep );

    //! Update the right hand side and the matrices
    /*!
     *  @param TimeStep The time step.
     */
    void _updateSystem( const Real& TimeStep );

    //! Update the P1 flux vector from U: M_Fluxi = F_h(Un) i=1,2 (works only for P1Seg elements)
    void _updateFlux();

    //! Call _updateFlux and update the P0 derivative of flux vector from U:
    /*!
     *  M_diffFluxij = dF_h/dU(Un) i,j=1,2
     *  M_diffFluxij(elem) = 1/2 [ dF/dU(U(node1(elem))) + dF/dU(U(node2(elem))) ]
     *
     *  (mean value of the two extremal values of dF/dU)
     *  BEWARE: works only for P1Seg elements
     */
    void _updateFluxDer();

    //! Update the P1 source vector from U: M_Sourcei = S_h(Un) i=1,2 (works only for P1Seg elements)
    void _updateSource();

    //! Call _updateSource and update the P0 derivative of source vector from U:
    /*!
     *  M_diffSrcij = dS_h/dU(Un) i,j=1,2
     *  M_diffSrcij(elem) = 1/2 [ dS/dU(U(node1(elem))) + dS/dU(U(node2(elem))) ]
     *
     *  (mean value of the two extremal values of dS/dU)
     *  BEWARE: works only for P1Seg elements
     */
    void _updateSourceDer();

    //! Update the matrices
    /*!
     *  M_massMatrixDiffSrcij, M_stiffMatrixDiffFluxij
     *  M_gradMatrixDiffFluxij, and M_divMatrixDiffSrcij (i,j=1,2)
     *
     *  from the values of diffFlux(Un) and diffSrc(Un)
     *  that are computed with _updateMatrixCoefficients.
     *
     *  call of  _updateMatrixCoefficients,
     *  _updateElemMatrices and _assemble_matrices.
     */
    void _updateMatrices();

    //! Update the coefficients (from the flux, source functions and their derivatives)
    void _updateMatrixCoefficients( const UInt& ii, const UInt& jj, const UInt& iedge);

    //! Update the element matrices with the current element
    void _updateElemMatrices();

    //! Assemble the matrices
    int _assemble_matrices(const UInt& ii, const UInt& jj );

    //! Update the vectors to take into account Dirichlet BC.
    /*!
     *  Modify the vectors to take into account
     *  the Dirichlet boundary conditions
     *  (works for P1Seg and canonic numbering!)
     */
    void _updateBCDirichletVector();

    //! Update the matrices to take into account Dirichlet BC.
    /*!
     *  Modify the matrix to take into account
     *  the Dirichlet boundary conditions
     *  (works for P1Seg and canonic numbering!)
     */
    void _updateBCDirichletMatrix( Matrix_Type& mat );

    //! Apply the inertial Flux correction:
    /*!
     *  We use a finite element scheme for the correction term:
     *  given the solution of Taylor-Galerkin scheme, solve
     *  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  ( m / rho ) * ( dQtildeh(n+1)/dz, dphi/dz )  //! stiff * Qtilde(U)
     *  = ( m / rho ) *       ( dQhath(n+1)/dz, dphi/dz )  //! stiff * Qhat(U)
     *
     *  m = rho_w h0 / ( 2 sqrt(pi) sqrt(A0) )
     */
    Vector_Type                   _correct_flux_inertial    (const Vector_Type &);

    //! Apply the viscoelastic Flux correction:
    /*!
     *  We use a finite element scheme for the correction term:
     *  given the solution of Taylor-Galerkin scheme, solve
     *  ( 1/(dt*Ah(n+1)) Qtildeh(n+1), phi) +        //! 1/A * massFactor^{-1} * Un+1
     *  ( gamma / rho ) *     ( dQtildeh(n+1)/dz, dphi/dz )  //! stiff * Qtilde(U)
     *  = - ( gamma / rho ) * ( dQhath(n+1)/dz, dphi/dz )  //! stiff * Qhat(U)
     *
     *  gamma = gamma_tilde / ( 2 sqrt(pi) )
     */
    Vector_Type                   _correct_flux_viscoelastic(const Vector_Type&, const Real& TimeStep );

    //! Apply the longitudinal Flux correction:
    /*!
     *  We use a finite element scheme for the correction term:
     *  given the solution of Taylor-Galerkin scheme, solve
     *  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  = ( 1/Ah(n+1) Qtildeh(n), phi) +             //! 1/A * massFactor^{-1} * Un+1
     *  + ( a / rho ) *       ( d3Ahath(n+1)/dz3, phi )  //! mass * d3Ahat(U)/dz
     */
    Vector_Type                   _correct_flux_longitudinal( );

    //! L2 Projection of the second derivative of Q over P1 space.
    //ScalVec                       _compute_d2Q_dx2( const ScalVec& );

    //@}

    Physics_PtrType                    M_Physics;
    Flux_PtrType                       M_Flux;
    Source_PtrType                     M_Source;
    FESpace_PtrType                    M_FESpace;
    Comm_PtrType                       M_Comm;
    Displayer                          M_Displayer;

    UInt                               M_leftNodeId;
    UInt                               M_leftInternalNodeId;
    UInt                               M_rightNodeId;
    UInt                               M_rightInternalNodeId;

    //! coefficient in front of the corresponding M_elmat*
    Real                               M_coeffMass;
    Real                               M_coeffStiff;
    Real                               M_coeffGrad;
    Real                               M_coeffDiv;

    boost::shared_ptr< ElemMat >       M_elmatMass;  //!< element mass matrix
    boost::shared_ptr< ElemMat >       M_elmatStiff; //!< element stiffness matrix
    boost::shared_ptr< ElemMat >       M_elmatGrad;  //!< element gradient matrix
    boost::shared_ptr< ElemMat >       M_elmatDiv;   //!< element divergence matrix

    //! Unknowns at present time step
    /*!
      U is a vector of ScalVec:
      U[0] = A, U[1] = Q
      U[2] = W1, U[3] = W2
      Other components of U may contain additional variables such as the pressure
    */
    Solution_PtrType                   M_U_thistime;

    //! Unknowns at previous time step (see savesol() )
    Solution_Type                      M_U_prevtime;
    Solution_Type                      M_U_2prevtime;

    //! Right hand sides of the linear system i: "mass * M_Ui = M_rhsi"
    std::vector<Vector_Type>           M_rhs;

    //! Flux F(U) (in P1)
    std::vector<Vector_Type>           M_FluxVector;

    //! diffFlux = dF(U)/dU (in P0)
    std::vector<ScalVec>               M_diffFlux;

    //! Source term S (in P1)
    std::vector<Vector_Type>           M_SourceVector;

    //! diffSrc = dSource(U)/dU (in P0)
    std::vector<ScalVec>               M_diffSrc;

    //! tridiagonal mass matrix
    Matrix_PtrType                     M_massMatrix;

    //! tridiagonal mass matrices multiplied by diffSrcij
    std::vector<Matrix_PtrType >       M_massMatrixDiffSrc;

    //! tridiagonal stiffness matrices multiplied by diffFluxij
    std::vector<Matrix_PtrType >       M_stiffMatrixDiffFlux;

    //! tridiagonal gradient matrix
    Matrix_PtrType                     M_gradMatrix;

    //! tridiagonal gradient matrices multiplied by diffFluxij
    std::vector<Matrix_PtrType >       M_gradMatrixDiffFlux;

    //! tridiagonal divergence matrices multiplied by diffSrcij
    std::vector<Matrix_PtrType >       M_divMatrixDiffSrc;

    //! The linear solver
    boost::shared_ptr<LinearSolver_Type>     M_linearSolver;

    //! ostringstream buffers to store the solution before writing it to file
    std::map<std::string, boost::shared_ptr<std::ostringstream> > M_post_process_buffer;

    //! position of the put pointer to ostringstream buffers
    std::map<std::string, long>  M_post_process_buffer_offset;

    //! Trick to use strings in C++ "switch" construct
    std::map<std::string, OneDInitializeVar> M_oneDstring2initializeVarMap;

    Container2D_Type                        M_bcDirLeft;  //! first -> U1, second ->U2
    Container2D_Type                        M_bcDirRight; //

    // maps associating each unknown to a string, a numeric index, a function pointer
    std::map< std::string, UInt>                                        M_variable_index_map;
    std::map< std::string, postproc_funptr>                             M_variable_filter_map;
};

}

#endif // ONEDMODELSOLVER_H
