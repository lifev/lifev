/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
/*!
  \file oneDModelSolver.hpp
  \author Vincent Martin
  \author Tiziano Passerini
  \author Lucia Mirabella

  \date 10/2006
  \version 1.1

  \author Gilles Fourestey
  \date 08/2009
  \version 2.0

  \brief This file contains a solver class for the 1D model.

  ---------------------------------------------------
  1 dimensional hyperbolic equation:
  ---------------------------------------------------
  dU/dt  +  dF(U)/dz  +  S(U) = 0
  ---------------------------------------------------

  with U = [U1,U2]^T in R^2.

  The non linear flux function F(U) and source function S(U)
  are quite independant of this solver : they are taken into
  account only via two classes that define a vectorial function
  and its derivatives.

  More precisely:
  two functions (M_fluxFun and M_sourceFun) have to be defined
  separately to allow the update (_updateFluxDer, _updateSourceDer)
  of the corresponding vectors (M_Fluxi, M_diffFluxij,
  M_Sourcei, M_diffSrcij). I also separated the treatment of the
  parameters that exist in the functions.
  Normally, one should be able to create a template to allow
  the user to select between different problems (linear, non-linear,
  etc.).


  -----------------------------------------------------
  Solver based on the 2nd order Taylor-Galerkin scheme.
  -----------------------------------------------------
  (see for instance Formaggia-Veneziani, Mox report no 21 june 2003)

  ---------------------------------------------------
  Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T )
  ---------------------------------------------------

  (Un+1, phi) =                             //! massFactor^{-1} * Un+1
  (Un, phi)                         //!            mass * U
  + dt     * (       Fh(Un), dphi/dz )         //!            grad * F(U)
  - dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     //! gradDiffFlux(U) * S(U)
  + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //!   divDiffSrc(U) * F(U)
  - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //!stiffDiffFlux(U) * F(U)
  - dt     * (       Sh(Un), phi )             //!            mass * S(U)
  + dt^2/2 * (diffSh(Un) Sh(Un), phi )         //!  massDiffSrc(U) * S(U)
  ---------------------------------------------------

  Approximation of the unknowns and non-linearities:

  Let's define:
  a) (phi_i)_{i in nodes} is the basis of P1 (the "hat" functions)
  b) (1_{i+1/2})_{i+1/2 in elements} is the basis of P0 (constant per
  element). The vertices of the element "i+1/2" are the nodes "i" and "i+1".

  Then:

  c) Uh    is in P1 : U = sum_{i in nodes} U_i phi_i

  d) Fh(U) is in P1 : F(U) = sum_{i in nodes} F(U_i) phi_i
  e) diffFh(U) is in P0 :
  diffFlux(U) = sum_{i+1/2 in elements} 1/2 { dF/dU(U_i) + dF/dU(U_i+1) } 1_{i+1/2}
  (means of the two extremal values of the cell)

  We note that d) allows us to define easily
  f) dF/dz(U) = sum_{i in nodes} F(U_i) d(phi_i)/dz

  g) Sh(U) is in P1 : S(U) = sum_{i in nodes} S(U_i) phi_i
  h) diffSh(U) is in P0 :
  diffSrc(U) = sum_{i+1/2 in elements} 1/2 { dS/dU(U_i) + dS/dU(U_i+1) } 1_{i+1/2}
  (means of the two extremal values of the cell)


  -----------------------------------------------------
  The option taken here is to define the different tridiagonal matrix
  operators (div, grad, mass, stiff) and reconstruct them at each time
  step (as they depend on diffFlux and diffSrc). They are thus rebuilt
  at the element level and reassembled.
  Afterwards, there remains to do only some tridiagonal matrix vector
  products to obtain the right hand side.

  This procedure might appear a bit memory consuming (there are 18
  tridiagonal matrices stored), but it has the advantage of being
  very clear. If it is too costly, it should be quite easy to improve
  it.

  -----------------------------------------------------
*/

#ifndef _ONEDMODELSOLVER_H_
#define _ONEDMODELSOLVER_H_

#include <string>
#include <sstream>

//
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/SolverAmesos.hpp>

#include <life/lifefem/assemb.hpp>
#include <life/lifecore/chrono.hpp>

//#include <life/lifearray/tridiagMatrix.hpp>
#include <life/lifealg/triDiagCholesky.hpp>
#include <life/lifealg/triDiagLU.hpp>


#include <lifemc/lifesolver/vectorFunction1D.hpp>

#include <life/lifefem/FESpace.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/function.hpp>

#include <lifemc/lifesolver/oneDNonLinModelParam.hpp>
#include <lifemc/lifesolver/oneDBCHandler.hpp>
#include <lifemc/lifesolver/dataOneDModel.hpp>

namespace ublas = boost::numeric::ublas;

namespace LifeV
{
enum OneDInitializeVar { OneDInitPressure,
                         OneDInitArea,
                         OneDInitFlux,
                         OneDInitReimann1,
                         OneDInitReimann2 };




/*!
  \class OneDModelSolver

  This class contains a solver class for the 1D model.

*/


//

template<typename Params, typename Flux, typename Source>
class OneDModelSolver
    //        :public OneDModelHandler
{
public:

    /*! \name Typedefs
     */
    //@{
    typedef Params                            param_type;
    typedef Flux                              flux_type;
    typedef Source                            source_type;

    typedef DataOneDModel                     data_type;

    typedef data_type::mesh_raw_type          Mesh;

    typedef data_type::Vec2D                  Vec2D;
//     typedef data_type::ScalVec                ScalVec;
//     typedef data_type::ScalVec_vector         ScalVec_vector;
//     typedef data_type::ScalVec_vector_range   ScalVec_vector_range;

//     typedef ublas::bounded_array<Real, 2>        Vec2D;
//     //typedef ublas::vector<Real>                  Vec2D;
//     //! Vector containing the nodal values of the scalar unknown
//     typedef ublas::vector<double>                ScalVec;
//     //! ublas::vector of Vectors
//     typedef ublas::vector<ScalVec>               ScalVec_vector;
//     //! ublas::vector_range of ScalVec_vector
//     typedef ublas::vector_range<ScalVec_vector>  ScalVec_vector_range;


    typedef SolverAmesos                      solver_type;
    typedef typename solver_type::vector_type vector_type;

    typedef typename solver_type::matrix_type matrix_type;
    typedef boost::shared_ptr<matrix_type>    matrix_ptrtype;



    //! 2D vector
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> Function;
    //@}

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param onedparam Variable containing parameters for OneD model
    */

    OneDModelSolver( const data_type&                   dataType,
                     const Params&                      parameters,
                     FESpace<Mesh, EpetraMap>&          FESpace,
                     Epetra_Comm&                       comm); //,

    ~OneDModelSolver(){}

    //! set up
    void setup(const GetPot& datafile, const std::string& section = "");

    //! return the solution at current time step (Area)
    const vector_type& U1_thistime()       const { return M_U_thistime[0]; }

    //! return the solution at current time step (Flux)
    const vector_type& U2_thistime()       const { return M_U_thistime[1]; }

    //! return the reimann invariant W1 at current time step
    const vector_type& W1_thistime()       const { return M_U_thistime[2]; }

    //! return the reimann invarant W2 at current time step
    const vector_type& W2_thistime()       const { return M_U_thistime[3]; }

    //! return the solution at current time step (U)
    const std::vector<vector_type>& U_thistime() const { return M_U_thistime; }

    //! return a const reference to parameter class
    const Params& oneDParams()         const { return M_oneDParams; }

    //! return the BC handler
    OneDBCHandler<Flux>& bcH() { return M_bcH;}

    //! Sets initial condition for the unknowns
    /*!
      \param var if var == "physical" then A = u10 and Q = u20
      if var == "reimann" then W1 = u10 and W2 = u20
    */
    void initialize(const Real& u10, const Real& u20,
                    const std::string& var = "physical");

    //! Sets initial condition for the unknowns
    /*!
      Give vectors (all the nodal values)
    */
    void initialize(const Vector& u10, const Vector& u20);

    //! Sets initial condition for the unknown
    //! (incremental approach): the initial time is t0, the time step dt
    void initialize(const Function& c0, Real t0, Real dt);

    //! Sets initial condition for the unknowns from file
    void initialize(const std::string & vname);

    //! Sets initial condition for the unknown
    void initialize(const Real& u20);

    //! Initialize with non constant (step) data
    void initialize(const GetPot& data_file, const std::string section = "");

    //! Save and recover solution at current time step
    void savesol();
    void loadsol();

    //! Update the right hand side  for time advancing
    void timeAdvance( const Real& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate( const Real& time , const int& count);

    //! Save results on file
    void postProcess( const Real& time );

    //! get the Dirichlet boundary conditions (left)
    Vec2D BCValuesLeft() const;

    //! get the value at neighboring node (left)
    Vec2D BCValuesInternalLeft() const;

    //! get the Dirichlet boundary conditions (right)
    Vec2D BCValuesRight() const;

    //! get the value at neighboring node (right)
    Vec2D BCValuesInternalRight() const;

    //! return the selection solution (P, A, Q, W1, W2)
    Real value(std::string var, UInt pos) const;


    //! set the Dirichlet boundary conditions (left)
    void setBCValuesLeft( const Real& bcL1, const Real& bcL2 );

    //! set the Dirichlet boundary conditions (right)
    void setBCValuesRight( const Real& bcR1, const Real& bcR2 );

    //! set left bctype to internal node
    void setBCLeft_internalnode();

    //! set right bctype to internal node
    void setBCRight_internalnode();

    //! get the flux function
    Flux const& FluxFun() const;
    //! get the source function
    Source const& SourceFun() const;

    //! get the left edge
    Mesh::EdgeType LeftEdge() const  {return M_leftEdge;}
    //! get the right edge
    Mesh::EdgeType RightEdge() const {return M_rightEdge;};

    //! get the left node identifier
    UInt LeftNodeId() const;
    //! get the left internal node (neighboring node)
    UInt LeftInternalNodeId() const;
    //! get the right node identifier
    UInt RightNodeId() const;
    //! get the right internal node (neighboring node)
    UInt RightInternalNodeId() const;

    //! simple cfl computation (correct for constant mesh)
    void CheckCFL() const;

    //! plotting
    void gplot();

    //! writer in plotmtv format
    void output_to_plotmtv(std::string fname, Real time_val,
                           //        const std::vector< Point1D >& ptlist,
                           const ScalVec& U,
                           std::string vname );

    //! store solutions on matlab readable files
    void output_to_matlab( std::string fname,
                           Real time_val,
                           const vector_type& U,
                           std::string vname );

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

    //! Create matlab file for postprocessing
    void create_movie_file();

    //! Print to screen informations on the solver class
    void showMe( std::ostream& c = std::cout, UInt verbose = false );

    //! timestep getters
    const Real timestep() const {return M_data.timestep();}

    //! left external point position
    const Real xLeft()    const {return M_data.xLeft();}

    //! left external point position
    const Real xRight()   const {return M_data.xRight();}

private:

    //! the data
    const data_type&                   M_data;

    //! the parameters
    Params                             M_oneDParams;

    //! The FESpace
    FESpace<Mesh, EpetraMap>&          M_FESpace;

    //! The FESpace dof size
    UInt                               M_dimDof;

    // const LinearSimpleParam& M_oneDParam;
    EpetraMap                          M_localMap;

    //! the flux function
    Flux                               M_fluxFun;

    //! the source function
    Source                             M_sourceFun ;

    const UInt                         M_leftNodeId;
    const UInt                         M_leftInternalNodeId;
    const UInt                         M_rightNodeId;
    const UInt                         M_rightInternalNodeId;

    //! boundary edges
    const Mesh::EdgeType               M_leftEdge;
    const Mesh::EdgeType               M_rightEdge;

    //! coefficient in front of the corresponding M_elmat*
    Real                               M_coeffMass;
    Real                               M_coeffStiff;
    Real                               M_coeffGrad;
    Real                               M_coeffDiv;

    ElemMat                            M_elmatMass;  //!< element mass matrix
    ElemMat                            M_elmatStiff; //!< element stiffness matrix
    ElemMat                            M_elmatGrad;  //!< element gradient matrix
    ElemMat                            M_elmatDiv;   //!< element divergence matrix

    //  ElemVec M_elvec; // Elementary right hand side

    //! Unknowns at present time step
    /*!
      U is a vector of ScalVec:
      U[0] = A, U[1] = Q
      U[2] = W1, U[3] = W2
      Other components of U may contain additional variables such as the pressure
    */

    std::vector<vector_type>               M_U_thistime;

    //! Unknowns at previous time step (see savesol() )
    std::vector<vector_type>               M_U_prevtime;
    std::vector<vector_type>               M_U_2prevtime;


    //
//     std::vector<ScalVec>               M_U_thistime;

//     //! Unknowns at previous time step (see savesol() )
//     std::vector<ScalVec>               M_U_prevtime;
//     std::vector<ScalVec>               M_U_2prevtime;

    //! Right hand sides of the linear system i: "mass * M_Ui = M_rhsi"
    std::vector<vector_type>           M_rhs;

    //! Solution
    std::vector<vector_type>           M_sol;

    //! Flux F(U) (in P1)
    std::vector<vector_type>           M_Flux;

    //! diffFlux = dF(U)/dU (in P0)
    std::vector<ScalVec>               M_diffFlux;

    //! diffSrc = dSource(U)/dU (in P0)
    std::vector<ScalVec>               M_diffSrc;

    //! Source term S (in P1)
    std::vector<vector_type>           M_Source;

    //! tridiagonal mass matrix
    matrix_type                        M_massMatrix;

    //! factorized tridiagonal mass matrix !NOT NEEDED ANYMOORE
    //matrix_type                        M_factorMassMatrix;

    //! tridiagonal mass matrices multiplied by diffSrcij
    std::vector<matrix_ptrtype >       M_massMatrixDiffSrc;

    //! tridiagonal stiffness matrices multiplied by diffFluxij
    std::vector<matrix_ptrtype >       M_stiffMatrixDiffFlux;

    //! tridiagonal gradient matrix
    matrix_type                        M_gradMatrix;

    //! tridiagonal gradient matrices multiplied by diffFluxij
    std::vector<matrix_ptrtype >       M_gradMatrixDiffFlux;

    //! tridiagonal divergence matrices multiplied by diffSrcij
    std::vector<matrix_ptrtype >       M_divMatrixDiffSrc;

    //! The linear solver
    solver_type                        M_linearSolver;

    //! Update the coefficients
    //! (from the flux, source functions and their derivatives)
    void _updateMatrixCoefficients(const UInt& ii, const UInt& jj,
                                   const UInt& iedge);

    //! Update the element matrices with the current element
    void _updateElemMatrices();

    //! assemble the matrices
    int _assemble_matrices(const UInt& ii, const UInt& jj );

    /*! update the matrices
      M_massMatrixDiffSrcij, M_stiffMatrixDiffFluxij
      M_gradMatrixDiffFluxij, and M_divMatrixDiffSrcij (i,j=1,2)

      from the values of diffFlux(Un) and diffSrc(Un)
      that are computed with _updateMatrixCoefficients.

      call of  _updateMatrixCoefficients,
      _updateElemMatrices and _assemble_matrices.
    */
    void _updateMatrices();

    /*! modify the matrix to take into account
      the Dirichlet boundary conditions
      (works for P1Seg and canonic numbering!)
    */
    void _updateBCDirichletMatrix( matrix_type& mat );

    /*! modify the vector to take into account
      the Dirichlet boundary conditions
      (works for P1Seg and canonic numbering!)
    */
    void _updateBCDirichletVector();

    /*! compute the M_bcDirLeft and M_bcDirRight
      from the external boundary condition
      and the compatibility condition.
      use the extrapolation of the pseudo-characteristics

      (from the solution at time n, obtain the
      value of the linearized characteristics
      at time n+1.)
      Used as Compatibility condition.
    */
    void _computeBC( const Real& time_val );

    //! Axpy product for 2D vectors (pairs)
    //! Axpy(alpha, x, beta, y) -> y = a*A*x + beta*y
    void Axpy(const Vec2D& line1, const Vec2D& line2,
              const Real& alpha,  const Vec2D& x,
              const Real& beta,   Vec2D& y) const;

    //! update the P1 flux vector from U: M_Fluxi = F_h(Un) i=1,2
    void                          _updateFlux();
    //! update the P1 source vector from U: M_Sourcei = S_h(Un) i=1,2
    void                          _updateSource();

    //! call _updateFlux and update the P0 derivative of flux vector from U:
    //! M_diffFluxij = dF_h/dU(Un) i,j=1,2
    void                          _updateFluxDer();
    //! call _updateSource and update the P0 derivative of source vector from U:
    //! M_diffSrcij = dS_h/dU(Un) i,j=1,2
    void                          _updateSourceDer();

    //! to manage non elastic wall behaviour
    vector_type                   _correct_flux_inertial    (const vector_type &);
    vector_type                   _correct_flux_viscoelastic(const vector_type &);
    vector_type                   _correct_flux_longitudinal( );

    ScalVec                       _compute_d2Q_dx2( const ScalVec& );

    //! ostringstream buffers to store the solution before writing it to file
    std::map<std::string, boost::shared_ptr<std::ostringstream> >
                                 M_post_process_buffer;

    //! position of the put pointer to ostringstream buffers
    std::map<std::string, long>  M_post_process_buffer_offset;

    //! Class to manage boundary conditions
    OneDBCHandler<Flux>          M_bcH;

    //! Trick to use strings in C++ "switch" construct
    std::map<std::string, OneDInitializeVar> M_oneDstring2initializeVarMap;

    //
    Vec2D                        M_bcDirLeft;  //! first -> U1, second ->U2
    Vec2D                        M_bcDirRight; //

    //! boolean: show CFL value during computations?
    bool                         M_CFL;

    //! boolean: use alternate solver (UpWind)?
    bool                         M_UW;

    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool                         M_inertial_wall;
    bool                         M_viscoelastic_wall;
    bool                         M_linearize_string_model;
    bool                         M_linearize_equations;
    bool                         M_longitudinal_wall;

    //! boolean: compute second spatial derivative of flux?
    bool                         M_flux_second_der;

    //! approximation of pressure temporal derivative: how many time steps?
    int                          M_dP_dt_steps;

    //! plotting using graceplot

    //GracePlot M_GracePlot; //!< for plotting

    //! member function pointer: to "automagically" manage postprocess
    typedef void (LifeV::OneDModelSolver<Params, Flux, Source>::*postproc_funptr)( std::string,
        Real, const vector_type& U, std::string );

    // maps associating each unknown to a string, a numeric index, a function pointer
    std::map< std::string, std::string>                                 M_variable_string_map;
    std::map< std::string, UInt>                                        M_variable_index_map;
    std::map< std::string, postproc_funptr>                             M_variable_filter_map;
    typedef std::map< std::string, UInt>::iterator                      M_variable_index_iter;
    typedef typename std::map< std::string, postproc_funptr >::iterator M_variable_filter_iter;
    typedef std::map< std::string, std::string>::iterator               M_variable_string_iter;


};


// ***********************************
// IMPLEMENTATION
// ***********************************


template< class Params, class Flux, class Source >
OneDModelSolver<Params, Flux, Source>::
OneDModelSolver( const data_type&                   dataType,
                 const Params&                      parameters,
                 FESpace<Mesh, EpetraMap>&          FESpace,
                 Epetra_Comm&                       comm):
//,OneDModelSolver(const GetPot& data_file):
        //        OneDModelHandler        (data_file),
        M_data                  (dataType),
        M_oneDParams            (parameters),
        M_FESpace               (FESpace),
        M_dimDof                (FESpace.dim()),
        M_localMap              (FESpace.map()),
        M_fluxFun               (M_oneDParams),
        M_sourceFun             (M_oneDParams),
        //! id of left and right bc nodes
        M_leftNodeId            ( 1 ),
        M_leftInternalNodeId    ( M_leftNodeId  + 1 ),
        M_rightNodeId           ( M_FESpace.dim()   ),
        M_rightInternalNodeId   ( M_rightNodeId - 1 ),
        //! boundary edges
        M_leftEdge              ( M_FESpace.mesh()->edgeList( 1 ) ),
        M_rightEdge             ( M_FESpace.mesh()->edgeList( M_FESpace.dim() - 1 ) ),
        // elementary matrices
        M_elmatMass             (M_FESpace.fe().nbNode,1,1),
        M_elmatStiff            (M_FESpace.fe().nbNode,1,1),
        M_elmatGrad             (M_FESpace.fe().nbNode,1,1),
        M_elmatDiv              (M_FESpace.fe().nbNode,1,1),
        // vectorial unknowns and rhs
        M_U_thistime            (5, vector_type(M_localMap)), // size should be at least 4
        M_U_prevtime            (2, vector_type(M_localMap)), // probably useless - could be components in U_thistime
        M_U_2prevtime           (2, vector_type(M_localMap)),
//         M_U_thistime            (0 /*, vector_type(M_localMap)*/), // size should be at least 4
//         M_U_prevtime            (2 /*, vector_type(M_localMap)*/), // probably useless - could be components in U_thistime
//         M_U_2prevtime           (2 /*, vector_type(M_localMap)*/),
        //
        M_rhs                   (2, vector_type(M_localMap)),
        // vectors and matrices of the non-linear function
        M_Flux                  (2, vector_type(M_localMap)),
        M_diffFlux              (4 /*, vector_type(M_localMap)*/),
        M_diffSrc               (4 /*, vector_type(M_localMap)*/),
        //
        M_Source                (2, vector_type(M_localMap)),
        // mass matrix (to be inverted)
        M_massMatrix            (M_localMap),
        //M_factorMassMatrix      (M_localMap),
        //M_tridiagSlv            (M_localMap),
        // matrices used to build the rhs
        M_gradMatrix            (M_localMap),
        // The linear solver
        M_linearSolver          ( comm ),
        // Handle boundary conditions
        M_bcH                   (M_U_thistime, /*M_postproc_variable,*/
                                 M_fluxFun,
                                 M_dimDof),
        M_bcDirLeft             (2),
        M_bcDirRight            (2)
{
}

/*! Set up */
template< class Params, class Flux, class Source >
void OneDModelSolver<Params, Flux, Source>::setup(const GetPot& data_file, const std::string& section)
{
    M_CFL                    = data_file((section + "miscellaneous/showCFL").data(),                        0);
    M_UW                     = data_file((section + "miscellaneous/alternate_solver").data(),               0);
    M_inertial_wall          = data_file((section + "miscellaneous/inertial_wall").data(),                  0);
    M_viscoelastic_wall      = data_file((section + "miscellaneous/viscoelastic_wall").data(),              0);
    M_linearize_string_model = data_file((section + "miscellaneous/linearize_string_model").data(),         1);
    M_linearize_equations    = data_file((section + "miscellaneous/linearize_equations").data(),            0);
    M_longitudinal_wall      = data_file((section + "miscellaneous/longitudinal_wall").data(),              0);
    M_flux_second_der        = data_file((section + "miscellaneous/compute_flux_second_derivative").data(), 0);
    M_dP_dt_steps            = data_file((section + "miscellaneous/pressure_derivative_steps").data(),      1);


        // These maps allow a more readable definition of the variables
    // the model is describing


    M_variable_string_map["A"]              = "Area";
    M_variable_string_map["Q"]              = "Flow Rate";
    M_variable_string_map["P"]              = "Pressure";
    M_variable_string_map["W1"]             = "First Reimann Invariant";
    M_variable_string_map["W2"]             = "Second Reimann Invariant";
    M_variable_string_map["Q_visc"]         = "Flow Rate correction (viscoelastic)";
    M_variable_string_map["Q_inertial"]     = "Flow Rate correction (inertial)";
    M_variable_string_map["Q_longitudinal"] = "Flow Rate correction (longitudinal)";
    M_variable_string_map["d2Q_dx2"]        = "Flow Rate second derivative";

    // These maps allow to use the switch... case construct
    // with string variables
    M_oneDstring2initializeVarMap["A"]      = OneDInitArea;
    M_oneDstring2initializeVarMap["Q"]      = OneDInitFlux;
    M_oneDstring2initializeVarMap["W1"]     = OneDInitReimann1;
    M_oneDstring2initializeVarMap["W2"]     = OneDInitReimann2;
    M_oneDstring2initializeVarMap["P"]      = OneDInitPressure;

    Debug( 6310 ) << "[OneDModelSolver] O-  Nb of unknowns: " << M_FESpace.dim()     << "\n";
    Debug( 6310 ) << "[OneDModelSolver] O-  Computing the constant matrices... \n";
    Debug( 6310 ) << "[OneDModelSolver] O-  Adopting a"
                  << std::string( M_viscoelastic_wall ? " viscoelastic " : "n elastic " )
                  << "model for vessel wall... \n";

    Chrono chrono;
    chrono.start();

    //! Vector and Matrices initialization
    //M_massMatrix.zero();
    // M_factorMassMatrix.zero();
    //M_factorMassMatrix.zero();
    //M_gradMatrix.zero();

    // insert variables!
    UInt nvar(0);
    M_variable_index_map.insert( make_pair("A",  nvar++ ) );
    M_variable_index_map.insert( make_pair("Q",  nvar++ ) );
    M_variable_index_map.insert( make_pair("W1", nvar++ ) );
    M_variable_index_map.insert( make_pair("W2", nvar++ ) );

    M_variable_index_map.insert( make_pair("P",  nvar++ ) );



    //! correction flux with viscoelastic term
    if( M_viscoelastic_wall )
        {
            // correction on the flux due to the viscoelastic term
            M_variable_index_map.insert( make_pair("Q_visc",  nvar++ ) );
            // viscoelastic contribution to the pressure
            M_variable_index_map.insert( make_pair("P_visc",  nvar++ ) );
            // elastic contribution to the pressure
            M_variable_index_map.insert( make_pair("P_elast", nvar++ ) );
            // time derivative of the section area
            M_variable_index_map.insert( make_pair("dA_dt",   nvar++ ) );
        }

    //! flux second derivative
    if( M_flux_second_der )
        M_variable_index_map.insert( make_pair("d2Q_dx2", nvar++ ) );

    //! correction flux with inertial term
    if( M_inertial_wall )
        M_variable_index_map.insert( make_pair("Q_inert", nvar++ ) );

    //! correction flux with longitudinal term
    if( M_longitudinal_wall )
        M_variable_index_map.insert( make_pair("Q_long", nvar++ ) );

    //! activate the export filters
    // matlab postprocessing

    M_variable_filter_map.insert( make_pair(".m",
                                            &LifeV::OneDModelSolver<Params, Flux, Source>::output_to_matlab) );

    // plotmtv postprocessing
    //    M_variable_filter_map.insert( make_pair(".mtv",
    //  &LifeV::OneDModelSolver<Params, Flux, Source>::output_to_plotmtv) );


    M_U_thistime.resize(nvar, vector_type(M_localMap));

    for (UInt ii = 0; ii < nvar; ++ii)
        M_U_thistime[ii] *= 0.;


    // initialize matrices
    std::fill( M_diffFlux.begin(),    M_diffFlux.end(),    ublas::zero_vector<double>(M_FESpace.dim()) );
    std::fill( M_diffSrc.begin() ,    M_diffSrc.end() ,    ublas::zero_vector<double>(M_FESpace.dim()) );

    M_massMatrixDiffSrc.resize(4);
    M_stiffMatrixDiffFlux.resize(4);
    M_gradMatrixDiffFlux.resize(4);
    M_divMatrixDiffSrc.resize(4);

    for( UInt i = 0; i < 4; ++i )
        {
            // using push_back because of no default constructor available
            //            M_massMatrixDiffSrc.  push_back( new matrix_type( M_localMap ));
            M_massMatrixDiffSrc[i].reset(new matrix_type( M_localMap ));
            //M_stiffMatrixDiffFlux.push_back( new matrix_type( M_localMap ));
            M_stiffMatrixDiffFlux[i].reset( new matrix_type( M_localMap ));
            //M_gradMatrixDiffFlux. push_back( new matrix_type( M_localMap ));
            M_gradMatrixDiffFlux[i].reset( new matrix_type( M_localMap ));
            //M_divMatrixDiffSrc.   push_back( new matrix_type( M_localMap ));
            M_divMatrixDiffSrc[i].reset( new matrix_type( M_localMap ));
        }

    //-------------------------------------------
    //! update first the constant matrices (cst w.r. to time iter)

    //! set the coeff to 1.
    M_coeffMass = 1.;
    M_coeffGrad = 1.;

    //! Elementary computation and matrix assembling
    //! Loop on elements


    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++)
        {
            //! set the elementary matrices to 0.
            M_elmatMass.zero();
            M_elmatGrad.zero();

            //! update the current element

            M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));

            //! update the mass and grad matrices
            mass( M_coeffMass, M_elmatMass, M_FESpace.fe(),0, 0 );
            grad( 0 , - M_coeffGrad, M_elmatGrad, M_FESpace.fe(), M_FESpace.fe(), 0, 0 );
            //! assemble the mass and grad matrices
            assemb_mat( M_massMatrix, M_elmatMass, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
            //assemb_mat( M_factorMassMatrix, M_elmatMass, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
            assemb_mat( M_gradMatrix, M_elmatGrad, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
        } //! end loop on elements

                //M_massMatrix.GlobalAssemble();
    //M_factorMassMatrix = M_massMatrix;


    //! Dirichlet boundary conditions set in the mass matrix
    M_massMatrix.GlobalAssemble();
    M_gradMatrix.GlobalAssemble();

    M_linearSolver.setUpPrec( data_file, section + "/prec");
    M_linearSolver.setDataFromGetPot( data_file, section + "/solver" );

    chrono.stop();

    Debug( 6310 ) << "[OneDModelSolver] \tdone in " << chrono.diff() << " s.\n";


}


/*! Update the coefficients
  (from the flux, source functions and their derivatives)
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateMatrixCoefficients(const UInt& ii,
    const UInt& jj , const UInt& iedge)
{
    Real dFluxdUelem(0), dSrcdUelem(0);

    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    dFluxdUelem = M_diffFlux[ 2*(ii-1) + jj-1 ]( iedge - 1 ); //! iedge starts from 1...
    dSrcdUelem  = M_diffSrc [ 2*(ii-1) + jj-1 ]( iedge - 1 );

    M_coeffGrad  = dFluxdUelem; //! term gradDiffFlux(U) [* S(U)]
    M_coeffDiv   = dSrcdUelem;  //! term  divDiffSrc(U) [* F(U)]
    M_coeffStiff = dFluxdUelem; //! term stiffDiffFlux(U) [* F(U)]
    M_coeffMass  = dSrcdUelem;  //! term  massDiffSrc(U) [* S(U)]
}

//! Update the element matrices with the updated
//! current element and updated coefficients
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateElemMatrices()
{
    //! set the elementary matrices to 0.
    M_elmatMass.zero();
    M_elmatStiff.zero();
    M_elmatGrad.zero();
    M_elmatDiv.zero();

    //! update the mass matrix
    mass( M_coeffMass, M_elmatMass, M_FESpace.fe(),0, 0 );
    //  M_elmatMass.showMe( std::cout );

    //! update the stiffness matrix
    stiff( M_coeffStiff, M_elmatStiff, M_FESpace.fe(),0 ,0 );
    // std::cout << "Elem Stiff matrix :" << std::endl;
    // M_elmatStiff.showMe( std::cout );

    /*! update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0 , - M_coeffGrad, M_elmatGrad, M_FESpace.fe(), M_FESpace.fe(), 0, 0 );
    //  std::cout << "Elem Grad matrix :" << std::endl;
    //  M_elmatGrad.showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this M_elmatDiv is not necessary
      as it is the transpose of the M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0 , - M_coeffDiv, M_elmatDiv, M_FESpace.fe(), M_FESpace.fe(), 0, 0 );
    //  std::cout << "Elem Div matrix :" << std::endl;
    //  M_elmatDiv.showMe( std::cout );
}


//! assemble the matrices
template< class Params, class Flux, class Source >
int
OneDModelSolver<Params, Flux, Source>::
_assemble_matrices(const UInt& ii, const UInt& jj )
{
    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    //! assemble the mass matrix
    assemb_mat( *M_massMatrixDiffSrc[ 2*(ii-1) + jj-1 ]  , M_elmatMass, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

    //! assemble the stiffness matrix
    assemb_mat( *M_stiffMatrixDiffFlux[ 2*(ii-1) + jj-1 ], M_elmatStiff, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

    //! assemble the gradient matrix
    assemb_mat( *M_gradMatrixDiffFlux[ 2*(ii-1) + jj-1 ] , M_elmatGrad, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

    //! assemble the divergence matrix
    assemb_mat( *M_divMatrixDiffSrc[ 2*(ii-1) + jj-1 ]   , M_elmatDiv, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

    return 0;

    ERROR_MSG("Invalid values for the _assemble_matrices method.");
    return 1;
}

/*! update the matrices
  M_massMatrixDiffSrcij, M_stiffMatrixDiffFluxij
  M_gradMatrixDiffFluxij, and M_divMatrixDiffSrcij (i,j=1,2)

  from the values of diffFlux(Un) and diffSrc(Un)
  that are computed with _updateMatrixCoefficients.
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateMatrices()
{
    //--------------------------------------------------------
    // Chrono chrono;
    // chrono.start();
    // std::cout << "o-loop over the matrices INIT... ";

    //! Matrices initialization
    for( UInt i = 0; i < 4; ++i )
        {
//             *M_massMatrixDiffSrc[i]   *= 0.;
//             *M_stiffMatrixDiffFlux[i] *= 0.;
//             *M_gradMatrixDiffFlux[i]  *= 0.;
//             *M_divMatrixDiffSrc[i]    *= 0.;
            M_massMatrixDiffSrc[i].reset(new matrix_type( M_localMap ));
            //M_stiffMatrixDiffFlux.push_back( new matrix_type( M_localMap ));
            M_stiffMatrixDiffFlux[i].reset( new matrix_type( M_localMap ));
            //M_gradMatrixDiffFlux. push_back( new matrix_type( M_localMap ));
            M_gradMatrixDiffFlux[i].reset( new matrix_type( M_localMap ));
            //M_divMatrixDiffSrc.   push_back( new matrix_type( M_localMap ));
            M_divMatrixDiffSrc[i].reset( new matrix_type( M_localMap ));



        }

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices... ";
    */

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++){

        //! update the current element
        M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));
        //  std::cout << M_FESpace.fe().currentId() << std::endl;

        for(UInt ii = 1; ii <= 2; ii ++) {
            for(UInt jj = 1; jj <= 2; jj ++) {

                //! update the M_coeff*
                _updateMatrixCoefficients( ii , jj, iedge);

                //! update the M_elmat*
                _updateElemMatrices();

                //! assemble the global matrices
                _assemble_matrices( ii, jj );
            }
        }
    } //! end loop on elements

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices BC DIR... ";
    */

}


/*! modify the matrix to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::
_updateBCDirichletMatrix( matrix_type& mat )
{
    UInt firstDof = M_leftNodeId;
    UInt lastDof  = M_rightNodeId;

    //! unsymmetric treatment (LU must be used!)
    //! modify the first row

    mat.diagonalize(firstDof - 1, 1, 0);

//     mat.set_mat_inc(firstDof, firstDof    , 1. );
//     mat.set_mat_inc(firstDof, firstDof + 1, 0. );

    //! modify the last row

    mat.diagonalize(lastDof - 1, 1, 0);
//     mat.set_mat_inc(lastDof,  lastDof      , 1.);
//     mat.set_mat_inc(lastDof,  lastDof  - 1 , 0.);

    //! symmetric treatment (cholesky can be used)
    //! modify the first row
//     mat.Diag()( firstDof )    = 1.;
//     mat.UpDiag()( firstDof )  = 0.;
//     mat.LowDiag()( firstDof ) = 0.; //!and second row

//     //! modify the last row
//     mat.Diag()( lastDof )      = 1.;
//     mat.UpDiag()( lastDof-1 )  = 0.; //!and penultimate row
//     mat.LowDiag()( lastDof-1 ) = 0.;
}

/*! modify the vector to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)

  \param val_left  : Dirichlet value inserted to the left
  \param val_right : Dirichlet value inserted to the right
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::
_updateBCDirichletVector()
{
    UInt firstDof = M_leftNodeId;
    UInt lastDof  = M_rightNodeId;

    Debug( 6310 ) << "[_updateBCDirichletVector] \t firstDof = " << firstDof
                  << " lastDof = " << lastDof << ";\n";

    Debug( 6310 ) << "[_updateBCDirichletVector] \t bcDirLeft[0] = " << M_bcDirLeft[0]
                  << " bcDirLeft[1] = " << M_bcDirLeft[1] << ";\n";

    Debug( 6310 ) << "[_updateBCDirichletVector] \t bcDirRight[0] = " << M_bcDirRight[0]
                  << " bcDirRight[1] = " << M_bcDirRight[1] << ";\n";

    //! unsymmetric treatment (LU must be used!)
    //! first row modified
    M_rhs[0]( firstDof ) = M_bcDirLeft[0];
    M_rhs[1]( firstDof ) = M_bcDirLeft[1];

    //! last row modified
    M_rhs[0]( lastDof ) = M_bcDirRight[0];
    M_rhs[1]( lastDof ) = M_bcDirRight[1];


    //! symmetric treatment (cholesky can be used)
//     for(UInt i=0; i<2; ++i) {
//         //! first row modified (Dirichlet)
//         M_rhs[i]( firstDof ) = M_bcDirLeft[i];
//         //! second row modified (for symmetry)
//         M_rhs[i]( firstDof + 1 ) += - M_massMatrix.LowDiag()( firstDof ) * M_bcDirLeft[i];
//         //! last row modified (Dirichlet)
//         M_rhs[i]( lastDof ) = M_bcDirRight[i];
//         //! penultimate row modified (for symmetry)
//         M_rhs[i]( lastDof - 1 ) += - M_massMatrix.UpDiag()( lastDof - 1 ) * M_bcDirRight[i];
    // }

}


//! compute the M_bcDirLeft and M_bcDirRight and set them to the new values
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_computeBC( const Real& time_val )
{
    M_bcH.applyBC(time_val, M_bcDirLeft, M_bcDirRight );
}


//! get the flux function
template< class Params, class Flux, class Source >
const Flux&
OneDModelSolver<Params, Flux, Source>::FluxFun() const
{
    return M_fluxFun;
}


//! get the source function
template< class Params, class Flux, class Source >
const Source&
OneDModelSolver<Params, Flux, Source>::SourceFun() const
{
    return M_sourceFun;
}


// //! get the left edge
// template< class Params, class Flux, class Source >
// Mesh::EdgeType
// OneDModelSolver<Params, Flux, Source>::LeftEdge() const
// {
//     return M_leftEdge;
// }


// //! get the right edge
// template< class Params, class Flux, class Source >
// Mesh::EdgeType
// OneDModelSolver<Params, Flux, Source>::RightEdge() const
// {
//     return M_rightEdge;
// }


//! get the left node
template< class Params, class Flux, class Source >
UInt
OneDModelSolver<Params, Flux, Source>::LeftNodeId() const
{
    return M_leftNodeId;
}


//! get the left internal node (neighboring node)
template< class Params, class Flux, class Source >
UInt
OneDModelSolver<Params, Flux, Source>::LeftInternalNodeId() const
{
    return M_leftInternalNodeId;
}


//! get the right node
template< class Params, class Flux, class Source >
UInt
OneDModelSolver<Params, Flux, Source>::RightNodeId() const
{
    return M_rightNodeId;
}


//

template< class Params, class Flux, class Source >
Real
OneDModelSolver<Params, Flux, Source>::value(std::string var, UInt pos) const
{
    //    M_variable_index_iter iter;

    std::map< std::string, UInt>::const_iterator it = M_variable_index_map.find(var);

    UInt offset = it->second;
    return M_U_thistime[offset](pos);


//     for( iter = M_variable_index_map.begin(); iter != M_variable_index_map.end(); ++iter)
//     {
//         if (iter->first == var)
//         {
//         }
//     }

//     return -1;
}


//! get the right internal node (neighboring node)
template< class Params, class Flux, class Source >
UInt
OneDModelSolver<Params, Flux, Source>::RightInternalNodeId() const
{
    return M_rightInternalNodeId;
}


//! get the Dirichlet boundary conditions (left)
template< class Params, class Flux, class Source >
typename OneDModelSolver<Params, Flux, Source>::Vec2D
OneDModelSolver<Params, Flux, Source>::BCValuesLeft() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )

    temp[0] = M_U_thistime[0]( LeftNodeId() );
    temp[1] = M_U_thistime[1]( LeftNodeId() );

    return temp;
}







//! get the value at neighboring node (left)
template< class Params, class Flux, class Source >
typename OneDModelSolver<Params, Flux, Source>::Vec2D
OneDModelSolver<Params, Flux, Source>::BCValuesInternalLeft() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = M_U_thistime[0]( LeftInternalNodeId() );
    temp[1] = M_U_thistime[1]( LeftInternalNodeId() );
    return temp;
}


//! get the Dirichlet boundary conditions (right)
template< class Params, class Flux, class Source >
typename OneDModelSolver<Params, Flux, Source>::Vec2D
OneDModelSolver<Params, Flux, Source>::BCValuesRight() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = M_U_thistime[0]( RightNodeId() );
    temp[1] = M_U_thistime[1]( RightNodeId() );
    return temp;
}


//! get the value at neighboring node (right)
template< class Params, class Flux, class Source >
typename OneDModelSolver<Params, Flux, Source>::Vec2D
OneDModelSolver<Params, Flux, Source>::BCValuesInternalRight() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = M_U_thistime[0]( RightInternalNodeId() );
    temp[1] = M_U_thistime[1]( RightInternalNodeId() );
    return temp;
}


//! set the Dirichlet boundary conditions (right)
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
{
    M_bcDirRight[0] = bcR1; // M_bcDirRight.first  = bcR1;
    M_bcDirRight[1] = bcR2; // M_bcDirRight.second = bcR2;
}


//! set the Dirichlet boundary conditions (left)
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
{
    M_bcDirLeft[0] = bcL1; // M_bcDirLeft.first  = bcL1;
    M_bcDirLeft[1] = bcL2;
}


template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::setBCLeft_internalnode()
{
    M_bcH.setBCLeft_internalnode();
}


template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::setBCRight_internalnode()
{
    M_bcH.setBCRight_internalnode();
}


//! simple cfl computation (correct for constant mesh)
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::CheckCFL() const
{

    Debug(6310) << "[CheckCFL] checking the CFL ... ";

    Real CFL = 0.;

    //! length of the first edge (arbitrary as they are all supposed equal).
    Real deltaX;

    Real deltaX_min = M_FESpace.mesh()->edgeList( 1 ).point(2).x() - M_FESpace.mesh()->edgeList( 1 ).point(1).x();

    Real Ainode, Qinode;

    Real lambda1_max = 0.;
    Real lambda2_max = 0.;

    Real eigval1, eigval2;
    Real tmp11, tmp12, tmp21, tmp22;

    for ( UInt inode = 1; inode <= M_FESpace.dim() ; inode++ )
        {
            Ainode = M_U_thistime[0]( inode );
            Qinode = M_U_thistime[1]( inode );

            //! compute the eigenvalues at node
            M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
                                                    eigval1, eigval2,
                                                    tmp11, tmp12,
                                                    tmp21, tmp22,
                                                    inode);

            lambda1_max = std::max<Real>( std::fabs(eigval1), lambda1_max );
            lambda2_max = std::max<Real>( std::fabs(eigval2), lambda2_max );

        }

    for ( UInt inode = 1; inode < M_FESpace.dim() ; inode++ )
        {
            deltaX     = M_FESpace.mesh()->edgeList( inode ).point(2).x() - M_FESpace.mesh()->edgeList( inode ).point(1).x();
            deltaX_min = std::min<Real>( std::fabs(deltaX), deltaX_min );
        }


    CFL = M_data.timestep() / deltaX_min * std::max<Real>( lambda1_max , lambda2_max );

    //if ( M_CFL )
    /*
      std::cout << "Old CFL = " << M_data.() /
      (M_FESpace.mesh()->edgeList( 1 ).pt2().x() - M_FESpace.mesh()->edgeList( 1 ).pt1().x())
      * std::max<Real>( lambda1_max , lambda2_max )
      << std::endl;
    */
    //    }

    if( CFL > 0.5774 )
        std::cout << "\n[CheckCFL] CFL not respected in " << M_data.PostFile()
                  << ": CFL = " << CFL << std::endl;

    //ASSERT( CFL < 0.5774 , "CFL not respected" );

    Debug(6310) << "[CheckCFL] ok.";
}


/*! Axpy product for 2D vectors (pairs)
  Axpy(alpha, x, beta, y) -> y = a*A*x + beta*y

  A is given by two pairs corresponding to the 2 lines.
  A = [line1;
  line2 ]
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::Axpy(const Vec2D& line1, const Vec2D& line2,
                                            const Real& alpha,  const Vec2D& x,
                                            const Real& beta,   Vec2D& y) const
{
    ASSERT_PRE( line1.size() == 2 && line2.size(),
                "Axpy works only for 2x2 matrices");
    ASSERT_PRE( x.size() == 2 && y.size(),
                "Axpy works only with 2D vectors");
    y[0] = alpha * ( line1[0] * x[0] + line1[1] * x[1]) + beta * y[0];
    y[1] = alpha * ( line2[0] * x[0] + line2[1] * x[1]) + beta * y[1];
}


//! update the P1 flux vector from U: M_Fluxi = F_h(Un) i=1,2
//! BEWARE: works only for P1Seg elements
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateFlux()
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace.dim() ; ii++ )
        {
            Aii = M_U_thistime[0]( ii );
            Qii = M_U_thistime[1]( ii );
            M_Flux[0]( ii ) = M_fluxFun( Aii, Qii, 1, ii - 1 );
            M_Flux[1]( ii ) = M_fluxFun( Aii, Qii, 2, ii - 1 );
        }
}


/*! call _updateFlux and update the P0 derivative of flux vector from U:
  M_diffFluxij = dF_h/dU(Un) i,j=1,2

  M_diffFluxij(elem) = 1/2 [ dF/dU(U(node1(elem))) + dF/dU(U(node2(elem))) ]

  (mean value of the two extremal values of dF/dU)

  BEWARE: works only for P1Seg elements
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateFluxDer()
{
    //! first update the Flux vector
    _updateFlux();

    //! then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem = 0; ielem < M_FESpace.dim() - 1; ielem++ ) {
        //! for P1Seg and appropriate mesh only!
        ii    = ielem;      //! left node of current element
        iip1  = ielem + 1;  //! right node of current element

        Aii   = M_U_thistime[0]( ielem + 1 );
        Qii   = M_U_thistime[1]( ielem + 1 );

        Aiip1 = M_U_thistime[0]( ielem + 2 );
        Qiip1 = M_U_thistime[1]( ielem + 2 );

        for( UInt ii=1; ii<3; ++ii )
            {
                for( UInt jj=1; jj<3; ++jj )
                    {
                        tmp  = M_fluxFun.diff(   Aii,   Qii, ii, jj, ii );
                        tmp += M_fluxFun.diff( Aiip1, Qiip1, ii, jj, iip1 );

                        M_diffFlux[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5*tmp;
                    }
            }
    }
}

//! update the P1 source vector from U: M_Sourcei = S_h(Un) i=1,2
//! BEWARE: works only for P1Seg elements
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateSource( )
{
    Real Aii, Qii;

    for ( UInt ii = 1; ii <= M_FESpace.dim() ; ii++ )
        {
            Aii = M_U_thistime[0]( ii );
            Qii = M_U_thistime[1]( ii );
            for(UInt k=0; k<2; ++k)
                M_Source[k]( ii ) = M_sourceFun( Aii, Qii, k+1, ii - 1);
    }
}

/*! call _updateSource and update the P0 derivative of source vector from U:
  M_diffSrcij = dS_h/dU(Un) i,j=1,2

  M_diffSrcij(elem) = 1/2 [ dS/dU(U(node1(elem))) + dS/dU(U(node2(elem))) ]

  (mean value of the two extremal values of dS/dU)

  BEWARE: works only for P1Seg elements
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::_updateSourceDer( )
{
    //! first update the Source vector
    _updateSource();

    //! then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < M_FESpace.dim() - 1 ; ielem++ )
        {
            //! for P1Seg and appropriate mesh only!
            ii = ielem;        //! left node of current element
            iip1 = ielem + 1;  //! right node of current element

            Aii   = M_U_thistime[0]( ielem + 1);
            Qii   = M_U_thistime[1]( ielem + 1);
            Aiip1 = M_U_thistime[0]( ielem + 2 );
            Qiip1 = M_U_thistime[1]( ielem + 2 );

            for( UInt ii=1; ii<3; ++ii )
                {
                    for( UInt jj=1; jj<3; ++jj )
                        {
                            tmp =  M_sourceFun.diff(   Aii,   Qii, ii, jj, ii );
                            tmp += M_sourceFun.diff( Aiip1, Qiip1, ii, jj, iip1 );
                            M_diffSrc[ 2*(ii - 1) + jj - 1 ]( ielem ) = 0.5 * tmp;
                    }
                }
        }
}

/*
  We use a finite element scheme for the correction term:
  given the solution of Taylor-Galerkin scheme, solve
  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
  ( m / rho ) * ( dQtildeh(n+1)/dz, dphi/dz )  //! stiff * Qtilde(U)
  = ( m / rho ) *       ( dQhath(n+1)/dz, dphi/dz )  //! stiff * Qhat(U)
  ---------------------------------------------------
  m = rho_w h0 / ( 2 sqrt(pi) sqrt(A0) )
*/
template< class Params, class Flux, class Source >
vector_type
OneDModelSolver<Params, Flux, Source>::_correct_flux_inertial( const vector_type& flux )
{
    matrix_type _matrixLHS(M_localMap);
    matrix_type _stiffRHS (M_localMap);

    ElemMat _elmatMassLHS  (M_FESpace.fe().nbNode, 1, 1);
    ElemMat _elmatStiffLHS (M_FESpace.fe().nbNode, 1, 1);
    ElemMat _elmatStiffRHS (M_FESpace.fe().nbNode, 1, 1);

    vector_type _rhs(M_localMap);

    Real _coeffMass;
    Real _coeffStiff;

    Real m, meanA0;

    //    std::ostringstream output;

    _matrixLHS *= 0.;
     _stiffRHS  *= 0.;

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS. zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        _coeffMass  = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass /= 2;
        _coeffMass  = 1./_coeffMass;

        meanA0  = M_oneDParams.Area0(iedge - 1) + M_oneDParams.Area0(iedge);
        meanA0 /= 2;

        m = M_oneDParams.DensityWall()*M_oneDParams.Thickness()/
            ( 2*std::sqrt(4*std::atan(1))*std::sqrt(meanA0) );

        _coeffStiff = m/M_oneDParams.DensityRho();

        //! update the current element
        M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));

        mass (   _coeffMass,  _elmatMassLHS,  M_FESpace.fe(), 0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace.fe(), 0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace.fe(), 0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

        Debug( 6310 ) << "\n\tm = "           << m
                      << "\n\t_coeffMass = "  << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } //! end loop on elements

    // update rhs
    //_stiffRHS.Axpy( 1., flux , 0., _rhs );

    _rhs = _stiffRHS*flux;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    //  _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    //  _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    //_tridiagsolver.Factor( _matrixLHS );

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    //_tridiagsolver.Solve( _matrixLHS, _rhs );

    vector_type _sol(_rhs);

    M_linearSolver.setMatrix(_matrixLHS);
    //@int numIter = M_linearSolver.solveSystem( _rhs, _sol, _matrixLHS, true);

    //std::cout <<" iterations number :  " << numIter << std::endl;

    return _sol;
}

/*
  L2 Projection of the second derivative of Q over P1 space
*/
template< class Params, class Flux, class Source >
ScalVec
OneDModelSolver<Params, Flux, Source>::_compute_d2Q_dx2( const ScalVec& flux )
{
    matrix_type _massLHS (M_FESpace.dim());
    matrix_type _stiffRHS(M_FESpace.dim());

    TriDiagCholesky< Real, matrix_type, Vector > _tridiagsolver(M_FESpace.dim());

    ElemMat _elmatMassLHS (M_FESpace.fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace.fe().nbNode,1,1);

    ScalVec _rhs(M_FESpace.dim());

    _massLHS  *= 0.;
    _stiffRHS *= 0.;

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS *= 0.;
        _elmatStiffRHS *= 0.;

        //! update the current element
        M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));
        //std::cout << M_FESpace.fe().currentId() << std::endl;

        mass( 1., _elmatMassLHS, M_FESpace.fe(),0, 0 );
        stiff( -1., _elmatStiffRHS, M_FESpace.fe(),0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

    } //! end loop on elements

    // update rhs
    _rhs = _stiffRHS*flux;

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.Diag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.UpDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.LowDiag():\n";
//     for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
//         Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
//     Debug( 6310 ) << "\n";

//     _tridiagsolver.Factor( _massLHS );

//     Debug( 6310 ) << "\n[_compute_d2Q_dx2] solving with rhs:\n";
//     for(int i=0; i<_massLHS.OrderMatrix(); ++i)
//         Debug( 6310 ) << "\t" << _rhs(i);
//     Debug( 6310 ) << "\n";

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}


/*
  We use a finite element scheme for the correction term:
  given the solution of Taylor-Galerkin scheme, solve
  ( 1/(dt*Ah(n+1)) Qtildeh(n+1), phi) +        //! 1/A * massFactor^{-1} * Un+1
  ( gamma / rho ) *     ( dQtildeh(n+1)/dz, dphi/dz )  //! stiff * Qtilde(U)
  = - ( gamma / rho ) * ( dQhath(n+1)/dz, dphi/dz )  //! stiff * Qhat(U)
  ---------------------------------------------------
  gamma = gamma_tilde / ( 2 sqrt(pi) )
*/
template< class Params, class Flux, class Source >
vector_type
OneDModelSolver<Params, Flux, Source>::_correct_flux_viscoelastic( const vector_type& flux )
{
    matrix_type _matrixLHS(M_localMap);
    matrix_type _stiffRHS (M_localMap);

    //TriDiagCholesky< Real, matrix_type, Vector > _tridiagsolver(M_FESpace.dim());

    ElemMat _elmatMassLHS  (M_FESpace.fe().nbNode,1,1);
    ElemMat _elmatStiffLHS (M_FESpace.fe().nbNode,1,1);
    ElemMat _elmatStiffRHS (M_FESpace.fe().nbNode,1,1);

    vector_type _rhs(M_FESpace.map());

    Real _coeffMass;
    Real _coeffStiff;

    Real gamma, meanA0(1.);

    _matrixLHS *= 0.;
    _stiffRHS  *= 0.;

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        // this comes from the exact derivation of generalized string model
        // + Voigt viscoelasticity
        _coeffMass = M_rhs[0]( iedge - 1 ) + M_rhs[0]( iedge );
        _coeffMass *= 0.5;
        _coeffMass = 1./std::sqrt(_coeffMass);

        if(M_linearize_string_model) {
            // this is to recover the linearized version (_coeffMass = 1/A)
            _coeffMass *= _coeffMass;

            meanA0  = M_oneDParams.Area0( iedge - 1 ) + M_oneDParams.Area0( iedge );
            meanA0 *= 0.5;

            if(M_linearize_equations) {
                // when using linearized equations, A \simeq A0
                _coeffMass = 1./ meanA0;
            }
        }
        gamma = M_oneDParams.Gamma() / ( 2 * std::sqrt(4*std::atan(1)) );
        gamma *= 1 / std::sqrt( meanA0 );
        _coeffStiff = M_data.timestep() * gamma / M_oneDParams.DensityRho();

        //! update the current element
        M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));
        //std::cout << M_FESpace.fe().currentId() << std::endl;

        mass (   _coeffMass,   _elmatMassLHS, M_FESpace.fe(),0, 0 );
        stiff(   _coeffStiff, _elmatStiffLHS, M_FESpace.fe(),0, 0 );
        stiff( - _coeffStiff, _elmatStiffRHS, M_FESpace.fe(),0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

        Debug( 6310 ) << "\n\tgamma = " << gamma
                      << "\n\t_coeffMass = " << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } //! end loop on elements

    // update rhs
    // rhs = _stiffRHS * rhs
    _rhs = _stiffRHS*flux;

    // NOTE I should add to rhs the value of boundary integral ("neumann-like")
    // BUT it's useless since I'm going to impose dirichlet conditions on all boundaries

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    vector_type _sol(_rhs);

    M_linearSolver.setMatrix(_matrixLHS);
    //int numIter = M_linearSolver.solveSystem( _rhs, _sol, _matrixLHS, true);

    return _sol;
}


/*
  We use a finite element scheme for the correction term:
  given the solution of Taylor-Galerkin scheme, solve
  ( 1/Ah(n+1) Qtildeh(n+1), phi) +             //! 1/A * massFactor^{-1} * Un+1
  = ( 1/Ah(n+1) Qtildeh(n), phi) +             //! 1/A * massFactor^{-1} * Un+1
  + ( a / rho ) *       ( d3Ahath(n+1)/dz3, phi )  //! mass * d3Ahat(U)/dz
  ---------------------------------------------------
  a = ?
*/
template< class Params, class Flux, class Source >
vector_type
OneDModelSolver<Params, Flux, Source>::_correct_flux_longitudinal( )
{
    matrix_type _massLHS(M_localMap);
    matrix_type _massRHS(M_localMap);

    //TriDiagCholesky< Real, matrix_type, Vector > _tridiagsolver(M_FESpace.dim());

    ElemMat _elmatMassLHS (M_FESpace.fe().nbNode,1,1);
    ElemMat _elmatMassRHS (M_FESpace.fe().nbNode,1,1);

    ScalVec _rhs(M_FESpace.dim());
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    ScalVec _g(M_FESpace.dim());
    ScalVec _f(M_FESpace.dim());

    //          _g = M_rhs[0];
    for( UInt i=0; i<M_FESpace.dim(); ++i )
        _g(i) = std::sqrt(M_rhs[0](i)) - std::sqrt(M_oneDParams.Area0(i));

    UInt inode;

    Real _coeffMassLHS;
    Real _coeffMassRHS;

    Real _a;
    //    std::ostringstream output;
    _massLHS *= 0.;
    _massRHS *= 0.;

    Real _h( M_oneDParams.Length() / static_cast<Real>(M_oneDParams.ParamSize() - 1) );

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= M_FESpace.mesh()->numEdges(); iedge++){

        inode = iedge - 1;

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        _coeffMassLHS = M_rhs[0]( inode ) + M_rhs[0]( inode+1 );
        _coeffMassLHS /= 2;
        _coeffMassLHS = 1./_coeffMassLHS;

        _a = M_oneDParams.CoeffA() / std::sqrt(4*std::atan(1));
        _coeffMassRHS = M_data.timestep() * _a / M_oneDParams.DensityRho();

        // backward differentiation when near to the left boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi) + A(xi+3) - 3A(xi+2) + 3A(xi+1) )

        // forward differentiation when near to the right boundary
        // d3 A (xi)/ dz3 = (1/h^3) * ( A(xi) - A(xi-3) + 3A(xi-2) - 3A(xi-1) )

        // central differentiation otherwise
        // d3 A (xi)/ dz3 = (1/h^3) * ( -A(xi-2) + 2A(xi-1) - 2A(xi+1) + A(xi+2) )

        Debug( 6310 ) << "\ninode = " << inode << "\n";
        if(inode<2) { // backward differentiation
            _f( inode ) = -_g( inode ) + _g( inode+3 )
                - 3*_g( inode+2 ) + 3*_g( inode+1 );
            Debug( 6310 ) << "\n\tbackward differentiation = " << _coeffMassLHS << "\n"; }
        else if(inode>M_FESpace.mesh()->numEdges()-2) { // forward differentiation
            _f( inode ) = _g( inode ) - _g( inode-3 )
                + 3*_g( inode-2 ) - 3*_g( inode-1 );
            Debug( 6310 ) << "\n\forward differentiation = " << _coeffMassLHS << "\n"; }
        else  { // central differentiation
            _f( inode ) = -_g( inode - 2 ) + 2*_g( inode-1 )
                - 2*_g( inode+1 ) + _g( inode+2 );
            Debug( 6310 ) << "\n\tcentral differentiation = " << _coeffMassLHS << "\n"; }

        _f(inode) *= 1/(2*std::pow(_h,3));

        //! update the current element
        M_FESpace.fe().updateFirstDerivQuadPt(M_FESpace.mesh()->edgeList(iedge));

        mass( _coeffMassLHS, _elmatMassLHS, M_FESpace.fe(),0, 0 );
        mass( _coeffMassRHS, _elmatMassRHS, M_FESpace.fe(),0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );
        assemb_mat( _massRHS, _elmatMassRHS, M_FESpace.fe(), M_FESpace.dof() , 0, 0 );

        Debug( 6310 ) << "\n\t_coeffMassLHS = " << _coeffMassLHS << "\n";
        Debug( 6310 ) << "\n\t_coeffMassRHS = " << _coeffMassRHS << "\n";

    } //! end loop on elements

    // update rhs
    //    _massLHS.Axpy( 1., M_U_thistime[4] , 0., _rhs );
    _rhs = _massRHS*_f;

    UInt firstDof = 1;
    UInt lastDof  = _rhs.size();

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

    //@_tridiagsolver.Factor( _massLHS );

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    //@_tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}


//! Initialize from Reimann invariants
//! Initialize with constant initial conditions
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const Real& u10, const Real& u20,
                                                  const std::string& var )
{

    Debug( 6310 ) << "[OneDModelSolver::initialize] O- Initialize: var " << var << "\n";

    if( var == "physical")
      {
          //
          Debug( 6310 ) << "[OneDModelSolver::initialize] O- Imposing real values ... \n";
          Debug( 6310 ) << "[OneDModelSolver::initialize] O- A0 = " << u10 << "\n";
          Debug( 6310 ) << "[OneDModelSolver::initialize] O- Q0 = " << u20 << "\n";

          M_U_thistime[0] = vector_type( M_localMap );
          //M_U_thistime[0][LeftNodeId()] = u10;


          M_U_thistime[1] = vector_type( M_localMap );
          //M_U_thistime[1][LeftNodeId()] = u20;

//           std::cout << "LeftNodeId() " << LeftNodeId() << std::endl;

//           std::cout << "A0 = " << M_U_thistime[0][1] << " " << BCValuesLeft()[0] << std::endl;
//           std::cout << "Q0 = " << M_U_thistime[0][2] << " " << BCValuesLeft()[1] << std::endl;

          for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
          {
              M_U_thistime[0][ielem + 1] = u10;
              M_U_thistime[1][ielem + 1] = u20;

              //              std::cout << ielem << " " << M_U_thistime[0][ielem + 1] << std::endl;
              M_oneDParams.W_from_U( M_U_thistime[2][ielem + 1], M_U_thistime[3][ielem + 1],
                                     M_U_thistime[0][ielem + 1], M_U_thistime[1][ielem + 1],
                                     ielem );
          }

          Debug( 6310 ) << "[OneDModelSolver::initialize] O- ok\n";
      }
    else if( var == "reimann" )
        {
            //        for( UInt i=0; i<2; ++i )
            M_U_thistime[2] = vector_type( M_localMap );
            M_U_thistime[2] = u10;
            M_U_thistime[3] = vector_type( M_localMap );
            M_U_thistime[3] = u20;

            for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
                M_oneDParams.U_from_W( M_U_thistime[0][ielem + 1], M_U_thistime[1][ielem + 1],
                                       M_U_thistime[2][ielem + 1], M_U_thistime[3][ielem + 1],
                                       ielem + 1 ); // WARNING the +1 is not debugged yet (GF 12/2009)
        }
    else
        {
            std::cout << "[initialize] trying to initialize " << var << " variables!" << std::endl;
            abort();
        }


    for( UInt i = 0; i < 2; ++i )
        {
            M_U_prevtime [i] = M_U_thistime[i];
            M_U_2prevtime[i] = M_U_prevtime[i];
        }

    ScalVec pressures(4*M_FESpace.dim());

    for (UInt ielem = 0; ielem < M_FESpace.dim() ; ++ielem )
      {
          subrange(pressures, 4*ielem, 4+4*ielem) =
              M_oneDParams.pressure( M_U_thistime [0][ielem + 1],
                                     M_U_prevtime [0][ielem + 1],
                                     M_U_2prevtime[0][ielem + 1],
                                     M_data.timestep(),
                                     ielem,
                                     M_dP_dt_steps,
                                     M_viscoelastic_wall,
                                     M_linearize_string_model );
          M_U_thistime[4][ielem + 1] = pressures(4*ielem);
      }

    if(M_viscoelastic_wall)
      {
        for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
	  {
	    M_U_thistime[M_variable_index_map.find("P_elast")->second][ielem + 1]
	      = pressures(1+4*ielem);
            M_U_thistime[M_variable_index_map.find("P_visc")->second][ielem + 1]
	      = pressures(2+4*ielem);
            M_U_thistime[M_variable_index_map.find("dA_dt")->second][ielem + 1]
	      = pressures(3+4*ielem);
        }
    }

    //! Prepare bc Handler
    M_bcH.setDefaultBC( M_FESpace, M_sourceFun, M_data.timestep() );


    //! create matlab scripts
    create_movie_file();

    openFileBuffers();

    output2FileBuffers( 0. );

    //postProcess( 0. );

}

//! Initialize with vectors containing all the nodal values
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const Vector& u10, const Vector& u20)
{

    M_U_thistime[0] = u10;
    M_U_thistime[1] = u20;

    for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ )
    {
//         M_U_thistime[0][ielem] = u10;
//         M_U_thistime[1][ielem] = u20;

        M_oneDParams.W_from_U( M_U_thistime[2][ielem], M_U_thistime[3][ielem],
                               M_U_thistime[0][ielem], M_U_thistime[1][ielem], ielem );
    }

    for( UInt i=0; i<2; ++i )
        {
            M_U_prevtime[i] = M_U_thistime[i];
            M_U_2prevtime[i] = M_U_prevtime[i];
        }

    Vector pressures(4 * M_FESpace.dim());

    for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ ) {
        subrange(pressures, 4*ielem, 4+4*ielem) =
            M_oneDParams.pressure( M_U_thistime[0][ielem],
                                   M_U_prevtime[0][ielem], M_U_2prevtime[0][ielem],
                                   M_data.timestep(), ielem,
                                   M_dP_dt_steps, M_viscoelastic_wall,
                                   M_linearize_string_model );

        M_U_thistime[4][ielem] = pressures(4*ielem);
    }

    if(M_viscoelastic_wall) {
        for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ ) {
            M_U_thistime[M_variable_index_map.find("P_elast")->second][ielem] =
              pressures(1+4*ielem);
            M_U_thistime[M_variable_index_map.find("P_visc")->second][ielem] =
              pressures(2+4*ielem);
            M_U_thistime[M_variable_index_map.find("dA_dt")->second][ielem] =
              pressures(3+4*ielem);
        }
    }

    //! Prepare bc Handler
    M_bcH.setDefaultBC( M_FESpace, M_sourceFun, M_data.timestep() );

    //! create matlab scripts
    //create_movie_file();

    //openFileBuffers();

    //output2FileBuffers( 0. );

    postProcess( 0. );

}

//! Initialize when initial conditions concentration
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const Function& /*c0*/,
    Real /*t0*/, Real /*dt*/)
{
    ERROR_MSG("Not yet implemented");
}

// ! Initialize when initial values for the concentration are read from file
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const std::string & vname)
{
    ERROR_MSG("Not yet implemented");
    /*
    //! create matlab scripts
    create_movie_file();
    */
    std::fstream Resfile(vname.c_str(),std::ios::in | std::ios::binary);
    if (Resfile.fail()) {
        std::cerr<<" Error in initialize: File not found or locked"<<std::endl;
        abort();
    }
    Resfile.read((char*)&M_U_thistime[0](1),M_U_thistime[0].size()*sizeof(Real));
    Resfile.close();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

}

//! Initialize only Flux (Area read from OneDNonLinParam)
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const Real& u20)
{
    M_U_thistime[1] = vector_type( M_localMap );
    M_U_thistime[1] = u20;


    ScalVec pressures(4*M_dimDof);

    for (UInt ielem = 0; ielem <= M_FESpace.dim() ; ielem++ )
        {
            M_U_thistime[0][ielem]=M_oneDParams.Area0(ielem);
            for( UInt i=0; i<2; ++i )
                {
                    M_U_prevtime [i][ielem]=M_U_thistime[i][ielem];
                    M_U_2prevtime[i][ielem]=M_U_prevtime[i][ielem];
                }
            M_oneDParams.W_from_U( M_U_thistime[2][ielem], M_U_thistime[3][ielem],
                                   M_U_thistime[0][ielem], M_U_thistime[1][ielem], ielem );

            //      for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ ) {
            subrange(pressures, 4*ielem, 4+4*ielem) =
                M_oneDParams.pressure( M_U_thistime[0][ielem],
                                       M_U_prevtime[0][ielem], M_U_2prevtime[0][ielem],
                                       M_data.timestep(), ielem,
                                       M_dP_dt_steps, M_viscoelastic_wall,
                                       M_linearize_string_model );

            M_U_thistime[4][ielem] = pressures(4*ielem);
        }

    if(M_viscoelastic_wall)
        {
            for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ )
                {
                    M_U_thistime[M_variable_index_map.find("P_elast")->second][ielem] =
                        pressures(1+4*ielem);
                    M_U_thistime[M_variable_index_map.find("P_visc")->second][ielem] =
                        pressures(2+4*ielem);
                    M_U_thistime[M_variable_index_map.find("dA_dt")->second][ielem] =
                        pressures(3+4*ielem);
                }
        }

    //! Prepare bc Handler
    M_bcH.setDefaultBC( M_FESpace.mesh(), M_sourceFun, M_data.timestep() );

    //! create matlab scripts
    create_movie_file();

    openFileBuffers();

    output2FileBuffers( 0. );

    //postProcess( 0. );

}


/*
  Initialize variables with non constant values in space.
  A step profile is approximated by smoothing the edges with a gaussian "bell".
  Read options from data file: section [initialize]
  var = { A, Q, P, W1, W2} the variable to be initialized
  rest_value = rest value
  firstnode =,
  lastnode = domain of the step
  multiplier = step value (= multiplier * value)
  amplitude = width of the gaussian bell
  X = rest_value ( 1 + multiplier * exp( - ( z - node ) / ( 2 * width^2 ) ) )
  where node is alternatively firstnode or lastnode
*/
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::initialize(const GetPot& data_file, const std::string section)
{
    // the discontinuity is comprised between firstnode and lastnode
    UInt firstnode( data_file((section+"/initialize/firstnode").data(),1) );
    UInt lastnode( data_file((section+"/initialize/lastnode").data(),2) );

    ASSERT_PRE( (firstnode <= lastnode) && (lastnode <= M_FESpace.dim()),
                "[initialize] outside tube boundaries" );

    // read initialization type from data file (see OneDModelSolver::initialize)
    std::string init_var( data_file((section+"/initialize/var").data(),"P") );
    Real multiplier( data_file((section+"/initialize/multiplier").data(),1.) );

    std::cout << "init_var = " << init_var << std::endl;
    // tell me what I am doing
    Debug( 6310 ) << "[initialize] 0- Initializing with values:\n";
    Debug( 6310 ) << "[initialize]\t\tinitialize var = " << init_var << "\n";
    Debug( 6310 ) << "[initialize]\t\tfirstnode      = " << firstnode << "\n";
    Debug( 6310 ) << "[initialize]\t\tlastnode       = " << lastnode << "\n";
    Debug( 6310 ) << "[initialize]\t\tmultiplier     = " << multiplier << "\n";

    Real value1, value2, width( data_file((section+"/initialize/width").data(),5.) );

    Vector exponent(M_FESpace.dim());


    //M_localMap.getMap(Unique)->Print(std::cout);

    for (UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
        {
            exponent[inode - 1] *= 0.;

            // first half of a gaussian signal, centered in firstnode;
            // second half of a gaussian signal, centered in lastnode;
            // width represents the total duration of the gaussian signal
            // (rise + decay)
            if( (inode < firstnode) || (inode > lastnode) ) {
                exponent[inode - 1]  = - std::pow( double(int( inode - firstnode )), 2 );
                exponent[inode - 1] /= 2*std::pow( width, 2 );
            }
        }


    switch( M_oneDstring2initializeVarMap[init_var] )
        {
            // case 1, 2: initialize physical variables to desired value
        case OneDInitPressure:
            //std::cout << "OneDInitPressure" << std::endl;
            // this is a pressure value! has to be converted in area value
            Debug( 6310 ) << "[initialize] 0- OneDInitPressure\n";

            value1 = data_file((section+"/initialize/rest_value").data(),0.);
            // HYPOTHESIS: when initializing pressure, flux is imposed constant = 0
            value2 = 0;

            Debug( 6310 ) << "[initialize] pressure " << value1 << "\n";

            //ScalarVector( M_FESpace.dim(), value2 );
            M_U_thistime[1] = vector_type(M_localMap);
            M_U_thistime[1] = value2;

            Debug( 6310 ) << "[initialize] Q done\n";

            for (UInt inode = M_leftNodeId; inode <= M_rightNodeId ; ++inode )
                {
                    // reusing value2 as help variable
                    value2 = value1*( 1 + multiplier*std::exp( exponent[inode - 1] ) );
                    M_U_thistime[0][inode] = M_oneDParams.A_from_P( value2 );

                    Debug( 6310 ) << "[initialize] A(" << inode <<") done\n";
                    M_oneDParams.W_from_U( M_U_thistime[2][inode], M_U_thistime[3][inode],
                                           M_U_thistime[0][inode], M_U_thistime[1][inode],
                                           inode - 1);
                    Debug( 6310 ) << "[initialize] W_i(" << inode <<") done\n";
                }
            break;

        case OneDInitArea:
            Debug( 6310 ) << "[initialize] 0- OneDInitArea\n";

            //std::cout << "OneDInitArea" << std::endl;
            value1 = data_file((section+"/initialize/value").data(), 0.);
            value2 = 0;

            //ScalarVector( M_U_thistime[1].size(), value2 );
            //M_U_thistime[0] = vector_type(M_localMap);
            M_U_thistime[1] = value2;


            for (UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
                {
                    M_U_thistime[0][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode - 1] ) );

                    M_oneDParams.W_from_U( M_U_thistime[2][inode], M_U_thistime[3][inode],
                                           M_U_thistime[0][inode], M_U_thistime[1][inode],
                                           inode - 1 );
                }
            break;

        case OneDInitFlux:
            // HYPOTHESIS: when initializing flux, area is equal to Area0
            Debug( 6310 ) << "[initialize] 0- OneDInitFlux\n";


            value1 = M_oneDParams.Area0(0); // this if Area0 is constant
            value2 = data_file((section+"/initialize/value").data(), 0.);

            for (UInt inode = M_leftNodeId; inode <= M_rightNodeId ; ++inode )
                {
                    M_U_thistime[0][inode] = M_oneDParams.Area0(inode - 1);
                    M_U_thistime[1][inode] = value2*( 1 + multiplier*std::exp( exponent[inode - 1] ) );

                    M_oneDParams.W_from_U( M_U_thistime[2][inode], M_U_thistime[3][inode],
                                           M_U_thistime[0][inode], M_U_thistime[1][inode],
                                           inode - 1 );
                }
            break;

        case OneDInitReimann1:
            Debug( 6310 ) << "[initialize] 0- OneDInitReimann1\n";
            value1 = data_file((section+"initialize/value").data(), 0.);
            value2 = -value1;
            std::cout << "[initialize] WARNING! Initializing W2 = - W1"
                      << " (assuming Q = 0)" << std::endl;

            for (UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
                {
                    M_U_thistime[2][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode - 1] ) );
                    M_U_thistime[3][inode] = value2 *
                      ( 1 + multiplier * std::exp( exponent[inode - 1] ) );
                    M_oneDParams.U_from_W( M_U_thistime[0][inode],
                                           M_U_thistime[1][inode],
                                           M_U_thistime[2][inode],
                                           M_U_thistime[3][inode],
                                           inode - 1);
                }

            break;

        case OneDInitReimann2:
            Debug( 6310 ) << "[initialize] 0- OneDInitReimann2\n";
            value1 = data_file("initialize/value",0.);
            value2 = - value1;
            std::cout << "[initialize] WARNING! Initializing W1 = - W2"
                      << " (assuming Q = 0)" << std::endl;

            for (UInt inode = M_leftNodeId; inode <= M_rightNodeId; ++inode )
                {
                    M_U_thistime[3][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode - 1] ) );
                    M_U_thistime[2][inode] = value2 *
                      ( 1 + multiplier * std::exp( exponent[inode - 1] ) );
                    M_oneDParams.U_from_W( M_U_thistime[0][inode],
                                           M_U_thistime[1][inode],
                                           M_U_thistime[2][inode],
                                           M_U_thistime[3][inode],
                                           inode );
                }


            break;

        default:
            ERROR_MSG("No such initializing option.");

        }


    Debug( 6310 ) << "[initialize]\t\tvalue1         = " << value1 << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue1_step    = " << value1 * multiplier << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue2         = " << value2 << "\n";

    for( UInt i = 0; i < 2; ++i )
        {
            M_U_prevtime [i] = M_U_thistime[i];
            M_U_2prevtime[i] = M_U_prevtime[i];
        }

    Vector pressures(4*M_FESpace.dim());


    for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
        {
            subrange(pressures, 4*ielem, 4 + 4*ielem) =
                M_oneDParams.pressure( M_U_thistime [0][ielem + 1],
                                       M_U_prevtime [0][ielem + 1],
                                       M_U_2prevtime[0][ielem + 1],
                                       M_data.timestep(),
                                       ielem,
                                       M_dP_dt_steps,
                                       M_viscoelastic_wall,
                                       M_linearize_string_model );

            M_U_thistime[4][ielem + 1] = pressures(4*ielem);
        }


    if(M_viscoelastic_wall)
        {
            for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
                {
                    M_U_thistime[M_variable_index_map.find("P_elast")->second][ielem + 1] =
                        pressures(1 + 4*ielem);
                    M_U_thistime[M_variable_index_map.find("P_visc")->second][ielem + 1]  =
                        pressures(2 + 4*ielem);
                    M_U_thistime[M_variable_index_map.find("dA_dt")->second][ielem + 1]   =
                        pressures(3 + 4*ielem);
                }
        }


    //! Prepare bc Handler
    M_bcH.setDefaultBC( M_FESpace, M_sourceFun, M_data.timestep() );

    //! resetting the file buffers

    //resetFileBuffers();

    //! create matlab scripts
    //create_movie_file();

    //! Prepare the buffers
    openFileBuffers();

    //! Write down initial condition
    //output2FileBuffers( 0. );
    postProcess( 0. );


}


// Remember the values at previous timestep
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::savesol()
{
    //    for( UInt i=0; i<2; ++i )
    M_U_prevtime[0] = M_U_thistime[0];
    M_U_prevtime[1] = M_U_thistime[1];
}

// Recover the solution at previous timestep
// BUT keep unaltered the boundary values
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::loadsol()
{
    Vec2D U_leftbd( 2 ), U_rightbd( 2 );
    for( UInt i=0; i<2; ++i ) {
        U_leftbd[i] = M_U_thistime[i]( M_leftNodeId );
        U_rightbd[i] = M_U_thistime[i]( M_rightNodeId );
    }

    //    for( UInt i=0; i<2; ++i )
    M_U_thistime[0] = M_U_prevtime[0];
    M_U_thistime[1] = M_U_prevtime[1];

    for( UInt i=0; i<2; ++i ) {
        M_U_thistime[i]( M_leftNodeId ) = U_leftbd[i];
        M_U_thistime[i]( M_rightNodeId ) = U_rightbd[i];
    }

    for (UInt inode=M_leftNodeId; inode <= M_rightNodeId ; ++inode )
        M_oneDParams.W_from_U( M_U_thistime[2][inode], M_U_thistime[3][inode],
                               M_U_thistime[0][inode], M_U_thistime[1][inode],
                               inode - 1 );

}


//! Update the right hand side for time advancing
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::timeAdvance( const Real& time_val )
{
    Chrono chrono;
    Chrono chronoall;

    chronoall.start();

    Real dt2over2 = M_data.timestep() * M_data.timestep() * 0.5;

    Debug( 6310 ) << "[timeAdvance]\t o-  updates of flux and sources... " << "\n";
    chrono.start();

    //! output cfl
    CheckCFL();
    chrono.stop();
    Debug( 6310 ) << "[timeAdvance]\t\t CFL checked in    " << chrono.diff() << " s.\n";
    //! update the vector containing the values of the flux at the nodes
    //! and its jacobian

    chrono.start();
    _updateFluxDer();
    chrono.stop();

    Debug( 6310 ) << "[timeAdvance]\t\t Flux updated in   " << chrono.diff() << " s.\n";


    //! update the vector containing the values of the source term at the nodes
    //! and its jacobian
    chrono.start();
    _updateSourceDer();
    chrono.stop();

    Debug( 6310 ) << "[timeAdvance]\t\t Source updated in " << chrono.diff() << " s.\n";
    Debug( 6310 ) << "[timeAdvance]\t o-  updates of matrices... ";

    chrono.start();
    //! update the matrices for the non-linear terms
    _updateMatrices();
    chrono.stop();

    Debug( 6310 ) << "done in " << chrono.diff() << " s.\n";

    /*!
      ---------------------------------------------------
      Taylor-Galerkin scheme: (explicit, U = [U1,U2]^T )
      ---------------------------------------------------

      (Un+1, phi) =                             //! massFactor^{-1} * Un+1
      (Un, phi)                         //!            mass * U
      + dt     * (       Fh(Un), dphi/dz )         //!            grad * F(U)
      - dt^2/2 * (diffFh(Un) Sh(Un), dphi/dz )     //! gradDiffFlux(U) * S(U)
      + dt^2/2 * (diffSh(Un) dFh/dz(Un), phi )     //!   divDiffSrc(U) * F(U)
      - dt^2/2 * (diffFh(Un) dFh/dz(Un), dphi/dz ) //!stiffDiffFlux(U) * F(U)
      - dt     * (       Sh(Un), phi )             //!            mass * S(U)
      + dt^2/2 * (diffSh(Un) Sh(Un), phi )         //!  massDiffSrc(U) * S(U)
      ---------------------------------------------------
    */

    Debug( 6310 ) << "[timeAdvance]\t o-  Matrix vector products... " << "\n";
    chrono.start();

    //! Reminder of the function Axpy:
    //! Axpy(alpha, x, beta, y) -> y = alpha*A*x + beta*y

    //!--------------------------------------------------------
    //! 1/, 2/ compute M_rhs[0], M_rhs[1] (systems in U1=A, U2=Q)
    //!--------------------------------------------------------

    //M_gradMatrix.GlobalAssemble();

    for( UInt i = 0; i < 2; ++i )
        {
            //! rhs = mass * Un
            //M_massMatrix.Axpy( 1., M_U_thistime[i] , 0., M_rhs[i] );
            M_rhs[i] = M_massMatrix*M_U_thistime[i];
            //M_gradMatrix.Axpy( M_data.timestep(), M_Flux[i] , 1., M_rhs[i] );
            M_rhs[i] += M_gradMatrix*(M_data.timestep()*M_Flux[i]);


            //! rhs = rhs - dt * mass * S(Un)
            M_rhs[i] += M_massMatrix*(- M_data.timestep()*M_Source[i]);

            for( UInt j=0; j<2; ++j )
                {
                    //! rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
                    M_gradMatrixDiffFlux[2*i + j]->GlobalAssemble();
                    M_rhs[i] += *M_gradMatrixDiffFlux[2*i + j]*(- dt2over2*M_Source[j]);

                    //! rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
                    //M_divMatrixDiffSrc[2*i+j].Axpy( dt2over2, M_Flux[j] , 1., M_rhs[i] );
                    M_divMatrixDiffSrc[2*i + j]->GlobalAssemble();
                    M_rhs[i] += *M_divMatrixDiffSrc[2*i + j]*(dt2over2*M_Flux[j]);

                    //! rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
                    //M_stiffMatrixDiffFlux[2*i+j].Axpy( -dt2over2, M_Flux[j] , 1., M_rhs[i] );
                    M_stiffMatrixDiffFlux[2*i + j]->GlobalAssemble();
                    M_rhs[i] += *M_stiffMatrixDiffFlux[2*i + j]*(- dt2over2*M_Flux[j]);

                    //! rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
                    //M_massMatrixDiffSrc[2*i+j].Axpy( dt2over2, M_Source[j] , 1., M_rhs[i] );
                    M_massMatrixDiffSrc[2*i + j]->GlobalAssemble();
                    M_rhs[i] += *M_massMatrixDiffSrc[2*i + j]*(dt2over2*M_Source[j]);
                }
        }

    Debug( 6310 ) << "[timeAdvance] \tcomputed rhs\n";
    //!---------------------------------------------------
    //! 3/ take into account the BOUNDARY CONDITIONS
    //!---------------------------------------------------
    //! compute the values for the boundary conditions
    Debug( 6310 ) << "[timeAdvance] \tcompute BC\n";

//     std::cout <<  M_U_thistime[0](101) << std::endl;
//     std::cout <<  M_U_thistime[1](101) << std::endl;
//     std::cout <<  M_U_thistime[2](101) << std::endl;
//     std::cout <<  M_U_thistime[3](101) << std::endl;

    _computeBC( time_val );
    //! take into account the bc
    Debug( 6310 ) << "[timeAdvance] \tcompute BC dirichlet vector\n";
    _updateBCDirichletVector();

    Debug( 6310 ) << "[timeAdvance] \trhs0 norm2 = " << M_rhs[0].Norm2() << "\n";
    Debug( 6310 ) << "[timeAdvance] \trhs1 norm2 = " << M_rhs[1].Norm2() << "\n";

    // *******************************************************
    chronoall.stop();
    Debug( 6310 ) << "[timeAdvance] ovall computation time " << chronoall.diff() << " s.\n";


}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::iterate( const Real& time_val , const int& count)
{
    Debug( 6310 ) << "[iterate] o-  Solving the system... t = " << time_val
                  << ", iter = " << count  << "... \n";
    Chrono chrono;

    Chrono chrono1;
    Chrono chrono2;
    Chrono chrono3;
    Chrono chrono4;
    Chrono chrono5;


    if( M_UW )
        {
            Real Ainode, Qinode;

            Real lambda1_plus  = 0.;
            Real lambda2_plus  = 0.;

            Real lambda1_minus = 0.;
            Real lambda2_minus = 0.;

            Real eigval1, eigval2;
            Real tmp11, tmp12, tmp21, tmp22;

            Real deltaX = M_FESpace.mesh()->edgeList( 1 ).point(1).x() - M_FESpace.mesh()->edgeList( 1 ).point(2).x();

            //! working on riemann invariants
            vector_type W1_UW(M_localMap);
            vector_type W2_UW(M_localMap);

            //! converting boundary conditions on physical variables
            //! in boundary conditions on characteristic variables
            M_oneDParams.W_from_U( W1_UW(0), W2_UW(0),
                                   M_rhs[0][0], M_rhs[1][0], 0 );
            M_oneDParams.W_from_U( W1_UW(M_FESpace.dim()-1), W2_UW(M_FESpace.dim()-1),
                                   M_rhs[0][M_FESpace.dim()-1], M_rhs[1][M_FESpace.dim()-1],
                                   M_FESpace.dim() - 1 );

            for ( UInt ii=1; ii < (M_FESpace.dim()-1) ; ii++ ) {
                //! compute the eigenvalues at node
                Ainode = M_U_thistime[0]( ii );
                Qinode = M_U_thistime[1]( ii );
                M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
                                                        eigval1, eigval2,
                                                        tmp11, tmp12,
                                                        tmp21, tmp22,
                                                        ii );

                lambda1_plus = std::max<Real>( eigval1, 0. );
                lambda1_minus = std::min<Real>( eigval1, 0. );
                lambda2_plus = std::max<Real>( eigval2, 0. );
                lambda2_minus = std::min<Real>( eigval2, 0. );
                //! update the solution for the next time step
                W1_UW[ii] = M_U_thistime[2][ii]
                    - (M_data.timestep() / deltaX) * lambda1_plus * ( M_U_thistime[2][ii] -
                                                                      M_U_thistime[2][ii-1])
                    - (M_data.timestep() / deltaX) * lambda1_minus * ( M_U_thistime[2][ii+1] -
                                                                       M_U_thistime[2][ii])
                    - M_data.timestep() * ( tmp11 * M_Source[0][ii] + tmp12 * M_Source[1][ii] );
                W2_UW[ii] = M_U_thistime[3][ii]
                    - (M_data.timestep() / deltaX) * lambda2_plus * ( M_U_thistime[3][ii] -
                                                                      M_U_thistime[3][ii-1])
                    - (M_data.timestep() / deltaX) * lambda2_minus * ( M_U_thistime[3][ii+1] -
                                                                       M_U_thistime[3][ii])
                    - M_data.timestep() * ( tmp21 * M_Source[0][ii] + tmp22 * M_Source[1][ii] );
            }

            M_U_thistime[2] = W1_UW;
            M_U_thistime[3] = W2_UW;

            for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ )
                M_oneDParams.U_from_W( M_U_thistime[0][ielem], M_U_thistime[1][ielem],
                                       M_U_thistime[2][ielem], M_U_thistime[3][ielem],
                                       ielem );

        }
    else
        {
            chrono.start();
            //! cholesky or lapack lu solve
            //! solve the system: rhs1 = massFactor^{-1} * rhs1

            chrono1.start();
            vector_type sol0(M_rhs[0]);

            //matrix_ptrtype matrFull( new matrix_type( M_localMap, M_factorMassMatrix.getMeanNumEntries()));
            matrix_ptrtype matrFull( new matrix_type( M_massMatrix ));
            //M_massMatrix.spy("mass");
            _updateBCDirichletMatrix( *matrFull );
            chrono1.stop();
            //M_factorMassMatrix.GlobalAssemble();
            //*matrFull = M_massMatrix;
            //matrFull->GlobalAssemble();

            //matrFull->spy("massmatr");

            chrono2.start();
            //            matrFull->spy("matr");
            M_linearSolver.setMatrix(*matrFull);
            M_linearSolver.setReusePreconditioner(false);
            M_linearSolver.solveSystem( M_rhs[0], sol0, matrFull );
            //std::cout << "sol0 norm2 = " << sol0.Norm2() << std::endl;
            chrono2.stop();

            chrono3.start();
            M_rhs[0] = sol0;
            vector_type sol1(M_rhs[1]);

            //std::cout << "rhs1 norm2 = " << M_rhs[1].Norm2() << std::endl;
//             //! solve the system: rhs2 = massFactor^{-1} * rhs2
//            M_linearSolver.setMatrix
            M_linearSolver.setReusePreconditioner(false);
            M_linearSolver.solveSystem( M_rhs[1], sol1, matrFull );
            //std::cout << "sol1 norm2 = " << sol1.Norm2() << std::endl;


            M_rhs[1] = sol1;

            chrono3.stop();
            //! correct flux with inertial term
            chrono4.start();
            if( M_inertial_wall )
                {
                    M_U_thistime[M_variable_index_map.find("Q_inert")->second] =
                        _correct_flux_inertial( M_rhs[1] );
                    M_rhs[1] += M_U_thistime[M_variable_index_map.find("Q_inert")->second];
                }
            //! correct flux with viscoelastic term
            if( M_viscoelastic_wall )
                {
                    M_U_thistime[M_variable_index_map.find("Q_visc")->second] =
                        _correct_flux_viscoelastic( M_rhs[1] );
                    M_rhs[1] += M_U_thistime[M_variable_index_map.find("Q_visc")->second];
                }
            //! compute L2 projection of d2Q_dx2
            //    if( M_flux_second_der )
            //      M_d2_U2_dx2 = _compute_d2Q_dx2( M_rhs[1] );

            //! correct flux with longitudinal term
            if( M_longitudinal_wall )
                {
//                     M_U_thistime[M_variable_index_map.find("Q_long")->second] =
//                         _correct_flux_longitudinal(  );
                    M_rhs[1] += M_U_thistime[M_variable_index_map.find("Q_long")->second];
                }

            //! store solution at previous timesteps & update the solution for the next time step
            for( UInt i=0; i<2; ++i )
                {
                    M_U_2prevtime[i] = M_U_prevtime[i];
                    M_U_prevtime [i] = M_U_thistime[i];
                    M_U_thistime [i] = M_rhs[i];
                }

            //      std::cout << std::endl;

            chrono4.stop();
            chrono5.start();
            Vector pressures(4 * M_FESpace.dim());
            for (UInt ielem = 0; ielem < M_FESpace.dim() ; ielem++ )
            {
                M_oneDParams.W_from_U( M_U_thistime[2][ielem + 1], M_U_thistime[3][ielem + 1],
                                       M_U_thistime[0][ielem + 1], M_U_thistime[1][ielem + 1],
                                       ielem);
                //        for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ ) {
                subrange(pressures, 4*ielem, 4 + 4*ielem) =
                    M_oneDParams.pressure( M_U_thistime [0][ielem + 1],
                                           M_U_prevtime [0][ielem + 1],
                                           M_U_2prevtime[0][ielem + 1],
                                           M_data.timestep(), ielem,
                                           M_dP_dt_steps, M_viscoelastic_wall,
                                           M_linearize_string_model );

                M_U_thistime[4][ielem + 1] = pressures(4*ielem);
            }

            if(M_viscoelastic_wall)
                {
                    for (UInt ielem=0; ielem <= M_FESpace.dim() ; ielem++ )
                        {
                            M_U_thistime[M_variable_index_map.find("P_elast")->second][ielem] =
                                pressures(1+4*ielem);
                            M_U_thistime[M_variable_index_map.find("P_visc")->second][ielem] =
                                pressures(2+4*ielem);
                            M_U_thistime[M_variable_index_map.find("dA_dt")->second][ielem] =
                                pressures(3+4*ielem);
                        }
                }

            chrono5.stop();
            chrono.stop();

        }



    Debug( 6310 ) << "chrono1 " << chrono1.diff() << "s.\n";
    Debug( 6310 ) << "chrono2 " << chrono2.diff() << "s.\n";
    Debug( 6310 ) << "chrono3 " << chrono3.diff() << "s.\n";
    Debug( 6310 ) << "chrono4 " << chrono4.diff() << "s.\n";
    Debug( 6310 ) << "chrono5 " << chrono5.diff() << "s.\n";
    Debug( 6310 ) << "chrono  " << chrono.diff() << " s.\n";

    Debug( 6310 ) << "[iterate] \tdone in " << chrono.diff() << " s.\n";


    //output2FileBuffers( time_val );

    //    postProcess( time_val );
}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::openFileBuffers()
{
    std::string file_output;

    boost::shared_ptr<std::ostringstream> buf;

    M_variable_index_iter iter_variable;
    M_variable_filter_iter iter_suffix;

    for( iter_variable = M_variable_index_map.begin();
         iter_variable != M_variable_index_map.end(); ++iter_variable ) {

        for( iter_suffix = M_variable_filter_map.begin();
             iter_suffix != M_variable_filter_map.end(); ++iter_suffix ) {

            file_output = M_data.PostDirectory() + "/" + M_data.PostFile() + iter_variable->first + iter_suffix->first;
            Debug( 6310 ) << "[openFileBuffers] setting output for file " << file_output << "\n";

            buf.reset( new std::ostringstream() );
            buf->setf(std::ios_base::scientific);
            buf->precision(5);
            buf->width(13);

            M_post_process_buffer.insert( std::map<std::string,
                                          boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );

            M_post_process_buffer_offset.insert( std::map<std::string,
                                                 long>::value_type( file_output, buf->tellp() ) );
        }
    }
}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::output2FileBuffers( const Real& time_val )
{

    Debug( 6310 ) << "[output2FileBuffers] begin \n";

    if( !( static_cast<int>( std::floor( (time_val/M_data.timestep()) + .5 ) ) % M_data.verbose() ) )
        {

            Debug( 6310 ) << "[output2FileBuffers] writting output \n";

            std::string file_output;

            M_variable_index_iter  iter_variable;
            M_variable_filter_iter iter_suffix;

            for( iter_variable = M_variable_index_map.begin(); iter_variable != M_variable_index_map.end(); ++iter_variable )
                {
                    for( iter_suffix = M_variable_filter_map.begin(); iter_suffix != M_variable_filter_map.end(); ++iter_suffix )
                        {

                            file_output = M_data.PostDirectory() + "/" + M_data.PostFile() +
                                iter_variable->first + iter_suffix->first;

                            Debug( 6310 ) << "[output2FileBuffers] setting output for file "
                                          << file_output << ", writing variable "
                                          << iter_variable->first << "\n";
                            Debug( 6310 ) << "[output2FileBuffers] variable size = "
                                          << M_U_thistime[iter_variable->second].size()
                                          << " variable norm2 = " << M_U_thistime[iter_variable->second].Norm2() << "\n";

                            (this->*(iter_suffix->second))
                                ( file_output, time_val, M_U_thistime[iter_variable->second], iter_variable->first );
                        }
                }
        }

}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::closeFileBuffers()
{
    // as I have a boost::shared_ptr, I expect the objects to be deallocated
    // now that the pointers are destroyed
    M_post_process_buffer.erase( M_post_process_buffer.begin(),
        M_post_process_buffer.end() );
    M_post_process_buffer_offset.erase( M_post_process_buffer_offset.begin(),
        M_post_process_buffer_offset.end() );
}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::postProcess( const Real& time_val )
{

    std::string str;

    std::ofstream outfile;

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator it;
    std::map< std::string, UInt>::iterator iter;

//     if (time_val==0.)
//         {
//             // the code is entering this for sure, as
//             // initialize invokes postProcess( 0. )

//             std::string file_output;

//             Real deltax( M_oneDParams.Length() /
//                          static_cast<Real>(M_oneDParams.ParamSize() - 1) );

//             file_output = "dati.m";

//             outfile.open(file_output.c_str(), std::ios::app );
//             outfile << "z" << M_data.PostFile()
//                     << " = (" << M_data.xLeft()
//                     << ":" << deltax
//                 << ":" << M_data.xRight() << ");\n"
//                     << std::endl;
//             outfile.close();

//         }

//     for (it = M_post_process_buffer.begin(); it != M_post_process_buffer.end(); ++it)
//         std::cout << it->first << std::endl;


    Debug( 6310 ) << "[postProcess] o- Dumping solutions on files (1d)!" << "/n";
    // dump solutions on files (buffers must be active!)

    //    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end(); ++iter, ++count)
    int count = 0;
    it  = M_post_process_buffer.begin();


    for( iter = M_variable_index_map.begin(); iter != M_variable_index_map.end(); ++iter, ++it)
    {


        std::string varname  = iter->first;

        int offset           = iter->second;

        std::string filename = it->first;
        Debug( 6310 ) << "Writing file " << filename << " with var " << offset << "/n";

        outfile.open( filename.c_str(), std::ios::app );
        //outfile << "# time = " << time_val << std::endl;

        for (int ii = 0; ii < M_dimDof; ++ii)
        {
            outfile << M_U_thistime[offset](ii + 1) << " ";
            //outfile << (*(*iter).second).str();
        }
        outfile << std::endl;
        outfile.close();

    }
    resetFileBuffers();

};



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::resetFileBuffers( )

{
    closeFileBuffers();
    openFileBuffers();
}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::seekpFileBuffers( )
{

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end();
         iter++ ){
        (*iter).second->seekp( M_post_process_buffer_offset[(*iter).first] );

    }

}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::tellpFileBuffers( )
{

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = M_post_process_buffer.begin(); iter != M_post_process_buffer.end();
         iter++ ){
        M_post_process_buffer_offset[(*iter).first] = (*iter).second->tellp();

    }

}



template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::gplot( )
{
    //M_GracePlot.Plot( M_FESpace.mesh()->pointList(), M_U_thistime[0] );
}

//! output for Plotmtv.
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::output_to_plotmtv(std::string fname, Real time_val,
                                                        const ScalVec& U,
                                                        std::string vname )
{

    boost::shared_ptr<std::ostringstream> buf;

    buf = M_post_process_buffer[fname];

    (*buf) << "$ DATA = CURVE2D\n % xlabel='z'\n"
           << "% toplabel='Section,time=" << time_val << "'\n % ylabel='" << vname << "'\n"
           << std::flush;

    for(UInt ii = 0; ii < U.size(); ii++){
        (*buf) << M_FESpace.mesh()->pointList()[ii].x() << " " << U[ii] << "\n" << std::flush;
    }


}

//! output for Matlab.
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::output_to_matlab( std::string  fname,
                                                         Real         /*time_val*/,
                                                         const        vector_type& U,
                                                         std::string  /*vname*/ )
{

#if 0

    fstream filestr;

    filestr.open (fname.c_str(), fstream::in | fstream::out | fstream::app);

    //filestr << "# time " << time_val << std::endl;
    for ( int ii = LeftNodeId(); ii <= RightNodeId() ; ++ii )
        {
            //        (*buf) << U[ii] << "; ";
            //std::cout << U[ii] << " " << std::endl;
            filestr << U[ii] << " ";
            //test << U[ii] << std::endl;
        }
    filestr << "\n";
    filestr.close();

#else


    boost::shared_ptr<std::ostringstream> buf;

    //std::cout << "buffer = " << fname << std::endl;
    buf = M_post_process_buffer[fname];
    //buf = &std::cout;
//     std::ostringstream test;
//     test.open("test.txt");

    //    (*buf) << vname << M_data.PostFile()
    //           << "( " << (static_cast<int>( std::floor( time_val/M_data.timestep() + 0.5 ) )
    //              /  M_data.verbose() )+1
    //           << ", : ) = [ " << std::flush;

    //std::cout << "Norm2 U = " << U.Norm2() << std::endl;

    for ( int ii = LeftNodeId(); ii <= RightNodeId() ; ++ii )
        {
            //        (*buf) << U[ii] << "; ";
            //std::cout << U[ii] << std::endl;
            (*buf) << U[ii] << " ";
            //std::cout << buf->str();
            //test << U[ii] << std::endl;
        }
    //    (*buf) << "]';\n" << std::endl;
    (*buf) << "\n";
    //std::cout << (*buf) << std::endl;
#endif

}

//! Create a Matlab script to visualize output matlab files
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::create_movie_file()
{
    std::ofstream outfile;
    std::string file_output;
    file_output = M_data.PostDirectory() + "/" + "areamovie"+M_data.PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Area"<<M_data.PostFile()<<";\n"
            <<"[n,m]=size(A" << M_data.PostFile() << ");\n"
            <<"Max=max(max(A" << M_data.PostFile() << "));\n"
            <<"Min=min(min(A" << M_data.PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(A" << M_data.PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<<M_data.timestep()<<"*"<<M_data.verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";

    outfile.close();
    file_output = M_data.PostDirectory() + "/" + "portatamovie"+M_data.PostFile()+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Portata"<<M_data.PostFile()<<";\n"
            <<"[n,m]=size(Q" << M_data.PostFile() << ");\n"
            <<"Max=max(max(Q" << M_data.PostFile() << "));\n"
            <<"Min=min(min(Q" << M_data.PostFile() << "));\n"
            <<"for i=1:n\n"
            <<"  plot(Q" << M_data.PostFile() << "(i,:));\n"
            <<"  title((i-1)*"<< M_data.timestep() <<"*"<<M_data.verbose()<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";
    outfile.close();

}


//! Print to screen information on the Solver class
template< class Params, class Flux, class Source >
void
OneDModelSolver<Params, Flux, Source>::showMe(std::ostream& c, UInt verbose)
{
    c << "\n--- One Dimensional Model Data\n";
    //this->showMe(c);
//     c << "\n--- One Dimensional Model Handler\n";
//     this->showMeHandler(c, verbose);
//     c << "\n--- One Dimensional Model Parameters\n";
    this->oneDParams().showMe(c);
    c << "--- End of One Dimensional Model\n";
}


}

#endif
