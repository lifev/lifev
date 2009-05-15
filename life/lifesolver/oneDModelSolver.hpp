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
  two functions (_M_fluxFun and _M_sourceFun) have to be defined
  separately to allow the update (_updateFluxDer, _updateSourceDer)
  of the corresponding vectors (_M_Fluxi, _M_diffFluxij,
  _M_Sourcei, _M_diffSrcij). I also separated the treatment of the
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

#include <life/lifesolver/oneDModelHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>

#include <life/lifefem/assembGeneric.hpp>
#include <life/lifecore/chrono.hpp>

#include <life/lifearray/tridiagMatrix.hpp>
#include <life/lifealg/triDiagCholesky.hpp>
#include <life/lifealg/triDiagLU.hpp>
#include <life/lifesolver/oneDNonLinModelParam.hpp>
#include <life/lifesolver/vectorFunction1D.hpp>

#include <life/lifefem/oneDBCHandler.hpp>

#include <boost/numeric/ublas/vector_proxy.hpp>

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
template<typename PARAM, typename FLUX, typename SOURCE>
class OneDModelSolver:
        public OneDModelHandler
{
public:

    /*! \name Typedefs
     */
    //@{
    typedef PARAM param_type;

    typedef FLUX flux_type;

    typedef SOURCE source_type;
    //@}

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param onedparam Variable containing parameters for OneD model
    */
    OneDModelSolver(const GetPot& data_file); //,

    ~OneDModelSolver(){}

    //! return the solution at current time step (Area)
    const ScalVec& U1_thistime() const { return _M_U_thistime[0];}

    //! return the solution at current time step (Flux)
    const ScalVec& U2_thistime() const { return _M_U_thistime[1];}

    //! return the reimann invariant W1 at current time step
    const ScalVec& W1_thistime() const { return _M_U_thistime[2];}

    //! return the reimann invarant W2 at current time step
    const ScalVec& W2_thistime() const { return _M_U_thistime[3];}

    //! return the solution at current time step (U)
    const ScalVec_vector& U_thistime() const { return _M_U_thistime;}

    //! return a const reference to parameter class
    const PARAM& oneDParam() const { return _M_oneDParam; }

    //! return the BC handler
    OneDBCHandler<FLUX>& bcH() { return _M_bcH;}

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
    void initialize(GetPot& data_file);

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

    //! set the Dirichlet boundary conditions (left)
    void setBCValuesLeft( const Real& bcL1, const Real& bcL2 );

    //! set the Dirichlet boundary conditions (right)
    void setBCValuesRight( const Real& bcR1, const Real& bcR2 );

    //! set left bctype to internal node
    void setBCLeft_internalnode();

    //! set right bctype to internal node
    void setBCRight_internalnode();

    //! get the flux function
    FLUX const& FluxFun() const;
    //! get the source function
    SOURCE const& SourceFun() const;

    //! get the left edge
    Edge1D LeftEdge() const;
    //! get the right edge
    Edge1D RightEdge() const;

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
    void output_to_matlab( std::string fname, Real time_val,
                           const ScalVec& U,
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
    void showMe( std::ostream& c, UInt verbose = false );

private:

    //! the parameters
    PARAM _M_oneDParam;
    // const LinearSimpleParam& _M_oneDParam;

    //! the flux function
    FLUX _M_fluxFun;
    //! the source function
    SOURCE _M_sourceFun ;

    const UInt _M_leftNodeId;
    const UInt _M_leftInternalNodeId;
    const UInt _M_rightNodeId;
    const UInt _M_rightInternalNodeId;

    //! boundary edges
    const Edge1D _M_leftEdge;
    const Edge1D _M_rightEdge;

    //! coefficient in front of the corresponding _M_elmat*
    Real _M_coeffMass;
    Real _M_coeffStiff;
    Real _M_coeffGrad;
    Real _M_coeffDiv;

    ElemMat _M_elmatMass;  //!< element mass matrix
    ElemMat _M_elmatStiff; //!< element stiffness matrix
    ElemMat _M_elmatGrad;  //!< element gradient matrix
    ElemMat _M_elmatDiv;   //!< element divergence matrix

    //  ElemVec _M_elvec; // Elementary right hand side

    //! Unknowns at present time step
    /*!
      U is a vector of ScalVec:
      U[0] = A, U[1] = Q
      U[2] = W1, U[3] = W2
      Other components of U may contain additional variables such as the pressure
    */
    ScalVec_vector _M_U_thistime;

    //! Unknowns at previous time step (see savesol() )
    ScalVec_vector _M_U_prevtime;
    ScalVec_vector _M_U_2prevtime;

    //! Right hand sides of the linear system i: "mass * _M_Ui = _M_rhsi"
    ScalVec_vector _M_rhs;

    //! Flux F(U) (in P1)
    std::vector<ScalVec > _M_Flux;

    //! diffFlux = dF(U)/dU (in P0)
    std::vector<ScalVec > _M_diffFlux;

    //! Source term S (in P1)
    std::vector<ScalVec > _M_Source;
    //! diffSrc = dSource(U)/dU (in P0)
    std::vector<ScalVec > _M_diffSrc;



    //! tridiagonal mass matrix
    TriDiagMatrix<Real> _M_massMatrix;

    //! factorized tridiagonal mass matrix
    TriDiagMatrix<Real> _M_factorMassMatrix;

    //!@{  TO BE Templatized !!
    //! cholesky factorization
    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _M_tridiagSlv;
    //! lapack LU factorization
    // TriDiagLU< Real, TriDiagMatrix<Real>, Vector > _M_tridiagSlv;
    //!@} end of  TO BE Templatized !!


    //! tridiagonal mass matrices multiplied by diffSrcij
    std::vector<TriDiagMatrix<Real> > _M_massMatrixDiffSrc;

    //! tridiagonal stiffness matrices multiplied by diffFluxij
    std::vector<TriDiagMatrix<Real> > _M_stiffMatrixDiffFlux;

    //! tridiagonal gradient matrix
    TriDiagMatrix<Real> _M_gradMatrix;
    //! tridiagonal gradient matrices multiplied by diffFluxij
    std::vector<TriDiagMatrix<Real> > _M_gradMatrixDiffFlux;

    //! tridiagonal divergence matrices multiplied by diffSrcij
    std::vector<TriDiagMatrix<Real> > _M_divMatrixDiffSrc;


    //! Update the coefficients
    //! (from the flux, source functions and their derivatives)
    void _updateMatrixCoefficients(const UInt& ii, const UInt& jj,
                                   const UInt& iedge);

    //! Update the element matrices with the current element
    void _updateElemMatrices();

    //! assemble the matrices
    int _assemble_matrices(const UInt& ii, const UInt& jj );

    /*! update the matrices
      _M_massMatrixDiffSrcij, _M_stiffMatrixDiffFluxij
      _M_gradMatrixDiffFluxij, and _M_divMatrixDiffSrcij (i,j=1,2)

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
    void _updateBCDirichletMatrix( TriDiagMatrix<Real>& mat );

    /*! modify the vector to take into account
      the Dirichlet boundary conditions
      (works for P1Seg and canonic numbering!)
    */
    void _updateBCDirichletVector();

    /*! compute the _M_bcDirLeft and _M_bcDirRight
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

    //! update the P1 flux vector from U: _M_Fluxi = F_h(Un) i=1,2
    void _updateFlux();
    //! update the P1 source vector from U: _M_Sourcei = S_h(Un) i=1,2
    void _updateSource();

    //! call _updateFlux and update the P0 derivative of flux vector from U:
    //! _M_diffFluxij = dF_h/dU(Un) i,j=1,2
    void _updateFluxDer();
    //! call _updateSource and update the P0 derivative of source vector from U:
    //! _M_diffSrcij = dS_h/dU(Un) i,j=1,2
    void _updateSourceDer();

    //! to manage non elastic wall behaviour
    ScalVec _correct_flux_inertial(const ScalVec &);
    ScalVec _correct_flux_viscoelastic(const ScalVec &);
    ScalVec _correct_flux_longitudinal( );

    ScalVec _compute_d2Q_dx2( const ScalVec& );

    //! ostringstream buffers to store the solution before writing it to file
    std::map<std::string, boost::shared_ptr<std::ostringstream> >
    _M_post_process_buffer;

    //! position of the put pointer to ostringstream buffers
    std::map<std::string, long> _M_post_process_buffer_offset;

    //! Class to manage boundary conditions
    OneDBCHandler<FLUX> _M_bcH;

    //! Trick to use strings in C++ "switch" construct
    std::map<std::string, OneDInitializeVar> _M_oneDstring2initializeVarMap;

    //! boolean: show CFL value during computations?
    bool _M_CFL;

    //! boolean: use alternate solver (UpWind)?
    bool _M_UW;

    //! boolean: activate inertial/ viscoelastic/ longitudinal term in pressure-area relationship?
    bool _M_inertial_wall;
    bool _M_viscoelastic_wall;
    bool _M_linearize_string_model;
    bool _M_linearize_equations;
    bool _M_longitudinal_wall;

    //! boolean: compute second spatial derivative of flux?
    bool _M_flux_second_der;

    //! approximation of pressure temporal derivative: how many time steps?
    int _M_dP_dt_steps;

    //! member function pointer: to "automagically" manage postprocess
    typedef void (LifeV::OneDModelSolver<PARAM, FLUX, SOURCE>::*postproc_funptr)( std::string,
        Real, const ScalVec& U, std::string );

    // maps associating each unknown to a string, a numeric index, a function pointer
    std::map< std::string, std::string> _M_variable_string_map;
    std::map< std::string, UInt> _M_variable_index_map;
    std::map< std::string, postproc_funptr> _M_variable_filter_map;
    typedef std::map< std::string, UInt>::iterator _M_variable_index_iter;
    typedef typename std::map< std::string, postproc_funptr >::iterator _M_variable_filter_iter;
    typedef std::map< std::string, std::string>::iterator _M_variable_string_iter;
};


// ***********************************
// IMPLEMENTATION
// ***********************************


template< class PARAM, class FLUX, class SOURCE >
OneDModelSolver<PARAM, FLUX, SOURCE>::OneDModelSolver(const GetPot& data_file):
    OneDModelHandler(data_file),
    _M_oneDParam(data_file),
    _M_fluxFun(_M_oneDParam),
    _M_sourceFun(_M_oneDParam),
    //! id of left and right bc nodes
    _M_leftNodeId( 0 ),
    _M_leftInternalNodeId( _M_leftNodeId + 1 ),
    _M_rightNodeId( _M_dimDof - 1 ),
    _M_rightInternalNodeId( _M_rightNodeId  - 1 ),
    //! boundary edges
    _M_leftEdge( _M_mesh.edgeList( 1 ) ),
    _M_rightEdge( _M_mesh.edgeList( _M_nb_elem ) ),
    // elementary matrices
    _M_elmatMass (_M_fe.nbNode,1,1),
    _M_elmatStiff(_M_fe.nbNode,1,1),
    _M_elmatGrad (_M_fe.nbNode,1,1),
    _M_elmatDiv  (_M_fe.nbNode,1,1),
    // vectorial unknowns and rhs
    _M_U_thistime(), // size should be at least 4
    _M_U_prevtime(2), // probably useless - could be components in U_thistime
    _M_U_2prevtime(2),
    _M_rhs(2),
    // vectors and matrices of the non-linear function
    _M_Flux(2),
    _M_diffFlux(4),
    _M_Source(2),
    _M_diffSrc(4),
    // mass matrix (to be inverted)
    _M_massMatrix(_M_dimDof),
    _M_factorMassMatrix(_M_dimDof),
    _M_tridiagSlv(_M_dimDof),
    // matrices used to build the rhs
    _M_gradMatrix(_M_dimDof),
    // Handle boundary conditions
    _M_bcH(_M_U_thistime, /*_M_postproc_variable,*/ _M_fluxFun, _M_dimDof),
    _M_CFL(data_file("miscellaneous/showCFL",0)),
    _M_UW(data_file("miscellaneous/alternate_solver",0)),
    _M_inertial_wall(data_file("miscellaneous/inertial_wall",0)),
    _M_viscoelastic_wall(data_file("miscellaneous/viscoelastic_wall",0)),
    _M_linearize_string_model(data_file("miscellaneous/linearize_string_model",1)),
    _M_linearize_equations(data_file("miscellaneous/linearize_equations",0)),
    _M_longitudinal_wall(data_file("miscellaneous/longitudinal_wall",0)),
    _M_flux_second_der(data_file("miscellaneous/compute_flux_second_derivative",0)),
    _M_dP_dt_steps(data_file("miscellaneous/pressure_derivative_steps",2))
{
    // These maps allow a more readable definition of the variables
    // the model is describing
    _M_variable_string_map["A"] = "Area";
    _M_variable_string_map["Q"] = "Flow Rate";
    _M_variable_string_map["P"] = "Pressure";
    _M_variable_string_map["W1"] = "First Reimann Invariant";
    _M_variable_string_map["W2"] = "Second Reimann Invariant";
    _M_variable_string_map["Q_visc"] = "Flow Rate correction (viscoelastic)";
    _M_variable_string_map["Q_inertial"] = "Flow Rate correction (inertial)";
    _M_variable_string_map["Q_longitudinal"] = "Flow Rate correction (longitudinal)";
    _M_variable_string_map["d2Q_dx2"] = "Flow Rate second derivative";

    // These maps allow to use the switch... case construct
    // with string variables
    _M_oneDstring2initializeVarMap["A"] = OneDInitArea;
    _M_oneDstring2initializeVarMap["Q"] = OneDInitFlux;
    _M_oneDstring2initializeVarMap["W1"] = OneDInitReimann1;
    _M_oneDstring2initializeVarMap["W2"] = OneDInitReimann2;
    _M_oneDstring2initializeVarMap["P"] = OneDInitPressure;

    Debug( 6310 ) << "[OneDModelSolver] O-  Nb of unknowns: " << _M_dimDof     << "\n";
    Debug( 6310 ) << "[OneDModelSolver] O-  Computing the constant matrices... \n";
    Debug( 6310 ) << "[OneDModelSolver] O-  Adopting a"
                  << std::string( _M_viscoelastic_wall ? " viscoelastic " : "n elastic " )
                  << "model for vessel wall... \n";

    Chrono chrono;
    chrono.start();

    //! Vector and Matrices initialization
    _M_massMatrix.zero();
    // _M_factorMassMatrix.zero();
    _M_factorMassMatrix.zero();
    _M_gradMatrix.zero();

    // insert variables!
    UInt nvar(0);
    _M_variable_index_map.insert( make_pair("A", nvar++ ) );
    _M_variable_index_map.insert( make_pair("Q", nvar++ ) );
    _M_variable_index_map.insert( make_pair("W1", nvar++ ) );
    _M_variable_index_map.insert( make_pair("W2", nvar++ ) );

    _M_variable_index_map.insert( make_pair("P", nvar++ ) );

    //! correction flux with viscoelastic term
    if( _M_viscoelastic_wall ) {
        // correction on the flux due to the viscoelastic term
        _M_variable_index_map.insert( make_pair("Q_visc", nvar++ ) );
        // viscoelastic contribution to the pressure
        _M_variable_index_map.insert( make_pair("P_visc", nvar++ ) );
        // elastic contribution to the pressure
        _M_variable_index_map.insert( make_pair("P_elast", nvar++ ) );
        // time derivative of the section area
        _M_variable_index_map.insert( make_pair("dA_dt", nvar++ ) );
    }

    //! flux second derivative
    if( _M_flux_second_der )
        _M_variable_index_map.insert( make_pair("d2Q_dx2", nvar++ ) );

    //! correction flux with inertial term
    if( _M_inertial_wall )
        _M_variable_index_map.insert( make_pair("Q_inert", nvar++ ) );

    //! correction flux with longitudinal term
    if( _M_longitudinal_wall )
        _M_variable_index_map.insert( make_pair("Q_long", nvar++ ) );

    //! activate the export filters
    // matlab postprocessing
    _M_variable_filter_map.insert( make_pair(".m",
        &LifeV::OneDModelSolver<PARAM, FLUX, SOURCE>::output_to_matlab) );
    // plotmtv postprocessing
    //    _M_variable_filter_map.insert( make_pair(".mtv",
    //  &LifeV::OneDModelSolver<PARAM, FLUX, SOURCE>::output_to_plotmtv) );

    _M_U_thistime.resize(nvar);

    // initialize vectors
    std::fill( _M_U_thistime.begin(), _M_U_thistime.end(), ZeroVector(_M_dimDof) );
    std::fill( _M_U_prevtime.begin(), _M_U_prevtime.end(), ZeroVector(_M_dimDof) );
    std::fill( _M_U_2prevtime.begin(), _M_U_2prevtime.end(), ZeroVector(_M_dimDof) );
    std::fill( _M_rhs.begin(), _M_rhs.end(), ZeroVector(_M_dimDof) );
    std::fill( _M_Flux.begin(), _M_Flux.end(), ZeroVector(_M_dimDof) );
    std::fill( _M_Source.begin(), _M_Source.end(), ZeroVector(_M_dimDof) );

    // initialize matrices
    std::fill( _M_diffFlux.begin(), _M_diffFlux.end(), ZeroVector(_M_nb_elem) );
    std::fill( _M_diffSrc.begin(), _M_diffSrc.end(), ZeroVector(_M_nb_elem) );
    for( UInt i=0; i<4; ++i ) {
        // using push_back because of no default constructor available
        _M_massMatrixDiffSrc.push_back( TriDiagMatrix<Real>(_M_dimDof) );
        _M_stiffMatrixDiffFlux.push_back( TriDiagMatrix<Real>(_M_dimDof) );
        _M_gradMatrixDiffFlux.push_back( TriDiagMatrix<Real>(_M_dimDof) );
        _M_divMatrixDiffSrc.push_back( TriDiagMatrix<Real>(_M_dimDof) );
    }

    //-------------------------------------------
    //! update first the constant matrices (cst w.r. to time iter)

    //! set the coeff to 1.
    _M_coeffMass = 1.;
    _M_coeffGrad = 1.;

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _M_elmatMass.zero();
        _M_elmatGrad.zero();

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
        //std::cout << _M_fe.currentId() << std::endl;

        //! update the mass and grad matrices
        mass( _M_coeffMass, _M_elmatMass, _M_fe,0, 0 );
        grad( 0 , - _M_coeffGrad, _M_elmatGrad, _M_fe, _M_fe, 0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _M_massMatrix, _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );
        assemb_mat( _M_gradMatrix, _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

    } //! end loop on elements

    _M_factorMassMatrix = _M_massMatrix;

    //! Dirichlet boundary conditions set in the mass matrix
    _updateBCDirichletMatrix( _M_factorMassMatrix );

    //! cholesky or lapack lu factorization of the mass matrix
    _M_tridiagSlv.Factor( _M_factorMassMatrix );

    chrono.stop();

    Debug( 6310 ) << "[OneDModelSolver] \tdone in " << chrono.diff() << " s.\n";

    // for debugging purposes
    //    std::ostringstream output, dummy_sstr;
    //    std::string dummy_file("second_derivative");
    //    dummy_sstr << _M_mesh.numEdges();
    //    dummy_file += dummy_sstr.str() + ".m";
    //    std::ofstream outfile;
    //
    //    outfile.open(dummy_file.c_str(), std::ios::app);
    //
    //    outfile << "Qhat = [ ];\n";
    //    outfile << "rhs = [ ];\n";
    //
    //    outfile.close();

}

/*! Update the coefficients
  (from the flux, source functions and their derivatives)
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateMatrixCoefficients(const UInt& ii,
    const UInt& jj , const UInt& iedge)
{
    Real dFluxdUelem(0), dSrcdUelem(0);

    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    dFluxdUelem = _M_diffFlux[ 2*(ii-1) + jj-1 ]( iedge - 1 ); //! iedge starts from 1...
    dSrcdUelem  = _M_diffSrc[ 2*(ii-1) + jj-1 ]( iedge - 1 );

    _M_coeffGrad  = dFluxdUelem; //! term gradDiffFlux(U) [* S(U)]
    _M_coeffDiv   = dSrcdUelem;  //! term  divDiffSrc(U) [* F(U)]
    _M_coeffStiff = dFluxdUelem; //! term stiffDiffFlux(U) [* F(U)]
    _M_coeffMass  = dSrcdUelem;  //! term  massDiffSrc(U) [* S(U)]
}

//! Update the element matrices with the updated
//! current element and updated coefficients
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateElemMatrices()
{
    //! set the elementary matrices to 0.
    _M_elmatMass.zero();
    _M_elmatStiff.zero();
    _M_elmatGrad.zero();
    _M_elmatDiv.zero();

    //! update the mass matrix
    mass( _M_coeffMass, _M_elmatMass, _M_fe,0, 0 );
    //  std::cout << "Elem Mass matrix :" << std::endl;
    //  _M_elmatMass.showMe( std::cout );

    //! update the stiffness matrix
    stiff( _M_coeffStiff, _M_elmatStiff, _M_fe,0 ,0 );
    // std::cout << "Elem Stiff matrix :" << std::endl;
    // _M_elmatStiff.showMe( std::cout );

    /*! update the gradient matrix
      gradient operator:
      grad_{ij} = \int_{fe} coeff \phi_j \frac{d \phi_i}{d x}

      BEWARE :
      \param 0: the first argument "0" corresponds to the first
      and only coordinate (1D!), and HERE it starts from 0... (Damm'!)

      \param - _M_coeffGrad: the sign "-" in the second argument
      is added to correspond to the described operator.
      (There is a minus in the elemOper implementation).
    */
    grad( 0 , - _M_coeffGrad, _M_elmatGrad, _M_fe, _M_fe, 0, 0 );
    //  std::cout << "Elem Grad matrix :" << std::endl;
    //  _M_elmatGrad.showMe( std::cout );

    /*! update the divergence matrix
      divergence operator: (transpose of the gradient)
      div_{ij} = \int_{fe} coeff \frac{d \phi_j}{d x} \phi_i

      \note formally this _M_elmatDiv is not necessary
      as it is the transpose of the _M_elmatGrad.
      But for the sake of clarity, I prefer to keep it. (low cost!)

      BEWARE : same remarks as grad (see above).
    */
    div( 0 , - _M_coeffDiv, _M_elmatDiv, _M_fe, _M_fe, 0, 0 );
    //  std::cout << "Elem Div matrix :" << std::endl;
    //  _M_elmatDiv.showMe( std::cout );
}


//! assemble the matrices
template< class PARAM, class FLUX, class SOURCE >
int
OneDModelSolver<PARAM, FLUX, SOURCE>::
_assemble_matrices(const UInt& ii, const UInt& jj )
{
    ASSERT_BD( 0 < ii && ii < 3 && 0 < jj && jj < 3 );

    //! assemble the mass matrix
    assemb_mat( _M_massMatrixDiffSrc[ 2*(ii-1) + jj-1 ], _M_elmatMass, _M_fe, _M_dof1D , 0, 0 );

    //! assemble the stiffness matrix
    assemb_mat( _M_stiffMatrixDiffFlux[ 2*(ii-1) + jj-1 ], _M_elmatStiff, _M_fe, _M_dof1D , 0, 0 );

    //! assemble the gradient matrix
    assemb_mat( _M_gradMatrixDiffFlux[ 2*(ii-1) + jj-1 ], _M_elmatGrad, _M_fe, _M_dof1D , 0, 0 );

    //! assemble the divergence matrix
    assemb_mat( _M_divMatrixDiffSrc[ 2*(ii-1) + jj-1 ], _M_elmatDiv, _M_fe, _M_dof1D , 0, 0 );

    return 0;

    ERROR_MSG("Invalid values for the _assemble_matrices method.");
    return 1;
}

/*! update the matrices
  _M_massMatrixDiffSrcij, _M_stiffMatrixDiffFluxij
  _M_gradMatrixDiffFluxij, and _M_divMatrixDiffSrcij (i,j=1,2)

  from the values of diffFlux(Un) and diffSrc(Un)
  that are computed with _updateMatrixCoefficients.
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateMatrices()
{
    //--------------------------------------------------------
    // Chrono chrono;
    // chrono.start();
    // std::cout << "o-loop over the matrices INIT... ";

    //! Matrices initialization
    for( UInt i=0; i<4; ++i ) {
        _M_massMatrixDiffSrc[i].zero();
        _M_stiffMatrixDiffFlux[i].zero();
        _M_gradMatrixDiffFlux[i].zero();
        _M_divMatrixDiffSrc[i].zero();
    }

    /*
      chrono.stop();
      std::cout << "done in " << chrono.diff() << " s." << std::endl;
      chrono.start();
      std::cout << "o-loop over the matrices... ";
    */

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
        //  std::cout << _M_fe.currentId() << std::endl;

        for(UInt ii = 1; ii <= 2; ii ++) {
            for(UInt jj = 1; jj <= 2; jj ++) {

                //! update the _M_coeff*
                _updateMatrixCoefficients( ii , jj, iedge);

                //! update the _M_elmat*
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

    /*
    //! useless ??????
    //! taking into account the dirichlet bc
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc11 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc12 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc21 );
    _updateBCDirichletMatrix( _M_massMatrixDiffSrc22 );

    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux11 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux12 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux21 );
    _updateBCDirichletMatrix( _M_stiffMatrixDiffFlux22 );

    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux11 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux12 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux21 );
    _updateBCDirichletMatrix( _M_gradMatrixDiffFlux22 );

    _updateBCDirichletMatrix( _M_divMatrixDiffSrc11 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc12 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc21 );
    _updateBCDirichletMatrix( _M_divMatrixDiffSrc22 );
    // chrono.stop();
    // std::cout << "done in " << chrono.diff() << " s." << std::endl;
    */
}


/*! modify the matrix to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::
_updateBCDirichletMatrix( TriDiagMatrix<Real>& mat )
{
    UInt firstDof = 0;
    UInt lastDof  = mat.OrderMatrix()-1;

    /*
    //! unsymmetric treatment (LU must be used!)
    //! modify the first row
    mat.Diag()( firstDof )   = 1.;
    mat.UpDiag()( firstDof ) = 0.;

    //! modify the last row
    mat.Diag()( lastDof )      = 1.;
    mat.LowDiag()( lastDof-1 ) = 0.;
    */

    //! symmetric treatment (cholesky can be used)
    //! modify the first row
    mat.Diag()( firstDof )    = 1.;
    mat.UpDiag()( firstDof )  = 0.;
    mat.LowDiag()( firstDof ) = 0.; //!and second row

    //! modify the last row
    mat.Diag()( lastDof )      = 1.;
    mat.UpDiag()( lastDof-1 )  = 0.; //!and penultimate row
    mat.LowDiag()( lastDof-1 ) = 0.;
}

/*! modify the vector to take into account
  the Dirichlet boundary conditions
  (works for P1Seg and canonic numbering!)

  \param val_left  : Dirichlet value inserted to the left
  \param val_right : Dirichlet value inserted to the right
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::
_updateBCDirichletVector()
{
    UInt firstDof = _M_leftNodeId;
    UInt lastDof  = _M_rightNodeId;

    Debug( 6310 ) << "\n[_updateBCDirichletVector] firstDof = " << firstDof
                  << "\nlastDof = " << lastDof << ";\n";
    /*
    //! unsymmetric treatment (LU must be used!)
    //! first row modified
    _M_rhs[0]( firstDof ) = _M_bcDirLeft.first;
    _M_rhs[1]( firstDof ) = _M_bcDirLeft.second;

    //! last row modified
    _M_rhs[0]( lastDof ) = _M_bcDirRight.first;
    _M_rhs[1]( lastDof ) = _M_bcDirRight.second;
    */

    //! symmetric treatment (cholesky can be used)
    for(UInt i=0; i<2; ++i) {
        //! first row modified (Dirichlet)
        _M_rhs[i]( firstDof ) = _M_bcDirLeft[i];
        //! second row modified (for symmetry)
        _M_rhs[i]( firstDof + 1 ) += - _M_massMatrix.LowDiag()( firstDof ) * _M_bcDirLeft[i];
        //! last row modified (Dirichlet)
        _M_rhs[i]( lastDof ) = _M_bcDirRight[i];
        //! penultimate row modified (for symmetry)
        _M_rhs[i]( lastDof - 1 ) += - _M_massMatrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight[i];
    }

}


//! compute the _M_bcDirLeft and _M_bcDirRight and set them to the new values
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_computeBC( const Real& time_val )
{
    _M_bcH.applyBC(time_val, _M_bcDirLeft, _M_bcDirRight );
}


//! get the flux function
template< class PARAM, class FLUX, class SOURCE >
const FLUX&
OneDModelSolver<PARAM, FLUX, SOURCE>::FluxFun() const
{
    return _M_fluxFun;
}


//! get the source function
template< class PARAM, class FLUX, class SOURCE >
const SOURCE&
OneDModelSolver<PARAM, FLUX, SOURCE>::SourceFun() const
{
    return _M_sourceFun;
}


//! get the left edge
template< class PARAM, class FLUX, class SOURCE >
Edge1D
OneDModelSolver<PARAM, FLUX, SOURCE>::LeftEdge() const
{
    return _M_leftEdge;
}


//! get the right edge
template< class PARAM, class FLUX, class SOURCE >
Edge1D
OneDModelSolver<PARAM, FLUX, SOURCE>::RightEdge() const
{
    return _M_rightEdge;
}


//! get the left node
template< class PARAM, class FLUX, class SOURCE >
UInt
OneDModelSolver<PARAM, FLUX, SOURCE>::LeftNodeId() const
{
    return _M_leftNodeId;
}


//! get the left internal node (neighboring node)
template< class PARAM, class FLUX, class SOURCE >
UInt
OneDModelSolver<PARAM, FLUX, SOURCE>::LeftInternalNodeId() const
{
    return _M_leftInternalNodeId;
}


//! get the right node
template< class PARAM, class FLUX, class SOURCE >
UInt
OneDModelSolver<PARAM, FLUX, SOURCE>::RightNodeId() const
{
    return _M_rightNodeId;
}


//! get the right internal node (neighboring node)
template< class PARAM, class FLUX, class SOURCE >
UInt
OneDModelSolver<PARAM, FLUX, SOURCE>::RightInternalNodeId() const
{
    return _M_rightInternalNodeId;
}


//! get the Dirichlet boundary conditions (left)
template< class PARAM, class FLUX, class SOURCE >
typename OneDModelSolver<PARAM, FLUX, SOURCE>::Vec2D
OneDModelSolver<PARAM, FLUX, SOURCE>::BCValuesLeft() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = _M_U_thistime[0]( LeftNodeId() );
    temp[1] = _M_U_thistime[1]( LeftNodeId() );
    return temp;
}


//! get the value at neighboring node (left)
template< class PARAM, class FLUX, class SOURCE >
typename OneDModelSolver<PARAM, FLUX, SOURCE>::Vec2D
OneDModelSolver<PARAM, FLUX, SOURCE>::BCValuesInternalLeft() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = _M_U_thistime[0]( LeftInternalNodeId() );
    temp[1] = _M_U_thistime[1]( LeftInternalNodeId() );
    return temp;
}


//! get the Dirichlet boundary conditions (right)
template< class PARAM, class FLUX, class SOURCE >
typename OneDModelSolver<PARAM, FLUX, SOURCE>::Vec2D
OneDModelSolver<PARAM, FLUX, SOURCE>::BCValuesRight() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = _M_U_thistime[0]( RightNodeId() );
    temp[1] = _M_U_thistime[1]( RightNodeId() );
    return temp;
}


//! get the value at neighboring node (right)
template< class PARAM, class FLUX, class SOURCE >
typename OneDModelSolver<PARAM, FLUX, SOURCE>::Vec2D
OneDModelSolver<PARAM, FLUX, SOURCE>::BCValuesInternalRight() const
{
    Vec2D temp(2);
    //    for( UInt i=0; i<2; ++i )
    temp[0] = _M_U_thistime[0]( RightInternalNodeId() );
    temp[1] = _M_U_thistime[1]( RightInternalNodeId() );
    return temp;
}


//! set the Dirichlet boundary conditions (right)
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::setBCValuesRight( const Real& bcR1, const Real& bcR2 )
{
    _M_bcDirRight[0] = bcR1; // _M_bcDirRight.first  = bcR1;
    _M_bcDirRight[1] = bcR2; // _M_bcDirRight.second = bcR2;
}


//! set the Dirichlet boundary conditions (left)
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::setBCValuesLeft( const Real& bcL1, const Real& bcL2 )
{
    _M_bcDirLeft[0] = bcL1; // _M_bcDirLeft.first  = bcL1;
    _M_bcDirLeft[1] = bcL2;
}


template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::setBCLeft_internalnode()
{
    _M_bcH.setBCLeft_internalnode();
}


template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::setBCRight_internalnode()
{
    _M_bcH.setBCRight_internalnode();
}


//! simple cfl computation (correct for constant mesh)
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::CheckCFL() const
{
    Real CFL = 0.;

    //! length of the first edge (arbitrary as they are all supposed equal).
    Real deltaX,
        deltaX_min = _M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x();
    Real Ainode, Qinode;

    Real lambda1_max = 0.;
    Real lambda2_max = 0.;
    Real eigval1, eigval2;
    Real tmp11, tmp12, tmp21, tmp22;

    for ( UInt inode=0; inode < _M_dimDof ; inode++ ) {

        Ainode = _M_U_thistime[0]( inode );
        Qinode = _M_U_thistime[1]( inode );

        //! compute the eigenvalues at node
        _M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
                                                 eigval1, eigval2,
                                                 tmp11, tmp12,
                                                 tmp21, tmp22,
                                                 inode );


        lambda1_max = std::max<Real>( std::fabs(eigval1), lambda1_max );
        lambda2_max = std::max<Real>( std::fabs(eigval2), lambda2_max );
    }

    for ( UInt inode=1; inode < _M_dimDof ; inode++ ) {

        deltaX = _M_mesh.edgeList( inode ).pt2().x() - _M_mesh.edgeList( inode ).pt1().x();

        deltaX_min = std::min<Real>( std::fabs(deltaX), deltaX_min );
    }

    CFL = _M_time_step / deltaX_min * std::max<Real>( lambda1_max , lambda2_max );

    if ( _M_CFL )
        std::cout << "CFL = " << CFL << std::endl;

    /*
      std::cout << "Old CFL = " << _M_time_step /
      (_M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x())
      * std::max<Real>( lambda1_max , lambda2_max )
      << std::endl;
    */
    //    }
    if( CFL > 0.5774 )
        std::cout << "\n[CheckCFL] CFL not respected in " << _M_post_file
                  << ": CFL = " << CFL << std::endl;
    //      ASSERT( CFL < 0.5774 , "CFL not respected" );
}


/*! Axpy product for 2D vectors (pairs)
  Axpy(alpha, x, beta, y) -> y = a*A*x + beta*y

  A is given by two pairs corresponding to the 2 lines.
  A = [line1;
  line2 ]
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::Axpy(const Vec2D& line1, const Vec2D& line2,
                                           const Real& alpha,  const Vec2D& x,
                                           const Real& beta,   Vec2D& y) const
{
    ASSERT_PRE( line1.size() == 2 && line2.size(),
                "Axpy works only for 2x2 matrices");
    ASSERT_PRE( x.size() == 2 && y.size(),
                "Axpy works only with 2D vectors");
    y[0]  = alpha * ( line1[0] * x[0] + line1[1] * x[1] ) + beta * y[0];
    y[1] = alpha * ( line2[0] * x[0] + line2[1] * x[1]) + beta * y[1];
}


//! update the P1 flux vector from U: _M_Fluxi = F_h(Un) i=1,2
//! BEWARE: works only for P1Seg elements
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateFlux()
{
    Real Aii, Qii;

    for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
        Aii = _M_U_thistime[0]( ii );
        Qii = _M_U_thistime[1]( ii );
        _M_Flux[0]( ii ) = _M_fluxFun( Aii, Qii, 1, ii );
        _M_Flux[1]( ii ) = _M_fluxFun( Aii, Qii, 2, ii );
    }
}


/*! call _updateFlux and update the P0 derivative of flux vector from U:
  _M_diffFluxij = dF_h/dU(Un) i,j=1,2

  _M_diffFluxij(elem) = 1/2 [ dF/dU(U(node1(elem))) + dF/dU(U(node2(elem))) ]

  (mean value of the two extremal values of dF/dU)

  BEWARE: works only for P1Seg elements
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateFluxDer()
{
    //! first update the Flux vector
    _updateFlux();

    //! then update the derivative of the Flux vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < _M_nb_elem ; ielem++ ) {
        //! for P1Seg and appropriate mesh only!
        ii = ielem;        //! left node of current element
        iip1 = ielem + 1;  //! right node of current element

        Aii = _M_U_thistime[0]( ielem );
        Qii = _M_U_thistime[1]( ielem );
        Aiip1 = _M_U_thistime[0]( ielem + 1 );
        Qiip1 = _M_U_thistime[1]( ielem + 1 );

        for( UInt ii=1; ii<3; ++ii ) {
            for( UInt jj=1; jj<3; ++jj ) {
                tmp = _M_fluxFun.diff(   Aii,   Qii, ii, jj, ii );
                tmp += _M_fluxFun.diff( Aiip1, Qiip1, ii, jj, iip1 );
                _M_diffFlux[ 2*(ii-1) + jj-1 ]( ielem ) = 0.5 * tmp;
            }
        }
    }
}

//! update the P1 source vector from U: _M_Sourcei = S_h(Un) i=1,2
//! BEWARE: works only for P1Seg elements
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateSource( )
{
    Real Aii, Qii;

    for ( UInt ii=0; ii < _M_dimDof ; ii++ ) {
        Aii = _M_U_thistime[0]( ii );
        Qii = _M_U_thistime[1]( ii );
        for(UInt k=0; k<2; ++k)
            _M_Source[k]( ii ) = _M_sourceFun( Aii, Qii, k+1, ii );
    }
}

/*! call _updateSource and update the P0 derivative of source vector from U:
  _M_diffSrcij = dS_h/dU(Un) i,j=1,2

  _M_diffSrcij(elem) = 1/2 [ dS/dU(U(node1(elem))) + dS/dU(U(node2(elem))) ]

  (mean value of the two extremal values of dS/dU)

  BEWARE: works only for P1Seg elements
*/
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::_updateSourceDer( )
{
    //! first update the Source vector
    _updateSource();

    //! then update the derivative of the Source vector
    Real tmp;
    Real Aii, Qii;
    Real Aiip1, Qiip1;
    UInt ii, iip1;

    for ( UInt ielem=0; ielem < _M_nb_elem ; ielem++ ) {
        //! for P1Seg and appropriate mesh only!
        ii = ielem;        //! left node of current element
        iip1 = ielem + 1;  //! right node of current element

        Aii = _M_U_thistime[0]( ielem );
        Qii = _M_U_thistime[1]( ielem );
        Aiip1 = _M_U_thistime[0]( ielem + 1 );
        Qiip1 = _M_U_thistime[1]( ielem + 1 );

        for( UInt ii=1; ii<3; ++ii ) {
            for( UInt jj=1; jj<3; ++jj ) {
                tmp =  _M_sourceFun.diff(   Aii,   Qii, ii, jj, ii );
                tmp += _M_sourceFun.diff( Aiip1, Qiip1, ii, jj, iip1 );
                _M_diffSrc[ 2*(ii-1) + jj-1 ]( ielem ) = 0.5 * tmp;
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
template< class PARAM, class FLUX, class SOURCE >
ScalUnknown<Vector>
OneDModelSolver<PARAM, FLUX, SOURCE>::_correct_flux_inertial( const ScalUnknown<Vector>& flux )
{
    TriDiagMatrix<Real> _matrixLHS(_M_dimDof);
    TriDiagMatrix<Real> _stiffRHS(_M_dimDof);

    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _tridiagsolver(_M_dimDof);

    ElemMat _elmatMassLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatStiffLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatStiffRHS (_M_fe.nbNode,1,1);

    ScalUnknown<Vector> _rhs(_M_dimDof);

    Real _coeffMass;
    Real _coeffStiff;

    Real m, meanA0;

    //    std::ostringstream output;

    _matrixLHS.zero();
    _stiffRHS.zero();

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        _coeffMass = _M_rhs[0]( iedge - 1 ) + _M_rhs[0]( iedge );
        _coeffMass /= 2;
        _coeffMass = 1./_coeffMass;

        meanA0 = _M_oneDParam.Area0(iedge-1) + _M_oneDParam.Area0(iedge);
        meanA0 /= 2;
        m = _M_oneDParam.DensityWall() * _M_oneDParam.Thickness() /
            ( 2 * std::sqrt(4*std::atan(1)) * std::sqrt(meanA0) );
        _coeffStiff = m / _M_oneDParam.DensityRho();

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
        //std::cout << _M_fe.currentId() << std::endl;

        mass( _coeffMass, _elmatMassLHS, _M_fe,0, 0 );
        stiff( _coeffStiff, _elmatStiffLHS, _M_fe,0, 0 );
        stiff( -_coeffStiff, _elmatStiffRHS, _M_fe,0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, _M_fe, _M_dof1D , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, _M_fe, _M_dof1D , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, _M_fe, _M_dof1D , 0, 0 );

        Debug( 6310 ) << "\n\tm = " << m
                      << "\n\t_coeffMass = " << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } //! end loop on elements

    // update rhs
    _stiffRHS.Axpy( 1., flux , 0., _rhs );

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * _M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    Debug( 6310 ) << "\n[_correct_flux_inertial] _matrix.Diag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.Diag()(i);
    Debug( 6310 ) << "\n[_correct_flux_inertial] _matrix.UpDiag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_correct_flux_inertial] _matrix.LowDiag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    _tridiagsolver.Factor( _matrixLHS );

    Debug( 6310 ) << "\n[_correct_flux_inertial] factorized _matrix.Diag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.Diag()(i);
    Debug( 6310 ) << "\n[_correct_flux_inertial] factorized _matrix.UpDiag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_correct_flux_inertial] factorized _matrix.LowDiag():\n";
    for(int i=0; i<_matrixLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _matrixLHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    Debug( 6310 ) << "\n[_correct_flux_inertial] solveing with rhs:\n";
    for(int i=0; i<_matrixLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _rhs(i);
    Debug( 6310 ) << "\n";

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _matrixLHS, _rhs );

    return _rhs;
}

/*
  L2 Projection of the second derivative of Q over P1 space
*/
template< class PARAM, class FLUX, class SOURCE >
ScalUnknown<Vector>
OneDModelSolver<PARAM, FLUX, SOURCE>::_compute_d2Q_dx2( const ScalUnknown<Vector>& flux )
{
    TriDiagMatrix<Real> _massLHS(_M_dimDof);
    TriDiagMatrix<Real> _stiffRHS(_M_dimDof);

    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _tridiagsolver(_M_dimDof);

    ElemMat _elmatMassLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatStiffRHS (_M_fe.nbNode,1,1);

    ScalUnknown<Vector> _rhs(_M_dimDof);

    _massLHS.zero();
    _stiffRHS.zero();

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffRHS.zero();

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
        //std::cout << _M_fe.currentId() << std::endl;

        mass( 1., _elmatMassLHS, _M_fe,0, 0 );
        stiff( -1., _elmatStiffRHS, _M_fe,0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, _M_fe, _M_dof1D , 0, 0 );
        assemb_mat( _stiffRHS, _elmatStiffRHS, _M_fe, _M_dof1D , 0, 0 );

    } //! end loop on elements

    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _stiffRHS.Diag():\n";
    for(int i=0; i<_stiffRHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _stiffRHS.Diag()(i);
    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _stiffRHS.UpDiag():\n";
    for(int i=0; i<_stiffRHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _stiffRHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _stiffRHS.LowDiag():\n";
    for(int i=0; i<_stiffRHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _stiffRHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    // update rhs
    _stiffRHS.Axpy( 1., flux , 0., _rhs );

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * _M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.Diag():\n";
    for(int i=0; i<_massLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.UpDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_compute_d2Q_dx2] _massLHS.LowDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    _tridiagsolver.Factor( _massLHS );

    Debug( 6310 ) << "\n[_compute_d2Q_dx2] solving with rhs:\n";
    for(int i=0; i<_massLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _rhs(i);
    Debug( 6310 ) << "\n";

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
template< class PARAM, class FLUX, class SOURCE >
ScalUnknown<Vector>
OneDModelSolver<PARAM, FLUX, SOURCE>::_correct_flux_viscoelastic( const ScalUnknown<Vector>& flux )
{
    TriDiagMatrix<Real> _matrixLHS(_M_dimDof);
    TriDiagMatrix<Real> _stiffRHS(_M_dimDof);

    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _tridiagsolver(_M_dimDof);

    ElemMat _elmatMassLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatStiffLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatStiffRHS (_M_fe.nbNode,1,1);

    ScalUnknown<Vector> _rhs(_M_dimDof);

    Real _coeffMass;
    Real _coeffStiff;

    Real gamma, meanA0(1.);

    _matrixLHS.zero();
    _stiffRHS.zero();

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatStiffLHS.zero();
        _elmatStiffRHS.zero();

        // this comes from the exact derivation of generalized string model
        // + Voigt viscoelasticity
        _coeffMass = _M_rhs[0]( iedge - 1 ) + _M_rhs[0]( iedge );
        _coeffMass *= 0.5;
        _coeffMass = 1./ std::sqrt(_coeffMass);

        if(_M_linearize_string_model) {
            // this is to recover the linearized version (_coeffMass = 1/A)
            _coeffMass *= _coeffMass;

            meanA0 = _M_oneDParam.Area0( iedge - 1 ) + _M_oneDParam.Area0( iedge );
            meanA0 *= 0.5;

            if(_M_linearize_equations) {
                // when using linearized equations, A \simeq A0
                _coeffMass = 1./ meanA0;
            }
        }
        gamma = _M_oneDParam.Gamma() / ( 2 * std::sqrt(4*std::atan(1)) );
        gamma *= 1 / std::sqrt( meanA0 );
        _coeffStiff = _M_time_step * gamma / _M_oneDParam.DensityRho();

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));
        //std::cout << _M_fe.currentId() << std::endl;

        mass( _coeffMass, _elmatMassLHS, _M_fe,0, 0 );
        stiff( _coeffStiff, _elmatStiffLHS, _M_fe,0, 0 );
        stiff( -_coeffStiff, _elmatStiffRHS, _M_fe,0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _matrixLHS, _elmatMassLHS, _M_fe, _M_dof1D , 0, 0 );
        assemb_mat( _matrixLHS, _elmatStiffLHS, _M_fe, _M_dof1D , 0, 0 );

        assemb_mat( _stiffRHS, _elmatStiffRHS, _M_fe, _M_dof1D , 0, 0 );

        Debug( 6310 ) << "\n\tgamma = " << gamma
                      << "\n\t_coeffMass = " << _coeffMass
                      << "\n\t_coeffStiff = " << _coeffStiff << "\n";

    } //! end loop on elements

    // update rhs
    // rhs = _stiffRHS * rhs
    _stiffRHS.Axpy( 1., flux , 0., _rhs );

    // NOTE I should add to rhs the value of boundary integral ("neumann-like")
    // BUT it's useless since I'm going to impose dirichlet conditions on all boundaries

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    // _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * _M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    // but nothing has to be added because of homogeneous dirichlet bc
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.second;

    _updateBCDirichletMatrix(_matrixLHS);

    _tridiagsolver.Factor( _matrixLHS );

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _matrixLHS, _rhs );

    return _rhs;
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
template< class PARAM, class FLUX, class SOURCE >
ScalUnknown<Vector>
OneDModelSolver<PARAM, FLUX, SOURCE>::_correct_flux_longitudinal( )
{
    TriDiagMatrix<Real> _massLHS(_M_dimDof);
    TriDiagMatrix<Real> _massRHS(_M_dimDof);

    TriDiagCholesky< Real, TriDiagMatrix<Real>, Vector > _tridiagsolver(_M_dimDof);

    ElemMat _elmatMassLHS (_M_fe.nbNode,1,1);
    ElemMat _elmatMassRHS (_M_fe.nbNode,1,1);

    ScalUnknown<Vector> _rhs(_M_dimDof);
    // let g = sqrt(A) - sqrt(A0)
    // now f = _d3g_dz3
    ScalUnknown<Vector> _g(_M_dimDof);
    ScalUnknown<Vector> _f(_M_dimDof);

    //          _g = _M_rhs[0];
    for( UInt i=0; i<_M_dimDof; ++i )
        _g(i) = std::sqrt(_M_rhs[0](i)) - std::sqrt(_M_oneDParam.Area0(i));

    UInt inode;

    Real _coeffMassLHS;
    Real _coeffMassRHS;

    Real _a;
    //    std::ostringstream output;
    _massLHS.zero();
    _massRHS.zero();

    Real _h( _M_oneDParam.Length() / static_cast<Real>(_M_oneDParam.ParamSize() - 1) );

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for(UInt iedge = 1; iedge <= _M_mesh.numEdges(); iedge++){

        inode = iedge - 1;

        //! set the elementary matrices to 0.
        _elmatMassLHS.zero();
        _elmatMassRHS.zero();

        // coeff (1/A) (average values over the element)
        _coeffMassLHS = _M_rhs[0]( inode ) + _M_rhs[0]( inode+1 );
        _coeffMassLHS /= 2;
        _coeffMassLHS = 1./_coeffMassLHS;

        _a = _M_oneDParam.CoeffA() / std::sqrt(4*std::atan(1));
        _coeffMassRHS = _M_time_step * _a / _M_oneDParam.DensityRho();

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
        else if(inode>_M_mesh.numEdges()-2) { // forward differentiation
            _f( inode ) = _g( inode ) - _g( inode-3 )
                + 3*_g( inode-2 ) - 3*_g( inode-1 );
            Debug( 6310 ) << "\n\forward differentiation = " << _coeffMassLHS << "\n"; }
        else  { // central differentiation
            _f( inode ) = -_g( inode - 2 ) + 2*_g( inode-1 )
                - 2*_g( inode+1 ) + _g( inode+2 );
            Debug( 6310 ) << "\n\tcentral differentiation = " << _coeffMassLHS << "\n"; }

        _f(inode) *= 1/(2*std::pow(_h,3));

        //! update the current element
        _M_fe.updateFirstDerivQuadPt(_M_mesh.edgeList(iedge));

        mass( _coeffMassLHS, _elmatMassLHS, _M_fe,0, 0 );
        mass( _coeffMassRHS, _elmatMassRHS, _M_fe,0, 0 );

        //! assemble the mass and grad matrices
        assemb_mat( _massLHS, _elmatMassLHS, _M_fe, _M_dof1D , 0, 0 );
        assemb_mat( _massRHS, _elmatMassRHS, _M_fe, _M_dof1D , 0, 0 );

        Debug( 6310 ) << "\n\t_coeffMassLHS = " << _coeffMassLHS << "\n";
        Debug( 6310 ) << "\n\t_coeffMassRHS = " << _coeffMassRHS << "\n";

    } //! end loop on elements

    // update rhs
    //    _massLHS.Axpy( 1., _M_U_thistime[4] , 0., _rhs );
    _massRHS.Axpy( 1., _f , 1., _rhs );

    UInt firstDof = 0;
    UInt lastDof  = _rhs.size()-1;

    //! symmetric treatment (cholesky can be used)
    //! first row modified (Dirichlet)
    _rhs( firstDof ) = 0.;
    //! second row modified (for symmetry)
    //    _rhs( firstDof + 1 ) += - _matrix.LowDiag()( firstDof ) * _M_bcDirLeft.second;

    //! last row modified (Dirichlet)
    _rhs( lastDof ) = 0.;
    //! penultimate row modified (for symmetry)
    //    _rhs( lastDof - 1 ) += - _matrix.UpDiag()( lastDof - 1 ) * _M_bcDirRight.second;

    _updateBCDirichletMatrix(_massLHS);

    Debug( 6310 ) << "\n[_correct_flux_longitudinal] _matrix.Diag():\n";
    for(int i=0; i<_massLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
    Debug( 6310 ) << "\n[_correct_flux_longitudinal] _matrix.UpDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_correct_flux_longitudinal] _matrix.LowDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    _tridiagsolver.Factor( _massLHS );

    Debug( 6310 ) << "\n[_correct_flux_longitudinal] factorized _matrix.Diag():\n";
    for(int i=0; i<_massLHS.OrderMatrix(); ++i)
        Debug( 6310 ) << "\t" << _massLHS.Diag()(i);
    Debug( 6310 ) << "\n[_correct_flux_longitudinal] factorized _matrix.UpDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.UpDiag()(i);
    Debug( 6310 ) << "\n[_correct_flux_longitudinal] factorized _matrix.LowDiag():\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _massLHS.LowDiag()(i);
    Debug( 6310 ) << "\n";

    Debug( 6310 ) << "\n[_correct_flux_longitudinal] solveing with rhs:\n";
    for(int i=0; i<_massLHS.OrderMatrix()-1; ++i)
        Debug( 6310 ) << "\t" << _rhs(i);
    Debug( 6310 ) << "\n";

    //! cholesky or lapack lu solve
    //! solve the system: rhs1 = massFactor^{-1} * rhs1
    _tridiagsolver.Solve( _massLHS, _rhs );

    return _rhs;
}


//! Initialize from Reimann invariants
//! Initialize with constant initial conditions
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(const Real& u10, const Real& u20,
    const std::string& var )
{

    if( var == "physical")
        {
            //        for( UInt i=0; i<2; ++i )
            _M_U_thistime[0] = ScalarVector( _M_dimDof, u10 );
            _M_U_thistime[1] = ScalarVector( _M_dimDof, u20 );

            for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
                _M_oneDParam.W_from_U( _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                                       _M_U_thistime[0][ielem], _M_U_thistime[1][ielem],
                                       ielem );
        }
    else if( var == "reimann" )
        {
            //        for( UInt i=0; i<2; ++i )
            _M_U_thistime[2] = ScalarVector( _M_dimDof, u10 );
            _M_U_thistime[3] = ScalarVector( _M_dimDof, u20 );

            for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
                _M_oneDParam.U_from_W( _M_U_thistime[0][ielem], _M_U_thistime[1][ielem],
                                       _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                                       ielem );
        }
    else
        {
            std::cout << "[initialize] trying to initialize " << var << " variables!"
            << std::endl;
            abort();
        }

    for( UInt i=0; i<2; ++i )
        {
            _M_U_prevtime[i] = _M_U_thistime[i];
            _M_U_2prevtime[i] = _M_U_prevtime[i];
        }

    Vector pressures(4 * _M_dimDof);
    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
        subrange(pressures, 4*ielem, 4+4*ielem) =
            _M_oneDParam.pressure( _M_U_thistime[0][ielem],
                                   _M_U_prevtime[0][ielem], _M_U_2prevtime[0][ielem],
                                   _M_time_step, ielem,
                                   _M_dP_dt_steps, _M_viscoelastic_wall,
                                   _M_linearize_string_model );

        _M_U_thistime[4][ielem] = pressures(4*ielem);
    }

    if(_M_viscoelastic_wall) {
        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            _M_U_thistime[_M_variable_index_map.find("P_elast")->second][ielem]
                                                                         = pressures(1+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("P_visc")->second][ielem]
                                                                        = pressures(2+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("dA_dt")->second][ielem]
                                                                       = pressures(3+4*ielem);
        }
    }
    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step );

    //! create matlab scripts
    create_movie_file();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

}

//! Initialize with vectors containing all the nodal values
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(const Vector& u10, const Vector& u20)
{
    _M_U_thistime[0] = u10;
    _M_U_thistime[1] = u20;

    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
        _M_oneDParam.W_from_U( _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                               _M_U_thistime[0][ielem], _M_U_thistime[1][ielem], ielem );

    for( UInt i=0; i<2; ++i )
        {
            _M_U_prevtime[i] = _M_U_thistime[i];
            _M_U_2prevtime[i] = _M_U_prevtime[i];
        }

    Vector pressures(4 * _M_dimDof);
    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
        subrange(pressures, 4*ielem, 4+4*ielem) =
            _M_oneDParam.pressure( _M_U_thistime[0][ielem],
                                   _M_U_prevtime[0][ielem], _M_U_2prevtime[0][ielem],
                                   _M_time_step, ielem,
                                   _M_dP_dt_steps, _M_viscoelastic_wall,
                                   _M_linearize_string_model );

        _M_U_thistime[4][ielem] = pressures(4*ielem);
    }

    if(_M_viscoelastic_wall) {
        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            _M_U_thistime[_M_variable_index_map.find("P_elast")->second][ielem] =
              pressures(1+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("P_visc")->second][ielem] =
              pressures(2+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("dA_dt")->second][ielem] =
              pressures(3+4*ielem);
        }
    }
    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step );

    //! create matlab scripts
    create_movie_file();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

}

//! Initialize when initial conditions concentration
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(const Function& /*c0*/,
    Real /*t0*/, Real /*dt*/)
{
    ERROR_MSG("Not yet implemented");
}

// ! Initialize when initial values for the concentration are read from file
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(const std::string & vname)
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
    Resfile.read((char*)&_M_U_thistime[0](1),_M_U_thistime[0].size()*sizeof(Real));
    Resfile.close();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

}

//! Initialize only Flux (Area read from OneDNonLinParam)
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(const Real& u20)
{
    _M_U_thistime[1] = ScalarVector( _M_dimDof, u20 );

    Vector pressures(4 * _M_dimDof);

    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
        _M_U_thistime[0][ielem]=_M_oneDParam.Area0(ielem);
        for( UInt i=0; i<2; ++i )
            {
                _M_U_prevtime[i][ielem]=_M_U_thistime[i][ielem];
                _M_U_2prevtime[i][ielem]=_M_U_prevtime[i][ielem];
            }
        _M_oneDParam.W_from_U( _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                               _M_U_thistime[0][ielem], _M_U_thistime[1][ielem], ielem );

        //      for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
        subrange(pressures, 4*ielem, 4+4*ielem) =
            _M_oneDParam.pressure( _M_U_thistime[0][ielem],
                                   _M_U_prevtime[0][ielem], _M_U_2prevtime[0][ielem],
                                   _M_time_step, ielem,
                                   _M_dP_dt_steps, _M_viscoelastic_wall,
                                   _M_linearize_string_model );

        _M_U_thistime[4][ielem] = pressures(4*ielem);
    }

    if(_M_viscoelastic_wall) {
        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            _M_U_thistime[_M_variable_index_map.find("P_elast")->second][ielem] =
              pressures(1+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("P_visc")->second][ielem] =
              pressures(2+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("dA_dt")->second][ielem] =
              pressures(3+4*ielem);
        }
    }

    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step );

    //! create matlab scripts
    create_movie_file();

    openFileBuffers();

    output2FileBuffers( 0. );

    postProcess( 0. );

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
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::initialize(GetPot& data_file)
{
    // the discontinuity is comprised between firstnode and lastnode
    UInt firstnode( data_file("initialize/firstnode",1) );
    UInt lastnode( data_file("initialize/lastnode",2) );

    ASSERT_PRE( (firstnode <= lastnode) && (lastnode <= _M_dimDof),
                "[initialize] outside tube boundaries" );

    // read initialization type from data file (see OneDModelSolver::initialize)
    std::string init_var( data_file("initialize/var","P") );
    Real multiplier( data_file("initialize/multiplier",1.) );

    // tell me what I am doing
    Debug( 6310 ) << "[initialize] 0- Initializing with values:\n";
    Debug( 6310 ) << "[initialize]\t\tinitialize var = " << init_var << "\n";
    Debug( 6310 ) << "[initialize]\t\tfirstnode = " << firstnode << "\n";
    Debug( 6310 ) << "[initialize]\t\tlastnode = " << lastnode << "\n";
    Debug( 6310 ) << "[initialize]\t\tmultiplier = " << multiplier << "\n";

    Real value1, value2, width( data_file("initialize/width",5.) );

    Vector exponent(_M_dimDof);
    for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
        {
            exponent[inode] *= 0.;

            // first half of a gaussian signal, centered in firstnode;
            // second half of a gaussian signal, centered in lastnode;
            // width represents the total duration of the gaussian signal
            // (rise + decay)
            if( (inode < firstnode) || (inode > lastnode) ) {
                exponent[inode] = - std::pow( double(int( inode - firstnode )), 2 );
                exponent[inode] /= 2 * std::pow( width, 2 );
            }
        }

    switch( _M_oneDstring2initializeVarMap[init_var] )
        {
            // case 1, 2: initialize physical variables to desired value
        case OneDInitPressure:
            // this is a pressure value! has to be converted in area value
            value1 = data_file("initialize/rest_value",0.);
            // HYPOTHESIS: when initializing pressure, flux is imposed constant = 0
            value2 = 0;

            Debug( 6310 ) << "[initialize] pressure " << value1 << "\n";

            _M_U_thistime[1] = ScalarVector( _M_dimDof, value2 );
            Debug( 6310 ) << "[initialize] Q done\n";

            for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
                {
                    // reusing value2 as help variable
                    value2 = value1 * ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_U_thistime[0][inode] = _M_oneDParam.A_from_P( value2 );
                    Debug( 6310 ) << "[initialize] A(" << inode <<") done\n";
                    _M_oneDParam.W_from_U( _M_U_thistime[2][inode], _M_U_thistime[3][inode],
                                           _M_U_thistime[0][inode], _M_U_thistime[1][inode],
                                           inode );
                    Debug( 6310 ) << "[initialize] W_i(" << inode <<") done\n";
                }
            break;

        case OneDInitArea:
            value1 = data_file("initialize/value",0.);
            value2 = 0;

            _M_U_thistime[1] = ScalarVector( _M_U_thistime[1].size(), value2 );

            for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
                {
                    _M_U_thistime[0][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_oneDParam.W_from_U( _M_U_thistime[2][inode], _M_U_thistime[3][inode],
                                           _M_U_thistime[0][inode], _M_U_thistime[1][inode],
                                           inode );
                }
            break;

        case OneDInitFlux:
            // HYPOTHESIS: when initializing flux, area is equal to Area0
            //      value1 = _M_oneDParam.Area0(0); // this if Area0 is constant
            value2 = data_file("initialize/value",0.);

            for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
                {
                    _M_U_thistime[0][inode] = _M_oneDParam.Area0(inode);
                    _M_U_thistime[1][inode] = value2 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_oneDParam.W_from_U( _M_U_thistime[2][inode], _M_U_thistime[3][inode],
                                           _M_U_thistime[0][inode], _M_U_thistime[1][inode],
                                           inode );
                }
            break;

        case OneDInitReimann1:
            value1 = data_file("initialize/value",0.);
            value2 = -value1;
            std::cout << "[initialize] WARNING! Initializing W2 = - W1"
                      << " (assuming Q = 0)" << std::endl;

            for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
                {
                    _M_U_thistime[2][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_U_thistime[3][inode] = value2 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_oneDParam.U_from_W( _M_U_thistime[0][inode],
                                           _M_U_thistime[1][inode],
                                           _M_U_thistime[2][inode],
                                           _M_U_thistime[3][inode],
                                           inode );
                }

            break;

        case OneDInitReimann2:
            value1 = data_file("initialize/value",0.);
            value2 = -value1;
            std::cout << "[initialize] WARNING! Initializing W1 = - W2"
                      << " (assuming Q = 0)" << std::endl;

            for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
                {
                    _M_U_thistime[3][inode] = value1 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_U_thistime[2][inode] = value2 *
                      ( 1 + multiplier * std::exp( exponent[inode] ) );
                    _M_oneDParam.U_from_W( _M_U_thistime[0][inode],
                                           _M_U_thistime[1][inode],
                                           _M_U_thistime[2][inode],
                                           _M_U_thistime[3][inode],
                                           inode );
                }


            break;

        default:
            ERROR_MSG("No such initializing option.");

        }

    Debug( 6310 ) << "[initialize]\t\tvalue1 = " << value1 << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue1_step = " << value1 * multiplier << "\n";
    Debug( 6310 ) << "[initialize]\t\tvalue2 = " << value2 << "\n";

    for( UInt i=0; i<2; ++i )
        {
            _M_U_prevtime[i] = _M_U_thistime[i];
            _M_U_2prevtime[i] = _M_U_prevtime[i];
        }

    Vector pressures(4 * _M_dimDof);
    for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
        subrange(pressures, 4*ielem, 4+4*ielem) =
            _M_oneDParam.pressure( _M_U_thistime[0][ielem],
                                   _M_U_prevtime[0][ielem], _M_U_2prevtime[0][ielem],
                                   _M_time_step, ielem,
                                   _M_dP_dt_steps, _M_viscoelastic_wall,
                                   _M_linearize_string_model );

        _M_U_thistime[4][ielem] = pressures(4*ielem);
    }

    if(_M_viscoelastic_wall) {
        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            _M_U_thistime[_M_variable_index_map.find("P_elast")->second][ielem] =
              pressures(1+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("P_visc")->second][ielem] =
              pressures(2+4*ielem);
            _M_U_thistime[_M_variable_index_map.find("dA_dt")->second][ielem] =
              pressures(3+4*ielem);
        }
    }

    //! Prepare bc Handler
    _M_bcH.setDefaultBC( _M_mesh, _M_sourceFun, _M_time_step );

    //! create matlab scripts
    create_movie_file();

    //! Prepare the buffers
    openFileBuffers();

    //! Write down initial condition
    output2FileBuffers( 0. );
    postProcess( 0. );

}


// Remember the values at previous timestep
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::savesol()
{
    //    for( UInt i=0; i<2; ++i )
    _M_U_prevtime[0] = _M_U_thistime[0];
    _M_U_prevtime[1] = _M_U_thistime[1];
}

// Recover the solution at previous timestep
// BUT keep unaltered the boundary values
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::loadsol()
{
    Vec2D U_leftbd( 2 ), U_rightbd( 2 );
    for( UInt i=0; i<2; ++i ) {
        U_leftbd[i] = _M_U_thistime[i]( _M_leftNodeId );
        U_rightbd[i] = _M_U_thistime[i]( _M_rightNodeId );
    }

    //    for( UInt i=0; i<2; ++i )
    _M_U_thistime[0] = _M_U_prevtime[0];
    _M_U_thistime[1] = _M_U_prevtime[1];

    for( UInt i=0; i<2; ++i ) {
        _M_U_thistime[i]( _M_leftNodeId ) = U_leftbd[i];
        _M_U_thistime[i]( _M_rightNodeId ) = U_rightbd[i];
    }

    for (UInt inode=_M_leftNodeId; inode <= _M_rightNodeId ; ++inode )
        _M_oneDParam.W_from_U( _M_U_thistime[2][inode], _M_U_thistime[3][inode],
                               _M_U_thistime[0][inode], _M_U_thistime[1][inode],
                               inode );

}


//! Update the right hand side for time advancing
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::timeAdvance( const Real& time_val )
{
    Chrono chrono;

    Real dt2over2 = _M_time_step * _M_time_step * 0.5;

    Debug( 6310 ) << "[timeAdvance] o-  updates of flux and sources... " << "\n";
    chrono.start();

    //! output cfl
    CheckCFL();

    //! update the vector containing the values of the flux at the nodes
    //! and its jacobian
    _updateFluxDer();
    //! update the vector containing the values of the source term at the nodes
    //! and its jacobian
    _updateSourceDer();
    chrono.stop();
    Debug( 6310 ) << "[timeAdvance] \tdone in " << chrono.diff() << " s.\n";


    Debug( 6310 ) << "[timeAdvance] o-  updates of matrices... " << "\n";
    chrono.start();
    //! update the matrices for the non-linear terms
    _updateMatrices();
    chrono.stop();
    Debug( 6310 ) << "[timeAdvance] \tdone in " << chrono.diff() << " s.\n";

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

    Debug( 6310 ) << "[timeAdvance] o-  Matrix vector products... " << "\n";
    chrono.start();

    //! Reminder of the function Axpy:
    //! Axpy(alpha, x, beta, y) -> y = alpha*A*x + beta*y

    //!--------------------------------------------------------
    //! 1/, 2/ compute _M_rhs[0], _M_rhs[1] (systems in U1=A, U2=Q)
    //!--------------------------------------------------------
    for( UInt i=0; i<2; ++i ) {
        //! rhs = mass * Un
        _M_massMatrix.Axpy( 1., _M_U_thistime[i] , 0., _M_rhs[i] );

        //! rhs = rhs + dt * grad * F(Un)
        _M_gradMatrix.Axpy( _M_time_step, _M_Flux[i] , 1., _M_rhs[i] );

        //! rhs = rhs - dt * mass * S(Un)
        _M_massMatrix.Axpy( - _M_time_step, _M_Source[i] , 1., _M_rhs[i] );

        for( UInt j=0; j<2; ++j ) {
            //! rhs = rhs - dt^2/2 * gradDiffFlux * S(Un)
            _M_gradMatrixDiffFlux[2*i+j].Axpy( -dt2over2, _M_Source[j] , 1., _M_rhs[i] );

            //! rhs = rhs + dt^2/2 * divDiffSrc * F(Un)
            _M_divMatrixDiffSrc[2*i+j].Axpy( dt2over2, _M_Flux[j] , 1., _M_rhs[i] );

            //! rhs = rhs - dt^2/2 * stiffDiffFlux * F(Un)
            _M_stiffMatrixDiffFlux[2*i+j].Axpy( -dt2over2, _M_Flux[j] , 1., _M_rhs[i] );

            //! rhs = rhs + dt^2/2 * massDiffSrc * S(Un)
            _M_massMatrixDiffSrc[2*i+j].Axpy( dt2over2, _M_Source[j] , 1., _M_rhs[i] );

            Debug( 6310 ) << "[timeAdvance] \tfilling position " << 2*i+j << "\n";
        }
    }

    Debug( 6310 ) << "[timeAdvance] \tcomputed rhs\n";
    //!---------------------------------------------------
    //! 3/ take into account the BOUNDARY CONDITIONS
    //!---------------------------------------------------
    //! compute the values for the boundary conditions
    _computeBC( time_val );
    //! take into account the bc
    _updateBCDirichletVector();

    // *******************************************************
    chrono.stop();
    Debug( 6310 ) << "[timeAdvance] \trhs computed in " << chrono.diff() << " s.\n";


}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::iterate( const Real& time_val , const int& count)
{
    Debug( 6310 ) << "[iterate] o-  Solving the system... t = " << time_val
                  << ", iter = " << count  << "... \n";
    Chrono chrono;
    chrono.start();

    if( _M_UW ){
        Real Ainode, Qinode;

        Real lambda1_plus = 0.;
        Real lambda2_plus = 0.;
        Real lambda1_minus = 0.;
        Real lambda2_minus = 0.;
        Real eigval1, eigval2;
        Real tmp11, tmp12, tmp21, tmp22;

        Real deltaX = _M_mesh.edgeList( 1 ).pt2().x() - _M_mesh.edgeList( 1 ).pt1().x();

        //! working on riemann invariants
        ScalUnknown<Vector> W1_UW(_M_dimDof);
        ScalUnknown<Vector> W2_UW(_M_dimDof);

        //! converting boundary conditions on physical variables
        //! in boundary conditions on characteristic variables
        _M_oneDParam.W_from_U( W1_UW(0), W2_UW(0),
                               _M_rhs[0][0], _M_rhs[1][0], 0 );
        _M_oneDParam.W_from_U( W1_UW(_M_dimDof-1), W2_UW(_M_dimDof-1),
                               _M_rhs[0][_M_dimDof-1], _M_rhs[1][_M_dimDof-1],
                               _M_dimDof-1 );

        for ( UInt ii=1; ii < (_M_dimDof-1) ; ii++ ) {
            //! compute the eigenvalues at node
            Ainode = _M_U_thistime[0]( ii );
            Qinode = _M_U_thistime[1]( ii );
            _M_fluxFun.jacobian_EigenValues_Vectors( Ainode, Qinode,
                                                     eigval1, eigval2,
                                                     tmp11, tmp12,
                                                     tmp21, tmp22,
                                                     ii );

            lambda1_plus = std::max<Real>( eigval1, 0. );
            lambda1_minus = std::min<Real>( eigval1, 0. );
            lambda2_plus = std::max<Real>( eigval2, 0. );
            lambda2_minus = std::min<Real>( eigval2, 0. );
            //! update the solution for the next time step
            W1_UW[ii] = _M_U_thistime[2][ii]
                - (_M_time_step / deltaX) * lambda1_plus * ( _M_U_thistime[2][ii] -
                    _M_U_thistime[2][ii-1])
                - (_M_time_step / deltaX) * lambda1_minus * ( _M_U_thistime[2][ii+1] -
                    _M_U_thistime[2][ii])
                - _M_time_step * ( tmp11 * _M_Source[0][ii] + tmp12 * _M_Source[1][ii] );
            W2_UW[ii] = _M_U_thistime[3][ii]
                - (_M_time_step / deltaX) * lambda2_plus * ( _M_U_thistime[3][ii] -
                    _M_U_thistime[3][ii-1])
                - (_M_time_step / deltaX) * lambda2_minus * ( _M_U_thistime[3][ii+1] -
                    _M_U_thistime[3][ii])
                - _M_time_step * ( tmp21 * _M_Source[0][ii] + tmp22 * _M_Source[1][ii] );
        }

        _M_U_thistime[2] = W1_UW;
        _M_U_thistime[3] = W2_UW;

        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ )
            _M_oneDParam.U_from_W( _M_U_thistime[0][ielem], _M_U_thistime[1][ielem],
                                   _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                                   ielem );

    }
    else{
        //! cholesky or lapack lu solve
        //! solve the system: rhs1 = massFactor^{-1} * rhs1
        _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs[0] );

        //! solve the system: rhs2 = massFactor^{-1} * rhs2
        _M_tridiagSlv.Solve( _M_factorMassMatrix, _M_rhs[1] );

        //! correct flux with inertial term
        if( _M_inertial_wall ) {
            _M_U_thistime[_M_variable_index_map.find("Q_inert")->second] =
              _correct_flux_inertial( _M_rhs[1] );
            _M_rhs[1] += _M_U_thistime[_M_variable_index_map.find("Q_inert")->second];
        }
        //! correct flux with viscoelastic term
        if( _M_viscoelastic_wall ) {
            _M_U_thistime[_M_variable_index_map.find("Q_visc")->second] =
              _correct_flux_viscoelastic( _M_rhs[1] );
            _M_rhs[1] += _M_U_thistime[_M_variable_index_map.find("Q_visc")->second];
        }
        //! compute L2 projection of d2Q_dx2
        //    if( _M_flux_second_der )
        //      _M_d2_U2_dx2 = _compute_d2Q_dx2( _M_rhs[1] );

        //! correct flux with longitudinal term
        if( _M_longitudinal_wall ) {
            _M_U_thistime[_M_variable_index_map.find("Q_long")->second] =
              _correct_flux_longitudinal(  );
            _M_rhs[1] += _M_U_thistime[_M_variable_index_map.find("Q_long")->second];
        }

        //! store solution at previous timesteps & update the solution for the next time step
        for( UInt i=0; i<2; ++i )
            {
                _M_U_2prevtime[i] = _M_U_prevtime[i];
                _M_U_prevtime[i] = _M_U_thistime[i];
                _M_U_thistime[i] = _M_rhs[i];
            }

        //      std::cout << std::endl;

        Vector pressures(4 * _M_dimDof);
        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            _M_oneDParam.W_from_U( _M_U_thistime[2][ielem], _M_U_thistime[3][ielem],
                                   _M_U_thistime[0][ielem], _M_U_thistime[1][ielem],
                                   ielem );
            //        for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
            subrange(pressures, 4*ielem, 4+4*ielem) =
                _M_oneDParam.pressure( _M_U_thistime[0][ielem],
                                       _M_U_prevtime[0][ielem], _M_U_2prevtime[0][ielem],
                                       _M_time_step, ielem,
                                       _M_dP_dt_steps, _M_viscoelastic_wall,
                                       _M_linearize_string_model );

            _M_U_thistime[4][ielem] = pressures(4*ielem);
            //        std::cout << 4*ielem << "\t";
        }
        //      std::cout << std::endl;

        if(_M_viscoelastic_wall) {
            for (UInt ielem=0; ielem <= _M_nb_elem ; ielem++ ) {
                _M_U_thistime[_M_variable_index_map.find("P_elast")->second][ielem] =
                  pressures(1+4*ielem);
                _M_U_thistime[_M_variable_index_map.find("P_visc")->second][ielem] =
                  pressures(2+4*ielem);
                _M_U_thistime[_M_variable_index_map.find("dA_dt")->second][ielem] =
                  pressures(3+4*ielem);
            }
        }

    }


    chrono.stop();
    Debug( 6310 ) << "[iterate] \tdone in " << chrono.diff() << " s.\n";

    output2FileBuffers( time_val );

}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::openFileBuffers()
{
    std::string file_output;

    boost::shared_ptr<std::ostringstream> buf;

    _M_variable_index_iter iter_variable;
    _M_variable_filter_iter iter_suffix;

    for( iter_variable = _M_variable_index_map.begin();
         iter_variable != _M_variable_index_map.end(); ++iter_variable ) {

        for( iter_suffix = _M_variable_filter_map.begin();
             iter_suffix != _M_variable_filter_map.end(); ++iter_suffix ) {

            file_output = _M_post_dir + "/" + _M_post_file +
              iter_variable->first + iter_suffix->first;
            Debug( 6310 ) << "\n[openFileBuffers] setting output for file "
            << file_output << "\n";

            buf.reset( new std::ostringstream("") );
            buf->setf(std::ios_base::scientific);
            buf->precision(5);
            buf->width(13);
            _M_post_process_buffer.insert( std::map<std::string,
                boost::shared_ptr<std::ostringstream> >::value_type( file_output, buf ) );
            _M_post_process_buffer_offset.insert( std::map<std::string,
                long>::value_type( file_output, buf->tellp() ) );
        }
    }
}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::output2FileBuffers( const Real& time_val )
{

    if( !( static_cast<int>( std::floor( (time_val/_M_time_step) + .5 ) ) %
        _M_verbose ) ){

        std::string file_output;

        _M_variable_index_iter iter_variable;
        _M_variable_filter_iter iter_suffix;

        for( iter_variable = _M_variable_index_map.begin();
             iter_variable != _M_variable_index_map.end(); ++iter_variable ) {

            for( iter_suffix = _M_variable_filter_map.begin();
                 iter_suffix != _M_variable_filter_map.end(); ++iter_suffix ) {

                file_output = _M_post_dir + "/" + _M_post_file +
                  iter_variable->first + iter_suffix->first;
                Debug( 6310 ) << "\n[output2FileBuffers] setting output for file "
                              << file_output << ", writing variable "
                              << iter_variable->first << ";\n" << "variable size = "
                              << _M_U_thistime[iter_variable->second].size();

                (this->*(iter_suffix->second))( file_output, time_val,
                    _M_U_thistime[iter_variable->second], iter_variable->first );
            }
        }
    }

}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::closeFileBuffers()
{
    // as I have a boost::shared_ptr, I expect the objects to be deallocated
    // now that the pointers are destroyed
    _M_post_process_buffer.erase( _M_post_process_buffer.begin(),
        _M_post_process_buffer.end() );
    _M_post_process_buffer_offset.erase( _M_post_process_buffer_offset.begin(),
        _M_post_process_buffer_offset.end() );
}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::postProcess( const Real& time_val )
{

    std::string str;

    std::ofstream outfile;

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    if (time_val==0.){
        // the code is entering this for sure, as
        // initialize invokes postProcess( 0. )

        std::string file_output;

        Real deltax( _M_oneDParam.Length() /
            static_cast<Real>(_M_oneDParam.ParamSize() - 1) );

        file_output = "dati.m";
        outfile.open(file_output.c_str(), std::ios::app );
        outfile << "z" << _M_post_file
                << " = (" << _M_x_left
                << ":" << deltax
                << ":" << _M_x_right << ");\n"
                << std::endl;
        outfile.close();

    }


    Debug( 6310 ) << "[postProcess] o- Dumping solutions on files (1d)!" << "\n";
    // dump solutions on files (buffers must be active!)
    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end();
         iter++ ){
        outfile.open( (*iter).first.c_str(), std::ios::app );
        outfile << (*(*iter).second).str();
        outfile.close();
    }


    resetFileBuffers();


};



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::resetFileBuffers( )
{
    closeFileBuffers();
    openFileBuffers();
}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::seekpFileBuffers( )
{

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end();
         iter++ ){
        (*iter).second->seekp( _M_post_process_buffer_offset[(*iter).first] );

    }

}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::tellpFileBuffers( )
{

    std::map< std::string, boost::shared_ptr<std::ostringstream> >::iterator iter;

    for( iter = _M_post_process_buffer.begin(); iter != _M_post_process_buffer.end();
         iter++ ){
        _M_post_process_buffer_offset[(*iter).first] = (*iter).second->tellp();

    }

}



template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::gplot( )
{
    _M_GracePlot.Plot( _M_mesh.pointList(), _M_U_thistime[0] );
}

//! output for Plotmtv.
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::output_to_plotmtv(std::string fname, Real time_val,
                                                        const ScalVec& U,
                                                        std::string vname )
{

    boost::shared_ptr<std::ostringstream> buf;

    buf = _M_post_process_buffer[fname];

    (*buf) << "$ DATA = CURVE2D\n % xlabel='z'\n"
           << "% toplabel='Section,time=" << time_val << "'\n % ylabel='" << vname << "'\n"
           << std::flush;

    for(UInt ii = 0; ii < U.size(); ii++){
        (*buf) << _M_mesh.pointList()[ii].x() << " " << U[ii] << "\n" << std::flush;
    }


}

//! output for Matlab.
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::output_to_matlab( std::string fname, Real /*time_val*/,
                                                        const ScalVec& U,
                                                        std::string /*vname*/ )
{

    boost::shared_ptr<std::ostringstream> buf;

    buf = _M_post_process_buffer[fname];

    //    (*buf) << vname << _M_post_file
    //           << "( " << (static_cast<int>( std::floor( time_val/_M_time_step + 0.5 ) )
    //              /  _M_verbose )+1
    //           << ", : ) = [ " << std::flush;

    for ( UInt ii=LeftNodeId(); ii <= RightNodeId() ; ++ii ) {
        //        (*buf) << U[ii] << "; ";
        (*buf) << U[ii] << " ";
    }
    //    (*buf) << "]';\n" << std::endl;
    (*buf) << std::endl;


}

//! Create a Matlab script to visualize output matlab files
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::create_movie_file()
{
    std::ofstream outfile;
    std::string file_output;
    file_output = _M_post_dir + "/" + "areamovie"+_M_post_file+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Area"<<_M_post_file<<";\n"
            <<"[n,m]=size(A" << _M_post_file << ");\n"
            <<"Max=max(max(A" << _M_post_file << "));\n"
            <<"Min=min(min(A" << _M_post_file << "));\n"
            <<"for i=1:n\n"
            <<"  plot(A" << _M_post_file << "(i,:));\n"
            <<"  title((i-1)*"<<_M_time_step<<"*"<<_M_verbose<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";

    outfile.close();
    file_output = _M_post_dir + "/" + "portatamovie"+_M_post_file+".m";
    outfile.open(file_output.c_str(), std::ios::out);
    outfile <<"Portata"<<_M_post_file<<";\n"
            <<"[n,m]=size(Q" << _M_post_file << ");\n"
            <<"Max=max(max(Q" << _M_post_file << "));\n"
            <<"Min=min(min(Q" << _M_post_file << "));\n"
            <<"for i=1:n\n"
            <<"  plot(Q" << _M_post_file << "(i,:));\n"
            <<"  title((i-1)*"<<_M_time_step<<"*"<<_M_verbose<<");\n"
            <<"  axis([0 m Min-(Min/100) Max+(Min/100)]);\n"
            <<"  %pause;\n"
            <<"  F(i) = getframe;\n"
            <<"end\n"
            <<"movie(F)";
    outfile.close();

}


//! Print to screen information on the Solver class
template< class PARAM, class FLUX, class SOURCE >
void
OneDModelSolver<PARAM, FLUX, SOURCE>::showMe(std::ostream& c, UInt verbose)
{
    c << "\n--- One Dimensional Model Data\n";
    this->showMeData(c);
    c << "\n--- One Dimensional Model Handler\n";
    this->showMeHandler(c, verbose);
    c << "\n--- One Dimensional Model Parameters\n";
    this->oneDParam().showMeData(c);
    c << "--- End of One Dimensional Model\n";

}


}

#endif
