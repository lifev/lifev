/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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
   \file main_lin.cpp
   \author Tiziano Passerini <tiziano@mathcs.emory.edu>
   \date 2009-04-12
   \brief Sets up a test case for the 1D  solver (for a linearized version
      of the Euler equations)

  A linearized version of the Euler equations an be derived as follows:
  -# we neglect the nonlinear term, therefore
     \f$ \alpha \frac{\partial}{\partial x} \left( \frac{Q^{2}}{A} \right) = 0 \f$;
  -# we linearize the coefficients with respect to \f$ A \f$, setting
     \f$ A(x, t) \simeq A_0 \f$.

  The resulting linear system of first order partial differential
  equations reads:
    \f[
      \frac{\partial A}{\partial t} + \frac{\partial Q}{\partial x} = 0
    \f]
    \f[
      \frac{\partial Q}{\partial t} + \frac{A_0}{\rho} \frac{\partial P}{\partial x}
      + \frac{K_R}{A_0} Q = 0
      \ .
    \f]

  We will consider a viscoelastic structural model for the relation between
  pressure and the cross-sectional area, therefore
  \f[
  \frac{\partial P}{\partial x} =
     a \frac{\partial A}{\partial x}
     + \gamma \frac{\partial^2 A}{\partial x \partial t}
  \ ,
  \f]
  with \f$ a = \frac{ \tilde{a} }{2 \sqrt{\pi A_0}} \f$,
  \f$ \tilde{a} = \frac{\beta}{A_0} \sqrt{\pi} \f$,
  \f$ \beta = \frac{\sqrt{\pi} h_0 E}{1-\xi^2} \f$,
  \f$ \gamma \f$ is called visco-elastic modulus and these
  parameters are assumed constant along x.

  The purely elastic case is recovered when \f$ \gamma = 0 \f$.

  Under the hypothesis of periodicity of the solution, we can look for
  solutions in the form of harmonic waves:
  \f{equation}
    \Re \bigl( A(x,t) \bigr) = \hat{A}( k ) \exp \bigl[ \Im (k) x
    \bigr] \cos \bigl( \omega(k) t - \Re (k) x \bigr)
  \f}
  \f{equation}
    \Re \bigl( Q(x,t) \bigr) = \exp \bigl[ \Im (k) x \bigr] \Bigl(
      \Re (\hat{Q}) \cos \bigl( \omega t - \Re (k) x \bigr) -
      \Im (\hat{Q}) \sin \bigl( \omega t - \Re (k) x \bigr) \Bigr)
    \ .
  \f}
  \f$ i \f$ being the imaginary unit, \f$ k \f$ the wave number,
  \f$ \omega \f$ the angular frequency,
  \f$ \hat{A} \f$ and \f$ \hat{Q} \f$ the area and flow rate wave amplitudes
  at \f$ (x, t) = (0, 0) \f$. We will consider
  \f$ \hat{A} \in \mathbf{R} \f$, \f$ \omega \in \mathbf{R} \f$,
  \f$ k \in \mathbf{C} \f$, \f$ \hat{Q} \in \mathbf{C} \f$,

  \sa T. Passerini, "Computational hemodynamics of the cerebral circulation:
  multiscale modeling from the circle of Willis to cerebral aneurysms", PhD Thesis,
  2009
 */
#include <life/lifecore/life.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifecore/GetPot.hpp>

#include <life/lifefem/sobolevNorms.hpp>

#include <life/lifesolver/oneDNonLinModelParam.hpp>
#include <life/lifesolver/vectorFunction1D.hpp>
#include <life/lifesolver/oneDModelSolver.hpp>

#include "ud_functions.hpp"

#include <sstream>

// *********************************
// Useful typedefs
// *********************************
typedef LifeV::LinearSimpleFluxFun1D<LifeV::LinearizedParam> Flux1D;
typedef LifeV::LinearSimpleSourceFun1D<LifeV::LinearizedParam> Source1D;
typedef LifeV::OneDModelSolver<LifeV::LinearizedParam, Flux1D, Source1D> onedsolver_type;
typedef boost::shared_ptr<LifeV::Analytical_Solution> analytical_sol_ptrtype;

// *********************************
// Forward declarations
// *********************************
/*!
  \brief Evaluate the relative error with respect to the analytical solution

  \param[in] onedsolver the 1D solver
  \param[in] sol_vec the vector of pointers to the analytic solutions
  \param[in] t the time
  \return a vector the same size as sol_vec, containing the relative errors
    of the computed solution wrt the respective analytical solutions
 */
LifeV::Vector computeRelativeError(onedsolver_type& onedsolver,
    std::vector<analytical_sol_ptrtype> const& sol_vec,
    LifeV::Real const& t);


int main(int argc, char** argv)
{
  using namespace LifeV;

  // *********************************
  // ***** Reading from data file
  // *********************************
  GetPot command_line(argc,argv);
  std::string data_file_name = command_line.follow("data_lin", 2, "-f","--file");
  GetPot data_file(data_file_name);

  // *********************************
  // Useful variables
  // *********************************
  // abscissa on which to evaluate the analytical solution
  Real _x_current( 0. );
  // relative error (one component for each available analytical solution)
  Vector relative_error;
  // timer
  Chrono chrono;
  // iteration counter
  int count(0);

  // parameters for the analytical solution
  Real kappa_Re( data_file("initialize/kappa_Re", 8.*std::atan(1.) ) );
  Real kappa_Im( data_file("initialize/kappa_Im", 8.*std::atan(1.) ) );
  Real omega( data_file("initialize/omega", 8.*std::atan(1.) ) );
  Real omega_over_kappa_Re( data_file("initialize/omega_over_kappa_Re", 0. ) );
  Real omega_over_kappa_Im( data_file("initialize/omega_over_kappa_Im", 0. ) );

  // for debugging purposes
  Debug() << "\n[main] omega: " << omega;
  Debug() << "\n[main] kappa_Re: " << kappa_Re;
  Debug() << "\n[main] kappa_Im: " << kappa_Im;
  Debug() << "\n[main] omega_over_kappa_Re: " << omega_over_kappa_Re;
  Debug() << "\n[main] omega_over_kappa_Im: " << omega_over_kappa_Im;
  Debug() << "\n";


  // Temporal parameters
  Real postprocess_dt = data_file("miscellaneous/postprocess_timestep",0.01);

  // Maximum relative error (see the exit condition for the test)
  Real _max_rel_err(0.);


  // *********************************
  // Build the (linear) 1D model
  // *********************************
  onedsolver_type onedm(data_file);
  onedm.showMe( std::cout );

  Real dt     = onedm.timestep();
  Real startT = onedm.inittime();
  Real T      = onedm.endtime();

  Debug() << "[main] startT T dt postprocess_dt "
  << startT << " " <<  T << " " << dt << " " << postprocess_dt << "\n";


  // *********************************
  // Boundary conditions
  // *********************************
  Real analytical_A_amplitude( onedm.oneDParam().Area0(0) );

  /*
    from the continuity equation it follows that
    \hat{Q} = \hat{A} * \omega / kappa
    so if we choose \hat{A} in \mathbb{R}, then the complex part of \hat{Q}
    is given by the term omega_over_kappa_Im
   */
  Real analytical_Q_amplitude_Re( analytical_A_amplitude
      * omega_over_kappa_Re);
  Real analytical_Q_amplitude_Im( analytical_A_amplitude
      * omega_over_kappa_Im);

  analytical_sol_ptrtype inflow;

  inflow.reset( new Analytical_Solution( analytical_Q_amplitude_Re,
      analytical_Q_amplitude_Im, kappa_Re, kappa_Im, omega ) );

  analytical_sol_ptrtype outflow;

  outflow.reset( new Analytical_Solution( analytical_Q_amplitude_Re,
      analytical_Q_amplitude_Im, kappa_Re, kappa_Im, omega ) );
  outflow->update_x( onedm.xRight() );

  onedm.bcH().setBC( inflow, "left", "first", "Q" );
  onedm.bcH().setBC( outflow, "right", "first", "Q" );

  onedm.showMeData();


  // *********************************
  // The analytical solution
  // *********************************
  // How many variables are we checking?
  UInt nVar(2);
  // we exploit the functions implemented in the BC classes
  std::vector<analytical_sol_ptrtype> analytical_sol_vec(nVar);
  // in this test we control area and flux
  analytical_sol_vec[0].reset( new Analytical_Solution( analytical_A_amplitude,
      0., kappa_Re, kappa_Im, omega ) );
  analytical_sol_vec[1].reset( new Analytical_Solution( analytical_Q_amplitude_Re,
      analytical_Q_amplitude_Im, kappa_Re, kappa_Im, omega ) );
  // prepare the output data structure
  relative_error.resize(nVar);

  // *********************************
  // Initialize with analytical solution
  // *********************************
  Real deltax( ( onedm.xRight() - onedm.xLeft() ) / onedm.nbElem() );
  Vector U10(onedm.nbElem()+1), U20(onedm.nbElem()+1);

  // evaluate the analytical solution on the mesh vertices
  for( UInt icomp=0; icomp < onedm.nbElem()+1; ++icomp )
    {
      _x_current = icomp * deltax;

      analytical_sol_vec[0]->update_x( _x_current );
      U10[icomp] = analytical_sol_vec[0]->evaluate(0.);

      analytical_sol_vec[1]->update_x( _x_current );
      U20[icomp] = analytical_sol_vec[1]->evaluate(0.);
    }

  onedm.initialize(U10, U20);

  // *********************************
  // Temporal loop
  // *********************************
  printf("\nTemporal loop:\n");

  for (Real time=startT+dt ; time <= T; time+=dt) {
    count++;
    chrono.start();

    onedm.timeAdvance( time );
    onedm.iterate( time , count );

    relative_error = computeRelativeError(onedm, analytical_sol_vec, time);

    if( !( static_cast<int>( std::floor( time/dt + 0.5 ) ) %
        static_cast<int>( std::floor( postprocess_dt/dt  + 0.5 ) ) ) )
      onedm.postProcess( time );

    if ( data_file( "miscellaneous/show_graceplot", 0 ) )
      onedm.gplot();

    chrono.stop();

    printf("\033[0GIter %d", count);
    printf("\033[10G, t = %f", time);
    printf("\033[25Gs (%f", chrono.diff());
    printf("\033[35Gs). Rel.Err: A = %f, ", relative_error[0]);
    printf("\033[60G, Q = %f.", relative_error[1]);

    // keep trace of the performances of the solver
    // store the biggest error
    _max_rel_err = (norm_inf(relative_error) > _max_rel_err) ?
        norm_inf(relative_error) : _max_rel_err;
  }

  printf("\nSimulation ended successfully.\n");

  // the threshold is obviously customizable
  return (_max_rel_err > 1);
}


LifeV::Vector computeRelativeError(onedsolver_type& onedsolver,
    std::vector<analytical_sol_ptrtype> const& sol_vec,
    LifeV::Real const& t)
{
  // How many variables are we considering?
  LifeV::UInt nVar(sol_vec.size());
  ASSERT_PRE(nVar <= onedsolver.U_thistime().size(),
      "You have more analytical solutions than unknowns!")
  // elementwise error
  std::vector<LifeV::Vector> _err_vec(nVar, LifeV::ZeroVector(onedsolver.nbElem()+1));
  // (P1) interpolate of the analytical solution
  std::vector<LifeV::Vector> _sol_interp_vec(nVar, LifeV::ZeroVector(onedsolver.nbElem()+1));
  // norms and errors
  LifeV::Vector _err_norm = LifeV::ZeroVector(nVar), _sol_norm = LifeV::ZeroVector(nVar),
  _rel_err = LifeV::ZeroVector(nVar);
  // abscissa on which to evaluate the analytical solution
  LifeV::Real _x_current( 0. );
  // delta x
  LifeV::Real deltax( ( onedsolver.xRight() - onedsolver.xLeft() )
      / onedsolver.nbElem() );

  // treat separately the first mesh vertex
  for( LifeV::UInt iVar = 0; iVar < nVar; ++iVar )
    {
      sol_vec[iVar]->update_x( 0. );
      _sol_interp_vec[iVar][0] = sol_vec[iVar]->evaluate(t);
      _err_vec[iVar][0] = _sol_interp_vec[iVar][0] - onedsolver.U_thistime()[iVar][0];
    }

  // cycle over the elements to evaluate the L2 norm of the error
  // (numbering starts from 1)
  for( LifeV::UInt iEl = 1; iEl <= onedsolver.Mesh().numEdges(); ++iEl )
    {
      // update the fe descriptor
      onedsolver.fe().updateJacQuadPt( onedsolver.Mesh().edgeList( iEl ) );
      // update the abscissa
      _x_current = iEl * deltax;
      // update the analytical solution for area
      for( LifeV::UInt iVar = 0; iVar < nVar; ++iVar )
        {
        sol_vec[iVar]->update_x( _x_current );
        _sol_interp_vec[iVar][iEl] = sol_vec[iVar]->evaluate(t);
        _err_vec[iVar][iEl] = _sol_interp_vec[iVar][iEl] - onedsolver.U_thistime()[iVar][iEl];
        // compute the norm on this element and accumulate
        _sol_norm[iVar] += elem_L2_2( _sol_interp_vec[iVar], onedsolver.fe(), onedsolver.dof() );
        // compute the error on this element and accumulate
        _err_norm[iVar] += elem_L2_2( _err_vec[iVar], onedsolver.fe(), onedsolver.dof() );
        }
    }

  // build the relative error vector
  for( LifeV::UInt iVar = 0; iVar < nVar; ++iVar )
    _rel_err[iVar] = sqrt(_err_norm[iVar]/_sol_norm[iVar]);

  return _rel_err;
}
