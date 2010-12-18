//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief File containing a class for handling temporal discretization with GetPot.
 *
 *  @date 01-10-2003
 *  @author Miguel Angel Fernandez
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _DATANEWTON_H_
#define _DATANEWTON_H_
#include <string>
#include <iostream>

#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{
/*!
  \class DataNewton

  Base class which holds data concerning Newton method

*/
class DataNewton
{
public:

  //! @name Constructor & Destructor
  //@{
  //! Constructor
  DataNewton( const GetPot& dfile, const std::string& section = "newton" );

  //! Virtual destructor
  virtual ~DataNewton();
  //@}

  //! @name Methods
  //@{

  //! Ouptut
  virtual void showMe( std::ostream& c = std::cout ) const;

  //@}

  //! @name Get Methods
  //@{

  //! The max number of interations
  const UInt getMaxiter() const;

  //! The absolute tolerance
  const Real getAbstol() const;

  //! The relative tolerance
  const Real getReltol() const;

  //! The maximum error tolerance for residual in linear solver.
  const Real getEtamax() const;

  //! The linesearch option
  const UInt getLinesearch() const;

  //@}

  /*const UInt __attribute__ ((__deprecated__)) maxiter()
  {
    return getMaxiter();
  }

  const UInt __attribute__ ((__deprecated__)) abstol()
  {
    return getAbstol();
  }

  const UInt __attribute__ ((__deprecated__)) reltol()
  {
    return getReltol();
  }

  const UInt __attribute__ ((__deprecated__)) etamax()
  {
    return getEtamax();
  }

  const UInt __attribute__ ((__deprecated__)) linesearch()
  {
    return getLinesearch();
    }*/

protected:

  UInt M_maxiter;      // max number of iterations
  Real M_abstol;       // the stopping criteria is abstol+reltol*norm(residual_0)
  Real M_reltol;       //
  Real M_etamax;       // Maximum error tolerance for residual in linear solver.
  // The linear solver terminates when the relative
  // linear residual is smaller than eta*| f(sol) |.
  // The value linear_rel_tol send for the relative tolerance
  // to the linear solver is therefore eta. eta is determined
  // by the modified Eisenstat-Walker formula if etamax > 0.
  // If eta_max < 0, then eta = |etamax| for the entire
  // iteration (e.g. etamax = -1e-6 ensures that the linear
  // tolerance would be always 1e-6).
  UInt M_linesearch;   // 0 (no linesearch) 1 (parabolic) 2 (cubic: recommended)
};
}
#endif
