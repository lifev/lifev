/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
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
/*!
  \file dataNewton.h
  \author M.A. Fernandez
  \date 10/2003 
  \version 1.0

  \brief File containing a class for handling temporal discretization with GetPot

*/
#ifndef _DATANEWTON_H_
#define _DATANEWTON_H_
#include <string>
#include <iostream>

#include "GetPot.hpp"
#include "lifeV.hpp"


/*! 
  \class DataNewton

  Base class which holds data concerning Newton method

*/
class DataNewton
{
 public:

  //! Constructor
  DataNewton(const GetPot& dfile, const std::string& section="newton");
  // The max number of interations

  //! The max number of interations
  UInt maxiter() const;

  //! The absolute tolerance
  Real abstol() const;

  //! The relative tolerance
  Real reltol() const;

  //! The maximum error tolerance for residual in linear solver.
  Real etamax() const;
  
  //! The linesearch option
  UInt linesearch() const;

  //! Ouptut
  virtual void showMe(std::ostream& c = std::cout) const;
 
  //! Virtual destructor
  virtual ~DataNewton();

 protected:
 
  UInt _maxiter;      // max number of iterations 
  Real _abstol;       // the stopping criteria is abstol+reltol*norm(residual_0)
  Real _reltol;       //
  Real _etamax;       // Maximum error tolerance for residual in linear solver.
                      // The linear solver terminates when the relative
		      // linear residual is smaller than eta*| f(sol) |.
		      // The value linear_rel_tol send for the relative tolerance
		      // to the linear solver is therefore eta. eta is determined
		      // by the modified Eisenstat-Walker formula if etamax > 0.
		      // If eta_max < 0, then eta = |etamax| for the entire
		      // iteration (e.g. etamax = -1e-6 ensures that the linear
		      // tolerance would be always 1e-6). 
  UInt _linesearch;   // 0 (no linesearch) 1 (parabolic) 2 (cubic: recommended)
};

#endif
