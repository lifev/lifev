/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
       Date: 2004-11-02

  Copyright (C) 2004 Politecnico di Milano

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
   \file interface2Vessels.hpp
   \author Vincent Martin <vincent.martin@mate.polimi.it>
   \date 2004-11-02
 */

#ifndef _INTERFACE2VESSELS_H_
#define _INTERFACE2VESSELS_H_

#include "oneDModelSolver.hpp"

namespace LifeV
{
/*!
  \class Interface2Vessels

  This class contains an interface to couple 2 vessels
  described by 1D models.

     tube alpha          tube beta
  |---------------|  |-----------------|


*/
class Interface2Vessels
{

    typedef pair< Real, Real > Vec2D;

public:

    Interface2Vessels( OneDModelSolver const& tube_alpha,
                       OneDModelSolver const& tube_beta );

    void updateInterface2Vessels( OneDModelSolver const& tube_alpha,
                                  OneDModelSolver const& tube_beta );
    Vector BcDir_alpha() const;
    Vector BcDir_beta() const;

    int computeInterface2TubesValues();

private:

    //! function computing the function and its gradient
    void f_jac( const Vector& x,
                Vector& f,
                Matrix& jac
                ) const;
    //! 2D dot product
    Real dot(const Vec2D& vec1, const Vec2D& vec2) const;

    Vec2D interpolLinear(const Real& point_bound, const Real& point_internal,
                         const Real& deltaT, const Real& eigenvalue,
                         const Vec2D& U_bound, const Vec2D& U_intern) const;

private:
    //! side : tube alpha
    Vec2D _M_Un_alpha_bd;
    Vec2D _M_Un_alpha_int;
    const Edge1D _M_edge_alpha;
    const UInt _M_dof_alpha;
    const NonLinearFluxFun1D&   _M_fluxFun_alpha;
    const NonLinearSourceFun1D& _M_sourceFun_alpha;
    //! result of the computation
    Vector _M_bcDir_alpha;

    //! side : tube beta
    Vec2D _M_Un_beta_bd;
    Vec2D _M_Un_beta_int;
    const Edge1D _M_edge_beta;
    const UInt _M_dof_beta;
    const NonLinearFluxFun1D&   _M_fluxFun_beta;
    const NonLinearSourceFun1D& _M_sourceFun_beta;
    //! result of the computation
    Vector _M_bcDir_beta;

    //! common values
    const Real _M_time_step;

};

}

#endif
