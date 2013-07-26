/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2005-04-19

  Copyright (C) 2005 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file cylinder.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2005-04-19
 */

#ifndef _RESISTANCEFSI_HPP
#define _RESISTANCEFSI_HPP

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
/*!
 * \class Outflow
 * \brief 2D/3D Cylinder Simulation class
 *
 *  @author Christophe Prud'homme
 *  @see
 */

namespace LifeV
{

class ResistanceBCs
{
public:

    ResistanceBCs();

    void initParameters      ( const int outflowFlag, const Real resistance, const Real hydrostatic, const std::string name );

    void renewParameters     ( OseenSolverShapeDerivative<RegionMesh<LinearTetra> >&  solver, const VectorEpetra& solution, const Real time );

    Real computeResistance   ( const Real time );

    Real fZero               ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/ );

    static Real outPressure0   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure1   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure2   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure3   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure4   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure5   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);
    static Real outPressure6   (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/);

    inline Real resistance( )
    {
        return M_resistance;
    }

    inline Real flow( )
    {
        return M_outflux;
    }

    inline std::string name( )
    {
        return M_name;
    }

    inline Real hydrostatic( )
    {
        return M_hydrostaticP;
    }

    inline Real outP( )
    {
        return M_outP;
    }

private:
    Real pi;

    Real M_outflux;
    Real M_resistance;
    Real M_hydrostaticP;
    Real M_outP;
    int  M_flag;
    std::string M_name;

    UInt conditionNumber;

    static std::vector<Real> outputVector;
};
}
#endif /* __RESISTANCE_H */
