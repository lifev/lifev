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
    @file
    @brief Contains methods to generate 1D meshes.

    @author -
    @contributor Mauro Perego <mperego@fsu.edu>
    @maintainer -

    @date 10-04-2011
 */

#ifndef REGIONMESH1DBUILDER_HPP
#define REGIONMESH1DBUILDER_HPP 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/RegionMesh1DStructured.hpp>

namespace LifeV
{

//! Build uniform mesh along the x axis.
/**
 *  @param mesh Reference to the mesh
 *  @param x_l Left end point
 *  @param x_r Right end point
 *  @param numberOfElements Number of elements inside the mesh.
 *
 *  Build 1D uniform mesh along the x axis, extending from x_l to x_r, with numberOfElements elements
 */
template <typename MC>
void LIFEV_DEPRECATED ( uniformMesh1D ( RegionMesh<LinearLine, MC>& mesh,
                                        const Real& x_l, const Real& x_r,
                                        const UInt& numberOfElements ) );

template <typename MC>
void uniformMesh1D ( RegionMesh<LinearLine, MC>& mesh,
                     const Real& x_l, const Real& x_r,
                     const UInt& numberOfElements )
{
    regularMesh1D ( mesh, 1, numberOfElements, false, x_r - x_l, x_l );
}

} // Namespace LifeV

#endif /* REGIONMESH1DBUILDER_HPP */
