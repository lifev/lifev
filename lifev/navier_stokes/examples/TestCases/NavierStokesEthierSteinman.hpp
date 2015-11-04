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
    @file NavierStokesEthierSteinman.hpp
    @brief Navier-Stokes cavity problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 19-08-2011
 */

#ifndef NAVIERSTOKESETHIERSTEINMAN_HPP
#define NAVIERSTOKESETHIERSTEINMAN_HPP

#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>

namespace LifeV
{

namespace
{
typedef RegionMesh< LinearTetra > esMesh_Type;
}

class NavierStokesEthierSteinman : public NavierStokesProblem< esMesh_Type >
{

private:

    // Walls
    static const Int LEFTWALL;
    static const Int RIGHTWALL;
    static const Int FRONTWALL;
    static const Int BACKWALL;
    static const Int TOPWALL;
    static const Int BOTTOMWALL;
    // Edges
    static const Int BOTTOMEDGE1;
    static const Int BOTTOMEDGE2;
    static const Int BOTTOMEDGE3;
    static const Int BOTTOMEDGE4;
    static const Int SIDEEDGE1;
    static const Int SIDEEDGE2;
    static const Int SIDEEDGE3;
    static const Int SIDEEDGE4;
    static const Int TOPEDGE1;
    static const Int TOPEDGE2;
    static const Int TOPEDGE3;
    static const Int TOPEDGE4;
    // Corners
    static const Int BOTTOMCORNER1;
    static const Int BOTTOMCORNER2;
    static const Int BOTTOMCORNER3;
    static const Int BOTTOMCORNER4;
    static const Int TOPCORNER1;
    static const Int TOPCORNER2;
    static const Int TOPCORNER3;
    static const Int TOPCORNER4;

public:

    //! @name Constructors, destructor
    //@{

    //! Empty constructor
    NavierStokesEthierSteinman();

    //! Destructor
    virtual ~NavierStokesEthierSteinman();

    //@}

    //! @name  Methods
    //@{

    //! Returns true if the problem has an exact solution
    bool hasExactSolution() const;

    //! Returns the value of the exact solution
    NavierStokesProblem< esMesh_Type >::function_Type xexact();

    //! Returns the value of the exact solution (velocity components only)
    NavierStokesProblem< esMesh_Type >::function_Type uexact();

    //! Returns the value of the derivative of the exact solution with respect to the time (velocity components only)
    NavierStokesProblem< esMesh_Type >::function_Type uderexact();

    //! Returns the value of the exact solution (pressure component only)
    NavierStokesProblem< esMesh_Type >::function_Type pexact();

    //! Display general information about the problem
    /*!
        @param output specify the output stream (std::cout by default)
     */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name  Set Methods
    //@{

    //! Setup the boundary conditions in the provided BCHandler
    /*!
        @param bcHandler shared pointer on a BCHandler object
     */
    void setBoundaryConditions ( std::shared_ptr<BCHandler> bcHandler ) const;

    //! Set the viscosity of the fluid
    /*!
        @param viscosity viscosity of the fluid
     */
    void setViscosity ( const Real& viscosity );

    //! Set the density of the fluid
    /*!
        @param density density of the fluid
     */
    void setDensity ( const Real& density );

    //@}

    //! @name  Get Methods
    //@{

    //! Getter for the problem mesh
    void mesh ( std::shared_ptr< esMesh_Type >& meshPart ) const;

    //! Getter for the boundary conditions in the provided BCHandler
    /*!
        @param bcHandler shared pointer on a BCHandler object
     */
    void boundaryConditions ( std::shared_ptr<BCHandler> bcHandler ) const;

    //! Returns the name of the problem
    std::string name() const;

    //@}

};

} // namespace LifeV

#endif /* NAVIERSTOKESETHIERSTEINMAN_HPP */
