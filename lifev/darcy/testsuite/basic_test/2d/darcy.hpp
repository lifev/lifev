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
/**
     @file darcy.hpp
    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
    @date 2012-06-13
*/

#ifndef __darcy_H
#define __darcy_H 1


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

/*!
    @class darcy_nonlinear
    @brief LifeV non-linear and transient Darcy test case

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    A 2D Darcy transient and non-linear test with Dirichlet boundary condition.
    Let us set \f$ \Omega = (0,1)^3 \f$ and \f$ I_T = \Omega \times [0,4)\f$, solve the problem in dual-mixed form
    \f[
        \left\{
        \begin{array}{l l l }
            \Lambda^{-1} \left( p \right) \sigma + \nabla p = f_v & \mathrm{in} & I_T,  \vspace{0.2cm} \\
            \displaystyle \phi \frac{\partial p}{\partial t} + \nabla \cdot \sigma + \pi p - f = 0
            & \mathrm{in} & I_T,  \vspace{0.2cm} \\
            p = g_D                              & \mathrm{on} & \partial \Omega \times [0,4),\vspace{0.2cm} \\
            p = p_0                              & \mathrm{in} & \Omega.
          \end{array}
        \right.
    \f]
    Furthermore the data are
    \f[
        \begin{array}{l l l}
            f_v(x,y) = \left( ty - tx, -t^2 y^2 \right)^T &
            \pi(x,y) = 1 &
            f(x,y) = tx^2 - 2t^2 - t - (6y + 2 t^2 y)(1 + t^4 x^4 + y^6 + 2 t^2 x^2 y^3)
            - (3y^2 + t^2 y^2)(6y^5 + 6t^2 x^2 y^2) + t^2 x^2 + y^3, \vspace{0.2cm} \\
            g_D(x,y,z) = t^2 x^2 + y^3, &
            \phi(x,y) = 0.5, &
            \Lambda (x,y,p) =
            \left(
            \begin{array}{c c}
                1 & 0 \\
                0 & 1+p^2
            \end{array}
            \right), \\
            p_0(x,y) = y^3. & &
        \end{array}
    \f]
    The analytical solutions are
    \f[
        p(x,y,z) = t^2 x^2 + y^3,
        \quad
        \sigma(x,y,z) =
        \left(
        \begin{array}{c}
            - 2t^2 x + t (y - x) \\
            - (3y^2 + t^2 y^2)( 1 + t^4 x^4 + y^6 +2 t^2 x^2 y^3)
        \end{array}
        \right).
    \f]
    Fixing \f$ \Delta t = 0.015625 \f$ and the tolerance of the fixed point scheme at
    \f$ 10^{-10} \f$, the computed space errors are
    <table border="1" align="center">
        <tr>
            <th> N </th>
            <th> \f$ \left\Vert p - p_h \right\Vert_{L^2} \f$ </th>
            <th> \f$ \left\Vert \sigma - \Pi_{P_0} \sigma_h \right\Vert_{L^2} \f$ </th>
            <th> N proc </th>
        </tr>
        <tr>
            <td> 10 </td>
            <td> 0.105324 </td>
            <td> 0.284968 </td>
            <td> 8 </td>
        </tr>
        <tr>
            <td> 20 </td>
            <td> 0.0375415 </td>
            <td> 0.140317 </td>
            <td> 8 </td>
        </tr>
        <tr>
            <td> 40 </td>
            <td> 0.0162273 </td>
            <td> 0.066338 </td>
            <td> 16 </td>
        </tr>
        <tr>
            <td> 80 </td>
            <td> 0.00773563 </td>
            <td> 0.0315562 </td>
            <td> 16 </td>
        </tr>
        <tr>
            <td> 160 </td>
            <td> 0.00381641 </td>
            <td> 0.015142 </td>
            <td> 32 </td>
        </tr>
    </table>
    where N is the number of subdivisions for each boundary.
    @image html darcy/2d.png "Example of the solution with N = 10."
*/
class darcy_nonlinear
//     :
//     public LifeV::Application
{
public:

    /* @name Constructors and destructor
     */
    //@{

    //! Constructor
    darcy_nonlinear ( int argc, char** argv );


    //! Destructor
    ~darcy_nonlinear () {}

    //@}

    /* @name  Methods
     */
    //@{

    //! To lunch the simulation
    LifeV::Real run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __darcy_H */
