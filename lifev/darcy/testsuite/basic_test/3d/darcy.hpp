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

#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

/*!
    @class darcy_linear
    @brief LifeV linear Darcy test case

    @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>

    Simple linear 3D Darcy test with Dirichlet, Neumann and Robin Boundary condition.
    <br>
    Solve the problem in dual-mixed form
    \f[
        \left\{
        \begin{array}{l l l }
            \Lambda^{-1} \sigma + \nabla p = f_v & \mathrm{in} & \Omega,  \vspace{0.2cm} \\
            \nabla \cdot \sigma + \pi p - f = 0  & \mathrm{in} & \Omega,  \vspace{0.2cm} \\
            p = g_D                              & \mathrm{on} & \Gamma_D,\vspace{0.2cm} \\
            \sigma \cdot n + h p = g_R           & \mathrm{on} & \Gamma_R, \vspace{0.2cm} \\
            \sigma \cdot n = g_N                 & \mathrm{on} & \Gamma_N,
          \end{array}
        \right.
    \f]
    where \f$ \Omega = (0,1)^3 \f$ with
    \f$ \Gamma_R = \left\{ z = 0 \right\} \f$,
    \f$ \Gamma_{N_1} = \left\{ z = 1 \right\} \f$,
    \f$ \Gamma_{N_2} = \left\{ y = 0 \right\} \f$ and
    \f$ \Gamma_D = \partial \Omega \setminus \left( \Gamma_R \cup \Gamma_{N_1}
    \cup \Gamma_{N_2} \right) \f$.
    Furthermore the data are
    \f[
        \begin{array}{l l l}
            f_v(x,y,z) = \left( x^3, 2y, 4z \right)^T &
            \pi(x,y,z) = xy - 0.5z &
            f(x,y,z) = 4x^2 - 4y^2 - 8xy + 6 + (xy - 0.5z)(x^2y^2 + 6x + 5z), \vspace{0.2cm} \\
            g_D(x,y,z) = x^2y^2 + 6x + 5z, &
            h(x,y,z) = 1, &
            g_R(x,y,z) =  5 - 4z + x^2y^2 + 6x + 5z, \vspace{0.2cm} \\
            g_{N_1}(x,y,z) = - 5 + 4z, &
            g_{N_2}(x,y,z) = 2xy^2 + 6 + 2x^2y - 2y - x^3, &
            \Lambda (x,y,z) =
            \left(
            \begin{array}{c c c}
                2 & 1 & 0 \\
                1 & 1 & 0 \\
                0 & 0 & 1
            \end{array}
            \right).
        \end{array}
    \f]
    The analytical solutions are
    \f[
        p(x,y,z) = x^2y^2 + 6x + 5z,
        \quad
        \sigma(x,y,z) =
        \left(
        \begin{array}{c}
            - 4xy^2 - 12 - 2x^2y + 2x^3 + 2y \\
            -2xy^2 - 6 - 2x^2y + 2y + x^3 \\
            - 5 + 4z
        \end{array}
        \right).
    \f]
    The computed errors are
    <table border="1" align="center">
        <tr>
            <th> N </th>
            <th> \f$ \left\Vert p - p_h \right\Vert_{L^2} \f$ </th>
            <th> \f$ \left\Vert \sigma - \Pi_{P_0} \sigma_h \right\Vert_{L^2} \f$ </th>
            <th> N proc </th>
        </tr>
        <tr>
            <td> 5 </td>
            <td> 0.338269 </td>
            <td> 0.411764 </td>
            <td> 3 </td>
        </tr>
        <tr>
            <td> 10 </td>
            <td> 0.166443 </td>
            <td> 0.204983 </td>
            <td> 3 </td>
        </tr>
        <tr>
            <td> 20 </td>
            <td> 0.0832153 </td>
            <td> 0.102579 </td>
            <td> 3 </td>
        </tr>
        <tr>
            <td> 40 </td>
            <td> 0.0416069 </td>
            <td> 0.051361 </td>
            <td> 32 </td>
        </tr>
        <tr>
            <td> 80 </td>
            <td> 0.0208047 </td>
            <td> 0.0258553 </td>
            <td> 96 </td>
        </tr>
    </table>
    where N is the number of subdivisions for each boundary.
    @image html darcy/3d.png "Example of the solution with N = 5."
*/
class darcy_linear
    //     :
    //     public LifeV::Application
{
public:

    /* @name Constructors and destructor
     */
    //@{

    //! Constructor
    darcy_linear ( int argc, char** argv );


    //! Destructor
    ~darcy_linear () {}

    //@}

    /* @name  Methods
     */
    //@{

    //! To lunch the simulation
    LifeV::Real run();

    //@}


private:
    struct Private;
    std::shared_ptr<Private> Members;
};

#endif /* __darcy_H */
