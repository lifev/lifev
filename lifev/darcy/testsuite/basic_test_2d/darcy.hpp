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
    @file
    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author Anna Scotti <anna.scotti@mail.polimi.it>

    @date 2012-03-30
*/


#ifndef __darcy_H
#define __darcy_H 1


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/darcy/solver/DarcySolver.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>

/*!
    @class darcy
    @brief LifeV Darcy test case

    @author Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
    @author Anna Scotti <anna.scotti@mail.polimi.it>

    Simple 2D Darcy test with Dirichlet, Neumann and Robin Boundary condition.
    <br>
    Solve the problem in dual-mixed form
    \f[
        \left\{
        \begin{array}{l l l }
            \Lambda^{-1} \sigma + \nabla p = f_v & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
            \nabla \cdot \sigma - f = 0          & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
            p = g_D                              & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
            \sigma \cdot n + h p = g_R           & \mathrm{on} & \Gamma_R\,, \vspace{0.2cm} \\
            \sigma \cdot n = g_{n,i}             & \mathrm{on} & \Gamma_{N,i} \,.
        \end{array}
        \right.
    \f]
    where \f$ \Omega = (0,1)^2 \f$ with \f$ \Gamma_R = \left\{ y = 0 \right\} \f$,
    \f$ \Gamma_{N,1} = \left\{ x = 1 \right\} \f$, \f$ \Gamma_{N,2} = \left\{ y = 1 \right\} \f$ and
    \f$ \Gamma_D = \left\{ x = 0 \right\} \f$.
    Furthermore the data are
    \f[
        \begin{array}{l l l}
            f(x,y) = -1 + 4y,
            &f_v(x,y) = \left( x-y,\, y^2 \right)^T,
            &g_D(x,y) = x^2 + xy - y^2, \\
            h(x,y) = 1,
            &g_R(x,y) = x^2 + xy - 3 y^2 + 3x - 2y,
            &g_{N,1}(x,y) = -3 x - 2 y + y^2, \\
            g_{N,2}(x,y) = -3 x + 2 y + 2 y^2,
            &K(x,y) =
            \left[
                \begin{array}{c c}
                    2 & 1 \\
                    1 & 2
                \end{array}
            \right]
        \end{array}
    \f]
    The analytical solutions are
    \f[
        p(x,y) = x^2 + xy - y^2,
        \quad
        \sigma(x,y) =
        \left(
            -3 x - 2 y + y^2,\, -3 x + 2 y + 2 y^2
        \right)^T.
    \f]
    The computed errors with three processors are
    <table border="1" align="center">
        <tr>
            <th> \f$ N \f$ </th>
            <th> \f$ \Vert p - p_h \Vert_{L^2} \f$ </th>
            <th> \f$ \Vert \sigma - \sigma_h \Vert_{L^2} \f$ </th>
        </tr>
        <tr>
            <td> 15 </td>
            <td> 0.0252484 </td>
            <td> 0.138088 </td>
        </tr>
        <tr>
            <td> 30 </td>
            <td> 0.012627 </td>
            <td> 0.0693529 </td>
        </tr>
        <tr>
            <td> 60 </td>
            <td> 0.00631386 </td>
            <td> 0.034748 </td>
        </tr>
        <tr>
            <td> 120 </td>
            <td> 0.00315697 </td>
            <td> 0.0174405 </td>
        </tr>
        <tr>
            <td> 240 </td>
            <td> 0.00157849 </td>
            <td> 0.00871514 </td>
        </tr>
    </table>
    where \f$ N \f$ is the number of subdivisions for each boundary.
*/
class darcy
//     :
//     public LifeV::Application
{
public:

    /*! @name Constructors and destructor
     */
    //@{

    darcy( int argc,
           char** argv );

    ~darcy()
    {}

    //@}

    /*! @name  Methods
     */
    //@{

    LifeV::Real run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __darcy_H */
