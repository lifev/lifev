/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-08-29

  Copyright (C) 2004 EPFL

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
   \file MatrixTest.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-08-29
*/
#ifndef __MatrixTest_H
#define __MatrixTest_H 1

#include <sparseArray.hpp>

namespace LifeV
{
class MatrixMass
{
public:
    typedef double value_type;
    typedef CSRMatr<CSRPatt,value_type> matrix_type;

    MatrixMass( int n );

    ~MatrixMass()
        {
            delete _M_mat;
            delete _M_pattern;
        }
    uint const * iaData()const { return _M_mat->Patt()->giveRawCSR_ia(); }
    uint const * jaData()const { return _M_mat->Patt()->giveRawCSR_ja(); }
    double* valueData() { return _M_mat->giveRawCSR_value(); }

    matrix_type const& matrix() const { return *_M_mat; }
    matrix_type &      matrix()       { return *_M_mat; }

private:
    matrix_type* _M_mat;
    CSRPatt* _M_pattern;
    std::vector<double> _M_val;
};

/*!
      \class MatrixConvectionDiffusion


    | T -I          |
    |-I  T -I       |
A = |   -I  T       |
    |        ...  -I|
    |           -I T|

    derived from the standard central difference discretization of the
     2-dimensional convection-diffusion operator (Laplacian u) + rho*(du/dx)
    on a unit square with zero Dirichlet boundary conditions.
    When rho*h/2 <= 1, the discrete convection-diffusion operator has real
    eigenvalues.  When rho*h/2 > 1, it has COMPLEX eigenvalues.

    */
class MatrixConvectionDiffusion
{
public:
    typedef double value_type;
    typedef CSRMatr<CSRPatt,value_type> matrix_type;

    MatrixConvectionDiffusion( int nx, value_type __rho = 0.0 );

    ~MatrixConvectionDiffusion()
        {
            delete _M_mat;
            delete _M_pattern;
        }
    uint const * iaData()const  { return _M_mat->Patt()->giveRawCSR_ia(); }
    uint const * jaData()const  { return _M_mat->Patt()->giveRawCSR_ja(); }
    double* valueData() { return _M_mat->giveRawCSR_value(); }

    matrix_type const& matrix() const { return *_M_mat; }
    matrix_type &      matrix()       { return *_M_mat; }

private:
    value_type _M_rho;
    matrix_type* _M_mat;
    CSRPatt* _M_pattern;
    std::vector<double> _M_val;
};

}
#endif /* __MatrixTest_H */
