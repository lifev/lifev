/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  This file is part of the LifeV library.

  Author: Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  Copyright (C) 2004 EPFL

  Distributed under the GPL(GNU Public License):
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/
/** \file MatrixTest.hpp
    
*/
#ifndef __MatrixTest_H
#define __MatrixTest_H 1

#include <values.hpp>

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


#endif /* __MatrixTest_H */
