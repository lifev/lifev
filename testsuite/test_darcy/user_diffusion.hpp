/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef _USERDIFFUSION_H
#define _USERDIFFUSION_H
#include "tab.hpp"

namespace LifeV
{
KNM<double> fibrous_permea(double visc,KNM<double> perm_gen,double x);
double gamma(double x);
double der_gamma(double x);

double permeability_sd009(const double& x, const double& y, const double& z);
double permeability_sd010(const double& x, const double& y, const double& z);
}
#endif


