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
#define N_PHY_UNK 1
#define DEBUG

#define CSR_ORD // ensuring Ordered CSR

#include <vector>
#include <algorithm>
#include "lifeV.hpp"
#include "readMesh3D.hpp"
#include "chrono.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "refFE.hpp"
#include "elemOper_ext.hpp"
#include "dof.hpp"
#include "markers.hpp"

#include "values.hpp"
#include "assemb.hpp"

namespace LifeV
{
//===================================
//
// DIFFERENTIAL OPERATORS WE WANT TO SOLVE FOR, WRAPPED IN A SUITABLE CLASS
//
typedef EOExpr<Real,Stiff>  EOStiff;


/* =========== Problem Dependent */
#define LENGTH 1.
#define OMEGA (M_PI/LENGTH)


// f
class SourceFct
{
public:
  inline Real operator()(Real x,Real y,Real z,int ic=0) const {
    return -6;
  }
};

class AnalyticSol
{
public:
  inline Real operator()(Real x,Real y,Real z,int ic=0) const {
    return x*x + y*y + z*z;
  }
};

}
