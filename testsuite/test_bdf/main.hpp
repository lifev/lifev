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
#include <life/lifecore/life.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/elemOper_ext.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifemesh/markers.hpp>

#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>

#include <life/lifealg/dataAztec.hpp>

namespace LifeV
{
//===================================
//
// DIFFERENTIAL OPERATORS WE WANT TO SOLVE FOR, WRAPPED IN A SUITABLE CLASS
//
typedef EOExpr<Real,Stiff>  EOStiff;
typedef EOExpr<Real,Mass>  EOMass;


// f
class SourceFct
{
public:
  inline Real operator()(Real x,Real y,Real z,Real t,int ic=0) const {
    return -6.;
  }
};

class AnalyticalSol
{
  // ic stands for a component index (unuseful in the scalar case)
public:
  inline Real operator()(Real t, Real x,Real y,Real z, UInt ic) const {
    return t*t*(x*x+y*y+z*z);
  }
  inline Real der_t(Real t, Real x,Real y,Real z, UInt ic) const {
    return 2*t*(x*x+y*y+z*z);
  }
  inline Real der_x(Real t, Real x,Real y,Real z, UInt ic) const {
    return 2*x*t*t;
  }
  inline Real der_y(Real t, Real x,Real y,Real z, UInt ic) const {
    return 2*y*t*t;
  }
  inline Real der_z(Real t, Real x,Real y,Real z, UInt ic) const {
    return 2*z*t*t;
  }
};


/////////////////
// AZTEC's STUFF
////////////////
void init_options(int options[], double params[])
{

  /*
   * Choose among AZTEC options (see User's Guide).
   */

 AZ_defaults(options, params);

 options[AZ_solver]   = AZ_gmres;
 options[AZ_scaling]  = AZ_none;
 options[AZ_conv]     = AZ_r0;
 options[AZ_output]   = 1; //AZ_warnings;
 options[AZ_pre_calc] = AZ_calc;
 options[AZ_max_iter] = 1550;
 options[AZ_poly_ord] = 5;
 options[AZ_overlap]  = AZ_none;
 options[AZ_kspace]   = 40;
 options[AZ_aux_vec]  = AZ_resid;
 options[AZ_precond] =  AZ_dom_decomp;
 options[AZ_subdomain_solve] = AZ_ilut;
 options[AZ_keep_info] = 0;
 params[AZ_tol]       = 1.00e-10;
 params[AZ_drop]      = 1.00e-4;
 params[AZ_ilut_fill] = 5;
 params[AZ_omega]     = 1.;

} /* init_options */
//////////////////////////////////////////////////

}
