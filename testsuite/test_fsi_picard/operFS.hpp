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
#ifndef _OPERFS
#define _OPERFS

#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "regionMesh3D_ALE.hpp"

namespace LifeV
{

using namespace std;

class operFS {

 public:

  operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
	 VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid);

  void setTime(const Real time) {
    _time = time;
  }

  void eval(Vector& dispNew, Vector& veloStruct, const Vector& disp,int status);

  NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& _fluid;

  VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& _solid;

 private:

  Real _time;

};


operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
	       VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid):
     _fluid(fluid),
     _solid(solid) {}


void operFS::eval(Vector& dispNew, Vector& velo, const Vector& disp, int status) {

  _solid.d() = disp;

  _fluid.updateMesh(_time);
  _fluid.iterate(_time);
  _solid.iterate();

  dispNew = _solid.d();
  velo    = _solid.w();

  cout << "                ::: norm(disp     ) = " << maxnorm(disp)    << endl;
  cout << "                ::: norm(dispNew  ) = " << maxnorm(dispNew) << endl;
  cout << "                ::: norm(velo     ) = " << maxnorm(velo)    << endl;
}

}


#endif
