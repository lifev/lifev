#include "NavierStokesAleSolverPC.hpp"
#include "VenantKirchhofSolver.hpp"
#include "norm.hpp"

#ifndef _OPERFS
#define _OPERFS


class operFS {

 public:

  operFS(NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& fluid, 
	 VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& solid);

  void eval(Vector& dispNew, Vector& veloStruct, const Vector& disp,int status);

  NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& _fluid;
  
  VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& _solid;
  
};


operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D<LinearTetra> >& fluid, 
	       VenantKirchhofSolver< RegionMesh3D<LinearTetra> >& solid):
     _fluid(fluid),
     _solid(solid) {}


void operFS::eval(Vector& dispNew, Vector& velo, const Vector& disp, int status) {

  _solid.d().vec() = disp;

  _fluid.updateMesh(); 
  _fluid.iterate();  
  _solid.iterate(); 
 
  dispNew = _solid.d().vec(); 
  velo    = _solid.w().vec(); 

  cout << "                ::: norm(disp     ) = " << maxnorm(disp)    << endl;
  cout << "                ::: norm(dispNew  ) = " << maxnorm(dispNew) << endl;
  cout << "                ::: norm(velo     ) = " << maxnorm(velo)    << endl;
}




#endif
