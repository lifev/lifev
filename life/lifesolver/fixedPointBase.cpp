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


#include "fixedPointBase.hpp"


namespace LifeV
{
fixedPoint::fixedPoint(GetPot &_dataFile):
    operFS(_dataFile)
{
}

fixedPoint::~fixedPoint()
{}


//
// Residual computation
//

void fixedPoint::eval(const Vector &_disp,
                      const int     _status)
{
    this->M_solid.d() = _disp;

    this->M_fluid.updateMesh(time());
    this->M_fluid.iterate   (time());

    this->M_fluid.postProcess();

    this->M_solid.setRecur(0);
    this->M_solid.iterate();

    this->M_solid.postProcess();
}


void fixedPoint::evalResidual(const Vector &_disp,
                              const int     _iter,
                              Vector       &_res)
{
    int status = 0;

    if(_iter == 0) status = 1;

    std::cout << "*** Residual computation g(x_" << _iter <<")";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    eval(M_dispStruct, status);

    M_dispStruct = this->M_solid.d();
    M_velo       = this->M_solid.w();

    std::cout << "                ::: norm(disp     ) = "
              << maxnorm(_disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = "
              << maxnorm(M_dispStruct) << std::endl;
    std::cout << "                ::: norm(velo     ) = "
              << maxnorm(M_velo) << std::endl;

    _res = M_dispStruct - _disp;
}


//
// Boundary conditions setup
//


void fixedPoint::setUpBC(function_type _bcf,
                         function_type _vel)
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid = this->M_solid.dDof().numTotalDof();
    UInt dim_fluid = this->M_fluid.uDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    BCVectorInterface g_wall(this->M_fluid.residual(),
                             dim_fluid,
                             M_dofFluidToStructure);
    //
    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    BCVectorInterface displ(this->M_solid.d(),
                            dim_solid,
                            M_dofStructureToFluidMesh);
    //
    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    BCVectorInterface u_wall(this->M_fluid.wInterpolated(),
                             dim_fluid,
                             M_dofMeshToFluid);
    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    BCFunctionBase bcf(_bcf);
    M_BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);
    M_BCh_mesh.addBC("Top",       3, Essential, Full, bcf,   3);
    M_BCh_mesh.addBC("Base",      2, Essential, Full, bcf,   3);
    M_BCh_mesh.addBC("Edges",    20, Essential, Full, bcf,   3);


    // Boundary conditions for the fluid velocity
    BCFunctionBase in_flow(_vel);
    M_BCh_u.addBC("Wall",   1,  Essential, Full, u_wall,  3);
    M_BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
    M_BCh_u.addBC("Edges",  20, Essential, Full, bcf,     3);

    // Boundary conditions for the solid displacement
    M_BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
    M_BCh_d.addBC("Top",       3, Essential, Full, bcf,  3);
    M_BCh_d.addBC("Base",      2, Essential, Full, bcf,  3);
}


//
// new step computation resolution
//


void  fixedPoint::solvePrec(const Vector  &_res,
                            const double   _linearRelTol,
                            Vector        &_muk)
{
    _muk = _res;
}


}
