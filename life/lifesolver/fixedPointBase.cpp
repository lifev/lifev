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
fixedPoint::fixedPoint( fluid_type& fluid,
                        solid_type& solid,
                        GetPot &_dataFile,
                        BCHandler &BCh_u,
                        BCHandler &BCh_d,
                        BCHandler &BCh_mesh):
    operFS(fluid, solid, _dataFile, BCh_u, BCh_d, BCh_mesh)
{
    M_defOmega =  _dataFile("problem/defOmega",0.01);
    std::cout << "Default aikten start value = " << M_defOmega
              << std::endl;
    setUpBC();
}

fixedPoint::~fixedPoint()
{}

void
fixedPoint::setup()
{
    // call operFS setup()
    super::setup();

    setUpBC();
}

void fixedPoint::eval(Vector& dispNew, Vector& velo, const Vector& disp, int status)
{
    if(status) M_nbEval = 0; // new time step
    M_nbEval++ ;

    M_solid->d() = disp;

    M_fluid->updateMesh(M_time);
    M_fluid->iterate(M_time);

    M_solid->setRecur(0);
    M_solid->iterate();

    dispNew = M_solid->d();
    velo    = M_solid->w();

    std::cout << "                ::: norm(disp     ) = " << norm_inf(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = " << norm_inf(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = " << norm_inf(velo) << std::endl;
}


// Residual evaluation
//
void fixedPoint::evalResidual(Vector &res, const Vector& disp, int iter)
{
    int status = 0;
    if(iter == 0) status = 1;
    std::cout << "*** Residual computation g(x_" << iter <<" )";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;
    eval(M_dispStruct, M_velo, disp, status);
    res = M_dispStruct - disp;
}



//
// Boundary conditions setup
//


void fixedPoint::setUpBC()
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid = this->M_solid->dDof().numTotalDof();
    UInt dim_fluid = this->M_fluid->uDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    BCVectorInterface g_wall(this->M_fluid->residual(),
                             dim_fluid,
                             M_dofFluidToStructure);
    //
    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    BCVectorInterface displ(this->M_solid->d(),
                            dim_solid,
                            M_dofStructureToFluidMesh);
    //
    // Passing data from structure to the fluid: solid velocity at the interface velocity
    //
    BCVectorInterface u_wall(this->M_fluid->wInterpolated(),
                             dim_fluid,
                             M_dofMeshToFluid);
    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    M_BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);

    // Boundary conditions for the solid displacement
    M_BCh_d.addBC("Interface", 1, Natural, Full, g_wall, 3);
}


//
// new step computation resolution
//


void  fixedPoint::solveJac(Vector        &_muk,
                           const Vector  &_res,
                           const double   _linearRelTol)
{
    _muk = _res;
}

//
// add fixedPoint to factory
//
namespace
{
operFS* createFP(){ return new fixedPoint(); }
static bool reg = FSIFactory::instance().registerProduct( "fixedPoint", &createFP );
}

}
