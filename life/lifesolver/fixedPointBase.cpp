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
fixedPoint::fixedPoint():
    super(),
    M_defOmega( 0.001 ),
    M_aitkFS()
{
    M_aitkFS.setDefault( M_defOmega );
}


fixedPoint::~fixedPoint()
{}


void
fixedPoint::setDataFromGetPot( GetPot const& data )
{
    // call the super class to setup the data from getpot file if needed
    super::setDataFromGetPot( data );

    M_defOmega = data("problem/defOmega", 0.001);

    Debug( 6205 ) << "fixedPoint::setDataFromGetPot(GetPot) OmegaS = " << M_defOmega << "\n";

    M_aitkFS.setDefault(M_defOmega, 0.001);
}


void
fixedPoint::setup()
{
    // call FSIOperator setup()
    super::setup();
    M_aitkFS.setup( 3*M_solid->dDof().numTotalDof() );
//    setUpBC();
}

void fixedPoint::eval(Vector& dispNew,
                      Vector& velo,
                      const Vector& disp,
                      int status)
{
    if(status) M_nbEval = 0; // new time step
    M_nbEval++ ;

    this->M_solid->d() = disp;

    this->M_fluid->updateMesh(M_time);
    this->M_fluid->iterate(M_time);

    this->M_solid->setRecur(0);
    this->M_solid->iterate();

    dispNew = M_solid->d();
    velo    = M_solid->w();

    std::cout << " ::: norm(disp     ) = " << norm_2(disp) << std::endl;
    std::cout << " ::: norm(dispNew  ) = " << norm_2(dispNew) << std::endl;
    std::cout << " ::: norm(velo     ) = " << norm_2(velo) << std::endl;

    std::cout << "Max ResidualF        = " << norm_inf(M_fluid->residual())
              << std::endl;
    std::cout << "Max ResidualS        = " << norm_inf(M_solid->residual())
              << std::endl;
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

    M_dispStructOld = disp;

    eval(M_dispStruct, M_velo, disp, status);

    res = M_dispStruct - disp;

//     transferOnInterface(res,
//                         M_solid->BC_solid(),
//                         "Interface",
//                         res);
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
                             M_dofFluidToStructure  );
    //
    // Passing data from structure to the fluid mesh: motion of the fluid domain
    //
    BCVectorInterface displ(this->M_solid->d(),
                            dim_solid,
                            M_dofStructureToFluidMesh );
    //========================================================================================
    //  BOUNDARY CONDITIONS
    //========================================================================================

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    M_BCh_mesh->addBC("Interface", 1, Essential, Full, displ, 3);

    // Boundary conditions for the solid displacement
    M_BCh_d->addBC("Interface", 1, Natural, Full, g_wall, 3);
}


//
// new step computation resolution
//


void  fixedPoint::solveJac(Vector        &_muk,
                           const Vector  &_res,
                           const double   /*_linearRelTol*/)
{
    if (M_nbEval == 1) M_aitkFS.restart();
    _muk = M_aitkFS.computeDeltaLambda(M_dispStructOld, -1.*_res);
}


//
// add fixedPoint to factory
//


namespace
{
FSIOperator* createFP(){ return new fixedPoint(); }
static bool reg = FSIFactory::instance().registerProduct( "fixedPoint", &createFP );
}

}
