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

#include "operFS.hpp"

namespace LifeV
{

// Constructors

operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
               VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
               BCHandler &BCh_du,
               BCHandler &BCh_dz,
               GetPot     &data_file):
    M_fluid       (fluid),
    M_solid       (solid),
    M_dispStruct  ( 3*M_solid.dDof().numTotalDof() ),
    M_velo        ( 3*M_solid.dDof().numTotalDof() ),
    M_dz          ( 3*M_solid.dDof().numTotalDof() ),
    M_rhs_dz      ( 3*M_solid.dDof().numTotalDof() ),
    M_residualS   ( M_solid.dDof().numTotalDof() ),
    M_residualF   ( M_fluid.uDof().numTotalDof() ),
    M_residualFSI ( M_fluid.uDof().numTotalDof() ),
    M_nbEval      (0),
    M_BCh_du      (BCh_du),
    M_BCh_dz      (BCh_dz),
    M_dataJacobian(this)
{
    M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
}

// Destructor
operFS::~operFS()
{
}

// Setters and getters




// Member functions



//! Computing the residual on the fluid/structure interface

void operFS::computeResidualFSI()
{
    M_residualFSI = 0.;

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                M_residualF[IDfluid - 1 + jDim*totalDofFluid] +
                M_residualS[IDsolid - 1 + jDim*totalDofSolid];
        }
    }
}


void operFS::setResidualFSI(double *_res)
{
    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                _res[IDsolid - 1 + jDim*totalDofSolid];
        }
    }
}

void operFS::setResidualFSI(const Vector _res)
{
    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                _res[IDsolid - 1 + jDim*totalDofSolid];
        }
    }
}


Vector operFS::getResidualFSIOnSolid()
{
    Vector vec = M_residualS;
    vec = 0.;

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            vec[IDsolid - 1 + jDim*totalDofSolid] =
                M_residualF[IDfluid - 1 + jDim*totalDofFluid] +
                M_residualS[IDsolid - 1 + jDim*totalDofSolid];
        }
    }

    return vec;
}

Vector operFS::getSolidInterfaceOnFluid(Vector &_vec)
{
    Vector vec(M_residualF.size());

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            vec[IDfluid - 1 + jDim*totalDofSolid] =
                _vec[IDsolid - 1 + jDim*totalDofSolid];
}
    }

    return vec;
}


Vector operFS::getFluidInterfaceOnSolid(Vector &_vec)
{
    Vector vec(M_residualS.size());

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

        for (UInt jDim = 0; jDim < nDimF; ++jDim)
        {
            vec[IDsolid - 1 + jDim*totalDofSolid] =
                _vec[IDfluid - 1 + jDim*totalDofSolid];
}
    }

    return vec;
}

//
// Residual evaluation
//


void operFS::eval(const Vector& disp,
                  int           status,
                  Vector&       dispNew,
                  Vector&       velo)
{
    if(status) M_nbEval = 0; // new time step
    M_nbEval++;

    UInt nDofInterface;
    nDofInterface = M_fluid.BC_fluid()[1].list_size();

    M_solid.d() = disp;

    M_solid._recur = 0;
    M_solid.iterate();

    M_fluid.updateMesh(M_time);
    M_fluid.iterate   (M_time);

    dispNew = M_solid.d();
    velo    = M_solid.w();

    std::cout << "                ::: norm(disp     ) = "
              << maxnorm(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = "
              << maxnorm(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = "
              << maxnorm(velo) << std::endl;
}

Vector operFS::evalResidual(Vector &disp,
                            int iter)
{
    Vector res(disp.size());

    int status = 0;

    if(iter == 0) status = 1;

    std::cout << "*** Residual computation g(x_" << iter <<" )";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    eval(disp, status, M_dispStruct, M_velo);

    M_residualS = M_solid.residual();
    M_residualF = M_fluid.residual();

    computeResidualFSI();

    res = getResidualFSIOnSolid();

    std::cout << "Max ResidualF   = " << maxnorm(M_residualF)
              << std::endl;
    std::cout << "Max ResidualS   = " << maxnorm(M_residualS)
              << std::endl;
    std::cout << "Max ResidualFSI = " << maxnorm(M_residualFSI)
              << std::endl;

    return res;
}

//
// Preconditionner computation using AZTEC
//

// void  operFS::updateJac(Vector& sol,int iter)
// {
// }


//

void  operFS::solvePrec(Vector &_uBarS)
{
    // data_org assigned "by hands": no parallel computation is performed

    UInt dim = _uBarS.size();


    for (int ii = 0; ii < (int) dim; ++ii)
        solid().d()[ii] =  _uBarS[ii];

    fluid().updateDispVelo();
    solveLinearFluid();
    solveLinearSolid();

    for (int i = 0; i < (int) dim; ++i)
        _uBarS[i] =  dz()[i] - _uBarS[i];

}

//
void  operFS::updateJac(Vector& sol,int iter) {
}



Vector  operFS::solvePrec(const Vector  &_res,
                          double        _linearRelTol)
{
    Vector muk(_res.size());

    M_linearRelTol = _linearRelTol;

    UInt precChoice = 0;

    std::cout << "  o-  Solving the preconditionned system... ";

    Chrono chrono;

    switch(precChoice)
    {
        case 0:
            // Dirichlet-Neumann preconditioner
            muk = invSfPrime(_res);
            break;
        case 1:
            // Dirichlet-Neumann preconditioner
            muk = invSsPrime(_res);
            break;
        case 2:
            // Dirichlet-Neumann preconditioner
        {
            Vector muF(_res.size());
            Vector muS(_res.size());

            muF = invSfPrime(_res);
            muS = invSsPrime(_res);

            muk = muS + muF;
        }
        break;
        default:
            // Newton preconditioner
            muk = invSfSsPrime(_res);
    }

    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    return muk;
}

//

void  operFS::solveLinearFluid()
{
    M_fluid.iterateLin(M_time, M_BCh_du);
}

//

void  operFS::solveLinearSolid()
{

    M_rhs_dz = 0.;
    M_dz     = 0.;

    if ( !M_BCh_dz.bdUpdateDone() )
        M_BCh_dz.bdUpdate(M_solid._mesh, M_solid._feBd,
                          M_solid._dof);

    bcManageVector(M_rhs_dz, M_solid._mesh, M_solid._dof,
                     M_BCh_dz, M_solid._feBd, 1., 1.);

    Real tol       = 1.e-10;

    M_BCh_dz.showMe(true, std::cout);
    std::cout << "rhs_dz norm = " << maxnorm(M_rhs_dz) << std::endl;
    M_solid._recur = 1;
    M_solid.solveJac(M_dz, M_rhs_dz, tol, M_BCh_dz);
    std::cout << "dz norm     = " << maxnorm(M_dz) << std::endl;
}



//
//
//


Vector operFS::invSfPrime(const Vector& _res)
{
    Vector step(_res.size());

    setResidualFSI(_res);
    solveLinearFluid();

    Vector deltaLambda = M_fluid.getDeltaLambda();

    step = getFluidInterfaceOnSolid(deltaLambda);
    return step;
}


//
//
//



Vector operFS::invSsPrime(const Vector& _res)
{
    Vector step(_res.size());

    setResidualFSI(_res);
    solveLinearSolid();

    for (int i = 0; i < (int) M_dz.size(); ++i)
        step[i] = dz()[i];

    return step;
}


//
//
//


Vector operFS::invSfSsPrime(const Vector& _res)
{
    UInt dimRes = _res.size();
    Vector step(dimRes);

    M_solverAztec.setTolerance(M_linearRelTol);

    M_solverAztec.setMatrixFree(dimRes,
                             &M_dataJacobian,
                             my_matvecSfSsPrime);

    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dimRes; ++i)
        step[i]=0.0;

    chrono.start();
    M_solverAztec.solve(step, _res);
    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    return step;
}


void my_matvecSfSsPrime(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
{
    // Extraction of data from J
    DataJacobian* my_data = static_cast< DataJacobian* >(AZ_get_matvec_data(J));

    UInt dim = my_data->M_pFS->dz().size();

    double xnorm =  AZ_gvector_norm(dim, -1, z, proc_config);
    std::cout << " ***** norm (z)= " << xnorm << std::endl<< std::endl;

    if ( xnorm == 0.0 )
        for (int i=0; i <(int)dim; ++i)
            Jz[i] =  0.0;
    else
    {
        my_data->M_pFS->setResidualFSI(z);

        // Here we have to replace steps
        //               S'_s^{-1}\circ S'_f  z
        // by
        //               S'_s \cdot z  +  S'_f \cdot z

        my_data->M_pFS->fluid().updateDispVelo();
        my_data->M_pFS->solveLinearFluid();
        my_data->M_pFS->solveLinearSolid();

        // Here we have to get rid of z[i]
        // Jz = to the computations just above.

        for (int i = 0; i < (int) dim; ++i)
            Jz[i] =  z[i] - my_data->M_pFS->dz()[i];
    }
    std::cout << " ***** norm (Jz)= "
              << AZ_gvector_norm(dim, -1, Jz, proc_config)
              << std::endl << std::endl;
}
}
