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
               BCHandler &BCh_dz):
//                   GetPot     &data_file):
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
//        M_solverAztec.setOptionsFromGetPot(data_file,"jacobian/aztec");
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

    M_residualFSI = 0.;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID IDfluid = BC_fluidInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BC_fluidInterface.pointerToBCVector());

        ID IDsolid = BCVInterface->
            dofInterface().getInterfaceDof(IDfluid);

//        std::cout << IDsolid << "->" << IDfluid << std::endl;
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


Vector operFS::setDispOnInterface(const Vector &_disp)
{
    Vector vec = _disp.size();
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
                _disp[IDsolid - 1 + jDim*totalDofSolid];
        }
    }

    return vec;
}



Vector operFS::getFluidInterfaceOnSolid(Vector &_vec)
{
    Vector vec(M_residualS.size());
    vec = 0;

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



Vector operFS::getFluidInterfaceOnSolid(Vector    &_vec,
                                        BCHandler &_BCh)
{
    Vector vec(M_residualF.size());
    vec = 0;

//     for (UInt iBC = 1;


    BCBase const &BC_fluidInterface = _BCh[0];
    BCBase const &BC_solidInterface = M_fluid.BC_fluid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;

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
            vec[IDsolid - 1 + jDim*totalDofFluid] =
                _vec[IDfluid - 1 + jDim*totalDofFluid];
        }
    }

    return vec;
}
void operFS::transferOnInterface(const Vector      &_vec1,
                                 const BCHandler   &_BC,
                                 const std::string &_BCName,
                                 Vector            &_vec2)
{
    int iBC = -1;
    for (UInt jBC = 0; jBC < _BC.size(); jBC++)
        if (_BC[jBC].name() == _BCName) iBC = jBC;
    if (iBC == -1)
    {
        std::cout << _BCName << " not found on BC ";
        return;
    }

    std::cout << "iBC " << _BCName << " = " << iBC << std::endl;
    std::cout << std::endl;
    // we know the BC number for each mesh, let's swap

    BCBase const &BCInterface = _BC[(UInt) iBC];

    UInt nDofInterface = BCInterface.list_size();

    std::cout << "nDofInterface = " << nDofInterface << std::endl;

    UInt nDim = BCInterface.numberOfComponents();

    UInt totalDof1 = _vec1.size()/ nDim;
    UInt totalDof2 = _vec2.size()/ nDim;

    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)
    {
        ID ID1 = BCInterface(iBC)->id();

        BCVectorInterface const *BCVInterface =
            static_cast <BCVectorInterface const *>
            (BCInterface.pointerToBCVector());

        ID ID2 = BCVInterface->
            dofInterface().getInterfaceDof(ID1);

//        std::cout << "ID1 = " << ID1 << " ID2 = " << ID2 << std::endl;
        for (UInt jDim = 0; jDim < nDim; ++jDim)
        {
            _vec2[ID2 - 1 + jDim*totalDof2] =
                _vec1[ID1 - 1 + jDim*totalDof1];
        }
    }
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


//    Vector sol(disp.size());

    M_solid.d() = setDispOnInterface(disp);
//    M_solid.d() = disp;


    M_solid._recur = 0;
    M_solid.iterate();

    M_fluid.updateMesh(M_time);
    M_fluid.iterate   (M_time);

    M_solid.postProcess();
    M_fluid.postProcess();

    dispNew = M_solid.d();
    velo    = M_solid.w();

    std::cout << "                ::: norm(disp     ) = "
              << maxnorm(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = "
              << maxnorm(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = "
              << maxnorm(velo) << std::endl;
}

void operFS::evalResidual(Vector &res,
                          Vector &disp,
                          int iter)
{
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
}

//
// Preconditionner computation using AZTEC
//

// void  operFS::updateJac(Vector& sol,int iter)
// {
// }


//

void  operFS::updateJac(Vector& sol,int iter) {
}



void  operFS::solveJac(Vector &muk,
                       const Vector  &_res,
                       double        _linearRelTol)
{
    UInt precChoice = 0;


//     Vector mup = muk;
//     invSfPrime(_res, _linearRelTol, muk);
//     invSsPrime(_res, _linearRelTol, mup);

//     for (int ii = 0; ii < muk.size(); ii++)
//         std::cout << muk[ii] << " " << mup[ii] << std::endl;

    switch(precChoice)
    {
        case 0:
            // Dirichlet-Neumann preconditioner
            invSfPrime(_res, _linearRelTol, muk);
            break;
        case 1:
            // Dirichlet-Neumann preconditioner
            invSsPrime(_res, _linearRelTol, muk);
            break;
        case 2:
            // Dirichlet-Neumann preconditioner
        {
            Vector muF(_res.size());
            Vector muS(_res.size());

            invSfPrime(_res, _linearRelTol, muF);
            invSsPrime(_res, _linearRelTol, muS);

            muk = muS + muF;
        }
        break;
        default:
            // Newton preconditioner
            invSfSsPrime(_res, _linearRelTol, muk);
    }

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

    std::cout << "rhs_dz norm = " << maxnorm(M_rhs_dz) << std::endl;
    M_solid._recur = 1;
    M_solid.updateJac(M_dz, 0);
    M_solid.solveJac(M_dz, M_rhs_dz, tol, M_BCh_dz);
    std::cout << "dz norm     = " << maxnorm(M_dz) << std::endl;
}



//
//
//


void  operFS::invSfPrime(const Vector& res,
                         double linear_rel_tol,
                         Vector& step)
{
    // step = S'_f^{-1} \cdot res

    setResidualFSI(-1*res);
//     transferOnInterface(res,
//                          M_BCh_dz,
//                          "Interface",
//                          M_residualFSI);
    this->M_fluid.updateDispVelo();
//    M_residualFSI = 0.;
    solveLinearFluid();

    Vector deltaLambda = this->M_fluid.getDeltaLambda();
//      for (int ii = 0; ii < deltaLambda.size(); ++ii)
//          std::cout << deltaLambda[ii] << std::endl;
    Vector deltaL = deltaLambda;
    deltaL = 0.;
    transferOnInterface(deltaLambda,
                        M_BCh_du,
                        "Wall",
//                         M_fluid.BC_fluid(),
//                         "Interface",
                        deltaL);

    std::cout << "deltaLambda norm = " << maxnorm(deltaLambda) << std::endl;
    std::cout << "ResidualFSI norm = " << maxnorm(M_residualFSI) << std::endl;

    transferOnInterface(deltaL,
                        M_fluid.BC_fluid(),
                        "Interface",
                        step);

//     for (int ii = 0; ii < M_residualFSI.size(); ++ii)
//         std::cout << M_residualFSI[ii] << " "
//                   << deltaLambda[ii] << " "
//                   << deltaL[ii] << std::endl;
    std::cout << "deltaLambda step = " << maxnorm(step) << std::endl;
}


//
//
//



void  operFS::invSsPrime(const Vector& res,
                         double linear_rel_tol,
                         Vector& step)
{
    setResidualFSI(res);
    solveLinearSolid();

    step = setDispOnInterface(M_dz);

//     for (int ii = 0; ii < step.size(); ++ii)
//         std::cout << step[ii] << std::endl;

//      for (int i = 0; i < (int) M_dz.size(); ++i)
//          step[i] = dz()[i];

//     return;
}


//
//
//


void  operFS::invSfSsPrime(const Vector& res,
                           double linear_rel_tol,
                           Vector& step)
{
    // AZTEC specifications for the second system
    int    data_org[AZ_COMM_SIZE];   // data organisation for J
    int    proc_config[AZ_PROC_SIZE];  // Processor information:
    int    options[AZ_OPTIONS_SIZE];   // Array used to select solver options.
    double params[AZ_PARAMS_SIZE];     // User selected solver paramters.
    double status[AZ_STATUS_SIZE];     // Information returned from AZ_solve()

    AZ_set_proc_config(proc_config, AZ_NOT_MPI);

    // data_org assigned "by hands": no parallel computation is performed
    UInt dim_res = res.size();
    data_org[AZ_N_internal]= dim_res;
    data_org[AZ_N_border]= 0;
    data_org[AZ_N_external]= 0;
    data_org[AZ_N_neigh]= 0;

    // Recovering AZTEC defaults options and params
    AZ_defaults(options,params);

    // Fixed Aztec options for this linear system
    options[AZ_solver]   = AZ_gmres;
    options[AZ_output]   = 1;
    options[AZ_poly_ord] = 5;
    options[AZ_kspace]   = 40;
    options[AZ_conv]     = AZ_rhs;
    params[AZ_tol]       = linear_rel_tol;

    //AZTEC matrix for the jacobian
    AZ_MATRIX *J;
    J = AZ_matrix_create(dim_res);

    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecSfSsPrime);

    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dim_res; ++i)
        step[i]=0.0;

    chrono.start();
    AZ_iterate(&step[0], &res[0], options, params, status, proc_config, J, NULL, NULL);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;


    AZ_matrix_destroy(&J);
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
