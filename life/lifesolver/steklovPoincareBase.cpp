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


#include "steklovPoincareBase.hpp"


namespace LifeV
{
Real fzeroSP(const Real& t,
           const Real& x,
           const Real& y,
           const Real& z,
           const ID& i)
{return 0.0;}

steklovPoincare::steklovPoincare(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
                                 VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
                                 GetPot &_dataFile,
                                 BCHandler &BCh_u,
                                 BCHandler &BCh_d,
                                 BCHandler &BCh_mesh):
    operFS(fluid, solid, _dataFile, BCh_u, BCh_d, BCh_mesh),
    M_BCh_du      (),
    M_BCh_dz      (),
    M_dz          ( 3*M_solid.dDof().numTotalDof() ),
    M_rhs_dz      ( 3*M_solid.dDof().numTotalDof() ),
    M_residualS   ( M_solid.dDof().numTotalDof() ),
    M_residualF   ( M_fluid.uDof().numTotalDof() ),
    M_residualFSI ( M_fluid.uDof().numTotalDof() ),
    M_aitkFS      ( 3*M_solid.dDof().numTotalDof() ),
    M_dataJacobian(this)
{
    M_precond   = _dataFile("problem/precond"  , 1);

    M_defOmegaS = _dataFile("problem/defOmegaS",0.005);
    M_defOmegaF = _dataFile("problem/defOmegaF",0.005);
    M_aitkFS.setDefault(M_defOmegaS, M_defOmegaF);

    // for the NN case, aitken is inside the step
    // computation routine
//     if (M_precond == 2)
//     {
//         M_defOmega  = -1.;
//         M_defOmegaS = _dataFile("problem/defOmegaS",0.005);
//         M_defOmegaF = _dataFile("problem/defOmegaF",0.005);
//         M_aitkFS.setDefault(M_defOmegaS, M_defOmegaF);
//     }
//     else
//     {
//         M_defOmega =  _dataFile("problem/defOmega",0.01);
//         std::cout << "Default aikten start value = " << M_defOmega
//                   << std::endl;
//     }
    setUpBC();
}

steklovPoincare::~steklovPoincare()
{}


//
// Residual computation
//

void steklovPoincare::eval(const Vector& disp,
                           int           status,
                           Vector&       dispNew,
                           Vector&       velo)
{
    if(status) M_nbEval = 0; // new time step
    M_nbEval++;

//    M_solid.d() = setDispOnInterface(disp);

//     transferOnInterface(disp,
//                         M_solid.BC_solid(),
//                         "Interface",
//                         M_solid.d());

    M_solid.d() = disp;

    M_fluid.updateMesh(M_time);
    M_fluid.iterate   (M_time);

    M_solid.setRecur(0);
    M_solid.iterate();
//    M_solid.solveLin(dispNew, disp, 1.e-8);
//    M_solid.postProcess();
//    M_fluid.postProcess();

    dispNew = M_solid.d();
    //M_solid.d() = dispNew;
    velo    = M_solid.w();

    std::cout << "                ::: norm(disp     ) = "
              << norm_2(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = "
              << norm_2(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = "
              << norm_2(velo) << std::endl;
}

void steklovPoincare::evalResidual(Vector &res,
                                   const Vector &disp,
                                   const int iter)
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

    std::cout << "max ResidualF   = " << norm_inf(M_residualF)
              << std::endl;
    std::cout << "max ResidualS   = " << norm_inf(M_residualS)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_residualFSI)
              << std::endl;

//     Vector muk = disp;
//     muk = ZeroVector( muk.size() );
//     invSsPrime(M_residualS, 1e-08, muk);
//     std::cout << "Norm_max d_disp = " << norm_inf(disp - muk) << std::endl;
//     muk = ZeroVector( muk.size() );
//     invSfPrime(M_residualF, 1e-08, muk);
//     std::cout << "Norm_max f_disp = " << norm_inf(disp - muk) << std::endl;
}

//
// Boundary conditions setup
//


// void steklovPoincare::setUpBC(function_type _bcf,
//                               function_type _vel)
void steklovPoincare::setUpBC()
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid = this->M_solid.dDof().numTotalDof();
    UInt dim_fluid = this->M_fluid.uDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================

    BCFunctionBase bcf(fzeroSP);

    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    BCVectorInterface g_wall(this->M_fluid.residual(),
                             dim_fluid,
                             M_dofFluidToStructure);
    //
    // Passing data from structure to the solid mesh: motion of the solid domain
    //
    BCVectorInterface d_wall(this->M_solid.d(),
                             dim_solid,
                             M_dofStructureToSolid);
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

    M_BCh_mesh.addBC("Interface", 1, Essential, Full, displ, 3);

    // Boundary conditions for the solid displacement

    M_BCh_d.addBC("Interface", 1, Essential, Full, d_wall, 3);

    //========================================================================================
    //  COUPLED FSI LINEARIZED OPERATORS
    //========================================================================================

    // Passing the residue to the linearized fluid: \sigma -> du
    //
    // rem: for now: no fluid.dwInterpolated().
    //      In the future this could be relevant

    BCVectorInterface du_wall(M_residualFSI,
                              dim_fluid,
                              M_dofMeshToFluid);
    // Passing the residual to the linearized structure: \sigma -> dz
    BCVectorInterface dg_wall(M_residualFSI,
                              dim_fluid,
                              M_dofFluidToStructure);

    // Boundary conditions for du

    M_BCh_du.addBC("Wall",   1,  Natural  , Full, du_wall,  3);
    M_BCh_du.addBC("Edges",  20, Essential, Full, bcf,      3);

    // Boundary conditions for dz

    M_BCh_dz.addBC("Interface", 1, Natural  , Full, dg_wall, 3);
    M_BCh_dz.addBC("Top",       3, Essential, Full, bcf,     3);
    M_BCh_dz.addBC("Base",      2, Essential, Full, bcf,     3);
    std::cout << "ok." << std::endl;
}


//
// new step computation resolution
//


void  steklovPoincare::solveJac(Vector &muk,
                                const Vector  &_res,
                                double        _linearRelTol)
{
    Vector muF(_res.size());
    Vector muS(_res.size());

    switch(M_precond)
    {
        case 0:
            // Neumann-Dirichlet preconditioner
            invSfPrime(-1*_res, _linearRelTol, muF);
            break;
        case 1:
            // Dirichlet-Neumann preconditioner
            invSsPrime(-1*_res, _linearRelTol, muS);
            break;
        case 2:
            // Neumann-Neumann preconditioner
        {
            invSsPrime(-1*_res, _linearRelTol, muS);
            invSfPrime(-1*_res, _linearRelTol, muF);

            std::cout << "norm_inf muS = " << norm_inf(muS) << std::endl;
            std::cout << "norm_inf muF = " << norm_inf(muF) << std::endl;

            // indeed, Here we should call Aitken (a different instance than
            // the one defined in nonLinRichardson) and
            // we should replace the defaultOmega for he nonLinRichardson
            // to -1, such that there is no relaxation at all there.
            //
            // Maybe we have to define precChoice as Memeber -> M_precChoice
            // and get the choice as input parameter with the constructor
            // and add a memeber
            // generalizedAitken aitkDN(...) (DN for Dirichlet Neumann)
            // then
            //
//            muk = M_aitkFS.computeDeltaLambda(getResidualFSIOnSolid(), muS, muF );
            // muk = muS + muF;
            // muk = .9*muS + .1*muF;

//            std::cout << "maxnorm muk = " << norm_inf(muk) << std::endl;
        }
        break;
        default:
            // Newton preconditioner
            invSfSsPrime(_res, _linearRelTol, muk);
    }

    if (M_nbEval == 1) M_aitkFS.restart();
    muk = M_aitkFS.computeDeltaLambda(M_dispStruct, muF, muS );
}



void steklovPoincare::solveLinearFluid()
{
    this->M_fluid.iterateLin(time(), M_BCh_du);
}

//

void steklovPoincare::solveLinearSolid()
{
    M_rhs_dz = ZeroVector( M_rhs_dz.size() );
    M_dz     = ZeroVector( M_dz.size() );

    if ( !M_BCh_dz.bdUpdateDone() )
        M_BCh_dz.bdUpdate(this->M_solid.mesh(),
                          this->M_solid.feBd(),
                          this->M_solid.dof());

    bcManageVector(M_rhs_dz,
                   this->M_solid.mesh(),
                   this->M_solid.dof(),
                   M_BCh_dz,
                   this->M_solid.feBd(),
                   1., 1.);

    Real tol       = 1.e-10;

    std::cout << "rhs_dz norm = " << norm_inf(M_rhs_dz) << std::endl;
    this->M_solid.setRecur(1);
    this->M_solid.updateJac(M_dz, 0);
    this->M_solid.solveJac(M_dz, M_rhs_dz, tol, M_BCh_dz);
    std::cout << "dz norm     = " << norm_inf(M_dz) << std::endl;
}



//
//
//


void  steklovPoincare::invSfPrime(const Vector& res,
                                  double linear_rel_tol,
                                  Vector& step)
{
    setResidualFSI(res);

    M_solid.d() = ZeroVector(M_solid.d().size());
    std::cout << "norm_inf residual FSI = " << norm_inf(M_residualFSI);
    this->M_fluid.updateDispVelo();
    solveLinearFluid();

    Vector deltaLambda = this->M_fluid.getDeltaLambda();

    transferOnInterface(deltaLambda,
                        M_fluid.BC_fluid(),
                        "Interface",
                        step);
    std::cout << "norm_2 deltaLambda = " << norm_2(deltaLambda) << std::endl;
    std::cout << "norm_2 step        = " << norm_2(step) << std::endl;
}


//
//
//



void  steklovPoincare::invSsPrime(const Vector& res,
                                  double linear_rel_tol,
                                  Vector& step)
{
    setResidualFSI(res);
    solveLinearSolid();

//    step = setDispOnInterface(M_dz);

//    step = M_dz;

    transferOnInterface(M_dz,
                        M_solid.BC_solid(),
                        "Interface",
                        step);

//     for (int ii = 0; ii < step.size(); ++ii)
//         std::cout << step[ii] << std::endl;

}


void steklovPoincare::invSfSsPrime(const Vector& res,
                                   double linear_rel_tol,
                                   Vector& step)
{
    UInt dimRes = res.size();

    M_solverAztec.setTolerance(M_linearRelTol);

    M_solverAztec.setMatrixFree(dimRes,
                                &M_dataJacobian,
                                my_matvecSfSsPrime);

    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dimRes; ++i)
        step[i]=0.0;

    chrono.start();
    M_solverAztec.solve(step, res);
    chrono.stop();

    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


void my_matvecSfSsPrime(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
{
    // Extraction of data from J
    DataJacobianSP* my_data = static_cast< DataJacobianSP* >(AZ_get_matvec_data(J));

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

void steklovPoincare::computeResidualFSI()
{
    M_residualFSI = ZeroVector( M_residualFSI.size() );

    UInt iBCf = M_fluid.BC_fluid().getBCbyName("Interface");
    UInt iBCd = M_solid.BC_solid().getBCbyName("Interface");

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[iBCf];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[iBCd];

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
                M_residualF[IDfluid - 1 + jDim*totalDofFluid] -
                M_residualS[IDsolid - 1 + jDim*totalDofSolid];
        }
    }
}


void steklovPoincare::setResidualFSI(double *_res)
{
    UInt iBCf = M_fluid.BC_fluid().getBCbyName("Interface");
    UInt iBCd = M_solid.BC_solid().getBCbyName("Interface");

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[iBCf];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[iBCd];

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

void steklovPoincare::setResidualFSI(const Vector _res)
{
    UInt iBCf = M_fluid.BC_fluid().getBCbyName("Interface");
    UInt iBCd = M_solid.BC_solid().getBCbyName("Interface");

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[iBCf];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[iBCd];

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


Vector steklovPoincare::getResidualFSIOnSolid()
{
    Vector vec = M_residualS;
    vec = ZeroVector( vec.size() );

    BCBase const &BC_fluidInterface = M_fluid.BC_fluid()[0];
    BCBase const &BC_solidInterface = M_solid.BC_solid()[0];

    UInt nDofInterface = BC_fluidInterface.list_size();

    UInt nDimF = BC_fluidInterface.numberOfComponents();
    UInt nDimS = BC_solidInterface.numberOfComponents();

    UInt totalDofFluid = M_residualF.size()/ nDimF;
    UInt totalDofSolid = M_residualS.size()/ nDimS;

    std::cout << "total dof fluid = " << totalDofFluid << std::endl;
    std::cout << "total dof solid = " << totalDofSolid << std::endl;

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
                M_residualF[IDfluid - 1 + jDim*totalDofFluid] -
                M_residualS[IDsolid - 1 + jDim*totalDofSolid];
        }
    }

    return vec;
}

Vector steklovPoincare::getSolidInterfaceOnFluid(Vector &_vec)
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


Vector steklovPoincare::getFluidInterfaceOnSolid(Vector &_vec)
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

void steklovPoincare::transferOnInterface(const Vector      &_vec1,
                                          const BCHandler   &_BC,
                                          const std::string &_BCName,
                                          Vector            &_vec2)
{
    int iBC = _BC.getBCbyName(_BCName);

    BCBase const &BCInterface = _BC[(UInt) iBC];

    UInt nDofInterface = BCInterface.list_size();

//     std::cout << "nDofInterface = " << nDofInterface << std::endl;

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

        for (UInt jDim = 0; jDim < nDim; ++jDim)
        {
            _vec2[ID2 - 1 + jDim*totalDof2] =
                _vec1[ID1 - 1 + jDim*totalDof1];
        }
    }
}

}
