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

steklovPoincare::steklovPoincare():
    super(),
    M_BCh_du( new BCHandler ),
    M_BCh_dz( new BCHandler ),
    M_BCh_du_inv( new BCHandler ),
    M_BCh_dz_inv( new BCHandler ),
    M_dzSolid(),
    M_dzFluid(),
    M_rhs_dz(),
    M_residualS(),
    M_residualF(),
    M_residualFSI(),
    M_defOmega( 0.005 ),
    M_defOmegaS( 0.005 ),
    M_defOmegaF( 0.005 ),
    M_aitkFS(),
    M_dataJacobian( this )
{
    this->setPreconditioner( DIRICHLET_NEUMANN );
    this->setDDNPreconditioner( DDN_DIRICHLET_NEUMANN );
    M_aitkFS.setDefault( M_defOmegaS, M_defOmegaF );
}

steklovPoincare::~steklovPoincare()
{}

void
steklovPoincare::setDataFromGetPot( GetPot const& data )
{
    // call the super class to setup the data from getpot file if needed
    super::setDataFromGetPot( data );

    M_defOmegaS = data("problem/defOmegaS",0.005);
    M_defOmegaF = data("problem/defOmegaF",0.005);

    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) OmegaS = " << M_defOmegaS << "\n";
    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) OmegaF = " << M_defOmegaF << "\n";

    M_aitkFS.setDefault(M_defOmegaS, M_defOmegaF);

    this->setPreconditioner   (  ( OperFSPreconditioner )data("problem/precond"     , DIRICHLET_NEUMANN ) );
    this->setDDNPreconditioner(  ( DDNPreconditioner    )data("problem/DDNprecond"  , DDN_DIRICHLET_NEUMANN ) );

    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) prec   = " << this->preconditioner() << "\n";
}
void
steklovPoincare::setup()
{
    // call operFS setup()
    super::setup();

    M_dzSolid.resize( 3*M_solid->dDof().numTotalDof() );
    M_dzFluid.resize( 3*M_fluid->uDof().numTotalDof() );
    M_rhs_dz.resize( 3*M_solid->dDof().numTotalDof() );
    M_residualS.resize( M_solid->dDof().numTotalDof() );
    M_residualF.resize( M_fluid->uDof().numTotalDof() );
    M_residualFSI.resize( M_fluid->uDof().numTotalDof() );

    M_aitkFS.setup( 3*M_solid->dDof().numTotalDof() );

}
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

    M_solid->d() = disp;

    M_fluid->updateMesh(M_time);
    M_fluid->iterate   (M_time);

    M_solid->setRecur(0);
    M_solid->iterate();

    dispNew = M_solid->d();
    velo    = M_solid->w();

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

    M_dispStructOld = disp;

    eval(disp, status, M_dispStruct, M_velo);

    M_residualS = M_solid->residual();
    M_residualF = M_fluid->residual();

    computeResidualFSI();
    res = getResidualFSIOnSolid();

    std::cout << "max ResidualF   = " << norm_inf(M_residualF)
              << std::endl;
    std::cout << "max ResidualS   = " << norm_inf(M_residualS)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_residualFSI)
              << std::endl;

//      Vector muk = disp;
//      muk = ZeroVector( muk.size() );
//      invSsPrime(M_residualS, 1e-08, muk);
//      std::cout << "Norm_max d_disp = " << norm_inf(disp - muk) << std::endl;
//     muk = ZeroVector( muk.size() );
//     invSfPrime(M_residualF, 1e-08, muk);
//     std::cout << "Norm_max f_disp = " << norm_inf(disp - muk) << std::endl;
}

//
// Boundary conditions setup
//


void steklovPoincare::setUpBC()
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid = this->M_solid->dDof().numTotalDof();
    UInt dim_fluid = this->M_fluid->uDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================

    BCFunctionBase bcf(fzeroSP);

    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    BCVectorInterface g_wall( this->M_fluid->residual(),
                              dim_fluid,
                              M_dofFluidToStructure );

    //
    // Passing data from structure to the solid mesh: motion of the solid domain
    //
    BCVectorInterface d_wall(this->M_solid->d(),
                             dim_solid,
                             M_dofStructureToSolid );
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

    M_BCh_mesh->addBC("Interface", 1, Essential, Full, displ, 3);

    // Boundary conditions for the solid displacement

    M_BCh_d->addBC("Interface", 1, Essential, Full, d_wall, 3);

    //========================================================================================
    //  COUPLED FSI LINEARIZED OPERATORS
    //========================================================================================

    // Passing the residue to the linearized fluid: \sigma -> du
    //
    // rem: for now: no fluid.dwInterpolated().
    //      In the future this could be relevant

    if (this->preconditioner() == NEWTON)
    {
        std::cout << "Steklov-Poincare:NEWTON boundary conditions" << std::endl;

        BCVectorInterface du_wall(M_fluid->dwInterpolated(),
                                  dim_fluid,
                                  M_dofMeshToFluid );

        M_BCh_du->addBC("Wall"     , 1, Essential , Full, du_wall, 3);

        BCVectorInterface du_wall_inv(M_residualFSI,
                                      dim_fluid,
                                      M_dofMeshToFluid);

        M_BCh_du_inv->addBC("Wall", 1, Natural  , Full, du_wall_inv, 3);


        BCVectorInterface dg_wall(M_solid->d(),
                                  dim_solid,
                                  M_dofStructureToSolid );

        M_BCh_dz->addBC("Interface", 1, Essential , Full, dg_wall, 3);

        BCVectorInterface dg_wall_inv(M_residualFSI,
                                      dim_fluid,
                                      M_dofFluidToStructure);

        M_BCh_dz_inv->addBC("Interface", 1, Natural  , Full, dg_wall_inv, 3);
    }
    else
    {
        BCVectorInterface du_wall(M_residualFSI,
                                  dim_fluid,
                                  M_dofMeshToFluid);
        // Passing the residual to the linearized structure: \sigma -> dz
        BCVectorInterface dg_wall(M_residualFSI,
                                  dim_fluid,
                                  M_dofFluidToStructure);
        M_BCh_du->addBC("Wall"     , 1, Natural  , Full, du_wall, 3);
        M_BCh_dz->addBC("Interface", 1, Natural  , Full, dg_wall, 3);

    }

    // Boundary conditions for du and inverse

    M_BCh_du->addBC("Edges",  20, Essential, Full, bcf,      3);
    M_BCh_du->addBC("InFlow", 2,  Natural,   Full, bcf,      3);
    M_BCh_du->addBC("OutFlow",3,  Natural,   Full, bcf,      3);

    M_BCh_du_inv->addBC("Edges",  20, Essential, Full, bcf,  3);
    M_BCh_du_inv->addBC("InFlow", 2,  Natural,   Full, bcf,  3);
    M_BCh_du_inv->addBC("OutFlow",3,  Natural,   Full, bcf,  3);

    // Boundary conditions for dz and inverse

    M_BCh_dz->addBC("Top",       3, Essential, Full, bcf,     3);
    M_BCh_dz->addBC("Base",      2, Essential, Full, bcf,     3);

    M_BCh_dz_inv->addBC("Top",       3, Essential, Full, bcf,     3);
    M_BCh_dz_inv->addBC("Base",      2, Essential, Full, bcf,     3);

    std::cout << "ok." << std::endl;

    std::cout << "BC U\n";
    M_BCh_u->showMe();
    std::cout << "BC dU Inv\n";
    M_BCh_du_inv->showMe();
    std::cout << "BC D\n";
    M_BCh_d->showMe();

    std::cout << "BC mesh\n";
    M_BCh_mesh->showMe();
#if 0
    UInt iBCd = M_solid->BC_solid().getBCbyName("Interface");
    UInt iBCf = M_fluid->BC_fluid().getBCbyName("Interface");
#endif
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
    Debug(  6205 ) << "steklovPoincare::solveJac _linearRelTol  : " << _linearRelTol << "\n";
    Debug(  6205 ) << "steklovPoincare::solveJac preconditioner : "  << this->preconditioner() << "\n";
    switch(this->preconditioner())
    {
        case NEUMANN_DIRICHLET:
            // Neumann-Dirichlet preconditioner
            invSfPrime(-1*_res, _linearRelTol, muF);
            break;
        case DIRICHLET_NEUMANN:
            // Dirichlet-Neumann preconditioner
            invSsPrime(-1*_res, _linearRelTol, muS);
            break;
        case NEUMANN_NEUMANN:
            // Neumann-Neumann preconditioner
        {
            invSsPrime(-1*_res, _linearRelTol, muS);
            invSfPrime(-1*_res, _linearRelTol, muF);

            std::cout << "norm_inf muS = " << norm_inf(muS) << std::endl;
            std::cout << "norm_inf muF = " << norm_inf(muF) << std::endl;
        }
        break;
        case NEWTON:
            // Dirichlet-Neumann preconditioner
            invSfSsPrime(_res, 0.001, muS);
            break;
        case NO_PRECONDITIONER:
        default:
        {
            std::ostringstream _ex;
            _ex << "The steklovPoincare operator needs a preconditioner : \n"
                << "NEUMANN_DIRICHLET, NEUMANN_NEUMANN, DIRICHLET_NEUMANN\n";
            throw std::logic_error( _ex.str() );
        }
    }

    if (this->preconditioner() != NEWTON)
    {
        if (M_nbEval == 1) M_aitkFS.restart();
        muk = M_aitkFS.computeDeltaLambda(M_dispStructOld, muF, muS );
    }
    else
    {
        muk = muS;
    }
}


void steklovPoincare::solveLinearFluid()
{
    this->M_fluid->iterateLin(time(), *M_BCh_du);
}


void steklovPoincare::solveInvLinearFluid()
{
    this->M_fluid->iterateLin(time(), *M_BCh_du_inv);
    this->M_dzFluid = M_fluid->getDeltaLambda();
}


//

void steklovPoincare::solveLinearSolid()
{
    M_rhs_dz  = ZeroVector( M_rhs_dz.size() );
    M_dzSolid = ZeroVector( M_dzSolid.size() );

    if ( !M_BCh_dz->bdUpdateDone() )
        M_BCh_dz->bdUpdate(this->M_solid->mesh(),
                           this->M_solid->feBd(),
                           this->M_solid->dof());

    bcManageVector(M_rhs_dz,
                   this->M_solid->mesh(),
                   this->M_solid->dof(),
                   *M_BCh_dz,
                   this->M_solid->feBd(),
                   1., 1.);

    Real tol       = 1.e-10;

    std::cout << "rhs_dz norm   = " << norm_inf(M_rhs_dz) << std::endl;
    this->M_solid->setRecur(1);
    this->M_solid->updateJac(M_dzSolid, 0);
    this->M_solid->solveLin(M_dzSolid, M_rhs_dz, tol, *M_BCh_dz);
    std::cout << "dz norm       = " << norm_inf(M_dzSolid) << std::endl;
    std::cout << "residual norm = " << norm_inf(M_solid->residual()) << std::endl;
}


void steklovPoincare::solveInvLinearSolid()
{
    M_rhs_dz  = ZeroVector( M_rhs_dz.size() );
    M_dzSolid = ZeroVector( M_dzSolid.size() );

    if ( !M_BCh_dz_inv->bdUpdateDone() )
        M_BCh_dz_inv->bdUpdate(this->M_solid->mesh(),
                               this->M_solid->feBd(),
                               this->M_solid->dof());

    bcManageVector(M_rhs_dz,
                   this->M_solid->mesh(),
                   this->M_solid->dof(),
                   *M_BCh_dz_inv,
                   this->M_solid->feBd(),
                   1., 1.);

    Real tol       = 1.e-10;

    std::cout << "rhs_dz norm   = " << norm_inf(M_rhs_dz) << std::endl;
    this->M_solid->setRecur(1);
    this->M_solid->updateJac(M_dzSolid, 0);
    this->M_solid->solveLin(M_dzSolid, M_rhs_dz, tol, *M_BCh_dz_inv);
    std::cout << "dz norm       = " << norm_inf(M_dzSolid) << std::endl;
    std::cout << "residual norm = " << norm_inf(M_solid->residual()) << std::endl;
}


//
//
//


void  steklovPoincare::invSfPrime(const Vector& res,
                                  double linear_rel_tol,
                                  Vector& step)
{
    setResidualFSI(res);

    M_solid->d() = ZeroVector(M_solid->d().size());
    std::cout << "norm_inf residual FSI = " << norm_inf(M_residualFSI);
//    this->M_fluid.updateDispVelo();
    solveLinearFluid();

    Vector deltaLambda = this->M_fluid->getDeltaLambda();

    transferOnInterface(deltaLambda,
                        M_fluid->BC_fluid(),
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

    transferOnInterface(M_dzSolid,
                        M_solid->BC_solid(),
                        "Interface",
                        step);
}


void steklovPoincare::invSfSsPrime(const Vector& _res,
                                   double _linearRelTol,
                                   Vector& _muk)
{
        // AZTEC specifications for the second system
    int    data_org[AZ_COMM_SIZE];   // data organisation for J
    int    proc_config[AZ_PROC_SIZE];  // Processor information:
    int    options[AZ_OPTIONS_SIZE];   // Array used to select solver options.
    double params[AZ_PARAMS_SIZE];     // User selected solver paramters.
    double status[AZ_STATUS_SIZE];     // Information returned from AZ_solve()

    AZ_set_proc_config(proc_config, AZ_NOT_MPI);

    // data_org assigned "by hands": no parallel computation is performed
    UInt dim_res = _res.size();
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
    params[AZ_tol]       = _linearRelTol;

    //AZTEC matrix for the jacobian
    AZ_MATRIX *J;
    J = AZ_matrix_create(dim_res);

    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:
    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecSfSsPrimePrec);
//    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecSfSsPrime);

    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dim_res; ++i)
        _muk[i]=0.0;

    chrono.start();
    AZ_iterate(&_muk[0], const_cast<double*>( &_res[0] ), options, params,
               status, proc_config, J, NULL, NULL);

    transferOnInterface(DDNprecond(_muk),
                        M_solid->BC_solid(),
                        "Interface",
                        _muk);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    AZ_matrix_destroy(&J);
}

void my_matvecSfSsPrime(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
{
    // Extraction of data from J
    steklovPoincare::DataJacobian* my_data = static_cast< steklovPoincare::DataJacobian* >(AZ_get_matvec_data(J));

    UInt dim = my_data->M_pFS->dzSolid().size();

    double xnorm =  AZ_gvector_norm(dim, -1, z, proc_config);
    std::cout << " ***** norm (z)  = " << xnorm << std::endl<< std::endl;

    Vector jz(dim);
    Vector zSolid(dim);

    if ( xnorm == 0.0 )
        for (int i=0; i <(int) dim; ++i)
            {
                Jz[i]  = 0.0;
            }
    else
    {
        zSolid =  my_data->M_pFS->getSolidInterfaceOnSolid(z);

        for (int i = 0; i <(int) dim; ++i)
        {
            my_data->M_pFS->solid().d()[i] =  zSolid[i];
        }

        my_data->M_pFS->fluid().updateDispVelo();
        my_data->M_pFS->solveLinearFluid();

        my_data->M_pFS->solveLinearSolid();

        my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
        my_data->M_pFS->setResidualF(my_data->M_pFS->fluid().residual());

        my_data->M_pFS->getResidualFSIOnSolid(jz);

        for (int i = 0; i < (int) dim; ++i)
            Jz[i] =  jz[i];
    }

    std::cout << " ***** norm (Jz) = "
              << AZ_gvector_norm(dim, -1, Jz, proc_config)
              << std::endl << std::endl;
}


void my_matvecSfSsPrimePrec(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
{
    // Extraction of data from J
    steklovPoincare::DataJacobian* my_data = static_cast< steklovPoincare::DataJacobian* >(AZ_get_matvec_data(J));

    UInt dim = my_data->M_pFS->dzSolid().size();

    double xnorm =  AZ_gvector_norm(dim, -1, z, proc_config);
    std::cout << " ***** norm (z)  = " << xnorm << std::endl<< std::endl;

    Vector jz(dim);
    Vector zSolid(dim);

    for (int ii = 0; ii < dim; ++ii)
        zSolid[ii] = z[ii];

    if ( xnorm == 0.0 )
        for (int i=0; i <(int) dim; ++i)
            {
                Jz[i]     = 0.0;
            }
    else
    {
        my_data->M_pFS->solid().d() = my_data->M_pFS->DDNprecond(zSolid);

        my_data->M_pFS->fluid().updateDispVelo();
        my_data->M_pFS->solveLinearFluid();

        my_data->M_pFS->solveLinearSolid();

        my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
        my_data->M_pFS->setResidualF(my_data->M_pFS->fluid().residual());

        my_data->M_pFS->getResidualFSIOnSolid(jz);

        for (int i = 0; i < (int) dim; ++i)
            Jz[i] =  jz[i];
    }

    std::cout << " ***** norm (Jz) = "
              << AZ_gvector_norm(dim, -1, Jz, proc_config)
              << std::endl << std::endl;
}


Vector steklovPoincare::DDNprecond(Vector const &_z)
{
    std::cout << "Domain Decompostion Newton Precond. using ";
    Vector Pz(_z.size());

    setResidualFSI(_z);

    this->M_dzSolid = ZeroVector(M_dzSolid.size());
    this->M_dzFluid = ZeroVector(M_dzFluid.size());

    switch(this->DDNpreconditioner())
    {
        case NEUMANN_DIRICHLET:
            // Neumann-Dirichlet preconditioner
            std::cout << " Neumann-Dirichlet Precond ... \n" << std::endl;
            solveInvLinearFluid();
            Pz = getFluidInterfaceOnSolid(M_dzFluid);
            break;
        case DIRICHLET_NEUMANN:
            // Dirichlet-Neumann preconditioner
            std::cout << " Dirichlet-Neumann Precond ... \n" << std::endl;
            solveInvLinearSolid();
            Pz = M_dzSolid;
            break;
        case NEUMANN_NEUMANN:
            // Neumann-Neumann preconditioner
        {
            std::cout << " Neumann-Neumann Precond ... \n" << std::endl;
            solveInvLinearFluid();
            solveInvLinearSolid();
            Pz = 0.9*M_dzSolid  +
                0.1*getFluidInterfaceOnSolid(M_dzFluid);
        }
        break;
        case NO_PRECONDITIONER:
        default:
        {
            std::cout << " no preconditioner ... " << std::endl;
            Pz = _z;
        }
    }

    return Pz;
}


void steklovPoincare::computeResidualFSI()
{
    M_residualFSI = ZeroVector( M_residualFSI.size() );

    FOR_EACH_INTERFACE_DOF( M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                            M_residualF[IDfluid - 1 + jDim*totalDofFluid] -
                            M_residualS[IDsolid - 1 + jDim*totalDofSolid] );
}


void steklovPoincare::setResidualFSI(double const* _res)
{
    FOR_EACH_INTERFACE_DOF( M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                            _res[IDsolid - 1 + jDim*totalDofSolid] );

}

void steklovPoincare::setResidualFSI( Vector const&  _res)
{
    FOR_EACH_INTERFACE_DOF( M_residualFSI[IDfluid - 1 + jDim*totalDofFluid] =
                            _res[IDsolid - 1 + jDim*totalDofSolid] );

}


Vector steklovPoincare::getResidualFSIOnSolid()
{
    Vector vec(M_residualS.size());

    FOR_EACH_INTERFACE_DOF( vec[IDsolid - 1 + jDim*totalDofSolid] =
                            M_residualF[IDfluid - 1 + jDim*totalDofFluid] -
                            M_residualS[IDsolid - 1 + jDim*totalDofSolid] );
    return vec;
}


void steklovPoincare::getResidualFSIOnSolid(Vector& _vec)
{
    FOR_EACH_INTERFACE_DOF( _vec[IDsolid - 1 + jDim*totalDofSolid] =
                            M_residualF[IDfluid - 1 + jDim*totalDofFluid] -
                            M_residualS[IDsolid - 1 + jDim*totalDofSolid] );
}


Vector steklovPoincare::getSolidInterfaceOnFluid(Vector const& _vec)
{
    Vector vec(M_residualF.size());

    FOR_EACH_INTERFACE_DOF( vec[IDfluid - 1 + jDim*totalDofFluid] =
                            _vec[IDsolid - 1 + jDim*totalDofSolid] );
    return vec;
}

Vector steklovPoincare::getSolidInterfaceOnSolid(double const* _vec)
{
    Vector vec(M_residualF.size());

    FOR_EACH_INTERFACE_DOF( vec[IDsolid - 1 + jDim*totalDofSolid] =
                            _vec[IDsolid - 1 + jDim*totalDofSolid] );
    return vec;
}


Vector steklovPoincare::getFluidInterfaceOnSolid(Vector const& _vec)
{
    Vector vec(M_residualS.size());

    FOR_EACH_INTERFACE_DOF( vec[IDsolid - 1 + jDim*totalDofSolid] =
                            _vec[IDfluid - 1 + jDim*totalDofSolid] );

    return vec;
}

//
// add steklovPoincare to factory
//
namespace
{
operFS* createSP(){ return new steklovPoincare(); }
static bool reg = FSIFactory::instance().registerProduct( "steklovPoincare", &createSP );
}
}
