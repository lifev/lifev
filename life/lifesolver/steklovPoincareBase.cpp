/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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


#include <life/lifesolver/steklovPoincareBase.hpp>
#include <life/lifesolver/reducedLinFluid.hpp>


namespace LifeV
{
Real fzeroSP(const Real& /*t*/,
             const Real& /*x*/,
             const Real& /*y*/,
             const Real& /*z*/,
             const ID& /*i*/)
{return 0.0;}

steklovPoincare::steklovPoincare():
        super                                ( ),
        dz                                   ( ),
        M_defOmega                           ( 0.005 ),
        M_defOmegaS                          ( 0.005 ),
        M_defOmegaF                          ( 0.005 ),
        M_aitkFS                             ( ),
{
    this->setPreconditioner( NEWTON );
    this->setDDNPreconditioner( DDN_DIRICHLET_NEUMANN );
//     M_aitkFS.setDefault( M_defOmegaS, M_defOmegaF );
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

    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) OmegaS        = " << M_defOmegaS << "\n";
    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) OmegaF        = " << M_defOmegaF << "\n";

    M_aitkFS.setDefault(M_defOmegaS, M_defOmegaF);

    this->setPreconditioner   (  ( Preconditioner )      data("problem/precond"     , NEWTON ) );
    this->setDDNPreconditioner(  ( DDNPreconditioner )data("problem/DDNprecond"  , DDN_DIRICHLET_NEUMANN ) );

    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) prec          = " << this->preconditioner() << "\n";
    Debug( 6205 ) << "steklovPoincare::setDataFromGetPot(GetPot) Newton-prec   = " << this->DDNpreconditioner() << "\n";
}

void
steklovPoincare::setup()
{

    // call FSIOperator setup()
    setLinearFluid(true);
    setLinearSolid(true);

    super::setup();

    M_dz.reset    (new vector_Type(*this->M_solidInterfaceMap));

    M_epetraOper.reset( new Epetra_SteklovPoincare(this));

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_Type(this->M_fluid->getMap()));
        M_beta.reset  (new vector_Type(this->M_fluid->getMap()));
    }


}



//
// Residual computation
//





void steklovPoincare::eval(const vector_Type& disp,
                           int                status,
                           vector_Type&       dispNew,
                           vector_Type&       velo)
{

    if (status)
    {
        this->M_fluid->resetPrec();
        //this->M_solid->resetPrec();
    }


    MPI_Barrier(MPI_COMM_WORLD);

    // possibly unsafe when using more cpus, since both has repeated maps
    *M_lambdaFluid       = _disp;
    *M_lambdaSolid       = _disp

                           *this->M_sigmaFluid *= 0.;

    if (this->isFluid())
    {
        this->M_meshMotion->iterate();

        vector_Type const meshDisplacement( M_meshMotion->displacement(), this->M_meshMotion->getRepeatedEpetraMap() );
        this->moveMesh(meshDisplacement);

        this->transferMeshMotionOnFluid(meshDisplacement,
                                        *this->M_veloFluidMesh);

        *this->M_veloFluidMesh    -= *M_dispFluidMeshOld;
        *this->M_veloFluidMesh    *= 1./(M_dataFluid->getTimeStep());

        // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
        // the value of the displacement for some points

        this->moveMesh(meshDisplacement);

        *M_beta *= 0.;

        vector_Type const meshDispDiff( M_meshMotion->dispDiff(), this->M_meshMotion->getRepeatedEpetraMap() );

        this->interpolateVelocity(meshDispDiff, *M_beta);

        *M_beta *= -1./M_dataFluid->getTimeStep();

        *M_beta  += *this->M_un;

        double alpha = 1./M_dataFluid->getTimeStep();

        *M_rhsNew   = *this->M_rhs;
        *M_rhsNew  *= alpha;

        this->M_fluid->updateSystem( alpha, *M_beta, *M_rhsNew );

        this->M_fluid->iterate( *M_BCh_u );

        this->transferFluidOnInterface(this->M_fluid->residual(), *this->M_sigmaFluid);

    }


    if ( true && this->isFluid() )
    {
        vector_Type vel  (this->fluid().velocityFESpace().map());
        vector_Type press(this->fluid().pressureFESpace().map());

        vel.subset(this->M_fluid->solution());
        press.subset(this->M_fluid->solution(), this->fluid().velocityFESpace().dim()*this->fluid().pressureFESpace().fieldDim());


        std::cout << "norm_inf( vel ) " << vel.NormInf() << std::endl;
        std::cout << "norm_inf( press ) " << press.NormInf() << std::endl;
        this->M_fluid->postProcess();
    }



    if (this->isSolid())
    {
        this->M_solid->iterate( *M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(),      *this->M_lambdaSolid);
        this->transferSolidOnInterface(this->M_solid->vel(),       *this->M_lambdaDotSolid);
        this->transferSolidOnInterface(this->M_solid->residual() , *this->M_sigmaSolid);
    }


    if ( true && this->isSolid() )
    {
        this->solid().postProcess();
    }

// possibly unsafe when using more cpus, since both has repeated maps


    MPI_Barrier(MPI_COMM_WORLD);


    std::cout << " ::: norm(disp     )    = " << _disp.NormInf() << std::endl;
    std::cout << " ::: norm(dispNew  )    = " << this->M_lambdaSolid->NormInf() << std::endl;
    std::cout << " ::: norm(velo     )    = " << this->M_lambdaDotSolid->NormInf() << std::endl;
    std::cout << " ::: max Residual Fluid = " << M_sigmaFluid->NormInf() << std::endl;
    std::cout << " ::: max Residual Solid = " << M_sigmaSolid->NormInf() << std::endl;

    if (this->isFluid())
        std::cout << "Max ResidualF        = " << M_fluid->residual().NormInf() << std::endl;
    if (this->isSolid())
    {
        std::cout << "NL2 DiplacementS     = " << M_solid->disp().Norm2() << std::endl;
        std::cout << "Max ResidualS        = " << M_solid->residual().NormInf() << std::endl;
    }

}

void steklovPoincare::evalResidual(Vector       &res,
                                   const Vector &disp,
                                   const int     iter)
{
    int status = 0;

//     M_solid->postProcess();
//     M_fluid->postProcess();

    if (iter == 0) status = 1;

    std::cout << "*** Residual computation g(x_" << iter <<" )";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    M_dispStructOld = disp;

    eval(disp, status, M_dispStruct, M_velo);

//    M_residualS = M_solid->residual();
    M_residualF = M_fluid->residual();

    computeResidualFSI();

    res = M_interfaceStress;

    std::cout << "max ResidualF   = " << norm_inf(M_residualF)
              << std::endl;
    std::cout << "max ResidualS   = " << norm_inf(M_residualS)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_interfaceStress)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_strongResidualFSI)
              << std::endl;
}



//
// new step computation resolution
//


void  steklovPoincare::solveJac(vector_Type        &muk,
                                const vector_Type  &_res,
                                double        _linearRelTol)
{
    vector_Type muF(_res.size());
    vector_Type muS(_res.size());

    Debug(  6215 ) << "steklovPoincare::solveJac _linearRelTol  : " << _linearRelTol << "\n";
    Debug(  6215 ) << "steklovPoincare::solveJac preconditioner : "  << this->preconditioner() << "\n";

    switch (this->preconditioner())
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
        invSfSsPrime(_res, _linearRelTol, muS);
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
        M_aitkFS.restart();
        muk = M_aitkFS.computeDeltaLambdaFSI(M_dispStructOld, muF, muS );
    }
    else
    {
        muk = muS;
    }

    M_reducedLinFluid->setComputedMatrix(false);
}


void steklovPoincare::solveLinearFluid()
{
    this->M_fluid->solveLinearSystem(time(), *M_BCh_du);
}


void steklovPoincare::solveInvLinearFluid()
{
    this->M_fluid->solveLinearSystem(time(), *M_BCh_du_inv);
    this->M_dzFluid = M_fluid->getDeltaLambda();
}


//

void steklovPoincare::solveLinearSolid()
{
    this->M_solid->setRecur(1);
    this->M_solid->updateJacobian(M_dzSolid, 0);
    this->M_solid->solveJacobian(0., M_BCh_dz);
    M_dzSolid = this->M_solid->ddisp();

    Debug(6215) << "dz norm       = " << norm_2(M_dzSolid) << "\n";
    std::cout << "S-  norm_inf residual = " << norm_inf(M_solid->residual()) << "\n";
}


void steklovPoincare::solveInvLinearSolid()
{
    this->M_solid->setRecur(1);
    this->M_solid->updateJacobian(M_dzSolid, 0);
    this->M_solid->solveJacobian(0., M_BCh_dz_inv);
    M_dzSolid = this->M_solid->ddisp();
    Debug(  6215 ) << "dz norm       = " << norm_2(M_dzSolid) << "\n";
    Debug(  6215 ) << "residual norm = " << norm_2(M_solid->residual()) << "\n";
}


//
//
//


void  steklovPoincare::invSfPrime(const vector_Type& res,
                                  double /*linear_rel_tol*/,
                                  vector_Type& step)
{
    setResidualFSI(res);

    M_solid->disp() = ZeroVector(M_solid->disp().size());

//    solveLinearFluid();

    M_reducedLinFluid->solveInvReducedLinearFluid();

    //Vector deltaLambda = this->M_fluid->getDeltaLambda();

    double dt  = this->M_fluid->getTimeStep();
    double rho = this->M_fluid->density();

    vector_Type deltaLambda = dt*dt/rho*M_reducedLinFluid->residual();

    transferOnInterface(deltaLambda,
                        M_fluid->bcHandler(),
                        "Interface",
                        step);

    std::cout << "norm_2 deltaLambda = " << norm_2(deltaLambda) << std::endl;
    std::cout << "norm_2 step        = " << norm_2(step) << std::endl;
}


//
//
//

void  steklovPoincare::invSsPrime(const vector_Type& res,
                                  double /*linear_rel_tol*/,
                                  vector_Type& step)
{
    setResidualFSI(res);
    solveLinearSolid();

    transferOnInterface(M_dzSolid,
                        M_solid->BCh_solid(),
                        "Interface",
                        step);
}


void steklovPoincare::invSfSsPrime(const vector_Type& _res,
                                   double _linearRelTol,
                                   vector_Type& _muk)
{
    // AZTEC specifications for the second system
    int    data_org   [AZ_COMM_SIZE];     // data organisation for J
    int    proc_config[AZ_PROC_SIZE];  // Processor information:
    int    options    [AZ_OPTIONS_SIZE];   // Array used to select solver options.
    double params     [AZ_PARAMS_SIZE];     // User selected solver paramters.
    double status     [AZ_STATUS_SIZE];     // Information returned from AZ_solve()

    AZ_set_proc_config(proc_config, AZ_NOT_MPI);

    // data_org assigned "by hands": no parallel computation is performed
    UInt dim_res = _res.size();
    data_org[AZ_N_internal] = dim_res;
    data_org[AZ_N_border]   = 0;
    data_org[AZ_N_external] = 0;
    data_org[AZ_N_neigh]    = 0;

    // Recovering AZTEC defaults options and params
    AZ_defaults(options,params);

    // Fixed Aztec options for this linear system
    options[AZ_solver]   = AZ_gmres;
    options[AZ_output]   = 1;
    options[AZ_poly_ord] = 5;
    options[AZ_kspace]   = 40;
    options[AZ_conv]     = AZ_rhs;
    params [AZ_tol]      = _linearRelTol;

    //AZTEC matrix for the jacobian
    AZ_MATRIX *J;
    J = AZ_matrix_create(dim_res);

    // data containing the matrices C, D, trD and H as pointers
    // are passed through A_ii and pILU_ii:

    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecSfSsPrime);

    std::cout << "  N-  Solving Jacobian system... ";
    Chrono chrono;

    vector_Type res = _res;

//     for (UInt ii = 0; ii < 20; ++ii)
//     {
//         res[ii] = 0;
//         res[ii + M_interfaceNbreDof] = 0;
//         res[ii + 2*M_interfaceNbreDof] = 0;
//         res[M_interfaceNbreDof - ii - 1] = 0;
//         res[2*M_interfaceNbreDof - ii - 1] = 0;
//         res[3*M_interfaceNbreDof - ii - 1] = 0;
//     }

    for (UInt i=0; i<dim_res; ++i)
    {
        _muk[i]=0.0;
    }
    chrono.start();


    std::cout << "norm 2 muk = " << norm_2(_muk) << std::endl;
    std::cout << "norm 2 res = " << norm_2(_res) << std::endl;

    AZ_iterate(&_muk[0], const_cast<double*>( &res[0] ), options, params,
               status, proc_config, J, NULL, NULL);

    _muk = DDNprecond(_muk);

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    AZ_matrix_destroy(&J);
}


void my_matvecSfSsPrime(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
{
    // Extraction of data from J
    steklovPoincare::DataJacobian* my_data = static_cast< steklovPoincare::DataJacobian* >(AZ_get_matvec_data(J));

    UInt dim = my_data->M_pFS->displacement().size();

    double xnorm =  AZ_gvector_norm(dim, -1, z, proc_config);
    std::cout << " ***** norm (z)  = " << xnorm << std::endl<< std::endl;

//    Vector jz(dim);
//    jz = ZeroVector(dim);

    vector_Type zSolid(dim);

    for (UInt ii = 0; ii < dim; ++ii)
    {
        zSolid[ii] = z[ii];
    }
    if ( xnorm == 0.0 )
        for (UInt i=0; i < dim; ++i)
        {
            Jz[i]     = 0.0;
        }
    else
    {
        if (!my_data->M_pFS->reducedFluid())
        {
            my_data->M_pFS->displacement() = my_data->M_pFS->DDNprecond(zSolid);

            if (my_data->M_pFS->fluidMpi())
            {
                my_data->M_pFS->fluid().updateDispVelo();
                my_data->M_pFS->solveLinearFluid();
                my_data->M_pFS->setResidualF(my_data->M_pFS->fluid().residual());
            }
            if (my_data->M_pFS->solidMpi())
            {
                my_data->M_pFS->solveLinearSolid();
                my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
            }
        }
        else
        {
            vector_Type da(dim);
            double dt = my_data->M_pFS->fluid().getTimeStep();
            double dti2 = 1.0/(dt*dt) ;

            vector_Type zSolidPrec(dim);
            zSolidPrec = my_data->M_pFS->DDNprecond(zSolid);

            if (my_data->M_pFS->fluidMpi())
            {
                da = - dti2*my_data->M_pFS->fluid().density()*zSolidPrec;

                if (my_data->M_pFS->nbEval() == 1) my_data->M_pFS->getReducedLinFluid()->setComputedMatrix(false);

                my_data->M_pFS->getReducedLinFluid()->setDacc(da);
                my_data->M_pFS->getReducedLinFluid()->solveReducedLinearFluid();
                my_data->M_pFS->setResidualF(my_data->M_pFS->getReducedLinFluid()->residual());
            }
            if (my_data->M_pFS->solidMpi())
            {
                my_data->M_pFS->displacement() = zSolidPrec;
                my_data->M_pFS->solveLinearSolid();

                my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
            }
        }
        my_data->M_pFS->computeResidualFSI();

        for (UInt i = 0; i < dim; ++i)
        {
            Jz[i] =  my_data->M_pFS->residual()[i];
        }
    }

    std::cout << " ***** norm (Jz) = "
              << AZ_gvector_norm(dim, -1, Jz, proc_config)
              << std::endl << std::endl;
}


vector_Type steklovPoincare::DDNprecond(vector_Type const &_z)
{
    std::cout << "DD-Newton Precond. using ";

    vector_Type Pz(_z.size());

    //setResidualFSI(_z);

    //M_interfaceStress = transferSolidOnInterface( _z );

    M_interfaceStress = _z;

    M_dzSolid = Zerovector_type(M_dzSolid.size());
    M_dzFluid = Zerovector_type(M_dzFluid.size());

    switch (this->DDNpreconditioner())
    {
    case NEUMANN_DIRICHLET:
        // Neumann-Dirichlet preconditioner
    {
        std::cout << " Neumann-Dirichlet Precond ... \n" << std::endl;

        M_reducedLinFluid->solveInvReducedLinearFluid();

        //Vector deltaLambda = this->M_fluid->getDeltaLambda();

        double dt  = this->M_fluid->getTimeStep();
        double rho = this->M_fluid->density();

        vector_Type deltaLambda = dt*dt/rho*M_reducedLinFluid->residual();

        transferOnInterface(deltaLambda,
                            M_fluid->bcHandler(),
                            "Interface",
                            Pz);

    }
    break;
    case DIRICHLET_NEUMANN:
        // Dirichlet-Neumann preconditioner
        std::cout << " Dirichlet-Neumann Precond ... \n" << std::flush << std::endl;
        if (this->mpi())
        {
            if (this->solidMpi())
            {
                solveInvLinearSolid();
                Pz = transferSolidOnInterface(M_dzSolid);
                MPI_Send(&Pz[0], Pz.size(), MPI_DOUBLE, FLUID, 0, MPI_COMM_WORLD);
            }
            if (this->fluidMpi())
            {
                MPI_Status err;
                MPI_Recv(&Pz[0], Pz.size(), MPI_DOUBLE, SOLID, 0, MPI_COMM_WORLD, &err);
            }
        }
        else
        {
            solveInvLinearSolid();
            Pz = transferSolidOnInterface(M_dzSolid);
        }
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

//     for (UInt ii = 0; ii < Pz.size(); ++ii)
//         std::cout << ii << " " << Pz[ii] << std::endl;

    std::cout << "ok" << std::endl;
    return Pz;
}


//
// Interface operators
//


void steklovPoincare::computeStrongResidualFSI()
{

    Chrono chrono, chronoloc;
    chrono.start();

    std::cout << "  SP- Computing strong residual ... " << std::endl;

    std::cout << "        builing the mass matrix ... " << std::flush;

    FSIOperator::dof_interface_type3D dofReducedFluidToMesh
    (new FSIOperator::dof_interface_type3D::element_type);

    //Real f(const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    BCFunctionBase f0( fzeroSP);

    BCHandler BCh;

    BCh.addBC("Interface", 1, Natural,   Full, f0, 3);
    BCh.addBC("interface", 1, Essential, Full, f0, 3);

    BCh.bdUpdate(this->fluid().mesh(),
                 this->fluid().feBd_u(),
                 this->fluid().uDof());

    BCBase const &BCbEss = BCh[1];

    // Number of local Dof (i.e. nodes) in this face
    UInt nDofF = this->fluid().feBd_u().nbNode;

    // Number of total scalar Dof
    UInt totalDof = this->fluid().uDof().numTotalDof();

    // Number of components involved in this boundary condition
    UInt nComp = BCbEss.numberOfComponents();

    for ( ID i = 1; i <= BCbEss.list_size(); ++i )
    {
        {
            // Global Dof
            UInt index = BCbEss( i ) ->id();// + ( BCbEss.component( j ) - 1 ) * totalDof;
            M_interfaceDisplacement[index] = i;
//             std::cout << index << " -> " << i << std::endl;
//            std::cout << M_dofReducedFluidToStructure->getInterfaceDof(51) << " -> " << i << std::endl;
        }
    }

//     typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type;

//     matrix_type mass    (3*nDofF, 3*nDofF);
//     Vector      rhs     (3*nDofF);
//     Vector      residual(3*nDofF);


//     BCBase const &BCb = BCh[0];

//     // Number of local Dof (i.e. nodes) in this face
//     UInt nDofF = this->fluid().feBd_u().nbNode;

//     // Number of total scalar Dof
//     UInt totalDof = this->fluid().uDof().numTotalDof();

//     // Number of components involved in this boundary condition
//     UInt nComp = BCb.numberOfComponents();

//     const IdentifierNatural* pId;

//     SimpleArray<UInt>  localToBoundary;
//     SimpleArray<UInt>  boundaryToLocal;

//     ID ibF, idDof, icDof, gDof;

//     typedef BoostMatrix<boost::numeric::ublas::column_major> matrix_type;

//     matrix_type mass(3*nDofF, 3*nDofF);
//     Vector      rhs;
//     Vector      residual;

//     SolverUMFPACK  solver;

//     solver.setMatrix(mass);
//     solver.solve(residual, rhs);

    chrono.stop();

    std::cout << "        total time                          " << chrono.diff() << " s." << std::endl;
}


//
// BC vector interface treatment
//

// Solid, Lin. Solid and inverses

void steklovPoincare::setSolidInterfaceDisp(Vector& disp,
                                            UInt type)
{
    M_bcvSolidInterfaceDisp->setup(disp,
                                   M_interfaceNbreDof,
                                   M_dofSolid,
                                   type);
}

void steklovPoincare::setSolidLinInterfaceDisp(Vector& disp,
                                               UInt type)
{
    M_bcvSolidLinInterfaceDisp->setup(disp,
                                      M_interfaceNbreDof,
                                      M_dofSolid,
                                      type);
}

void steklovPoincare::setSolidInvLinInterfaceStress(Vector &stress,
                                                    UInt type)
{
    M_bcvSolidInvLinInterfaceStress->setup(stress,
                                           M_interfaceNbreDof,
                                           M_dofSolid,
                                           type);
}


// void steklovPoincare::setSolidInterfaceStress(Vector& stress,
//                                                  UInt type)
// {
//     M_bcvSolidInterfaceDisplacement->setup(stress,
//                                            stress.size(),
//                                            M_dofSolid,
//                                            type);
// }
// void steklovPoincare::setSolidInterfaceDisp(Vector &disp,
//                                                UInt type)
// {
//     M_bcvSolidInterfaceDisplacement->setup(disp,
//                                            disp.size(),
//                                            M_dofSolid,
//                                            type);
// }


// Fluid, Lin. and inverses

void steklovPoincare::setFluidInterfaceDisp(Vector &disp,
                                            UInt type)
{
    M_bcvFluidInterfaceDisp->setup(disp,
                                   M_interfaceNbreDof,
                                   M_dofFluid,
                                   type);
}

void steklovPoincare::setFluidLinInterfaceVel(Vector &vel,
                                              UInt type)
{
    M_bcvFluidLinInterfaceVel->setup(vel,
                                     M_interfaceNbreDof,
                                     M_dofFluid,
                                     type);
}

// reduced fluid and inverse


void steklovPoincare::setReducedFluidInterfaceAcc(Vector &acc,
                                                  UInt type)
{
    M_bcvReducedFluidInterfaceAcc->setup(acc,
                                         M_interfaceNbreDof,
                                         M_dofFluid,
                                         type);
}

void steklovPoincare::setReducedFluidInvInterfaceAcc(Vector &acc,
                                                     UInt type)
{
    M_bcvReducedFluidInvInterfaceAcc->setup(acc,
                                            M_interfaceNbreDof,
                                            M_dofFluid,
                                            type);
}





//
// Moment projection on the interface
//

void steklovPoincare::computeResidualFSI()
{

    // Add the fluid and solid stress on the interface to get the residual

    std::map<ID, ID> fluidDofMap = M_dofFluid->locDofMap();
    std::map<ID, ID> solidDofMap = M_dofSolid->locDofMap();
    std::map<ID, ID> FSIDofMap   = M_dofStructureToHarmonicExtension->locDofMap();

    UInt totalDofFluid = M_fluid->uDof().numTotalDof();
    UInt totalDofSolid = M_solid->dDof().numTotalDof();


    int size = 3*M_interfaceNbreDof;

    M_interfaceStress = ZeroVector(size);

    double* residual;
    double* residualFS;

    residual   = (double*) malloc(size*sizeof(double));
    residualFS = (double*) malloc(size*sizeof(double));

    if (this->mpi())
    {
        if (this->fluidMpi())
        {
            FOR_EACH_DOF_INTERFACE(residualFS[localS + jDim*M_interfaceNbreDof - 1] =
                                       this->M_residualF[dofF + jDim*totalDofFluid - 1]);
        }
        if (this->solidMpi())
        {
            FOR_EACH_DOF_INTERFACE(residualFS[localS + jDim*M_interfaceNbreDof - 1] =
                                       - this->M_residualS[dofS + jDim*totalDofSolid - 1]);
        }

        MPI_Reduce( residualFS,
                    residual,
                    size,
                    MPI_DOUBLE,
                    MPI_SUM,
                    SOLID,
                    MPI_COMM_WORLD
                  );

        if (this->solidMpi())
        {
            for (UInt ii = 0; ii < size; ++ii)
            {
                M_interfaceStress[ii] = residual[ii];
            }
            MPI_Send(residual, size, MPI_DOUBLE, FLUID, 0, MPI_COMM_WORLD);
        }
        if (this->fluidMpi())
        {
            MPI_Status err;
            MPI_Recv(residual, size, MPI_DOUBLE, SOLID, 0, MPI_COMM_WORLD, &err);

            for (UInt ii = 0; ii < size; ++ii)
            {
                M_interfaceStress[ii] = residual[ii];
            }
        }
    }
    else
    {
//        std::cout << "transfer ..." << M_interfaceStress.size() << std::endl;

        Vector vecF = transferFluidOnInterface(M_residualF);
        Vector vecS = transferSolidOnInterface(M_residualS);

//         for (UInt ii = 0; ii < M_residualF.size(); ++ii)
//         {
//             std::cout << ii << " " << M_residualF[ii] << std::endl;
//         }

//         for (UInt ii = 0; ii < vecS.size(); ++ii)
//         {
//             std::cout << ii << " " << vecF[ii] << " " << vecS[ii] << std::endl;
//         }
//          for (UInt ii = 0; ii < 20; ++ii)
//          {
//              vecF[ii] = 0;
//              vecF[ii + M_interfaceNbreDof] = 0;
//              vecF[ii + 2*M_interfaceNbreDof] = 0;
//              vecF[M_interfaceNbreDof - ii - 1] = 0;
//              vecF[2*M_interfaceNbreDof - ii - 1] = 0;
//              vecF[3*M_interfaceNbreDof - ii - 1] = 0;
//          }
//         M_interfaceStress = transferFluidOnInterface(M_residualF) -
//             transferSolidOnInterface(M_residualS);
        M_interfaceStress = vecF - vecS;
//         FOR_EACH_DOF_INTERFACE(M_interfaceStress[localS + jDim*M_interfaceNbreDof - 1] =
//                     this->M_residualF[dofF + jDim*totalDofFluid - 1] -
//                     this->M_residualS[dofS + jDim*totalDofSolid - 1]);
    }

    for (UInt ii = 0; ii < 30; ++ii)
    {
        M_interfaceStress[ii] = 0;
        M_interfaceStress[ii + M_interfaceNbreDof] = 0;
        M_interfaceStress[ii + 2*M_interfaceNbreDof] = 0;
        M_interfaceStress[M_interfaceNbreDof - ii - 1] = 0;
        M_interfaceStress[2*M_interfaceNbreDof - ii - 1] = 0;
        M_interfaceStress[3*M_interfaceNbreDof - ii - 1] = 0;
    }

}

// void steklovPoincare::registerMyProducts( )
// {
// FSIFactory::instance().registerProduct( "steklovPoincare", &createSP );
// solid_raw_type::StructureSolverFactory::instance().registerProduct( "LinearVenantKirchhof", &createLinearStructure );
// solid_raw_type::StructureSolverFactory::instance().registerProduct( "NonLinearVenantKirchhof", &createNonLinearStructure );
// }

//
// add steklovPoincare to factory
//

namespace
{
FSIOperator* createSP() { return new steklovPoincare(); }
static bool reg = FSIOperator::FSIFactory::instance().registerProduct( "steklovPoincare", &createSP );
}


}
