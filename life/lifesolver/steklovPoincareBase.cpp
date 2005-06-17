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
    super(),
    M_dzSolid(),
    M_dzFluid(),
    M_rhs_dz(),
    M_residualS(),
    M_residualF(),
    M_residualFSI(),
    M_strongResidualFSI(),
    M_defOmega( 0.005 ),
    M_defOmegaS( 0.005 ),
    M_defOmegaF( 0.005 ),
    M_aitkFS(),
    M_dataJacobian( this )
{
    this->setPreconditioner( NEWTON );
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
    super::setup();

    M_dzSolid.resize( 3*M_solid->dDof().numTotalDof() );
    M_dzFluid.resize( 3*M_fluid->uDof().numTotalDof() );

    M_rhs_dz.resize( 3*M_solid->dDof().numTotalDof() );

    M_residualS.resize( M_solid->dDof().numTotalDof() );
    M_residualF.resize( M_fluid->uDof().numTotalDof() );
    M_residualFSI.resize( M_fluid->uDof().numTotalDof() );
    M_strongResidualFSI.resize( M_fluid->uDof().numTotalDof() );

    M_aitkFS.setup( 3*M_solid->dDof().numTotalDof() );

    M_reducedLinFluid.reset(new reducedLinFluid(this, M_fluid, M_solid));
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

    M_solid->disp() = disp;

    M_fluid->updateMesh(M_time);

//      M_solid->postProcess();
//      M_fluid->postProcess();

    M_fluid->iterate   (M_time);

    M_solid->setRecur(0);
    M_solid->iterate();

    dispNew = M_solid->disp();
    velo    = M_solid->w();

    std::cout << "                ::: norm(disp     ) = "
              << norm_inf(disp) << std::endl;
    std::cout << "                ::: norm(dispNew  ) = "
              << norm_inf(dispNew) << std::endl;
    std::cout << "                ::: norm(velo     ) = "
              << norm_inf(velo) << std::endl;
}

void steklovPoincare::evalResidual(Vector       &res,
                                   const Vector &disp,
                                   const int     iter)
{
    int status = 0;

//     M_solid->postProcess();
//     M_fluid->postProcess();

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

    computeStrongResidualFSI();

    std::cout << "max ResidualF   = " << norm_inf(M_residualF)
              << std::endl;
    std::cout << "max ResidualS   = " << norm_inf(M_residualS)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_residualFSI)
              << std::endl;
    std::cout << "max ResidualFSI = " << norm_inf(M_strongResidualFSI)
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

    setBC();

    if (this->preconditioner() == NEWTON)
    {
        setInterfaceNewtonBC();
    }
    else
    {
        setInterfaceBC();
    }
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

    Debug(  6215 ) << "steklovPoincare::solveJac _linearRelTol  : " << _linearRelTol << "\n";
    Debug(  6215 ) << "steklovPoincare::solveJac preconditioner : "  << this->preconditioner() << "\n";

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
        if (M_nbEval == 1) M_aitkFS.restart();
        muk = M_aitkFS.computeDeltaLambda(M_dispStructOld, muF, muS );
    }
    else
    {
        muk = muS;
    }

    M_reducedLinFluid->setComputedMatrix(false);
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
    this->M_solid->setRecur(1);
    this->M_solid->updateJacobian(M_dzSolid, 0);
    this->M_solid->solveJacobian(0., M_BCh_dz);
    M_dzSolid = this->M_solid->ddisp();
    Debug(6215) << "dz norm       = " << norm_inf(M_dzSolid) << "\n";
    std::cout << "S-  norm_inf residual = " << norm_inf(M_solid->residual()) << "\n";
}


void steklovPoincare::solveInvLinearSolid()
{
    this->M_solid->setRecur(1);
    this->M_solid->updateJacobian(M_dzSolid, 0);
    this->M_solid->solveJacobian(0., M_BCh_dz_inv);
    M_dzSolid = this->M_solid->ddisp();
    Debug(  6215 ) << "dz norm       = " << norm_inf(M_dzSolid) << "\n";
    Debug(  6215 ) << "residual norm = " << norm_inf(M_solid->residual()) << "\n";
}


//
//
//


void  steklovPoincare::invSfPrime(const Vector& res,
                                  double /*linear_rel_tol*/,
                                  Vector& step)
{
    setResidualFSI(res);

    M_solid->disp() = ZeroVector(M_solid->disp().size());

//    solveLinearFluid();

    M_reducedLinFluid->solveInvReducedLinearFluid();

    //Vector deltaLambda = this->M_fluid->getDeltaLambda();

    double dt  = this->M_fluid->timestep();
    double rho = this->M_fluid->density();

    Vector deltaLambda = dt*dt/rho*M_reducedLinFluid->residual();

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



void  steklovPoincare::invSsPrime(const Vector& res,
                                  double /*linear_rel_tol*/,
                                  Vector& step)
{
    setResidualFSI(res);
    solveLinearSolid();

    transferOnInterface(M_dzSolid,
                        M_solid->BCh_solid(),
                        "Interface",
                        step);
}


void steklovPoincare::invSfSsPrime(const Vector& _res,
                                   double _linearRelTol,
                                   Vector& _muk)
{
    // AZTEC specifications for the second system
    int    data_org[AZ_COMM_SIZE];     // data organisation for J
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

    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecSfSsPrime);

    std::cout << "  N-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dim_res; ++i)
        _muk[i]=0.0;

    chrono.start();
    AZ_iterate(&_muk[0], const_cast<double*>( &_res[0] ), options, params,
               status, proc_config, J, NULL, NULL);

    transferOnInterface(DDNprecond(_muk),
                        M_solid->BCh_solid(),
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
    jz = ZeroVector(dim);

    Vector zSolid(dim);

    for (int ii = 0; ii < (int) dim; ++ii)
        zSolid[ii] = z[ii];

    if ( xnorm == 0.0 )
        for (int i=0; i <(int) dim; ++i)
            {
                Jz[i]     = 0.0;
            }
    else
    {
        if (!my_data->M_pFS->reducedFluid())
        {
            my_data->M_pFS->solid().disp() = my_data->M_pFS->DDNprecond(zSolid);

            my_data->M_pFS->fluid().updateDispVelo();
            my_data->M_pFS->solveLinearFluid();

            my_data->M_pFS->solveLinearSolid();

            my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
            my_data->M_pFS->setResidualF(my_data->M_pFS->fluid().residual());
        }
        else
        {
            Vector da(dim);
            double dt = my_data->M_pFS->fluid().timestep();
            double dti2 = 1.0/(dt*dt) ;

            Vector zSolidPrec(dim);
            zSolidPrec = my_data->M_pFS->DDNprecond(zSolid);
            da = - dti2*my_data->M_pFS->fluid().density()*zSolidPrec;

            if (my_data->M_pFS->nbEval() == 1) my_data->M_pFS->getReducedLinFluid()->setComputedMatrix(false);

            my_data->M_pFS->getReducedLinFluid()->setDacc(da);
            my_data->M_pFS->getReducedLinFluid()->solveReducedLinearFluid();

            my_data->M_pFS->solid().disp() = zSolidPrec;
            my_data->M_pFS->solveLinearSolid();

            my_data->M_pFS->setResidualS(my_data->M_pFS->solid().residual());
            my_data->M_pFS->setResidualF(my_data->M_pFS->getReducedLinFluid()->residual());
        }

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
    std::cout << "DD-Newton Precond. using ";

    Vector Pz(_z.size());

    setResidualFSI(_z);

    this->M_dzSolid = ZeroVector(M_dzSolid.size());
    this->M_dzFluid = ZeroVector(M_dzFluid.size());

    switch(this->DDNpreconditioner())
    {
        case NEUMANN_DIRICHLET:
            // Neumann-Dirichlet preconditioner
            {
                std::cout << " Neumann-Dirichlet Precond ... \n" << std::endl;
//                 solveInvLinearFluid();
//                 Vector deltaLambda = this->M_fluid->getDeltaLambda();
//                 transferOnInterface(deltaLambda,
//                                     M_fluid->BCh_fluid(),
//                                     "Interface",
//                                     Pz);

                M_reducedLinFluid->solveInvReducedLinearFluid();

                //Vector deltaLambda = this->M_fluid->getDeltaLambda();

                double dt  = this->M_fluid->timestep();
                double rho = this->M_fluid->density();

                Vector deltaLambda = dt*dt/rho*M_reducedLinFluid->residual();

                transferOnInterface(deltaLambda,
                                    M_fluid->bcHandler(),
                                    "Interface",
                                    Pz);

            }
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

    std::cout << "ok" << std::endl;
    return Pz;
}


//
// Boundary conditions setup
//

void steklovPoincare::setBC()
{
//     UInt dim_solid        = this->M_solid->dDof().numTotalDof();
//     UInt dim_fluid        = this->M_fluid->uDof().numTotalDof();

    // Boundary conditions for du and inverse

    BCFunctionBase bcf(fzeroSP);

    setStructureDispToSolid            (this->M_solid->disp());
    setStructureDispToHarmonicExtension(this->M_solid->disp());

    M_BCh_mesh->addBC("Interface", 1, Essential, Full,
                      *bcvStructureDispToHarmonicExtension(), 3);
    M_BCh_d->addBC("Interface", 1, Essential, Full,
                   *bcvStructureDispToSolid(), 3);


//     M_BCh_mesh->bdUpdate(this->M_fluid->mesh(),
//                          this->M_fluid->feBd_u(),
//                          this->M_fluid->uDof());
    //    COUPLED FSI LINEARIZED OPERATORS
    //
    // Passing the residue to the linearized fluid: \sigma -> du
    //
    // rem: for now: no fluid.dwInterpolated().
    //      In the future this could be relevant

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

    // solid acceleration
    // Boundary conditions for dp

     if (reducedFluid())
     {
         setDerStructureAccToReducedFluid(M_reducedLinFluid->dacc(), 2);
         M_BCh_dp->addBC("Wall",        1, Natural,   Scalar, //da_wall);
                         *bcvDerStructureAccToReducedFluid());
         M_BCh_dp->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
         M_BCh_dp->addBC("InFlow",      2, Essential, Scalar, bcf);
         M_BCh_dp->addBC("OutFlow",     3, Essential, Scalar, bcf);

         setDerReducedFluidLoadToStructure(M_strongResidualFSI);
         M_BCh_dp_inv->addBC("Wall",        1, Essential, Scalar,//dr_wall);
                             *bcvDerReducedFluidLoadToStructure());
         M_BCh_dp_inv->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
         M_BCh_dp_inv->addBC("InFlow",      2, Essential, Scalar, bcf);
         M_BCh_dp_inv->addBC("OutFlow",     3, Essential, Scalar, bcf);

         M_reducedLinFluid->setUpBC(M_BCh_dp);
         M_reducedLinFluid->setUpInvBC(M_BCh_dp_inv);
     }
}


void steklovPoincare::setInterfaceBC()
{
    setDerFluidLoadToFluid(residualFSI());
    M_BCh_du->addBC("Wall"     , 1, Natural  , Full,
                    *bcvDerFluidLoadToFluid(), 3);

    setDerFluidLoadToStructure(residualFSI());
    M_BCh_dz->addBC("Interface", 1, Natural  , Full,
                    *bcvDerFluidLoadToStructure(), 3);
}



void steklovPoincare::setInterfaceNewtonBC()
{
    std::cout << "Steklov-Poincare: NEWTON boundary conditions" << std::flush << std::endl;

    setDerHarmonicExtensionVelToFluid(this->M_fluid->dwInterpolated());
    M_BCh_du->addBC("Wall",   1,  Essential, Full,
                    *bcvDerHarmonicExtensionVelToFluid(), 3);

    setDerStructureDispToSolid(M_solid->disp());
    M_BCh_dz->addBC("Interface", 1, Essential , Full,
                    *bcvDerStructureDispToSolid(), 3);

    //! inverse operators

    setDerFluidLoadToFluid(residualFSI());
    M_BCh_du_inv->addBC("Wall", 1, Natural  , Full,
                        *bcvDerFluidLoadToFluid(), 3);


    setDerFluidLoadToStructure(residualFSI());
    M_BCh_dz_inv->addBC("Interface", 1, Natural  , Full,
                        *bcvDerFluidLoadToStructure(), 3);
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

    FSIOperator::dof_interface_type dofReducedFluidToMesh
        (new FSIOperator::dof_interface_type::element_type);

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

    M_volToSurf.resize(totalDof);


    for ( ID i = 1; i <= BCbEss.list_size(); ++i )
    {
        // Loop on components involved in this boundary condition
//        for ( ID j = 1; j <= nComp; ++j )
//
        {
            // Global Dof
            UInt index = BCbEss( i ) ->id();// + ( BCbEss.component( j ) - 1 ) * totalDof; 
            M_volToSurf[index] = i;
            std::cout << index << " -> " << i << std::endl;
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
    vec = ZeroVector( vec.size() );

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
FSIOperator* createSP(){ return new steklovPoincare(); }
static bool reg = FSIFactory::instance().registerProduct( "steklovPoincare", &createSP );
}


}
