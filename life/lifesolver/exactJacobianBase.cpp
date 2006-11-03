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


#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/reducedLinFluid.hpp>


namespace LifeV
{

Real fzeroEJ(const Real& /*t*/,
             const Real& /*x*/,
             const Real& /*y*/,
             const Real& /*z*/,
             const ID& /*i*/)
{return 0.0;}


exactJacobian::~exactJacobian()
{}

// void
// exactJacobian::setDataFromGetPot( GetPot const& data )
// {
//     // call the super class to setup the data from getpot file if needed
//     super::setDataFromGetPot( data );
// }


void
exactJacobian::setup()
{
    super::setup();

    M_dz.resize(3*M_solid->dDof().numTotalDof());
    M_rhs_dz.resize(3*M_solid->dDof().numTotalDof());
    M_reducedLinFluid.reset(new reducedLinFluid(this, M_fluid, M_solid));
}
//
// Residual computation
//

void exactJacobian::eval(const Vector &_disp,
                         const int     _status)
{
    if(_status) M_nbEval = 0; // new time step
    M_nbEval++ ;

    this->M_solid->disp() = _disp;

    this->M_fluid->updateMesh(time());
    this->M_fluid->iterate   (time());

    this->M_solid->setRecur(0);
    this->M_solid->iterate();

//      this->M_fluid->postProcess();
//      this->M_solid->postProcess();
}

void exactJacobian::evalResidual(Vector &_res,
                                 const Vector &_disp,
                                 const int     _iter)
{
    int status = 0;

    if(_iter == 0) status = 1;

    std::cout << "*** Exact Jacobian: Residual computation g(x_" << _iter <<")";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    eval(_disp, status);

    M_dispStruct = this->M_solid->disp();
    M_velo       = this->M_solid->w();

    std::cout << " ::: norm(disp     ) = " << norm_inf(_disp)  << std::endl;
    std::cout << " ::: norm(dispNew  ) = " << norm_inf(M_dispStruct) << std::endl;
    std::cout << " ::: norm(velo     ) = " << norm_inf(M_velo) << std::endl;

    std::cout << "Max ResidualF   = " << norm_inf(M_fluid->residual())
              << std::endl;
    std::cout << "Max ResidualS   = " << norm_inf(M_solid->residual())
              << std::endl;

    _res = _disp - M_dispStruct;
}


//
// Boundary conditions setup
//


void exactJacobian::setUpBC()
{
    std::cout << "Boundary Conditions setup ... ";

    UInt dim_solid        = this->M_solid->dDof().numTotalDof();
    UInt dim_fluid        = this->M_fluid->uDof().numTotalDof();
    UInt dim_reducedfluid = this->M_fluid->pDof().numTotalDof();

    //========================================================================================
    //  DATA INTERFACING BETWEEN BOTH SOLVERS
    //========================================================================================
    //
    // Passing data from the fluid to the structure: fluid load at the interface
    //
    setFluidLoadToStructure(this->M_fluid->residual());
    //
    // Passing data from structure to the harmonic Extension: motion of the fluid domain
    //
    setStructureDispToHarmonicExtension(this->M_solid->disp());
    //========================================================================================
    //  Interface BOUNDARY CONDITIONS
    //========================================================================================

    BCFunctionBase bcf(fzeroEJ);

    // Boundary conditions for the harmonic extension of the
    // interface solid displacement
    M_BCh_mesh->addBC("Interface", 1, Essential, Full,
                      *bcvStructureDispToHarmonicExtension(), 3);

//     M_BCh_mesh->bdUpdate(this->M_fluid->mesh(),
//                          this->M_fluid->feBd_u(),
//                          this->M_fluid->uDof());

    // Boundary conditions for the solid displacement
    M_BCh_d->addBC("Interface", 1, Natural,   Full,
                   *bcvFluidLoadToStructure(), 3);

    //========================================================================================
    //  Linear operators BOUNDARY CONDITIONS
    //========================================================================================

    if (!reducedFluid())
    {
        // Passing data from the harmonic extension to the fluid
        setDerHarmonicExtensionVelToFluid(this->M_fluid->dwInterpolated());
        // Passing data from fluid to the structure: du -> dz
        setDerFluidLoadToStructure(this->M_fluid->residual());
        // Boundary conditions for du
        M_BCh_du->addBC("Wall",   1,  Essential, Full,
                        *bcvDerHarmonicExtensionVelToFluid(), 3);
        M_BCh_du->addBC("Edges",  20, Essential, Full, bcf,      3);
        M_BCh_du->addBC("InFlow", 2,  Natural,   Full, bcf,      3);
        M_BCh_du->addBC("OutFlow",3,  Natural,   Full, bcf,      3);


        // Boundary conditions for dz
        M_BCh_dz->addBC("Interface", 1, Natural,   Full,
                        *bcvDerFluidLoadToStructure(), 3);
        M_BCh_dz->addBC("Top",       3, Essential, Full, bcf,  3);
        M_BCh_dz->addBC("Base",      2, Essential, Full, bcf,  3);
    }
    else
    {
        setDerStructureAccToReducedFluid(M_reducedLinFluid->dacc(), 2);
        M_BCh_dp->addBC("Wall",        1, Natural,   Scalar, //da_wall);
                        *bcvDerStructureAccToReducedFluid());
        M_BCh_dp->addBC("Wall_Edges", 20, Essential, Scalar, bcf);
        M_BCh_dp->addBC("InFlow",      2, Essential, Scalar, bcf);
        M_BCh_dp->addBC("OutFlow",     3, Essential, Scalar, bcf);

        M_reducedLinFluid->setUpBC(M_BCh_dp);

        setDerReducedFluidLoadToStructure(M_reducedLinFluid->minusdp(), 1);
        M_BCh_dz->addBC("Interface", 1, Natural,   Full, //dg_wall, 3);
                        *bcvDerReducedFluidLoadToStructure(), 2);
        M_BCh_dz->addBC("Top",       3, Essential, Full, bcf,     3);
        M_BCh_dz->addBC("Base",      2, Essential, Full, bcf,     3);
    }
}


//
// new step computation resolution
//


void  exactJacobian::solveJac(Vector         &_muk,
                              const Vector  &_res,
                              const double   _linearRelTol)
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
    AZ_set_MATFREE(J, &M_dataJacobian, my_matvecJacobianEJ);

    std::cout << "  o-  Solving Jacobian system... ";
    Chrono chrono;

    for (UInt i=0;i<dim_res; ++i)
        _muk[i]=0.0;

    chrono.start();
    AZ_iterate(&_muk[0], const_cast<double*>( &_res[0] ), options, params, status, proc_config, J, NULL, NULL);
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

    AZ_matrix_destroy(&J);
}


void  exactJacobian::solveLinearFluid()
{
    this->M_fluid->iterateLin(time(), *M_BCh_du);
}


//


void  exactJacobian::solveLinearSolid()
{
    M_rhs_dz = ZeroVector( M_rhs_dz.size() );
    M_dz     = ZeroVector( M_dz.size() );

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
//    std::cout << "rhs_dz norm = " << norm_2(M_rhs_dz) << std::endl;
    this->M_solid->setRecur(1);
    this->M_solid->solveJac(M_dz, M_rhs_dz, tol, M_BCh_dz);
//    std::cout << "dz norm     = " << norm_inf(M_dz) << std::endl;
}


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[])
{
    // Extraction of data from J
    exactJacobian::dataJacobian* my_data = static_cast< exactJacobian::dataJacobian* >(AZ_get_matvec_data(J));

    UInt dim = my_data->M_pFS->dz().size();

    double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
    std::cout << " ***** norm (z)= " << xnorm << std::endl << std::endl;

    Vector zSolid(dim);

    for (int ii = 0; ii < (int) dim; ++ii)
        zSolid[ii] = z[ii];

    if ( xnorm == 0.0 )
    {
        for (int i=0; i <(int)dim; ++i)
            Jz[i] =  0.0;
    }
    else
    {
        if (!my_data->M_pFS->reducedFluid())
        {
            for (int i=0; i <(int)dim; ++i)
            {
                my_data->M_pFS->solid().disp()[i] =  z[i];
            }
            my_data->M_pFS->fluid().updateDispVelo();
            my_data->M_pFS->solveLinearFluid();
            my_data->M_pFS->solveLinearSolid();
        }
        else
        {
            Vector da(dim);
            double dt   = my_data->M_pFS->fluid().dt();
            double dti2 = 1.0/( dt*dt);

            da = - my_data->M_pFS->fluid().density()*dti2*zSolid;

            my_data->M_pFS->getReducedLinFluid()->setDacc(da);
            my_data->M_pFS->getReducedLinFluid()->solveReducedLinearFluid();
            my_data->M_pFS->solveLinearSolid();
        }
        for (int i=0; i <(int)dim; ++i)
            Jz[i] =  z[i] - my_data->M_pFS->dz()[i];
    }

    std::cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim, -1, Jz, proc_config)
              << std::endl << std::endl;
}


//
// add steklovPoincare to factory
//
namespace
{
FSIOperator* createEJ(){ return new exactJacobian(); }
static bool reg = FSIFactory::instance().registerProduct( "exactJacobian", &createEJ );
}

}
