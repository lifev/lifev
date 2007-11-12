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
//#include <life/lifesolver/reducedLinFluid.hpp>


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

    M_dz.reset    (new vector_type(*this->M_solidInterfaceMap));
    M_rhs_dz.reset(new vector_type(*this->M_solidInterfaceMap));

    if ( this->isFluid() )
    {
        M_rhsNew.reset(new vector_type(this->M_fluid->getMap()));
        M_beta.reset(new vector_type(this->M_fluid->getMap()));
    }
//    M_reducedLinFluid.reset(new reducedLinFluid(this, M_fluid, M_solid));
}
//
// Residual computation
//

void exactJacobian::eval(vector_type&       dispNew,
                         vector_type&       velo,
                         const vector_type &disp,
                         const int          status)
{

    if(status)
        {
            M_nbEval = 0; // new time step
            this->M_fluid->resetPrec();
            this->M_solid->resetPrec();
        }

    M_nbEval++ ;


    MPI_Barrier(MPI_COMM_WORLD);

    // possibly unsafe when using more cpus, since both has repeated maps
    *M_lambdaFluid       = disp;


    *this->M_sigmaFluid *= 0.;

    if (this->isFluid())
    {
        this->M_meshMotion->iterate();

        vector_type const meshDisplacement( M_meshMotion->displacement(), this->M_meshMotion->getRepeatedEpetraMap() );

        this->transferMeshMotionOnFluid(meshDisplacement,
                                        *this->M_veloFluidMesh);

        *this->M_veloFluidMesh    -= *M_dispFluidMeshOld;
        *this->M_veloFluidMesh    *= 1./(M_dataFluid->timestep());

        // copying displacement to a repeated indeces displacement, otherwise the mesh wont know
        // the value of the displacement for some points

        this->moveMesh(meshDisplacement);

        *M_beta *= 0.;

        vector_type const meshDispDiff( M_meshMotion->dispDiff(), this->M_meshMotion->getRepeatedEpetraMap() );

        this->interpolateVelocity(meshDispDiff, *M_beta);

        *M_beta *= -1./M_dataFluid->timestep();

        *M_beta  += *this->M_un;

        double alpha = 1./M_dataFluid->timestep();

        *M_rhsNew   = *this->M_rhs;
        *M_rhsNew  *= alpha;

        this->M_fluid->updateSystem( alpha, *M_beta, *M_rhsNew );

        this->M_fluid->iterate( *M_BCh_u );

        this->transferFluidOnInterface(this->M_fluid->residual(), *this->M_sigmaFluid);

    }

    MPI_Barrier(MPI_COMM_WORLD);

     if ( false && this->isFluid() )
         {
//             *M_velAndPressure = this->M_fluid->solution();
//             M_ensightFluid->postProcess( M_nbEval );
         }



    // possibly unsafe when using more cpus, since both has repeated maps
    *this->M_sigmaSolid = *this->M_sigmaFluid;

    MPI_Barrier(MPI_COMM_WORLD);


    if (this->isSolid())
    {
        this->M_solid->iterate( *M_BCh_d );
        this->transferSolidOnInterface(this->M_solid->disp(), *this->M_lambdaSolid);
        this->transferSolidOnInterface(this->M_solid->vel() , *this->M_lambdaDotSolid);
        this->transferSolidOnInterface(this->M_solid->residual() , *this->M_sigmaSolid);
    }

    // possibly unsafe when using more cpus, since both has repeated maps


//     dispNew = *this->M_lambdaSolid;
//     velo    = *this->M_lambdaDotSolid;

    MPI_Barrier(MPI_COMM_WORLD);


    std::cout << " ::: norm(disp     )    = " << disp.NormInf() << std::endl;
    std::cout << " ::: norm(dispNew  )    = " << dispNew.NormInf() << std::endl;
    std::cout << " ::: norm(velo     )    = " << velo.NormInf() << std::endl;
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

void exactJacobian::evalResidual(vector_type&       res,
                                 const vector_type& disp,
                                 const int     iter)
{
    int status = 0;

    if(iter == 0) status = 1;

    std::cout << "*** Exact Jacobian: Residual computation g(x_" << iter <<")";
    if (status) std::cout << " [NEW TIME STEP] ";
    std::cout << std::endl;

    eval(*M_lambdaSolid, *M_lambdaDotSolid, disp, status);

    res  = *M_lambdaSolid;
    res -= disp;
}


//
// new step computation resolution
//


void  exactJacobian::solveJac(vector_type         &_muk,
                              const vector_type  &_res,
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
//    this->M_fluid->updateLinearSystem(
    this->M_fluidLin->iterate(*M_BCh_du);
}


//


void  exactJacobian::solveLinearSolid()
{
//     M_rhs_dz = Zerovector_type( M_rhs_dz.size() );
//     M_dz     = Zerovector_type( M_dz.size() );

//     if ( !M_BCh_dz->bdUpdateDone() )
//         M_BCh_dz->bdUpdate(this->M_solid->mesh(),
//                            this->M_solid->feBd(),
//                            this->M_solid->dDof());

//     bcManageVector(M_rhs_dz,
//                    this->M_solid->mesh(),
//                    this->M_solid->dDof(),
//                    *M_BCh_dz,
//                    this->M_solid->feBd(),
//                    1., 1.);

//     Real tol       = 1.e-10;
//     std::cout << "rhs_dz norm = " << norm_2(M_rhs_dz) << std::endl;
//     this->M_solid->setRecur(1);
//     this->M_solid->solveJac(M_dz, M_rhs_dz, tol, M_BCh_dz);
//     std::cout << "dz norm     = " << norm_inf(M_dz) << std::endl;
}


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[])
{
    // Extraction of data from J
//     exactJacobian::dataJacobian* my_data = static_cast< exactJacobian::dataJacobian* >(AZ_get_matvec_data(J));

// //     UInt dim = my_data->M_pFS->dz().size();

// //     double xnorm =  AZ_gvector_norm(dim, -1, z, proc_config);
// //     std::cout << " ***** norm (z)= " << xnorm << std::endl << std::endl;

//     exactJacoian::vector_type zSolid(z);

// //     for (int ii = 0; ii < (int) dim; ++ii)
// //         zSolid[ii] = z[ii];

//     if ( xnorm == 0.0 )
//     {
//         for (int i=0; i <(int)dim; ++i)
//             Jz[i] =  0.0;
//     }
//     else
//     {
//         for (int i = 0; i <( int) dim; ++i)
//         {
//             my_data->M_pFS->solid().disp()[i] =  z[i];
//         }
// //         my_data->M_pFS->fluid().updateDispVelo();
//         my_data->M_pFS->solveLinearFluid();
//         my_data->M_pFS->solveLinearSolid();
//     }

// //     for (int i=0; i <(int)dim; ++i)
// //         Jz[i] =  z[i] - my_data->M_pFS->dz()[i];

//     std::cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim, -1, Jz, proc_config)
//               << std::endl << std::endl;
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
