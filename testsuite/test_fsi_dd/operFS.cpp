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
                   BC_Handler &BCh_du,
                   BC_Handler &BCh_dz):
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

// Member functions
    
    
    void operFS::eval(const Vector& disp,
                      int           status,
                      Vector&       dispNew,
                      Vector&       velo)
    {
        if(status) M_nbEval = 0; // new time step
        M_nbEval++;

        std::cout << "dim solid.d = " << M_solid.d().size()
                  << std::endl;
        std::cout << "dim disp    = " << disp.size()
                  << std::endl;
        
        M_solid.d() = disp;

        M_fluid.updateMesh(M_time);
        M_fluid.iterate   (M_time);

        M_solid._recur = 0;

        M_solid.iterate();

        dispNew = M_solid.d();
        velo    = M_solid.w();

        std::cout << "                ::: norm(disp     ) = "
                  << maxnorm(disp) << std::endl;
        std::cout << "                ::: norm(dispNew  ) = "
                  << maxnorm(dispNew) << std::endl;
        std::cout << "                ::: norm(velo     ) = "
                  << maxnorm(velo) << std::endl;
    }

//! Computing the residual on the fluid/structure interface
    void operFS::computeResidualFSI()
    {
        for (UInt ifsi = 0; ifsi < M_residualFSI.size(); ++ifsi)
        {
        }
    }
        
    
//
// Residual evaluation
//
    
    void operFS::evalResidual(Vector &disp,
                              int iter,
                              Vector &res)
    {
        int status = 0;

        if(iter == 0) status = 1;

        std::cout << "*** Residual computation g(x_" << iter <<" )";
        if (status) std::cout << " [NEW TIME STEP] ";
        std::cout << std::endl;

        eval(disp, status, M_dispStruct, M_velo);
        
        res = disp - M_dispStruct;

        M_residualS = M_solid.residual();

        M_residualF = M_fluid.residual();

        //disp = M_dispStruct;
    }


    

    

    
//    
// Preconditionner computation using AZTEC
//
    
    void  operFS::updatePrec(Vector& sol,int iter)
    {
    }


//

    void  operFS::solvePrec(Vector &_uBarS)
    {
        // data_org assigned "by hands": no parallel computation is performed

        UInt dim = _uBarS.size();

        std::cout << "  o-  Solving the preconditionned system... ";
        Chrono chrono;

        for (int ii = 0; ii < (int) dim; ++ii)
            solid().d()[ii] =  _uBarS[ii];

        fluid().updateDispVelo();
        solveLinearFluid();
        solveLinearSolid();

        for (int i = 0; i < (int) dim; ++i)
            _uBarS[i] =  dz()[i] - _uBarS[i];

        chrono.stop();

        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }
 

    void  operFS::solvePrec(Vector& step,
                            const Vector& res,
                            double& linear_rel_tol)
    {
        Real omega(0.1);

        generalizedAitken<Vector,Real> aitk(res.size(),omega);

        UInt precChoice = 1;

        switch(precChoice)
        {
            case 0:
                // Dirichlet-Neumann preconditioner
                invSfPrime(res, linear_rel_tol, step);
                break;
            case 1:
                // Dirichlet-Neumann preconditioner
                invSfPrime(res, linear_rel_tol, step);
                break;
            case 2:
                // Dirichlet-Neumann preconditioner
            //     Vector mu_f(res.size());
//                 Vector mu_s(res.size());

//                 invSfPrime(res, linear_rel_tol, mu_f);
//                 invSsPrime(res, linear_rel_tol, mu_s);

//                 step = aitk.computeDeltaLambda(mu_f,mu_s);
                break;
            default:
                // Newton preconditioner
                invSfSsPrime(res, linear_rel_tol, step);
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
        M_rhs_dz = 0.0;

        if ( !M_BCh_dz.bdUpdateDone() )
            M_BCh_dz.bdUpdate(M_solid._mesh, M_solid._feBd,
                              M_solid._dof);

        bc_manage_vector(M_rhs_dz, M_solid._mesh, M_solid._dof,
                         M_BCh_dz, M_solid._feBd, 1.0, 1.0);

        Real tol       = 1.e-10;

        M_solid._recur = 1;
        M_solid.solveJac(M_dz, M_rhs_dz, tol);
    }


    void  operFS::invSfPrime(const Vector& res,
                             double linear_rel_tol,
                             Vector& step)
    {
        // step = S'_f^{-1} \cdot res

      

    }

    void  operFS::invSsPrime(const Vector& res,
                             double linear_rel_tol,
                             Vector& step)
    {
        // step  = S'_s^{-1} \cdot res


    }

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
        AZ_set_MATFREE(J, &M_dataJacobian, my_matvecJacobian);

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


    void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX *J, int proc_config[])
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
            for (int i=0; i <(int)dim; ++i)
                my_data->M_pFS->solid().d()[i] =  z[i];

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
 
