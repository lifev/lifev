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
    operFS::operFS(NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> >& fluid,
                   VenantKirchhofSolver< RegionMesh3D_ALE<LinearTetra> >& solid,
                   BC_Handler& BCh_du,
                   BC_Handler& BCh_dz):
        M_fluid       (fluid),
        M_solid       (solid),
        M_dispStruct  ( 3*M_solid.dDof().numTotalDof() ),
        M_velo        ( 3*M_solid.dDof().numTotalDof() ),
        M_dz          ( 3*M_solid.dDof().numTotalDof() ),
        M_rhs_dz      ( 3*M_solid.dDof().numTotalDof() ),
        M_nbEval      (0),
        M_BCh_du      (BCh_du),
        M_BCh_dz      (BCh_dz),
        M_dataJacobian(this)
    {
// constructors already called, nothing to do
    }


    void operFS::eval(Vector& dispNew,
                      Vector& velo,
                      const Vector& disp,
                      int status)
    {
        if(status) M_nbEval = 0; // new time step
        M_nbEval++ ;

        M_solid.d() = disp;

        M_fluid.updateMesh(M_time);
        M_fluid.iterate   (M_time);

        M_solid._recur=0;
        
        M_solid.iterate();

        dispNew = M_solid.d();
        velo    = M_solid.w();

        cout << "                ::: norm(disp     ) = "
             << maxnorm(disp) << endl;
        cout << "                ::: norm(dispNew  ) = "
             << maxnorm(dispNew) << endl;
        cout << "                ::: norm(velo     ) = "
             << maxnorm(velo) << endl;

    }


// Residual evaluation
//
    void operFS::evalResidual(Vector& res,
                              const Vector& disp,
                              int iter)
    {
        int status = 0;
        if(iter == 0) status = 1;
        cout << "*** Residual computation g(x_" << iter <<" )";
        if (status) cout << " [NEW TIME STEP] ";
        cout << endl;
        eval(M_dispStruct,M_velo,disp,status);
        res = disp - M_dispStruct;

    }


//
    void  operFS::updateJac(Vector& sol,int iter)
    {
    }


//
    void  operFS::solveJac(Vector& step, const Vector& res, double& linear_rel_tol) {

        // AZTEC specifications for the second system
        int    data_org[AZ_COMM_SIZE];   // data organisation for J
        int    proc_config[AZ_PROC_SIZE];  // Processor information:
        int    options[AZ_OPTIONS_SIZE];   // Array used to select solver options.
        double params[AZ_PARAMS_SIZE];     // User selected solver paramters.
        double status[AZ_STATUS_SIZE];     // Information returned from AZ_solve()

        AZ_set_proc_config(proc_config, AZ_NOT_MPI);

        // data_org assigned "by hands": no parallel computation is performed
        UInt dim_res = res.size();

        data_org[AZ_N_internal] = dim_res;
        data_org[AZ_N_border  ] = 0;
        data_org[AZ_N_external] = 0;
        data_org[AZ_N_neigh   ] = 0;

        // Recovering AZTEC defaults options and params

        AZ_defaults(options,params);

        // Fixed Aztec options for this linear system
        options[AZ_solver  ] = AZ_gmres;
        options[AZ_output  ] = 1;
        options[AZ_poly_ord] = 5;
        options[AZ_kspace  ] = 40;
        options[AZ_conv    ] = AZ_rhs;
        params [AZ_tol     ] = linear_rel_tol;

        //AZTEC matrix for the jacobian

        AZ_MATRIX *J;
        
        J = AZ_matrix_create(dim_res);

        // data containing the matrices C, D, trD and H as pointers
        // are passed through A_ii and pILU_ii:
        AZ_set_MATFREE(J, &M_dataJacobian, my_matvecJacobian);

        cout << "  o-  Solving Jacobian system... ";
        Chrono chrono;

        for (UInt i = 0; i < dim_res; ++i)
            step[i] = 0.0;

        chrono.start();
        
        AZ_iterate(&step[0], &res[0], options, params,
                   status, proc_config, J, NULL, NULL);

        chrono.stop();

        cout << "done in " << chrono.diff() << " s." << endl;


        AZ_matrix_destroy(&J);
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


    UInt operFS::nbEval()
    {
        return M_nbEval;
    }

    void operFS::setTime(const Real& time)
    {
        M_time = time;
    }

    void my_matvecJacobian(double *z, double *Jz, AZ_MATRIX* J, int proc_config[])
    {
        // Extraction of data from J
        DataJacobian* my_data = static_cast< DataJacobian* >(AZ_get_matvec_data(J));

        UInt dim = my_data->_pFS->M_dz.size();

        double xnorm =  AZ_gvector_norm(dim,-1,z,proc_config);
        cout << " ***** norm (z)= " << xnorm << endl<< endl;

        if ( xnorm == 0.0 ) {
            for (int i=0; i <(int)dim; ++i)
                Jz[i] =  0.0;
        }
        else {
            for (int i=0; i <(int)dim; ++i) {
                my_data->_pFS->M_solid.d()[i] =  z[i];
            }
            my_data->_pFS->M_fluid.updateDispVelo();
            my_data->_pFS->solveLinearFluid();
            my_data->_pFS->solveLinearSolid();
            for (int i = 0; i < (int) dim; ++i)
                Jz[i] =  z[i] - my_data->_pFS->M_dz[i];
        }
        cout << " ***** norm (Jz)= " << AZ_gvector_norm(dim,-1,Jz,proc_config)<< endl<< endl;
    }
}
