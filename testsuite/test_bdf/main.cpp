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
/* ========================================================

  Simple Fourier test with Dirichlet Boundary condition

  Solve the problem

           \partial_t u - \nu(t) \Delta u + sigma(t) u = f(t)

	                u = g on the boundary
                        u(t=0) = u0 initial condition

			on a cube
\nu, \sigma and \source can be function of time
(which implies that the matrix needs to be reassembled each time)

 Purpose: Test BDF of different order
*/

#include <GetPot.hpp>


#include "main.hpp"
#include "ud_functions.hpp"
#include "bc_manage.hpp"
#include "elemMat.hpp"
#include "elemOper.hpp"
#include "bdf.hpp"
#include "openDX_wrtrs.hpp"
#include "vtk_wrtrs.hpp"
#include "sobolevNorms.hpp"

#undef  OPER_TEMPLATE
#define P2
#undef INRIA

int main() {
    using namespace LifeV;
    using namespace std;
    try
    {
        Chrono chrono;

        // ===================================================
        // Boundary conditions definition
        // ===================================================

        BCFunctionBase gv1(g1); // Functor storing the user definded function g
        BCFunctionBase gv2(g2); // Functor storing the user definded function g
        BCHandler BCh(2); // We impose two boundary conditions

        BCh.addBC("Inlet",  10, Essential, Scalar, gv1);
        BCh.addBC("Outlet",  20, Essential, Scalar, gv2);
        // Ouput
        BCh.showMe();

        // ===================================================
        // Finite element stuff
        // ===================================================
        const GeoMap& geoMap  = geoLinearTetra;
        const QuadRule& qr    = quadRuleTetra64pt;

        const GeoMap& geoMapBd = geoLinearTria;
        const QuadRule& qrBd   = quadRuleTria3pt;

#ifndef P2
        // P1 elements
        const RefFE& refFE    = feTetraP1;
        const RefFE& refBdFE   = feTriaP1;
#else
        //P2 elements
        const RefFE& refFE    = feTetraP2;
        const RefFE& refBdFE   = feTriaP2;
#endif

        // ===================================================
        // Mesh stuff
        // ===================================================
#ifdef P2
        RegionMesh3D<LinearTetra> aMesh;
        //  RegionMesh3D<QuadraticTetra> aMesh;
#else
        RegionMesh3D<LinearTetra> aMesh;
#endif

        GetPot datafile( "data" );
        long int  m=1;
        std::string mesh_type = datafile( "mesh_type", "INRIA" );
        string mesh_dir = datafile( "mesh_dir", "." );
        string fname=mesh_dir+datafile( "mesh_file", "cube_6007.mesh" );

        if ( mesh_type == "INRIA" )
        {
            readINRIAMeshFile(aMesh,fname,m);
        }
        else if ( mesh_type == "MESH++" )
        {
            readMppFile(aMesh,fname,m);
        }
        else
        {
            std::cerr << "wrong mesh type. It can be either MESH++ or INRIA" << std::endl;
            return EXIT_FAILURE;
        }

        // Avaliable meshes
        //string fname=mesh_dir+"cube_48.m++";
        aMesh.showMe();

        cout<< "Now building local Edges/faces Stuff"<<endl<<endl;
        aMesh.updateElementEdges();
        aMesh.updateElementFaces();
        aMesh.showMe();

        // ===================================================
        // Current FE classes for the problem under study with
        // mapping and quadrature rules
        // ===================================================

        CurrentFE fe(refFE,geoMap,qr);
        CurrentBdFE feBd(refBdFE,geoMapBd,qrBd);

        // ===============================================
        // Update of the Dof for the particular FE problem
        // and for the boundary conditions
        // ===============================================

        Dof dof(refFE);
        dof.update(aMesh);

        BCh.bdUpdate( aMesh,  feBd, dof );

        UInt dim = dof.numTotalDof();

        dof.showMe();

        // initialization of vector of unknowns and rhs
        ScalUnknown<Vector> U(dim), F(dim);
        U=0.0;
        F=0.0;

        // ==========================================
        // Definition of the time integration stuff
        // ==========================================
        Real Tfin = datafile( "bdf/endtime", 10.0 );
        Real delta_t = datafile( "bdf/timestep", 0.5 );
        Real t0 = 0.;
        UInt ord_bdf = datafile( "bdf/order", 3 );;
        Bdf bdf(ord_bdf);
        Real coeff=bdf.coeff_der(0)/delta_t;

        bdf.showMe();

        // ==========================================
        // Pattern construction and time independent matrix assembling
        // ==========================================
        cout << "dim                    = " << dim     << endl << endl;

        // pattern for stiff operator
        MSRPatt pattA(dof);

        cout << "Values" << endl;

        // A: Fourier operator = alpha/dt Mass + stiff
        MSRMatr<double> A(pattA);
        // M : mass matrix
        MSRMatr<double> M(pattA);
        M.zeros();

        cout << "*** Matrix computation           : "<<endl;
        chrono.start();
        //
        SourceFct sourceFct;
#ifdef OPER_TEMPLATE
        // assembling of A: stiff operator
        Stiff Ostiff(&fe);
        EOStiff stiff(Ostiff);
        Mass Omass(&fe);
        EOMass mass(Omass);
        assemble_symm(mass,aMesh,fe,dof,sourceFct,M,F,t0); //mass matrix
        //  assemble(coeff*mass+stiff,aMesh,fe,dof,sourceFct,A,F);
#else
        ElemMat elmat(fe.nbNode,1,1);
        ElemVec elvec(fe.nbNode,1);
        for(UInt i = 1; i<=aMesh.numVolumes(); i++){
            fe.updateJacQuadPt(aMesh.volumeList(i));
            elmat.zero();
            mass(1.,elmat,fe);
            assemb_mat(M,elmat,fe,dof,0,0);
        }
#endif

        // ==============================
        // Resolution of the linear system
        // ==============================
        int    proc_config[AZ_PROC_SIZE];// Processor information:
        //  proc_config[AZ_node] = node name
        //  proc_config[AZ_N_procs] = # of nodes
        int    options[AZ_OPTIONS_SIZE]; // Array used to select solver options.
        double params[AZ_PARAMS_SIZE];   // User selected solver paramters.
        int    *data_org;                // Array to specify data layout
        double status[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
        // indicating success or failure.
        // altre dichiarazioni per AZTEC
        int    *update,                  // vector elements updated on this node.
        *external;                // vector elements needed by this node.
        int    *update_index;            // ordering of update[] and external[]
        int    *extern_index;            // locally on this processor.
        //  int    *bindx;                 // Sparse matrix to be solved is stored
        //  double *val;                   // in these MSR arrays.
        int    N_update;                 // # of unknowns updated on this node
        //
        cout << "*** Linear System Solving (AZTEC)" << endl;
        AZ_set_proc_config(proc_config, AZ_NOT_MPI );
        //   cout << AZ_PROC_SIZE << " " << AZ_node << " " << AZ_N_procs << endl;
        //   for (UInt ii=0; ii<AZ_PROC_SIZE; ++ii)
        //     cout << proc_config[ii] << endl;
        AZ_read_update(&N_update, &update, proc_config, U.size(), 1, AZ_linear);

        AZ_defaults(options,params);

        AZ_transform(proc_config, &external,
                     (int *)pattA.giveRaw_bindx(), A.giveRaw_value(),
                     update, &update_index,
                     &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
                     AZ_MSR_MATRIX);

        chrono.start();
        init_options(options,params);

        // ==========================
        // output writing with openDX
        // ==========================

        string NameFile= "outputU.dx";
        wr_opendx_header(NameFile,aMesh,dof,fe,"P1");


        // =====================================
        // TIME LOOP
        // =====================================

        bdf.initialize_unk(u0,aMesh,refFE,fe,dof,t0,delta_t,1);
        bdf.showMe();

        for (Real t=t0+delta_t;t<=Tfin;t+=delta_t)
        {
            cout << "Now we are at time " << t << endl;

            A.zeros();
            F=0;
            // ======================================================================
            // Update of the right hand sied with the solution at the previous steps
            // ======================================================================

            Real visc=nu(t);
            Real s=sigma(t);
#ifdef OPER_TEMPLATE
            // assembling of A
            assemble_symm((coeff+s)*mass+visc*stiff,aMesh,fe,dof,sourceFct,A,F,t);
#else
            for(UInt i = 1; i<=aMesh.numVolumes(); i++){
                fe.updateFirstDerivQuadPt(aMesh.volumeList(i));
                elmat.zero();
                elvec.zero();
                mass(coeff+s,elmat,fe);
                stiff(visc,elmat,fe);
                source(sourceFct,elvec,fe,t,0);
                assemb_mat(A,elmat,fe,dof,0,0);
                assemb_vec(F,elvec,fe,dof,0);
            }
#endif
            chrono.stop();
            cout << "A has been constructed" << endl;
            cout << chrono.diff() << "s." << endl;

            // Handling of the right hand side
            F += M*bdf.time_der(delta_t);

            // ====================================
            // Treatment of the Boundary conditions
            // ====================================
            cout << "*** BC Management: "<<endl;

            Real tgv=1.;

            chrono.start();
            bc_manage(A,F,aMesh,dof,BCh,feBd,tgv,t);

            chrono.stop();
            cout << chrono.diff() << "s." << endl;
            chrono.start();
            AZ_solve(U.giveVec(), F.giveVec(), options, params, NULL,
                     (int *)pattA.giveRaw_bindx(), NULL, NULL, NULL,
                     A.giveRaw_value(), data_org,
                     status, proc_config);
            //
            chrono.stop();
            cout << "*** Solution computed in " << chrono.diff() << "s." << endl;

            wr_opendx_scalar(NameFile,"scalar",U);

            //  wr_vtk_ascii_header("scal.vtk","Scalar output",aMesh, dof, fe);
            //wr_vtk_ascii_scalar("scal.vtk","scalar",U.giveVec(),U.size());

            /*comparison with the analytical solution possible
              in case of using the same example as OFICIAL_TEST/Test_ESSENTIAL/
              taking its source function and boundary conditions.
            */
            // =======================================
            // comparison with the analytical solution
            // =======================================

            AnalyticalSol analyticSol;

            Real normL2=0., normL2diff=0., normL2sol=0.;
            Real normH1=0., normH1diff=0., normH1sol=0.;

            for(UInt i=1; i<=aMesh.numVolumes(); ++i){
                //
                fe.updateFirstDeriv(aMesh.volumeList(i));

                normL2     += elem_L2_2(U,fe,dof);
                normL2sol  += elem_L2_2(analyticSol,fe,t,U.nbcomp());// U.nbcomp()=1 for a scalar problem
                normL2diff += elem_L2_diff_2(U,analyticSol,fe,dof,t,U.nbcomp());

                normH1     += elem_H1_2(U,fe,dof);
                normH1sol  += elem_H1_2(analyticSol,fe,t,U.nbcomp());
                normH1diff += elem_H1_diff_2(U,analyticSol,fe,dof,t,U.nbcomp());
            }

            normL2     = sqrt(normL2);
            normL2sol  = sqrt(normL2sol);
            normL2diff = sqrt(normL2diff);

            normH1     = sqrt(normH1);
            normH1sol  = sqrt(normH1sol);
            normH1diff = sqrt(normH1diff);

            cout << "|| U       ||_{L^2}                   = " << normL2 << endl;
            cout << "|| sol     ||_{L^2}                   = " << normL2sol << endl;
            cout << "|| U - sol ||_{L^2}                   = " << normL2diff<< endl;
            cout << "|| U - sol ||_{L^2} / || sol ||_{L^2} = " << normL2diff/normL2sol
                 << endl;

            cout << "|| U       ||_{H^1}                   = " << normH1 << endl;
            cout << "|| sol     ||_{H^1}                   = " << normH1sol << endl;
            cout << "|| U - sol ||_{H^1}                   = " << normH1diff<< endl;
            cout << "|| U - sol ||_{H^1} / || sol ||_{H^1} = " << normH1diff/normH1sol
                 << endl;



            bdf.shift_right(U);


        } // END OF TIME LOOP
    }
    catch( std::exception const& __e )
    {
        std::cerr << "An exception was caught in LifeV\n"
        << "reason: \n" << __e.what() << "\n";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
