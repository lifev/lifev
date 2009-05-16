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
/* ========================================================

  Simple Laplacian test with Natural Boundary condition

  Solve the problem

               - \Delta u = f in [0,1]^3

                    du/dn = h on {z = 1}

                u = g on the other faces


   Data:

      f(x,y,z)    = -6
      h(x,y,z)    = 2 * z
      g(x,y,z)    = x^2+y^2+z^2

   the exact solution is u(x,y,z) = x^2+y^2+z^2

   ========================================================
*/

#include <cassert>

#include <life/lifecore/GetPot.hpp>

#include "main.hpp"
#include "ud_functions.hpp"
#include <life/lifefem/bcManage.hpp>



int main()
{
    using namespace LifeV;
    Chrono chrono;


    // ===================================================
    // Boundary conditions definition
    // ===================================================


    BCFunctionBase gv(g); // Functor storing the user definded function g

    BCFunctionBase hv(h); // Functor storing the user definded function h

    BCHandler BCh(2); // We impose two boundary conditions


    BCh.addBC("Inlet",  10, Essential, Scalar, gv);
    BCh.addBC("Outlet", 20, Natural, Scalar, hv);


    // Ouput
    BCh.showMe();




    // ===================================================
    // Finite element staff
    // ===================================================
    const GeoMap& geoMap  = geoLinearTetra;
    const QuadRule& qr    = quadRuleTetra4pt;

    const GeoMap& geoMapBd = geoMap.boundaryMap();
    const QuadRule& qrBd   = quadRuleTria3pt;

    // P1 elements
    //const RefFE& refFE    = feTetraP1;
    //const RefFE& refBdFE   = feTriaP1;

    // P2 elements
    const RefFE& refFE    = feTetraP2;
    const RefFE& refBdFE   = refFE.boundaryFE();


    // ===================================================
    // Mesh staff
    // ===================================================
    RegionMesh3D<LinearTetra> aMesh;

    GetPot datafile( "data" );
    std::string mesh_dir = datafile( "mesh_dir", "." );//../data/mesh/mesh++/";
    std::string fname=mesh_dir+datafile( "mesh_file", "cube_48_neumann.m++" );

    long int  m=1;
    readMppFile(aMesh,fname,m);
    aMesh.showMe();

    Debug( 10010 )<< "Now building local Edges/faces Stuff"<<"\n";
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
    U = ZeroVector( dim );
    F = ZeroVector( dim );

    // ==========================================
    // Pattern construction and matrix assembling
    // ==========================================
    Debug( 10010 ) << "dim                    = " << dim     << "\n" << "\n";

    // pattern for stiff operator
    MSRPatt pattA(dof);

    Debug( 10010 ) << "Values" << "\n";

    // A: stiff operator
    MSRMatr<double> A(pattA);

    Debug( 10010 ) << "*** Matrix computation           : "<<"\n";
    chrono.start();
    //
    Stiff Ostiff(&fe);
    EOStiff stiff(Ostiff);

    SourceFct source;

    // assembling of A: stiff operator
    assemble(stiff,aMesh,fe,dof,source,A,F);
    Debug( 10010 ) << "A has been constructed" << "\n";

    chrono.stop();
    Debug( 10010 ) << chrono.diff() << "s." << "\n";

    // ====================================
    // Treatment of the Boundary conditions
    // ====================================

    // BC manage for the velocity
    Debug( 10010 ) << "*** BC Management: "<<"\n";

    Real tgv=1.;

    chrono.start();
    bcManage(A,F,aMesh,dof,BCh,feBd,tgv,0.0);

    chrono.stop();
    //Debug( 10010 ) << chrono.diff() << "s." << "\n";

    // ==============================
    // Reolution of the linear system
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
    Debug( 10010 ) << "*** Linear System Solving (AZTEC)" << "\n";
    AZ_set_proc_config(proc_config, AZ_NOT_MPI );
    //   Debug( 10010 ) << AZ_PROC_SIZE << " " << AZ_node << " " << AZ_N_procs << "\n";
    //   for (UInt ii=0; ii<AZ_PROC_SIZE; ++ii)
    //     Debug( 10010 ) << proc_config[ii] << "\n";

    AZ_read_update(&N_update, &update, proc_config, U.size(), 1, AZ_linear);

    AZ_defaults(options,params);

    AZ_transform(proc_config, &external,
                 (int *)pattA.giveRaw_bindx(), A.giveRaw_value(),
                 update, &update_index,
                 &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
                 AZ_MSR_MATRIX);

    chrono.start();
    init_options(options,params);
    options[AZ_output]=AZ_warnings;

    AZ_solve(U.giveVec(), F.giveVec(), options, params, NULL,
             (int *)pattA.giveRaw_bindx(), NULL, NULL, NULL,
             A.giveRaw_value(), data_org,
             status, proc_config);
    //
    chrono.stop();
    Debug( 10010 ) << "*** Solution computed in " << chrono.diff() << "s." << "\n";

    //

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
        normL2sol  += elem_L2_2(analyticSol,fe);
        normL2diff += elem_L2_diff_2(U,analyticSol,fe,dof);

        normH1     += elem_H1_2(U,fe,dof);
        normH1sol  += elem_H1_2(analyticSol,fe);
        normH1diff += elem_H1_diff_2(U,analyticSol,fe,dof);
    }

    normL2     = sqrt(normL2);
    normL2sol  = sqrt(normL2sol);
    normL2diff = sqrt(normL2diff);

    normH1     = sqrt(normH1);
    normH1sol  = sqrt(normH1sol);
    normH1diff = sqrt(normH1diff);

    std::cout << "|| U       ||_{L^2}                   = " << normL2 << std::endl;
    std::cout << "|| sol     ||_{L^2}                   = " << normL2sol << std::endl;
    std::cout << "|| U - sol ||_{L^2}                   = " << normL2diff<< std::endl;
    std::cout << "|| U - sol ||_{L^2} / || sol ||_{L^2} = " << normL2diff/normL2sol
         << std::endl;

    std::cout << "|| U       ||_{H^1}                   = " << normH1 << std::endl;
    std::cout << "|| sol     ||_{H^1}                   = " << normH1sol << std::endl;
    std::cout << "|| U - sol ||_{H^1}                   = " << normH1diff<< std::endl;
    std::cout << "|| U - sol ||_{H^1} / || sol ||_{H^1} = " << normH1diff/normH1sol
         << std::endl;

    return 0;
}
