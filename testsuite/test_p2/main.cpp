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

  Simple Laplacian test with Dirichlet Boundary condition

  Solve the problem

               - \Delta u = f

	                u = g on the boundary


			on a cube
*/
#include <GetPot.hpp>

#include "main.hpp"
#include "ud_functions.hpp"
#include "bcManage.hpp"
#include "elemMat.hpp"
#include "elemOper.hpp"
#include "openDX_wrtrs.hpp"
#include "vtk_wrtrs.hpp"
//#include "sobolevNorms.hpp"

#undef  OPER_TEMPLATE
#define P2
#undef INRIA

int main() {
    using namespace LifeV;
    using namespace std;

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
    const QuadRule& qr    = quadRuleTetra5pt;

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
    //  RegionMesh3D<LinearTetra> aMesh;
    RegionMesh3D<QuadraticTetra> aMesh;
#else
    RegionMesh3D<LinearTetra> aMesh;
#endif

    long int  m=1;
    GetPot datafile( "data" );
    std::string mesh_type = datafile( "mesh_type", "INRIA" );
    if ( mesh_type == "INRIA" )
    {
        string mesh_dir = datafile( "mesh_dir", "." );
        string fname=mesh_dir+datafile( "mesh_file", "cube_6007.mesh" );
        readINRIAMeshFile(aMesh,fname,m);
    }
    else if ( mesh_type == "MESH++" )
    {
        string mesh_dir = datafile( "mesh_dir", "." );
        string fname=mesh_dir+datafile( "mesh_file", "cube_48.m++" );
        //  string fname=mesh_dir+"cube_48.m++";
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
    U=ZeroVector( dim );
    F=ZeroVector( dim );

    // ==========================================
    // Pattern construction and matrix assembling
    // ==========================================
    cout << "dim                    = " << dim     << endl << endl;

    // pattern for stiff operator
    MSRPatt pattA(dof);

    cout << "Values" << endl;

    // A: stiff operator
    MSRMatr<double> A(pattA);

    cout << "*** Matrix computation           : "<<endl;
    chrono.start();
    //
    SourceFct sourceFct;
#ifdef OPER_TEMPLATE
    // assembling of A: stiff operator
    Stiff Ostiff(&fe);
    EOStiff stiff(Ostiff);
    assemble(mu*stiff+sigma*mass,aMesh,fe,dof,sourceFct,A,F);
#else
    ElemMat elmat(fe.nbNode,1,1);
    ElemVec elvec(fe.nbNode,1);
    for(UInt i = 1; i<=aMesh.numVolumes(); i++){
        fe.updateFirstDerivQuadPt(aMesh.volumeList(i));
        elmat.zero();
        elvec.zero();
        stiff(1.,elmat,fe);
        source(sourceFct,elvec,fe,0);
        assemb_mat(A,elmat,fe,dof,0,0);
        assemb_vec(F,elvec,fe,dof,0);
    }
#endif
    chrono.stop();
    cout << "A has been constructed" << endl;
    cout << chrono.diff() << "s." << endl;

    // ====================================
    // Treatment of the Boundary conditions
    // ====================================

    // BC manage for the velocity
    cout << "*** BC Management: "<<endl;

    Real tgv=1.;

    chrono.start();
    bcManage(A,F,aMesh,dof,BCh,feBd,tgv,0.0);

    chrono.stop();
    cout << chrono.diff() << "s." << endl;

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

    AZ_solve(U.giveVec(), F.giveVec(), options, params, NULL,
             (int *)pattA.giveRaw_bindx(), NULL, NULL, NULL,
             A.giveRaw_value(), data_org,
             status, proc_config);
    //
    chrono.stop();
    cout << "*** Solution computed in " << chrono.diff() << "s." << endl;

    // ==========================
    // output writing with openDX
    // ==========================
    /*
      string NameFile= "outputU.dx";
      wr_opendx_header(NameFile,aMesh,dof,fe,"P2");
      wr_opendx_scalar(NameFile,"scalar",U);
    */
    wr_vtk_ascii_header("scal.vtk","Scalar output",aMesh, dof, fe);
    wr_vtk_ascii_scalar("scal.vtk","scalar",U.giveVec(),U.size());

    /* comparison with the analytical solution possible
       in case of using the same example as OFICIAL_TEST/Test_ESSENTIAL/
       taking its source function and boundary conditions.

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
    */
    for (UInt jj=0;jj<dim;jj++)
        cout << U(jj) << " ** ";

    cout << " " << endl;

    return 0;
}
