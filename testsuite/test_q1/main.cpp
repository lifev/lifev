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
	       with hexa

*/

#include <GetPot.hpp>


#include "main.hpp"
#include "ud_functions.hpp"
#include "bc_manage.hpp"
#include "vtk_wrtrs.hpp"
#include "elemMat.hpp"
#include "elemOper.hpp"

#undef OPER_TEMPLATE

int main() {
    using namespace LifeV;
    using namespace std;

    Chrono chrono;

    // ===================================================
    // Boundary conditions definition
    // ===================================================

    BCFunctionBase gv1(g1); // Functor storing the user definded function g1
    BCFunctionBase gv2(g2); // Functor storing the user definded function g2
    BCFunctionBase gv3(g3); // Functor storing the user definded function g3

    BCHandler BCh(5); // We impose five boundary conditions

    BCh.addBC("Wall",  2, Essential, Scalar, gv2);
    BCh.addBC("BInlet",  4, Essential, Scalar, gv1);
    BCh.addBC("Inlet",  1, Essential, Scalar, gv1);
    BCh.addBC("BOutlet",  5, Essential, Scalar, gv3);
    BCh.addBC("outlet",  3, Essential, Scalar, gv3);

    // Ouput
    BCh.showMe();
    // ===================================================
    // Finite element stuff
    // ===================================================
    const GeoMap& geoMap  = geoBilinearHexa;
    const QuadRule& qr    = quadRuleHexa8pt;

    const GeoMap& geoMapBd = geoBilinearQuad;
    const QuadRule& qrBd   = quadRuleQuad4pt;

    // Q1 elements
    const RefFE& refFE    = feHexaQ1;
    const RefFE& refBdFE   = feQuadQ1;

    // ===================================================
    // Mesh stuff
    // ===================================================
    RegionMesh3D<LinearHexa> aMesh;
    long int  m=1;

    GetPot datafile( "data" );
    string dirname=datafile( "mesh_dir","." );//"../data/mesh/inria/";
    string fname=dirname+datafile( "mesh_file","cylhexa.mesh" );//dirname+"cylhexa.mesh";
    readINRIAMeshFile(aMesh,fname,m);


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
    assemble(stiff,aMesh,fe,dof,sourceFct,A,F);
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
    bc_manage(A,F,aMesh,dof,BCh,feBd,tgv,0.0);

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

    //


    //  cout << "The approximation at the nodes ( node, u(node) ): " << endl;
    //  for (UInt i=0; i < U.size() ; ++i)
    //    cout << i+1 << " " << U[i] << endl;
    wr_vtk_ascii_header("post.vtk","Title",aMesh, dof, fe);
    wr_vtk_ascii_scalar("post.vtk","scal",U.giveVec(),U.size());

    return 0;
}
