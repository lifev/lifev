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
#define INRIA

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

    GetPot datafile( "data" );
    long int  m=1;
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
        string fname=mesh_dir+datafile( "mesh_file", "cube_6007.m++" );
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
    U = ZeroVector( dim );
    F = ZeroVector( dim );

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
    //cout << chrono.diff() << "s." << endl;
    cout << "NNz= "<< A.Patt()->nNz()<<endl;
    //A.spy("prima");
    A*=0;
    cout << "NNz= "<< A.Patt()->nNz()<<endl;
    A.spy("dopo");
    return 0;
}
