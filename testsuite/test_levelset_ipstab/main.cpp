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

/*
  \Simple levelset advection problem with IP stabilization
  \After Burman-Hansbo, Comput. Methods Appl. Mech. Engrg. 193 (2004)
  \pp. 1437-1453
*/

#include <GetPot.hpp>
#include <SolverAztec.hpp>
#include <LevelSetSolver.hpp>

#include "main.hpp"
#include "ud_functions.hpp"
#include "bcManage.hpp"
#include "elemMat.hpp"
#include "elemOper.hpp"
#include "openDX_wrtrs.hpp"
#include "vtk_wrtrs.hpp"
#include "bdf.hpp"

#undef INRIA
#define USE_AZTEC_SOLVER

int main() {
    using namespace LifeV;

    Chrono chrono;

    // Boundary conditions definition

    /*
      No boundary condition is needed in this testcase since the normal
      component of velocity is zero all over the boundary and hence no inlet
      boundary exists. This part of the code is however only commented out 
      since it may serve as a template.

      BCFunctionBase gv1(g1); // Functor storing the user definded function g
      BCFunctionBase gv2(g2); // Functor storing the user definded function g
      BCHandler BCh(2); // We impose two boundary conditions

      BCh.addBC("Inlet",  10, Essential, Scalar, gv1);
      BCh.addBC("Outlet",  20, Essential, Scalar, gv2);
    */
    BCHandler BCh;

    // Finite element stuff

    const GeoMap& geoMap = geoLinearTetra;
    const GeoMap& geoMapBd = geoMap.boundaryMap();//geoLinearTria;
    
    const QuadRule& qr = quadRuleTetra5pt;
    const QuadRule& qrBd = quadRuleTria3pt;

    const RefFE& refFE = feTetraP1;
    const RefFE& refBdFE = refFE.boundaryFE();//feTriaP1;

    // Mesh stuff

    typedef RegionMesh3D<LinearTetra> meshType;
    meshType mesh;

    long int  m=1;
    GetPot datafile( "data" );
    std::string mesh_type = datafile( "hyp/discretization/mesh_type", "INRIA" );
    std::cout << mesh_type << " pippo " << std::endl;
    if ( mesh_type == "INRIA" )
    {
        std::string mesh_dir = datafile( "hyp/discretization/mesh_dir", "." );
        std::string fname=mesh_dir+datafile( "hyp/discretization/mesh_file", "cube_6000.mesh" );
        readINRIAMeshFile(mesh,fname,m);
    }
    else if ( mesh_type == ".m++" )
    {
        std::string mesh_dir = datafile( "hyp/discretization/mesh_dir", "." );
        std::string fname=mesh_dir+datafile( "hyp/discretization/mesh_file", "cube_6000.m++" );
        readMppFile(mesh,fname,m);
    }
    else
    {
        std::cerr << "wrong mesh type. It can be either MESH++ or INRIA" << std::endl;
        return EXIT_FAILURE;
    }

    mesh.updateElementFaces();
    mesh.updateElementEdges();
    
    // Build ALL faces (interior faces are necessary for IP stabilization)

    std::cout << " o-> Building face list" << std::endl;

    UInt numIFaces = mesh.numFaces() - mesh.numBFaces();
    UInt numBFaces = mesh.numBFaces();
    buildFaces(mesh, std::cout, std::cerr, numBFaces, numIFaces, true, true, false);

    // Current FE classes for the problem under study with mapping and
    // quadrature rules

    CurrentFE fe(refFE, geoMap, qr);
    CurrentBdFE feBd(refBdFE, geoMapBd, qrBd);

    // Update of the Dof for the particular FE problem and for the boundary
    // conditions
 
    Dof dof(refFE);
    dof.update(mesh);

    /*
      No boundary condition is needed in this testcase since the normal
      component of velocity is zero all over the boundary and hence no inlet
      boundary exists. This part of the code is however only commented out 
      since it may serve as a template.

      BCh.bdUpdate( mesh,  feBd, dof );
    */

    UInt dim = dof.numTotalDof();

    // Velocity field projection

    bool analyticalBeta = datafile("levelset/parameters/analytical_beta", false);

    SolverAztec solverMass;
    solverMass.setOptionsFromGetPot(datafile, "levelset/solver-mass");

    ElemMat elmatM(fe.nbNode, 1, 1);
    PhysVectUnknown<Vector> betaVec(dim);
    betaVec = ZeroVector(NDIM * dim);

    if (!analyticalBeta) {
        std::cout << "O-> Projecting velocity field onto P1 fe space" << std::endl;

        MSRPatt pattM_NDIM(dof, mesh, NDIM);
        MSRMatr<Real> M_NDIM(pattM_NDIM);

        for(int ic = 0; ic < NDIM; ic++)
            for(UInt i = 1; i <= mesh.numVolumes(); i++){
                fe.updateJac(mesh.volumeList(i));

                elmatM.zero();

                mass(1., elmatM, fe, ic, ic);
                assemb_mat(M_NDIM, elmatM, fe, dof, ic, ic);
            }

    

        vortex analyticalVelocityField;
        projectVelocityField(mesh, fe, dof, betaVec, analyticalVelocityField, NDIM, M_NDIM, solverMass);

    }

    // Boundary conditions handling

    /*
      No boundary condition is needed in this testcase since the normal
      component of velocity is zero all over the boundary and hence no inlet
      boundary exists. This part of the code is however only commented out 
      since it may serve as a template.

      std::cout << "O-> BC Management " << std::endl;

      Real tgv=1.;

      bcManage(A,F,mesh,dof,BCh,feBd,tgv,0.0);
    */

    std::cout << "O-> Solving the problem" << std::endl;

    SolverAztec solver;

    LevelSetSolver<meshType> lss(datafile, refFE, qr, qrBd, BCh, fe, dof, betaVec);
    lss.initialize(sphere);

    const LevelSetSolver<meshType>::lsfunction_type& U = lss.lsfunction();
    std::string outputFileRoot = "./results/ls_ip";

    // Save initial conditions

    wr_opendx_header(outputFileRoot + "0000.dx", mesh, dof); 
    wr_opendx_scalar(outputFileRoot + "0000.dx", "levelset_ipstab", U);

    // Retrieve time advancement parameters

    Real t0 = datafile("hyp/bdf/t0", 0.);
    Real delta_t = datafile("hyp/bdf/delta_t", 0.01);
    Real T = datafile("hyp/bdf/T", 0.04);

    UInt save_every = datafile("levelset/parameters/save_every", 1);
    UInt reini_every = datafile("levelset/parameters/reinit_every", 10);

    UInt current_step = 1;
    UInt steps_after_last_reini = 1;
    UInt steps_after_last_save = 1;
    
    for(Real t = t0; t < T; t += delta_t) {
        std::cout << " o-> Step: " << current_step << ", t = " << t << std::endl;
        lss.timeAdvance();
        
        if(steps_after_last_save == save_every) {
            std::ostringstream number;
            number.width(4);
            number.fill('0');
            number << current_step;

            wr_opendx_header(outputFileRoot + number.str() + ".dx", mesh, dof); 
            wr_opendx_scalar(outputFileRoot + number.str() + ".dx", "levelset_ipstab", U);
        } else
            steps_after_last_save++;

        if(steps_after_last_reini == reini_every) {
            std::cout << "  - reinitializing signed distance function" << std::endl;
            lss.directReinitialization();
        }
        else
            steps_after_last_reini++;

        current_step++;
    }
    
    std::cout << "O-> Done" << std::endl;

    return 0;
}
