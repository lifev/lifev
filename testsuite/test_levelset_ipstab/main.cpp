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

#include <main.hpp>

#include <importer.hpp>

#include <bcManage.hpp>

#include <elemMat.hpp>
#include <elemOper.hpp>

#include <openDX_wrtrs.hpp>
#include <vtk_wrtrs.hpp>
#include <bdf.hpp>

int main() {
    using namespace LifeV;

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
    const GeoMap& geoMapBd = geoLinearTria;

    GetPot datafile( "data" );

    std::string __fe_name = datafile("levelset/discretization/FEM", "P1");
    const RefFE& refFE = __fe_name == "P1" ? feTetraP1 : feTetraP2;
    const RefFE& refBdFE = __fe_name == "P1" ? feTriaP1 : feTriaP2;

    const QuadRule& qr = __fe_name == "P1" ? quadRuleTetra5pt : quadRuleTetra15pt;
    const QuadRule& qrBd = __fe_name == "P1" ? quadRuleTria3pt : quadRuleTria7pt;


    // Mesh stuff

    std::cout << "** LS test ** Importing mesh" << std::endl;

    typedef RegionMesh3D<LinearTetra> meshType;
    meshType __mesh;

    std::string __mesh_file = datafile("levelset/discretization/mesh_dir", "../data/mesh++/");
    __mesh_file += datafile("levelset/discretization/mesh_file", "cube_48.m++");

    importer __importer(__mesh_file, MESHPP);
    __importer.import(__mesh, 1);

    std::cout << "** LS test ** Regenerating edge and face information" << std::endl;

    __mesh.updateElementFaces();
    __mesh.updateElementEdges();

    // Build ALL faces (interior faces are necessary for IP stabilization)

    std::cout << "** LS test ** Building face list" << std::endl;

    UInt numIFaces = __mesh.numFaces() - __mesh.numBFaces();
    UInt numBFaces = __mesh.numBFaces();
    buildFaces(__mesh, std::cout, std::cerr, numBFaces, numIFaces,
               true, true, false);

    // Current FE classes for the problem under study with mapping and
    // quadrature rules

    std::cout << "** LS test ** Defining current finite element" << std::endl;

    CurrentFE __fe(refFE, geoMap, qr);
    CurrentBdFE feBd(refBdFE, geoMapBd, qrBd);
    
    // Update of the Dof for the particular FE problem and for the boundary
    // conditions

    std::cout << "** LS test ** Updating DOF table" << std::endl;

    Dof __dof(refFE);
    __dof.update(__mesh);

    // No boundary condition is needed in this test case since the normal
    // component of velocity is zero all over the boundary and hence no inlet
    // boundary exists. This part of the code is however only commented out
    // since it may serve as a template.

    //BCh.bdUpdate( __mesh,  feBd, __dof );

    UInt dim = __dof.numTotalDof();

    // Velocity field projection

    bool analyticalBeta = datafile("levelset/parameters/analytical_beta",
                                   false);

    PhysVectUnknown<Vector> betaVec(dim);
    betaVec = ZeroVector(NDIM * dim);

    if (!analyticalBeta) {
        std::cout << "** LS test ** Projecting velocity field onto fe space"
                  << std::endl;

        SolverAztec solverMass;
        solverMass.setOptionsFromGetPot(datafile, "levelset/solver-mass");

        ElemMat elmatM(__fe.nbNode, NDIM, NDIM);

        MSRPatt pattM_NDIM(__dof, __mesh, NDIM);
        MSRMatr<Real> M_NDIM(pattM_NDIM);

        std::cout << "** LS test ** Computing projection RHS" << std::endl;
        for(int ic = 0; ic < NDIM; ic++)
            for(UInt i = 1; i <= __mesh.numVolumes(); i++){
                __fe.updateJac( __mesh.volumeList(i) );

                elmatM.zero();

                mass(1., elmatM, __fe, ic, ic);

                assemb_mat(M_NDIM, elmatM, __fe, __dof, ic, ic);
            }

        std::cout << "** LS test ** Solving mass matrix system" << std::endl;
        vortex analyticalVelocityField;
        projectVelocityField(__mesh, __fe, __dof, betaVec, analyticalVelocityField,
                             NDIM, M_NDIM, solverMass);

    }

    // Export velocity field projection to OpenDX format

    wr_opendx_header("results/beta.dx", __mesh, __dof, __fe, __fe_name);
    wr_opendx_vector("results/beta.dx", "u", betaVec, 3);

    // Boundary conditions handling

    // Real tgv=1.;
    // bcManage(A,F,__mesh,__dof,BCh,feBd,tgv,0.0);

    std::cout << "** LS test ** Solving the problem" << std::endl;

    // Retrieve time advancement parameters

    Real t0 = datafile("levelset/bdf/t0", 0.);
    Real delta_t = datafile("levelset/bdf/delta_t", 0.01);
    Real T = datafile("levelset/bdf/T", 0.04);

    // Solver initialization

    LevelSetSolver<meshType> lss(__mesh, datafile, "levelset", refFE,
                                 qr, qrBd, BCh, __fe, __dof, betaVec);
    lss.initialize(sphere, t0, delta_t);
    lss.setVerboseMode();

    const LevelSetSolver<meshType>::lsfunction_type& U = lss.lsfunction();
    std::string outputFileRoot = "./results/ls_ip";

    // Save initial conditions

    wr_opendx_header(outputFileRoot + "0000.dx", __mesh, __dof, lss.fe(), __fe_name );
    wr_opendx_scalar(outputFileRoot + "0000.dx", "levelset_ipstab", U);

    // Reinitialize and save re-initialized IC

    lss.directReinitialization();
    wr_opendx_header(outputFileRoot + "0000r.dx", __mesh, __dof, lss.fe(), __fe_name );
    wr_opendx_scalar(outputFileRoot + "0000r.dx", "levelset_ipstab", U);

    // Retrieve parameters from data file

    UInt save_every = datafile("levelset/parameters/save_every", 1);
    UInt reini_every = datafile("levelset/parameters/reinit_every", 10);

    // Initialize ostream to save mass information

    std::ofstream __ofile("results/mass.txt");

    // Initialize counters

    UInt current_step = 1;
    UInt steps_after_last_reini = 1;
    UInt steps_after_last_save = 1;
    
    for(Real t = t0; t < T; t += delta_t) {
        std::cout << "** LS test ** Step: " << current_step << ", t = " << t
                  << std::endl;
        lss.timeAdvance();
        __ofile << t << "\t" << lss.computeMass( LevelSetSolver<meshType>::fluid1 ) << std::endl;
        std::cout << lss;
        lss.spy("./results/");

        std::ostringstream number;
        number.width(4);
        number.fill('0');
        number << current_step;

        if(steps_after_last_save == save_every) {
            wr_opendx_header(outputFileRoot + number.str() + ".dx", __mesh, __dof, lss.fe(), __fe_name );
            wr_opendx_scalar(outputFileRoot + number.str() + ".dx",
                             "levelset_ipstab", U);

            steps_after_last_save = 1;
        } else
            steps_after_last_save++;

        if(steps_after_last_reini == reini_every) {
            std::cout << "** LS test ** Reinitializing signed distance function"
                      << std::endl;
            lss.directReinitialization();

            steps_after_last_reini = 1;

            wr_opendx_header(outputFileRoot + number.str() + "r.dx",
                             __mesh, __dof, lss.fe(), __fe_name );
            wr_opendx_scalar(outputFileRoot + number.str() + "r.dx",
                             "levelset_ipstab", U);
        }
        else
            steps_after_last_reini++;

        current_step++;
    }
    
    __ofile.close();

    std::cout << "** LS test ** Done" << std::endl;
    
    return 0;
}
