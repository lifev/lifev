/* -*- mode: c++ -*-
   This program is part of the LifeV library
   
   Author(s): Daniele A. Di Pietro <dipietro@unibg.it>
   Date: 2004-10
   
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

#include <GetPot.hpp>
#include "main.hpp"
#include "SolverAztec.hpp"

int main(){
  using namespace LifeV;
  using namespace std;

  bool built;

  Vortex velocity;
  
  //============================================================================
  // Plotting stuff
  //============================================================================
  
  // Define plot nodes
  KNM<Real> PN(4, 3);
  
  for(UInt i = 1; i < 4; i++)
    for(UInt j = 0; j < 3; j++)
      PN(i, j) = (j==i - 1) ? 1. : 0.;
  PN(0, 0) = 0.; PN(0, 1) = 0.; PN(0, 2) = 0.;
  cout << "** Plot nodes' coordinates on the reference element" << endl;
  for(UInt i = 0; i < 4; i++)
  {
    for(UInt j = 0; j < 3; j++)
    {
      cout << PN(i, j) << "\t";
    }
    cout << endl;
  }
  
  KNM<UInt> Conn(1, 4);
  Conn(0, 0) = 0; Conn(0, 1) = 1; Conn(0, 2) = 2; Conn(0, 3) = 3;
  
  //============================================================================
  // Boundary conditions
  //============================================================================

  BCFunctionBase gv1(g1);
  BCHandler BCh(1);
   
  BCh.addBC("Inlet", 10, Essential, Scalar, gv1);

  //============================================================================
  // Finite element stuff
  //============================================================================

  const GeoMap & geoMap      = geoLinearTetra;
  const QuadRule & qr        = quadRuleTetra64pt;

  const GeoMap & geoMapBd    = geoLinearTria;
  const QuadRule & qrBd      = quadRuleTria4pt;

  const GeoMapDG & geoMapDG  = geoLinearTetraDG;

  const RefFEDG & refFEDG    = feDGTetraP1;

  const RefFE & refBdFE      = feTriaP1;

  //============================================================================
  // Mesh stuff
  //============================================================================

  RegionMesh3D<LinearTetra> mesh;

  long int  m=1;
  GetPot datafile( "data" );
  std::string mesh_type = datafile( "mesh/mesh_type", "INRIA" );
  if ( mesh_type == "INRIA" )
    {
      string mesh_dir = datafile( "mesh/mesh_dir", "." );
      string fname = mesh_dir + datafile( "mesh/mesh_file", "cube_6007.mesh" );
      readINRIAMeshFile(mesh,fname,m);
      cout << "Successfully opened mesh file " << fname << endl;
    }
  else if ( mesh_type == "MESH++" )
    {
      string mesh_dir = datafile( "mesh/mesh_dir", "." );
      string fname = mesh_dir + datafile( "mesh/mesh_file", "cube_48.m++" );      
      readMppFile(mesh,fname,m);
      cout << "Successfully opened mesh file " << fname << endl;
    }
  else
    {
      std::cerr << "wrong mesh type. It can be either MESH++ or INRIA" << std::endl;
      return EXIT_FAILURE;
    }
  
  UInt numIFaces = mesh.numFaces() - mesh.numBFaces();
  UInt numBFaces = mesh.numBFaces();

  built =  buildFaces(mesh, cout, cerr, numBFaces, numIFaces, true, true, false);

  //============================================================================
  // Current discontinuous FE, IF, BF classes for the problem
  //============================================================================

  CurrentFEDG feDG(refFEDG, geoMap, qr);
  cout << "** CurrentFEDG created " << endl;

  CurrentIFDG ifDG(refFEDG, refBdFE, geoMapBd, geoMapDG, qrBd);
  cout << "** CurrentIFDG created " << endl;

  CurrentBFDG bfDG(refFEDG, refBdFE, geoMapBd, geoMapDG, qrBd);
  cout << "** CurrentBFDG created " << endl;

  //============================================================================
  // Update of the Dof and DofByFace
  //============================================================================

  DofDG dof(elPattern_P1_DG_3D);
  dof.update(mesh);
  cout << "** DOF created and updated" << endl;

  DofByFace dofByFace(facePattern_P1_DG_3D);
  dofByFace.update(mesh, dof);
  
  cout << "** DOF by face created and updated" << endl;

  //============================================================================
  // Initialization of vector unknowns and rhs
  //============================================================================

  UInt dim = dof.numTotalDof();

  SourceFct sourceFct;

  ScalUnknown<Vector> U(dim), F(dim), rhs(dim), K1(dim), K2(dim);
  U = 0.0;
  F = 0.0;

  cout << "** Vector unknown and rhs initializated" << endl;

  //============================================================================
  // Pattern construction and matrix assembling
  //============================================================================

  MSRPatt pattA(dof, dofByFace, "dg");
  MSRPatt pattM(dof);

  MSRMatr<double> A(pattA);
  MSRMatr<double> M(pattM);

  cout << "** Pattern for the problem matrices created" << endl;

  //============================================================================
  // Differential operators we want to solve for, wrapped in a suitable class
  //============================================================================

  AdvecDG<Vortex> OadvecDG(&feDG);

  EOAdvecDG advecDG(OadvecDG);

  AdvecIFUW1DG<Vortex> OadvecIFUW1DG(&ifDG);
  AdvecIFUW2DG<Vortex> OadvecIFUW2DG(&ifDG);

  EOAdvecIFUW1DG advecIFUW1DG(OadvecIFUW1DG);
  EOAdvecIFUW2DG advecIFUW2DG(OadvecIFUW2DG);

  AdvecBFUWDG<Vortex> OadvecBFUWDG(&bfDG);
  EOAdvecBFUWDG advecBFUWDG(OadvecBFUWDG);

  //============================================================================
  // Assembling problem matrix
  //============================================================================

  assemble_AdvecDG(advecDG, advecIFUW1DG + advecIFUW2DG, advecBFUWDG, mesh,
  		   BCh, velocity, feDG, ifDG, bfDG, dof, dofByFace, sourceFct, A, M, F);

		   
  cout << "** Finished to assemble matrices A and M" << endl;
  
  A.spy( "stiffness.m" );
  M.spy( "mass.m" );
  
  //============================================================================
  // Initial conditions
  //============================================================================
  ID globalDof;
  Real xDof, yDof, zDof;
  Sphere IC;

  for(UInt j = 1; j <= dof.numElements(); j++){
    for(UInt i = 1; i <= dof.numLocalDof(); i++){
      globalDof = dof.localToGlobal(j, i);

      xDof = mesh.volumeList(j).point(i).x();
      yDof = mesh.volumeList(j).point(i).y();
      zDof = mesh.volumeList(j).point(i).z();

      U(globalDof - 1) = IC(xDof, yDof, zDof);
    }
  }
  //==LINESTOREMOVE
  std::ofstream ofile( "u0.txt", std::ios::app );
  for(UInt ig = 0; ig < dof.numTotalDof(); ig++)
    ofile << U(ig) << endl;
  ofile.close();
  //==ENDOFLINESTOREMOVE
  elem_wise_wrtr("sphere_initial.dx", mesh, dof, feDG, PN, Conn, U);
  
  //============================================================================
  // Advancing in time with a second order RK method
  //============================================================================  
  
  Real h = datafile("time_advancement/h", 0.005);
  Real t0 = datafile("time_advancement/t0", 0.);
  Real T = datafile("time_advancement/T", 0.01);
  
  LifeV::SolverAztec solver;
  solver.setOptionsFromGetPot(datafile);
  solver.setMatrix(M);
  
  for (Real t = t0; t <= T; t += h) {
    cout << "t = " << t << endl;
    
    // RK step 1
    rhs = A * U;
    solver.solve(K1, rhs);
    K1 = U - h * K1;
    
    // RK step 2
    rhs = A * K1;
    solver.solve(K2, rhs);
    K2 = K1 - h * K2;
    
    // RK Final solution
    U = 0.5 * U + 0.5 * K2;
  }
  
  
  elem_wise_wrtr("sphere_final.dx", mesh, dof, feDG, PN, Conn, U);
  cout << "** Solution written on file sphere.dx" << endl;
  return 0;
}
