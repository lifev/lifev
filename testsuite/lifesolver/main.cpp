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

  Shows how ParabolicSolver can be invoked, to solve
  by iterations the problem

         \nu(t,x,y,z,u)\partial u/\partial t - \nabla\cdot(k(t,x,y,z,u)u)+
                                 \sigma(t,x,y,z,u)u=f(t,x,z,y,t)

  with on boundary
     u = guD(t,x,y,z,u)
  or
         \paritial u/\partial n=guN(t,x,y,z,u)

  on a cube
*/


#include <life/lifecore/GetPot.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefilters/openDX_wrtrs.hpp>
#include <life/lifefilters/vtk_wrtrs.hpp>
#include <vector>
#include <algorithm>
#include <life/lifecore/life.hpp>

#include <life/lifefilters/vtk_wrtrs.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifemesh/markers.hpp>
#include <life/lifealg/dataAztec.hpp>


#include <life/lifesolver/timeSolver.hpp>
#include <life/lifesolver/parabolicSolver.hpp>
#include "ud_functions.hpp"

#undef  OPER_TEMPLATE
//#define P1
#define P2
#undef INRIA

int
main ()
{
  using namespace LifeV;
  using namespace std;


  Real timeStep;
  Chrono chrono;
  UInt i;

  // ===================================================
  // Boundary conditions definition
  //   (UDep means depending on solution)
  // ===================================================

  BCFunctionUDepBase gvu1 (gu1);
  BCFunctionUDepBase gvu3 (gu3);
  BCFunctionUDepBase gvu2 (gu2);
  BCHandler BCh (3);

  BCh.addBC ("Lateral", 1, Natural, Scalar, gvu1);
  BCh.addBC ("Left", 3, Essential, Scalar, gvu3);
  BCh.addBC ("Right", 2, Essential, Scalar, gvu2);





  // Ouput
  BCh.showMe ();

  // ===================================================
  // Finite element stuff
  // ===================================================
  const GeoMap & geoMap = geoLinearTetra;
  const QuadRule & qr = quadRuleTetra5pt;

  const GeoMap & geoMapBd = geoLinearTria;
  const QuadRule & qrBd = quadRuleTria3pt;

#ifndef P2
  // P1 elements
  const RefFE & refFE = feTetraP1;
  const RefFE & refBdFE = feTriaP1;
  RegionMesh3D < LinearTetra > mesh;
  typedef RegionMesh3D < LinearTetra > mesh_type;
#else
  //P2 elements
  const RefFE & refFE = feTetraP2;
  const RefFE & refBdFE = feTriaP2;
  RegionMesh3D < QuadraticTetra > mesh;
  typedef RegionMesh3D < QuadraticTetra > mesh_type;
#endif

  // ===================================================
  // Mesh stuff
  // ===================================================

  long int m = 1;

  GetPot datafile ("data");
  std::string mesh_format = datafile ("mesh_format", "NETGEN");
  string mesh_dir = datafile ("mesh_dir", ".");
  string fname = mesh_dir + datafile ("mesh_file", ".");
  if (mesh_format == "INRIA")
    {
      readINRIAMeshFile (mesh, fname, m);
    }
  else if (mesh_format == "MESH++")
    {
      readMppFile (mesh, fname, m);
    }
  else if (mesh_format == "NETGEN")
    {
      readNetgenMesh (mesh, fname, m);
    }
  else
    {
      std::
    cerr << "wrong mesh type. It can be either MESH++ or INRIA or NETGEN"
    << std::endl;
      return EXIT_FAILURE;
    }

  timeStep=datafile ("timeStep", 0.5);
  cout << "timeStep: " << timeStep << endl;

  // Avaliable meshes
//  mesh.showMe ();

  cout << "Now building local Edges/faces Stuff" << endl << endl;
  mesh.updateElementEdges ();
  mesh.updateElementFaces ();
//  mesh.showMe ();

  // ===================================================
  // Current FE classes for the problem under study with
  // mapping and quadrature rules
  // ===================================================

  CurrentFE fe (refFE, geoMap, qr);
  CurrentBdFE feBd (refBdFE, geoMapBd, qrBd);

  // ===============================================
  // Update of the Dof for the particular FE problem
  // and for the boundary conditions
  // ===============================================

  Dof dof (refFE);
  dof.update (mesh);

  BCh.bdUpdate (mesh, feBd, dof);

  UInt dim = dof.numTotalDof ();

  dof.showMe ();

  // initialization of vector of unknowns and rhs
  ScalUnknown < Vector > U (dim), U0 (dim);
  U0 = ZeroVector (dim);

  //init solution at time 0: U0
  for (i = 0; i < dim; i++)
    U0[i] = 0;

  // ==========================================
  // Pattern construction and matrix assembling
  // ==========================================
  cout << "dim                    = " << dim << endl << endl;

  chrono.start ();
  ParabolicSolver < mesh_type > pSolver (nu, mu, mesh, fe, feBd, dof, BCh,
                     fct, U0, timeStep, 8);
  //here it solves, note that first iteration is the most heavy, next
  //are enought tuned to converge easily
  pSolver.automaticSolver (0.01, 12);
  U = pSolver.getU (pSolver.getStep ());
  chrono.stop ();
  cout << "U computed in " << chrono.diff () << "second\n";


  //you want an dump of last U computed
  for (UInt jj = 0; jj < dim; jj++)
    cout << U (jj) << " ** ";

  cout << endl<<"saving "<<pSolver.getStep()<<" solutions on dir sol/"<<endl;

  //I save all the solutions in directory sol
  pSolver.saveSolution ();

  return 0;
}


