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
#define DISCGAL
#define PI 3.141592653589793283

#include <vector>
#include <algorithm>
#include <cmath>

#include "ud_functions.hpp"

#include <life/lifefem/geoMap.hpp>

#include <life/lifefem/refEleDG.hpp>
#include <life/lifefem/refFEDG.hpp>

#include <life/lifefem/currentFEDG.hpp>
#include <life/lifefem/currentIFDG.hpp>
#include <life/lifefem/currentBFDG.hpp>

#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/readMesh3D.hpp>

#include <life/lifefem/dofDG.hpp>
#include <life/lifefem/dofByFace.hpp>
//#include "elementAdjacency.hpp"

#include <life/lifearray/pattern.hpp>

#include <life/lifefem/values.hpp> // Contiene la definizione della classe MSRMatr

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>

#include <life/lifefem/elemOper_ext.hpp>

#include <life/lifefem/assemb.hpp>

namespace LifeV
{
//============================================================================
// Velocity function
//============================================================================

class VectorFunction {
 public:
  virtual Real operator()(Real x, Real y, Real z, UInt comp) {return 0;};
};

class Vortex: public VectorFunction {
 public:
  Real operator()(Real& x, Real& y, Real& z, int& comp) {
    Real s = 0.;

    switch(comp){
    case 0:
      s = sin(PI * x) * sin(PI * x) * (sin(2 * PI * z) - sin(2 * PI * y));
      break;
    case 1:
      s = sin(PI * y) * sin(PI * y) * (sin(2 * PI * x) - sin(2 * PI * z));
      break;
    case 2:
      s = sin(PI * z) * sin(PI * z) * (sin(2 * PI * y) - sin(2 * PI * x));
      break;
    }

    return s;
  }
};

//============================================================================
// Source function
//============================================================================

class SourceFct
{
 public:
  inline Real operator()(Real x,Real y,Real z,int ic=0) const {
    return 0.;
  }
};

//============================================================================
// Operators related to advection terms
//============================================================================

template<typename Velocity>
class AdvecDG{
 public: 
  AdvecDG(CurrentFEDG* fe):_fe(fe){};
 
  CurrentFEDG* fe_ptr() {return _fe;};

  Real operator()(int i, int j, int iq, Real x, Real y, Real z, int K, int H, int ic, int jc)
    {
      Real s = 0.;

      for(int icoor = 0; icoor < (int)_fe -> nbCoor; icoor++){
    s += (- _fe -> phi(j, iq) * _fe -> phiDer(i, icoor, iq) * u(x, y, z, icoor));
      }

      return s;
    }
 private:
  CurrentFEDG* _fe;
  Velocity u;
};

template<typename Velocity>
class AdvecIFUW1DG {
 public: 
  AdvecIFUW1DG(CurrentIFDG* fe):_fe(fe){};
 
  CurrentIFDG* fe_ptr() {return _fe;};

  Real operator()(int i, int j, int iq, Real x, Real y, Real z, int K, int H, int ic, int jc)
    {
      Real s = 0.;

      JumpIF jump(_fe);
      AvgIF avg(_fe);

      for(int icoor = 0; icoor < _fe -> nbCoorAd; icoor++){
    s += u(x, y, z, icoor) * jump(i, icoor, iq, H) * avg(j, iq, K); 
      }

      return s;
    }

 private:
  CurrentIFDG* _fe;
  Velocity u;
};

template<typename Velocity>
class AdvecIFUW2DG{
 public: 
  AdvecIFUW2DG(CurrentIFDG* fe):_fe(fe){};
 
  CurrentIFDG* fe_ptr() {return _fe;};

  Real operator()(int i, int j, int iq, Real x, Real y, Real z, int K, int H, int ic, int jc)
    {
      Real s = 0.;
      Real provv = 0.;
      
      JumpIF jump(_fe);

      for(int icoor = 0; icoor < _fe -> nbCoorAd; icoor++)
    provv += u(x, y, z, icoor) * _fe -> normal(icoor, iq);

      provv = 0.5 * fabs(provv);

      for(int icoor = 0; icoor < _fe -> nbCoorAd; icoor++)
    s += provv * jump(i, icoor, iq, H) * jump(j, icoor, iq, K);

      return s;
    }

 private:
  CurrentIFDG* _fe;
  Velocity u;
};

template<typename Velocity>
class AdvecBFUWDG{
 public: 
  AdvecBFUWDG(CurrentBFDG* fe):_fe(fe){};
 
  CurrentBFDG* fe_ptr() {return _fe;};

  Real operator()(int i, int j, int iq, Real x, Real y, Real z, int K, int H, int ic, int jc)
    {
      Real s = 0.;
      Real provv = 0.;
      
      JumpBF jump(_fe);

      for(int icoor = 0; icoor < _fe -> nbCoorAd; icoor++)
    provv += u(x, y, z, icoor) * _fe -> normal(icoor, iq);

      if(provv >= 0){
    for(int icoor = 0; icoor < _fe -> nbCoorAd; icoor++)
      s += provv * jump(i, icoor, iq, H) * jump(j, icoor, iq, K);
      }

      return s;
    }

 private:
  CurrentBFDG* _fe;
  Velocity u;
};

//============================================================================
// Finite element stuff
//============================================================================

const GeoMapDG geoLinearTetraDG("Linear mapping on a tri (DG)", TETRA, 4, 3,
                fct_P1_3D, derfct_P1_3D, der2fct_P1_3D,
                refcoor_P1_3D,
                allQuadRuleTetra, TRIANGLE,
                4, 3,
                refCoorFaces_P1_DG_3D,
                allQuadRuleTria, geoLinearTria, &geoLinearTria);
               
RefEleDG MyRefEleDG("My Ref Ele DG", TETRA, 4, 3, 
            fct_P1_3D, derfct_P1_3D, der2fct_P1_3D, refcoor_P1_3D, 
            allQuadRuleTetra, TRIANGLE, 4, 3, refCoorFaces_P1_DG_3D, 
            allQuadRuleTria, geoLinearTria);

const LocalDofPattern elPattern_P1_DG_3D(4, 1, 0, 0, 0, STANDARD_PATTERN); 

const LocalDofPattern facePattern_P1_DG_3D(8, 1, 0, 0, 0, STANDARD_PATTERN);

const RefFEDG feDGTetraP1("DG P1 on a tetrahedra", FE_DG_P1_3D, TETRA,
              1, 0, 0, 0, 4, 3,
              fct_P1_3D, derfct_P1_3D, der2fct_P1_3D, 
              refcoor_P1_3D, 
              allQuadRuleTetra, elPattern_P1_DG_3D, &feTriaP1,
              TRIANGLE,
              4, 3,
              refCoorFaces_P1_DG_3D,
              allQuadRuleTria, facePattern_P1_DG_3D, geoLinearTria);
              
//============================================================================
// Operators' stuff
//============================================================================

typedef EOExpr<Real, AdvecDG<Vortex> > EOAdvecDG;

typedef EOExpr<Real, AdvecIFUW1DG<Vortex> > EOAdvecIFUW1DG;

typedef EOExpr<Real, AdvecIFUW2DG<Vortex> > EOAdvecIFUW2DG;

typedef EOExpr<Real, AdvecBFUWDG<Vortex> > EOAdvecBFUWDG;

//============================================================================
// OpenDX writer suitable for DG-FEM
//============================================================================
template<class MESH, class CURRENTFE, class DOF, class VECTOR>
void elem_wise_wrtr(std::string fname, const MESH& mesh, const DOF& dof, CURRENTFE& fe, 
                    const KNM<Real>& PlotNodes, const KNM<UInt>& Connections, const VECTOR& U)
{
  // WARNING: Connections must start from 0
  
  // Open output file
  std::ofstream ofile( fname.c_str(), std::ios::app );
  
  ASSERT( ofile, "Error: Output file cannot be open" );
  
  // Number of coordinates
  ASSERT_PRE(fe.nbCoor == PlotNodes.M(), "At present time works in 3D only");
  
  // Number of local and total plot nodes, tetrahedra
  UInt nE = mesh.numVolumes();
  const int nbPlotNodes = PlotNodes.N();
  int nbTotalPlotNodes = nbPlotNodes * nE;
  int nbConn = Connections.N();
  
  // Number of (local) dofs
  const int nbDof = fe.refFE.nbDof;
  
  // Multi-index array storing shape functions' values on plot nodes
  KNM<Real> PhiOnPlotNodes(nbDof, nbPlotNodes);
  
  // Variables to store plot nodes' coordinates on the reference element
  Real xi, eta, zeta;
  
  // Variables to store plot nodes' physical coordinates
  Real x, y, z;
  
  for(int i = 0; i < PhiOnPlotNodes.N(); i++)
    for(int j = 0; j < PhiOnPlotNodes.M(); j++)
    {
      xi = PlotNodes(j, 0); eta = PlotNodes(j, 1); zeta = PlotNodes(j, 2);
      PhiOnPlotNodes(i, j) = fe.refFE.phi(i, xi, eta, zeta);
    }
  
  // ** Write positions **
  ofile << "# POSITIONS" << std::endl;
  ofile << "object 1 class array type float rank 1 shape 3 items " << nbTotalPlotNodes << " data follows" << std::endl;
  for(UInt iE = 1; iE <= nE; iE++)
  {
    for(int iPN = 0; iPN < nbPlotNodes; iPN++)
    {
      xi = PlotNodes(iPN, 0); eta = PlotNodes(iPN, 1); zeta = PlotNodes(iPN, 2);
      
      // Update element
      fe.update( mesh.volumeList(iE) );
      
      // Plot node's physical coordinates
      fe.coorMap(x, y, z, xi, eta, zeta);
      
      ofile << (float)x << " " << (float)y << " " << (float)z << std::endl;
    }
  }
  
  // ** Write connections **
  ofile << "# CONNECTIONS" << std::endl;
  int nbTotalConn = nE * nbConn;
  ofile << "object 2 class array type int rank 1 shape 4 items " << nbTotalConn << " data follows" << std::endl;

  UInt nPN = 0;
  for(UInt iE = 1; iE <= nE; iE++)
  {
    for(int iConn = 0; iConn < nbConn; iConn++)
    {
      for(int i = 0; i < 4; i++)
      {
        ofile << Connections(iConn, i) + nPN << " ";
      }
      ofile << std::endl;
    }
    nPN += nbPlotNodes;
  } 
  ofile << "attribute \"element type\" string \"tetrahedra\"" << std::endl;
  ofile << "attribute \"ref\" string \"positions\"" << std::endl;
  
  // ** Write data **
  Real UOnPlotNode;
  ofile << "# DATA" << std::endl;
  ofile << "object 3 class array type float rank 0 items " << nbTotalPlotNodes << " data follows" << std::endl;
  UInt ig;
  
  for(UInt iE = 1; iE <= nE; iE++)
  {
    for(int iPN = 0; iPN < nbPlotNodes; iPN++)
    {
      // Solution value on plot node
      UOnPlotNode = 0.;      
      for(int iDof = 0; iDof < nbDof; iDof++)
      {
        ig = dof.localToGlobal(iE, iDof + 1) - 1;
        UOnPlotNode += PhiOnPlotNodes(iDof, iPN)  * U(ig);
      }
      ofile << (float)UOnPlotNode << std::endl;
    }
  }
  
  ofile << "attribute \"dep\" string \"positions\"" << std::endl;
  ofile << "# FIELD" << std::endl;
  ofile << "object \"irregular positions irregular connections\" class field" << std::endl;
  ofile << "component \"positions\" value 1" << std::endl;
  ofile << "component \"connections\" value 2" << std::endl;
  ofile << "component \"data\" value 3" << std::endl;
  ofile << "end" << std::endl;
  
  ofile.close();
  }

}
