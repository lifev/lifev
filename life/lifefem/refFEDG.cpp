/*-*- mode: c++ -*-
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <life/lifefem/refFEDG.hpp>
#include <set>

namespace LifeV
{
  RefFEDG::RefFEDG(std::string _name, int _type,
		   ReferenceShapes _shape,
		   int /*_nbDofPerVertex*/, int /*_nbDofPerEdge*/, int /*_nbDofPerFace*/, int /*_nbDofPerVolume*/, int _nbDof,
		   int _nbCoor,
		   const Fct* phi, const Fct* dPhi, const Fct* d2Phi, 
		   const Real* refCoor, 
		   const SetOfQuadRule& sqr, const LocalDofPattern& _elPattern, const RefFE* boundaryFE,
		   ReferenceShapes _shapeFaces, 
		   int _nbFaces, int _nbGeoNodeFaces, 
		   const Real* refCoorFaces, 
		   const SetOfQuadRule& sqrFaces, const LocalDofPattern& _facePattern, const GeoMap _geoMap):
    RefEleDG( _name, _shape, _nbDof, _nbCoor, phi, dPhi, d2Phi, refCoor, sqr, _shapeFaces, _nbFaces, _nbGeoNodeFaces, refCoorFaces, sqrFaces, _geoMap),
    _boundaryFE(boundaryFE), 
    elPattern(_elPattern),
    facePattern(_facePattern),
    type(_type)
  {
    CONSTRUCTOR("RefFEDG");
  }

  RefFEDG::~RefFEDG()
  {
    DESTRUCTOR("RefFEDG");
  }

  std::ostream& operator << (std::ostream& f, const RefFEDG& fe)
  {
    f << "--------------------------------------------------------------------------------" << std::endl;
    f << "Reference Discontinuous Finite Element: " << fe.name << std::endl;
    f << "--------------------------------------------------------------------------------" << std::endl;

    f << std::endl << "*** Shape : " << fe.shape << std::endl;

    f << std::endl << "*** Local coordinate of the nodes :\n";
    for(int i=0;i<fe.nbDof;i++){
      f << fe.xi(i) << " " << fe.eta(i) << " " << fe.zeta(i) << std::endl;
    }
    f << std::endl << "*** Local coordinates of the face vertices : \n";
    for(int iFace = 0; iFace < fe.nbFaces; iFace++){
      f << "    - Coordinates of face : " << iFace << std::endl;
      for(int i = 0; i < fe.nbGeoNodeFaces; i++){
	f << fe.xiFace(iFace,i) << " " << fe.etaFace(iFace,i) << " " << fe.zetaFace(iFace,i) << std::endl;
      } // for i
    } // for iFace
   
    f << std::endl << "*** Pattern :\n";
    //f << "*** Pattern :\n";
    //for(int i=0;i<fe.nbPattern();i++) f << "(" << fe.patternFirst(i) << "," << fe.patternSecond(i) << ") \n";
    //f << std::endl << fe.elPattern;
    //f << std::endl << fe.facePattern;

    for(int k=0;k<fe._sqr->nbQuadRule;k++){
      const QuadRule& qr = fe._sqr->quadRule(k);
      f << std::endl << "*** Quadrature rule : " << qr.name << std::endl;
      for(int ig=0;ig<qr.nbQuadPt;ig++){
	f << "    - Quadrature point : " << ig << std::endl;
	//      f << "     number and values of basis functions = " << fe.phiQuadPt(ig,qr) << std::endl;
	for(int i=0;i<fe.nbDof;i++){
	  f << "      Basif fct " << i << std::endl;
	  f << "         Value = " << fe.phi(i,ig,qr) << std::endl;
	  f << "         Derivatives = " ;
	  for(int icoor=0;icoor<fe.nbCoor;icoor++) f << " " << fe.dPhi(i,icoor,ig,qr);
	  f << std::endl;
	  f << "         Second derivatives = " ;
	  for(int icoor=0;icoor<fe.nbCoor;icoor++)
	    for(int jcoor=0;jcoor<fe.nbCoor;jcoor++)
	      f << " " << fe.d2Phi(i,icoor,jcoor,ig,qr);
	  f << std::endl;
	}
      }
    }

    // This part is for discontinuous elements only
    f << std::endl << "*** Face shape: " << fe.shapeFaces << std::endl;
    f << std::endl << "*** Reference coordinates for the Faces  :" << std::endl;
    for(int iFace = 0; iFace < fe.nbFaces; iFace++){
      f << std::endl << "*** Face #" << iFace << std::endl;
      for(int i = 0; i < fe.nbGeoNodeFaces; i++){
	f << fe.xiFace(iFace,i) << " " << fe.etaFace(iFace,i) << " " << fe.zetaFace(iFace,i) << std::endl;
      }
  
      for(int k = 0; k < fe._sqrFaces -> nbQuadRule; k++){
	const QuadRule& qrFaces = fe._sqrFaces -> quadRule(k);
	f << std::endl << "*** Quadrature rule : " << qrFaces.name << std::endl;
	for(int ig = 0; ig < qrFaces.nbQuadPt; ig++){
	  f << "   - Quadrature point : " << ig << std::endl;
	  for(int i = 0; i < fe.nbDof; i++){
	    f << "      Basis fct " << i << std::endl;
	    f << "         Value = " << fe.phiFace(iFace, i, ig, qrFaces) << std::endl;
	    f << "         Derivatives = ";
	    for(int icoor = 0; icoor < fe.nbCoor; icoor++)
	      f << " " << fe.dPhiFace(iFace, i, icoor, ig, qrFaces);
	    f << std::endl;
	    f << "         Second derivatives = ";
	    for(int icoor = 0; icoor < fe.nbCoor; icoor++)
	      for(int jcoor = 0; jcoor < fe.nbCoor; jcoor++)
		f << " " << fe.d2PhiFace(iFace, i, icoor, jcoor, ig, qrFaces);
	    f << std::endl;
	  } // for i
	} // for ig
      } // for k
    } // for iFace
    return f;
  }
}
