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
#include <life/lifesolver/NavierStokesAleSolverCT.hpp>

namespace LifeV
{

  using namespace std;
  
  
  void compute_divuq(Real alpha, ElemVec& uLoc,  ElemVec& elvec, const CurrentFE& fe_u, const CurrentFE& fe_p, int iblock  )
  {
    int i, j, ic, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    
    for (i = 0; i < fe_p.nbNode; i++) {
      s = 0;
      for (iq = 0; iq < fe_p.nbQuadPt; ++iq )
	for (j = 0; j < fe_u.nbNode; ++j)
	  for (ic = 0; ic < (int)nDimensions; ++ic) 
	    s += uLoc[ic*fe_u.nbNode+j]*fe_u.phiDer(j,ic,iq)*fe_p.phi(i,iq)* fe_p.weightDet( iq );
      
      vec( i ) += s*alpha;
    }
    
  }

  void compute_gradpv(Real alpha, ElemVec& pLoc,  ElemVec& elvec, const CurrentFE& fe_p, const CurrentFE& fe_u, int iblock )
  {
    int i, j, iq;
    ElemVec::vector_view vec = elvec.block( iblock );
    Real s;
    
    for ( i = 0;i < fe_u.nbNode;i++ )
      {
        s = 0;
        for (iq = 0; iq < fe_u.nbQuadPt; ++iq )
	  for (j = 0; j < fe_p.nbNode; ++j)
	    s += pLoc[j]*fe_p.phiDer(j,iblock,iq)*fe_u.phi(i,iq)*fe_u.weightDet( iq );
	vec( i ) += s*alpha;
      }
  }
}
