/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
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
#ifndef _SOBOLEVNORMS_H_INCLUDED
#define _SOBOLEVNORMS_H_INCLUDED
#include "lifeV.hpp"
#include "dof.hpp"
#include "currentFE.hpp"

//! returns the square of the L2 norm of u on the current element
template<typename VectorType,typename DOF>
Real elem_L2_2(const VectorType & u,const CurrentFE& fe,const DOF& dof)
{
  int i,inod,ig;
  UInt eleID=fe.currentId();
  Real s=0,u_ig;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    u_ig=0.;
    for(i=0;i<fe.nbNode;i++){
      inod = dof.localToGlobal(eleID,i+1)-1;
      u_ig += u(inod) * fe.phi(i,ig);
    }
    s += u_ig * u_ig * fe.weightDet(ig);
  }
  return s;
}

//! version for vectorial problem
template<typename VectorType,typename DOF>
Real elem_L2_2(const VectorType & u,const CurrentFE& fe,const DOF& dof,
	       const int nbcomp)
{
  int i,inod,ig;
  UInt eleID=fe.currentId();
  UInt ic;
  Real s=0,u_ig;
  for (ic=0; ic< nbcomp; ic++){
    for(ig=0;ig<fe.nbQuadPt;ig++){
      u_ig=0.;
      for(i=0;i<fe.nbNode;i++){
	inod = dof.localToGlobal(eleID,i+1)-1 + ic*dof.numTotalDof();
	u_ig += u(inod) * fe.phi(i,ig);
      }
      s += u_ig * u_ig * fe.weightDet(ig);
    }
  }
  return s;
}

//! returns the square of the L2 norm of fct on the current element  
template<typename UsrFct>
Real elem_L2_2(const UsrFct& fct,const CurrentFE& fe)
{
  Real s=0.,f,x,y,z;
  for(int ig=0;ig<fe.nbQuadPt;ig++){
    fe.coorQuadPt(x,y,z,ig);
    f = fct(x,y,z);
    s += f * f * fe.weightDet(ig);
  }
  return s;
}

//! for time dependent+vectorial.
template<typename UsrFct>
Real elem_L2_2(const UsrFct& fct,const CurrentFE& fe, const Real t,
	       const UInt nbcomp)
{
  int ig;
  UInt ic;
  Real s=0.,f,x,y,z;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    fe.coorQuadPt(x,y,z,ig);
    for (ic=0; ic<nbcomp; ic++){
      f = fct(x,y,z,t,ic);
      s += f * f * fe.weightDet(ig);
    }
  }
  return s;
}

//! returns the square of the H1 norm of u on the current element
template<typename VectorType,typename DOF>
Real elem_H1_2(const VectorType & u,const CurrentFE& fe,const DOF& dof)
{
  int i,inod,ig;
  UInt eleID=fe.currentId();
  Real s=0,u_ig,dx_u_ig,dy_u_ig,dz_u_ig;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    u_ig=0.;
    dx_u_ig = 0.;
    dy_u_ig = 0.;
    dz_u_ig = 0.;
    for(i=0;i<fe.nbNode;i++){
      inod = dof.localToGlobal(eleID,i+1)-1;
      u_ig += u(inod) * fe.phi(i,ig);
      dx_u_ig += u(inod)*fe.phiDer(i,0,ig);
      dy_u_ig += u(inod)*fe.phiDer(i,1,ig);
      dz_u_ig += u(inod)*fe.phiDer(i,2,ig);
    }
    s += (u_ig * u_ig + dx_u_ig * dx_u_ig + dy_u_ig * dy_u_ig +
	  dz_u_ig * dz_u_ig) * fe.weightDet(ig);
  }
  return s;
}

//! returns the square of the H1 norm of fct on the current element
template<typename UsrFct>
Real elem_H1_2(const UsrFct& fct,const CurrentFE& fe)
{
  Real s=0.,s1=0.,f,fx,fy,fz,x,y,z;
  for(int ig=0;ig<fe.nbQuadPt;ig++){
    fe.coorQuadPt(x,y,z,ig);
    f = fct(x,y,z);
    fx = fct.der_x(x,y,z);
    fy = fct.der_y(x,y,z);
    fz = fct.der_z(x,y,z);
    s += f * f * fe.weightDet(ig);
    s1 += (fx*fx + fy*fy + fz*fz) * fe.weightDet(ig);
  }
  return s+s1;
}

//! returns the square of the H1 norm of fct on the current element (time-dependent case)
template<typename UsrFct>
Real elem_H1_2(const UsrFct& fct,const CurrentFE& fe,const Real t, const UInt nbcomp)
{
  Real s=0.,s1=0.,f,fx,fy,fz,x,y,z;
  for(int ig=0;ig<fe.nbQuadPt;ig++){
    fe.coorQuadPt(x,y,z,ig);
    for (UInt ic=0;ic<nbcomp;ic++){
    f = fct(x,y,z,t,ic);
    fx = fct.der_x(x,y,z,t,ic);
    fy = fct.der_y(x,y,z,t,ic);
    fz = fct.der_z(x,y,z,t,ic);
    s += f * f * fe.weightDet(ig);
    s1 += (fx*fx + fy*fy + fz*fz) * fe.weightDet(ig);
   }
  }
  return s+s1;
}


//! returns the square of the L2 norm of (u-fct) on the current element  
template<typename VectorType,typename UsrFct,typename DOF>
Real elem_L2_diff_2(const VectorType & u,const UsrFct& fct,const CurrentFE& fe,
		    const DOF& dof)
{
  int i,inod,ig;
  UInt eleID=fe.currentId();
  Real s=0.,u_ig,u_minus_f,x,y,z;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    u_ig=0.;
    for(i=0;i<fe.nbNode;i++){
      inod = dof.localToGlobal(eleID,i+1)-1;
      u_ig += u(inod) * fe.phi(i,ig);
    }
    fe.coorQuadPt(x,y,z,ig);
    u_minus_f = u_ig - fct(x,y,z);
    s += u_minus_f * u_minus_f * fe.weightDet(ig);
  }
  return s;
}

//! for time dependent+vectorial
template<typename VectorType,typename UsrFct,typename DOF>
Real elem_L2_diff_2(const VectorType & u,const UsrFct& fct,const CurrentFE& fe,
		    const DOF& dof, const Real t, const int nbcomp)
{
  // returns the square of the L2 norm of (u-fct) on the current element  
  int i,inod,ig;
  UInt eleID=fe.currentId();
  int ic;
  Real s=0.,u_ig,u_minus_f,x,y,z;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    for (ic= 0; ic <nbcomp; ic++){
      u_ig=0.;
      for(i=0;i<fe.nbNode;i++){
	inod = dof.localToGlobal(eleID,i+1)-1+ic*dof.numTotalDof();
	u_ig += u(inod) * fe.phi(i,ig);
      }
      fe.coorQuadPt(x,y,z,ig);
      u_minus_f = u_ig - fct(x,y,z,t,ic);
      s += u_minus_f * u_minus_f * fe.weightDet(ig);
    }
  }
  return s;
}

//! returns the square of the H1 norm of (u-fct) on the current element  
template<typename VectorType,typename UsrFct,typename DOF>
Real elem_H1_diff_2(const VectorType & u,const UsrFct& fct,const CurrentFE& fe,
		    const DOF& dof)
{
  int i,inod,ig;
  UInt eleID=fe.currentId();
  Real s=0.,u_ig,u_minus_f,x,y,z,dx_u_ig,dy_u_ig,dz_u_ig,u_minus_fx,u_minus_fy,u_minus_fz;
  for(ig=0;ig<fe.nbQuadPt;ig++){
    u_ig=0.;
    dx_u_ig = 0.;
    dy_u_ig = 0.;
    dz_u_ig = 0.;

    for(i=0;i<fe.nbNode;i++){
      inod = dof.localToGlobal(eleID,i+1)-1;
      u_ig += u(inod) * fe.phi(i,ig);
      dx_u_ig += u(inod)*fe.phiDer(i,0,ig);
      dy_u_ig += u(inod)*fe.phiDer(i,1,ig);
      dz_u_ig += u(inod)*fe.phiDer(i,2,ig);
    }
    fe.coorQuadPt(x,y,z,ig);
    u_minus_f = u_ig - fct(x,y,z);
    u_minus_fx = dx_u_ig - fct.der_x(x,y,z);
    u_minus_fy = dy_u_ig - fct.der_y(x,y,z);
    u_minus_fz = dz_u_ig - fct.der_z(x,y,z);

    s += (u_minus_f * u_minus_f +
	  u_minus_fx * u_minus_fx +
	  u_minus_fy * u_minus_fy +
	  u_minus_fz * u_minus_fz )
      * fe.weightDet(ig);
  }
  return s;
}

//! returns the square of the H1 norm of (u-fct) on the current element  (time-dependent case) 
template<typename VectorType,typename UsrFct,typename DOF>
Real elem_H1_diff_2(const VectorType & u,const UsrFct& fct,const CurrentFE& fe,
		    const DOF& dof, const Real t, const UInt nbcomp)
{
  int i,inod,ig;
  UInt ic;
  UInt eleID=fe.currentId();
  Real s=0.,u_ig,u_minus_f,x,y,z,dx_u_ig,dy_u_ig,dz_u_ig,u_minus_fx,u_minus_fy,u_minus_fz;
  for(ig=0;ig<fe.nbQuadPt;ig++){
   for (ic= 0; ic <nbcomp; ic++){
    u_ig=0.;
    dx_u_ig = 0.;
    dy_u_ig = 0.;
    dz_u_ig = 0.;

    for(i=0;i<fe.nbNode;i++){
      inod = dof.localToGlobal(eleID,i+1)-1+ic*dof.numTotalDof();
      u_ig += u(inod) * fe.phi(i,ig);
      dx_u_ig += u(inod)*fe.phiDer(i,0,ig);
      dy_u_ig += u(inod)*fe.phiDer(i,1,ig);
      dz_u_ig += u(inod)*fe.phiDer(i,2,ig);
    }
    fe.coorQuadPt(x,y,z,ig);
    u_minus_f = u_ig - fct(x,y,z,t,ic);
    u_minus_fx = dx_u_ig - fct.der_x(x,y,z,t,ic);
    u_minus_fy = dy_u_ig - fct.der_y(x,y,z,t,ic);
    u_minus_fz = dz_u_ig - fct.der_z(x,y,z,t,ic);

    s += (u_minus_f * u_minus_f +
	  u_minus_fx * u_minus_fx +
	  u_minus_fy * u_minus_fy +
	  u_minus_fz * u_minus_fz )
      * fe.weightDet(ig);
   }
  }
  return s;
}

#endif
