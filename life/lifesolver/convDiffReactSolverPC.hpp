/*
  This file is part of the LifeV library
  Copyright (C) 2004 EPFL, INRIA and Politechnico di Milano

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
/*!
  \file convDiffReactSolverPC.h
  \author M. Prosi
  \date 03/2004
  \version 1.0

  \brief This file contains a solver class for the Convection-Diffusion-Reaction equation
*/

#ifndef _CONVDIFFREACTSOLVERPC_H_
#define _CONVDIFFREACTSOLVERPC_H_

#include "convDiffReactHandler.hpp"
#include "elemMat.hpp"
#include "elemVec.hpp"
#include "elemOper.hpp"
#include "values.hpp"
#include "pattern.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"
#include "algebraic_facto.hpp"
#include "bcCond.hpp"
#include "chrono.hpp"
#include "dataAztec.hpp"
#include "bdf.hpp"
#include "openDX_wrtrs.hpp"
#include <string>

namespace LifeV
{
/*!
  \class convDiffReactSolverPC

   This class implements a solver for the Convection-Diffusion-Reaction equation

*/
template<typename Mesh>
class ConvDiffReactSolverPC:
public ConvDiffReactHandler<Mesh> {

 public:

  typedef  typename  ConvDiffReactHandler<Mesh>::Function Function;

  //! Constructor
  /*!
    \param data_file GetPot data file
    \param refFE_c reference FE for the concentration
    \param Qr_c volumic quadrature rule for the concentration
    \param bdQr_c surface quadrature rule for the concentration
    \param BCh_c boundary conditions for the concentration
  */
  ConvDiffReactSolverPC(const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
	    const QuadRule& bdQr_c, BC_Handler& BCh_c);

  //! Update the right  hand side  for time advancing
  /*!
    \param source volumic source
    \param time present time
  */
  void timeAdvance(const Function source, const Real& time);

  //! Update convective term, bc treatment and solve the linearized ns system
  void iterate(const Real& time);

  //! Projection of the velocity on grid of concentration discretization
  template <typename RegionMesh3D>
  void getvel(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u, const Real& time);

  //! Calculate the local coordinates of concentration gridpoints in the
  //! velocity grid (is needed for the Projection)
  template <typename RegionMesh3D>
  void getcoord(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u);

  //! Calculate the volume of a tetrahedra given by its corner nodes
  Real calcvol(Real x[4], Real y[4], Real z[4]);

  //! tests if point (xp, yp, zp) is in the tetrahedra (x[4], y[4], z[4]) and returns
  //! interpolation coefficents (1-b1-b2-b3, b1, b2, b3)
  int test(Real x[4], Real y[4], Real z[4], Real & xp, Real & yp, Real & zp, Real & b1, Real & b2, Real & b3);

protected:

   // inherited from parent class
   typedef typename ConvDiffReactHandler<Mesh>::intpolcoord intpolcoord;

 private:

  //! Pattern of M
  MSRPatt _pattM;

  //! Matrix C:  1/dt*Cmass + D*Cstiff operator + r*Cmass
  MSRMatr<double> _DR;

  //! Matrix C:  1/dt*Cmass + D*Cstiff operator + Convective_transport term + r*Cmass
  MSRMatr<double> _CDR;

  //! Matrix C_u: Cmass
  MSRMatr<double> _M_c;

  //! Elementary matrices and vectors
  ElemMat _elmatC; //Concentration stiffnes
  ElemMat _elmatM_c; //Concentration mass
  ElemMat _elmatMR_c; // Concentration mass and reaction
  ElemVec _elvec; // Elementary right hand side
  ElemVec _elvec_u; // Elementary velocity for convection term

  //! Right  hand  side for the concentration
  ScalUnknown<Vector> _f_c;
  ScalUnknown<Vector> _f_c_noBC;

  DataAztec _dataAztec_o;
};

//
//                                         IMPLEMENTATION
//
template<typename Mesh> ConvDiffReactSolverPC<Mesh>::
ConvDiffReactSolverPC(const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
		      const QuadRule& bdQr_c, BC_Handler& BCh_c):
  ConvDiffReactHandler<Mesh>(data_file,refFE_c,Qr_c, bdQr_c, BCh_c),
     _pattM(_dof_c),
     _DR(_pattM),
     _CDR(_pattM),
     _M_c(_pattM),
     _elmatC(_fe_c.nbNode,1,1),
     _elmatM_c(_fe_c.nbNode,1,1),
     _elmatMR_c(_fe_c.nbNode,1,1),
     _elvec(_fe_c.nbNode,1),
     _elvec_u(_fe_c.nbNode,nDimensions),
     _f_c(_dim_c),
     _dataAztec_o(data_file,"masstransport/aztec_o"){

  std::cout << endl;
  std::cout << "O-  Concentration unknowns: " << _dim_c     << std::endl;
  std::cout << "O-  Computing mass and stiffness matrices... ";

  Chrono chrono;
  chrono.start();

  // Matrices initialization
  _DR.zeros();
  _CDR.zeros();
  _M_c.zeros();

  //inverse of dt:
  Real dti=1./_dt;

  // *******************************************************
  // Coefficient of the mass term at time t^{n+1}
  Real first_coeff = _bdf.coeff_der(0);
  std::cout << endl;
  std::cout << "Bdf CDR first coeff " << first_coeff << std::endl;

  _bdf.showMe();

  // Elementary computation and matrix assembling

  for(UInt i = 1; i <= _mesh.numVolumes(); i++){          // Loop on elements

    _fe_c.updateFirstDerivQuadPt(_mesh.volumeList(i));

    _elmatC.zero();
    _elmatM_c.zero();
    _elmatMR_c.zero();

    stiff(_diffusivity,_elmatC,_fe_c);

   if(_stationary){
      mass((-_react),_elmatMR_c,_fe_c);}
   else{
      mass((first_coeff*dti-_react),_elmatMR_c,_fe_c);
      mass(dti,_elmatM_c,_fe_c);}

    _elmatC.mat() += _elmatMR_c.mat();

    // stiffness + mass + reaction term
    assemb_mat(_DR,_elmatC,_fe_c,_dof_c);

    // mass
    assemb_mat(_M_c,_elmatM_c,_fe_c,_dof_c);

  }

  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;

}

template<typename Mesh>
void ConvDiffReactSolverPC<Mesh>::
timeAdvance(const Function source, const Real& time) {

  std::cout << "  o-  Updating mass term on right hand side (concentration)... ";

  Chrono chrono;
  chrono.start();

  // Right hand side for the velocity at time
  _f_c=0.;

  // loop on volumes: assembling source term
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){
     _elvec.zero();
     _fe_c.update(_mesh.volumeList(i));

      compute_vec(source,_elvec,_fe_c,time,0); // compute local vector
      assemb_vec(_f_c,_elvec,_fe_c,_dof_c,0); // assemble local vector into global one
  }

  // *******************************************************
  _f_c += _M_c*_bdf.time_der(); //_M_u is the mass matrix divided by the time step
  _f_c_noBC=_f_c;
  chrono.stop();
  std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template<typename Mesh>
void ConvDiffReactSolverPC<Mesh>::
iterate(const Real& time) {

  UInt nc_u=_u_c.nbcomp();

  Chrono  chrono;

  // CDR = DR + convective term (C)
  chrono.start();
  _CDR=_DR;
  _f_c=_f_c_noBC;
  chrono.stop();


  std::cout << "  o-  Diffusion-Reaction matrix was copied in " << chrono.diff() << "s." << std::endl;
  std::cout << "  o-  Updating convective transport... ";

  chrono.start();

  // loop on volumes
  for(UInt i=1; i<=_mesh.numVolumes(); ++i){

    _fe_c.updateFirstDeriv(_mesh.volumeList(i)); // as updateFirstDer

    _elmatC.zero();

    UInt eleID = _fe_c.currentId();

// ********** copy global velocity vector to local velocity vector *************
// ********** assuming velocity is given on concentration mesh *****************

    for (UInt k=0 ; k<(UInt)_fe_c.nbNode ; k++){
       UInt  iloc = _fe_c.patternFirst(k);
       for (UInt ic=0; ic<nc_u; ++ic){
	  UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;
	  _elvec_u[iloc+ic*_fe_c.nbNode] = _u_c(ig);
       }
     }

    grad(0,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);
    grad(1,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);
    grad(2,_elvec_u,_elmatC,_fe_c,_fe_c,_fe_c);


// *************************************** Upwind ******************************

      Real VLoc_infty=0.;
      Real VLoc_mean=0.;
      Real VLoc_c=0.;
      for (UInt ih_c=0 ; ih_c<(UInt)_fe_c.nbNode ; ih_c++){
         UInt  iloc = _fe_c.patternFirst(ih_c);
	 for (UInt ic=0; ic<nc_u;++ic){
	   UInt ig=_dof_c.localToGlobal(eleID,iloc+1)-1+ic*_dim_c;
           _elvec_u[iloc+ic*_fe_c.nbNode] = _u_c(ig);
	   VLoc_c+=_u_c(ig)*_u_c(ig);}
	 VLoc_c=sqrt(VLoc_c);
	 VLoc_mean += VLoc_c;
	if (VLoc_c>VLoc_infty) VLoc_infty=VLoc_c;
      }
      VLoc_mean=VLoc_mean/_fe_c.nbNode;

      Real coef_stab, Pe_loc;
//      coef_stab=_fe_c.diameter()*VLoc_infty; // Alessandro - method

      Pe_loc=VLoc_infty*_fe_c.diameter()/(2.0*_diffusivity);

//      coef_stab=(1.0/tanh(Pe_loc))-(1.0/Pe_loc); // classical approach

      if(Pe_loc < -3.0)
	 coef_stab= -1.0;
      else {
	 if(Pe_loc > 3.0)
	    coef_stab=1.0;
	 else
	    coef_stab=Pe_loc/3.0;}

// ******************************* STREAMLINEUPWIND ****************************
     stiff_sd(coef_stab/(VLoc_mean*VLoc_mean),_elvec_u,_elmatC,_fe_c,_fe_c);

// ************************* Assembling ****************************************

     assemb_mat(_CDR,_elmatC,_fe_c,_dof_c);

  }
  chrono.stop();
  std::cout << "done in " << chrono.diff() << "s." << std::endl;

  // for BC treatment (done at each time-step)
  Real tgv=1.e02;

  std::cout << "  o-  Applying boundary conditions... ";
  chrono.start();
  // BC manage for the concentration
  if ( !_BCh_c.bdUpdateDone() )
    _BCh_c.bdUpdate(_mesh, _feBd_c, _dof_c);
  bc_manage(_CDR, _f_c, _mesh, _dof_c, _BCh_c, _feBd_c, tgv, time);
  chrono.stop();

  std::cout << "done in " << chrono.diff() << "s." << std::endl;


  int    proc_config_o[AZ_PROC_SIZE];// Processor information:
  //  proc_config[AZ_node] = node name
  //  proc_config[AZ_N_procs] = # of nodes
  int    options_o[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params_o[AZ_PARAMS_SIZE];   // User selected solver paramters.
  int    *data_org_o;                // Array to specify data layout
  double status_o[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                   // indicating success or failure.
  // other declarations for AZTEC
  int    *update_o,                  // vector elements updated on this node.
         *external_o;                // vector elements needed by this node.
  int    *update_index_o;            // ordering of update[] and external[]
  int    *extern_index_o;            // locally on this processor.
  //  int    *bindx;                 // Sparse matrix to be solved is stored
  //  double *val;                   // in these MSR arrays.
  int    N_update_o;                 // # of unknowns updated on this node

  AZ_set_proc_config(proc_config_o, AZ_NOT_MPI );

    AZ_read_update(&N_update_o, &update_o, proc_config_o, _dim_c, 1, AZ_linear);
    AZ_defaults(options_o,params_o);
    _dataAztec_o.aztecOptionsFromDataFile(options_o,params_o);
    AZ_transform(proc_config_o, &external_o,
	       (int *)_pattM.giveRaw_bindx(), _CDR.giveRaw_value(),
	       update_o, &update_index_o,
	       &extern_index_o, &data_org_o, N_update_o, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);

    chrono.start();
//    init_options_c(options_o,params_o);

    AZ_solve(_c.giveVec(),_f_c.giveVec(), options_o, params_o, NULL,
	   (int *)_pattM.giveRaw_bindx(), NULL, NULL, NULL,
	   _CDR.giveRaw_value(), data_org_o,status_o, proc_config_o);
  //
    chrono.stop();
    std::cout << "*** Solution (Concentration) computed in " << chrono.diff() << "s." << std::endl;
  _bdf.shift_right(_c);

}



template<typename Mesh> template<typename RegionMesh3D>
void ConvDiffReactSolverPC<Mesh>::
getvel(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u, const Real& time){

  //   for (UInt j=0; j<3; j++){
  //       for(UInt i=0; i< _dim_c; i++){
  //	  _u_c(i+j*_dim_c)=u(i+j*u.size()/3);
  //       }}

  for(UInt i=0; i < _mesh.numVertices(); i++){


    if(_u_to_c[i].ele == 0){
    // Dirichlet boundary for the velocity -> get velocity for boundary function
      Real xp=_mesh.point(i+1).x();
      Real yp=_mesh.point(i+1).y();
      Real zp=_mesh.point(i+1).z();
      for (ID jj=0; jj<3; ++jj)
	_u_c(i+jj*_u_c.size()/3)=BCh_u[(int)_u_to_c[i].b[0]](time,xp,yp,zp,BCh_u[(int)_u_to_c[i].b[0]].component(jj+1));

    }
    else{
    // Velocity interpolation
      _u_c(i)=_u_to_c[i].b[0]*u(umesh.volume(_u_to_c[i].ele).point(1).id()-1)+
                    _u_to_c[i].b[1]*u(umesh.volume(_u_to_c[i].ele).point(2).id()-1)+
                    _u_to_c[i].b[2]*u(umesh.volume(_u_to_c[i].ele).point(3).id()-1)+
	            _u_to_c[i].b[3]*u(umesh.volume(_u_to_c[i].ele).point(4).id()-1);

      _u_c(i+_dim_c)=_u_to_c[i].b[0]*u(umesh.volume(_u_to_c[i].ele).point(1).id()-1+u.size()/3)+
                                        _u_to_c[i].b[1]*u(umesh.volume(_u_to_c[i].ele).point(2).id()-1+u.size()/3)+
                                        _u_to_c[i].b[2]*u(umesh.volume(_u_to_c[i].ele).point(3).id()-1+u.size()/3)+
	                                _u_to_c[i].b[3]*u(umesh.volume(_u_to_c[i].ele).point(4).id()-1+u.size()/3);

      _u_c(i+2*_dim_c)=_u_to_c[i].b[0]*u(umesh.volume(_u_to_c[i].ele).point(1).id()-1+2*u.size()/3)+
                                          _u_to_c[i].b[1]*u(umesh.volume(_u_to_c[i].ele).point(2).id()-1+2*u.size()/3)+
                                          _u_to_c[i].b[2]*u(umesh.volume(_u_to_c[i].ele).point(3).id()-1+2*u.size()/3)+
	                                  _u_to_c[i].b[3]*u(umesh.volume(_u_to_c[i].ele).point(4).id()-1+2*u.size()/3);
    }


  }

}

template<typename Mesh> template<typename RegionMesh3D>
void ConvDiffReactSolverPC<Mesh>::
getcoord(RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BC_Handler& BCh_u){

  Real b1, b2, b3;
  intpolcoord localcoord;

  Chrono  chrono;
  chrono.start();

  _u_c=-100.0;

  Real x[4], y[4], z[4], xt[4], yt[4], zt[4];
  UInt vid, i1, i2, i3, v1, v2, v3, v4;
  LinearTetra ele;

  SimpleVect<GeoElement3D<LinearTetra> >::iterator iv = umesh.volumeList.begin();

  for(UInt i=0; i < _mesh.numVertices(); i++){

    x[0]=_mesh.point(i+1).x();
    y[0]=_mesh.point(i+1).y();
    z[0]=_mesh.point(i+1).z();

    if(_mesh.point(i+1).boundary()){

      UInt l;

      for(UInt k=0; k < BCh_u.size(); k++){
	if(BCh_u[k].flag() == _mesh.point(i+1).marker()){ l=k;}}

      switch( BCh_u[l].type() ) {
      case Essential:
	localcoord.b[0]=(Real)l;
	localcoord.b[1]=(Real)l;
	localcoord.b[2]=(Real)l;
	localcoord.b[3]=(Real)l;
	localcoord.ele=0;
	_u_to_c.push_back(localcoord);
	break;
      default:
	 for(Int found = 0;found == 0;){
	    vid = iv -> id();
	    Real volume[umesh.numLocalFaces()],minvolume=100.0;
	    UInt jk = UInt( -1 ); // initialized to max UInt+1 (trick for debugging if needed )
	    for(UInt jj=1; jj <= umesh.numLocalFaces();jj++){
	       i1=ele.fToP(jj,1);
	       i2=ele.fToP(jj,2);
	       i3=ele.fToP(jj,3);

	       i1=(iv->point(i1)).id();
	       i2=(iv->point(i2)).id();
	       i3=(iv->point(i3)).id();

	       x[1]=umesh.point(i1).x();
	       x[2]=umesh.point(i2).x();
	       x[3]=umesh.point(i3).x();
	       y[1]=umesh.point(i1).y();
	       y[2]=umesh.point(i2).y();
	       y[3]=umesh.point(i3).y();
	       z[1]=umesh.point(i1).z();
	       z[2]=umesh.point(i2).z();
	       z[3]=umesh.point(i3).z();
	       volume[jj]=calcvol(x,y,z);
	       if(volume[jj] < minvolume){
	       minvolume = volume[jj];
	       jk=jj;}}
	 if(minvolume < -0.0000000001){
	    if(umesh.faceList(umesh.localFaceId(vid,jk)).ad_first() == vid){
	       iv=iv+umesh.faceList(umesh.localFaceId(vid,jk)).ad_second()-vid;}
	    else{
	       iv=iv+umesh.faceList(umesh.localFaceId(vid,jk)).ad_first()-vid;}
	 }
	 else{
	    found = 1;
	    v1=(iv->point(1)).id();
	    v2=(iv->point(2)).id();
	    v3=(iv->point(3)).id();
	    v4=(iv->point(4)).id();
	    xt[0]=umesh.point(v1).x();
	    yt[0]=umesh.point(v1).y();
	    zt[0]=umesh.point(v1).z();
	    xt[1]=umesh.point(v2).x();
	    yt[1]=umesh.point(v2).y();
	    zt[1]=umesh.point(v2).z();
	    xt[2]=umesh.point(v3).x();
	    yt[2]=umesh.point(v3).y();
	    zt[2]=umesh.point(v3).z();
	    xt[3]=umesh.point(v4).x();
	    yt[3]=umesh.point(v4).y();
	    zt[3]=umesh.point(v4).z();

            if(test(xt, yt, zt, x[0], y[0], z[0], b1, b2, b3)){
	       localcoord.b[0]=1.0-b1-b2-b3;
	       localcoord.b[1]=b1;
	       localcoord.b[2]=b2;
	       localcoord.b[3]=b3;
	       localcoord.ele=vid;
	       _u_to_c.push_back(localcoord);}
	    else{
	       localcoord.b[0]=1.0-b1-b2-b3;
	       localcoord.b[1]=b1;
	       localcoord.b[2]=b2;
	       localcoord.b[3]=b3;
	       localcoord.ele=vid;
	       _u_to_c.push_back(localcoord);
	    }}}
	 break;
      }
    }
    else{
       for(Int found = 0;found == 0;){
	  vid = iv -> id();
	  Real volume[umesh.numLocalFaces()],minvolume=100.0;
	  UInt jk = UInt( -1 );
	  for(UInt jj=1; jj <= umesh.numLocalFaces();jj++){
	     i1=ele.fToP(jj,1);
	     i2=ele.fToP(jj,2);
	     i3=ele.fToP(jj,3);

	     i1=(iv->point(i1)).id();
	     i2=(iv->point(i2)).id();
	     i3=(iv->point(i3)).id();

	     x[1]=umesh.point(i1).x();
	     x[2]=umesh.point(i2).x();
	     x[3]=umesh.point(i3).x();
	     y[1]=umesh.point(i1).y();
	     y[2]=umesh.point(i2).y();
	     y[3]=umesh.point(i3).y();
	     z[1]=umesh.point(i1).z();
	     z[2]=umesh.point(i2).z();
	     z[3]=umesh.point(i3).z();
	     volume[jj]=calcvol(x,y,z);
	     if(volume[jj] < minvolume){
		minvolume = volume[jj];
		jk=jj;}}
	  if(minvolume < -0.00000001){
	     if(umesh.faceList(umesh.localFaceId(vid,jk)).ad_first() == vid){
		iv=iv+umesh.faceList(umesh.localFaceId(vid,jk)).ad_second()-vid;}
	     else{
		iv=iv+umesh.faceList(umesh.localFaceId(vid,jk)).ad_first()-vid;}
	  }
	  else{
	     found = 1;
	     v1=(iv->point(1)).id();
	     v2=(iv->point(2)).id();
	     v3=(iv->point(3)).id();
	     v4=(iv->point(4)).id();
	     xt[0]=umesh.point(v1).x();
	     yt[0]=umesh.point(v1).y();
	     zt[0]=umesh.point(v1).z();
	     xt[1]=umesh.point(v2).x();
	     yt[1]=umesh.point(v2).y();
	     zt[1]=umesh.point(v2).z();
	     xt[2]=umesh.point(v3).x();
	     yt[2]=umesh.point(v3).y();
	     zt[2]=umesh.point(v3).z();
	     xt[3]=umesh.point(v4).x();
	     yt[3]=umesh.point(v4).y();
	     zt[3]=umesh.point(v4).z();

	     if(test(xt, yt, zt, x[0], y[0], z[0], b1, b2, b3)){
		localcoord.b[0]=1.0-b1-b2-b3;
		localcoord.b[1]=b1;
		localcoord.b[2]=b2;
		localcoord.b[3]=b3;
	       localcoord.ele=vid;
	       _u_to_c.push_back(localcoord);}
	     else{
		localcoord.b[0]=1.0-b1-b2-b3;
		localcoord.b[1]=b1;
		localcoord.b[2]=b2;
		localcoord.b[3]=b3;
		localcoord.ele=vid;
		_u_to_c.push_back(localcoord);
	     }
	  }
       }
    }

  }

  chrono.stop();
  std::cout << " Calculation of the projection coordinates " << chrono.diff() << "s." << std::endl;

}

template<typename Mesh>
Real ConvDiffReactSolverPC<Mesh>::
calcvol(Real x[4], Real y[4], Real z[4]){
     Real volume=0.0;
     volume +=((x[1]-x[2])*(y[1]-y[3])*(z[1]-z[0]));
     volume +=((y[1]-y[2])*(z[1]-z[3])*(x[1]-x[0]));
     volume +=((z[1]-z[2])*(x[1]-x[3])*(y[1]-y[0]));
     volume -=((x[1]-x[0])*(y[1]-y[3])*(z[1]-z[2]));
     volume -=((y[1]-y[0])*(z[1]-z[3])*(x[1]-x[2]));
     volume -=((z[1]-z[0])*(x[1]-x[3])*(y[1]-y[2]));
     return volume;
}


template<typename Mesh>
int ConvDiffReactSolverPC<Mesh>::
test(Real x[4], Real y[4], Real z[4], Real & xp, Real & yp, Real & zp, Real & b1, Real & b2, Real & b3){

    Real a11, a12, a13, a21, a22, a23, a31, a32, a33, zw;
    int found;

    a11 = x[1]-x[0];
    a12 = x[2]-x[0];
    a13 = x[3]-x[0];
    a21 = y[1]-y[0];
    a22 = y[2]-y[0];
    a23 = y[3]-y[0];
    a31 = z[1]-z[0];
    a32 = z[2]-z[0];
    a33 = z[3]-z[0];

    b1  = xp-x[0];
    b2  = yp-y[0];
    b3  = zp-z[0];

//  Solve the equation system (without loop - faster)

        if(abs(a11) < max(abs(a21),abs(a31))){
          if(abs(a21) > abs(a31)){
	     zw  = a11;
	     a11 = a21;
	     a21 = zw;
	     zw  = a12;
	     a12 = a22;
	     a22 = zw;
	     zw  = a13;
	     a13 = a23;
	     a23 = zw;
	     zw  = b1;
	     b1  = b2;
	     b2  = zw;}
	  else{
	     zw  = a11;
	     a11 = a31;
	     a31 = zw;
	     zw  = a12;
	     a12 = a32;
	     a32 = zw;
	     zw  = a13;
	     a13 = a33;
	     a33 = zw;
	     zw  = b1;
	     b1  = b3;
	     b3  = zw;}}

        zw  = a21 / a11;
        a22 = a22 - zw * a12;
        a23 = a23 - zw * a13;
        b2  = b2  - zw * b1;
        zw  = a31 / a11;
        a32 = a32 - zw * a12;
        a33 = a33 - zw * a13;
        b3  = b3  - zw * b1;

        if(abs(a32) > abs(a22)){
	   zw  = a22;
	   a22 = a32;
	   a32 = zw;
	   zw  = a23;
	   a23 = a33;
	   a33 = zw;
	   zw  = b2;
	   b2  = b3;
	   b3  = zw;}

        zw = a32 / a22;
        a33 = a33 - zw * a23;

        b3 = b3 - zw * b2;
        b3 = b3 / a33;
        b2 = (b2 - a23 * b3) / a22;
        b1 = (b1 - a12 * b2 - a13 * b3) / a11;

	if((b1 >= 0.0) && (b2 >= 0.0) && (b3 >= 0.0) && (b1+b2+b3 <= 1.0))
	  found = 1;
	else
	  found = 0;

	return found;
}
}

#endif
