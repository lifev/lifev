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
#ifndef _DARCY_SOLVER_H
#define _DARCY_SOLVER_H
#include "darcyHandler.hpp"

/*!
  \brief A mixed hybrid Darcy solver
  \file darcySolver.h
  \author J.-F. Gerbeau and V. Martin
  \date 11/2002
*/

class DarcySolver:
  public DarcyHandler
{
public:
  DarcySolver(const GetPot& data_file);
  MSRPatt msrPattern;
  MSRMatr<double> mat;
  ScalUnknown<Vector> globalTP;  //!< Trace of Pressure (TP)
  ScalUnknown<Vector> globalF;   //!< Source term
  ScalUnknown<Vector> globalP;   //!< Pressure (P)
  ScalUnknown<Vector> globalFlux;//!< Flux (U)

  ElemVec elvecHyb;//!< Element vector for the Rhs (lives in RTkHyb(fe))
  ElemVec elvecSource;//!< Source element vector (lives in Qk(fe))
  ElemVec elvecFlux;  //!< Flux element vector (lives in RTk(fe))
  ElemMat elmatMix;//!< Mixed hybrid Matrix: [A B C] in [A B C; Bt 0 0; Ct 0 0]
  ElemMat elmatHyb;//!< Hybridization Matrix

  KNM<Real> BtB; //!< tmp array (NBP x NBP)
  KNM<Real> CtC; //!< tmp array (NBL x NBL): storage of the final hybrid array.
  KNM<Real> BtC; //!< tmp array (NBP x NBL)
  /*!
    signLocalFace(ivol,ilocface)
    =  1  if ivol is the first adjacent of face ilocface
    = -1                 second
    
    Remark: this array could be in regionMesh.
  */
  KNM<Real> signLocalFace; 

  SourceFct sourceFct;

  void computeHybridMatrix(); //!< compute the matrix for TP
  void applyBC(); //!< apply the b.c. for the TP problem
  void solveDarcy();//!< solve the linear system for TP
  void computePresFlux();//!< Compute P and U (once TP known)
  void postProcessPressureQ0();//!< postprocess P constant per element
  void postProcessPressureQ1();//!< projection of P (Q0) on Q1 and postproc.
  void postProcessVelocityQ1();//!< projection of U (RT0) on Q1 and postproc.
};

#endif
