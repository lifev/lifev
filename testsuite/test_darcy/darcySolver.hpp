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
