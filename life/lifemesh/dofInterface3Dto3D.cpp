#include <cmath>
#include "dofInterface3Dto3D.hpp"


//! Constructor for interfacing Dof of the same type (RefFE)
/*!
  \param refFe the reference FE used in both meshes
  \param dof1 the Dof object of the mesh in which we want to make the computations
  \param dof2 the Dof object of the mesh which provides de data at the interface
*/
DofInterface3Dto3D::DofInterface3Dto3D(const RefFE& refFE, const Dof& dof1, 
	     const Dof& dof2):_refFE1(refFE),_dof1(dof1),_refFE2(refFE),_dof2(dof2),_dof(refFE) {}

//! Constructor for interfacing Dof of diferent type (RefFE)
/*!
  \param refFe1 the reference FE used in the mesh in which we want to make the computations
  \param dof1 the Dof object of the mesh in which we want to make the computations
  \param refFe2 the reference FE used in the mesh which provides de data at the interface
  \param dof2 the Dof object of the mesh which provides de data at the interface
*/
DofInterface3Dto3D::DofInterface3Dto3D(const RefFE& refFE1, const Dof& dof1, const RefFE refFE2,
	     const Dof& dof2):_refFE1(refFE1),_dof1(dof1),_refFE2(refFE2),_dof2(dof2),_dof(refFE1) {}


//! Returns true if the vectors v1 and v2 are equal with respect to the tolerance tol  
bool coincide(const KN_<Real>& v1, const KN_<Real>& v2, const Real& tol) {
 
  Real sumDif=0.;

  for (UInt i=0; i<nDimensions; ++i)
    sumDif += abs(v1[i]-v2[i]);

  if (sumDif <= tol )
    return true; 
  else 
    return false;
}

 
//! Returns true if points (x1,y1,z1) and (x2,y2,z2) are equal with respect to the tolerance tol 
bool coincide(const Real& x1, const Real& y1, const Real& z1, const Real& x2, const Real& y2, const Real& z2, const Real& tol) {

  Real sumDif = abs(x1-x2) +  abs(y1-y2) + abs(z1-z2);

  if (sumDif <= tol )
    return true; 
  else 
    return false;
}
