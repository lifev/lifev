#include "currentFEDG.hpp"

namespace LifeV
{

CurrentFEDG::CurrentFEDG(const RefFEDG& _refFE, const GeoMap& _geoMap, const QuadRule& _qr):
  nbGeoNode(_geoMap.nbDof), nbNode(_refFE.nbDof), nbCoor(_refFE.nbCoor),
  nbQuadPt(_qr.nbQuadPt), 
  nbDiag( _refFE.elPattern.nbDiag() ), nbUpper( _refFE.elPattern.nbUpper() ) , nbPattern( _refFE.elPattern.nbPattern() ),
  point(nbGeoNode,nbCoor),
  refFE(_refFE), geoMap(_geoMap), qr(_qr),
  phi(nbNode,nbQuadPt), dPhiRef(nbNode,nbCoor,nbQuadPt), dPhiRef2(nbNode,nbCoor,nbCoor,nbQuadPt),
  phiDer(nbNode,nbCoor,nbQuadPt), jacobian(nbCoor,nbCoor,nbQuadPt), tInvJac(nbCoor,nbCoor,nbQuadPt),
  phiGeo(nbGeoNode,nbQuadPt), dPhiGeo(nbGeoNode,nbCoor,nbQuadPt),
  weightDet(nbQuadPt), detJac(nbQuadPt),quadPt(nbQuadPt,3),
  mass(nbNode, nbNode),
  invMass(nbNode, nbNode),
  phiDer2(nbNode,nbCoor,nbCoor,nbQuadPt)
{
  CONSTRUCTOR("CurrentFEDG");

  for(int ig = 0; ig < nbQuadPt; ig++){
    for(int i = 0; i < nbNode; i++){
      phi(i, ig) = refFE.phi(i, ig, qr);
      for(int icoor = 0; icoor < nbCoor; icoor++){
	dPhiRef(i, icoor, ig) = refFE.dPhi(i, icoor, ig, qr);
        for(int jcoor = 0; jcoor < nbCoor; jcoor++)
          dPhiRef2(i, icoor, jcoor, ig) = refFE.d2Phi(i, icoor, jcoor, ig, qr); 
      } // for icoor
    } // for i

    for(int k = 0; k < nbGeoNode; k++){
      phiGeo(k, ig) = geoMap.phi(k, ig, qr);
      for(int icoor = 0; icoor < nbCoor; icoor++){
	dPhiGeo(k, icoor, ig) = geoMap.dPhi(k, icoor, ig, qr);
      } // for icoor
    } // for k
  }// for ig

#ifdef TEST_PRE
  _hasJac = false;
  _hasFirstDeriv = false;
  _hasSecondDeriv = false;
  _hasQuadPtCoor = false;
#endif
}

//==============================================================================

void CurrentFEDG::coorMap(Real& x,Real& y,Real& z,
			const Real & xi,const Real & eta, const Real & zeta) const
{
  x = y = z = 0.;
  for(int i = 0; i < nbGeoNode; i++){
    x += point(i,0) * geoMap.phi(i, xi, eta, zeta);
    y += point(i,1) * geoMap.phi(i, xi, eta, zeta);
#if defined(THREEDIM)
    z += point(i,2) * geoMap.phi(i, xi, eta, zeta);
#endif
  }
}

//==============================================================================

Real CurrentFEDG::measure() const
{
  ASSERT_PRE(_hasJac,"The det of jacobian must have been computed")
  Real meas = 0.;
  for(int ig = 0; ig < nbQuadPt; ig++) meas += weightDet(ig);
  return meas;
}

Real CurrentFEDG::diameter() const
{
  int i, j, icoor;
  Real s, h = 0.;
  for(i = 0; i < nbGeoNode-1; i++){
    for(j = i + 1; j < nbGeoNode; j++){
      s = 0.;
      for(icoor = 0; icoor < 3; icoor++){
	s += fabs(point(i, icoor) - point(j, icoor));
      }
      if(s > h) h = s;
    }
  }
  return h;
}

//==============================================================================

void CurrentFEDG::_comp_quad_point_coor()
{
  for(int ig = 0;ig < nbQuadPt;ig++){
    coorMap(quadPt(ig,0), quadPt(ig,1), quadPt(ig,2),
	    qr.quadPointCoor(ig,0), qr.quadPointCoor(ig,1), qr.quadPointCoor(ig,2));
  }
}

//==============================================================================

void CurrentFEDG::_comp_jacobian()
{
  Real fctDer;

  // derivatives of geo map:

  for(int ig = 0; ig < nbQuadPt; ig++){
    for(int icoor = 0; icoor < nbCoor; icoor++){
      for(int jcoor = 0; jcoor < nbCoor; jcoor++){
	fctDer = 0.;
	for(int j = 0; j < nbGeoNode; j++){
	  fctDer += point(j, icoor) * dPhiGeo(j, jcoor, ig);
	}
	jacobian(icoor, jcoor, ig) = fctDer;
      }
    }
  }

  // determinant on integrations points

#if defined(TWODIM)
  // *** 2D code *** 
  Real a,b,c,d;
  for(int ig=0;ig<nbQuadPt;ig++){
    a = jacobian(0,0,ig);
    b = jacobian(0,1,ig);
    c = jacobian(1,0,ig);
    d = jacobian(1,1,ig);
    detJac(ig) = a*d - b*c;
    weightDet(ig) = detJac(ig) * qr.weight(ig);
  }
#elif defined(THREEDIM)
  // *** 3D code *** 
  Real a,b,c,d,e,f,g,h,i,ei,fh,bi,ch,bf,ce;
  for(int ig=0;ig<nbQuadPt;ig++){
    a = jacobian(0,0,ig);
    b = jacobian(0,1,ig);
    c = jacobian(0,2,ig);
    d = jacobian(1,0,ig);
    e = jacobian(1,1,ig);
    f = jacobian(1,2,ig);
    g = jacobian(2,0,ig);
    h = jacobian(2,1,ig);
    i = jacobian(2,2,ig);
    ei=e*i;
    fh=f*h;
    bi=b*i;
    ch=c*h;
    bf=b*f;
    ce=c*e;
    detJac(ig) = a*(ei-fh) + d*(ch-bi) + g*(bf-ce);
    weightDet(ig) = detJac(ig) * qr.weight(ig);
  }
#endif
}

//==============================================================================

void CurrentFEDG::_comp_inv_jacobian()
{
    Real fctDer;
    // derivatives of geo map:
    for(int ig=0;ig<nbQuadPt;ig++){
      for(int icoor=0;icoor<nbCoor;icoor++){
	for(int jcoor=0;jcoor<nbCoor;jcoor++){
	  fctDer = 0.;
	  for(int j=0;j<nbGeoNode;j++){
	    fctDer += point(j,icoor)*dPhiGeo(j,jcoor,ig);
	  }
	  jacobian(icoor,jcoor,ig) = fctDer;
	}
      }
    }

    // determinant on integrations points an inverse tranpose jacobian

#if defined(TWODIM)
    // *** 2D code *** 
    Real a,b,c,d,det;
    for(int ig=0;ig<nbQuadPt;ig++){
      a = jacobian(0,0,ig);
      b = jacobian(0,1,ig);
      c = jacobian(1,0,ig);
      d = jacobian(1,1,ig);
      det = a*d - b*c;
      detJac(ig) = det;
      weightDet(ig) = detJac(ig) * qr.weight(ig);
      tInvJac(0,0,ig) = d/det ;
      tInvJac(0,1,ig) =-c/det ;   
      tInvJac(1,0,ig) =-b/det ;
      tInvJac(1,1,ig) = a/det ;
    }
#elif defined(THREEDIM)
    // *** 3D code *** 
    Real a,b,c,d,e,f,g,h,i,ei,fh,bi,ch,bf,ce,det;
    for(int ig=0;ig<nbQuadPt;ig++){
      a = jacobian(0,0,ig);
      b = jacobian(0,1,ig);
      c = jacobian(0,2,ig);
      d = jacobian(1,0,ig);
      e = jacobian(1,1,ig);
      f = jacobian(1,2,ig);
      g = jacobian(2,0,ig);
      h = jacobian(2,1,ig);
      i = jacobian(2,2,ig);
      ei=e*i;
      fh=f*h;
      bi=b*i;
      ch=c*h;
      bf=b*f;
      ce=c*e;
      det = a*(ei-fh) + d*(ch-bi) + g*(bf-ce);
      detJac(ig) = det;
      weightDet(ig) = detJac(ig) * qr.weight(ig);
      tInvJac(0,0,ig) = (  ei - fh )/det ;
      tInvJac(0,1,ig) = (-d*i + f*g)/det ;
      tInvJac(0,2,ig) = ( d*h - e*g)/det ;
      
      tInvJac(1,0,ig) = ( -bi + ch )/det ;
      tInvJac(1,1,ig) = ( a*i - c*g)/det ;
      tInvJac(1,2,ig) = (-a*h + b*g)/det ;
      
      tInvJac(2,0,ig) = (  bf - ce )/det ;
      tInvJac(2,1,ig) = (-a*f + c*d)/det ;
      tInvJac(2,2,ig) = ( a*e - b*d)/det ;
    }
#endif
}

//==============================================================================

void CurrentFEDG::_comp_mass()
{
#ifdef TEST_PRE
  ASSERT_PRE(_hasJac, "Mass matrix computation needs updated Jacobian")
#endif

    Real s;

    for(int i = 0; i < nbNode; i++){
      for(int j = 0; j < nbNode; j++){
	s = 0.;

	for(int ig = 0; ig < nbQuadPt; ig++){
	  s += phi(i, ig) * phi(j, ig) * weightDet(ig);
	}
	mass(i, j) = s;
      }
    }
}

//==============================================================================

void CurrentFEDG::_comp_inv_mass()
{
#ifdef TEST_PRE
  ASSERT_PRE(_hasJac, "Mass matrix computation needs updated Jacobian")
#endif
    Real s;
  
  // Compute mass matrix

    for(int i = 0; i < nbNode; i++){
      for(int j = 0; j < nbNode; j++){
	s = 0.;
	for(int ig = 0; ig < nbQuadPt; ig++){
	  s += phi(i, ig) * phi(j, ig) * weightDet(ig);
	}

	mass(i, j) = s;
      }
    }

    // Leverrier(mass, invMass);
}

//==============================================================================
// Useful stuff
//==============================================================================

template<typename P>
void Leverrier(KNM<P>& A, KNM<P>& invA)
{
  // A is assumed to be square matrices of size A.N()

  ASSERT_PRE((A.N() == A.M()), "Matrix must be square")

  const int n = A.N();

  P alpha_k = 0.;

  KNM<Real> B_k_minus_1(n, n);
  KNM<Real> B_k(n, n);
  KNM<Real> Temp(n, n);

  // B_0 = I
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      if(i == j){
	B_k_minus_1(i, j) = 1.;
      }else{
	B_k_minus_1(i, j) = 0.;
      }
    }
  }

  for(int k = 1; k < n; k++){
    MatMult(A, B_k_minus_1, Temp);

    alpha_k = Trace(Temp) / k;

    // B_k = - A * B_k_minus_1 + alpha_k * I
    for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
	B_k(i, j) = - Temp(i, j);
	if(i == j){
	  B_k(i, j) += alpha_k;
	}
      }
    }
    B_k_minus_1 = B_k;
  }

  MatMult(A, B_k, Temp);
  P alpha_n = Trace(Temp) / n;

  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      invA(i, j) = B_k(i, j)/alpha_n;
    }
  }

}

// Multiplies A by B and stores everything in C
template<typename MATRIX>
void MatMult(MATRIX& A, MATRIX& B, MATRIX& C){
  ASSERT_PRE((A.M() == B.N()), "Columns of first matrix must equal rows of second matrix")

    C = 0.;

    for(int i = 0; i < A.N(); i++){
      for(int j = 0; j < B.M(); j++){
	for(int k = 0; k < A.M(); k++){
	  C(i, j) += A(i, k) * B(k, j);
	}
      }
    }
}

// Returns the trace of matrix A
template<typename P>
P Trace(KNM<P>& A){
  ASSERT_PRE(A.N() == A.M(), "Can only compute the trace of a square matrix")

  P s = 0.;

  for(int i = 0; i < A.N(); i++)
    s += A(i, i);

  return s;
}
}
