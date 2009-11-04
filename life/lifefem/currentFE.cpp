/* -*- mode: c++ -*-
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
#include <life/lifefem/currentFE.hpp>

namespace LifeV
{
CurrentFE::CurrentFE( const RefFE& _refFE, const GeoMap& _geoMap, const QuadRule& _qr )
        :
        nbGeoNode( _geoMap.nbDof ),
        nbNode( _refFE.nbDof ),
        nbCoor( _refFE.nbCoor ),
        nbQuadPt( _qr.nbQuadPt ),
        nbDiag( _refFE.nbDiag() ),
        nbUpper( _refFE.nbUpper() ) ,
        nbPattern( _refFE.nbPattern() ),
        point( nbGeoNode, nbCoor ),
        refFE( _refFE ),
        geoMap( _geoMap ),
        qr( _qr ),
        phi( nbNode, nbQuadPt ),
        dPhiRef( nbNode, nbCoor, nbQuadPt ),
        dPhiRef2( nbNode, nbCoor, nbCoor, nbQuadPt ),
        phiDer( nbNode, nbCoor, nbQuadPt ),
        jacobian( nbCoor, nbCoor, nbQuadPt ),
        tInvJac( nbCoor, nbCoor, nbQuadPt ),
        phiGeo( nbGeoNode, nbQuadPt ),
        dPhiGeo( nbGeoNode, nbCoor, nbQuadPt ),
        weightDet( nbQuadPt ),
        detJac( nbQuadPt ),
        quadPt( nbQuadPt, 3 ),
        phiDer2( nbNode, nbCoor, nbCoor, nbQuadPt )
{
    CONSTRUCTOR( "CurrentFE" );
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int i = 0;i < nbNode;i++ )
        {
            phi( i, ig ) = refFE.phi( i, ig, qr );
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                dPhiRef( i, icoor, ig ) = refFE.dPhi( i, icoor, ig, qr );
                for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
                    dPhiRef2( i, icoor, jcoor, ig ) = refFE.d2Phi( i, icoor, jcoor, ig, qr );
            }
        }
        for ( int k = 0;k < nbGeoNode;k++ )
        {
            phiGeo( k, ig ) = geoMap.phi( k, ig, qr );
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                dPhiGeo( k, icoor, ig ) = geoMap.dPhi( k, icoor, ig, qr );
            }
        }
    }

#ifdef TEST_PRE
    _hasJac = false;
    _hasFirstDeriv = false;
    _hasSecondDeriv = false;
    _hasQuadPtCoor = false;
#endif
}
//----------------------------------------------------------------------
void CurrentFE::coorMap( Real& x, Real& y, Real& z,
                         const Real & xi, const Real & eta, const Real & zeta ) const
{
    x = y = z = 0.;

    /* former code using THREEDIM (commented out VM 07/04)
    // remove it permanently if the switch seems ok.
    for(int i=0;i<nbGeoNode;i++){
      x += point(i,0) * geoMap.phi(i,xi,eta,zeta);
      y += point(i,1) * geoMap.phi(i,xi,eta,zeta);
    #if defined(THREEDIM)
      z += point(i,2) * geoMap.phi(i,xi,eta,zeta);
    #endif
    }
    */

    //! I did NOT fully test this (VM)
    switch ( nbCoor )
    {
    case 1:    //! 1D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 ) * geoMap.phi( i, xi, eta, zeta );
        }
        break;
    case 2:    //! 2D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 ) * geoMap.phi( i, xi, eta, zeta );
            y += point( i, 1 ) * geoMap.phi( i, xi, eta, zeta );
        }
        break;
    case 3:    //! 3D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 ) * geoMap.phi( i, xi, eta, zeta );
            y += point( i, 1 ) * geoMap.phi( i, xi, eta, zeta );
            z += point( i, 2 ) * geoMap.phi( i, xi, eta, zeta );
        }
        break;
    default:
        ERROR_MSG( "Dimension (nbCoor): only 1, 2 or 3!" );
    }
}

  // Added by S. Quinodoz
void CurrentFE::coorBackMap(const Real& x, const Real& y, const Real& z,
			      Real & xi, Real & eta, Real& zeta) const
{
  // If the Mapping is not P1, this simple "back mapping" does not
  // work properly, as it suppose that the mapping is P1. If another
  // mapping (higher order, or for different shape like Q1) is needed
  // a Newton algorithm has to be implemented to solve the problem
  // M(x)-p = 0
  // with M the mapping, p the given point and x its coordinates in
  // the reference frame. However, this is more expensive and inexact.
  
  // ASSERT(geoMap.name.compare("P1")==0," The 'coorBackMap' function should NOT be used in this case! ");

  // In the P1 mapping case, we want to solve the problem 
  // P = \lamda_1 V_1 + \lambda_2 V_2 + \lambda_3 V_3
  // where V_i is the ith vertice of the cell because we have then
  // \hat{P} = \lambda_1 \hat{V}_1 + \lambda_2 \hat{V}_2 + \lambda_3 \hat{V}_3
  // where \hat{P} is our point in the reference frame
  //
  // For "simplicity", we use Cramer's formula (no need for an heavy solver)

  
  if (nbCoor ==1)
    {

      // 1D Case

      const Real a11=point(1,0)-point(0,0);
      const Real b1=x-point(0,0);
      
      Real lambda=b1/a11;
      
      xi  = lambda * 1; // 1= refFE.xi(1)-refFE.xi(0);
      eta = 0;
      zeta= 0;

    } else if ( nbCoor ==2)
    {
   
      // 2D Case
    
      const Real a11=point(1,0)-point(0,0);
      const Real a12=point(2,0)-point(0,0);
      const Real a21=point(1,1)-point(0,1);
      const Real a22=point(2,1)-point(0,1);
      const Real b1 =x-point(0,0);
      const Real b2 =y-point(0,1);
      
      Real D= a11*a22 - a21*a12;
      Real lambda_1= b1*a22 - b2*a12;
      Real lambda_2= a11*b2 - a21*b1;
      
      lambda_1/=D;
      lambda_2/=D;
      
      xi  = lambda_1 * 1 + lambda_2 * 0; // 1=refFE.xi(1)  -refFE.xi(0)  ; 0=refFE.xi(2)  -refFE.xi(0);
      eta = lambda_1 * 0 + lambda_2 * 1; // 0=refFE.eta(1) -refFE.eta(0) ; 1=refFE.eta(2) -refFE.eta(0);
      zeta= 0;
      
    }else if (nbCoor == 3)
    {
      
      // 3D Case
      
      const Real a11=point(1,0)-point(0,0);
      const Real a12=point(2,0)-point(0,0);
      const Real a13=point(3,0)-point(0,0);
      const Real a21=point(1,1)-point(0,1);
      const Real a22=point(2,1)-point(0,1);
      const Real a23=point(3,1)-point(0,1);
      const Real a31=point(1,2)-point(0,2);
      const Real a32=point(2,2)-point(0,2);
      const Real a33=point(3,2)-point(0,2);
      const Real b1 =x-point(0,0);
      const Real b2 =y-point(0,1);
      const Real b3 =z-point(0,2);
      
      const Real D= a11*a22*a33 + a31*a12*a23 + a21*a32*a13 - a11*a32*a23 - a31*a22*a13 - a21*a12*a33;
      Real lambda_1= b1*a22*a33 + b3*a12*a23 + b2*a32*a13 - b1*a32*a23 - b3*a22*a13 - b2*a12*a33;
      Real lambda_2= a11*b2*a33 + a31*b1*a23 + a21*b3*a13 - a11*b3*a23 - a31*b2*a13 - a21*b1*a33;
      Real lambda_3= a11*a22*b3 + a31*a12*b2 + a21*a32*b1 - a11*a32*b2 - a31*a22*b1 - a21*a12*b3;
      
      lambda_1/=D;
      lambda_2/=D;
      lambda_3/=D;
      
      xi  = lambda_1 * 1 + lambda_2 * 0 + lambda_3 * 0;
      eta = lambda_1 * 0 + lambda_2 * 1 + lambda_3 * 0;
      zeta= lambda_1 * 0 + lambda_2 * 0 + lambda_3 * 1;

    } else
    {   
      ERROR_MSG("Impossible dimension to invert coordinates");
    };
  
};

// Compute the Jacobian at the given point
// this means that jacobian(P1,P2,P3,i,j)
// computes d x_i / d zeta_j
// where x are the global coordinates
//       zeta the reference coordinates

Real CurrentFE::pointJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z, 
			      int comp_x, int comp_zeta) const
{
  Real jac(0);

  for ( int i = 0; i < nbGeoNode;i++ )
  {
    jac += point( i, comp_x ) * geoMap.dPhi( i , comp_zeta , hat_x, hat_y, hat_z );
  };
  
  return jac;
};

Real CurrentFE::pointInverseJacobian(const Real& hat_x, const Real& hat_y, const Real& hat_z,
				     int comp_x, int comp_zeta) const
{
  if ( nbCoor ==1 )
    {
 
      return 1/ pointJacobian(hat_x,hat_y,hat_z,comp_x,comp_zeta);
    } else if (nbCoor ==2)
    {
      
      Real a11= pointJacobian(hat_x,hat_y,hat_z,0,0);
      Real a12= pointJacobian(hat_x,hat_y,hat_z,0,1);
      Real a21= pointJacobian(hat_x,hat_y,hat_z,1,0);
      Real a22= pointJacobian(hat_x,hat_y,hat_z,1,1);

      int total(comp_x+comp_zeta);
      int mysign(1);
      if (total % 2 == 1)
	mysign=-1;

      return mysign*pointJacobian(hat_x,hat_y,hat_z,comp_zeta,comp_x)/(a11*a22-a12*a21);
      
    } else if (nbCoor ==3)
    {

      Real a11= pointJacobian(hat_x,hat_y,hat_z,0,0);
      Real a12= pointJacobian(hat_x,hat_y,hat_z,0,1);
      Real a13= pointJacobian(hat_x,hat_y,hat_z,0,2);
      Real a21= pointJacobian(hat_x,hat_y,hat_z,1,0);
      Real a22= pointJacobian(hat_x,hat_y,hat_z,1,1);
      Real a23= pointJacobian(hat_x,hat_y,hat_z,1,2);
      Real a31= pointJacobian(hat_x,hat_y,hat_z,2,0);
      Real a32= pointJacobian(hat_x,hat_y,hat_z,2,1);
      Real a33= pointJacobian(hat_x,hat_y,hat_z,2,2);
      
      //std::cout << a11 << "  " << a12 << "  " << a13 << std::endl;
      //std::cout << a21 << "  " << a22 << "  " << a23 << std::endl;
      //std::cout << a31 << "  " << a32 << "  " << a33 << std::endl;
      
      Real det= a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a13*a22*a31 - a12*a21*a33;
      
      std::vector<std::vector< Real > > cof(3,std::vector<Real>(3,0));
      
      cof[0][0]=  (a22*a33 - a23*a32);
      cof[0][1]= -(a21*a33 - a31*a23);
      cof[0][2]=  (a21*a32 - a31*a22);
      cof[1][0]= -(a12*a33 - a32*a13);
      cof[1][1]=  (a11*a33 - a13*a31);
      cof[1][2]= -(a11*a32 - a31*a12);
      cof[2][0]=  (a12*a23 - a13*a22);
      cof[2][1]= -(a11*a23 - a13*a21);
      cof[2][2]=  (a11*a22 - a12*a21);

      //std::cout << cof[0][0] << " ; " <<  cof[0][1] << " ; " <<  cof[0][2] << std::endl;
      //std::cout << cof[1][0] << " ; " <<  cof[1][1] << " ; " <<  cof[1][2] << std::endl;
      //std::cout << cof[2][0] << " ; " <<  cof[2][1] << " ; " <<  cof[2][2] << std::endl;
      //std::cout << det << std::endl;
     
      // inverse need the TRANSPOSE!
      return cof[comp_x][comp_zeta]/det;
    } else 
    {
        ERROR_MSG( "Dimension (nbCoor): only 1, 2 or 3!" );
    };
    return 0.;
};

//========Barycenter===============================
void CurrentFE::barycenter( Real& x, Real& y, Real& z )
{
    x = y = z = 0.;
    /*
    // former code using THREEDIM (commented out VM 07/04)
    // remove it permanently if the switch seems ok.
    for(int i=0;i<nbGeoNode;i++){
      x += point(i,0);
      y += point(i,1);
    #if defined(THREEDIM)
      z += point(i,2);
    #endif
    }
    */
    switch ( nbCoor )
    {
    case 1:    //! 1D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 );
        }
        break;
    case 2:    //! 2D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 );
            y += point( i, 1 );
        }
        break;
    case 3:    //! 3D
        for ( int i = 0;i < nbGeoNode;i++ )
        {
            x += point( i, 0 );
            y += point( i, 1 );
            z += point( i, 2 );
        }
        break;
    default:
        ERROR_MSG( "Dimension (nbCoor): only 1, 2 or 3!" );
    }

    x /= nbGeoNode;
    y /= nbGeoNode;
    z /= nbGeoNode;
}
//=====End of Barycenter=============================
//----------------------------------------------------------------------
Real CurrentFE::measure() const
{
    ASSERT_PRE( _hasJac, "The det of jacobian must have been computed" )
    Real meas = 0.;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
        meas += weightDet( ig );
    return meas;
}

//----------------------------------------------------------------------
Real CurrentFE::diameter() const
{
	int i, j, icoor;
    Real s, h = 0.;
    for ( i = 0;i < nbGeoNode - 1;i++ )
    {
        for ( j = i + 1;j < nbGeoNode;j++ )
        {
            s = 0.;
            for ( icoor = 0;icoor < (int)nDimensions; icoor++ )
            {
                s += fabs( point( i, icoor ) - point( j, icoor ) );
            }
            if ( s > h )
                h = s;
        }
    }
    return h;
}


//----------------------------------------------------------------------
void CurrentFE::_comp_quad_point_coor()
{
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        coorMap( quadPt( ig, 0 ), quadPt( ig, 1 ), quadPt( ig, 2 ),
                 qr.quadPointCoor( ig, 0 ), qr.quadPointCoor( ig, 1 ),
                 qr.quadPointCoor( ig, 2 ) );
    }
}

//----------------------------------------------------------------------
void CurrentFE::_comp_jacobian()
{
    Real fctDer;
    // derivatives of geo map:
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int icoor = 0;icoor < nbCoor;icoor++ )
        {
            for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
            {
                fctDer = 0.;
                for ( int j = 0;j < nbGeoNode;j++ )
                {
                    fctDer += point( j, icoor ) * dPhiGeo( j, jcoor, ig );
                }
                jacobian( icoor, jcoor, ig ) = fctDer;
            }
        }
    }
}
//----------------------------------------------------------------------
void CurrentFE::_comp_jacobian_and_det()
{
    //!first compute the jacobian matrix
    _comp_jacobian();

    /*
    // former code using THREEDIM (commented out VM 07/04)
    // remove it permanently if the switch seems ok.

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
    */
    Real a, b, c, d, e, f, g, h, i, ei, fh, bi, ch, bf, ce;
    // determinant on integrations points
    switch ( nbCoor )
    {
    case 1:    //! 1D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            detJac( ig ) = jacobian( 0, 0, ig );
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
        }
        break;
    case 2:    //! 2D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            a = jacobian( 0, 0, ig );
            b = jacobian( 0, 1, ig );
            c = jacobian( 1, 0, ig );
            d = jacobian( 1, 1, ig );
            detJac( ig ) = a * d - b * c;
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
        }
        break;
    case 3:    //! 3D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            a = jacobian( 0, 0, ig );
            b = jacobian( 0, 1, ig );
            c = jacobian( 0, 2, ig );
            d = jacobian( 1, 0, ig );
            e = jacobian( 1, 1, ig );
            f = jacobian( 1, 2, ig );
            g = jacobian( 2, 0, ig );
            h = jacobian( 2, 1, ig );
            i = jacobian( 2, 2, ig );
            ei = e * i;
            fh = f * h;
            bi = b * i;
            ch = c * h;
            bf = b * f;
            ce = c * e;
            detJac( ig ) = a * ( ei - fh ) + d * ( ch - bi ) + g * ( bf - ce );
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
        }
        break;
    default:
        ERROR_MSG( "Dimension (nbCoor): only 1, 2 or 3!" );
    }

}
//----------------------------------------------------------------------
void CurrentFE::_comp_inv_jacobian_and_det()
{
    //!first compute the jacobian matrix
    _comp_jacobian();

    /*
    // former code using THREEDIM (commented out VM 07/04)
    // remove it permanently if the switch seems ok.

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
    */

    Real a, b, c, d, e, f, g, h, i, ei, fh, bi, ch, bf, ce, det;
    //! determinant on integrations points an inverse transpose jacobian
    switch ( nbCoor )
    {
    case 1:    //! 1D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            det = jacobian( 0, 0, ig );
            detJac( ig ) = det;
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
            tInvJac( 0, 0, ig ) = 1 / det ;
        }
        break;
    case 2:    //! 2D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            a = jacobian( 0, 0, ig );
            b = jacobian( 0, 1, ig );
            c = jacobian( 1, 0, ig );
            d = jacobian( 1, 1, ig );
            det = a * d - b * c;
            detJac( ig ) = det;
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
            tInvJac( 0, 0, ig ) = d / det ;
            tInvJac( 0, 1, ig ) = -c / det ;
            tInvJac( 1, 0, ig ) = -b / det ;
            tInvJac( 1, 1, ig ) = a / det ;
        }
        break;
    case 3:    //! 3D
        for ( int ig = 0;ig < nbQuadPt;ig++ )
        {
            a = jacobian( 0, 0, ig );
            b = jacobian( 0, 1, ig );
            c = jacobian( 0, 2, ig );
            d = jacobian( 1, 0, ig );
            e = jacobian( 1, 1, ig );
            f = jacobian( 1, 2, ig );
            g = jacobian( 2, 0, ig );
            h = jacobian( 2, 1, ig );
            i = jacobian( 2, 2, ig );
            ei = e * i;
            fh = f * h;
            bi = b * i;
            ch = c * h;
            bf = b * f;
            ce = c * e;
            det = a * ( ei - fh ) + d * ( ch - bi ) + g * ( bf - ce );
            detJac( ig ) = det;
            weightDet( ig ) = detJac( ig ) * qr.weight( ig );
            tInvJac( 0, 0, ig ) = ( ei - fh ) / det ;
            tInvJac( 0, 1, ig ) = ( -d * i + f * g ) / det ;
            tInvJac( 0, 2, ig ) = ( d * h - e * g ) / det ;

            tInvJac( 1, 0, ig ) = ( -bi + ch ) / det ;
            tInvJac( 1, 1, ig ) = ( a * i - c * g ) / det ;
            tInvJac( 1, 2, ig ) = ( -a * h + b * g ) / det ;

            tInvJac( 2, 0, ig ) = ( bf - ce ) / det ;
            tInvJac( 2, 1, ig ) = ( -a * f + c * d ) / det ;
            tInvJac( 2, 2, ig ) = ( a * e - b * d ) / det ;
        }
        break;
    default:
        ERROR_MSG( "Dimension (nbCoor): only 1, 2 or 3!" );
    }
}

//----------------------------------------------------------------------
void CurrentFE::_comp_phiDer()
{
    Real x;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int j = 0;j < nbNode;j++ )
        {
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                x = 0.;
                for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
                {
                    x += tInvJac( icoor, jcoor, ig ) * dPhiRef( j, jcoor, ig ) ;
                }
                phiDer( j, icoor, ig ) = x;
            }
        }
    }
}

//----------------------------------------------------------------------
void CurrentFE::_comp_phiDer2()
{
    Real x;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int j = 0;j < nbNode;j++ )
        {
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
                {
                    x = 0.;
                    for ( int k1 = 0;k1 < nbCoor;k1++ )
                    {
                        for ( int k2 = 0;k2 < nbCoor;k2++ )
                        {
                            x += tInvJac( icoor, k1, ig ) * dPhiRef2( j, k1, k2, ig ) * tInvJac( jcoor, k2, ig );
                        }
                    }
                    phiDer2( j, icoor, jcoor, ig ) = x;
                }
            }
        }
    }
}
//----------------------------------------------------------------------
void CurrentFE::_comp_phiDerDer2()
{
    Real x1, x2;
    for ( int ig = 0;ig < nbQuadPt;ig++ )
    {
        for ( int j = 0;j < nbNode;j++ )
        {
            for ( int icoor = 0;icoor < nbCoor;icoor++ )
            {
                x1 = 0.;
                for ( int jcoor = 0;jcoor < nbCoor;jcoor++ )
                {
            x1 += tInvJac(icoor,jcoor,ig)*dPhiRef(j,jcoor,ig);
                    x2 = 0.;
                    for ( int k1 = 0;k1 < nbCoor;k1++ )
                    {
                        for ( int k2 = 0;k2 < nbCoor;k2++ )
                        {
                            x2 += tInvJac( icoor, k1, ig ) * dPhiRef2( j, k1, k2, ig ) * tInvJac( jcoor, k2, ig );
                        }
                    }
                    phiDer2( j, icoor, jcoor, ig ) = x2;
                }
                phiDer( j, icoor, ig ) = x1;
            }
        }
    }
}
}
