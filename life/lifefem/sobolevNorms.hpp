/*
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
#ifndef _SOBOLEVNORMS_H_INCLUDED
#define _SOBOLEVNORMS_H_INCLUDED

#include <boost/function.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/currentFE.hpp>

namespace LifeV
{
//! returns the square of the L2 norm of u on the current element
template <typename VectorType, typename DofType>
Real
elem_L2_2( const VectorType & u, const CurrentFE& fe, const DofType& dof )
{
    int i, inod, ig;
    UInt eleID = fe.currentLocalId();
    Real s = 0, u_ig;
    for ( ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        u_ig = 0.;
        for ( i = 0;i < fe.nbFENode();i++ )
        {
            inod = dof.localToGlobal( eleID, i + 1 );
            u_ig += u( inod ) * fe.phi( i, ig );
        }
        s += u_ig * u_ig * fe.weightDet( ig );
    }
    return s;
}

//! version for vectorial problem
template <typename VectorType>
Real
elem_L2_2( const VectorType & u, const CurrentFE& fe, const Dof& dof,
           const int nbcomp )
{
    int i, inod, ig;
    UInt eleID = fe.currentLocalId();
    UInt ic;
    Real s = 0, u_ig;
    for ( ic = 0; ic < (UInt)nbcomp; ic++ )
    {
        for ( ig = 0;ig < fe.nbQuadPt();ig++ )
        {
            u_ig = 0.;
            for ( i = 0;i < fe.nbFENode();i++ )
            {
                inod = dof.localToGlobal( eleID, i + 1 ) + ic * dof.numTotalDof();
                u_ig += u( inod ) * fe.phi( i, ig );
            }
            s += u_ig * u_ig * fe.weightDet( ig );
        }
    }
    return s;
}

//! returns the square of the L2 norm of fct on the current element
inline Real
elem_L2_2( boost::function<Real( Real,Real,Real )> fct,
           const CurrentFE& fe )
{
    Real s = 0., f, x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        f = fct( x, y, z );
        s += f * f * fe.weightDet( ig );
    }
    return s;
}

//! for time dependent+vectorial.
inline Real
elem_L2_2( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
           const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    int ig;
    UInt ic;
    Real s = 0., f, x, y, z;
    for ( ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( ic = 0; ic < nbcomp; ic++ )
        {
            f = fct( t, x, y, z, ic + 1 );
            s += f * f * fe.weightDet( ig );
        }
    }
    return s;
}

//! returns the square of the L2 norm of fct on the current element
inline Real
elem_f_L2_2( boost::function<Real( Real,Real,Real )> fct,
           const CurrentFE& fe )
{
    Real s = 0., f, x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        f = fct( x, y, z );
        s += f * f * fe.weightDet( ig );
    }
    return s;
}

//! for time dependent+vectorial.
inline Real
elem_f_L2_2( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
         const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    int ig;
    UInt ic;
    Real s = 0., f, x, y, z;
    for ( ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( ic = 0; ic < nbcomp; ic++ )
        {
            f = fct( t, x, y, z, ic + 1 );
            s += f * f * fe.weightDet( ig );
        }
    }
    return s;
}

//! returns the square of the H1 norm of u on the current element

template <typename VectorType>
Real
elem_H1_2( const VectorType & u, const CurrentFE& fe, const Dof& dof, const int nbcomp=1 )
{
    UInt eleID = fe.currentLocalId();
    Real s = 0;
    for (int ic = 0; ic < nbcomp; ic++ )
    {
    	for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
		{
    		Real u_ig(0.);
    		Vector grad_u_ig = ZeroVector(fe.nbCoor());
    		for ( int i = 0;i < fe.nbFENode();i++ )
    		{
    			int inod = dof.localToGlobal( eleID, i + 1 ) + ic * dof.numTotalDof();
    			u_ig += u( inod ) * fe.phi( i, ig );
    			for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
    				grad_u_ig(icoor) += u( inod ) * fe.phiDer( i, icoor, ig );
    		}
    		Real s_tmp = u_ig * u_ig;
    		for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
    			s_tmp += pow(grad_u_ig(icoor),2);
    		s += s_tmp*fe.weightDet( ig );
		}
    }
    return s;
}

//! returns the square of the H1 norm of fct on the current element
template<typename UsrFct>
Real
elem_H1_2( const UsrFct& fct, const CurrentFE& fe )
{
    Real s(0), x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        Real s_tmp = pow(fct( x, y, z ),2);
        for(int icoor = 0; icoor < fe.nbCoor(); icoor++)
        	s_tmp += pow(fct.grad(icoor+1, x,y,z),2);
        s += s_tmp * fe.weightDet( ig );
    }
    return s;
}

//! returns the square of the H1 norm of fct on the current element (time-dependent case)
template <typename UsrFct>
Real elem_H1_2( const UsrFct& fct, const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    Real s(0), x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( UInt ic = 0;ic < nbcomp;ic++ )
        {
            Real s_tmp = pow(fct(t, x, y, z, ic+1 ),2);
            for(int icoor = 0; icoor < fe.nbCoor(); icoor++)
            	s_tmp += pow(fct.grad(icoor+1, t,x,y,z, ic+1),2);
            s += s_tmp * fe.weightDet( ig );
        }
    }
    return s;
}


//! returns the square of the L2 norm of (u-fct) on the current element
template <typename VectorType>
Real elem_L2_diff_2( VectorType & u,
                     boost::function<Real( Real, Real, Real )> fct,
                     const CurrentFE& fe,
                     const Dof& dof )
{
    UInt eleID = fe.currentLocalId();
    Real s(0), x, y, z;
    for (int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        Real u_ig(0);
        for (int i = 0;i < fe.nbFENode();i++ )
        {
            int inod = dof.localToGlobal( eleID, i + 1 );
            u_ig += u( inod ) * fe.phi( i, ig );
        }
        fe.coorQuadPt( x, y, z, ig );
        Real diff_ig = u_ig - fct( x, y, z );
        s += diff_ig * diff_ig * fe.weightDet( ig );
    }
    return s;
}

//! returns the square of the L2 norm of (u-fct) on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elem_L2_diff_2( VectorType & u,
                     boost::function<Real( Real, Real, Real, Real, UInt )> fct,
                     const CurrentFE& fe,
                     const Dof& dof, const Real t, const int nbcomp )
{
    // returns the square of the L2 norm of (u-fct) on the current element
    UInt eleID = fe.currentLocalId();
    Real s(0), x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        for (int ic = 0; ic < nbcomp; ic++ )
        {
            Real u_ig(0);
            for (int i = 0;i < fe.nbFENode();i++ )
            {
                int inod = dof.localToGlobal( eleID, i + 1 ) + ic * dof.numTotalDof();
                u_ig += u( inod ) * fe.phi( i, ig );
            }
            fe.coorQuadPt( x, y, z, ig );
            Real diff_ig = u_ig - fct( t, x, y, z, ic + 1 );
            s += diff_ig * diff_ig * fe.weightDet( ig );
        }
    }
    return s;
}

//! returns the square of the H1 norm of (u-fct) on the current element
template <typename VectorType, typename UsrFct>
Real elem_H1_diff_2( const VectorType & u, const UsrFct& fct, const CurrentFE& fe,
                     const Dof& dof )
{
    UInt eleID = fe.currentLocalId();
    Real s(0), x, y, z;
    for (int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        Real u_ig = 0.;
        Vector grad_u_ig = ZeroVector(fe.nbCoor());

        for (int i = 0;i < fe.nbFENode();i++ )
        {
            int inod = dof.localToGlobal( eleID, i + 1 );
            u_ig += u( inod ) * fe.phi( i, ig );
            for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
            	grad_u_ig(icoor) += u( inod ) * fe.phiDer( i, icoor, ig );
        }
        fe.coorQuadPt( x, y, z, ig );
        Real diff_ig = u_ig - fct( x, y, z );
        Vector grad_diff_ig = grad_u_ig;
        for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
        	grad_diff_ig(icoor) -= fct.grad(icoor+1, x, y, z);
        Real s_tmp = diff_ig*diff_ig;
        for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
        	s_tmp += pow(grad_diff_ig(icoor),2);
        s += s_tmp* fe.weightDet( ig );
    }
    return s;
}

//! returns the square of the H1 norm of (u-fct) on the current element  (time-dependent case)
template <typename VectorType, typename UsrFct>
Real elem_H1_diff_2( const VectorType & u, const UsrFct& fct, const CurrentFE& fe,
                     const Dof& dof, const Real t, const UInt nbcomp )
{
	UInt eleID = fe.currentLocalId();
	Real s(0), x, y, z;
	for (int ig = 0;ig < fe.nbQuadPt();ig++ )
	{
        for (UInt ic = 0; ic < nbcomp; ic++ )
        {
        	Real u_ig = 0.;
        	Vector grad_u_ig = ZeroVector(fe.nbCoor());

        	for (int i = 0;i < fe.nbFENode();i++ )
        	{
        		int inod = dof.localToGlobal( eleID, i + 1 ) + ic * dof.numTotalDof();
        		u_ig += u( inod ) * fe.phi( i, ig );
        		for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
        			grad_u_ig(icoor) += u( inod ) * fe.phiDer( i, icoor, ig );
        	}
        	fe.coorQuadPt( x, y, z, ig );
            Real diff_ig = u_ig - fct(t, x, y, z, ic+1);
            Vector grad_diff_ig = grad_u_ig;
            for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
            	grad_diff_ig(icoor) -= fct.grad(icoor+1, t, x, y, z, ic+1);
            Real s_tmp = diff_ig*diff_ig;
            for (int icoor = 0; icoor < fe.nbCoor(); icoor++)
            	s_tmp += pow(grad_diff_ig(icoor),2);
            s += s_tmp* fe.weightDet( ig );
        }
	}
    return s;
}

//! returns the integral of (u-fct) on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elem_integral_diff( VectorType & u,
                         boost::function<Real( Real, Real, Real, Real, UInt )> fct,
                         const CurrentFE& fe,
                         const Dof& dof, const Real t, const int nbcomp )
{
    UInt eleID = fe.currentLocalId();
    Real s(0), x, y, z;
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        Real u_ig(0);
        for ( int i = 0;i < fe.nbFENode();i++ )
        {
            int inod = dof.localToGlobal( eleID, i + 1 )+ (nbcomp-1) * dof.numTotalDof();
            u_ig += u( inod ) * fe.phi( i, ig );
        }
        Real diff_ig = u_ig - fct( t, x, y, z, nbcomp );
        s += diff_ig * fe.weightDet( ig );
    }
    return s;
}

//! returns the integral of u on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elem_integral( VectorType & u,
                    const CurrentFE& fe,
                    const Dof& dof, const int nbcomp )
{
    UInt eleID = fe.currentLocalId();
    Real s(0);
    for ( int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        Real u_ig(0);
        for ( int i = 0;i < fe.nbFENode();i++ )
        {
            int inod = dof.localToGlobal( eleID, i + 1 ) + (nbcomp-1) * dof.numTotalDof();
            u_ig += u( inod ) * fe.phi( i, ig );
        }
        s += u_ig * fe.weightDet( ig );
    }
    return s;
}

//! returns the integral of fct on the current element
//! for time dependent+vectorial
inline Real
elem_integral( boost::function<Real( Real, Real, Real,
                                       Real, UInt )> fct,
               const CurrentFE& fe, const Real t, const int nbcomp )
{
    Real s(0), x, y, z;
    for (int ig = 0;ig < fe.nbQuadPt();ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        s += fct( t, x, y, z, nbcomp ) * fe.weightDet( ig );
    }
    return s;
}

} // namespace LifeV
#endif
