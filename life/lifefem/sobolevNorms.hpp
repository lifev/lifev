//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief File containing procedures for computing norm and errors.

    @author

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#ifndef _SOBOLEVNORMS_H_INCLUDED
#define _SOBOLEVNORMS_H_INCLUDED

#include <boost/function.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/currentFE.hpp>

namespace LifeV
{

inline Real
elemL22( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
         const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    UInt ic;
    Real s = 0., f, x, y, z;
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ++ig )
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

//! returns the square of the L2 norm of u on the current element
// TO BE DELETED
template <typename VectorType, typename DofType>
Real
elem_L2_2( const VectorType & u, const CurrentFE& fe, const DofType& dof )
{
    Int inod;
    UInt eleID(fe.currentLocalId());
    Real s(0);
    Real uQuadPt;

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt = 0.;
        for ( UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
        {
            inod = dof.localToGlobal( eleID, iDof + 1 );
            uQuadPt += u( inod ) * fe.phi( iDof, iQuadPt );
        }
        s += uQuadPt * uQuadPt * fe.weightDet( iQuadPt );
    }
    return s;
}

//! version for vectorial problem
template <typename VectorType>
Real
elem_L2_2( const VectorType & u, const CurrentFE& fe, const Dof& dof,
           const UInt nbComp )
{
    Int dofID;
    UInt eleID (fe.currentLocalId());
    Real sum(0.0);
    Real uQuadPt;

    for ( UInt iComp (0); iComp < nbComp; ++iComp )
    {
        for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
        {
            uQuadPt = 0.;
            for ( UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
            {
                dofID = dof.localToGlobal( eleID, iDof + 1 ) + iComp * dof.numTotalDof();
                uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
            }
            sum += uQuadPt * uQuadPt * fe.weightDet( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the L2 norm of fct on the current element
inline Real
elem_L2_2( boost::function<Real( Real,Real,Real )> fct,
           const CurrentFE& fe )
{
    Real sum(0.0);
    Real f;
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt = 0; iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );
        f = fct( x, y, z );
        sum += f * f * fe.weightDet( iQuadPt );
    }
    return sum;
}

//! for time dependent+vectorial.
inline Real
elem_L2_2( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
           const CurrentFE& fe, const Real t, const UInt nbComp )
{
    Real sum(0.0);
    Real f;
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt = 0; iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );
        for ( UInt iComp(0); iComp < nbComp; ++iComp )
        {
            f = fct( t, x, y, z, iComp + 1 );
            sum += f * f * fe.weightDet( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the L2 norm of fct on the current element
// TO BE DELETED (duplicated)
inline Real
elem_f_L2_2( boost::function<Real( Real,Real,Real )> fct,
             const CurrentFE& fe )
{
    Real s = 0., f, x, y, z;
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        f = fct( x, y, z );
        s += f * f * fe.weightDet( ig );
    }
    return s;
}

//! for time dependent+vectorial.
// TO BE DELETED (duplicated)
inline Real
elem_f_L2_2( boost::function<Real( Real, Real, Real, Real, UInt )> fct,
             const CurrentFE& fe, const Real t, const UInt nbcomp )
{
    Real s = 0., f, x, y, z;
    for ( UInt ig = 0; ig < fe.nbQuadPt(); ig++ )
    {
        fe.coorQuadPt( x, y, z, ig );
        for ( UInt ic = 0; ic < nbcomp; ic++ )
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
elem_H1_2( const VectorType & u, const CurrentFE& fe, const Dof& dof, const int nbComp=1 )
{
    UInt eleID (fe.currentLocalId());
    Real sum (0.0);
    Real sum2(0.0);
    std::vector<Real> graduQuadPt(fe.nbCoor(),0.0);
    Real uQuadPt(0.0);

    for (UInt iComp (0); iComp < nbComp; ++iComp )
    {
        for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
        {

            uQuadPt = 0.0;
            for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
            {
                graduQuadPt[iCoor] =0.0;
            }

            for ( UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobal( eleID, iDof + 1 ) + iComp * dof.numTotalDof();
                uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
                for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
                {
                    graduQuadPt[iCoor] += u( dofID ) * fe.dphi( iDof, iCoor, iQuadPt );
                }
            }

            sum2 = uQuadPt * uQuadPt;
            for (UInt icoor = 0; icoor < fe.nbCoor(); icoor++)
            {
                sum2 += graduQuadPt[icoor] * graduQuadPt[icoor];
            }
            sum += sum2*fe.weightDet( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the H1 norm of fct on the current element
template<typename UsrFct>
Real
elem_H1_2( const UsrFct& fct, const CurrentFE& fe )
{
    Real sum(0.0);
    Real sum2(0.0);
    Real x;
    Real y;
    Real z;

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );
        sum2 = fct( x, y, z ) * fct( x, y, z );
        for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
        {
            sum2 += fct.grad(iCoor+1, x,y,z) *  fct.grad(iCoor+1, x,y,z);
        }
        sum += sum2 * fe.weightDet( iQuadPt );
    }
    return sum;
}

//! returns the square of the H1 norm of fct on the current element (time-dependent case)
template <typename UsrFct>
Real elem_H1_2( const UsrFct& fct, const CurrentFE& fe, const Real t, const UInt nbComp )
{
    Real sum(0.0);
    Real sum2(0.0);
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );

        for ( UInt iComp = 0; iComp < nbComp; ++iComp )
        {
            sum2 = fct(t, x, y, z, iComp+1 ) *  fct(t, x, y, z, iComp+1 );

            for (UInt iCoor = 0; iCoor < fe.nbCoor(); ++iCoor)
            {
                sum2 += fct.grad(iCoor+1, t,x,y,z, iComp+1) * fct.grad(iCoor+1, t,x,y,z, iComp+1);
            }
            sum += sum2 * fe.weightDet( iQuadPt );
        }
    }
    return sum;
}


//! returns the square of the L2 norm of (u-fct) on the current element
template <typename VectorType>
Real elem_L2_diff_2( VectorType & u,
                     boost::function<Real( Real, Real, Real )> fct,
                     const CurrentFE& fe,
                     const Dof& dof )
{
    UInt eleID (fe.currentLocalId());
    Real sum(0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt(0.0);
    Real diffQuadPt(0.0);

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt=0.0;
        for (UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobal( eleID, iDof + 1 );
            uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
        }
        fe.coorQuadPt( x, y, z, iQuadPt );
        diffQuadPt = uQuadPt - fct( x, y, z );
        sum += diffQuadPt * diffQuadPt * fe.weightDet( iQuadPt );
    }
    return sum;
}

//! returns the square of the L2 norm of (u-fct) on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elem_L2_diff_2( VectorType & u,
                     boost::function<Real( Real, Real, Real, Real, UInt )> fct,
                     const CurrentFE& fe,
                     const Dof& dof, const Real t, const UInt nbComp )
{
    // returns the square of the L2 norm of (u-fct) on the current element

    UInt eleID ( fe.currentLocalId() );
    Real sum(0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt;
    Real diffQuadPt;

    for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        for (UInt iComp = 0; iComp < nbComp; ++iComp )
        {
            uQuadPt=0.0;
            for (UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobal( eleID, iDof + 1 ) + iComp * dof.numTotalDof();
                uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
            }
            fe.coorQuadPt( x, y, z, iQuadPt );
            diffQuadPt = uQuadPt - fct( t, x, y, z, iComp + 1 );
            sum += diffQuadPt * diffQuadPt * fe.weightDet( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the H1 norm of (u-fct) on the current element
template <typename VectorType, typename UsrFct>
Real elem_H1_diff_2( const VectorType & u, const UsrFct& fct, const CurrentFE& fe,
                     const Dof& dof )
{
    UInt eleID (fe.currentLocalId());
    Real sum(0.0);
    Real sum2(0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt(0.0);
    Real diffQuadPt(0.0);

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt=0.0;
        Vector graduQuadPt = ZeroVector(fe.nbCoor());


        for (UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobal( eleID, iDof + 1 );
            uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
            for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
            {
                graduQuadPt(iCoor) += u( dofID ) * fe.dphi( iDof, iCoor, iQuadPt );
            }
        }

        fe.coorQuadPt( x, y, z, iQuadPt );
        diffQuadPt = uQuadPt - fct( x, y, z );

        Vector diffGradQuadPt = graduQuadPt;
        for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
        {
            diffGradQuadPt(iCoor) -= fct.grad(iCoor+1, x, y, z);
        }
        Real sum2 = diffQuadPt*diffQuadPt;
        for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
        {
            sum2 += diffGradQuadPt(iCoor) * diffGradQuadPt(iCoor);
        }
        sum += sum2* fe.weightDet( iQuadPt );
    }
    return sum;
}

//! returns the square of the H1 norm of (u-fct) on the current element  (time-dependent case)
template <typename VectorType, typename UsrFct>
Real elem_H1_diff_2( const VectorType & u, const UsrFct& fct, const CurrentFE& fe,
                     const Dof& dof, const Real t, const UInt nbComp )
{
    UInt eleID = fe.currentLocalId();
    Real sum(0.0);
    Real sum2(0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt(0.0);
    Real diffQuadPt(0.0);

    for (UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        for (UInt iComp(0); iComp < nbComp; ++iComp )
        {
            uQuadPt = 0.0;
            Vector graduQuadPt = ZeroVector(fe.nbCoor());

            for (UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobal( eleID, iDof + 1 ) + iComp * dof.numTotalDof();
                uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
                for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
                {
                    graduQuadPt(iCoor) += u( dofID ) * fe.dphi( iDof, iCoor, iQuadPt );
                }
            }

            fe.coorQuadPt( x, y, z, iQuadPt );

            diffQuadPt = uQuadPt - fct(t, x, y, z, iComp+1);

            Vector diffGradQuadPt = graduQuadPt;
            for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
            {
                diffGradQuadPt(iCoor) -= fct.grad(iCoor+1, t, x, y, z, iComp+1);
            }

            sum2 = diffQuadPt*diffQuadPt;
            for (UInt iCoor(0); iCoor < fe.nbCoor(); ++iCoor)
            {
                sum2 += diffGradQuadPt(iCoor) * diffGradQuadPt(iCoor);
            }
            sum += sum2* fe.weightDet( iQuadPt );
        }
    }
    return sum;
}

//! returns the integral of (u-fct) on the current element
//! for time dependent+vectorial
//! THIS IS WRONG if nbComp != 1 !!!
template <typename VectorType>
Real elem_integral_diff( VectorType & u,
                         boost::function<Real( Real, Real, Real, Real, UInt )> fct,
                         const CurrentFE& fe,
                         const Dof& dof, const Real t, const UInt nbComp )
{
    UInt eleID (fe.currentLocalId());
    Real sum(0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt;
    Real diffQuadPt;

    for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );
        uQuadPt =0.0;
        for ( UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobal( eleID, iDof + 1 )+ (nbComp-1) * dof.numTotalDof();
            uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
        }
        diffQuadPt = uQuadPt - fct( t, x, y, z, nbComp );
        sum += diffQuadPt * fe.weightDet( iQuadPt );
    }
    return sum;
}

//! returns the integral of u on the current element
//! for time dependent+vectorial
//! WRONG if nbComp is not 1 !!
template <typename VectorType>
Real elem_integral( VectorType & u,
                    const CurrentFE& fe,
                    const Dof& dof, const UInt nbComp )
{
    UInt eleID (fe.currentLocalId());
    Real sum(0.0);
    Real uQuadPt(0.0);

    for ( UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt = 0.0;

        for ( UInt iDof(0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobal( eleID, iDof + 1 ) + (nbComp-1) * dof.numTotalDof();
            uQuadPt += u( dofID ) * fe.phi( iDof, iQuadPt );
        }
        sum += uQuadPt * fe.weightDet( iQuadPt );
    }
    return sum;
}

//! returns the integral of fct on the current element
//! for time dependent+vectorial
//! WRONG if nbComp is not 1 !!
inline Real
elem_integral( boost::function<Real( Real, Real, Real,
                                     Real, UInt )> fct,
               const CurrentFE& fe, const Real t, const UInt nbComp )
{
    Real sum(0.0);
    Real x;
    Real y;
    Real z;

    for (UInt iQuadPt(0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        fe.coorQuadPt( x, y, z, iQuadPt );
        sum += fct( t, x, y, z, nbComp ) * fe.weightDet( iQuadPt );
    }
    return sum;
}

} // namespace LifeV
#endif
