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
    @contributor Mauro Perego <perego.mauro@gmail.com>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

 */

#ifndef _SOBOLEVNORMS_H_INCLUDED
#define _SOBOLEVNORMS_H_INCLUDED

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/fem/CurrentFE.hpp>

namespace LifeV
{
//! @name Public typedefs
//@{
typedef boost::numeric::ublas::vector<Real> Vector;
typedef boost::numeric::ublas::zero_vector<Real> ZeroVector;
//@}

//! version for vectorial problem
template <typename VectorType>
Real
elementaryL2NormSquare ( const VectorType& u, const CurrentFE& fe, const DOF& dof,
                         const UInt nbComp )
{
    Int dofID;
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real uQuadPt;

    for ( UInt iComp (0); iComp < nbComp; ++iComp )
    {
        for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
        {
            uQuadPt = 0.;
            for ( UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                dofID = dof.localToGlobalMap ( eleID, iDof ) + iComp * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
            }
            sum += uQuadPt * uQuadPt * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the L2 norm of fct on the current element
inline Real
elementaryFctL2NormSquare ( boost::function<Real ( Real, Real, Real ) > fct,
                            const CurrentFE& fe )
{
    Real sum (0.0);
    Real f;
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt = 0; iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);
        f = fct ( x, y, z );
        sum += f * f * fe.weightDet ( iQuadPt );
    }
    return sum;
}


//! for time dependent+vectorial.
inline Real
elementaryFctL2NormSquare ( boost::function<Real ( Real, Real, Real, Real, UInt ) > fct,
                            const CurrentFE& fe, const Real t, const UInt nbComp )
{
    Real sum (0.0);
    Real f;
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt = 0; iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);
        for ( UInt iComp (0); iComp < nbComp; ++iComp )
        {
            f = fct ( t, x, y, z, iComp );
            sum += f * f * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the H1 norm of u on the current element
template <typename VectorType>
Real
elementaryH1NormSquare ( const VectorType& u, const CurrentFE& fe, const DOF& dof, const UInt nbComp = 1 )
{
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real sum2 (0.0);
    std::vector<Real> graduQuadPt (fe.nbCoor(), 0.0);
    Real uQuadPt (0.0);

    for (UInt iComp (0); iComp < nbComp; ++iComp )
    {
        for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
        {

            uQuadPt = 0.0;
            for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
            {
                graduQuadPt[iCoor] = 0.0;
            }

            for ( UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobalMap ( eleID, iDof ) + iComp * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
                for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
                {
                    graduQuadPt[iCoor] += u ( dofID ) * fe.dphi ( iDof, iCoor, iQuadPt );
                }
            }

            sum2 = uQuadPt * uQuadPt;
            for (UInt icoor = 0; icoor < fe.nbCoor(); icoor++)
            {
                sum2 += graduQuadPt[icoor] * graduQuadPt[icoor];
            }
            sum += sum2 * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the H1 norm of fct on the current element
template<typename FunctionType>
Real
elementaryFctH1NormSquare ( const FunctionType& fct, const CurrentFE& fe )
{
    Real sum (0.0);
    Real sum2 (0.0);
    Real x;
    Real y;
    Real z;

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);

        sum2 = std::pow (fct ( x, y, z ), 2);
        for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
        {
            sum2 += std::pow (fct.grad (iCoor, x, y, z), 2);
        }
        sum += sum2 * fe.weightDet ( iQuadPt );
    }
    return sum;
}

//! returns the square of the H1 norm of fct on the current element (time-dependent case)
template <typename FunctionType>
Real elementaryFctH1NormSquare ( const FunctionType& fct, const CurrentFE& fe, const Real t, const UInt nbComp )
{
    Real sum (0.0);
    Real sum2 (0.0);
    Real x;
    Real y;
    Real z;
    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);

        for ( UInt iComp = 0; iComp < nbComp; ++iComp )
        {
            sum2 = std::pow (fct (t, x, y, z, iComp ), 2);

            for (UInt iCoor = 0; iCoor < fe.nbCoor(); ++iCoor)
            {
                sum2 += std::pow (fct.grad (iCoor, t, x, y, z, iComp), 2);
            }
            sum += sum2 * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the L2 norm of (u-fct) on the current element
template <typename VectorType>
Real elementaryDifferenceL2NormSquare ( VectorType& u,
                                        boost::function<Real ( Real, Real, Real ) > fct,
                                        const CurrentFE& fe,
                                        const DOF& dof )
{
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt (0.0);
    Real diffQuadPt (0.0);

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt = 0.0;
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);
        for (UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobalMap ( eleID, iDof );
            uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
        }
        diffQuadPt = uQuadPt - fct ( x, y, z );
        sum += diffQuadPt * diffQuadPt * fe.weightDet ( iQuadPt );
    }
    return sum;
}

//! returns the square of the L2 norm of (u-fct) on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elementaryDifferenceL2NormSquare ( VectorType& u,
                                        boost::function<Real ( Real, Real, Real, Real, UInt ) > fct,
                                        const CurrentFE& fe,
                                        const DOF& dof, const Real t, const UInt nbComp )
{
    // returns the square of the L2 norm of (u-fct) on the current element

    UInt eleID ( fe.currentLocalId() );
    Real sum (0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt;
    Real diffQuadPt;

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);
        for (UInt iComp = 0; iComp < nbComp; ++iComp )
        {
            uQuadPt = 0.0;
            for (UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobalMap ( eleID, iDof ) + iComp * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
            }
            diffQuadPt = uQuadPt - fct ( t, x, y, z, iComp );
            sum += diffQuadPt * diffQuadPt * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the square of the H1 norm of (u-fct) on the current element
template <typename VectorType, typename UsrFct>
Real elementaryDifferenceH1NormSquare ( const VectorType& u, const UsrFct& fct, const CurrentFE& fe,
                                        const DOF& dof )
{
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt (0.0);
    Real diffQuadPt (0.0);

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt = 0.0;
        Vector graduQuadPt = ZeroVector (fe.nbCoor() );

        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);


        for (UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
        {
            UInt dofID = dof.localToGlobalMap ( eleID, iDof );
            uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
            for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
            {
                graduQuadPt (iCoor) += u ( dofID ) * fe.dphi ( iDof, iCoor, iQuadPt );
            }
        }

        diffQuadPt = uQuadPt - fct ( x, y, z );

        Vector diffGradQuadPt = graduQuadPt;
        for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
        {
            diffGradQuadPt (iCoor) -= fct.grad (iCoor, x, y, z);
        }
        Real sum2 = diffQuadPt * diffQuadPt;
        for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
        {
            sum2 += diffGradQuadPt (iCoor) * diffGradQuadPt (iCoor);
        }
        sum += sum2 * fe.weightDet ( iQuadPt );
    }
    return sum;
}

//! returns the square of the H1 norm of (u-fct) on the current element  (time-dependent case)
template <typename VectorType, typename UsrFct>
Real elementaryDifferenceH1NormSquare ( const VectorType& u, const UsrFct& fct, const CurrentFE& fe,
                                        const DOF& dof, const Real t, const UInt nbComp )
{
    UInt eleID = fe.currentLocalId();
    Real sum (0.0);
    Real sum2 (0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt (0.0);
    Real diffQuadPt (0.0);

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);

        for (UInt iComp (0); iComp < nbComp; ++iComp )
        {
            uQuadPt = 0.0;
            Vector graduQuadPt = ZeroVector (fe.nbCoor() );

            for (UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobalMap ( eleID, iDof ) + iComp * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
                for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
                {
                    graduQuadPt (iCoor) += u ( dofID ) * fe.dphi ( iDof, iCoor, iQuadPt );
                }
            }

            diffQuadPt = uQuadPt - fct (t, x, y, z, iComp);

            Vector diffGradQuadPt = graduQuadPt;
            for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
            {
                diffGradQuadPt (iCoor) -= fct.grad (iCoor, t, x, y, z, iComp);
            }

            sum2 = diffQuadPt * diffQuadPt;
            for (UInt iCoor (0); iCoor < fe.nbCoor(); ++iCoor)
            {
                sum2 += diffGradQuadPt (iCoor) * diffGradQuadPt (iCoor);
            }
            sum += sum2 * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}

//! returns the integral of (u-fct) of u on the current element
//! for time dependent+vectorial
template <typename VectorType>
Real elementaryDifferenceIntegral ( VectorType& u,
                                    boost::function<Real ( Real, Real, Real, Real, UInt ) > fct,
                                    const CurrentFE& fe,
                                    const DOF& dof, const Real t, const UInt nbComp = 1)
{
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real x;
    Real y;
    Real z;
    Real uQuadPt;
    Real diffQuadPt (0.0);

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);
        uQuadPt = 0.0;
        for (ID component = 0; component < nbComp; component++)
        {
            for ( UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobalMap ( eleID, iDof ) + component * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
            }
            diffQuadPt = uQuadPt - fct ( t, x, y, z, component );
        }
        sum += diffQuadPt * fe.weightDet ( iQuadPt );
    }
    return sum;
}

//! returns the integral of u on the current element
template <typename VectorType>
Real elementaryIntegral ( VectorType& u,
                          const CurrentFE& fe,
                          const DOF& dof, const UInt nbComp = 1)
{
    UInt eleID (fe.currentLocalId() );
    Real sum (0.0);
    Real uQuadPt (0.0);

    for ( UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        uQuadPt = 0.0;
        for (ID component = 0; component < nbComp; component++)
        {
            for ( UInt iDof (0); iDof < fe.nbFEDof(); ++iDof )
            {
                UInt dofID = dof.localToGlobalMap ( eleID, iDof ) + component * dof.numTotalDof();
                uQuadPt += u ( dofID ) * fe.phi ( iDof, iQuadPt );
            }
        }
        sum += uQuadPt * fe.weightDet ( iQuadPt );
    }
    return sum;
}

//! returns the integral of fct on the current element
inline Real
elementaryFctIntegral ( boost::function < Real ( Real, Real, Real,
                                                 Real, UInt ) > fct,
                        const CurrentFE& fe, const Real t, const UInt nbComp = 1)
{
    Real sum (0.0);
    Real x;
    Real y;
    Real z;

    for (UInt iQuadPt (0); iQuadPt < fe.nbQuadPt(); ++iQuadPt )
    {
        x = fe.quadNode (iQuadPt, 0);
        y = fe.quadNode (iQuadPt, 1);
        z = fe.quadNode (iQuadPt, 2);

        for (ID component = 0; component < nbComp; component++)
        {
            sum += fct ( t, x, y, z,  component ) * fe.weightDet ( iQuadPt );
        }
    }
    return sum;
}


} // namespace LifeV
#endif
