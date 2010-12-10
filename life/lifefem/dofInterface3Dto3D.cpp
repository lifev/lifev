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
    @brief Class for interfacing dofs between two 3D meshes, implementation

    @author M.A. Fernandez
    @date 00-11-2002

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This file contains the class which may be used to update and hold the connections between the dof
    on two matching meshes.
 */

#include <life/lifefem/dofInterface3Dto3D.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

DofInterface3Dto3D::DofInterface3Dto3D( const RefFE& refFE, const Dof& dof1,
                                        const Dof& dof2 )
        :
        M_refFE1( &refFE ),
        M_dof1( &dof1 ),
        M_refFE2( &refFE ),
        M_dof2( &dof2 ),
        M_dof( new Dof( refFE ) )
{}

DofInterface3Dto3D::DofInterface3Dto3D( const RefFE& refFE1, const Dof& dof1, const RefFE& refFE2,
                                        const Dof& dof2 )
        :
        M_refFE1( &refFE1 ),
        M_dof1( &dof1 ),
        M_refFE2( &refFE2 ),
        M_dof2( &dof2 ),
        M_dof( new Dof( refFE1 ) )
{}

// ===================================================
// Methods
// ===================================================

void
DofInterface3Dto3D::setup( const RefFE& refFE, const Dof& dof1, const Dof& dof2 )
{
    M_refFE1 = &refFE;
    M_dof1 = &dof1;
    M_refFE2 = &refFE;
    M_dof2 = &dof2;
    M_dof = boost::shared_ptr<Dof>( new Dof( refFE ) );
}

void
DofInterface3Dto3D::setup( const RefFE& refFE1, const Dof& dof1, const RefFE& refFE2, const Dof& dof2 )
{
    M_refFE1 = &refFE1;
    M_dof1 = &dof1;
    M_refFE2 = &refFE2;
    M_dof2 = &dof2;
    M_dof = boost::shared_ptr<Dof>( new Dof( refFE1 ) );
}

// ===================================================
// Helpers
// ===================================================

bool coincide( const KN_<Real>& v1, const KN_<Real>& v2, const Real& tol )
{

    Real normDiff(0.0);

    for ( UInt i(0); i < nDimensions; ++i )
    {
        normDiff += fabs( v1[ i ] - v2[ i ] );
    }

    if ( normDiff <= tol )
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool coincide( const Real& x1, const Real& y1, const Real& z1, const Real& x2, const Real& y2, const Real& z2, const Real& tol )
{

    Real normDiff (fabs( x1 - x2 ) + fabs( y1 - y2 ) + fabs( z1 - z2 ));

    if ( normDiff <= tol )
    {
        return true;
    }
    else
    {
        return false;
    }
}

}
