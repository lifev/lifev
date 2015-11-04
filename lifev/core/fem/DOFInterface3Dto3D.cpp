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

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

DOFInterface3Dto3D::DOFInterface3Dto3D ( const ReferenceFE& refFE, const DOF& dof1,
                                         const DOF& dof2 )
    :
    M_refFE1 ( &refFE ),
    M_dof1 ( &dof1 ),
    M_refFE2 ( &refFE ),
    M_dof2 ( &dof2 ),
    M_dof ( new DOF ( refFE ) )
{}

DOFInterface3Dto3D::DOFInterface3Dto3D ( const ReferenceFE& refFE1, const DOF& dof1, const ReferenceFE& refFE2,
                                         const DOF& dof2 )
    :
    M_refFE1 ( &refFE1 ),
    M_dof1 ( &dof1 ),
    M_refFE2 ( &refFE2 ),
    M_dof2 ( &dof2 ),
    M_dof ( new DOF ( refFE1 ) )
{}

DOFInterface3Dto3D::DOFInterface3Dto3D ( const ReferenceFE& refFE, const DOF& dof )
    :
    M_refFE1 ( & refFE ),
    M_dof1 ( &dof ),
    M_refFE2 ( & refFE ),
    M_dof2 ( &dof )
{}

// ===================================================
// Methods
// ===================================================

void
DOFInterface3Dto3D::setup ( const ReferenceFE& refFE, const DOF& dof1, const DOF& dof2 )
{
    M_refFE1 = &refFE;
    M_dof1 = &dof1;
    M_refFE2 = &refFE;
    M_dof2 = &dof2;
    M_dof = std::shared_ptr<DOF> ( new DOF ( refFE ) );
}

void
DOFInterface3Dto3D::setup ( const ReferenceFE& refFE1, const DOF& dof1, const ReferenceFE& refFE2, const DOF& dof2 )
{
    M_refFE1 = &refFE1;
    M_dof1 = &dof1;
    M_refFE2 = &refFE2;
    M_dof2 = &dof2;
    M_dof = std::shared_ptr<DOF> ( new DOF ( refFE1 ) );
}

// ===================================================
// Helpers
// ===================================================

bool
coincide ( const std::vector<Real>& p1, const std::vector<Real>& p2, const Real& tol )
{
    Real normDiff = fabs ( p1[ 0 ] - p2[ 0 ]) + fabs ( p1[ 1 ] - p2[ 1 ]) + fabs ( p1[ 2 ] - p2[ 2 ]);

    return ( normDiff <= tol );
}

}
