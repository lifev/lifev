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
    @brief Degrees of freedom, the class that provides the localtoglobal table

    @author M.A. Fernandez & Luca Formaggia
    @date 00-07-2002

    @contributor Vincent Martin
                 Mohamed Belhadj
                 Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#include <life/lifefem/dof.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

Dof::Dof( const LocalDofPattern& fePattern, UInt offset ) : M_elementDofPattern( fePattern ), M_offset( offset ), M_totalDof( 0 ),
        M_numElement( 0 ), M_nbLocalVertex( 0 ), M_nbLocalEdge( 0 ), M_nbLocalFace( 0 ), M_localToGlobal(),
        M_nbFace(0),M_localToGlobalByFace(),M_globalToLocalByFace()
{
    //Getting the face
    switch ( fePattern.nbLocalDof() )
    {
    case 2:
        // No M_faceToPoint (it is 1D)
        M_numLocalDofByFace = 1;
        break;
    case 4:
        M_faceToPoint = LinearTetra::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 5:
        M_faceToPoint = LinearTetraBubble::faceToPoint;
        M_numLocalDofByFace = 3;
        break;
    case 10:
        M_faceToPoint = QuadraticTetra::faceToPoint;
        M_numLocalDofByFace = 6;
        break;
    case 8:
        M_faceToPoint = LinearHexa::faceToPoint;
        M_numLocalDofByFace = 4;
        break;
    case 27:
        M_faceToPoint = QuadraticHexa::faceToPoint;
        M_numLocalDofByFace = 27;
        break;
    default:
        std::cout << "Warning: This refFE is not available for the dof by face." << std::endl;
        M_numLocalDofByFace = 0;
        break;
    }

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = 0;
}

Dof::Dof( const Dof & dof2 ) : M_elementDofPattern( dof2.M_elementDofPattern ), M_offset( dof2.M_offset ),
        M_totalDof( dof2.M_totalDof ), M_numElement( dof2.M_numElement ),
        M_nbLocalVertex( dof2.M_nbLocalVertex ), M_nbLocalEdge( dof2.M_nbLocalEdge ), M_nbLocalFace( dof2.M_nbLocalFace ),
        M_localToGlobal( dof2.M_localToGlobal ),M_nbFace(dof2.M_nbFace),
        M_localToGlobalByFace(dof2.M_localToGlobalByFace),M_globalToLocalByFace(dof2.M_globalToLocalByFace),
        M_faceToPoint(dof2.M_faceToPoint),M_numLocalDofByFace(dof2.M_numLocalDofByFace)
{
    if ( &dof2 == this )
        return ;

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = dof2.M_dofPositionByEntity[ i ];
}

// ===================================================
// Methods
// ===================================================

ID Dof::localToGlobalByFace(const ID& faceId, const ID& localDof, bool& exist ) const
{
    ASSERT_PRE( (M_numLocalDofByFace>0) , "This data are not available for this reference element");
    std::map<ID,ID>::const_iterator mapIt(M_globalToLocalByFace.find(faceId) );

    if (mapIt != M_globalToLocalByFace.end())
    {
        exist = true;
        return M_localToGlobalByFace[(*mapIt).second][localDof-1];
    }
    else
    {
        exist = false;
        return 0;
    }
}

void Dof::showMe( std::ostream & out, bool verbose ) const
{
    out << " Degree of Freedom (Dof) Object" << std::endl;
    out << " Total Dof Stored             " << M_totalDof << std::endl;
    out << " With offset (min. Dof Id) =  " << M_offset << std::endl;
    out << " Dof's on Vertices  from " << M_dofPositionByEntity[ 0 ] << " , to:" << M_dofPositionByEntity[ 1 ] - 1 << std::endl;
    out << " Dof's on Edges     from " << M_dofPositionByEntity[ 1 ] << " , to:" << M_dofPositionByEntity[ 2 ] - 1 << std::endl;
    out << " Dof's on Faces     from " << M_dofPositionByEntity[ 2 ] << " , to:" << M_dofPositionByEntity[ 3 ] - 1 << std::endl;
    out << " Dof's on Volumes   from " << M_dofPositionByEntity[ 3 ] << " , to:" << M_dofPositionByEntity[ 4 ] - 1 << std::endl;
    if ( verbose )
    {
        out << "************************************************************" << std::endl;
        out << "           Local to Global DOF table" << std::endl;
        out << "************************************************************" << std::endl;
        out << "Element Id   Loc. N.   Global N.  #  Element Id  Loc. N. Global N. " << std::endl;


        for ( UInt i = 0; i < M_numElement; ++i )
        {
            for ( UInt j = 0; j < numLocalDof(); ++j )
            {
                out.width( 10 );
                out << i + 1;
                out.width( 10 );
                out << j + 1;
                out.width( 10 );
                out << localToGlobal( i + 1, j + 1 );
                out << " # ";
                if ( (i*numLocalDof()+j) % 2 != 0 )
                    out << std::endl;
            }

        }
        out << std::endl;

    }

}

void Dof::showMeByFace(std::ostream& out, bool verbose) const
{
    out << "--------------------------------------------------------------------------------" << std::endl;
    out << " Degree of freedom by face object " << std::endl;
    out << "--------------------------------------------------------------------------------" << std::endl;

    out << " Offset (min Dof Id) = " << M_offset << std::endl;
    out << " Number of local dof per face = " << M_numLocalDofByFace << std::endl;

    if (verbose)
    {
        out << "*********************************************************************************" << std::endl;
        out << " Local-to-global DOF table (DOF grouped by internal face)" << std::endl;
        out << "*********************************************************************************" << std::endl;
        out << "=================================================================================" << std::endl;
        out << "Face ID     Local DOF   Global DOF  " << std::endl;
        out << "=================================================================================" << std::endl;

        for (UInt i = 0; i < M_nbFace; ++i)
        {
            for (UInt j = 0; j < M_numLocalDofByFace; ++j)
            {
                out.width(12);
                out << i + 1;
                out.width(12);
                out << j + 1;
                out.width(12);
                out << localToGlobal(i+1, j+1);
                out << " # ";
                if (j % 2 != 0) out << std::endl;
            } // for j
        } //for i
    } // if verbose
}

//End of namespace LifeV
}
