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


#include <life/lifefem/DOF.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

DOF::DOF( const DOFLocalPattern& fePattern) : M_elementDofPattern( fePattern ), M_totalDof( 0 ),
        M_numElement( 0 ), M_nbLocalPeaks( 0 ), M_nbLocalRidges( 0 ), M_nbLocalFacets( 0 ), M_localToGlobal(),
        M_nbFacets(0),M_localToGlobalByFacet(),M_globalToLocalByFacet()
{
    //Getting the face
    switch ( fePattern.nbLocalDof() )
    {
    case 1: //P0 Q0
    	M_numLocalDofByFace = 0;
    	break;
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
        std::cout << "Warning: This refFE is not available for the dof by face."  <<  std::endl;
        M_numLocalDofByFace = 0;
        break;
    }

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = 0;
}

DOF::DOF( const DOF & dof2 ) : M_elementDofPattern( dof2.M_elementDofPattern ), //, M_offset( dof2.M_offset ),
        M_totalDof( dof2.M_totalDof ), M_numElement( dof2.M_numElement ),
        M_nbLocalPeaks( dof2.M_nbLocalPeaks ), M_nbLocalRidges( dof2.M_nbLocalRidges ), M_nbLocalFacets( dof2.M_nbLocalFacets ),
        M_localToGlobal( dof2.M_localToGlobal ),M_nbFacets(dof2.M_nbFacets),
        M_localToGlobalByFacet(dof2.M_localToGlobalByFacet),M_globalToLocalByFacet(dof2.M_globalToLocalByFacet),
        M_facetToPoint(dof2.M_facetToPoint),M_numLocalDofByFacet(dof2.M_numLocalDofByFacet)
{
    if ( &dof2 == this )
        return ;

    for ( UInt i = 0; i < 5; ++i )
        M_dofPositionByEntity[ i ] = dof2.M_dofPositionByEntity[ i ];
}

// ===================================================
// Methods
// ===================================================

ID DOF::localToGlobalMapByBdFacet(const ID& facetId, const ID& localDof ) const
{
	ASSERT_PRE( (M_localToGlobalByBdFacet.size()>0) , "The local to global map by boundary facet is void");
	return M_localToGlobalByBdFacet[facetId][localDof];
}


void DOF::showMe( std::ostream & out, bool verbose ) const
{
    out << " Degree of Freedom (DOF) Object" << std::endl;
    out << " Total DOF Stored             " << M_totalDof << std::endl;
    out << " DOF's on Vertices  from " << M_dofPositionByEntity[ 0 ] << " , to:" << M_dofPositionByEntity[ 1 ] << std::endl;
    out << " DOF's on Edges     from " << M_dofPositionByEntity[ 1 ] << " , to:" << M_dofPositionByEntity[ 2 ] << std::endl;
    out << " DOF's on Faces     from " << M_dofPositionByEntity[ 2 ] << " , to:" << M_dofPositionByEntity[ 3 ] << std::endl;
    out << " DOF's on Volumes   from " << M_dofPositionByEntity[ 3 ] << " , to:" << M_dofPositionByEntity[ 4 ] << std::endl;
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
                out << i;
                out.width( 10 );
                out << j;
                out.width( 10 );
                out << localToGlobalMap( i, j );
                out << " # ";
                if ( (i*numLocalDof()+j) % 2 != 0 )
                    out << std::endl;
            }

        }
        out << std::endl;

    }

}

void DOF::showMeByFace(std::ostream& out, bool verbose) const
{
    out << "--------------------------------------------------------------------------------" << std::endl;
    out << " Degree of freedom by face object " << std::endl;
    out << "--------------------------------------------------------------------------------" << std::endl;

    out << " Number of local dof per face = " << M_numLocalDofByFacet << std::endl;

    if (verbose)
    {
        out << "*********************************************************************************" << std::endl;
        out << " Local-to-global DOF table (DOF grouped by internal face)" << std::endl;
        out << "*********************************************************************************" << std::endl;
        out << "=================================================================================" << std::endl;
        out << "Face ID     Local DOF   Global DOF  " << std::endl;
        out << "=================================================================================" << std::endl;

        for (UInt i = 0; i < M_nbFacets; ++i)
        {
            for (UInt j = 0; j < M_numLocalDofByFacet; ++j)
            {
                out.width(12);
                out << i;
                out.width(12);
                out << j;
                out.width(12);
                out << localToGlobalMap(i, j);
                out << " # ";
                if (j % 2 != 0) out << std::endl;
            } // for j
        } //for i
    } // if verbose
}

//End of namespace LifeV
}
