//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Implementation of the localDofPattern class.
 */

#include <life/lifefem/localDofPattern.hpp>

namespace LifeV
{

// This is the 3D constructor

LocalDofPattern::LocalDofPattern( const UInt& _nbLocalDof, const UInt& _nbDofPerVertex,
                                  const UInt& _nbDofPerEdge, const UInt& _nbDofPerFace,
                                  const UInt& _nbDofPerVolume, const DofPatternType& _patternType ) :
        M_dim(3), M_nbLocalDof( _nbLocalDof ), M_nbDofPerDimEntity(std::vector< UInt> (4)),
        M_patternType( _patternType )
{
    // Store the location of the dofs
    M_nbDofPerDimEntity[0]=_nbDofPerVertex;
    M_nbDofPerDimEntity[1]=_nbDofPerEdge;
    M_nbDofPerDimEntity[2]=_nbDofPerFace;
    M_nbDofPerDimEntity[3]=_nbDofPerVolume;

    // Decide the pattern depending on the type
    switch ( M_patternType )
    {
    case STANDARD_PATTERN:
    {
        setupStandardPattern();
        break;
    }
    default:
    {
        std::ostringstream _err_msg;
        _err_msg << "Pattern " << M_patternType << " not available for " << M_dim << "D. ";
        ERROR_MSG( _err_msg.str().c_str() );
    };
    };
}


// This is the 2D constructor

LocalDofPattern::LocalDofPattern( const UInt& _nbLocalDof, const UInt& _nbDofPerVertex,
                                  const UInt& _nbDofPerEdge, const UInt& _nbDofPerFace,
                                  const DofPatternType& _patternType ) :
        M_dim(2), M_nbLocalDof( _nbLocalDof ), M_nbDofPerDimEntity(std::vector< UInt> (3)),
        M_patternType( _patternType )
{
    // Store the location of the dofs
    M_nbDofPerDimEntity[0]=_nbDofPerVertex;
    M_nbDofPerDimEntity[1]=_nbDofPerEdge;
    M_nbDofPerDimEntity[2]=_nbDofPerFace;

    // Decide the pattern depending on the type
    switch ( M_patternType )
    {
    case STANDARD_PATTERN:
    {
        setupStandardPattern();
        break;
    }
    case P1ISOP2_TRIA_PATTERN:
    {
        setupP1isoP2TriaPattern();
        break;
    }
    default:
    {
        std::ostringstream _err_msg;
        _err_msg << "Pattern " << M_patternType << " not available for " << M_dim << "D. ";
        ERROR_MSG( _err_msg.str().c_str() );
    };
    };
}


// This is the 1D constructor

LocalDofPattern::LocalDofPattern( const UInt& _nbLocalDof, const UInt& _nbDofPerVertex,
                                  const UInt& _nbDofPerEdge, const DofPatternType& _patternType ) :
        M_dim(1), M_nbLocalDof( _nbLocalDof ), M_nbDofPerDimEntity(std::vector< UInt> (2)),
        M_patternType( _patternType )
{
    // Store the location of the dofs
    M_nbDofPerDimEntity[0]=_nbDofPerVertex;
    M_nbDofPerDimEntity[1]=_nbDofPerEdge;

    // Decide the pattern depending on the type
    switch ( M_patternType )
    {
    case STANDARD_PATTERN:
    {
        setupStandardPattern();
        break;
    }
    case P1ISOP2_SEG_PATTERN:
    {
        setupP1isoP2SegPattern();
        break;
    }
    default:
    {
        std::ostringstream _err_msg;
        _err_msg << "Pattern " << M_patternType << " not available for " << M_dim << "D. ";
        ERROR_MSG( _err_msg.str().c_str() );

    };
    }; // end of the switch
};


// The copy constructor

LocalDofPattern::LocalDofPattern( const LocalDofPattern& _localDofPattern) :
        M_dim(_localDofPattern.M_dim),
        M_nbLocalDof (_localDofPattern.M_nbLocalDof ),
        M_nbDofPerDimEntity (_localDofPattern.M_nbDofPerDimEntity),
        M_patternType (_localDofPattern.M_patternType ),
        M_pattern (_localDofPattern.M_pattern),
        M_nbPattern(_localDofPattern.M_nbPattern),
        M_nbDiag (_localDofPattern.M_nbDiag ),
        M_nbUpper (_localDofPattern.M_nbUpper )
{};


void LocalDofPattern::showMe( std::ostream& output) const
{
    output << " Size of the pattern : " << M_nbPattern << std::endl;
    output << " Diag: " << M_nbDiag << "  Upper: " << M_nbUpper << std::endl;
    output << " Pattern type: " << M_patternType << std::endl;

    for (UInt iter(0); iter< M_nbPattern; ++iter)
    {
        output << iter << " : " << M_pattern[iter].first << " - " << M_pattern[iter].second << std::endl;
    };
}


void LocalDofPattern::setupStandardPattern()
{

    // This is the standard pattern with all the
    // degrees of freedom coupled together.

    // Number of couplings
    M_nbPattern = M_nbLocalDof * M_nbLocalDof;
    M_nbDiag = M_nbLocalDof;
    M_nbUpper = M_nbLocalDof * ( M_nbLocalDof - 1 ) / 2;

    // initialization of the pattern
    M_pattern = std::vector< std::pair< UInt, UInt > > (M_nbPattern);

    // First, put the diagonal entries
    for ( UInt i = 0; i < M_nbLocalDof; i++ )
    {
        M_pattern[i] = std::pair<UInt,UInt> (i,i);
    }

    // Upper diagonal entries
    int ip ( M_nbLocalDof );
    for ( UInt i = 0; i < M_nbLocalDof - 1; i++ )
    {
        for ( UInt j = i + 1; j < M_nbLocalDof; j++ )
        {
            M_pattern[ip] = std::pair<UInt,UInt> (i,j);
            ip++;
        }
    }

    // Lower diagonal entries
    for ( UInt i = 1; i < M_nbLocalDof; i++ )
    {
        for ( UInt j = 0; j < i; j++ )
        {
            M_pattern[ip] = std::pair<UInt,UInt> (i,j);
            ip++;
        }
    };
};


void LocalDofPattern::setupP1isoP2SegPattern()
{
    // Some check to ensure consistency
    ASSERT(M_nbDofPerDimEntity[0] == 1, " Inconsistent P1 iso P2 (Vertices)");
    ASSERT(M_nbDofPerDimEntity[1] == 1, " Inconsistent P1 iso P2 (Edges)");
    ASSERT(M_nbDofPerDimEntity[2] == 0, " Inconsistent P1 iso P2 (Faces)");

    M_nbPattern = 7;
    M_nbDiag = 3;
    M_nbUpper = 2;
    M_pattern = std::vector< std::pair < UInt,UInt > > (M_nbPattern);

    // Diagonal entries
    for (UInt diag(0); diag<3; ++diag)
    {
        M_pattern[diag] = std::pair<UInt,UInt> (diag,diag);
    };

    // Upper diagonal entries
    M_pattern[3]  = std::pair<UInt,UInt> (0,2);
    M_pattern[4]  = std::pair<UInt,UInt> (1,2);

    // Lower diagonal entries
    M_pattern[5] = std::pair<UInt,UInt> (2,0);
    M_pattern[6] = std::pair<UInt,UInt> (2,1);
}


void LocalDofPattern::setupP1isoP2TriaPattern()
{
    // Some check to ensure consistency
    ASSERT(M_nbDofPerDimEntity[0] == 1, " Inconsistent P1 iso P2 (Vertices)");
    ASSERT(M_nbDofPerDimEntity[1] == 1, " Inconsistent P1 iso P2 (Edges)");
    ASSERT(M_nbDofPerDimEntity[2] == 0, " Inconsistent P1 iso P2 (Faces)");

    M_nbPattern = 24;
    M_nbDiag = 6;
    M_nbUpper = 9;
    M_pattern = std::vector< std::pair < UInt,UInt > > (M_nbPattern);

    // Diagonal entries
    for (UInt diag(0); diag<6; ++diag)
    {
        M_pattern[diag] = std::pair<UInt,UInt> (diag,diag);
    };

    // Upper diagonal entries
    M_pattern[6]  = std::pair<UInt,UInt> (0,3);
    M_pattern[7]  = std::pair<UInt,UInt> (0,5);
    M_pattern[8]  = std::pair<UInt,UInt> (1,3);
    M_pattern[9]  = std::pair<UInt,UInt> (1,4);
    M_pattern[10] = std::pair<UInt,UInt> (2,4);
    M_pattern[11] = std::pair<UInt,UInt> (2,5);
    M_pattern[12] = std::pair<UInt,UInt> (3,4);
    M_pattern[13] = std::pair<UInt,UInt> (3,5);
    M_pattern[14] = std::pair<UInt,UInt> (4,5);

    // Lower diagonal entries
    M_pattern[15] = std::pair<UInt,UInt> (3,0);
    M_pattern[16] = std::pair<UInt,UInt> (3,1);
    M_pattern[17] = std::pair<UInt,UInt> (4,1);
    M_pattern[18] = std::pair<UInt,UInt> (4,2);
    M_pattern[19] = std::pair<UInt,UInt> (4,3);
    M_pattern[20] = std::pair<UInt,UInt> (5,0);
    M_pattern[21] = std::pair<UInt,UInt> (5,2);
    M_pattern[22] = std::pair<UInt,UInt> (5,3);
    M_pattern[23] = std::pair<UInt,UInt> (5,4);
}

} // end of the namespace LifeV
