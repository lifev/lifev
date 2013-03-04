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
    @brief This file contains the definition of the DOFLocalPattern class.

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#include <lifev/core/fem/DOFLocalPattern.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

// constructor
DOFLocalPattern::DOFLocalPattern ( const UInt& nbLocalDof, const UInt& nbDofPerVertex,
                                   const UInt& nbDofPerEdge, const UInt& nbDofPerFace,
                                   const UInt& nbDofPerVolume, const DofPatternType& patternType, UInt nbCoor  ) :
    M_dim (nbCoor), M_nbLocalDof ( nbLocalDof ), M_nbDofPerDimEntity (std::vector< UInt> (4) ),
    M_patternType ( patternType )
{
    // Store the location of the dofs
    M_nbDofPerDimEntity[0] = nbDofPerVertex;
    M_nbDofPerDimEntity[1] = nbDofPerEdge;
    M_nbDofPerDimEntity[2] = nbDofPerFace;
    M_nbDofPerDimEntity[3] = nbDofPerVolume;

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
            std::ostringstream errorMessage;
            errorMessage << "Pattern " << M_patternType << " not available for " << M_dim << "D. ";
            ERROR_MSG ( errorMessage.str().c_str() );
        }
    }
}



// The copy constructor

DOFLocalPattern::DOFLocalPattern ( const DOFLocalPattern& localDofPattern) :
    M_dim (localDofPattern.M_dim),
    M_nbLocalDof (localDofPattern.M_nbLocalDof ),
    M_nbDofPerDimEntity (localDofPattern.M_nbDofPerDimEntity),
    M_patternType (localDofPattern.M_patternType ),
    M_pattern (localDofPattern.M_pattern),
    M_nbPattern (localDofPattern.M_nbPattern),
    M_nbDiag (localDofPattern.M_nbDiag ),
    M_nbUpper (localDofPattern.M_nbUpper )
{}

// ===================================================
// Methods
// ===================================================

void DOFLocalPattern::showMe ( std::ostream& output) const
{
    output << " Size of the pattern : " << M_nbPattern << std::endl;
    output << " Diag: " << M_nbDiag << "  Upper: " << M_nbUpper << std::endl;
    output << " Pattern type: " << M_patternType << std::endl;

    for (UInt iter (0); iter < M_nbPattern; ++iter)
    {
        output << iter << " : " << M_pattern[iter].first << " - " << M_pattern[iter].second << std::endl;
    }
}

// ===================================================
// Private Methods
// ===================================================

void DOFLocalPattern::setupStandardPattern()
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
        M_pattern[i] = std::pair<UInt, UInt> (i, i);
    }

    // Upper diagonal entries
    int ip ( M_nbLocalDof );
    for ( UInt i = 0; i < M_nbLocalDof - 1; i++ )
    {
        for ( UInt j = i + 1; j < M_nbLocalDof; j++ )
        {
            M_pattern[ip] = std::pair<UInt, UInt> (i, j);
            ip++;
        }
    }

    // Lower diagonal entries
    for ( UInt i = 1; i < M_nbLocalDof; i++ )
    {
        for ( UInt j = 0; j < i; j++ )
        {
            M_pattern[ip] = std::pair<UInt, UInt> (i, j);
            ip++;
        }
    }
}


void DOFLocalPattern::setupP1isoP2SegPattern()
{
    // Some check to ensure consistency
    ASSERT (M_nbDofPerDimEntity[0] == 1, " Inconsistent P1 iso P2 (Vertices)");
    ASSERT (M_nbDofPerDimEntity[1] == 1, " Inconsistent P1 iso P2 (Edges)");
    ASSERT (M_nbDofPerDimEntity[2] == 0, " Inconsistent P1 iso P2 (Faces)");

    M_nbPattern = 7;
    M_nbDiag = 3;
    M_nbUpper = 2;
    M_pattern = std::vector< std::pair < UInt, UInt > > (M_nbPattern);

    // Diagonal entries
    for (UInt diag (0); diag < 3; ++diag)
    {
        M_pattern[diag] = std::pair<UInt, UInt> (diag, diag);
    }

    // Upper diagonal entries
    M_pattern[3]  = std::pair<UInt, UInt> (0, 2);
    M_pattern[4]  = std::pair<UInt, UInt> (1, 2);

    // Lower diagonal entries
    M_pattern[5] = std::pair<UInt, UInt> (2, 0);
    M_pattern[6] = std::pair<UInt, UInt> (2, 1);
}


void DOFLocalPattern::setupP1isoP2TriaPattern()
{
    // Some check to ensure consistency
    ASSERT (M_nbDofPerDimEntity[0] == 1, " Inconsistent P1 iso P2 (Vertices)");
    ASSERT (M_nbDofPerDimEntity[1] == 1, " Inconsistent P1 iso P2 (Edges)");
    ASSERT (M_nbDofPerDimEntity[2] == 0, " Inconsistent P1 iso P2 (Faces)");

    M_nbPattern = 24;
    M_nbDiag = 6;
    M_nbUpper = 9;
    M_pattern = std::vector< std::pair < UInt, UInt > > (M_nbPattern);

    // Diagonal entries
    for (UInt diag (0); diag < 6; ++diag)
    {
        M_pattern[diag] = std::pair<UInt, UInt> (diag, diag);
    };

    // Upper diagonal entries
    M_pattern[6]  = std::pair<UInt, UInt> (0, 3);
    M_pattern[7]  = std::pair<UInt, UInt> (0, 5);
    M_pattern[8]  = std::pair<UInt, UInt> (1, 3);
    M_pattern[9]  = std::pair<UInt, UInt> (1, 4);
    M_pattern[10] = std::pair<UInt, UInt> (2, 4);
    M_pattern[11] = std::pair<UInt, UInt> (2, 5);
    M_pattern[12] = std::pair<UInt, UInt> (3, 4);
    M_pattern[13] = std::pair<UInt, UInt> (3, 5);
    M_pattern[14] = std::pair<UInt, UInt> (4, 5);

    // Lower diagonal entries
    M_pattern[15] = std::pair<UInt, UInt> (3, 0);
    M_pattern[16] = std::pair<UInt, UInt> (3, 1);
    M_pattern[17] = std::pair<UInt, UInt> (4, 1);
    M_pattern[18] = std::pair<UInt, UInt> (4, 2);
    M_pattern[19] = std::pair<UInt, UInt> (4, 3);
    M_pattern[20] = std::pair<UInt, UInt> (5, 0);
    M_pattern[21] = std::pair<UInt, UInt> (5, 2);
    M_pattern[22] = std::pair<UInt, UInt> (5, 3);
    M_pattern[23] = std::pair<UInt, UInt> (5, 4);
}

} // end of the namespace LifeV
