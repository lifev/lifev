/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#include <life/lifefem/localDofPattern.hpp>

namespace LifeV
{
LocalDofPattern::LocalDofPattern( int _nbLocalDof, int _nbDofPerVertex, int _nbDofPerEdge, int _nbDofPerFace,
                                  int _nbDofPerVolume, PatternType _patternType ) :
        nbLocalDof( _nbLocalDof ), nbDofPerVertex( _nbDofPerVertex ), nbDofPerEdge( _nbDofPerEdge ),
        nbDofPerFace( _nbDofPerFace ), nbDofPerVolume( _nbDofPerVolume ),
        patternType( _patternType )
        // nbdof is total number of dof (not used for the initialization of a member)
{
    // Pattern
    switch ( patternType )
    {
    case STANDARD_PATTERN:
        {
            // full element matrix
            _nbPattern = nbLocalDof * nbLocalDof;
            _nbDiag = nbLocalDof;
            _nbUpper = nbLocalDof * ( nbLocalDof - 1 ) / 2;
            _patternFirst = new int[ _nbPattern ];
            _patternSecond = new int[ _nbPattern ];
            for ( int i = 0;i < nbLocalDof;i++ )
                _patternFirst[ i ] = _patternSecond[ i ] = i;
            int ip = nbLocalDof;
            for ( int i = 0;i < nbLocalDof - 1;i++ )
            {
                for ( int j = i + 1;j < nbLocalDof;j++ )
                {
                    _patternFirst[ ip ] = i;
                    _patternSecond[ ip ] = j;
                    ip++;
                }
            }
            for ( int i = 1;i < nbLocalDof;i++ )
            {
                for ( int j = 0;j < i;j++ )
                {
                    _patternFirst[ ip ] = i;
                    _patternSecond[ ip ] = j;
                    ip++;
                }
            }
            break;
        }
    case P1ISOP2_TRIA_PATTERN:
        {
            _nbPattern = 24;
            _nbDiag = 6;
            _nbUpper = 9;
            typedef std::pair<int, int> COO;
            COO pattern_p1isop2_tria[ 24 ] = {
                                                 COO( 0, 0 ), COO( 1, 1 ), COO( 2, 2 ), COO( 3, 3 ), COO( 4, 4 ), COO( 5, 5 ),  // diagonal entries
                                                 COO( 0, 3 ), COO( 0, 5 ), COO( 1, 3 ), COO( 1, 4 ), COO( 2, 4 ), COO( 2, 5 ), COO( 3, 4 ), COO( 3, 5 ), COO( 4, 5 ),  //upper entries
                                                 COO( 3, 0 ), COO( 3, 1 ), COO( 4, 1 ), COO( 4, 2 ), COO( 4, 3 ), COO( 5, 0 ), COO( 5, 2 ), COO( 5, 3 ), COO( 5, 4 )  // lower entries
                                             };
            for ( int i = 0;i < 24;i++ )
            {
                _patternFirst[ i ] = pattern_p1isop2_tria[ i ].first;
                _patternSecond[ i ] = pattern_p1isop2_tria[ i ].second;
            }
        }
    default:
    		std::ostringstream _err_msg;
    		_err_msg << "Unknown pattern " << patternType << "I cannot build local pattern!";
        ERROR_MSG( _err_msg.str().c_str() );
    }
}
}
