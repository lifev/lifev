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
#include <life/lifefem/dof.hpp>

namespace LifeV
{
//! Constructor
Dof::Dof( const LocalDofPattern& _fe, UInt off ) : fe( _fe ), _offset( off ), _totalDof( 0 ),
        _nEl( 0 ), nlv( 0 ), nle( 0 ), nlf( 0 ), _ltg()
{
    for ( UInt i = 0; i < 5; ++i )
        _ncount[ i ] = 0;
}

//! Copy constructor
Dof::Dof( const Dof & dof2 ) : fe( dof2.fe ), _offset( dof2._offset ),
        _totalDof( dof2._totalDof ), _nEl( dof2._nEl ),
        nlv( dof2.nlv ), nle( dof2.nle ), nlf( dof2.nlf ),
        _ltg( dof2._ltg )
{
    if ( &dof2 == this )
        return ;

    //fe=dof2.fe;
    //  _offset=dof2._offset;
    //  _ltg=dof2._ltg;
    //  _totalDof=dof2._totalDof;
    for ( UInt i = 0; i < 5; ++i )
        _ncount[ i ] = dof2._ncount[ i ];
}

//! Ouput
void Dof::showMe( std::ostream & out, bool verbose ) const
{
    out << " Degree of Freedom (Dof) Object" << std::endl;
    out << " Total Dof Stored             " << _totalDof << std::endl;
    out << " With offset (min. Dof Id) =  " << _offset << std::endl;
    out << " Dof's on Vertices  from " << _ncount[ 0 ] << " , to:" << _ncount[ 1 ] - 1 << std::endl;
    out << " Dof's on Edges     from " << _ncount[ 1 ] << " , to:" << _ncount[ 2 ] - 1 << std::endl;
    out << " Dof's on Faces     from " << _ncount[ 2 ] << " , to:" << _ncount[ 3 ] - 1 << std::endl;
    out << " Dof's on Volumes   from " << _ncount[ 3 ] << " , to:" << _ncount[ 4 ] - 1 << std::endl;
    if ( verbose )
    {
        out << "************************************************************" << std::endl;
        out << "           Local to Global DOF table" << std::endl;
        out << "************************************************************" << std::endl;
        out << "Element Id   Loc. N.   Global N.  #  Element Id  Loc. N. Global N. " << std::endl;


        for ( UInt i = 0; i < _nEl;++i )
        {
            for ( UInt j = 0; j < numLocalDof();++j )
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

        }out << std::endl;

    }

}
}
