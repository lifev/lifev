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
#include <life/lifearray/vecUnknown.hpp>

namespace LifeV
{
//! the case of VectorBlock type
PhysVectUnknown<VectorBlock>::
PhysVectUnknown( UInt const Ndof ) : super( Ndof, nDimensions )
{}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl( const ScalUnknown<VectorBlock> &RhScU1,
               const ScalUnknown<VectorBlock> &RhScU2 ) :
        super( RhScU1.size() + RhScU2.size() ), _size( RhScU1.size() + RhScU2.size() ),
        _nbcomp( RhScU1.nbcomp() + RhScU2.nbcomp() )
{
    for ( UInt i = 0; i < RhScU1.size(); ++i )
        numBlock( i ) = RhScU1.numBlock( i );
    for ( UInt i = 0; i < RhScU2.size(); ++i )
        numBlock( i + RhScU1.size() ) = RhScU2.numBlock( i );
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl( const PhysVectUnknown<VectorBlock> &RhPhVU,
               const ScalUnknown<VectorBlock> &RhScU ) :
        super( RhPhVU.size() + RhScU.size() ), _size( RhPhVU.size() + RhScU.size() ),
        _nbcomp( RhPhVU.nbcomp() + RhScU.nbcomp() )
{
    for ( UInt i = 0; i < RhPhVU.size(); ++i )
        numBlock( i ) = RhPhVU.numBlock( i );
    for ( UInt i = 0; i < RhScU.size(); ++i )
        numBlock( RhPhVU.size() + i ) = RhScU.numBlock( i );
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl( const ScalUnknown<VectorBlock> &RhScU,
               const PhysVectUnknown<VectorBlock> &RhPhVU ) :
        super( RhPhVU.size() + RhScU.size() ), _size( RhPhVU.size() + RhScU.size() ),
        _nbcomp( RhPhVU.nbcomp() + RhScU.nbcomp() )
{
    for ( UInt i = 0; i < RhScU.size(); ++i )
        numBlock( i ) = RhScU.numBlock( i );
    for ( UInt i = 0; i < RhPhVU.size(); ++i )
        numBlock( RhScU.size() + i ) = RhPhVU.numBlock( i );
}

//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl( const PhysVectUnknown<VectorBlock> &RhPhVU,
               const GenericVecHdl<VectorBlock> &RhGenVec ) :
        super( RhPhVU.size() + RhGenVec.size() ),
        _size( RhPhVU.size() + RhGenVec.size() ),
        _nbcomp( RhPhVU.nbcomp() + RhGenVec.nbcomp() )
{
    for ( UInt i = 0;i < RhPhVU.size();++i )
        numBlock( i ) = RhPhVU.numBlock( i );
    for ( UInt i = 0;i < RhGenVec.size();++i )
        numBlock( RhPhVU.size() + i ) = RhGenVec.numBlock( i );
}
//! the case of VectorBlock type
GenericVecHdl<VectorBlock>::
GenericVecHdl( const GenericVecHdl<VectorBlock> &RhGenVec,
               const PhysVectUnknown<VectorBlock> &RhPhVU ) :
        super( RhPhVU.size() + RhGenVec.size() ),
        _size( RhPhVU.size() + RhGenVec.size() ),
        _nbcomp( RhPhVU.nbcomp() + RhGenVec.nbcomp() )
{
    for ( UInt i = 0;i < RhGenVec.size();++i )
        numBlock( i ) = RhGenVec.numBlock( i );
    for ( UInt i = 0;i < RhPhVU.size();++i )
        numBlock( i + RhGenVec.size() ) = RhPhVU.numBlock( i );
}
}
