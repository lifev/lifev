/*-*- mode: c++ -*-
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
#include <life/lifefem/refHybridFE.hpp>
#include <set>

namespace LifeV
{

RefHybridFE::RefHybridFE( const UInt& nbdfe, const StaticBdFE* bdfelist,
                          std::string _name, int _type, ReferenceShapes _shape,
                          int _nbDofPerVertex, int _nbDofPerEdge,
                          int _nbDofPerFace, int _nbDofPerVolume,
                          int _nbDof, int _nbCoor, const Real* refCoor,
                          PatternType _patternType ) :
        LocalDofPattern( _nbDof, _nbDofPerVertex, _nbDofPerEdge, _nbDofPerFace, _nbDofPerVolume, _patternType ),
        _nBdFE( nbdfe ), _bdfeList( bdfelist ), _refCoor( refCoor ),
        name( _name ), type( _type ), shape( _shape ),
        nbDof( _nbDof ), nbCoor( _nbCoor )
{
    CONSTRUCTOR( "RefHybridFE" );

    //! simple consistency test: (to be removed some day)
    Int nbdofbdfe = 0;
    for ( UInt nf = 0 ; nf < _nBdFE ; nf ++ )
        nbdofbdfe += _bdfeList[ nf ].nbNode;
    ASSERT_PRE( nbdofbdfe == nbDof , \
                "There should be as many dof in the refHybrid element as the sum of all dof in the boundary elements." );
}

RefHybridFE::~RefHybridFE()
{
    DESTRUCTOR( "RefHybridFE" );
}

//! extracting a BdFE from the boundary elements list. //to be checked
const StaticBdFE& RefHybridFE::operator[] ( const Index_t& i ) const
{
    ASSERT_BD( i < ( Index_t ) _nBdFE );
    return _bdfeList[ i ];
}

void RefHybridFE::check() const
{
    Real sumphi, sumdphi;
    int nbdofbdfe = 0;
    std::cout << "*** Check " << name << std::endl;
    for ( UInt nf = 0 ; nf < _nBdFE ; nf ++ )
    {
        const StaticBdFE& bdfe = _bdfeList[ nf ];
        const QuadRule& qr = bdfe.qr;
        std::cout << std::endl << "    " << bdfe.refFE.name << ", whose identity is " << bdfe.currentId() << std::endl;
        std::cout << "    " << qr.name << std::endl;
        nbdofbdfe += bdfe.nbNode;

        for ( int i = 0 ; i < bdfe.nbNode ; i ++ )
        {
            sumphi = sumdphi = 0.;
            for ( int ig = 0 ; ig < qr.nbQuadPt ; ig ++ )
            {
                sumphi += bdfe.phi( i, ig ) * bdfe.weightMeas( ig );
                for ( int icoor = 0 ; icoor < bdfe.nbCoor ; icoor ++ )
                    sumdphi += bdfe.dPhiRef( i, icoor, ig ) * bdfe.weightMeas( ig );
            }
            std::cout << "      integral_Face phi_i                                   = " << sumphi << std::endl;
            std::cout << "      sum_{icoor} integral dphi_i/dx_icoor                  = " << sumdphi << std::endl;
        }
    }
    std::cout << "     number of dof (should be : " << nbDof << " )             = " << nbdofbdfe << std::endl;
    //! simple consistency test:
    ASSERT_PRE( nbdofbdfe == nbDof , \
                "There should be as many dof in the refHybrid element as the sum of all dof in the boundary elements." );
    std::cout << std::endl;
}

std::ostream& operator << ( std::ostream& f, const RefHybridFE& fe )
{
    f << "-------------------------\n";
    f << "Reference Finite Element: " << fe.name << std::endl;
    f << "-------------------------\n";
    f << "*** Shape : " << fe.shape << std::endl;
    f << "*** Local coordinate of the nodes :\n";
    for ( int i = 0;i < fe.nbDof;i++ )
    {
        f << fe.xi( i ) << " " << fe.eta( i ) << " " << fe.zeta( i ) << std::endl;
    }
    f << "*** Pattern :\n";
    for ( int i = 0;i < fe.nbPattern();i++ )
        f << "(" << fe.patternFirst( i ) << "," << fe.patternSecond( i ) << ") \n";
    return f;
}
}
