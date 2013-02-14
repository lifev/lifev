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
    @brief Vector for elementary assembly

    @contributor Matteo Astorino <matteo.astorino@epfl.ch>
    @mantainer Matteo Astorino <matteo.astorino@epfl.ch>

 */

#ifndef ELEMVEC_H
#define ELEMVEC_H

#include <iomanip>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/RNM.hpp>

namespace LifeV
{
class VectorElemental
    :
public KN<Real>
    //public Tab1d
{
public:
    //typedef Tab1d super;
    typedef KN<Real> super;
    typedef KN_<Real> vector_view;


    VectorElemental ( int nNode1, int nbr1 );
    VectorElemental ( int nNode1, int nbr1,
                      int nNode2, int nbr2 );
    VectorElemental ( int nNode1, int nbr1,
                      int nNode2, int nbr2,
                      int nNode3, int nbr3 );

    VectorElemental& operator= ( super const& __v )
    {
        if ( this == &__v )
        {
            return * this;
        }
        super::operator= ( ( super const& ) __v );
        return *this;
    }

    super& vec()
    {
        return * this;
    };
    const super& vec() const
    {
        return * this;
    };
    int nBlockRow() const
    {
        return _nBlockRow;
    }
    //Tab1dView block( int i )
    vector_view block ( int i )
    {
        return ( *this ) ( SubArray ( _nRow[ i ], _firstRow[ i ] ) );
        //return Tab1dView( *this, TabRange( _firstRow[i], _nRow[ i ] ) );
    }
    //  inline void zero(){_vec=Tab1d(_nBlockRow*_vec.N(),0.0);};
    void zero()
    {
        //( *this ) = ZeroVector( this->size() );
        super& __super = ( super& ) * this;
        __super = super ( this->N(), 0.0 );
    }
    void showMe ( std::ostream& c = std::cout );
private:
    int _nBlockRow; // number of block rows
    std::vector<int> _nRow; // _nRow[i]=nb of rows in the i-th block row
    std::vector<int> _firstRow; //_firstRow[i]=index of first row of i-th block row

};
}


#endif

