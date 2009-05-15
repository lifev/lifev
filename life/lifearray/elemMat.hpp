/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#ifndef _ELEMMAT_H_INCLUDED
#define _ELEMMAT_H_INCLUDED
#include <vector>
/*
#include "comprow_double.hpp"   // SparseLib++ matrix for "assemble_SparseLib"
*/

#include <life/lifearray/vecUnknown.hpp>

namespace LifeV
{
class ElemMat
{

public:
    typedef KNM<Real> matrix_type;
    typedef KNM_<Real> matrix_view;
    //typedef Tab2d matrix_type;

    ~ElemMat();
    ElemMat( UInt nNode1, UInt nbr1, UInt nbc1 ); // constructor for 1 finite element
    ElemMat( UInt nNode1, UInt nbr1, UInt nbc1,
             UInt nNode2, UInt nbr2, UInt nbc2 ); // constructor for 2 finite elements
    ElemMat( UInt nNode1, UInt nbr1, UInt nbc1,
             UInt nNode2, UInt nbr2, UInt nbc2,
             UInt nNode3, UInt nbr3, UInt nbc3 ); // constructor for 3 finite elements
    matrix_type& mat()
        {
            return _mat;
        }
    UInt nBlockRow() const
        {
            return _nBlockRow;
        }
    UInt nBlockCol() const
        {
            return _nBlockCol;
        }

    //Tab2dView block( UInt i, UInt j )
    matrix_view block( UInt i, UInt j )
        {
            return _mat( SubArray( _nRow[ i ], _firstRow[ i ] ),
                         SubArray( _nCol[ j ], _firstCol[ j ] ) );
            //Tab2dView __mr (_mat, TabRange(_firstRow[i], _nRow[i]), TabRange(_firstCol[j], _nCol[j]));
            //return __mr;
        }
    void zero()
        {
            //_mat = ZeroMatrix( _mat.size1(), _mat.size2() );
            _mat = 0.0;
        };
    void showMe( std::ostream& c = std::cout );

    void   operator *= (Real coef)
    {
        this->_mat *= coef;
    }
private:

    matrix_type _mat;

    UInt _nBlockRow; // number of block rows
    UInt _nBlockCol; // number of block columns
    std::vector<UInt> _nRow; // _nRow[i]=nb of rows in the i-th block row
    std::vector<UInt> _firstRow; //_firstRow[i]=index of first row of i-th block row
    std::vector<UInt> _nCol; // _nCol[i]=nb of col in the i-th block col
    std::vector<UInt> _firstCol; //_firstCol[i]=index of first col of i-th block col
};
}
#endif


