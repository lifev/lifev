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
/* --------------------------------------------------------------------------*
/                                                                            /
/      ...                                                                   /
/                                                                            /
/                                                                            /
/ MATRIX VALUES                                                              /
/                                                                            /
/ #Version 0.0 Experimental:   21/06/2000 Alessandro Veneziani               /
/                                                                            /
/  I am sorry for my C++ at a very beginner level !!!!                       /
/                                                                            /
/ #Purpose: Provides the basic classes for sparse matrix handling            /
/ including common formats (CSR,MSR),                                        /
/ assembling in the LifeV environment and utilities                          /
/                                                                            /
/ #Note: It will be suitably modified and linked with                        /
/ yAMLib (by Bertolazzi, Formaggia, Manzini, Veneziani)                      /
/ and its sparse environment                                                 /
/                                                                            /
/ #attempt to define the VBR matrix class in order to handle vectorial       /
/  problems. 26/10/01, Alain Gauthier.                                       /
/                                                                            /
/ #modification of CSR matrice: _value can be a vector of matrices, each     /
/  matrix being a RNM class. Attempt of 15/11/01.                            /
/                                                                            /
/ #class used in IML++: construction of a DiagPreconditioner class.          /
/  21/11/01.                                                                 /
/                                                                            /
/ #defintest_fe/Neumann/vition of MixedMatr class: holds the values of block matrix           /
/  stored under a MixedPattern format. 30/04/02.                             /
/                                                                            /
/---------------------------------------------------------------------------*/
#ifndef _ASSEMBLE_MATRIX_HH
#define _ASSEMBLE_MATRIX_HH
#include <life/lifecore/life.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include  <life/lifearray/sparseArray.hpp>
#include <life/lifefem/dof.hpp>

#ifndef _LIFEV_HH_
//more correct version
typedef size_t UInt;
//original version
typedef std::vector<UInt>::iterator UIIter;
#endif

#ifndef _VEC_UNKNOWN_HH
#include <life/lifearray/vecUnknown.hpp>
#endif

namespace LifeV
{

/*---------------------------------------------------------------------------------/
 / I M P L E M E N T A T I O N S                                                    /
 /---------------------------------------------------------------------------------*/
template <typename FE1, typename GeoMap, typename MatrixType>
void assemble_first( MatrixType& M, ElemMat& elmat,
                     const FE1& fe1, const GeoMap& geo, const Dof& dof )
{
    if ( elmat.nBlockRow() != 1 || elmat.nBlockCol() != 1 )
    {
        std::cout << "assemble for vector elem mat not yet implemented\n";
        exit( 1 );
    }
    Tab2dView mat = elmat.block( 0, 0 );
    int i, j;
    UInt ig, jg, eleID = geo.currentID();
    for ( UInt k = 0 ; k < FE1::nbPattern ; k++ )
    {
        i = FE1::patternFirst( k );
        j = FE1::patternSecond( k );
        ig = dof.localToGlobal( eleID, i + 1 ) - 1;  // damned 1-base vs 0-base !
        jg = dof.localToGlobal( eleID, j + 1 ) - 1;  // damned 1-base vs 0-base !
        M.set_mat_inc( ig, jg, mat( i, j ) );
    }
}

template <typename Vector, typename FE1, typename GeoMap>
void assemble_vec( Vector& V, ElemVec& elvec, const FE1& fe1, const GeoMap& geo,
                   const Dof& dof )
{
    if ( elvec.nBlockRow() != 1 )
    {
        std::cout << "assemble for vector elem vec not yet implemented\n";
        exit( 1 );
    }
    Tab1dView vec = elvec.block( 0 );
    UInt i;
    UInt ig, eleID = geo.currentID();
    for ( i = 0 ; i < FE1::nbNode ; i++ )
    {
        ig = dof.localToGlobal( eleID, i + 1 ) - 1;
        V[ ig ] += vec( i );
    }
}
}
#endif
