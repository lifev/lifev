/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef _MONOLITHIC_HPP
#define _MONOLITHIC_HPP


#include <life/lifesolver/FSIOperator.hpp>
//#include <Epetra_IntVector.h>

namespace LifeV
{

class Monolithic : public FSIOperator
{
public:


    typedef FSIOperator                            super;
    typedef FSIOperator::fluid_type::value_type::matrix_type   matrix_type;
    typedef boost::shared_ptr<matrix_type>                    matrix_ptrtype;


    // constructors
    Monolithic();

    // destructor

    ~Monolithic();



    void   evalResidual(vector_type&        res,
                        const vector_type& _disp,
                        const int          _iter);

    void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const double       _linearRelTol);

    void updateSystem(const vector_type& displacement);

    void setDataFromGetPot( GetPot const& data );

    void monolithicToInterface(    vector_type& lambdaSolid, const vector_type& disp) ;

    void monolithicToSolid(const vector_type& disp, vector_type& dispSolid);

    void setDispSolid();

    void monolithicToFluid(const vector_type& disp, vector_type& dispFluid);

    void iterateMesh(const vector_type& disp);
    void coupling(matrix_ptrtype & bigMatrix, vector_ptrtype rhs = vector_ptrtype());
    void buildSystem();
    void addDiagonalEntries(Real& entry, matrix_ptrtype & bigMatrix);
    void setup();
    EpetraMap& monolithicMap(){return M_monolithicMap;}

private:

    int                                               M_updateEvery;
    bool                                              firstIter;
    EpetraMap                                         M_monolithicMap;
    UInt                                              M_interface;
    bool                                              M_semiImplicit;
    EpetraMap                                         M_interfaceMap;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
    bool                                              M_isDiagonalBlockPrec;
    boost::shared_ptr<vector_type>                    M_beta;

    //    boost::shared_ptr<Epetra_IntVector>               M_numerationInterfaceInt;


};

}
#endif
