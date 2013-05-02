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
    @brief Composed preconditioner for a three blocks coupled problem

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Aug 2010

    Solves a Dirichlet--Neumann--Dirichlet system for the three blocks
 */

#ifndef COMPOSEDDND_H
#define COMPOSEDDND_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>

namespace LifeV
{

//! MonolithicBlockComposedDND - Modular preconditioner for (e.g.) geometry implicit monolithic FSI, three factors splitting
/*!
    @author Paolo Crosetto

Class implementing a modular preconditioner for FSI with the fluid geometry implicit. The preconditioner si split into
three factor, which can be recomputed every time or reused.
 */
class MonolithicBlockComposedDND : public MonolithicBlockComposedDN
{
public:

    //! @name Public Types
    //@{

    typedef MonolithicBlockComposedDN super_Type;

    //@}
    //! @name Constructors and destructor
    //@{

    MonolithicBlockComposedDND ( const std::vector<Int>& flag, const std::vector<Int>& order ) :
        super_Type (flag, order),
        M_swapped (false)
    {
    }

    ~MonolithicBlockComposedDND() {}

    //@}
    //!@name Public Methods
    //@{

    void blockAssembling( );

    //@}



    static MonolithicBlock* createComposedDNDGI()
    {
        const Int order[] = {  MonolithicBlockComposed::mesh, MonolithicBlockComposed::solid, MonolithicBlockComposed::fluid };
        const Int couplingsDNGI2[] = { 0, 7, 0 };
        const std::vector<Int> couplingVectorDNGI2 (couplingsDNGI2, couplingsDNGI2 + 3);
        const std::vector<Int> orderVector (order, order + 3);
        return new MonolithicBlockComposedDND ( couplingVectorDNGI2, orderVector );
    }

    static MonolithicBlock* createComposedDND2GI()
    {
        const Int order[] = { MonolithicBlockComposed::mesh, MonolithicBlockComposed::fluid , MonolithicBlockComposed::solid};
        const Int couplingsDN2GI2[] = { 8, 6, 0 };
        const std::vector<Int> couplingVectorDN2GI2 (couplingsDN2GI2, couplingsDN2GI2 + 3);
        const std::vector<Int> orderVector (order, order + 3);
        return new MonolithicBlockComposedDND ( couplingVectorDN2GI2, orderVector );
    }


private:

    //!@name Protected Members
    //@{

    bool M_swapped;

    //@}

};

} // Namespace LifeV

#endif /* COMPOSEDDNN_H */
