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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Aug 2010

    A more detailed description of the file (if necessary)
 */

#ifndef COMPOSEDDND_H
#define COMPOSEDDND_H 1

#include <life/lifecore/life.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>

namespace LifeV
{

//! ComposedDND - Modular preconditioner for geometry implicit monolithic FSI, three factors splitting
/*!
    @author Paolo Crosetto

Class implementing a modular preconditioner for FSI with the fluid geometry implicit. The preconditioner si split into
three factor, which can be recomputed every time or reused.
 */
class ComposedDND : public ComposedDN
{
public:

    typedef ComposedDN super;

    //! @name Constructors and destructor
    //@{

    ComposedDND( const std::vector<Int>& flag, const std::vector<Block>& order ):
            super(flag, order),
            M_swapped(false)
    {
    }

    ~ComposedDND() {}

    //@}


    void blockAssembling( );

private:

    bool M_swapped;

};

} // Namespace LifeV

#endif /* COMPOSEDDNN_H */
