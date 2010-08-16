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
  \include ../../testsuite/test_monolithic/fluidstructure.dox
    @file
    @ This file contains the implementation of a composed preconditioner for FSI

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 08 Jun 2010

 */

#ifndef COMPOSEDDN2_H
#define COMPOSEDDN2_H 1

#include <lifemc/lifesolver/ComposedDN.hpp>

namespace LifeV {

//! ComposedDN2 - Class handling a composed preconditioner for FSI with two or three factors.
/*!
    @author Paolo Crosetto
    @see  \ref CDFQ

    The preconditioner implemented here is a variant of the one in ComposedDN. See \ref CDFQ for a complete explanation.
    The only peculiarity of this class with respect to ComposedDN is the different reordering of the factor. For this
    reason probably they will be merged together in the future.
 */
class ComposedDN2 : public ComposedDN
{
public:
    typedef ComposedDN super;

    //! @name Constructor
    //@{
    ComposedDN2( const std::vector<Int>& flags ):
        super( flags )
    {
    }
    //@}

    //! @name Public Methods
    //@{

    void coupler(map_shared_ptrtype& map,
                 const std::map<ID, ID>& locDofMap,
                 const vector_ptrtype& numerationInterface,
                 const Real& timeStep);
    //@}

private:

};

} // Namespace LifeV

#endif /* COMPOSEDDN2_H */
