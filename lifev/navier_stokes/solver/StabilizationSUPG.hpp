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
    @brief SUPG Stabilization

    @author Davide Forti <davide.forti@epfl.ch>
    @mantainer Davide Forti <davide.forti@epfl.ch>
    @date 04-02-2014

    Implementation of SUPG stabilization for inf-sup incompatible finite elements for the Navier-Stokes equations. For the moment it has been coded only for semi-implicit
    approximation of the convective term.
 */

#ifndef STABILIZATIONSUPG_HPP
#define STABILIZATIONSUPG_HPP

#include <lifev/core/util/LifeChrono.hpp>

namespace LifeV
{

//! StabilizationIP Class
/*!
 * @brief SUPG Stabilization
 * @author Davide Forti <davide.forti@epfl.ch>
 *
 * Implementation of SUPG. stabilization for inf-sup incompatible finite elements
 * for the Navier-Stokes equations. <br>
 *
 *
 */

template<typename MeshType, typename MapType, UInt SpaceDim>
class StabilizationSUPG
{
public:

    //! @name Public Types
    //@{

    //@}

    //! @name Constructor and Destructor
    //@{

    //! Default Constructor
    StabilizationSUPG();

    ~StabilizationSUPG() {};
    //@}


private:


}; // class StabilizationSUPG


//=============================================================================
// Constructor
//=============================================================================

template<typename MeshType, typename MapType, UInt SpaceDim>
StabilizationSUPG<MeshType, MapType, SpaceDim>::StabilizationSUPG()
{
}

//=============================================================================
// Method
//=============================================================================


//=============================================================================
// Setters method
//=============================================================================


} // namespace LifeV

#endif /* STABILIZATIONSUPG_HPP */
