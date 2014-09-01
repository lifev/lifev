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
/**
   @file multirate.hpp
   @author L. Oldani <luca.oldani@mail.polimi.it>
   @date 2013-12-27
 */


#ifndef __multirate_H
#define __multirate_H 1


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/solver/MultirateSolver.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#include <lifev/core/solver/HyperbolicData.hpp>
#include <lifev/core/fem/HyperbolicFlux.hpp>
#include <lifev/core/fem/HyperbolicNumericalFlux.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>


/*!
    @class multirate_2d
    @brief LifeV multirate test case

    @author L. Oldani <luca.oldani@mail.polimi.it>
*/

class multirate_2d
//     :
//     public LifeV::Application
{
public:

    /*! @name Constructors and destructor
     */
    //@{

    //! Constructor
    multirate_2d ( int argc, char* argv[] );

    //! Destructor
    ~multirate_2d () {}

    //@}

    /*! @name  Methods
     */
    //@{

    //! To lunch the simulation
    LifeV::Real run();

    //@}


private:
    struct Private;
    boost::shared_ptr<Private> Members;
};

#endif /* __multirate_H */
