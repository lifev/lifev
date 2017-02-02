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
   @file hyperbolic.hpp
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
 */


#ifndef __hyperbolic_H
#define __hyperbolic_H 1


// ===================================================
//! Includes
// ===================================================

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/array/GhostHandler.hpp>
#include <lifev/core/solver/HyperbolicSolver.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>

/*!
 @class hyperbolic
 @brief LifeV hyperbolic test case

 @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
 */
class hyperbolic
    //     :
    //     public LifeV::Application
{
public:

    /*! @name Constructors and destructor
     */
    //@{

    hyperbolic ( int argc,
                 char** argv );

    ~hyperbolic()
    {}

    //@}

    /*! @name  Methods
     */
    //@{

    LifeV::Real run();

    //@}


private:
    struct Private;
    std::shared_ptr<Private> Members;
};

#endif /* __hyperbolic_H */
