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
   @file ETA_Blocks2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
*/

#ifndef __ETA_ABLOCKS2DTEST_H
#define __ETA_ABLOCKS2DTEST_H 1


// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh2DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>

#include <boost/shared_ptr.hpp>

#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/util/LifeChrono.hpp>


// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

/*!
    @class ETA_Blocks2DTest
    @brief Simple ETA test to compute a block matrix associated to the Stokes problem

    @author L. Pasquale <luca.pasquale@mail.polimi.it>

    This test computes the matrix and teh RHS vector associated to the Stokes problem: \n
    \f$
        \left\{\begin{array}{ll}
        - \delta u + \nabla p = f \\
        \nabla \cdot u=0
    \end{array}\right. \mathrm{in} (-1,1)\mathrm{x}(-1,1)
    \f$ \n
    with \f$u\f$ a scalar field approximated with P2 functions and p approximated with P1 \n
*/
class ETA_Blocks2DTest
    //     :
    //     public LifeV::Application
{
public:

    /* @name Constructors and destructor
     */
    //@{

    //! Constructor
    ETA_Blocks2DTest ();


    //! Destructor
    ~ETA_Blocks2DTest () {}

    //@}

    /* @name  Methods
     */
    //@{

    //! To lunch the simulation
    std::vector<LifeV::Real> run();

    //@}


private:
    boost::shared_ptr<Epetra_Comm>   M_comm;
    bool M_verbose;

};

#endif /* __ETA_ABLOCKS2DTEST_H */
