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
   @file ETA_ADR2DTest.cpp
   @author L. Pasquale <luca.pasquale@mail.polimi.it>
   @date 2012-11-20
*/

#ifndef __ETA_ADR2DTEST_H
#define __ETA_ADR2DTEST_H 1


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

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// We include two additional files which are useful to make the
// assembly using the "classical way". We also include a file to
// monitor the different timings.
// ---------------------------------------------------------------

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>

#include <lifev/core/util/LifeChrono.hpp>


// ---------------------------------------------------------------
// We work in the LifeV namespace and define the mesh, matrix and
// vector types that we will need several times.
// ---------------------------------------------------------------

/*!
    @class ETA_ADR2DTest
    @brief Simple ETA test to compare ETA to ADRAssembly in a 2D space

    @author L. Pasquale <luca.pasquale@mail.polimi.it>

    This test computes the matrix associated to the advection-diffusion-reaction problem: \n
    \f$
        -\alpha \delta u + \beta \cdot \nabla u +\gamma u = f
    \f$ \n
    with \f$\alpha=1\f$, \f$\beta=\left(\begin{array}{c}1\\0\end{array}\right)\f$, \f$\gamma=2\f$, \f$f=2\f$ and \f$u\f$ a scalar field \n
    The domain is the square (-1,1)x(-1,1). \n
    The matrix and the vector for the RHS are first assembled using the class ADRAssembler and then the ExpressionAssembly::integrate function of the ETA module. \n
    if the two matrices (and vectors) are assembled correctly, they should be almost exactly the same. \n
*/
class ETA_ADR2DTest
//     :
//     public LifeV::Application
{
public:

    /* @name Constructors and destructor
     */
    //@{

    //! Constructor
    ETA_ADR2DTest ();


    //! Destructor
    ~ETA_ADR2DTest () {}

    //@}

    /* @name  Methods
     */
    //@{

    //! To lunch the simulation
    LifeV::Real run();

    //@}


private:
    boost::shared_ptr<Epetra_Comm>   M_comm;
    bool M_verbose;

};

#endif /* __ETA_ADR2DTEST_H */
