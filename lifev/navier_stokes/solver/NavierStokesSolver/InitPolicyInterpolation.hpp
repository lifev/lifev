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
    @file InitPolicyInterpolation class
    @brief This class is a strategy to initialize a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 11-12-2012
 */

#ifndef INITPOLICYINTERPOLATION_HPP
#define INITPOLICYINTERPOLATION_HPP

#include <iostream>
#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>

namespace LifeV
{

struct InitPolicyInterpolation
{
    typedef VectorEpetra                             vector_Type;
    typedef boost::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef MapEpetra                                map_Type;
    typedef boost::shared_ptr<map_Type>              mapPtr_Type;
    typedef RegionMesh<LinearTetra>                  mesh_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef boost::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef boost::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef boost::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;
    typedef boost::shared_ptr< NavierStokesProblem > NSProblemPtr_Type;

    void setupInit ( Teuchos::ParameterList& list );

    void initSimulation ( bcContainerPtr_Type bchandler,
                          vectorPtr_Type solution );

    virtual fespacePtr_Type uFESpace() const = 0;
    virtual fespacePtr_Type pFESpace() const = 0;
    virtual Real timestep() const = 0;
    virtual Real initialTime() const = 0;
    virtual bdfPtr_Type bdf() const = 0;
    virtual NSProblemPtr_Type problem() const = 0;
};

} // namespace LifeV

#endif /* INITPOLICYINTERPOLATION_HPP */
