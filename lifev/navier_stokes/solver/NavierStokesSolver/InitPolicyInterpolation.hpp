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
#include <string>
#include <boost/shared_ptr.hpp>


#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>


#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>

namespace LifeV
{

template< class mesh_Type >
struct InitPolicyInterpolation
{
    typedef VectorEpetra                             vector_Type;
    typedef std::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef MapEpetra                                map_Type;
    typedef std::shared_ptr<map_Type>              mapPtr_Type;
    typedef FESpace< mesh_Type, map_Type >           fespace_Type;
    typedef std::shared_ptr< fespace_Type >        fespacePtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef std::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef std::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;
    typedef std::shared_ptr< NavierStokesProblem<mesh_Type> > NSProblemPtr_Type;

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

template< class mesh_Type >
void
InitPolicyInterpolation< mesh_Type >::
setupInit ( Teuchos::ParameterList& /*list*/ )
{

}

template< class mesh_Type >
void
InitPolicyInterpolation< mesh_Type >::
initSimulation ( bcContainerPtr_Type /*bchandler*/,
                 vectorPtr_Type solution )
{
    ASSERT ( problem()->hasExactSolution(), "Error: You cannot use the interpolation method if the problem has not an analytical solution." );

    Real currentTime = initialTime() - timestep() * bdf()->order();
    UInt pressureOffset = uFESpace()->fieldDim() * uFESpace()->dof().numTotalDof();

    vectorPtr_Type velocity;
    velocity.reset ( new vector_Type ( uFESpace()->map(), Unique ) );

    vectorPtr_Type pressure;
    pressure.reset ( new vector_Type ( pFESpace()->map(), Unique ) );

    *solution = 0.0;
    uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
    pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
    solution->add ( *velocity );
    solution->add ( *pressure, pressureOffset );
    bdf()->setInitialCondition ( *solution );

    currentTime += timestep();
    for ( ; currentTime <=  initialTime() + timestep() / 2.0; currentTime += timestep() )
    {
        *solution = 0.0;
        *velocity = 0.0;
        *pressure = 0.0;

        uFESpace()->interpolate ( problem()->uexact(), *velocity, currentTime );
        pFESpace()->interpolate ( problem()->pexact(), *pressure, currentTime );
        solution->add ( *velocity );
        solution->add ( *pressure, pressureOffset );

        // Updating bdf
        bdf()->shiftRight ( *solution );
    }
}

} // namespace LifeV

#endif /* INITPOLICYINTERPOLATION_HPP */
