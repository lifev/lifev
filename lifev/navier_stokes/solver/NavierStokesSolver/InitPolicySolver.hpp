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
    @file InitPolicySolver class
    @brief This class is a strategy to initialize a Navier-Stokes problem

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 11-12-2012
 */

#ifndef INITPOLICYSOLVER_HPP
#define INITPOLICYSOLVER_HPP

#include <iostream>
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
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/navier_stokes/solver/NavierStokesSolver/NavierStokesProblem.hpp>


namespace LifeV
{

template< class mesh_Type, class TimeIterationPolicy >
struct InitPolicySolver : public virtual TimeIterationPolicy
{
    typedef VectorEpetra                             vector_Type;
    typedef std::shared_ptr<VectorEpetra>          vectorPtr_Type;
    typedef MeshPartitioner< mesh_Type >             meshPartitioner_Type;
    typedef MapEpetra                                map_Type;
    typedef std::shared_ptr<map_Type>              mapPtr_Type;
    typedef TimeAdvanceBDF<vector_Type>              bdf_Type;
    typedef std::shared_ptr< bdf_Type >            bdfPtr_Type;
    typedef BCHandler                                bcContainer_Type;
    typedef std::shared_ptr<bcContainer_Type>      bcContainerPtr_Type;
    typedef std::shared_ptr< NavierStokesProblem<mesh_Type> > NSProblemPtr_Type;

    InitPolicySolver() {}
    virtual ~InitPolicySolver() {}

    void setupInit ( Teuchos::ParameterList& list );

    void initSimulation ( bcContainerPtr_Type bchandler,
                          vectorPtr_Type solution );

    virtual Displayer displayer() = 0;
    virtual Real timestep() const = 0;
    virtual Real initialTime() const = 0;
    virtual bdfPtr_Type bdf() const = 0;
};

template< class mesh_Type, class TimeIterationPolicy >
void
InitPolicySolver< mesh_Type, TimeIterationPolicy >::
setupInit ( Teuchos::ParameterList& list )
{
    TimeIterationPolicy::initTimeIteration ( list );
}

template< class mesh_Type, class TimeIterationPolicy >
void
InitPolicySolver< mesh_Type, TimeIterationPolicy >::
initSimulation ( bcContainerPtr_Type bchandler,
                 vectorPtr_Type solution )
{
    displayer().leaderPrint ( "\n[Initializing the problem]\n" );

    LifeChrono nsTimeLoopChrono;
    nsTimeLoopChrono.start();

    Real currentTime = initialTime() - timestep() * ( bdf()->order() - 1 );

    *solution = 0.0;
    bdf()->setInitialCondition ( *solution );

    while ( currentTime <= initialTime() + timestep() / 2.0 )
    {
        LifeChrono iterChrono;
        iterChrono.start();

        displayer().leaderPrint ( "\n[t = ", currentTime, " s.]\n" );

        // if( !M_usePreviousSolutionAsGuess ) *solution = 0;

        TimeIterationPolicy::iterate ( solution,
                                       bchandler,
                                       currentTime );

        // Updating the BDF scheme
        bdf()->shiftRight ( *solution );

        iterChrono.stop();
        displayer().leaderPrintMax ( "Iteration time: ", iterChrono.diff(), " s.\n" );

        currentTime += timestep();

#ifdef HAVE_MPI
        MPI_Barrier ( MPI_COMM_WORLD );
#endif
    }

    nsTimeLoopChrono.stop();
    displayer().leaderPrintMax ("Time of the temporal loop: ", nsTimeLoopChrono.diff(), " s.\n");
}

} // namespace LifeV

#endif /* INITPOLICYSTOKES_HPP */
