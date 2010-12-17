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
    @brief

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 08-10-2010

 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include <mpi.h>

#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/hdf5exporter.hpp>

#include <life/lifefem/bcManage.hpp>

#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#include <life/lifemesh/regionMesh3D.hpp>

#include <life/lifesolver/LevelSetSolver.hpp>
#include <life/lifemesh/dataMesh.hpp>

using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}


Real betaFct( const Real& /* t */, const Real& /* x */, const Real& /* y */, const Real& /* z */, const ID& i )
{
    if (i == 1)
        return 1;
    return 0;
}

Real initLSFct( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */)
{
    return 	sqrt((x+0.7)*(x+0.7) + y*y + z*z ) - 0.3;
}

Real exactSolution( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /* i */)
{
    return 	sqrt((x+0.7-t)*(x+0.7-t) + y*y + z*z ) - 0.3;
}


typedef RegionMesh3D<LinearTetra> mesh_type;
typedef EpetraMatrix<Real> matrix_type;
typedef EpetraVector vector_type;

int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    boost::shared_ptr<Epetra_Comm> Comm(new Epetra_SerialComm);
#endif

    const bool verbose(Comm->MyPID()==0);

// Read first the data needed

    if (verbose) std::cout << " -- Reading the data ... " << std::flush;
    GetPot dataFile( "data" );
    if (verbose) std::cout << " done ! " << std::endl;

    const UInt Nelements(dataFile("mesh/nelements",20));
    if (verbose) std::cout << " ---> Number of elements : " << Nelements << std::endl;

// Build and partition the mesh

    if (verbose) std::cout << " -- Building the mesh ... " << std::flush;
    boost::shared_ptr< mesh_type > fullMeshPtr(new RegionMesh3D<LinearTetra>);
    regularMesh3D( *fullMeshPtr, 1, Nelements, Nelements, Nelements, false,
                   2.0,   2.0,   2.0,
                   -1.0,  -1.0,  -1.0);
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Partitioning the mesh ... " << std::flush;
    partitionMesh< mesh_type >   meshPart(fullMeshPtr, Comm);
    if (verbose) std::cout << " done ! " << std::endl;

    if (verbose) std::cout << " -- Freeing the global mesh ... " << std::flush;
    fullMeshPtr.reset();
    if (verbose) std::cout << " done ! " << std::endl;

// Build the FESpaces

    if (verbose) std::cout << " -- Building FESpaces ... " << std::flush;
    std::string uOrder("P1");
    std::string bOrder("P1");
    boost::shared_ptr<FESpace< mesh_type, EpetraMap > > uFESpace( new FESpace< mesh_type, EpetraMap >(meshPart,uOrder, 1, Comm));
    boost::shared_ptr<FESpace< mesh_type, EpetraMap > > betaFESpace( new FESpace< mesh_type, EpetraMap >(meshPart,bOrder, 3, Comm));
    if (verbose) std::cout << " done ! " << std::endl;
    if (verbose) std::cout << " ---> Dofs: " << uFESpace->dof().numTotalDof() << std::endl;

// Read the data

    boost::shared_ptr<DataLevelSet> data_level_set(new DataLevelSet);
    data_level_set->setup(dataFile,"level-set");

// Build the solver

    LevelSetSolver<mesh_type> level_set(uFESpace,betaFESpace);
    level_set.setup(data_level_set);

    vector_type initLS(level_set.map());
    uFESpace->interpolate(initLSFct,initLS,0.0);
    level_set.initialize(initLS);

    level_set.setupLinearSolver(dataFile,"");

    vector_type beta(betaFESpace->map());
    betaFESpace->interpolate(betaFct,beta,0.0);

    BCHandler bchandler(26);
    BCFunctionBase BCu( exactSolution );
    bchandler.addBC("Dirichlet",1,Essential,Full,BCu,1);
    for (UInt i(2); i<=26; ++i)
    {
        bchandler.addBC("Dirichlet",i,Essential,Full,BCu,1);
    }

    Hdf5exporter<mesh_type> exporter ( dataFile, meshPart.meshPartition(), "solution", Comm->MyPID());
    exporter.setMultimesh(false);
    boost::shared_ptr<vector_type> solutionPtr (new vector_type(level_set.solution(),Repeated));
    exporter.addVariable( ExporterData::Scalar, "level-set", solutionPtr, UInt(0), uFESpace->dof().numTotalDof() );
    exporter.postProcess(0);

    Real current_time( data_level_set->dataTime()->getInitialTime() );
    Real dt( data_level_set->dataTime()->getTimeStep() );
    Real final_time( data_level_set->dataTime()->getEndTime() );

    while ( current_time < final_time)
    {
        current_time += dt;
        data_level_set->dataTime()->updateTime();
        if (verbose) std::cout << " We are now at time " << current_time << std::endl;

        if (verbose) std::cout << "[LS] Building system ... " << std::flush;
        level_set.updateSystem(beta,bchandler,current_time);
        if (verbose) std::cout << " done " << std::endl;
        if (verbose) std::cout << "[LS] Solving system ... " << std::flush;
        Chrono c;
        c.start();
        level_set.iterate();
        c.stop();
        if (verbose) std::cout << " Iterate done in " << c.diff() << std::endl;
        if (verbose) std::cout << "[LS] Reinitizialization ... " << std::flush;
        level_set.reinitializationDirect();
        if (verbose) std::cout << " done " << std::endl;
        if (verbose) std::cout << "[LS] Exporting ... " << std::flush;
        *solutionPtr = level_set.solution();
        exporter.postProcess(current_time);
        if (verbose) std::cout << " done " << std::endl;
    }

    Real N(level_set.solution().norm1());
    if (verbose) std::cout << "Final norm of the solution : " << N << std::endl;

    if ((N < 6900) || (N>7100))
    {
        return (EXIT_FAILURE);
    }

    exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}


