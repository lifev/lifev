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

@author Davide Forti <davide.forti@epfl.ch>
@date 05-28-2015
*/

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

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/interpolation/Interpolation.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

using namespace LifeV;

double fx (double x, double y, double z)
{
    return sin (x) + y + z;
}

double fy (double x, double y, double z)
{
    return sin (y) + x + z;
}

double fz (double x, double y, double z)
{
    return sin (z) + x + y ;
}

Real exact_sol (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
	switch (i)
	{
	case 0:
		return sin (x) + y + z;
		break;
	case 1:
		return sin (y) + x + z;
		break;
	case 2:
		return sin (z) + x + y ;
		break;
	default:
		ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
		return 0.;
		break;
	}
}

Real exact_sol_easy (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
	switch (i)
	{
	case 0:
		return 20.;
		break;
	case 1:
		return 30.;
		break;
	case 2:
		return -15.;
		break;
	default:
		ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
		return 0.;
		break;
	}
}

int main (int argc, char** argv )
{
    typedef VectorEpetra                          vector_Type;
    typedef boost::shared_ptr<vector_Type >       vectorPtr_Type;
    typedef RegionMesh<LinearTetra >              mesh_Type;
    typedef boost::shared_ptr< mesh_Type >        meshPtr_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > FESpace_Type;

    boost::shared_ptr<Epetra_Comm> Comm;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    Comm.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    Comm.reset ( new Epetra_SerialComm() );
#endif

    // DATAFILE
    GetPot command_line (argc, argv);
    GetPot dataFile ( command_line.follow ("data_rbf3d", 2, "-f", "--file" ) );

    // BELOS XML file
    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList_rbf3d.xml" );

    // LOADING MESHES
    MeshData Solid_mesh_data;
    meshPtr_Type Solid_mesh_ptr ( new mesh_Type ( Comm ) );
    Solid_mesh_data.setup (dataFile, "solid/space_discretization");
    readMesh (*Solid_mesh_ptr, Solid_mesh_data);

    MeshData Fluid_mesh_data;
    meshPtr_Type Fluid_mesh_ptr ( new mesh_Type ( Comm ) );
    Fluid_mesh_data.setup (dataFile, "fluid/space_discretization");
    readMesh (*Fluid_mesh_ptr, Fluid_mesh_data);

    // PARTITIONING MESHES
    MeshPartitioner<mesh_Type>   Solid_mesh_part;
    boost::shared_ptr<mesh_Type> Solid_localMesh;
    Solid_mesh_part.doPartition (Solid_mesh_ptr, Comm);
    Solid_localMesh = Solid_mesh_part.meshPartition();

    MeshPartitioner<mesh_Type>   Fluid_mesh_part;
    boost::shared_ptr<mesh_Type> Fluid_localMesh;
    Fluid_mesh_part.doPartition (Fluid_mesh_ptr, Comm);
    Fluid_localMesh = Fluid_mesh_part.meshPartition();

    // ORDINE SPAZI FE
    std::string orderFluid = dataFile("fluid/space_discretization/order","default");
    std::string orderStructure = dataFile("solid/space_discretization/order","default");

    // CREO SPAZI ELEMENTI FINITI
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Solid_fieldFESpace_whole (new FESpace<mesh_Type, MapEpetra> (Solid_mesh_ptr, orderStructure, 3, Comm) );
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Fluid_fieldFESpace_whole (new FESpace<mesh_Type, MapEpetra> (Fluid_mesh_ptr, orderFluid, 3, Comm) );

    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Solid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Solid_localMesh, orderStructure, 3, Comm) );
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Fluid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Fluid_localMesh, orderFluid, 3, Comm) );

    Interpolation intergrid;
    intergrid.setup(dataFile, belosList);
    intergrid.setFlag(1);

    // Set vectors for interpolation
    vectorPtr_Type Solid_vector (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Fluid_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    intergrid.setVectors(Solid_vector, Fluid_solution);

    // Build table of DOFs
    intergrid.buildTableDofs_known ( Solid_fieldFESpace_whole );
    intergrid.buildTableDofs_unknown ( Fluid_fieldFESpace_whole );

    // Build table of DOFs
    intergrid.identifyNodes_known ( Solid_fieldFESpace );
    intergrid.identifyNodes_unknown ( Fluid_fieldFESpace );

    // build maps
    intergrid.buildKnownInterfaceMap();
    intergrid.buildUnknownInterfaceMap();

    intergrid.buildOperators();
    
    vectorPtr_Type Solid_vector_one (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    Solid_fieldFESpace->interpolate ( static_cast<FESpace_Type::function_Type> ( exact_sol_easy ), *Solid_vector_one, 0.0 );

    // GET THE EXACT SOLUTION THAT WILL BE INTERPOLATED
    Solid_fieldFESpace->interpolate ( static_cast<FESpace_Type::function_Type> ( exact_sol ), *Solid_vector, 0.0 );

    // TESTING THE RESTRICTION AND EXPANDING METHODS
    vectorPtr_Type solidVectorOnGamma;
    vectorPtr_Type solidVectorOnOmega;

    intergrid.restrictOmegaToGamma_Known(Solid_vector, solidVectorOnGamma);
    intergrid.expandGammaToOmega_Known(solidVectorOnGamma, solidVectorOnOmega);

    // EXPORTING THE DEFINED FIELD
    ExporterHDF5<mesh_Type> Solid_exporter (dataFile, Solid_localMesh, "InputField", Comm->MyPID() );
    Solid_exporter.setMeshProcId (Solid_localMesh, Comm->MyPID() );
    Solid_exporter.exportPID (Solid_localMesh, Comm, true );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "f(x,y,z)", Solid_fieldFESpace, Solid_vector, UInt (0) );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "restrict-expand", Solid_fieldFESpace, solidVectorOnOmega, UInt (0) );
    Solid_exporter.postProcess (0);
    Solid_exporter.closeFile();

    // TESTING THE METHOD updateRhs
    intergrid.updateRhs (Solid_vector);

    // PERFORMING INTERPOLATION WITH A VECTOR OF ONE
    intergrid.interpolate();

    // GET THE SOLUTION
    intergrid.solution (Fluid_solution);
    
    Fluid_solution->spy("Fluid_solution");

    // GET THE SOLUTION ON GAMMA
    vectorPtr_Type Fluid_solutionOnGamma;
    intergrid.getSolutionOnGamma (Fluid_solutionOnGamma);

    Fluid_solutionOnGamma->spy("Fluid_solutionOnGamma");
    
    // EXPORTING THE SOLUTION
    ExporterHDF5<mesh_Type> Fluid_exporter (dataFile, Fluid_localMesh, "Results", Comm->MyPID() );
    Fluid_exporter.setMeshProcId (Fluid_localMesh, Comm->MyPID() );
    Fluid_exporter.exportPID (Fluid_localMesh, Comm, true );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Solution", Fluid_fieldFESpace, Fluid_solution, UInt (0) );
    Fluid_exporter.postProcess (0);
    Fluid_exporter.closeFile();
    
//    std::cout << "aaaaaaaaaaa";
//
//    for (int i = 0; i < 3000000; ++i )
//    	for (int j = 0; j < 3000000; ++j )
//    		int a = i + j;




#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    return (EXIT_SUCCESS);

}
