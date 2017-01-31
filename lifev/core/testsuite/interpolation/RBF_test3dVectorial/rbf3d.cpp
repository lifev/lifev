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
#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>
#include <lifev/core/interpolation/RBFrescaledVectorial.hpp>
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

Real exact_sol_easy (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
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
    typedef std::shared_ptr<vector_Type >       vectorPtr_Type;
    typedef RegionMesh<LinearTetra >              mesh_Type;
    typedef std::shared_ptr< mesh_Type >        meshPtr_Type;
    typedef RBFInterpolation<mesh_Type>           interpolation_Type;
    typedef std::shared_ptr<interpolation_Type> interpolationPtr_Type;
    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra > FESpace_Type;

    std::shared_ptr<Epetra_Comm> Comm;
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
    std::shared_ptr<mesh_Type> Solid_localMesh;
    Solid_mesh_part.doPartition (Solid_mesh_ptr, Comm);
    Solid_localMesh = Solid_mesh_part.meshPartition();

    MeshPartitioner<mesh_Type>   Fluid_mesh_part;
    std::shared_ptr<mesh_Type> Fluid_localMesh;
    Fluid_mesh_part.doPartition (Fluid_mesh_ptr, Comm);
    Fluid_localMesh = Fluid_mesh_part.meshPartition();

    // CREATING A FE-SPACE FOR THE GRID ON WHICH WE ASSUME TO KNOW THE INTERFACE FIELD. DEFINING AN INTERFACE VECTOR TO BE INTERPOLATED.
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > Solid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Solid_localMesh, "P1", 3, Comm) );
    vectorPtr_Type Solid_vector (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    vectorPtr_Type Solid_vector_one (new vector_Type (Solid_fieldFESpace->map(), Unique) );
    Solid_fieldFESpace->interpolate ( static_cast<FESpace_Type::function_Type> ( exact_sol_easy ), *Solid_vector_one, 0.0 );

    // GET THE EXACT SOLUTION THAT WILL BE INTERPOLATED
    Solid_fieldFESpace->interpolate ( static_cast<FESpace_Type::function_Type> ( exact_sol ), *Solid_vector, 0.0 );

    // CREATING A FE-SPACE FOR THE GRID ON WHICH WE WANT TO INTERPOLATE THE DATA. INITIALIZING THE SOLUTION VECTOR.
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > Fluid_fieldFESpace (new FESpace<mesh_Type, MapEpetra> (Fluid_localMesh, "P1", 3, Comm) );
    vectorPtr_Type Fluid_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );

    // NUMBER OF FLAGS CONSIDERED: THE DOFS WHOSE FLAG IS 1 AND 20 ARE TAKEN INTO ACCOUNT
    // NOTE: from mesh to mesh the vector of flag has to contain one element with value -1
    int nFlags = 2;
    std::vector<int> flags (nFlags);
    flags[0] = 1;
    flags[1] = 20;
    
    // INITIALIZE THE INTERPOLANT
    interpolationPtr_Type RBFinterpolant;

    // READING FROM DATAFILE WHICH INTERPOLANT HAS TO BE USED
    RBFinterpolant.reset ( interpolation_Type::InterpolationFactory::instance().createObject (dataFile("interpolation/interpolation_Type","none")));

    // SET UP FOR THE INTERPOLANT: GLOBAL AND LOCAL MESHES PLUS THE VECTOR OF FLAGS
    RBFinterpolant->setup(Solid_mesh_ptr, Solid_localMesh, Fluid_mesh_ptr, Fluid_localMesh, flags);
    
    // PASSING THE VECTOR WITH THE DATA, THE ONE THAT WILL STORE THE SOLUTION, THE DATAFILE AND THE BELOS LIST
    RBFinterpolant->setupRBFData (Solid_vector_one, Fluid_solution, dataFile, belosList);

    // CREATING THE RBF OPERATORS
    RBFinterpolant->buildOperators();

    // TESTING THE RESTRICTION AND EXPANDING METHODS
    vectorPtr_Type solidVectorOnGamma;
    vectorPtr_Type solidVectorOnOmega;

    RBFinterpolant->restrictOmegaToGamma_Known(Solid_vector, solidVectorOnGamma);
    RBFinterpolant->expandGammaToOmega_Known(solidVectorOnGamma, solidVectorOnOmega);

    // EXPORTING THE DEFINED FIELD
    ExporterHDF5<mesh_Type> Solid_exporter (dataFile, Solid_localMesh, "InputField", Comm->MyPID() );
    Solid_exporter.setMeshProcId (Solid_localMesh, Comm->MyPID() );
    Solid_exporter.exportPID (Solid_localMesh, Comm, true );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "f(x,y,z)", Solid_fieldFESpace, Solid_vector, UInt (0) );
    Solid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "restrict-expand", Solid_fieldFESpace, solidVectorOnOmega, UInt (0) );
    Solid_exporter.postProcess (0);
    Solid_exporter.closeFile();

    // PERFORMING INTERPOLATION WITH A VECTOR OF ONE
    RBFinterpolant->interpolate();

    // GET THE SOLUTION
    RBFinterpolant->solution (Fluid_solution);
    
    // GET THE SOLUTION ON GAMMA
    vectorPtr_Type Fluid_solutionOnGamma;
    RBFinterpolant->getSolutionOnGamma (Fluid_solutionOnGamma);

    // TESTING THE METHOD updateRhs
    RBFinterpolant->updateRhs (Solid_vector);
    
    // PERFORMING INTERPOLATION WITH A VECTOR OF ONE
    RBFinterpolant->interpolate();
    
    // GET THE SOLUTION
    RBFinterpolant->solution (Fluid_solution);

    // COMPUTING THE ERROR
    vectorPtr_Type Fluid_exact_solution (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type myError (new vector_Type (Fluid_fieldFESpace->map(), Unique) );
    vectorPtr_Type rbfError (new vector_Type (Fluid_fieldFESpace->map(), Unique) );

    // NOTE: HERE WE DO NOT USE THE INTERPOLATE METHOD OF THE FESPACE CLASS SINCE WE WANT THAT THE EXACT SOLUTION
    // IS DEFINED JUST AT THE INTERFACE BETWEEN THE MESHES
    for ( UInt j = 0; j < 3; ++j)
        for ( UInt i = 0; i < Fluid_exact_solution->epetraVector().MyLength(); ++i )
            if ( Fluid_exact_solution->blockMap().GID (i) < Fluid_mesh_ptr->pointList.size())
                if ( Fluid_mesh_ptr->point ( Fluid_exact_solution->blockMap().GID (i) ).markerID() == flags[0] ||
                     Fluid_mesh_ptr->point ( Fluid_exact_solution->blockMap().GID (i) ).markerID() == flags[1] )
                    if ( Fluid_exact_solution->blockMap().LID ( static_cast<EpetraInt_Type> ( Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ) ) != -1)
                    {
                        switch (j)
                        {
                        case(0):
                            (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fx ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                            break;
                        case(1):
                            (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fy ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                            break;
                        case(2):
                            (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = fz ( Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).x(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).y(), Fluid_mesh_ptr->point (Fluid_exact_solution->blockMap().GID (i) ).z() );
                            break;
                        }

                        (*myError) [myError->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] = (*Fluid_exact_solution) [Fluid_exact_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ] - (*Fluid_solution) [Fluid_solution->blockMap().GID (i) + Fluid_exact_solution->size()/3*j ];
                    }

    // COMPUTING SOME ERRORS
    Real err_inf  = myError->normInf();

    // EXPORTING THE SOLUTION
    ExporterHDF5<mesh_Type> Fluid_exporter (dataFile, Fluid_localMesh, "Results", Comm->MyPID() );
    Fluid_exporter.setMeshProcId (Fluid_localMesh, Comm->MyPID() );
    Fluid_exporter.exportPID (Fluid_localMesh, Comm, true );
    Fluid_exporter.addVariable (ExporterData<mesh_Type>::VectorField, "Solution", Fluid_fieldFESpace, Fluid_solution, UInt (0) );
    Fluid_exporter.postProcess (0);
    Fluid_exporter.closeFile();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( std::abs(err_inf) < 1.0e-1 )
    {
        return ( EXIT_SUCCESS );
    }
    else
    {
        return ( EXIT_FAILURE );
    }

}
