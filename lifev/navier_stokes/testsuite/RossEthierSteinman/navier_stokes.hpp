/* -*- mode: c++ -*-
 
 This file is part of the LifeV library.
 Copyright (C) 2010 EPFL
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/*!
 @file navierStokes.hpp
 @author Davide Forti <davide.forti@epfl.ch>
 @date 2014-02-06
 */

#ifndef NAVIERSTOKES_H
#define NAVIERSTOKES_H 1

#include <lifev/navier_stokes/solver/OseenSolver.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/navier_stokes/solver/OseenData.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/navier_stokes/fem/TimeAdvanceBDFNavierStokes.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>
#include <lifev/core/filter/ExporterVTK.hpp>
#include <lifev/core/filter/ExporterEmpty.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/filter/PartitionIO.hpp>

// Include for boundary conditions
#include "functions.hpp"

using namespace LifeV;

class NavierStokes
{

public:

    typedef RegionMesh<LinearTetra>                       mesh_Type;
    
    typedef boost::shared_ptr<mesh_Type >                 meshPtr_Type;
    
    typedef FESpace< mesh_Type, MapEpetra >               feSpace_Type;
    
    typedef boost::shared_ptr<feSpace_Type>               feSpacePtr_Type;
    
    typedef OseenSolver< mesh_Type >                      fluid_Type;
    
    typedef VectorEpetra                                  vector_Type;
    
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
    
    typedef boost::shared_ptr< Exporter<mesh_Type > >     exporterPtr_Type;
    
    NavierStokes ( int argc,
                   char** argv,
                   boost::shared_ptr<Epetra_Comm> Comm,
                   const std::string defaultDataName = "data",
                   const std::string outputName = "result");
    
    ~NavierStokes(){}
    
    void run();
    
private:
    
    boost::shared_ptr<Epetra_Comm>       M_comm;
    boost::shared_ptr<GetPot>            M_dataFile;
    std::string                          M_outputFilename;
};


using namespace LifeV;

NavierStokes::NavierStokes ( int argc,
                            char** argv,
                            boost::shared_ptr<Epetra_Comm> Comm,
                            const std::string defaultDataName,
                            const std::string outputName):
M_comm (Comm)
{
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow (defaultDataName.c_str(), 2, "-f", "--file");
    M_dataFile.reset(new GetPot ( data_file_name ) );
    M_outputFilename = outputName;
}

void NavierStokes::run()
{
    //////////////////////////////////////////////////////////////
    // Flag to set variable for parallel output on the terminal //
    //////////////////////////////////////////////////////////////
    
    bool verbose = (M_comm->MyPID() == 0);
    
    //////////////////////
    // Loading the mesh //
    //////////////////////
    
    MeshData meshData;
    meshPtr_Type fullMeshPtr ( new mesh_Type ( M_comm ) );
    meshData.setup (*M_dataFile, "fluid/space_discretization");
    if (verbose){
        std::cout << "\n\n[Reading the mesh " << meshData.meshFile() << " ]\n\n";
    }
    readMesh (*fullMeshPtr, meshData);
    
//    meshPtr_Type fullMeshPtr ( new mesh_Type ( M_comm ) );
//    regularMesh3D( *fullMeshPtr, 1, 8, 8, 8, false, 2.0, 2.0, 2.0, -1.0,  -1.0,  -1.0);
    
    ///////////////////////////
    // Partitioning the mesh //
    ///////////////////////////
    
    meshPtr_Type localMeshPtr;
    if (verbose){
        std::cout << "\n[Partioning the mesh ]\n\n";
    }
    MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, M_comm);
    localMeshPtr = meshPart.meshPartition();
    fullMeshPtr.reset();
    
    ///////////////////////////
    // Finite element spaces //
    ///////////////////////////
    
    std::string uOrder = (*M_dataFile)("fluid/space_discretization/vel_order","P1");
    std::string pOrder = (*M_dataFile)("fluid/space_discretization/pres_order","P1");
    
    feSpacePtr_Type uFESpace;
    uFESpace.reset (new feSpace_Type (localMeshPtr, uOrder, 3, M_comm) );
    
    feSpacePtr_Type pFESpace;
    pFESpace.reset (new feSpace_Type (localMeshPtr, pOrder, 1, M_comm) );
    
    UInt totalPressDof = pFESpace->dof().numTotalDof();
    UInt pressureOffset =  uFESpace->fieldDim() * uFESpace->dof().numTotalDof();
    
    if (verbose)
        std::cout << "\n\n[Total Velocity Dof = " << pressureOffset << " ]\n";
    
    if (verbose)
        std::cout << "\n[Total Pressure Dof = " << totalPressDof << " ]\n";
    
    /////////////////////////////////////////
    // Boundary conditions for the problem //
    /////////////////////////////////////////
    
    BCFunctionBase u_dirichlet( exactVelocity );
    BCFunctionBase u_neumann( normalStress );
    BCHandler bcH;
    
    if (verbose)
        std::cout << "\n[Setting boundary conditions ] \n";
    
    bcH.addBC( "S", 5, Essential, Full, u_dirichlet, 3 );
    bcH.addBC( "S", 1, Essential, Full, u_dirichlet, 3 );
    bcH.addBC( "S", 3, Essential, Full, u_dirichlet, 3 );
    bcH.addBC( "S", 4, Essential, Full, u_dirichlet, 3 );
    bcH.addBC( "S", 6, Essential, Full, u_dirichlet, 3 );
    bcH.addBC( "S", 2, Natural,   Full, u_neumann,   3 );
    bcH.bcUpdate ( *localMeshPtr, uFESpace->feBd(), uFESpace->dof() );
    
    //////////////////////////
    // Setting the exporter //
    //////////////////////////
    
    std::string const exporterType =  (*M_dataFile)( "exporter/type", "hdf5");
    exporterPtr_Type exporter;
    
    if (verbose)
        std::cout << "\n[Setting the exporter ] \n";
    
    if (exporterType.compare ("hdf5") == 0){
        exporter.reset ( new ExporterHDF5<mesh_Type > ( *M_dataFile, M_outputFilename ) );
        exporter->setMeshProcId ( localMeshPtr, M_comm->MyPID() );
    }
    else if(exporterType.compare ("vtk") == 0){
        exporter.reset ( new ExporterVTK<mesh_Type > ( *M_dataFile, M_outputFilename ) );
        exporter->setMeshProcId ( localMeshPtr, M_comm->MyPID() );
    }    
    
    ///////////////////////////
    // Initialize the solver //
    ///////////////////////////
    
    if (verbose)
        std::cout << "\n[Initializing the solver ] \n";
    
    boost::shared_ptr<OseenData> oseenData (new OseenData() );
    oseenData->setup ( (*M_dataFile) );
    
    OseenSolver< mesh_Type > fluid (oseenData, *uFESpace, *pFESpace, M_comm);
    MapEpetra fullMap (fluid.getMap() );
    
    //////////////////////////////////////////////////
    // Vectors that will contain the exact solution //
    //////////////////////////////////////////////////
    
    vectorPtr_Type exactPressPtr ( new vector_Type(pFESpace->map(), Repeated) );
    vectorPtr_Type exactVelPtr   ( new vector_Type(uFESpace->map(), Repeated) );
    pFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactPressure ), *exactPressPtr, 0 );
    uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactVelocity ), *exactVelPtr, 0 );
    
    /////////////////////////////////////////////////////
    // Vector that will contain the numerical solution //
    /////////////////////////////////////////////////////
    
    vectorPtr_Type velAndPressure ( new vector_Type (fullMap, exporter->mapType() ) );
    
    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "u",       uFESpace, velAndPressure, UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "p",       pFESpace, velAndPressure, pressureOffset );
    exporter->addVariable ( ExporterData<mesh_Type>::VectorField, "u_exact", uFESpace, exactVelPtr,    UInt (0) );
    exporter->addVariable ( ExporterData<mesh_Type>::ScalarField, "p_exact", pFESpace, exactPressPtr,  UInt (0) );
    exporter->postProcess ( 0 );
    
    
    ///////////////////////////////
    // Initialize the simulation //
    ///////////////////////////////
    
    // Building the constant matrices of Navier-Stokes
    fluid.setUp (*M_dataFile);
    fluid.buildSystem();
    
    Real dt     = oseenData->dataTime()->timeStep();
    Real t0     = oseenData->dataTime()->initialTime();
    Real tFinal = oseenData->dataTime()->endTime();
    
    oseenData->dataTime()->setTime (t0);
    
    TimeAdvanceBDFNavierStokes<vector_Type> bdf;
    bdf.setup (oseenData->dataTimeAdvance()->orderBDF() );
    
    // For the moment just using interpolation
    std::vector<vectorPtr_Type> solutionStencil;
    solutionStencil.resize ( bdf.bdfVelocity().size() );
    for ( UInt i(0); i < bdf.bdfVelocity().size() ; ++i){
        *exactVelPtr    *= 0;
        *velAndPressure *= 0;
        uFESpace->interpolate ( static_cast<feSpace_Type::function_Type> ( exactVelocity ), *exactVelPtr, t0-i*dt );
        velAndPressure->subset(*exactVelPtr, uFESpace->map(), 0, 0);
        solutionStencil[ i ] = velAndPressure;
    }
    
    bdf.bdfVelocity().setInitialCondition (solutionStencil);
    
    Real alphaOverDt = bdf.bdfVelocity().coefficientFirstDerivative ( 0 ) / oseenData->dataTime()->timeStep();

    vector_Type u_star ( fullMap );                 // extrapolated velocity
    
    vector_Type rhs    ( fullMap );                 // right hand side
    
    int iter = 1;
    
    Real time = t0 + dt;
    
    for ( ; time <= tFinal + dt / 2.; time += dt, iter++){
        
        if (verbose)
            std::cout << "\n  [Solving for time = " << time << " ] \n";
        
        oseenData->dataTime()->setTime (time);
        
        // extrapolate the velocity
        bdf.bdfVelocity().extrapolation (u_star);
        
        // update the right hand side term
        bdf.bdfVelocity().updateRHSContribution(dt);

        // compute the right hand side, it contains just the term associated to the previous timesteps
        // CAUTION: the vector bdf.bdfVelocity().rhsContributionFirstDerivative() is already divided by the timestep
        rhs  = fluid.matrixMass() * bdf.bdfVelocity().rhsContributionFirstDerivative();
        
        // update the convective term of the Navier Stokes matrix and set the right hand side in the solver
        fluid.updateSystem ( alphaOverDt, u_star, rhs );
        
        // apply boundary conditions and solving the linear system
        fluid.iterate ( bcH );
        
        // shifting right the computed solution in the time-advance stencil
        bdf.bdfVelocity().shiftRight ( *fluid.solution() );
        
        // copy the solution for the postprocessing
        *velAndPressure = *fluid.solution();
        
        // evaluate the exact fluid velocity and pressure at the current time
        pFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( exactPressure ), *exactPressPtr, time );
        uFESpace->interpolate ( static_cast<typename feSpace_Type::function_Type> ( exactVelocity ), *exactVelPtr, time );
        
        // postprocessing
        exporter->postProcess ( time );
    
    }
}

#endif /* NAVIERSTOKES_H */
