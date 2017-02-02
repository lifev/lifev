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
    @brief Interpolate test

    @author Mauro Perego <mperego@fsu.edu>
    @contributor
    @maintainer Mauro Perego <mperego@fsu.edu>

    @date 06-30-2010

The program tests the interpolation methods between different finite elements (mostly between scalar continuous finite elements).
Also it test the interpolation of an analytical function into a finite element space.
 */


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif



#include <lifev/core/fem/ReferenceFE.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>
#include <lifev/core/fem/CurrentFE.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/fem/DOF.hpp>
#include <lifev/core/filter/MeshWriter.hpp>

#include <lifev/core/mesh/MeshData.hpp>


#include "test_interpolate.hpp"

using namespace LifeV;

int main (int argc, char** argv )
{
    typedef FESpace < RegionMesh<LinearTetra>, MapEpetra > FESpaceTetra_Type;
    typedef std::shared_ptr < FESpaceTetra_Type > FESpaceTetraPtr_Type;

    typedef FESpace < RegionMesh<LinearHexa>, MapEpetra > FESpaceHexa_Type;
    typedef std::shared_ptr < FESpaceHexa_Type > FESpaceHexaPtr_Type;


#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    std::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    UInt verbose = (Comm->MyPID() == 0);

    bool check (true);


    //definition of Array of Errors which will be used to check the correctness of the interpolate test.

    const Real errArrayBilinear[2] = { 0.0312819802, 0. };

    const Real errArrayQuadratic[12] = {    0.0136247667, 0.0005088372, 0.0005577494, 0.0005088372,
                                            0.0136172446, 0.0005088372, 0.0004270717, 0.0005088372,
                                            0.0136172446, 0.0005088372, 0.0004270717,           0.
                                       };

    const Real errArrayBubble[12] = {   0.0094702745, 3.584186e-10, 3.67611e-10,  3.584186e-10,
                                        0.0094702745, 3.584186e-10,           0., 3.584186e-10,
                                        0.0094702745, 3.584186e-10, 3.67611e-10,  3.584186e-10
                                    };

    const Real errArrayLinear[12] = {   0.010437463587, 0., 0., 0.,
                                        0.010437463587, 0., 0., 0.,
                                        0.010437463587, 0., 0., 0.
                                    };


    const std::string stringArrayP[12] = {  "P1  -> P0 ",  "P1  -> P1 ",  "P1  -> P1b", "P1  -> P2 ",
                                            "P1b -> P0 ",  "P1b -> P1 ",  "P1b -> P1b", "P1b -> P2 ",
                                            "P2  -> P0 ",  "P2  -> P1 ",  "P2  -> P1b", "P2  -> P2 "
                                         };

    const std::string stringArrayQ[2] = {"Q1  -> Q0 ", "Q1  -> Q1 "};


    // Import/Generate an hexahedral and  a Tetrahedral mesh.

    std::shared_ptr<RegionMesh<LinearTetra> > fullMeshTetraPtr ( new RegionMesh<LinearTetra> ( Comm ) );
    std::shared_ptr<RegionMesh<LinearHexa> > fullMeshHexaPtr ( new RegionMesh<LinearHexa> ( Comm ) );
    UInt nEl (10);
    GetPot dataFile ("./data");
    MeshData meshData (dataFile, "interpolate/space_discretization");
    readMesh (*fullMeshHexaPtr, meshData);
    regularMesh3D ( *fullMeshTetraPtr, 1, nEl, nEl, nEl);


    // Partition the meshes using ParMetis
    std::shared_ptr<RegionMesh<LinearHexa> > localMeshHexaPtr;
    {
        // Create the partitioner
        MeshPartitioner< RegionMesh<LinearHexa> >  meshPartHexa;

        // Partition the mesh using ParMetis
        meshPartHexa.doPartition ( fullMeshHexaPtr, Comm );

        // Get the local mesh
        localMeshHexaPtr = meshPartHexa.meshPartition();
    }
    std::shared_ptr<RegionMesh<LinearTetra> > localMeshTetraPtr;
    {
        // Create the partitioner
        MeshPartitioner< RegionMesh<LinearTetra> >  meshPartTetra;

        // Partition the mesh using ParMetis
        meshPartTetra.doPartition ( fullMeshTetraPtr, Comm );

        // Get the local mesh
        localMeshTetraPtr = meshPartTetra.meshPartition();
    }

    //Building finite element spaces

    //Finite element space of the first scalar field - P0
    FESpaceTetraPtr_Type feSpaceP0 ( new FESpaceTetra_Type ( localMeshTetraPtr, feTetraP0, quadRuleTetra4pt,
                                                             quadRuleTria4pt, 3, Comm ) );

    //Finite element space of the second scalar field - P1
    FESpaceTetraPtr_Type feSpaceP1 ( new FESpaceTetra_Type ( localMeshTetraPtr, feTetraP1, quadRuleTetra4pt,
                                                             quadRuleTria4pt, 3, Comm ) );

    // Finite element space of the second scalar field - P1bubble
    FESpaceTetraPtr_Type feSpaceP1Bubble ( new FESpaceTetra_Type ( localMeshTetraPtr, feTetraP1bubble, quadRuleTetra15pt,
                                                                   quadRuleTria4pt, 3, Comm ) );

    // Finite element space of the second scalar field - P2
    FESpaceTetraPtr_Type feSpaceP2 ( new FESpaceTetra_Type ( localMeshTetraPtr, feTetraP2, quadRuleTetra4pt,
                                                             quadRuleTria4pt, 3, Comm ) );

    // Finite element space of the first scalar field - Q0
    FESpaceHexaPtr_Type feSpaceQ0 ( new FESpaceHexa_Type ( localMeshHexaPtr, feHexaQ0, quadRuleHexa8pt,
                                                           quadRuleQuad4pt, 2, Comm ) );

    // Finite element space of the second scalar field - Q1
    FESpaceHexaPtr_Type feSpaceQ1 ( new FESpaceHexa_Type ( localMeshHexaPtr, feHexaQ1, quadRuleHexa8pt,
                                                           quadRuleQuad4pt, 2, Comm ) );

    //vectors containing the original and final  FE spaces
    //(FE vectors will be interpolated from original FE spaces into final FE Spaces)
    std::vector<FESpaceHexaPtr_Type > originalFeSpaceHexaVec (1),  finalFeSpaceHexaVec (2);
    std::vector<FESpaceTetraPtr_Type> originalFeSpaceTetraVec (3), finalFeSpaceTetraVec (4);

    originalFeSpaceHexaVec[0] = feSpaceQ1;
    finalFeSpaceHexaVec[0] = feSpaceQ0;
    finalFeSpaceHexaVec[1] = feSpaceQ1;

    originalFeSpaceTetraVec[0] = feSpaceP1;
    originalFeSpaceTetraVec[1] = feSpaceP1Bubble;
    originalFeSpaceTetraVec[2] = feSpaceP2;

    finalFeSpaceTetraVec[0] = feSpaceP0;
    finalFeSpaceTetraVec[1] = feSpaceP1;
    finalFeSpaceTetraVec[2] = feSpaceP1Bubble;
    finalFeSpaceTetraVec[3] = feSpaceP2;

    Real time = 0.5;

    if (verbose)
    {
        std::cout << "\nA bilinear function is interpolated into Q1 vector. \nThen this FE vector is interpolated into Q0 and Q1 vectors. \nThese are the errors with respect to the analytical solution.\n";
    }
    check = check_interpolate (originalFeSpaceHexaVec, finalFeSpaceHexaVec, Unique, bilinearFunction, errArrayBilinear, stringArrayQ, 1e-10, time, verbose);


    if (verbose)
        std::cout << "\nA linear function is interpolated into P1, P1b, P2 vectors. "
                  "\nThen these FE vectors are interpolated into the following finite elements \nand error with respect to the analytic solution are reported.\n";

    check &= check_interpolate (originalFeSpaceTetraVec, finalFeSpaceTetraVec, Repeated, linearFunction, errArrayLinear, stringArrayP, 1e-10, time, verbose);


    if (verbose)
        std::cout << "\nA quadratic function is interpolated into P1, P1b, P2 vectors. "
                  "\nThen these FE vectors are interpolated into the following finite elements \nand error with respect to the analytic solution are reported.\n";

    check &= check_interpolate (originalFeSpaceTetraVec, finalFeSpaceTetraVec, Unique, quadraticFunction, errArrayQuadratic, stringArrayP, 1e-10, time, verbose);


    if (verbose)
        std::cout << "\nA linear bubble function is interpolated into P1, P1b, P2 vectors. "
                  "\nThen these FE vectors are interpolated into the following finite elements \nand error with respect to the analytic solution are reported.\n";

    check &= check_interpolate (originalFeSpaceTetraVec, finalFeSpaceTetraVec, Repeated, linearBubbleFunction,  errArrayBubble, stringArrayP, 1e-10, time, verbose);


#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    if (verbose)
    {
        if (check)
        {
            std::cout << "\nTEST INTERPOLATE WAS SUCCESSFUL.\n\n";
        }
        else
        {
            std::cout << "\nTEST INTERPOLATE FAILED.\n\n";
        }
    }

    if (check)
    {
        return EXIT_SUCCESS;
    }
    else
    {
        return EXIT_FAILURE;
    }
}//end main


Real linearFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic)
{
    return  2 * x + -z + 4 + y * (ic + 1) + t;
}

//linear function containing a bubble. The bubble is defined on a particular mesh tetrahedra.
Real linearBubbleFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic)
{
    Real xx = std::max (x, 0.);
    Real yy = std::max (y - x, 0.);
    Real zz = std::max (z - y, 0.);
    return  2 * x + -z + 4 + y * ic + t + xx * yy * zz * std::max (0.1 - xx - yy - zz, 0.);
}

Real quadraticFunction (const Real& t, const Real& x, const Real& y, const Real& z, const ID& ic)
{
    return  2 * x * x + t * y * y - z * z + x * y - x * z + 2 * y * z - x + y * (ic + 1) - z - 8;
}

Real bilinearFunction (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& ic)
{
    return  x * y - x * z + 2 * y * z - x + y * (ic + 1) - z - 8;
}

