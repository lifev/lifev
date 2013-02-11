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
  @brief Offline mesh partitioning for FSI

  @date 24-08-2010
  @author Radu Popescu <radu.popescu@epfl.ch>

  @maintainer Radu Popescu <radu.popescu@epfl.ch>
*/

#ifndef MESH_PARTITIONER_OFFLINE_FSI_H
#define MESH_PARTITIONER_OFFLINE_FSI_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/DOFInterface3Dto3D.hpp>
#include <lifev/core/mesh/Marker.hpp>

namespace LifeV
{

//! MeshPartitionerOfflineFSI - Offline mesh partitioning for FSI
/*!
  @author Radu Popescu

  This class handles the offline partitioning for FSI simulations.
  It cuts the two unpartitioned meshes for the fluid and the solid
  into a set number of partitions for each. It also generates the
  map of the degrees of freedom at the interface between the two
  materials.

  Usage: create object
  call setup(...) method
  call execute() method
*/
template<typename MeshType>
class MeshPartitionerOfflineFSI
{
public:

    //! @name Public Types
    //@{
    typedef MeshType mesh_Type;
    typedef MeshType uncutMesh_Type;
    typedef boost::shared_ptr<uncutMesh_Type> uncutMeshPtr_Type;
    typedef std::vector<uncutMeshPtr_Type> uncutMeshVector_Type;
    typedef boost::shared_ptr<uncutMeshVector_Type> uncutMeshVectorPtr_Type;

    typedef std::vector<std::vector<Int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graphPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;

    typedef FESpace<uncutMesh_Type, MapEpetra> feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;
    typedef std::vector<feSpacePtr_Type> feSpaceVector_Type;
    typedef boost::shared_ptr<feSpaceVector_Type> feSpaceVectorPtr_Type;

    typedef DOFInterface3Dto3D interface_Type;
    typedef boost::shared_ptr<interface_Type> interfacePtr_Type;
    typedef std::vector<interfacePtr_Type> interfaceVector_Type;
    // The vector contains pointers to each fluid partition's interface with
    // the solid. The vector must be wrapped in a pointer so it can be stored
    // inside HDF5Filter3DMesh when doing output.
    typedef boost::shared_ptr<interfaceVector_Type> interfaceVectorPtr_Type;

    typedef MeshPartitioner<uncutMesh_Type> meshCutter_Type;
    typedef boost::scoped_ptr<meshCutter_Type> meshCutterPtr_Type;

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MeshPartitionerOfflineFSI() {}

    //! Destructor
    virtual ~MeshPartitionerOfflineFSI() {}
    //@}

    //! @name Methods
    //@{

    //! Setup the data members of the class after construction
    /*!
      This methods is called to configure the MeshPartitionerOfflineFSI object
      after it is constructed.
      @param uncutFluidMesh const boost::shared_ptr to the unpartitioned
      fluid mesh
      @param uncutSolidMesh const boost::shared_ptr to the unpartitioned
      solid mesh
      @param fluidPartitionNumber Int
      @param solidPartitionNumber Int
      @param velocityOrder std::string
      @param displacementOrder std::string
      @param fluidInterfaceFlag LifeV::MarkerIDStandardPolicy::markerID_Type (Int)
      @param solidInterfaceFlag LifeV::MarkerIDStandardPolicy::markerID_Type (Int)
      @param interfaceTolerance Real
      @param fluidInterfaceVertexFlag Int
      @param comm boost::shared_ptr to a Epetra_Comm object
    */
    void setup (const uncutMeshPtr_Type& uncutFluidMesh,
                const uncutMeshPtr_Type& uncutSolidMesh,
                const Int& fluidPartitionNumber,
                const Int& solidPartitionNumber,
                const std::string& velocityOrder,
                const std::string& displacementOrder,
                const markerID_Type& fluidInterfaceFlag,
                const markerID_Type& solidInterfaceFlag,
                const Real& interfaceTolerance,
                const Int& fluidInterfaceVertexFlag,
                const Int& solidInterfaceVertexFlag,
                const commPtr_Type& comm);

    //! Execute the partitioning and create the interface map
    /*!
      This is a wrapper method that calls the following private methods:
      runTheCutters(), createSpaces(), mapTheInterface()
    */
    void execute();

    //! Display general information about the content of the class
    /*!
      Displays the current value of the data members.
      @param output std::ostream - stream used for output
    */
    void showMe (std::ostream& output = std::cout) const;

    //@}

    //! @name Get Methods
    //@{
    const markerID_Type& fluidInterfaceFlag() const
    {
        return M_fluidInterfaceFlag;
    }
    const markerID_Type& solidInterfaceFlag() const
    {
        return M_solidInterfaceFlag;
    }
    const graphPtr_Type& fluidGraph() const
    {
        return M_fluidMeshCutter->elementDomains();
    }
    const graphPtr_Type& solidGraph() const
    {
        return M_solidMeshCutter->elementDomains();
    }
    const uncutMeshVectorPtr_Type& fluidPartitions() const
    {
        return M_fluidMeshCutter->meshPartitions();
    }

    const uncutMeshVectorPtr_Type& solidPartitions() const
    {
        return M_solidMeshCutter->meshPartitions();
    }
    /*
      const interfaceVectorPtr_Type& dofFluidToStructure()
      {
      return M_dofFluidToStructure;
      }
    */
    const interfaceVectorPtr_Type& dofStructureToHarmonicExtension() const
    {
        return M_dofStructureToHarmonicExtension;
    }
    /*
      const interfaceVectorPtr_Type& dofStructureToSolid()
      {
      return M_dofStructureToSolid;
      }

      const interfaceVectorPtr_Type& dofStructureToFluid()
      {
      return M_dofStructureToFluid;
      }

      const interfaceVectorPtr_Type& dofHarmonicExtensionToFluid()
      {
      return M_dofHarmonicExtensionToFluid;
      }
    */
    //@}

private:
    // !!! ATTENTION !!!
    // Copy constructor and assignment operator are disabled
    // because there is no need for such operations in this case.
    // As a safety measure, there are no definitions for these
    // functions, so any attempts to use them result in link errors.
    MeshPartitionerOfflineFSI (const MeshPartitionerOfflineFSI&);
    MeshPartitionerOfflineFSI& operator= (const MeshPartitionerOfflineFSI&);

    //! @name Private Methods
    //@{

    //! Create the MeshPartitioner objects for fluid and solid
    /*!
      Allocate the MeshPartitioner objects for the fluid and solid and
      execute the partitioning of the two materials
    */
    void runTheCutters();

    //! Create the finite element spaces
    /*!
      This method creates the finite element spaces necessary for
      the offline partitioning: one FE space for the uncut solid
      mesh and one FE space for each of the fluid mesh partitions.
    */
    void createSpaces();

    //! Create the interface map between fluid and solid
    /*!
      This method creates the interface map between the two
      regions. One interface map is created between each of the
      fluid partitions and the solid uncut mesh
    */
    void mapTheInterface();

    //@}

    commPtr_Type M_comm;

    Int M_fluidPartitionNumber;
    Int M_solidPartitionNumber;

    Real M_interfaceTolerance;

    markerID_Type M_fluidInterfaceFlag;
    markerID_Type M_solidInterfaceFlag;

    boost::scoped_ptr<const Int> M_fluidInterfaceVertexFlag;
    boost::scoped_ptr<const Int> M_solidInterfaceVertexFlag;

    std::string M_velocityOrder;
    std::string M_displacementOrder;

    uncutMeshPtr_Type M_uncutFluidMesh;
    uncutMeshPtr_Type M_uncutSolidMesh;

    meshCutterPtr_Type M_fluidMeshCutter;
    meshCutterPtr_Type M_solidMeshCutter;

    feSpaceVectorPtr_Type M_velocityFESpaces;
    feSpacePtr_Type M_displacementFESpace;

    interfaceVectorPtr_Type M_dofFluidToStructure;
    interfaceVectorPtr_Type M_dofStructureToHarmonicExtension;
    interfaceVectorPtr_Type M_dofStructureToSolid;
    interfaceVectorPtr_Type M_dofStructureToFluid;
    interfaceVectorPtr_Type M_dofHarmonicExtensionToFluid;
};

///////////////////////////
//    IMPLEMENTATION     //
///////////////////////////

// ===============================
// Pubic methods
// ===============================

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::setup (const uncutMeshPtr_Type& uncutFluidMesh,
                                                 const uncutMeshPtr_Type& uncutSolidMesh,
                                                 const Int& fluidPartitionNumber,
                                                 const Int& solidPartitionNumber,
                                                 const std::string& velocityOrder,
                                                 const std::string& displacementOrder,
                                                 const markerID_Type& fluidInterfaceFlag,
                                                 const markerID_Type& solidInterfaceFlag,
                                                 const Real& interfaceTolerance,
                                                 const Int& fluidInterfaceVertexFlag,
                                                 const Int& solidInterfaceVertexFlag,
                                                 const commPtr_Type& comm)
{
    M_comm = comm;

    // Make sure we are running in serial
    if (M_comm->NumProc() == 1)
    {
        M_fluidPartitionNumber = fluidPartitionNumber;
        M_solidPartitionNumber = solidPartitionNumber;

        M_velocityOrder = velocityOrder;
        M_displacementOrder = displacementOrder;

        M_fluidInterfaceFlag = fluidInterfaceFlag;
        M_solidInterfaceFlag = solidInterfaceFlag;
        M_interfaceTolerance = interfaceTolerance;

        if (fluidInterfaceVertexFlag > 0)
        {
            M_fluidInterfaceVertexFlag.reset (new const Int (fluidInterfaceVertexFlag) );
        }
        if (solidInterfaceVertexFlag > 0)
        {
            M_solidInterfaceVertexFlag.reset (new const Int (solidInterfaceVertexFlag) );
        }

        M_uncutFluidMesh = uncutFluidMesh;
        M_uncutSolidMesh = uncutSolidMesh;

        // Allocate and configure the mesh cutters
        M_fluidMeshCutter.reset (new meshCutter_Type);
        M_fluidMeshCutter->setup (M_fluidPartitionNumber, M_comm);
        M_fluidMeshCutter->attachUnpartitionedMesh (M_uncutFluidMesh);

        M_solidMeshCutter.reset (new meshCutter_Type);
        M_solidMeshCutter->setup (M_solidPartitionNumber, M_comm);
        M_solidMeshCutter->attachUnpartitionedMesh (M_uncutSolidMesh);
    }
    else
    {
        ERROR_MSG ("Offline FSI partitioning is designed to run"
                   " with a single process.");
    }
}

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::execute()
{
    // Cut the meshes ...
    runTheCutters();

    // ... create the appropriate FE spaces ...
    createSpaces();

    // ... then create the F-S interface
    mapTheInterface();
}

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::showMe (std::ostream& output) const
{
    output << std::endl;
    output << "=======================================" << std::endl;
    output << "Internal state of MeshPartitionerOfflineFSI" << std::endl;
    output << std::endl;

    output << "Number of fluid partitions: " << M_fluidPartitionNumber
           << std::endl;
    output << "Number of solid partitions: " << M_solidPartitionNumber
           << std::endl;
    output << "Velocity order: " << M_velocityOrder
           << std::endl;
    output << "Displacement order: " << M_displacementOrder
           << std::endl;
    output << "Fluid interface flag: " << M_fluidInterfaceFlag
           << std::endl;
    output << "Solid interface flag: " << M_solidInterfaceFlag
           << std::endl;
    output << "Interface tolerance: " << M_interfaceTolerance
           << std::endl;
    if (M_fluidInterfaceVertexFlag.get() )
    {
        output << "Fluid interface vertex flag: " << *M_fluidInterfaceVertexFlag
               << std::endl;
    }
    if (M_solidInterfaceVertexFlag.get() )
    {
        output << "Solid interface vertex flag: " << *M_solidInterfaceVertexFlag
               << std::endl;
    }
    output << "Fluid mesh is stored at: " << M_uncutFluidMesh.get()
           << std::endl;
    output << "Solid mesh is stored at: " << M_uncutSolidMesh.get()
           << std::endl;
    output << "=======================================" << std::endl;
    output << std::endl;
}

// ===============================
// Private methods
// ===============================

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::runTheCutters()
{
    std::cout << "\nPartitioning fluid mesh...\n" << std::endl;
    M_fluidMeshCutter->doPartitionGraph();
    M_fluidMeshCutter->doPartitionMesh();
    M_fluidMeshCutter->releaseUnpartitionedMesh();
    std::cout << "\nPartitioning solid mesh...\n" << std::endl;
    M_solidMeshCutter->doPartitionGraph();
    M_solidMeshCutter->doPartitionMesh();
    M_solidMeshCutter->releaseUnpartitionedMesh();
    std::cout << std::endl;
}

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::createSpaces()
{
    // Set the appropriate reference elements and quad rules
    const ReferenceFE*    refFE_vel (0);
    const QuadratureRule* qR_vel (0);
    const QuadratureRule* bdQr_vel (0);

    const ReferenceFE*    refFE_disp (0);
    const QuadratureRule* qR_disp (0);
    const QuadratureRule* bdQr_disp (0);

    if ( M_velocityOrder.compare ("P2") == 0 )
    {
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt;
        bdQr_vel  = &quadRuleTria3pt;
    }
    else
    {
        if ( M_velocityOrder.compare ("P1") == 0 )
        {
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;
            bdQr_vel  = &quadRuleTria3pt;
        }
        else
        {
            if ( M_velocityOrder.compare ("P1Bubble") == 0 )
            {
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;
                bdQr_vel  = &quadRuleTria3pt;
            }
            else
            {
                ERROR_MSG (M_velocityOrder + " velocity FE not implemented yet.");
            }
        }
    }

    if ( M_displacementOrder.compare ("P2") == 0 )
    {
        refFE_disp = &feTetraP2;
        qR_disp    = &quadRuleTetra15pt;
        bdQr_disp  = &quadRuleTria3pt;
    }
    else
    {
        if ( M_displacementOrder.compare ("P1") == 0 )
        {
            refFE_disp = &feTetraP1;
            qR_disp    = &quadRuleTetra4pt;
            bdQr_disp  = &quadRuleTria3pt;
        }
        else
        {
            ERROR_MSG (M_displacementOrder + " structure FE not implemented yet.");
        }
    }

    // Create finite element spaces for velocity and displacement
    std::cout << "Creating velocity finite element space... ";
    M_velocityFESpaces.reset (new feSpaceVector_Type);
    M_velocityFESpaces->resize (M_fluidPartitionNumber);
    for (Int i = 0; i < M_fluidPartitionNumber; ++i)
    {
        (*M_velocityFESpaces) [i].reset (new feSpace_Type (M_fluidMeshCutter->getPartition (i),
                                                           *refFE_vel,
                                                           *qR_vel,
                                                           *bdQr_vel,
                                                           3,
                                                           M_comm) );
    }
    std::cout << "done." << std::endl;
    std::cout << "Creating displacement finite element space... ";
    M_displacementFESpace.reset (new feSpace_Type (M_uncutSolidMesh,
                                                   *refFE_disp,
                                                   *qR_disp,
                                                   *bdQr_disp,
                                                   3,
                                                   M_comm) );
    std::cout << "done." << std::endl;
}

template<typename MeshType>
void MeshPartitionerOfflineFSI<MeshType>::mapTheInterface()
{
    // Create the DOFInterface3Dto3D objects for each fluid partition
    std::cout << std::endl;
    std::cout << "Creating the DOF interfaces..." << std::endl;
    std::cout << std::endl;
    std::cout << "Fluid interface flag is: " << M_fluidInterfaceFlag << std::endl;
    std::cout << "Solid interface flag is: " << M_solidInterfaceFlag << std::endl;
    std::cout << std::endl;

    M_dofStructureToHarmonicExtension.reset (new interfaceVector_Type);
    M_dofStructureToHarmonicExtension->resize (M_fluidPartitionNumber);

    for (Int i = 0; i < M_fluidPartitionNumber; ++i)
    {
        interfacePtr_Type& ifPtr = (*M_dofStructureToHarmonicExtension) [i];
        feSpacePtr_Type& velSpacePtr = (*M_velocityFESpaces) [i];
        ifPtr.reset (new interface_Type);
        ifPtr->setup (velSpacePtr->refFE(),
                      velSpacePtr->dof(),
                      M_displacementFESpace->refFE(),
                      M_displacementFESpace->dof() );
        ifPtr->update (* (velSpacePtr->mesh() ), M_fluidInterfaceFlag,
                       * (M_displacementFESpace->mesh() ), M_solidInterfaceFlag,
                       M_interfaceTolerance, M_fluidInterfaceVertexFlag.get() );
    }
    std::cout << "dofStructureToHarmonicExtension done." << std::endl;

    std::cout << "All done." << std::endl << std::endl;
}

} // Namespace LifeV

#endif // MESH_PARTITIONER_OFFLINE_FSI_H
