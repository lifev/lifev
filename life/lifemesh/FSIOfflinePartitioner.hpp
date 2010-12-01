//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
  @file
  @brief Offline mesh partitioning for FSI

  @author Radu Popescu <rpopescu@mathicsepc3.epfl.ch>
  @date 24 Aug 2010

  Class thet handles the fluid and solid mesh partitioning for FSI
  simulations and also generates the interface map between the
  two regions.
*/

#ifndef FSIOFFLINEPARTITIONER_H
#define FSIOFFLINEPARTITIONER_H 1

#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>

#include <life/lifemesh/markers_base.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/dofInterface3Dto3D.hpp>

namespace LifeV
{

//! FSIOfflinePartitioner - Offline mesh partitioning for FSI
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
template<typename Mesh>
class FSIOfflinePartitioner
{
public:

    //! @name Public Types
    //@{
    typedef Mesh uncutMesh_Type;
    typedef boost::shared_ptr<uncutMesh_Type> uncutMesh_ptrType;
    typedef std::vector<uncutMesh_ptrType> uncutMesh_vector_Type;
    typedef boost::shared_ptr<uncutMesh_vector_Type> uncutMesh_vector_ptrType;

    typedef std::vector<std::vector<int> > graph_Type;
    typedef boost::shared_ptr<graph_Type> graph_ptrType;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> comm_ptrType;

    typedef FESpace<uncutMesh_Type, EpetraMap> feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpace_ptrType;
    typedef std::vector<feSpace_ptrType> feSpace_vector_Type;
    typedef boost::shared_ptr<feSpace_vector_Type> feSpace_vector_ptrType;

    typedef DofInterface3Dto3D interface_Type;
    typedef boost::shared_ptr<interface_Type> interface_ptrType;
    typedef std::vector<interface_ptrType> interface_vector_Type;
    // The vector contains pointers to each fluid partition's interface with
    // the solid. The vector must be wrapped in a pointer so it can be stored
    // inside HDF5Filter3DMesh when doing output.
    typedef boost::shared_ptr<interface_vector_Type> interface_vector_ptrType;

    typedef partitionMesh<uncutMesh_Type> meshCutter_Type;
    typedef boost::scoped_ptr<meshCutter_Type> meshCutter_ptrType;

    typedef MarkerTraits_Base::EntityFlag entityFlag_Type;
    //@}

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    FSIOfflinePartitioner() {}

    //! Destructor
    virtual ~FSIOfflinePartitioner() {}

    //@}

    //! @name Methods
    //@{

    //! Setup the data members of the class after construction
    /*!
      This methods is called to configure the FSIOfflinePartitioner object
      after it is constructed.
      @param uncutFluidMesh const boost::shared_ptr to the unpartitioned
      fluid mesh
      @param uncutSolidMesh const boost::shared_ptr to the unpartitioned
      solid mesh
      @param fluidPartitionNumber int
      @param solidPartitionNumber int
      @param velocityOrder std::string
      @param displacementOrder std::string
      @param fluidInterfaceFlag LifeV::MarkerTraits_Base::EntityFlag (int)
      @param solidInterfaceFlag LifeV::MarkerTraits_Base::EntityFlag (int)
      @param interfaceTolerance double
      @param fluidInterfaceVertexFlag int
      @param comm boost::shared_ptr to a Epetra_Comm object
    */
    void setup(const uncutMesh_ptrType& uncutFluidMesh,
               const uncutMesh_ptrType& uncutSolidMesh,
               const int& fluidPartitionNumber,
               const int& solidPartitionNumber,
               const std::string& velocityOrder,
               const std::string& displacementOrder,
               const entityFlag_Type& fluidInterfaceFlag,
               const entityFlag_Type& solidInterfaceFlag,
               const Real& interfaceTolerance,
               const int& fluidInterfaceVertexFlag,
               const int& solidInterfaceVertexFlag,
               const comm_ptrType& comm);

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
    void showMe(std::ostream& output = std::cout) const;

    //@}

    //! @name Operators
    //@{

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{
    const int fluidInterfaceFlag() const
    {
        return M_fluidInterfaceFlag;
    }
    const int solidInterfaceFlag() const
    {
        return M_solidInterfaceFlag;
    }
    const graph_ptrType& fluidGraph() const
    {
        return M_fluidMeshCutter->graph();
    }
    const graph_ptrType& solidGraph() const
    {
        return M_solidMeshCutter->graph();
    }
    const uncutMesh_vector_ptrType& fluidPartitions() const
    {
        return M_fluidMeshCutter->meshAllPartitions();
    }

    const uncutMesh_vector_ptrType& solidPartitions() const
    {
        return M_solidMeshCutter->meshAllPartitions();
    }
    /*
        const interface_vector_ptrType& dofFluidToStructure()
        {
            return M_dofFluidToStructure;
        }
    */
    const interface_vector_ptrType& dofStructureToHarmonicExtension() const
    {
        return M_dofStructureToHarmonicExtension;
    }
    /*
        const interface_vector_ptrType& dofStructureToSolid()
        {
            return M_dofStructureToSolid;
        }

        const interface_vector_ptrType& dofStructureToFluid()
        {
            return M_dofStructureToFluid;
        }

        const interface_vector_ptrType& dofHarmonicExtensionToFluid()
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
    FSIOfflinePartitioner(const FSIOfflinePartitioner&);
    FSIOfflinePartitioner& operator=(const FSIOfflinePartitioner&);

    //! @name Private Methods
    //@{

    //! Create the partitionMesh objects for fluid and solid
    /*!
      Allocate the partitionMesh objects for the fluid and solid and
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

    comm_ptrType M_comm;

    int M_fluidPartitionNumber;
    int M_solidPartitionNumber;

    Real M_interfaceTolerance;

    entityFlag_Type M_fluidInterfaceFlag;
    entityFlag_Type M_solidInterfaceFlag;

    boost::scoped_ptr<const int> M_fluidInterfaceVertexFlag;
    boost::scoped_ptr<const int> M_solidInterfaceVertexFlag;

    std::string M_velocityOrder;
    std::string M_displacementOrder;

    uncutMesh_ptrType M_uncutFluidMesh;
    uncutMesh_ptrType M_uncutSolidMesh;

    meshCutter_ptrType M_fluidMeshCutter;
    meshCutter_ptrType M_solidMeshCutter;

    feSpace_vector_ptrType M_velocityFESpaces;
    feSpace_ptrType M_displacementFESpace;

    interface_vector_ptrType M_dofFluidToStructure;
    interface_vector_ptrType M_dofStructureToHarmonicExtension;
    interface_vector_ptrType M_dofStructureToSolid;
    interface_vector_ptrType M_dofStructureToFluid;
    interface_vector_ptrType M_dofHarmonicExtensionToFluid;
};

///////////////////////////
//    IMPLEMENTATION     //
///////////////////////////

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::setup(const uncutMesh_ptrType& uncutFluidMesh,
                                        const uncutMesh_ptrType& uncutSolidMesh,
                                        const int& fluidPartitionNumber,
                                        const int& solidPartitionNumber,
                                        const std::string& velocityOrder,
                                        const std::string& displacementOrder,
                                        const entityFlag_Type& fluidInterfaceFlag,
                                        const entityFlag_Type& solidInterfaceFlag,
                                        const Real& interfaceTolerance,
                                        const int& fluidInterfaceVertexFlag,
                                        const int& solidInterfaceVertexFlag,
                                        const comm_ptrType& comm)
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
            M_fluidInterfaceVertexFlag.reset(new const int(fluidInterfaceVertexFlag));
        }
        if (solidInterfaceVertexFlag > 0)
        {
            M_solidInterfaceVertexFlag.reset(new const int(solidInterfaceVertexFlag));
        }

        M_uncutFluidMesh = uncutFluidMesh;
        M_uncutSolidMesh = uncutSolidMesh;

        // Allocate and configure the mesh cutters
        M_fluidMeshCutter.reset(new meshCutter_Type);
        M_fluidMeshCutter->setup(M_fluidPartitionNumber, M_comm);
        M_fluidMeshCutter->attachUnpartitionedMesh(M_uncutFluidMesh);

        M_solidMeshCutter.reset(new meshCutter_Type);
        M_solidMeshCutter->setup(M_solidPartitionNumber, M_comm);
        M_solidMeshCutter->attachUnpartitionedMesh(M_uncutSolidMesh);
    }
    else
    {
        ERROR_MSG("Offline FSI partitioning is designed to run"
                  " with a single process.");
    }
}

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::execute()
{
    // Cut the meshes ...
    runTheCutters();

    // ... create the appropriate FE spaces ...
    createSpaces();

    // ... then create the F-S interface
    mapTheInterface();
}

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::showMe(std::ostream& output) const
{
    output << std::endl;
    output << "=======================================" << std::endl;
    output << "Internal state of FSIOfflinePartitioner" << std::endl;
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
    if (M_fluidInterfaceVertexFlag.get())
    {
        output << "Fluid interface vertex flag: " << *M_fluidInterfaceVertexFlag
        << std::endl;
    }
    if (M_solidInterfaceVertexFlag.get())
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

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::runTheCutters()
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

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::createSpaces()
{
    // Set the appropriate reference elements and quad rules
    const RefFE*    refFE_vel(0);
    const QuadRule* qR_vel(0);
    const QuadRule* bdQr_vel(0);

    const RefFE*    refFE_disp(0);
    const QuadRule* qR_disp(0);
    const QuadRule* bdQr_disp(0);

    if ( M_velocityOrder.compare("P2") == 0 )
    {
        refFE_vel = &feTetraP2;
        qR_vel    = &quadRuleTetra15pt;
        bdQr_vel  = &quadRuleTria3pt;
    }
    else
    {
        if ( M_velocityOrder.compare("P1") == 0 )
        {
            refFE_vel = &feTetraP1;
            qR_vel    = &quadRuleTetra4pt;
            bdQr_vel  = &quadRuleTria3pt;
        }
        else
        {
            if ( M_velocityOrder.compare("P1Bubble") == 0 )
            {
                refFE_vel = &feTetraP1bubble;
                qR_vel    = &quadRuleTetra64pt;
                bdQr_vel  = &quadRuleTria3pt;
            }
            else
            {
                ERROR_MSG(M_velocityOrder + " velocity FE not implemented yet.");
            }
        }
    }

    if ( M_displacementOrder.compare("P2") == 0 )
    {
        refFE_disp = &feTetraP2;
        qR_disp    = &quadRuleTetra15pt;
        bdQr_disp  = &quadRuleTria3pt;
    }
    else
    {
        if ( M_displacementOrder.compare("P1") == 0 )
        {
            refFE_disp = &feTetraP1;
            qR_disp    = &quadRuleTetra4pt;
            bdQr_disp  = &quadRuleTria3pt;
        }
        else
        {
            ERROR_MSG(M_displacementOrder + " structure FE not implemented yet.");
        }
    }

    // Create finite element spaces for velocity and displacement
    std::cout << "Creating velocity finite element space... ";
    M_velocityFESpaces.reset(new feSpace_vector_Type);
    M_velocityFESpaces->resize(M_fluidPartitionNumber);
    for (int i = 0; i < M_fluidPartitionNumber; ++i)
    {
        (*M_velocityFESpaces)[i].reset(new feSpace_Type(M_fluidMeshCutter->mesh(i),
                                                        *refFE_vel,
                                                        *qR_vel,
                                                        *bdQr_vel,
                                                        3,
                                                        M_comm));
    }
    std::cout << "done." << std::endl;
    std::cout << "Creating displacement finite element space... ";
    M_displacementFESpace.reset(new feSpace_Type(M_uncutSolidMesh,
                                                 *refFE_disp,
                                                 *qR_disp,
                                                 *bdQr_disp,
                                                 3,
                                                 M_comm));
    std::cout << "done." << std::endl;
}

template<typename Mesh>
void FSIOfflinePartitioner<Mesh>::mapTheInterface()
{
    // Create the dofInterface3Dto3D objects for each fluid partition
    std::cout << std::endl;
    std::cout << "Creating the DOF interfaces..." << std::endl;
    std::cout << std::endl;
    std::cout << "Fluid interface flag is: " << M_fluidInterfaceFlag << std::endl;
    std::cout << "Solid interface flag is: " << M_solidInterfaceFlag << std::endl;
    std::cout << std::endl;

    M_dofStructureToHarmonicExtension.reset(new interface_vector_Type);
    M_dofStructureToHarmonicExtension->resize(M_fluidPartitionNumber);

    for (int i = 0; i < M_fluidPartitionNumber; ++i)
    {
        interface_ptrType& ifPtr = (*M_dofStructureToHarmonicExtension)[i];
        feSpace_ptrType& velSpacePtr = (*M_velocityFESpaces)[i];
        ifPtr.reset(new interface_Type);
        ifPtr->setup(velSpacePtr->refFE(),
                     velSpacePtr->dof(),
                     M_displacementFESpace->refFE(),
                     M_displacementFESpace->dof());
        ifPtr->update(*(velSpacePtr->mesh()), M_fluidInterfaceFlag,
                      *(M_displacementFESpace->mesh()), M_solidInterfaceFlag,
                      M_interfaceTolerance, M_fluidInterfaceVertexFlag.get());
    }
    std::cout << "dofStructureToHarmonicExtension done." << std::endl;

    std::cout << "All done." << std::endl << std::endl;
}

} // Namespace LifeV

#endif // FSIOFFLINEPARTITIONER_H
