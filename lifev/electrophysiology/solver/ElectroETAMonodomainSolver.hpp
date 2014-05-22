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
 @brief Class for solving the Monodomain model in electrophysiology.

 @date 02 - 2013
 @author Simone Rossi <simone.rossi@epfl.ch>

 @last update 04 - 2014

 This class provides interfaces to solve the monodomain equation
 ( reaction diffusion equation ) using the ETA framework.
 The solution can be performed using three different methods:
 -operator splitting method (at this point available only with forward Euler
 for the reaction step and backward Euler for the diffusion step.
 Second order splitting is available but still experimental. );
 -Ionic Currents Interpolation (at this point only forward Euler);
 -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _ELECTROETAMONODOMAINSOLVER_H_
#define _ELECTROETAMONODOMAINSOLVER_H_


#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <Epetra_LocalMap.h>

#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>
#include <lifev/electrophysiology/stimulus/ElectroStimulus.hpp>

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>


namespace LifeV
{

//! monodomainSolver - Class featuring the solver for monodomain equations

template<typename Mesh, typename IonicModel>
class ElectroETAMonodomainSolver
{

    //!Monodomain Solver
    /*!
     The monodomain equation reads
     \f \Chi

     */

public:

    //! @name Type definitions
    //@{

    //! Mesh
    typedef Mesh                                                        mesh_Type;

    typedef boost::shared_ptr<mesh_Type>                                meshPtr_Type;

    //! Distributed vector // For parallel usage
    typedef VectorEpetra                                                vector_Type;

    typedef boost::shared_ptr<VectorEpetra>                             vectorPtr_Type;

    typedef vector_Type                                                 solution_Type;

    typedef vectorPtr_Type                                              solutionPtr_Type;

    typedef std::vector<vectorPtr_Type>                                 vectorOfPtr_Type;

    //! Distributed Matrix // For parallel usage
    typedef MatrixEpetra<Real>                                          matrix_Type;

    typedef boost::shared_ptr<matrix_Type>                              matrixPtr_Type;

    //! Communicator to exchange informations among processes
    typedef Epetra_Comm                                                 comm_Type;

    typedef boost::shared_ptr<comm_Type>                                commPtr_Type;

    //! Expression template  scalar finite element space
    //! To be used in the expression assembly namespace
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1>                       ETFESpace_Type;

    typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> >   ETFESpacePtr_Type;

    //! Expression template vectorial finite element space
    //! To be used in the expression assembly namespace
    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3>                       ETFESpaceVectorial_Type;

    typedef boost::shared_ptr<ETFESpaceVectorial_Type>                  ETFESpaceVectorialPtr_Type;

    //! Finite element space
    typedef FESpace<mesh_Type, MapEpetra>                               feSpace_Type;

    typedef boost::shared_ptr<feSpace_Type>                             feSpacePtr_Type;

    //! Linear Solver
    typedef LinearSolver                                                linearSolver_Type;

    typedef boost::shared_ptr<LinearSolver>                             linearSolverPtr_Type;

    //! Exporter to save the solution
    typedef Exporter<mesh_Type>                                         IOFile_Type;

    typedef boost::shared_ptr<IOFile_Type>                              IOFilePtr_Type;

    //! Exporter data
    typedef ExporterData<mesh_Type>                                     IOData_Type;

    typedef ExporterEnsight<mesh_Type>                                  ensightIOFile_Type;

#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type >                                   hdf5IOFile_Type;
#endif

    //! Preconditioner
    typedef LifeV::Preconditioner                                       basePrec_Type;

    typedef boost::shared_ptr<basePrec_Type>                            basePrecPtr_Type;

    //! MultiLevel Preconditioner
    typedef LifeV::PreconditionerML                                     prec_Type;

    typedef boost::shared_ptr<prec_Type>                                precPtr_Type;

    //! Ionic model
    typedef IonicModel                                                  ionicModel_Type;

    //! Base class of the ionic model
    typedef ElectroIonicModel                                           superIonicModel;

    typedef boost::shared_ptr<superIonicModel>                          ionicModelPtr_Type;

    //! xml list to read parameters
    typedef Teuchos::ParameterList                                      list_Type;

    //! boost function
    typedef boost::function < Real (const Real& t,
                                    const Real& x,
                                    const Real& y,
                                    const Real& z,
                                    const ID&   i) >                     function_Type;

    //! 3x3 matrix
    typedef MatrixSmall<3, 3>                                           matrixSmall_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    ElectroETAMonodomainSolver();

    //! Constructor
    /*!
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel>  chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    ElectroETAMonodomainSolver (GetPot& dataFile, ionicModelPtr_Type model,
                                meshPtr_Type meshPtr);

    //! Constructor
    /*!
     * @param meshName file name of the mesh
     * @param meshPath path to the mesh
     * @param datafile  GetPot file to setup the preconditioner
     * @param model shared pointer to the chosen ionic model
     */
    ElectroETAMonodomainSolver (std::string meshName, std::string meshPath,
                                GetPot& dataFile, ionicModelPtr_Type model);

    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */
    ElectroETAMonodomainSolver (std::string meshName, std::string meshPath,
                                GetPot& dataFile, ionicModelPtr_Type model, commPtr_Type comm);

    //! Copy Constructor
    /*!
     * @param ElectroETAmonodomainSolver object
     */
    ElectroETAMonodomainSolver (const ElectroETAMonodomainSolver& solver);

    //!Operator=()
    /*!
     * @param ElectroETAmonodomainSolver object
     */
    ElectroETAMonodomainSolver<Mesh, IonicModel>& operator= (
        const ElectroETAMonodomainSolver& solver);

    //! Destructor
    virtual ~ElectroETAMonodomainSolver()
    {
    }

    //@}

    //! @name Get Methods
    //@{

    //! get the surface to volume ratio
    /*!
     * Not used in the code ( implicit definition inside the diffusion tensor)
     */
    inline const Real& surfaceVolumeRatio() const
    {
        return M_surfaceVolumeRatio;
    }

    //! get the initial time (by default 0)
    inline const Real& initialTime() const
    {
        return M_initialTime;
    }

    //! get the final time
    inline const Real& timeStep() const
    {
        return M_timeStep;
    }

    //! get the time step
    inline const Real& endTime() const
    {
        return M_endTime;
    }

    //! get the diagonal diffusion tensor
    inline const VectorSmall<3>& diffusionTensor() const
    {
        return M_diffusionTensor;
    }

    //! get the order of the elements
    inline const std::string elementsOrder() const
    {
        return M_elementsOrder;
    }

    //! get the pointer to the ionic model
    inline const ionicModelPtr_Type ionicModelPtr() const
    {
        return M_ionicModelPtr;
    }

    //! get the pointer to the Epetra communicator
    inline const commPtr_Type commPtr() const
    {
        return M_commPtr;
    }

    //! get the pointer to the partitioned mesh
    inline const meshPtr_Type localMeshPtr() const
    {
        return M_localMeshPtr;
    }

    //! get the pointer to the partitioned mesh
    inline meshPtr_Type localMeshPtr()
    {
        return M_localMeshPtr;
    }

    //! get the pointer to the partitioned mesh
    inline const meshPtr_Type fullMeshPtr() const
    {
        return M_fullMeshPtr;
    }

    //! get the pointer to the ETA finite element space
    inline const ETFESpacePtr_Type ETFESpacePtr() const
    {
        return M_ETFESpacePtr;
    }

    //! get the pointer to the usual finite element space
    inline const feSpacePtr_Type feSpacePtr() const
    {
        return M_feSpacePtr;
    }

    //! get the pointer to the mass matrix
    inline const matrixPtr_Type massMatrixPtr() const
    {
        return M_massMatrixPtr;
    }

    //! get the pointer to the mass matrix
    inline const matrixPtr_Type fullMassMatrixPtr() const
    {
        return M_fullMassMatrixPtr;
    }

    //! get the pointer to the stiffness matrix
    inline const matrixPtr_Type stiffnessMatrixPtr() const
    {
        return M_stiffnessMatrixPtr;
    }

    //! get the pointer to the global matrix
    /*!
     *  \f[
     *  A = \frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     */
    inline const matrixPtr_Type globalMatrixPtr() const
    {
        return M_globalMatrixPtr;
    }

    //! get the pointer to the right hand side
    inline const vectorPtr_Type rhsPtr() const
    {
        return M_rhsPtr;
    }

    //! get the pointer to the unique version of the right hand side
    inline const vectorPtr_Type rhsPtrUnique() const
    {
        return M_rhsPtrUnique;
    }

    //! get the pointer to the transmembrane potential
    inline const vectorPtr_Type potentialPtr() const
    {
        return M_potentialPtr;
    }

    //! get the pointer to the fiber vector
    inline const vectorPtr_Type fiberPtr() const
    {
        return M_fiberPtr;
    }

    //! get the pointer to the fiber vector
    inline vectorPtr_Type fiberPtr()
    {
        return M_fiberPtr;
    }

    //! get the pointer to the applied current vector
    inline const vectorPtr_Type appliedCurrentPtr()
    {
        return M_ionicModelPtr->appliedCurrentPtr();
    }

    //! get the pointer to the linear solver
    inline const linearSolverPtr_Type linearSolverPtr() const
    {
        return M_linearSolverPtr;
    }

    //! get the pointer to the vector of pointers containing the transmembrane potential (at 0) and the gating variables
    inline const vectorOfPtr_Type& globalSolution() const
    {
        return M_globalSolution;
    }

    //! get the pointer to the vector of pointers containing the rhs for transmembrane potential (at 0) and the gating variables
    inline const vectorOfPtr_Type& globalRhs() const
    {
        return M_globalRhs;
    }

    //! getter for the boolean to know if we want a lumped matrix
    inline bool lumpedMassMatrix() const
    {
        return M_lumpedMassMatrix;
    }
    //@}

    //! @name Set Methods
    //@{

    //! set the surface to volume ratio
    /*!
     @param surfaceVolumeRatio surface to volume ratio
     */
    inline void setSurfaceVolumeRatio (const Real& surfaceVolumeRatio)
    {
        this->M_surfaceVolumeRatio = surfaceVolumeRatio;
    }

    //! set the starting time
    /*!
     @param initialTime initial time
     */
    inline void setInitialTime (const Real& initialTime)
    {
        this->M_initialTime = initialTime;
    }

    //! set the ending time
    /*!
     @param timeStep ending time
     */
    inline void setTimeStep (const Real& timeStep)
    {
        this->M_timeStep = timeStep;
    }

    //! set the time step
    /*!
     @param endTime time step
     */
    inline void setEndTime (const Real& endTime)
    {
        this->M_endTime = endTime;
    }

    //! set the diagonal diffusion tensor
    /*!
     @param  diffusionTensor diagonal diffusion tensor
     */
    inline void setDiffusionTensor (const VectorSmall<3>& diffusionTensor)
    {
        this->M_diffusionTensor = diffusionTensor;
    }

    //! set the pointer to the ionic model
    /*!
     @param ionicModelPtr pointer to the ionic model
     */
    inline void setIonicModelPtr (const ionicModelPtr_Type ionicModelPtr)
    {
        this->M_ionicModelPtr = ionicModelPtr;
    }

    //! set  ionic model
    /*!
     @param ionicModel  ionic model
     */
    inline void setIonicModel (const ionicModel_Type& ionicModel)
    {
        (* (M_ionicModelPtr) ) = ionicModel;
    }


    //! set the pointer to the Epetra communicator
    /*!
     @param commPtr pointer to the Epetra communicator
     */
    inline void setCommPtr (const commPtr_Type commPtr)
    {
        this->M_commPtr = commPtr;
    }

    //! set the Epetra communicator
    /*!
     @param comm Epetra communicator
     */    inline void setComm (const comm_Type& comm)
    {
        (* (M_commPtr) ) = comm;
    }

    //! set the pointer to the partitioned mesh
    /*!
     @param localMeshPtr pointer to the partitioned mesh
     */
    inline void setLocalMeshPtr (const meshPtr_Type localMeshPtr)
    {
        this->M_localMeshPtr = localMeshPtr;
    }

    //! set the partitioned mesh
    /*!
     @param localMesh  partitioned mesh
     */
    inline void setLocalMesh (const mesh_Type& localMesh)
    {
        (* (M_localMeshPtr) ) = localMesh;
    }


    //! set the pointer to the non partitioned mesh
    /*!
     @param fullMeshPtr pointer to the partitioned mesh
     */
    inline void setFullMeshPtr (const meshPtr_Type fullMeshPtr)
    {
        this->M_fullMeshPtr = fullMeshPtr;
    }
    //! set the non partitioned mesh
    /*!
     @param fullMesh pointer to the partitioned mesh
     */
    inline void setFullMesh (const mesh_Type& fullMesh)
    {
        (* (M_fullMeshPtr) ) = fullMesh;
    }


    //! set the pointer to the ETA fe space
    /*!
     @param  ETFESpacePtr pointer to the ETA fe space
     */
    inline void setETFESpacePtr (const ETFESpacePtr_Type ETFESpacePtr)
    {
        this->M_ETFESpacePtr = ETFESpacePtr;
    }

    //! set the scalar ETA fe space
    /*!
     @param  ETFESpace  scalar ETA fe space
     */
    inline void setETFESpace (const ETFESpace_Type& ETFESpace)
    {
        (* (M_ETFESpacePtr) ) = ETFESpace;
    }

    //! set the pointer to the usual fe space
    /*!
     @param feSpacePtr pointer to the usual fe space
     */
    inline void setFeSpacePtr (const feSpacePtr_Type feSpacePtr)
    {
        this->M_feSpacePtr = feSpacePtr;
    }

    //! set the fe space
    /*!
     @param feSpace the  fe space
     */
    inline void setFeSpace (const feSpace_Type& feSpace)
    {
        (* (M_feSpacePtr) ) = feSpace;
    }


    //! set the pointer to the  mass matrix
    /*!
     @param massMatrixPtr pointer to the mass matrix
     */
    inline void setMassMatrixPtr (const matrixPtr_Type massMatrixPtr)
    {
        this->M_massMatrixPtr = massMatrixPtr;
    }

    //! set the mass matrix
    /*!
     @param massMatrix  the mass matrix
     */
    inline void setMassMatrix (const matrix_Type& massMatrix)
    {
        (* (M_massMatrixPtr) ) = massMatrix;
    }

    //! set the pointer to the full mass matrix
    /*!
     @param massMatrixPtr pointer to the mass matrix
     */
    inline void setFullMassMatrixPtr (const matrixPtr_Type fullMassMatrixPtr)
    {
        this->M_fullMassMatrixPtr = fullMassMatrixPtr;
    }

    //! set the full mass matrix
    /*!
     @param massMatrix  the mass matrix
     */
    inline void setFullMassMatrix (const matrix_Type& fullMassMatrix)
    {
        (* (M_fullMassMatrixPtr) ) = fullMassMatrix;
    }

    //! set the pointer to the stiffness matrix
    /*!
     @param stiffnessMatrixPtr pointer to the stiffness matrix
     */
    inline void setStiffnessMatrixPtr (const matrixPtr_Type stiffnessMatrixPtr)
    {
        this->M_stiffnessMatrixPtr = stiffnessMatrixPtr;
    }

    //! set the stiffness matrix
    /*!
     @param stiffnessMatrix  the stiffness matrix
     */
    inline void setStiffnessMatrix (const matrix_Type& stiffnessMatrix)
    {
        (* (M_stiffnessMatrixPtr) ) = stiffnessMatrix;
    }


    //! set the pointer to the global matrix
    /*!
     @param globalMatrix pointer to the global matrix
     */
    inline void setGlobalMatrixPtr (const matrixPtr_Type globalMatrixPtr)
    {
        this->M_globalMatrixPtr = globalMatrixPtr;
    }

    //! set the global matrix
    /*!
     @param globalMatrix the global matrix
     */
    inline void setGlobalMatrix (const matrix_Type& globalMatrix)
    {
        (* (M_globalMatrixPtr) ) = globalMatrix;
    }


    //! set the pointer to the right hand side
    /*!
     @param rhsPtr pointer to the right hand side
     */
    inline void setRhsPtr (const vectorPtr_Type rhsPtr)
    {
        this->M_rhsPtr = rhsPtr;
    }

    //! set the right hand side
    /*!
     @param rhs  the right hand side
     */
    inline void setRhs (const vector_Type& rhs)
    {
        (* (M_rhsPtr) ) = rhs;
    }


    //! set the pointer to the unique version of the right hand side
    /*!
     @param rhsPtrUnique  pointer to the  unique version of the right hand side
     */
    inline void setRhsPtrUnique (const vectorPtr_Type rhsPtrUnique)
    {
        this->M_rhsPtrUnique = rhsPtrUnique;
        this->M_globalRhs.at (0) = M_rhsPtrUnique;
    }

    //! set the unique version of the right hand side
    /*!
     @param rhsPtrUnique  the  unique version of the right hand side
     */
    inline void setRhsUnique (const vector_Type& rhsPtrUnique)
    {
        (* (M_rhsPtrUnique) ) = rhsPtrUnique;
    }

    //! set the pointer to the potential
    /*!
     @param potentialPtr   pointer to the  potential
     */
    inline void setPotentialPtr (const vectorPtr_Type potentialPtr)
    {
        this->M_potentialPtr = potentialPtr;
        this->M_globalSolution.at (0) = M_potentialPtr;
    }

    //! set the potential
    /*!
     @param potential   the  potential
     */
    inline void setPotential (const vector_Type& potential)
    {
        (* (M_potentialPtr) ) = potential;
    }

    //! set the pointer to the applied current vector
    /*!
     @param appliedCurrentPtr  pointer to the applied current vector
     */
    inline void setAppliedCurrentPtr (const vectorPtr_Type appliedCurrentPtr)
    {
        M_ionicModelPtr->setAppliedCurrentPtr (appliedCurrentPtr);
    }

    //! set the applied current vector
    /*!
     @param appliedCurrent   the applied current vector
     */
    inline void setAppliedCurrent (const vector_Type& appliedCurrent)
    {
        M_ionicModelPtr->setAppliedCurrent (appliedCurrent);
    }

    //! set the pointer to the linear solver
    /*!
     @param linearSolverPtr pointer to the linear solver
     */
    inline void setLinearSolverPtr (const linearSolverPtr_Type linearSolverPtr)
    {
        this->M_linearSolverPtr = linearSolverPtr;
    }

    //! set the linear solver
    /*!
     @param linearSolver the linear solver
     */
    inline void setLinearSolver (const linearSolver_Type& linearSolver)
    {
        (* (M_linearSolverPtr) ) = linearSolver;
    }

    //! set the vector of pointers containing the transmembrane potential (at 0) and the gating variables
    /*!
     @param globalSolution vector of pointers containing the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalSolutionPtrs (const vectorOfPtr_Type& globalSolution)
    {
        this->M_globalSolution = globalSolution;
    }

    //! set the vectors of unknowns:  containing the transmembrane potential (at 0) and the gating variables
    /*!
     @param p vector of pointers containing the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalSolution (const vectorOfPtr_Type& globalSolution)
    {
        for (int j = 0; j < M_ionicModelPtr->Size(); j++)
        {
            (* (M_globalSolution.at (j) ) ) = (* (globalSolution.at (j) ) );
        }
    }

    //! set the pointer to the \[j\]-th  gating variable
    /*!
     @param gatingVariable  pointer to the gating variable
     @param j index of the gating variable
     */
    inline void setVariablePtr (const vectorPtr_Type gatingVariable, int j)
    {
        M_globalSolution.at (j) = gatingVariable;
    }

    //! set  the \[j\]-th  gating variable
    /*!
     @param gatingVariable  the gating variable
     @param j index of the gating variable
     */
    inline void setVariablePtr (const vector_Type& gatingVariable, int j)
    {
        * (M_globalSolution.at (j) ) = gatingVariable;
    }

    //! set the vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
    /*!
     @param globalRhs vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalRhsPtrs (const vectorOfPtr_Type& globalRhs)
    {
        this->M_globalRhs = globalRhs;
    }

    //! set the vectors containing the rhs for the transmembrane potential (at 0) and the gating variables
    /*!
     @param globalRhs vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalRhs (const vectorOfPtr_Type& globalRhs)
    {
        for (int j = 0; j < M_ionicModelPtr->Size(); j++)
        {
            (* (M_globalRhs.at (j) ) ) = (* (globalRhs.at (j) ) );
        }
    }

    //! set the pointer to the fiber direction vector
    /*!
     @param fiberPtr pointer to the fiber direction vector
     */
    inline void setFiberPtr (const vectorPtr_Type fiberPtr)
    {
        this->M_fiberPtr = fiberPtr;
    }

    //! set the fiber direction vector
    /*!
     @param fiber  the fiber direction vector
     */
    inline void setFiber (const vector_Type& fiber)
    {
        (* (M_fiberPtr) ) = fiber;
    }

    //! set the the choice of lumping
    /*!
     @param isLumped true if you want to lump the mass matrix
     */
    inline void setLumpedMassMatrix (bool isLumped)
    {
        M_lumpedMassMatrix = isLumped;
    }

    //@}

    //! @name Methods
    //@{

    //! setup method used in the constructor
    /*!
     @param dataFile needed to set up the preconditioner
     @param ionicSize number of equation in the ionic model
     */
    virtual void setup (GetPot& dataFile, short int ionicSize);

    //! setup method used in the constructor
    /*!
     @param meshFile filename of the mesh
     @param meshPath directory where we have the mesh
     @param dataFile needed to set up the preconditioner
     @param ionicSize number of equation in the ionic model
     */
    virtual void setup (std::string meshName, std::string meshPath, GetPot& dataFile,
                short int ionicSize);

    //! create mass matrix
    /*!
     * Computes the mass matrix calling different methods if the mass should be lumped
     */
    void setupMassMatrix(); //! create mass matrix

    //! create mass matrix
    /*!
     * Computes the lumped mass matrix by nodal integration.
     */
    void setupLumpedMassMatrix();

    //! create stiffness matrix
    /*!
      * Computes the stiffness matrix calling different methods
      */
    virtual void setupStiffnessMatrix();

    //! create stiffness matrix given a diagonal diffusion tensor
    /*!
     @param diffusion vector defining the conductivities
     */
    void setupStiffnessMatrix (VectorSmall<3> diffusion);

    //! setup the total matrix
    /*!
     *  \f[
     *  A = C_m\frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     *
     *  where \fC_m\f is the membrane capacitance, \f\Delta t\f the timestep,
     *  M is the mass matrix and K is the stiffness matrix which depends on
     *  the fiber direction.
     */
    void setupGlobalMatrix();

    //! setup the linear solver
    /*!
     * A file named MonodomainSolverParamList.xml must be in the execution folder
     * with the parameters to set the linear solver
     @param dataFile GetPot to setup the preconditioners
     */
    void setupLinearSolver (GetPot dataFile);

    //! Initialize the potential to the value k
    /*!
     @param k value to intialize the potential
     */
    void inline initializePotential (Real k = 0.0)
    {
        (*M_potentialPtr) = k;
    }

    //! Initialize the applied current to the value k
    /*!
     @param k value to intialize the applied current
     */
    void inline initializeAppliedCurrent (Real k = 0.0)
    {
        (* (M_ionicModelPtr->appliedCurrentPtr() ) ) = k;
    }

    //! creates a vector of pointers to store the solution
    /*!
     * The first pointer points to the vector of the transmembrane potential,
     * while the others point to the gating variables
     @param ionicSize number of unknowns in the ionic model
     */
    void setupGlobalSolution (short int ionicSize);

    //! creates a vector of pointers to store the rhs
    /*!
     * The first pointer points to the rhs of the transmembrane potential,
     * while the others point to the rhs of the gating variables
     @param ionicSize number of unknowns in the ionic model
     */
    void setupGlobalRhs (short int ionicSize);

    //! Set parameters from an xml file
    /*!
     @param list Teuchos parameter list with the monodomain parameters
     */
    void setParameters (list_Type list);

    //! partition the mesh
    /*!
     @param meshName filename of the mesh
     @param meshPath path to the folder where you have the mesh
     */
    void inline partitionMesh (std::string meshName, std::string meshPath)
    {
        MeshUtility::loadMesh (M_localMeshPtr, M_fullMeshPtr, meshName, meshPath);
    }

    //! given a boost function initialize the potential
    /*!
     @param f function defining the potential
     @param time time at which we want to evaluate the function
     */
    void inline setPotentialFromFunction (function_Type& f, Real time = 0.0)
    {
        M_feSpacePtr->interpolate (
            static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
            *M_potentialPtr, time);
    }

    //! given a boost function initialize the applied current
    /*!
     @param f function defining the applied current
     @param time time at which we want to evaluate the function
     */
    void inline setAppliedCurrentFromFunction (function_Type& f,
                                               Real time = 0.0)
    {
        M_ionicModelPtr->setAppliedCurrentFromFunction (f, M_feSpacePtr, time);
    }

    //! given a ElectroStimulus object initialize the applied current
    /*!
     @param stimulus pacing protocol
     @param time time at which we want to evaluate the stimulus
     */
    void inline setAppliedCurrentFromElectroStimulus (ElectroStimulus& stimulus,
                                                      Real time = 0.0)
    {

        M_ionicModelPtr->setAppliedCurrentFromElectroStimulus (stimulus, M_feSpacePtr, time);
    }

    //! Solves one reaction step using the forward Euler scheme and N subiterations
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^{n+k/N} + \dfrac{\Delta t}{N} I_{ion}(\mathbf{V}^{n+k/N}), \quad \text{for } k=0,\dots,N-1.
     * \f]
     */
    /*!
     @param int number of subiterations
     */
    void solveOneReactionStepFE (int subiterations = 1);

    //! Solves one reaction step using the forward Euler scheme
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^{n+k/N} + \dfrac{\Delta t}{N} M I_{ion}(\mathbf{V}^{n+k/N}), \quad \text{for } k=0,\dots,N-1.
     * \f]
     */
    /*!
     @param matrix_Type full mass matrix
     @param int number of subiterations
     */
    void solveOneReactionStepFE (matrix_Type& mass, int subiterations = 1);

    //! Solves one reaction step using the Rush-Larsen scheme
    /*!
     @param int number of subiterations
     */
    void solveOneReactionStepRL (int subiterations = 1);

    //! Update the rhs
    /*!
     * \f[
     * rhs \leftarrow C_m \frac{M}{\Delta t} \mathbf{V}^n
     * \f]
     */
    void inline updateRhs()
    {
        (*M_rhsPtrUnique) += (*M_massMatrixPtr) * (*M_potentialPtr)
                             * (M_ionicModelPtr -> membraneCapacitance() / (M_timeStep) );
    }

    //! Solves one diffusion step using the BDF2 scheme
    /*!
     * \f[
     * ( \frac{3}{\Delta t}M + A ) V^{n+1} = \frac{1}{\Delta t}M(4 V^n -V^{n-1})
     * \f]
     */
    /*!
     @param previousPotentialPtr potential at \f n-1 \f
     */
    void solveOneDiffusionStepBDF2 (vectorPtr_Type previousPotentialPtr);

    //! Solves one diffusion step using the backward Euler scheme
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( C_m \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    void solveOneDiffusionStepBE();

    //!Solve one full step with operator splitting
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    void solveOneSplittingStep();

    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep
    void solveSplitting();

    //!Solve one full step with operator splitting and export the solution
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^*.
     * \f]
     */
    /*!
     @param exporter where you want to save the solution
     @param  t time at which we save the solution
     */
    void solveOneSplittingStep (IOFile_Type& exporter, Real t);

    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    /*!
     @param exporter where you want to save the solution
     */
    void solveSplitting (IOFile_Type& exporter);

    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    /*!
     @param exporter where you want to save the solution
     @param  t time at which we save the solution
     */
    void solveSplitting (IOFile_Type& exporter, Real dt);

    //! add to a given exporter the pointer to the potential saved with name fileName
    /*!
     @param exporter where you want to save the solution
     @param  t time at which we save the solution
     */
    void setupPotentialExporter (IOFile_Type& exporter, std::string fileName = "Potential");

    //! add to a given exporter the pointer to the potential and to the gating variables saved with name fileName
    /*!
     @param exporter where you want to save the solution
     @param fileName name of the file we wish to export
     @param folder directory where to save the solution
     */
    void setupExporter (IOFile_Type& exporter, std::string fileName = "output",
                        std::string folder = "./");

    //! Generates a default fiber direction (0,1,0)
    void setupFibers();

    //! Generates the fiber direction given the three component of the vector (F_x,F_y,F_z)
    /*!
     @param fibers vector with the fiber direction
     */
    void setupFibers (VectorSmall<3> fibers);

    //! Imports the fiber direction from a hdf5 file
    /*!
     @param fibersFile name of the hdf5 file with the fibers
     */
    inline void setupFibers (std::string fibersFile, const std::string& filePath = "./")
    {
        ElectrophysiologyUtility::importFibers (M_fiberPtr, fibersFile, M_localMeshPtr, filePath);
    }

    //! Imports the fiber direction from a vtk file ( format = 0), or text file
    /*!
     @param fibersFile name of the file with the fibers
     @param directory folder in which we have the file for the fibers
     @param format format in which fibers are saved
    *
    * format 0 = fibers saved as (fx, fy, fz) in each row
    *
    * format 1 = fibers saved as fx in each row for all the mesh
    *                            fy in each row for all the mesh
    *                            fz in each row for all the mesh
     */
    inline void setupFibers (std::string fibersFile, std::string directory,
                             int format = 0)
    {
        ElectrophysiologyUtility::importFibersFromTextFile (M_fiberPtr, fibersFile,
                                                            directory, format);
    }

    //! Solves the gating variables with forward Euler
    void solveOneStepGatingVariablesFE();

    //! Solves the gating variables with Rush-Larsen scheme
    void solveOneStepGatingVariablesRL();

    //! Compute the rhs using state variable interpolation
    void computeRhsSVI();

    //! Compute the rhs using ionic current interpolation
    void computeRhsICI();

    //! Compute the rhs using ionic current interpolation
    /*!
     * This method is useful to solve ICI without lumping the mass matrix
     * in fron of the reaction term.
     * Lump the mass matrix, and pass as argument a full mass matrix
     */
    void computeRhsICIWithFullMass ();

    //!Solve one full step with ionic current interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     */
    virtual void solveOneICIStep();

    //! solves using ionic current interpolation
    /*!
     * This method is useful to solve ICI without lumping the mass matrix
     * in front of the reaction term.
     * Lump the mass matrix, and pass as argument a full mass matrix
     */
    void solveOneICIStepWithFullMass ();

    //!Solve one full step with state variable interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneSVIStep();

    //! solve system using ICI from M_initialTime to the M_endTime with time step M_timeStep
    void solveICI();

    //! solve system using SVI from M_initialTime to the M_endTime with time step M_timeStep
    void solveSVI();

    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveOneICIStep (IOFile_Type& exporter, Real t);

    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     * This method is useful to solve ICI without lumping the mass matrix
     * in front of the reaction term.
     * Lump the mass matrix, and pass as argument a full mass matrix
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveOneICIStepWithFullMass (IOFile_Type& exporter, Real t);

    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
     * \f]
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveOneSVIStep (IOFile_Type& exporter, Real t);

    //! solve system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    /*!
    @param exporter where we want to save the solution
     */
    void solveICI (IOFile_Type& exporter);

    //! solve system using SVI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    /*!
    @param exporter where we want to save the solution
     */
    void solveSVI (IOFile_Type& exporter);

    //!Solve the system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    /*!
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveICI (IOFile_Type& exporter, Real dt);

    //!Solve the system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    /*!
     * This method is useful to solve ICI without lumping the mass matrix
     * in front of the reaction term.
     * Lump the mass matrix, and pass as argument a full mass matrix
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveICIWithFullMass (IOFile_Type& exporter, Real dt);

    //!Solve the using SVI from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    /*!
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void solveSVI (IOFile_Type& exporter, Real dt);

    //! Generates a file where the fiber direction is saved
    /*!
    @param postDir directory in which we save the fibers
     */
    void exportFiberDirection (std::string postDir = "./");

    //! save the fiber direction into the given exporter
    /*!
    @param exporter where we want to save the solution
     */
    void exportFiberDirection (IOFile_Type& exporter);

    //! Save the solution in the exporter
    /*!
    @param exporter where we want to save the solution
    @param t time at wich we save the solution
     */
    void inline exportSolution (IOFile_Type& exporter, Real t)
    {
        exporter.postProcess (t);
    }

    //! Importer: lead a solution from an hdf5 file
    /*!
    @param prefix name of the hdf5 file to import
    @param postDir folder where the file to import is
    @param time  time at wich we want to import the solution
     */
    void importSolution (std::string prefix, std::string postDir, Real time = 0.0);

    //! Initialize the solution with resting values of the ionic model
    void inline setInitialConditions()
    {
        M_ionicModelPtr->initialize (M_globalSolution);
    }

    //! save the fiber direction into the given exporter
    /*!
     * The activationTimeVector is required to be initialized with negative values
     *
    @param activationTimeVector vector where we register the activation time
    @param time  time at which we are
    @param threshold value for which we consider activation
     */
    void registerActivationTime (vector_Type& activationTimeVector, Real time,
                                 Real threshold = 0.0);
    //! set the verbosity
    /*!
     *
    @param verbose show additional output
     */
    void setVerbosity (bool verbose);

    //@}

private:

    //! Set default parameters
    void setParameters();

    //! initialization in constructor
    void init();

    //! initialization in constructor
    /*!
    @param comm epetra communicator
     */
    void init (commPtr_Type comm);

    //! initialization in constructor
    /*!
    @param meshPtr pointer to the mesh
     */
    void init (meshPtr_Type meshPtr);

    //! initialization in constructor
    /*!
    @param model pointer to the ionic model
     */
    void init (ionicModelPtr_Type model);

protected:
    //surface to volume ration
    Real M_surfaceVolumeRatio;
    //ionic model
    ionicModelPtr_Type M_ionicModelPtr;
    //Epetra communicator
    commPtr_Type M_commPtr;
    //partitioned mesh
    meshPtr_Type M_localMeshPtr;
    //non partitioned mesh
    meshPtr_Type M_fullMeshPtr;
    //ET finite element space
    ETFESpacePtr_Type M_ETFESpacePtr;
    //finite element space
    feSpacePtr_Type M_feSpacePtr;
    //mass matrix
    matrixPtr_Type M_massMatrixPtr;
    //full mass matrix
    matrixPtr_Type M_fullMassMatrixPtr;
    //stiffness matrix
    matrixPtr_Type M_stiffnessMatrixPtr;
    //mass matrix / timestep + stiffness matrix
    matrixPtr_Type M_globalMatrixPtr;
    //initial time t_0
    Real M_initialTime;
    //final time t_f
    Real M_endTime;
    //timestep dt
    Real M_timeStep;
    // conductivity tensor
    VectorSmall<3> M_diffusionTensor;
    //rhs of the potential equation
    vectorPtr_Type M_rhsPtr;
    //rhs of the potential equation - unique version
    vectorPtr_Type M_rhsPtrUnique;
    //potential
    vectorPtr_Type M_potentialPtr;
    //linear solver
    linearSolverPtr_Type M_linearSolverPtr;
    //vector of pointers pointing to the potential
    // and to the other gating variables
    vectorOfPtr_Type M_globalSolution;
    //vector of pointers pointing to rhs of the potential
    // and to the ths of the other gating variables
    vectorOfPtr_Type M_globalRhs;
    //order of the elements (P1)
    std::string M_elementsOrder;
    //fiber field
    vectorPtr_Type M_fiberPtr;
    //Create the identity matrix I
    matrixSmall_Type M_identity;
    //using lumped mass matrix
    bool            M_lumpedMassMatrix;
    //verbosity
    bool            M_verbose;

};
// class MonodomainSolver

//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
//!Empty constructor
template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>::ElectroETAMonodomainSolver()
{
    M_verbose = false;
    setParameters();
    init();
}

//!constructor
template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>::ElectroETAMonodomainSolver (
    std::string meshName, std::string meshPath, GetPot& dataFile,
    ionicModelPtr_Type model)
{
    M_verbose = false;
    setParameters();
    init (model);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr->Size() );
}

//!constructor with communicator
template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>::ElectroETAMonodomainSolver (
    std::string meshName, std::string meshPath, GetPot& dataFile,
    ionicModelPtr_Type model, commPtr_Type comm) :
    M_ionicModelPtr (model), M_verbose (false)
{
    setParameters();
    init (comm);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>::ElectroETAMonodomainSolver (
    GetPot& dataFile, ionicModelPtr_Type model, meshPtr_Type meshPtr) :
    M_ionicModelPtr (model), M_verbose (false)
{
    setParameters();
    init (meshPtr);
    setup (dataFile, M_ionicModelPtr->Size() );
}

//! Copy constructor
template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>::ElectroETAMonodomainSolver (
    const ElectroETAMonodomainSolver& solver) :
    M_surfaceVolumeRatio (solver.M_surfaceVolumeRatio),
    M_ionicModelPtr (new IonicModel (*solver.M_ionicModelPtr) ),
    M_commPtr (solver.M_commPtr),
    M_localMeshPtr ( solver.M_localMeshPtr),
    M_fullMeshPtr (solver.M_fullMeshPtr),
    M_ETFESpacePtr ( new ETFESpace_Type (*solver.M_ETFESpacePtr) ),
    M_feSpacePtr ( new feSpace_Type (*solver.M_feSpacePtr) ),
    M_massMatrixPtr ( new matrix_Type (* (solver.M_massMatrixPtr) ) ),
    M_stiffnessMatrixPtr ( new matrix_Type (* (solver.M_stiffnessMatrixPtr) ) ),
    M_globalMatrixPtr ( new matrix_Type (* (solver.M_globalMatrixPtr) ) ),
    M_initialTime ( solver.M_initialTime),
    M_endTime (solver.M_endTime),
    M_timeStep ( solver.M_timeStep),
    M_diffusionTensor (solver.M_diffusionTensor),
    M_rhsPtr ( new vector_Type (* (solver.M_rhsPtr) ) ),
    M_rhsPtrUnique (  new vector_Type (* (M_rhsPtr), Unique) ),
    M_potentialPtr ( new vector_Type (solver.M_ETFESpacePtr->map() ) ),
    M_linearSolverPtr ( new LinearSolver (* (solver.M_linearSolverPtr) ) ),
    M_elementsOrder ( solver.M_elementsOrder),
    M_fiberPtr ( new vector_Type (* (solver.M_fiberPtr) ) ) ,
    M_lumpedMassMatrix (solver.M_lumpedMassMatrix),
    M_verbose (solver.M_verbose),
    M_identity (solver.M_identity)
{
    setupGlobalSolution (M_ionicModelPtr->Size() );
    setGlobalSolution (solver.M_globalSolution);
    setupGlobalRhs (M_ionicModelPtr->Size() );
    setGlobalRhs (solver.M_globalRhs);
}

//! Assignment operator
template<typename Mesh, typename IonicModel>
ElectroETAMonodomainSolver<Mesh, IonicModel>& ElectroETAMonodomainSolver < Mesh,
                           IonicModel >::operator= (const ElectroETAMonodomainSolver& solver)
{
    if (M_verbose && M_commPtr -> MyPID() == 0)
    {
        std::cout << "\n WARNING!!! ETA Monodomain Solver: you are using the assignment operator! This method is outdated.";
        std::cout << "\n WARNING!!! ETA Monodomain Solver: Don't count on it, at this moment. Feel free to update  yourself...";
    }
    M_surfaceVolumeRatio = solver.M_surfaceVolumeRatio;
    setIonicModel ( (*solver.M_ionicModelPtr) );
    M_commPtr = solver.M_commPtr;
    M_localMeshPtr = solver.M_localMeshPtr;
    M_fullMeshPtr = solver.M_fullMeshPtr;
    setETFESpace (* (solver.M_ETFESpacePtr) );
    setFeSpace (* (solver.M_feSpacePtr) );
    setMassMatrix (* (solver.M_massMatrixPtr) );
    setStiffnessMatrix (* (solver.M_stiffnessMatrixPtr) );
    setGlobalMatrix (* (solver.M_globalMatrixPtr) );
    M_initialTime = solver.M_initialTime;
    M_endTime = solver.M_endTime;
    M_timeStep = solver.M_timeStep;
    M_diffusionTensor = solver.M_diffusionTensor;
    setRhs (* (solver.M_rhsPtr) );
    setRhsUnique (* (solver.M_rhsPtrUnique) );
    setPotential (* (solver.M_potentialPtr) );

    setLinearSolver (* (solver.M_linearSolverPtr) );
    setGlobalSolution (solver.M_globalSolution);
    setGlobalRhs (solver.M_globalRhs);
    M_elementsOrder = solver.M_elementsOrder;
    setFiber (* (solver.M_fiberPtr) );
    M_verbose = solver.M_verbose;
    M_identity = solver.M_identity;

    return *this;
}

/********* SETUP METHODS */ //////
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupFibers()
{

    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Space3D (
        new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
                                           3, M_commPtr) );

    M_fiberPtr.reset (new vector_Type (Space3D->map() ) );

    int d1 = (*M_fiberPtr).epetraVector().MyLength() / 3;
    (*M_fiberPtr) *= 0;
    int j (0);
    for (int k (0); k < d1; k++)
    {
        j = (*M_fiberPtr).blockMap().GID (k + d1);
        (*M_fiberPtr) [j] = 1.0;
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupFibers (
    VectorSmall<3> fibers)
{
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Space3D (
        new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
                                           3, M_commPtr) );

    M_fiberPtr.reset (new vector_Type (Space3D->map() ) );

    ElectrophysiologyUtility::setupFibers (*M_fiberPtr, fibers);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::exportFiberDirection (
    std::string postDir)
{
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Space3D (
        new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
                                           3, M_commPtr) );

    ExporterHDF5 < mesh_Type > exp;
    exp.setMeshProcId (M_localMeshPtr, M_commPtr->MyPID() );
    exp.setPostDir (postDir);
    exp.setPrefix ("FiberDirection");
    exp.addVariable (ExporterData<mesh_Type>::VectorField, "fibers", Space3D,
                     M_fiberPtr, UInt (0) );
    exp.postProcess (0);
    exp.closeFile();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::exportFiberDirection (
    IOFile_Type& exporter)
{
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Space3D (
        new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
                                           3, M_commPtr) );

    exporter.addVariable (ExporterData<mesh_Type>::VectorField, "fibers",
                          Space3D, M_fiberPtr, UInt (0) );
    exporter.postProcess (0);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::importSolution (std::string prefix, std::string postDir, Real time)
{
    IOFilePtr_Type importer (new hdf5IOFile_Type() );
    importer->setPrefix (prefix);
    importer->setPostDir (postDir);

    importer->setMeshProcId (M_feSpacePtr->mesh(),
                             M_feSpacePtr->map().comm().MyPID() );

    importer->addVariable (IOData_Type::ScalarField, "Variable0", M_feSpacePtr,
                           M_potentialPtr, static_cast<UInt> (0) );
    for (int i (1); i < M_ionicModelPtr->Size(); i++)
        importer->addVariable (IOData_Type::ScalarField,
                               "Variable" + number2string (i), M_feSpacePtr,
                               M_globalSolution.at (i), static_cast<UInt> (0) );
    importer->importFromTime (time);
    importer->closeFile();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setup (GetPot& dataFile,
                                                          short int ionicSize)
{
    M_feSpacePtr.reset ( new feSpace_Type (M_localMeshPtr, M_elementsOrder, 1, M_commPtr) );

    M_ETFESpacePtr.reset ( new ETFESpace_Type (M_localMeshPtr, & (M_feSpacePtr -> refFE() ) , M_commPtr) );

    M_massMatrixPtr.reset (new matrix_Type (M_ETFESpacePtr->map() ) );

    M_stiffnessMatrixPtr.reset (new matrix_Type (M_ETFESpacePtr->map() ) );

    M_globalMatrixPtr.reset (new matrix_Type (M_ETFESpacePtr->map() ) );

    M_rhsPtr.reset (new vector_Type (M_ETFESpacePtr->map(), Repeated) );

    M_rhsPtrUnique.reset (new vector_Type (* (M_rhsPtr), Unique) );

    M_potentialPtr.reset (new vector_Type (M_ETFESpacePtr->map() ) );

    //***********************//
    //  Setup Linear Solver  //
    //***********************//
    setupLinearSolver (dataFile);

    //**************************//
    //  Setup Initial condition //
    //**************************//
    initializePotential (0.0);

    vector_Type Iapp (M_feSpacePtr->map() );
    Iapp *= 0.0;

    M_ionicModelPtr->setAppliedCurrent (Iapp);

    setupGlobalSolution (ionicSize);

    setupGlobalRhs (ionicSize);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setup (std::string meshName,
                                                          std::string meshPath,
                                                          GetPot& dataFile,
                                                          short int ionicSize)
{
    //partitioning the mesh
    partitionMesh (meshName, meshPath);
    setup (dataFile, ionicSize);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupMassMatrix()
{
    if (M_lumpedMassMatrix)
    {
        setupLumpedMassMatrix();
    }
    else
    {
        *M_massMatrixPtr *= 0.0;
        if (M_verbose && M_localMeshPtr->comm()->MyPID() == 0)
        {
            std::cout << "\nETA Monodomain Solver: Setting up mass matrix";
        }

        {
            using namespace ExpressionAssembly;

            integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(),
                       M_ETFESpacePtr, M_ETFESpacePtr, phi_i * phi_j)
                    >> M_massMatrixPtr;

        }
        M_massMatrixPtr->globalAssemble();
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupLumpedMassMatrix()
{

    M_lumpedMassMatrix = true;
    *M_massMatrixPtr *= 0.0;
    if (M_verbose && M_localMeshPtr->comm()->MyPID() == 0)
    {
        std::cout << "\nETA Monodomain Solver: Setting up lumped mass matrix";
    }
    {
        using namespace ExpressionAssembly;

        integrate (elements (M_localMeshPtr), quadRuleTetra4ptNodal,
                   M_ETFESpacePtr, M_ETFESpacePtr, phi_i * phi_j)
                >> M_massMatrixPtr;

    }
    M_massMatrixPtr->globalAssemble();
}


template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupStiffnessMatrix()
{
    setupStiffnessMatrix (M_diffusionTensor);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupStiffnessMatrix (
    VectorSmall<3> diffusion)
{
    if (M_verbose && M_localMeshPtr->comm()->MyPID() == 0)
    {
        std::cout
                << "\nETA Monodomain Solver: Setting up stiffness matrix (only fiber field)";
    }

    *M_stiffnessMatrixPtr *= 0.0;
    Real sigmal = M_diffusionTensor[0];
    Real sigmat = M_diffusionTensor[1];

    ETFESpaceVectorialPtr_Type spaceVectorial (
        new ETFESpaceVectorial_Type (M_localMeshPtr, & (M_feSpacePtr -> refFE() ), M_commPtr) );

    {
        using namespace ExpressionAssembly;

        auto I = value(M_identity);
        auto f0 = value (spaceVectorial, *M_fiberPtr);
        auto D = value (sigmat) * I + (value (sigmal) - value (sigmat) ) * outerProduct (f0, f0);

        integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                   M_ETFESpacePtr,
                   dot ( D * grad (phi_i), grad (phi_j) ) )
                >> M_stiffnessMatrixPtr;

    }

    M_stiffnessMatrixPtr->globalAssemble();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupGlobalMatrix()
{
    (*M_globalMatrixPtr) *= 0;
    (*M_globalMatrixPtr) = (*M_stiffnessMatrixPtr);
    (*M_globalMatrixPtr) *= 1.0 / M_surfaceVolumeRatio;
    (*M_globalMatrixPtr) += ( (*M_massMatrixPtr) * ( M_ionicModelPtr -> membraneCapacitance() / M_timeStep) );
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupLinearSolver (
    GetPot dataFile)
{
    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot (dataFile, "prec");
    precPtr.reset (precRawPtr);

    Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp ( new Teuchos::ParameterList);

    std::string xmlpath = dataFile ("electrophysiology/monodomain_xml_path", "./");
    std::string xmlfile = dataFile ("electrophysiology/monodomain_xml_file", "MonodomainSolverParamList.xml");

    solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);

    M_linearSolverPtr->setCommunicator (M_commPtr);
    M_linearSolverPtr->setParameters (*solverParamList);
    M_linearSolverPtr->setPreconditioner (precPtr);
    M_linearSolverPtr->setOperator (M_globalMatrixPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupGlobalSolution (
    short int ionicSize)
{
    M_globalSolution.push_back (M_potentialPtr);
    for (int k = 1; k < ionicSize; ++k)
    {
        M_globalSolution.push_back (
            * (new vectorPtr_Type (new VectorEpetra (M_ETFESpacePtr->map() ) ) ) );
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupGlobalRhs (
    short int ionicSize)
{
    M_globalRhs.push_back (M_rhsPtrUnique);
    for (int k = 1; k < ionicSize; ++k)
    {
        M_globalRhs.push_back (
            * (new vectorPtr_Type (new VectorEpetra (M_ETFESpacePtr->map() ) ) ) );
    }
}

/************** EXPORTER */    //////////////
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupPotentialExporter (
    IOFile_Type& exporter, std::string fileName)
{
    exporter.setMeshProcId (M_localMeshPtr, M_commPtr->MyPID() );
    exporter.setPrefix (fileName);
    exporter.exportPID (M_localMeshPtr, M_commPtr);
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "Potential",
                          M_feSpacePtr, M_potentialPtr, UInt (0) );
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setupExporter (
    IOFile_Type& exporter, std::string fileName, std::string folder)
{
    exporter.setMeshProcId (M_localMeshPtr, M_commPtr->MyPID() );
    exporter.setPrefix (fileName);
    exporter.exportPID (M_localMeshPtr, M_commPtr);
    exporter.setPostDir (folder);
    std::string variableName;
    for (int i = 0; i < M_ionicModelPtr->Size(); i++)
    {
        variableName = "Variable" + boost::lexical_cast<std::string> (i);
        exporter.addVariable (ExporterData<mesh_Type>::ScalarField, variableName,
                              M_feSpacePtr, M_globalSolution.at (i), UInt (0) );
    }
}

/********* SOLVING METHODS */    ////////////////////////
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneReactionStepFE (
    int subiterations)
{
    M_ionicModelPtr->superIonicModel::computeRhs (M_globalSolution, M_globalRhs);

    for (int i = 0; i < M_ionicModelPtr->Size(); i++)
    {
        if (i == 0)
            * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                           + ( (M_timeStep) / subiterations / M_ionicModelPtr -> membraneCapacitance() ) * (* (M_globalRhs.at (i) ) );
        else
            * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                           + ( (M_timeStep) / subiterations) * (* (M_globalRhs.at (i) ) );
    }
}


template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneReactionStepFE (matrix_Type& mass, int subiterations)
{
    M_ionicModelPtr->superIonicModel::computeRhs (M_globalSolution, M_globalRhs);

    for (int i = 0; i < M_ionicModelPtr->Size(); i++)
    {
        if (i == 0)
        {

            vector_Type aux ( M_potentialPtr -> map() );
            aux = mass.operator * ( (* (M_globalRhs.at (i) ) ) );
            * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                           + ( (M_timeStep) / subiterations / M_ionicModelPtr -> membraneCapacitance() ) * aux;
        }
        else
            * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                           + ( (M_timeStep) / subiterations) * (* (M_globalRhs.at (i) ) );
    }
}



template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneReactionStepRL (
    int subiterations)
{
    M_ionicModelPtr->superIonicModel::computeRhs (M_globalSolution, M_globalRhs);

    * (M_globalSolution.at (0) ) = * (M_globalSolution.at (0) )
                                   + ( (M_timeStep) / subiterations / M_ionicModelPtr -> membraneCapacitance() ) * (* (M_globalRhs.at (0) ) );

    M_ionicModelPtr->superIonicModel::computeGatingVariablesWithRushLarsen (
        M_globalSolution, M_timeStep / subiterations);
    int offset = M_ionicModelPtr->numberOfGatingVariables() + 1;
    for (int i = offset; i < M_ionicModelPtr->Size(); i++)
    {
        * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                       + ( (M_timeStep) / subiterations) * (* (M_globalRhs.at (i) ) );
    }

}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneDiffusionStepBDF2 (
    vectorPtr_Type previousPotentialPtr)
{
    matrixPtr_Type OperatorBDF2 (new matrix_Type (M_feSpacePtr->map() ) );
    vectorPtr_Type rhsBDF2 (new vector_Type (M_feSpacePtr->map() ) );
    (*OperatorBDF2) *= 0;
    (*OperatorBDF2) = (*M_stiffnessMatrixPtr);
    (*OperatorBDF2) += (*M_stiffnessMatrixPtr);
    (*OperatorBDF2) += ( (*M_massMatrixPtr) * (3.0 / M_timeStep) );
    (*rhsBDF2) *= 0;
    (*rhsBDF2) = 4.0 * (*M_rhsPtrUnique);
    (*rhsBDF2) -= (*M_massMatrixPtr) * (*previousPotentialPtr)
                  * (1.0 / M_timeStep);
    M_linearSolverPtr->setOperator (OperatorBDF2);
    M_linearSolverPtr->setRightHandSide (rhsBDF2);
    M_linearSolverPtr->solve (M_potentialPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneDiffusionStepBE()
{
    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
    M_linearSolverPtr->solve (M_potentialPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneSplittingStep()
{
    solveOneReactionStepFE();
    (*M_rhsPtrUnique) *= 0;
    updateRhs();
    solveOneDiffusionStepBE();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneSplittingStep (
    IOFile_Type& exporter, Real t)
{
    solveOneSplittingStep();
    exportSolution (exporter, t);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSplitting()
{
    for (Real t = M_initialTime; t < M_endTime;)
    {
        t = t + M_timeStep;
        solveOneSplittingStep();
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSplitting (
    IOFile_Type& exporter)
{
    if (M_endTime > M_timeStep)
    {
        for (Real t = M_initialTime; t < M_endTime;)
        {
            t = t + M_timeStep;
            solveOneSplittingStep (exporter, t);
        }
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSplitting (
    IOFile_Type& exporter, Real dt)
{
    assert (
        dt >= M_timeStep
        && "Cannot save the solution for step smaller than the timestep!");
    int iter ( (dt / M_timeStep) + 1e-9);
    int k (0);
    if (M_endTime > M_timeStep)
    {
        for (Real t = M_initialTime; t < M_endTime;)
        {

            t += M_timeStep;
            k++;
            if (k % iter == 0)
            {
                solveOneSplittingStep (exporter, t);
            }
            else
            {
                solveOneSplittingStep();
            }

        }
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneStepGatingVariablesFE()
{
    M_ionicModelPtr->superIonicModel::computeGatingRhs (M_globalSolution,
                                                        M_globalRhs);

    for (int i = 1; i < M_ionicModelPtr->Size(); i++)
    {
        * (M_globalSolution.at (i) ) = * (M_globalSolution.at (i) )
                                       + M_timeStep * (* (M_globalRhs.at (i) ) );
    }
}
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneStepGatingVariablesRL()
{

    M_ionicModelPtr->superIonicModel::computeGatingVariablesWithRushLarsen (
        M_globalSolution, M_timeStep);
    M_ionicModelPtr->superIonicModel::computeNonGatingRhs (M_globalSolution,
                                                           M_globalRhs);
    int offset = M_ionicModelPtr->numberOfGatingVariables() + 1;
    for (int i = offset; i < M_ionicModelPtr->Size(); i++)
    {
        * (M_globalRhs[i]) *= M_timeStep;
        * (M_globalSolution[i]) += * (M_globalRhs[i]);
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::computeRhsICI()
{
    M_ionicModelPtr->superIonicModel::computePotentialRhsICI (M_globalSolution,
                                                              M_globalRhs, (*M_massMatrixPtr) );
    updateRhs();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::computeRhsICIWithFullMass ()
{
    if (M_fullMassMatrixPtr)
    {
        M_ionicModelPtr->superIonicModel::computePotentialRhsICI (M_globalSolution,
                                                                  M_globalRhs, *M_fullMassMatrixPtr);
    }
    else
    {
        assert (0 && "fullMassMatrix Pointer was not set! Use the computeRhsICI() method!");
    }
    updateRhs();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::computeRhsSVI()
{
    if (M_verbose && M_commPtr -> MyPID() == 0)
    {
        std::cout << "\nETA Monodomain Solver: updating rhs with SVI";
    }
    M_ionicModelPtr->superIonicModel::computePotentialRhsSVI (M_globalSolution,
                                                              M_globalRhs, (*M_feSpacePtr) );
    updateRhs();
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneICIStep()
{
    computeRhsICI();
    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
    M_linearSolverPtr->solve (M_potentialPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneICIStepWithFullMass ()
{
    computeRhsICIWithFullMass ();
    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
    M_linearSolverPtr->solve (M_potentialPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneSVIStep()
{
    computeRhsSVI();
    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
    M_linearSolverPtr->solve (M_potentialPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveICI()
{
    for (Real t = M_initialTime; t < M_endTime;)
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneICIStep();
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSVI()
{
    for (Real t = M_initialTime; t < M_endTime;)
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneSVIStep();
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneICIStep (
    IOFile_Type& exporter, Real t)
{
    solveOneICIStep();
    exportSolution (exporter, t);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneICIStepWithFullMass (
    IOFile_Type& exporter, Real t)
{
    solveOneICIStepWithFullMass();
    exportSolution (exporter, t);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveOneSVIStep (
    IOFile_Type& exporter, Real t)
{
    solveOneSVIStep();
    exportSolution (exporter, t);
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveICI (
    IOFile_Type& exporter)
{

    Real dt = M_timeStep;
    for (Real t = M_initialTime; t < M_endTime;)
    {
        t += M_timeStep;
        if (t > M_endTime)
        {
            M_timeStep = M_endTime - (t - dt);
        }
        solveOneStepGatingVariablesFE();
        solveOneICIStep (exporter, t);
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSVI (
    IOFile_Type& exporter)
{
    for (Real t = M_initialTime; t < M_endTime;)
    {
        t = t + M_timeStep;
        solveOneStepGatingVariablesFE();
        solveOneSVIStep (exporter, t);
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveICI (
    IOFile_Type& exporter, Real dt)
{
    assert (
        dt >= M_timeStep
        && "Cannot save the solution for step smaller than the timestep!");
    int iter ( (dt / M_timeStep) + 1e-9);
    int k (0);

    if (M_endTime > M_timeStep)
    {
        for (Real t = M_initialTime; t < M_endTime;)
        {

            t += M_timeStep;
            if (t > M_endTime)
            {
                M_timeStep = M_endTime - (t - dt);
            }
            k++;
            solveOneStepGatingVariablesFE();
            if (k % iter == 0)
            {
                solveOneICIStep (exporter, t);
            }
            else
            {
                solveOneICIStep();
            }

        }
    }
}


template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveICIWithFullMass (
    IOFile_Type& exporter, Real dt)
{
    assert (
        dt >= M_timeStep
        && "Cannot save the solution for step smaller than the timestep!");
    if (!M_fullMassMatrixPtr)
    {
        solveICI (exporter, dt);
    }
    else
    {
        int iter ( (dt / M_timeStep) + 1e-9);
        int k (0);

        if (M_endTime > M_timeStep)
        {
            for (Real t = M_initialTime; t < M_endTime;)
            {

                t += M_timeStep;
                if (t > M_endTime)
                {
                    M_timeStep = M_endTime - (t - dt);
                }
                k++;
                solveOneStepGatingVariablesFE();
                if (k % iter == 0)
                {
                    solveOneICIStepWithFullMass (exporter, t);
                }
                else
                {
                    solveOneICIStepWithFullMass( );
                }

            }
        }
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::solveSVI (
    IOFile_Type& exporter, Real dt)
{
    assert (
        dt >= M_timeStep
        && "Cannot save the solution for step smaller than the timestep!");
    int iter ( (dt / M_timeStep) + 1e-9);
    int k (0);
    if (M_endTime > M_timeStep)
    {
        for (Real t = M_initialTime; t < M_endTime;)
        {

            t += M_timeStep;
            k++;
            solveOneStepGatingVariablesFE();
            if (k % iter == 0)
            {
                solveOneSVIStep (exporter, t);
            }
            else
            {
                solveOneSVIStep();
            }

        }
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::registerActivationTime (
    vector_Type& activationTimeVector, Real time, Real threshold)
{
    int n1 = M_potentialPtr->epetraVector().MyLength();
    int i (0);
    for (int l (0); l < n1; l++)
    {
        i = M_potentialPtr->blockMap().GID (l);
        if (activationTimeVector[i] < 0 && (* (M_potentialPtr) ) [i] > threshold)
        {
            activationTimeVector[i] = time;
        }

    }
}

/********   INITIALIZITION FOR CONSTRUCTOR ****/    //////
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::init()
{
    M_linearSolverPtr.reset (new LinearSolver() );
    M_globalSolution = vectorOfPtr_Type();
    M_globalRhs = vectorOfPtr_Type();

    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;

    M_commPtr.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
    M_localMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::init (
    ionicModelPtr_Type model)
{
    init();
    M_ionicModelPtr = model;
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::init (commPtr_Type comm)
{
    init();
    M_localMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_commPtr = comm;
}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::init (meshPtr_Type meshPtr)
{
    init();
    //TODO change the meshPtr to pass the fullMeshPtr
    M_localMeshPtr = meshPtr;
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_commPtr = meshPtr->comm();
}

/********* parameter initialization */    ////////
template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setParameters()
{
    M_surfaceVolumeRatio = 1400.0;
    M_diffusionTensor[0] = 1.0;
    M_diffusionTensor[1] = 1.0;
    M_diffusionTensor[2] = 1.0;
    M_initialTime = 0.0;
    M_endTime = 100.0;
    M_timeStep = 0.01;
    M_elementsOrder = "P1";
    M_lumpedMassMatrix = false;

}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setParameters (
    list_Type list)
{
    M_surfaceVolumeRatio = list.get ("surfaceVolumeRatio", 1400.0);
    M_diffusionTensor[0] = list.get ("longitudinalDiffusion", 1.0);
    M_diffusionTensor[1] = list.get ("transversalDiffusion", 1.0);
    M_diffusionTensor[2] = M_diffusionTensor[1];
    M_initialTime = list.get ("initialTime", 0.0);
    M_endTime = list.get ("endTime", 100.0);
    M_timeStep = list.get ("timeStep", 0.01);
    M_elementsOrder = list.get ("elementsOrder", "P1");
    M_lumpedMassMatrix = list.get ("LumpedMass", false);

}

template<typename Mesh, typename IonicModel>
void ElectroETAMonodomainSolver<Mesh, IonicModel>::setVerbosity (
    bool verbose)
{
    M_verbose = verbose;
}
} // namespace LifeV

#endif //_MONODOMAINSOLVER_H_
