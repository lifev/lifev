//HEADER
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
 @brief Class for solving the Bidomain equations in electrophysiology.

 @date 09-10-2013
 @author Toni Lassila <toni.lassila@epfl.ch>
 @author Matthias Lange <m.lange@sheffield.ac.uk>

 @last update 10-10-2013

 This class provides interfaces to solve the bidomain equations
 ( reaction diffusion equation ) using the ETA framework.
 The solution can be performed using three different methods:
 -operator splitting method (at this point available only with forward Euler
 for the reaction step and backward Euler for the diffusion step. );
 -Ionic Currents Interpolation (at this point only forward Euler);
 -State Variable interpolation (at this point only forward Euler).
 */

#ifndef _ELECTROETABIDOMAINSOLVER_H_
#define _ELECTROETABIDOMAINSOLVER_H_

#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <string>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

#include <Epetra_LocalMap.h>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/MatrixSmall.hpp>

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/VectorEpetraStructured.hpp>
#include <lifev/core/array/MatrixEpetraStructuredUtility.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/electrophysiology/solver/IonicModels/ElectroIonicModel.hpp>

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/electrophysiology/util/HeartUtility.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "Teuchos_XMLParameterListHelpers.hpp"

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/core/mesh/MeshLoadingUtility.hpp>

#include <lifev/core/algorithm/LinearSolver.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

#include <boost/typeof/typeof.hpp>

#include <lifev/core/fem/GradientRecovery.hpp>

namespace LifeV
{

//! BidomainSolver - Class featuring the usual solver for bidomain equations

template<typename Mesh, typename IonicModel>
class ElectroETABidomainSolver
{

    //!Bidomain Solver
    /*!
     The Bidomain equation reads


     */

public:

    //! @name Type definitions
    //@{

    typedef Mesh mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    typedef VectorEpetra vector_Type;
    typedef boost::shared_ptr<VectorEpetra> vectorPtr_Type;

    typedef VectorEpetraStructured blockVector_Type;
    typedef boost::shared_ptr<blockVector_Type> blockVectorPtr_Type;

    typedef std::vector<vectorPtr_Type> vectorOfPtr_Type;



    typedef MatrixEpetra<Real> matrix_Type;
    typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

    typedef MatrixEpetraStructured<Real> blockMatrix_Type;
    typedef boost::shared_ptr<blockMatrix_Type> blockMatrixPtr_Type;

    typedef Epetra_Comm comm_Type;
    typedef boost::shared_ptr<comm_Type> commPtr_Type;

    typedef ETFESpace<mesh_Type, MapEpetra, 3, 1> ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 1> > ETFESpacePtr_Type;

    typedef ETFESpace<mesh_Type, MapEpetra, 3, 3> ETFESpaceVectorial_Type;
    typedef boost::shared_ptr<ETFESpaceVectorial_Type> ETFESpaceVectorialPtr_Type;

    typedef FESpace<mesh_Type, MapEpetra> feSpace_Type;
    typedef boost::shared_ptr<feSpace_Type> feSpacePtr_Type;

    typedef LinearSolver linearSolver_Type;
    typedef boost::shared_ptr<LinearSolver> linearSolverPtr_Type;


    typedef Exporter<mesh_Type> IOFile_Type;
    typedef boost::shared_ptr<IOFile_Type> IOFilePtr_Type;
    typedef ExporterData<mesh_Type> IOData_Type;
    typedef ExporterEnsight<mesh_Type> ensightIOFile_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5< mesh_Type > hdf5IOFile_Type;
#endif

    typedef LifeV::Preconditioner basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type> basePrecPtr_Type;
    typedef LifeV::PreconditionerML prec_Type;
    typedef boost::shared_ptr<prec_Type> precPtr_Type;

    typedef IonicModel ionicModel_Type;
    typedef ElectroIonicModel superIonicModel;
    typedef boost::shared_ptr<ionicModel_Type> ionicModelPtr_Type;

    typedef Teuchos::ParameterList list_Type;

    typedef boost::function <
    Real (const Real& t, const Real& x, const Real& y, const Real& z,
          const ID& i) > function_Type;

    typedef MatrixSmall<3, 3> matrixSmall_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //!Empty Constructor
    /*!
     */
    ElectroETABidomainSolver();

    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot  datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     */

    ElectroETABidomainSolver (list_Type list, GetPot& dataFile,
                              ionicModelPtr_Type model);

    //! Constructor
    /*!
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel>  chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */
    ElectroETABidomainSolver (GetPot& dataFile, ionicModelPtr_Type model,
                              meshPtr_Type meshPtr);
    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */

    ElectroETABidomainSolver (list_Type list, GetPot& dataFile,
                              ionicModelPtr_Type model, commPtr_Type comm);

    //! Constructor
    /*!
     * @param Teuchos::ParameterList parameter list
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Mesh> Pointer to the partitioned mesh
     */

    ElectroETABidomainSolver (list_Type list, GetPot& dataFile,
                              ionicModelPtr_Type model, meshPtr_Type meshPtr);

    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     */

    ElectroETABidomainSolver (std::string meshName, std::string meshPath,
                              GetPot& dataFile, ionicModelPtr_Type model);



    //! Constructor
    /*!
     * @param string file name of the mesh
     * @param string path to the mesh
     * @param GetPot datafile (for preconditioner)
     * @param boost::shared_ptr<IonicModel> chosen ionic model pointer
     * @param boost::shared_ptr<Epetra_Comm> Epetra communicator
     */
    ElectroETABidomainSolver (std::string meshName, std::string meshPath,
                              GetPot& dataFile, ionicModelPtr_Type model, commPtr_Type comm);



    //! Copy Constructor
    /*!
     * @param ElectroETABidomainSolver object
     */
    ElectroETABidomainSolver (const ElectroETABidomainSolver& solver);

    //!Operator=()
    /*!
     * @param ElectroETABidomainSolver object
     */
    ElectroETABidomainSolver<Mesh, IonicModel>& operator= (
        const ElectroETABidomainSolver& solver);


    //! Destructor
    virtual ~ElectroETABidomainSolver()
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
    //! get the diagonal intra cellular diffusion tensor
    inline const VectorSmall<3>& diffusionTensorIntra() const
    {
        return M_diffusionTensorIntra;
    }
    //! get the diagonal extra cellular diffusion tensor
    inline const VectorSmall<3>& diffusionTensorExtra() const
    {
        return M_diffusionTensorExtra;
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
    inline const blockMatrixPtr_Type massMatrixPtr() const
    {
        return M_massMatrixPtr;
    }
    //! get the pointer to the stiffness matrix
    inline const blockMatrixPtr_Type stiffnessMatrixPtr() const
    {
        return M_stiffnessMatrixPtr;
    }

    //! get the pointer to the global matrix
    /*!
     *  \f[
     *  A = \frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     */
    inline const blockMatrixPtr_Type globalMatrixPtr() const
    {
        return M_globalMatrixPtr;
    }
    //! get the pointer to the right hand side
    inline const blockVectorPtr_Type rhsPtr() const
    {
        return M_rhsPtr;
    }
    //! get the pointer to the unique version of the right hand side
    inline const blockVectorPtr_Type rhsPtrUnique() const
    {
        return M_rhsPtrUnique;
    }
    //! get the pointer to the transmembrane potential
    inline const vectorPtr_Type potentialTransPtr() const
    {
        return M_potentialTransPtr;
    }
    //! get the pointer to the  potential
    inline blockVectorPtr_Type potentialGlobalPtr()
    {
        return M_potentialGlobalPtr;
    }
    //! get the pointer to the extra cellular potential
    inline const vectorPtr_Type potentialExtraPtr() const
    {
        return M_potentialExtraPtr;
    }

    //! get the pointer to the fiber vector
    inline const vectorPtr_Type fiberPtr() const
    {
        return M_fiberPtr;
    }
    //! get the pointer to the applied intra cellular current vector
    inline const vectorPtr_Type appliedCurrentIntraPtr()
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
    //! get the pointer to the transmembrane potential
    inline vectorPtr_Type displacementPtr() const
    {
        return M_displacementPtr;
    }

    inline ETFESpaceVectorialPtr_Type displacementETFESpacePtr() const
    {
        return M_displacementETFESpacePtr;
    }

    inline bool lumpedMassMatrix() const
    {
        return M_lumpedMassMatrix;
    }
    //@}

    //! @name Set Methods
    //@{


    //! set the surface to volume ratio
    /*!
     @param Real surface to volume ratio
     */
    inline void setSurfaceVolumeRatio (const Real& p)
    {
        this->M_surfaceVolumeRatio = p;
    }

    //! set the starting time
    /*!
     @param Real initial time
     */
    inline void setInitialTime (const Real& p)
    {
        this->M_initialTime = p;
    }
    //! set the ending time
    /*!
     @param Real ending time
     */
    inline void setTimeStep (const Real& p)
    {
        this->M_timeStep = p;
    }
    //! set the time step
    /*!
     @param Real time step
     */
    inline void setEndTime (const Real& p)
    {
        this->M_endTime = p;
    }
    //! set the intra cellular diagonal diffusion tensor
    /*!
     @param  VectorSmall<3> diagonal intra cellular diffusion tensor
     */
    inline void setDiffusionTensorIntra (const VectorSmall<3>& p)
    {
        this->M_diffusionTensorIntra = p;
    }
    //! set the extra cellular diagonal diffusion tensor
    /*!
     @param  VectorSmall<3> diagonal extra cellular diffusion tensor
     */
    inline void setDiffusionTensorExtra (const VectorSmall<3>& p)
    {
        this->M_diffusionTensorExtra = p;
    }

    //! set the pointer to the ionic model
    /*!
     @param boost::shared_ptr<IonicModel> pointer to the ionic model
     */
    inline void setIonicModelPtr (const ionicModelPtr_Type p)
    {
        this->M_ionicModelPtr = p;
    }
    //! set the pointer to the Epetra communicator
    /*!
     @param boost::shared_ptr<Epetra_Comm> pointer to the Epetra communicator
     */

    inline void setCommPtr (const commPtr_Type p)
    {
        this->M_commPtr = p;
    }
    //! set the pointer to the partitioned mesh
    /*!
     @param boost::shared_ptr<Mesh> pointer to the partitioned mesh
     */
    inline void setLocalMeshPtr (const meshPtr_Type p)
    {
        this->M_localMeshPtr = p;
    }
    //! set the pointer to the partitioned mesh
    /*!
     @param boost::shared_ptr<Mesh> pointer to the partitioned mesh
     */
    inline void setFullMeshPtr (const meshPtr_Type p)
    {
        this->M_fullMeshPtr = p;
    }
    //! set the pointer to the ETA fe space
    /*!
     @param  boost::shared_ptr<ETFESpace<Mesh,MapEpetra,3,1>> pointer to the ETA fe space
     */
    inline void setETFESpacePtr (const ETFESpacePtr_Type p)
    {
        this->M_ETFESpacePtr = p;
    }
    //! set the pointer to the usual fe space
    /*!
     @param boost::shared_ptr<IFESpace<Mesh,MapEpetra>> pointer to the usual fe space
     */
    inline void setFeSpacePtr (const feSpacePtr_Type p)
    {
        this->M_feSpacePtr = p;
    }

    //! set the pointer to the  mass matrix
    /*!
     @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the mass matrix
     */
    inline void setMassMatrixPtr (const blockMatrixPtr_Type p)
    {
        this->M_massMatrixPtr = p;
    }
    //! set the pointer to the stiffness matrix
    /*!
     @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the stiffness matrix
     */
    inline void setStiffnessMatrixPtr (const blockMatrixPtr_Type p)
    {
        this->M_stiffnessMatrixPtr = p;
    }
    //! set the pointer to the global matrix
    /*!
     @param boost::shared_ptr<MatrixEpetra<Real>> pointer to the global matrix
     */
    inline void setGlobalMatrixPtr (const blockMatrixPtr_Type p)
    {
        this->M_globalMatrixPtr = p;
    }
    //! set the pointer to the right hand side
    /*!
     @param boost::shared_ptr<VectorEpetra> pointer to the right hand side
     */
    inline void setRhsPtr (const blockVectorPtr_Type p)
    {
        this->M_rhsPtr = p;
    }
    //! set the pointer to the unique version of the right hand side
    /*!
     @param boost::shared_ptr<VectorEpetra>  pointer to the  unique version of the right hand side
     */
    inline void setRhsPtrUnique (const blockVectorPtr_Type p)
    {
        this->M_rhsPtrUnique = p;
    }

    inline void setPotentialExtraPtr (const vectorPtr_Type p)
    {
        this->M_potentialExtraPtr = p;
        this->M_potentialGlobalPtr->replace ( (*p), this->M_potentialGlobalPtr->block (1)->firstIndex);
    }

    inline void setPotentialTransPtr (const vectorPtr_Type p)
    {
        this->M_potentialTransPtr = p;
        this->M_potentialGlobalPtr->replace ( (*p), this->M_potentialGlobalPtr->block (0)->firstIndex() );
    }


    //! set the pointer to the intra cell applied current vector
    /*!
     @param boost::shared_ptr<VectorEpetra>  pointer to the intra cellular applied current vector
     */
    inline void setAppliedCurrentIntraPtr (const vectorPtr_Type p)
    {
        M_ionicModelPtr->setAppliedCurrentPtr (p);
    }

    //! set the pointer to the linear solver
    /*!
     @param boost::shared_ptr<LinearSolver> pointer to the linear solver
     */
    inline void setLinearSolverPtr (const linearSolverPtr_Type p)
    {
        this->M_linearSolverPtr = p;
    }
    //! set the vector of pointers containing the transmembrane potential (at 0) and the gating variables
    /*!
     @param std::vector<boost::shared_ptr<VectorEpetra>> vector of pointers containing the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalSolution (const vectorOfPtr_Type& p)
    {
        this->M_globalSolution = p;
    }
    //! set the vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
    /*!
     @param std::vector<boost::shared_ptr<VectorEpetra>> vector of pointers containing the rhs for the transmembrane potential (at 0) and the gating variables
     */
    inline void setGlobalRhs (const vectorOfPtr_Type& p)
    {
        this->M_globalRhs = p;
    }
    //! set the pointer to the fiber direction vector
    /*!
     @param boost::shared_ptr<VectorEpetra> pointer to the fiber direction vector
     */
    inline void setFiberPtr (const vectorPtr_Type p)
    {
        this->M_fiberPtr = p;
    }
    //! set the pointer to displacement of the tissue
    inline void setDisplacementPtr (const vectorPtr_Type p)
    {
        this->M_displacementPtr = p;
    }

    //@}

    //! @name Copy Methods
    //@{
    inline void setIonicModel (const ionicModel_Type& p)
    {
        (* (M_ionicModelPtr) ) = p;
    }
    inline void setComm (const comm_Type& p)
    {
        (* (M_commPtr) ) = p;
    }
    inline void setLocalMesh (const mesh_Type& p)
    {
        (* (M_localMeshPtr) ) = p;
    }
    inline void setFullMesh (const mesh_Type& p)
    {
        (* (M_fullMeshPtr) ) = p;
    }
    inline void setETFESpace (const ETFESpace_Type& p)
    {
        (* (M_ETFESpacePtr) ) = p;
    }
    inline void setFeSpace (const feSpace_Type& p)
    {
        (* (M_feSpacePtr) ) = p;
    }
    inline void setMassMatrix (const matrix_Type& p)
    {
        (* (M_massMatrixPtr) ) = p;
    }
    inline void setStiffnessMatrix (const matrix_Type& p)
    {
        (* (M_stiffnessMatrixPtr) ) = p;
    }
    inline void setGlobalMatrix (const matrix_Type& p)
    {
        (* (M_globalMatrixPtr) ) = p;
    }
    inline void setRhs (const vector_Type& p)
    {
        (* (M_rhsPtr) ) = p;
    }
    inline void setRhsUnique (const vector_Type& p)
    {
        (* (M_rhsPtrUnique) ) = p;
    }
    inline void setPotentialTrans (const vector_Type& p)
    {
        (* (M_potentialTransPtr) ) = p;
        this->M_potentialGlobalPtr->replace (p, this->M_potentialGlobalPtr->block (0)->firstIndex() );
    }
    inline void setPotentialExtra (const blockVector_Type& p)
    {
        (* (M_potentialExtraPtr) ) = p;
        this->M_potentialGlobalPtr->replace (p, this->M_potentialGlobalPtr->block (1)->firstIndex() );
    }
    inline void setFiber (const vector_Type& p)
    {
        (* (M_fiberPtr) ) = p;
    }
    inline void setDisplacement (const vector_Type& p)
    {
        (* (M_displacementPtr) ) = p;
    }
    inline void setAppliedCurrentIntra (const vector_Type& p)
    {
        M_ionicModelPtr->setAppliedCurrent (p);
    }
    inline void setLinearSolver (const linearSolver_Type& p)
    {
        (* (M_linearSolverPtr) ) = p;
    }
    inline void copyGlobalSolution (const vectorOfPtr_Type& p)
    {
        for (int j = 0; j < M_ionicModelPtr->Size(); j++)
        {
            (* (M_globalSolution.at (j) ) ) = (* (p.at (j) ) );
        }
    }

    inline void setVariablePtr (const vectorPtr_Type p, int j)
    {
        M_globalSolution.at (j) = p;
    }

    inline void copyGlobalRhs (const vectorOfPtr_Type& p)
    {
        for (int j = 0; j < M_ionicModelPtr->Size(); j++)
        {
            (* (M_globalRhs.at (j) ) ) = (* (p.at (j) ) );
        }
    }

    inline void setLumpedMassMatrix (bool p) const
    {
        M_lumpedMassMatrix = p;
    }


    //@}

    //! @name Methods
    //@{

    void setup (GetPot& dataFile, short int ionicSize);

    void setup (std::string meshName, std::string meshPath, GetPot& dataFile,
                short int ionicSize);

    //! create mass matrix
    void setupMassMatrix(); //! create mass matrix
    void setupMassMatrix (vector_Type& disp);
    //! create mass matrix
    void setupLumpedMassMatrix();
    void setupLumpedMassMatrix (vector_Type& disp);
    //! create stiffness matrix
    void setupStiffnessMatrix();
    //! create stiffness matrix in a moving domain
    void setupStiffnessMatrix (vectorPtr_Type disp);
    //! create stiffness matrix given a diagonal diffusion tensor
    void setupStiffnessMatrix (VectorSmall<3> diffusionIntra, VectorSmall<3> diffusionExtra);
    //! create stiffness matrix given the fiber direction and a diagonal diffusion tensor
    //void setupStiffnessMatrix (VectorEpetra& fiber, VectorSmall<3> diffusion);
    //! setup the total matrix
    /*!
     *  \f[
     *  A = \frac{M}{\Delta t} + K(\mathbf{f})
     *  \f]
     */
    void setupGlobalMatrix();
    //! setup the linear solver
    /*!
     * A file named BidomainSolverParamList.xml must be in the execution folder
     * with the parameters to set the linear solver
     */
    void setupLinearSolver (GetPot dataFile);
    //! setup the linear solver
    void setupLinearSolver (GetPot dataFile, list_Type list);
    //! Initialize the potentials to zero
    void inline initializePotential()
    {
        (*M_potentialTransPtr) *= 0;
        (*M_potentialExtraPtr) *= 0;
        (*M_potentialGlobalPtr) *= 0;
    }
    //! Initialize the potential to the value k
    void inline initializePotentialTrans (Real k)
    {
        (* (M_globalSolution.at (0) ) ) = k;
        this->M_potentialGlobalPtr->replace ( (* (M_globalSolution.at (0) ) ), this->M_potentialGlobalPtr->block (0)->firstIndex);

    }
    void inline initializePotentialExtra (Real k)
    {
        (*M_potentialExtraPtr) = k;
        this->M_potentialGlobalPtr->replace ( (*M_potentialExtraPtr), this->M_potentialGlobalPtr->block (1)->firstIndex);

    }

    //! Initialize the applied current to zero
    void inline initializeAppliedCurrentIntra()
    {
        (* (M_ionicModelPtr->appliedCurrentPtr() ) ) *= 0;
        //(*M_IappExtraPtr)*=0;
    }
    //! Initialize the intra cellular applied current to the value k
    void inline initializeAppliedCurrentIntra (Real k)
    {
        (* (M_ionicModelPtr->appliedCurrentPtr() ) ) = k;
    }
    //! creates a vector of pointers to store the solution
    /*!
     * The first pointer points to the vector of the transmembrane potential,
     * while the others point to the gating variables
     */
    void setupGlobalSolution (short int ionicSize);
    //! creates a vector of pointers to store the rhs
    /*!
     * The first pointer points to the rhs of the transmembrane potential,
     * while the others point to the rhs of the gating variables
     */
    void setupGlobalRhs (short int ionicSize);
    //! Set parameters from an xml file
    void setParameters (list_Type list);
    //! partition the mesh
    void inline partitionMesh (std::string meshName, std::string meshPath)
    {
        MeshUtility::loadMesh (M_localMeshPtr, M_fullMeshPtr, meshName,
                               meshPath);
    }
    //! given a boost function initialize the transmembrane potential
    void inline setPotentialFromFunctionTrans (function_Type& f, Real time = 0.0)
    {
        M_feSpacePtr->interpolate (
            static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
            * (M_globalSolution.at (0) ), time);
    }
    //! given a boost function initialize the extra cellular potential
    void inline setPotentialFromFunctionExtra (function_Type& f, Real time = 0.0)
    {
        M_feSpacePtr->interpolate (
            static_cast<FESpace<RegionMesh<LinearTetra>, MapEpetra>::function_Type> (f),
            *M_potentialExtraPtr, time);
    }
    //! given a boost function initialize the applied intra cellular current
    void inline setAppliedCurrentFromFunctionIntra (function_Type& f,
                                                    Real time = 0.0)
    {
        M_ionicModelPtr->setAppliedCurrentFromFunction (f, M_feSpacePtr, time);
    }
    //! Solves one reaction step using the forward Euler scheme
    /*!
     * \f[
     * \mathbf{V}^* = \mathbf{V}^n + \Delta t I_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneReactionStepFE (int subiterations = 1);
    void solveOneReactionStepFE (matrix_Type& mass, int subiterations = 1);
    void solveOneReactionStepRL (int subiterations = 1);

    //! Update the rhs
    /*!
     * \f[
     * rhs \leftarrow C_m \frac{M}{\Delta t} \mathbf{V}^n
     * \f]
     */
    void inline updateRhs()
    {
        blockVectorPtr_Type tmp ( new blockVector_Type ( *M_rhsPtrUnique ) );
        M_massMatrixPtr -> multiply (false, (*M_potentialGlobalPtr) * (M_ionicModelPtr -> membraneCapacitance() / (M_timeStep) ), *tmp );
        (*M_rhsPtrUnique) += (*tmp);
    }

    //! Solves one diffusion step using the BDF2 scheme
    /*!
     * \f[
     * ( \frac{3}{\Delta t}M + A ) V^{n+1} = \frac{1}{\Delta t}M(4 V^n -V^{n-1})
     * \f]
     */
    void solveOneDiffusionStepBDF2 (vectorPtr_Type previousPotentialGlobalPtr);

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
    void solveOneSplittingStep (IOFile_Type& exporter, Real t);
    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveSplitting (IOFile_Type& exporter);
    //!Solve the system with operator splitting from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    void solveSplitting (IOFile_Type& exporter, Real dt);
    //! add to a given exporter the pointer to the potential
    void setupPotentialExporter (IOFile_Type& exporter);
    //! add to a given exporter the pointer to the potential saved with name fileName
    void setupPotentialExporter (IOFile_Type& exporter, std::string fileName);
    //! add to a given exporter the pointer to the potential and to the gating variables saved with name fileName
    void setupExporter (IOFile_Type& exporter, std::string fileName = "output",
                        std::string folder = "./");

    //! Generates a default fiber direction (0,1,0)
    void setupFibers();
    //! Generates the fiber direction given the three component of the vector (F_x,F_y,F_z)
    void setupFibers (VectorSmall<3> fibers);
    //! Imports the fiber direction from a hdf5 file
    inline void setupFibers (std::string fibersFile)
    {
        ElectrophysiologyUtility::importFibers (M_fiberPtr, fibersFile, M_localMeshPtr);
    }
    //! Imports the fiber direction from a vtk file ( format = 0), or text file
    //
    inline void setupFibers (std::string fibersFile, std::string directory,
                             int format = 0)
    {
        ElectrophysiologyUtility::importFibersFromTextFile (M_fiberPtr, fibersFile,
                                                            directory, format);
    }
    //! Solves the gating variables with forward Euler
    void solveOneStepGatingVariablesFE();
    void solveOneStepGatingVariablesRL();
    //! Compute the rhs using state variable interpolation
    void computeRhsSVI();
    //! Compute the rhs using ionic current interpolation
    void computeRhsICI();
    void computeRhsICI (matrix_Type& mass);
    //!Solve one full step with ionic current interpolation
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+M\mathbf{I},
     * \f]
     * where $\mathbf{I}$ is the vector of the ionic currents $I_j = I_{ion}(V_j^n)$
     */
    void solveOneICIStep();
    void solveOneICIStep (matrix_Type& mass);
    //!Solve one full step with ionic current interpolation
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
     */
    void solveOneICIStep (IOFile_Type& exporter, Real t);
    //!Solve one full step with ionic current interpolation  and export the solution
    /*!
     * \f[
     * A\mathbf{V}^{n+1} = \left( \frac{M}{\Delta t} + K(\mathbf{f}) \right)\mathbf{V}^{n+1} =\frac{M}{\Delta t} \mathbf{V}^n+\mathbf{I}_{ion}(\mathbf{V}^n).
     * \f]
     */
    void solveOneSVIStep (IOFile_Type& exporter, Real t);

    //! solve system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveICI (IOFile_Type& exporter);
    //! solve system using SVI from M_initialTime to the M_endTime with time step M_timeStep and export the solution
    void solveSVI (IOFile_Type& exporter);
    //!Solve the system using ICI from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    void solveICI (IOFile_Type& exporter, Real dt);
    //!Solve the using SVI from M_initialTime to the M_endTime with time step M_timeStep and export the solution every dt
    void solveSVI (IOFile_Type& exporter, Real dt);
    //! Generates a file where the fiber direction is saved
    void exportFiberDirection (std::string postDir = "./");
    //! save the fiber direction into the given exporter
    void exportFiberDirection (IOFile_Type& exporter);
    //! Save the solution in the exporter
    void inline exportSolution (IOFile_Type& exporter, Real t)
    {
        exporter.postProcess (t);
    }
    //! Import solution
    void importSolution (GetPot& dataFile, std::string prefix,
                         std::string postDir, Real time);

    void inline setInitialConditions()
    {
        M_ionicModelPtr->initialize (M_globalSolution);
    }

    //! save the fiber direction into the given exporter
    void registerActivationTime (vector_Type& activationTimeVector, Real time,
                                 Real threshold = 0.0);

    //@}

private:

    void setParameters();
    void init();
    void init (commPtr_Type comm);
    void init (meshPtr_Type meshPtr);
    void init (ionicModelPtr_Type model);
    void init (commPtr_Type comm, ionicModelPtr_Type model);
    void init (meshPtr_Type meshPtr, ionicModelPtr_Type model);

    Real M_surfaceVolumeRatio;

    ionicModelPtr_Type M_ionicModelPtr;

    commPtr_Type M_commPtr;
    meshPtr_Type M_localMeshPtr;
    meshPtr_Type M_fullMeshPtr;
    ETFESpacePtr_Type M_ETFESpacePtr;
    feSpacePtr_Type M_feSpacePtr;

    blockMatrixPtr_Type M_massMatrixPtr;
    blockMatrixPtr_Type M_stiffnessMatrixPtr;
    blockMatrixPtr_Type M_globalMatrixPtr;

    Real M_initialTime;
    Real M_endTime;
    Real M_timeStep;

    VectorSmall<3> M_diffusionTensorIntra;
    VectorSmall<3> M_diffusionTensorExtra;

    blockVectorPtr_Type M_rhsPtr;
    blockVectorPtr_Type M_rhsPtrUnique;
    blockVectorPtr_Type M_potentialGlobalPtr;
    vectorPtr_Type M_potentialTransPtr;
    vectorPtr_Type M_potentialExtraPtr;

    linearSolverPtr_Type M_linearSolverPtr;

    vectorOfPtr_Type M_globalSolution;
    vectorOfPtr_Type M_globalRhs;

    std::string M_elementsOrder;

    vectorPtr_Type M_fiberPtr;

    vectorPtr_Type M_displacementPtr;

    //Create the identity for F
    matrixSmall_Type M_identity;

    ETFESpaceVectorialPtr_Type M_displacementETFESpacePtr;

    bool            M_lumpedMassMatrix;

protected:
    void importFibers (vectorPtr_Type M_fiberPtr, std::string fibersFile, meshPtr_Type M_localMeshPtr);
};



// class BidomainSolver

//
// IMPLEMENTATION
//
// ===================================================
//! Constructors
// ===================================================
template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver()
{
    setParameters();
    init();
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    list_Type list, GetPot& dataFile, ionicModelPtr_Type model)
{
    init (model);
    setParameters (list);
    setup (list.get ("meshName", "lid16.mesh"), list.get ("meshPath", "./"),
           dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    list_Type list, GetPot& dataFile, ionicModelPtr_Type model,
    commPtr_Type comm)
{
    setParameters (list);
    init (comm);
    setup (list.get ("meshName", "lid16.mesh"), list.get ("meshPath", "./"),
           dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    list_Type list, GetPot& dataFile, ionicModelPtr_Type model,
    meshPtr_Type meshPtr)
{
    setParameters (list);
    init (meshPtr);
    setup (dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    std::string meshName, std::string meshPath, GetPot& dataFile,
    ionicModelPtr_Type model)
{
    setParameters();
    init (model);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    std::string meshName, std::string meshPath, GetPot& dataFile,
    ionicModelPtr_Type model, commPtr_Type comm) :
    M_ionicModelPtr (model)
{
    setParameters();
    init (comm);
    setup (meshName, meshPath, dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    GetPot& dataFile, ionicModelPtr_Type model, meshPtr_Type meshPtr) :
    M_ionicModelPtr (model)
{
    setParameters();
    init (meshPtr);
    setup (dataFile, M_ionicModelPtr->Size() );
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>::ElectroETABidomainSolver (
    const ElectroETABidomainSolver& solver) :
    M_surfaceVolumeRatio (solver.M_surfaceVolumeRatio), M_ionicModelPtr (
        solver.M_ionicModelPtr), M_commPtr (solver.M_commPtr), M_localMeshPtr (
            solver.M_localMeshPtr), M_fullMeshPtr (solver.M_fullMeshPtr), M_ETFESpacePtr (
                solver.M_ETFESpacePtr), M_feSpacePtr (solver.M_feSpacePtr), M_massMatrixPtr (
                    new blockMatrix_Type (* (solver.M_massMatrixPtr) ) ), M_stiffnessMatrixPtr (
                        new blockMatrix_Type (* (solver.M_stiffnessMatrixPtr) ) ), M_globalMatrixPtr (
                            new blockMatrix_Type (* (solver.M_globalMatrixPtr) ) ), M_initialTime (
                                solver.M_initialTime), M_endTime (solver.M_endTime), M_timeStep (
                                    solver.M_timeStep), M_diffusionTensorIntra (solver.M_diffusionTensorIntra),
    M_diffusionTensorExtra (solver.M_diffusionTensorExtra), M_rhsPtr (
        new blockVector_Type (* (solver.M_rhsPtr) ) ), M_rhsPtrUnique (
            new blockVector_Type (* (M_rhsPtr), Unique) ), M_potentialGlobalPtr (
                new blockVector_Type (solver.M_ETFESpacePtr->map() | solver.M_ETFESpacePtr->map() ) ), M_linearSolverPtr (
                    new LinearSolver (* (solver.M_linearSolverPtr) ) ), M_elementsOrder (
                        solver.M_elementsOrder), M_fiberPtr (
                            new vector_Type (* (solver.M_fiberPtr) )
                        ) ,
    M_lumpedMassMatrix (false)
{
    if (M_commPtr -> MyPID() == 0)
    {
        std::cout << "\n WARNING!!! ETA Bidomain Solver: you are using the copy constructor! This method is outdated.";
        std::cout << "\n WARNING!!! ETA Bidomain Solver: Don't count on it, at this moment. Feel free to update  yourself...";
    }

    setupGlobalSolution (M_ionicModelPtr->Size() );
    copyGlobalSolution (solver.M_globalSolution);
    setupGlobalRhs (M_ionicModelPtr->Size() );
    copyGlobalRhs (solver.M_globalRhs);
}

template<typename Mesh, typename IonicModel>
ElectroETABidomainSolver<Mesh, IonicModel>& ElectroETABidomainSolver < Mesh,
                         IonicModel >::operator= (const ElectroETABidomainSolver& solver)
{
    if (M_commPtr -> MyPID() == 0)
    {
        std::cout << "\n WARNING!!! ETA Bidomain Solver: you are using the assignment operator! This method is outdated.";
        std::cout << "\n WARNING!!! ETA Bidomain Solver: Don't count on it, at this moment. Feel free to update  yourself...";
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
    M_diffusionTensorIntra = solver.M_diffusionTensorIntra;
    M_diffusionTensorExtra = solver.M_diffusionTensorExtra;
    setRhs (* (solver.M_rhsPtr) );
    setRhsUnique (* (solver.M_rhsPtrUnique) );
    setPotential (* (solver.M_potentialPtr) );

    setLinearSolver (* (solver.M_linearSolverPtr) );
    copyGlobalSolution (solver.M_globalSolution);
    copyGlobalRhs (solver.M_globalRhs);
    M_elementsOrder = solver.M_elementsOrder;
    setFiber (* (solver.M_fiberPtr) );

    return *this;
}


/********* SETUP METHODS */ //////
template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupFibers()
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
void ElectroETABidomainSolver<Mesh, IonicModel>::setupFibers (
    VectorSmall<3> fibers)
{
    boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > Space3D (
        new FESpace<mesh_Type, MapEpetra> (M_localMeshPtr, M_elementsOrder,
                                           3, M_commPtr) );

    M_fiberPtr.reset (new vector_Type (Space3D->map() ) );

    ElectrophysiologyUtility::setupFibers (*M_fiberPtr, fibers);
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::exportFiberDirection (
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
void ElectroETABidomainSolver<Mesh, IonicModel>::exportFiberDirection (
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
void ElectroETABidomainSolver<Mesh, IonicModel>::importSolution (
    GetPot& dataFile, std::string prefix, std::string postDir, Real time)
{
    const std::string exporterType = dataFile ("exporter/type", "ensight");

    IOFilePtr_Type importer (new hdf5IOFile_Type() );
    importer->setDataFromGetPot (dataFile);
    importer->setPrefix (prefix);
    importer->setPostDir (postDir);

    importer->setMeshProcId (M_feSpacePtr->mesh(),
                             M_feSpacePtr->map().comm().MyPID() );

    int i;
    for (i = 0; i < M_ionicModelPtr->Size(); i++)
        importer->addVariable (IOData_Type::ScalarField,
                               "Variable" + number2string (i), M_feSpacePtr,
                               M_globalSolution.at (i), static_cast<UInt> (0) );
    importer->addVariable (IOData_Type::ScalarField, "Variable" + number2string (i), M_feSpacePtr, M_potentialExtraPtr, static_cast<UInt> (0) );
    importer->importFromTime (time);
    importer->closeFile();

}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setup (GetPot& dataFile,
                                                        short int ionicSize)
{



    M_feSpacePtr.reset (
        new feSpace_Type (M_localMeshPtr, M_elementsOrder, 1, M_commPtr) );
    M_ETFESpacePtr.reset (
        new ETFESpace_Type (M_localMeshPtr, & (M_feSpacePtr -> refFE() ) , M_commPtr) );
    M_displacementETFESpacePtr.reset (
        new ETFESpaceVectorial_Type (M_localMeshPtr, & (M_feSpacePtr -> refFE() ), M_commPtr) );

    M_massMatrixPtr.reset (new blockMatrix_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map()  ) );
    M_stiffnessMatrixPtr.reset (new blockMatrix_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map()  ) );
    M_globalMatrixPtr.reset (new blockMatrix_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map() ) );

    M_rhsPtr.reset (new blockVector_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map(), Repeated) );
    M_rhsPtrUnique.reset (new blockVector_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map(), Unique) );
    M_potentialGlobalPtr.reset (new blockVector_Type (M_ETFESpacePtr->map() | M_ETFESpacePtr->map() ) );
    M_potentialTransPtr.reset (new vector_Type (M_ETFESpacePtr->map() ) );
    M_potentialExtraPtr.reset (new vector_Type (M_ETFESpacePtr->map() ) );


    //***********************//
    //  Setup Linear Solver  //
    //***********************//
    setupLinearSolver (dataFile);

    //**************************//
    //  Setup Initial condition //
    //**************************//
    initializePotential();
    vector_Type Iapp (M_feSpacePtr->map() );
    Iapp *= 0.0;
    M_ionicModelPtr->setAppliedCurrent (Iapp);

    setupGlobalSolution (ionicSize);
    setupGlobalRhs (ionicSize);
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setup (std::string meshName,
                                                        std::string meshPath, GetPot& dataFile, short int ionicSize)
{
    //partitioning the mesh
    partitionMesh (meshName, meshPath);
    setup (dataFile, ionicSize);
}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupMassMatrix()
{


    if (M_lumpedMassMatrix)
    {
        //      if (M_displacementPtr)
        //          setupLumpedMassMatrix(*M_displacementPtr);
        //      else
        setupLumpedMassMatrix();
    }
    else
    {
        if (M_displacementPtr)
        {
            setupMassMatrix (*M_displacementPtr);
        }
        else
        {
            *M_massMatrixPtr *= 0.0;
            if (M_localMeshPtr->comm()->MyPID() == 0)
            {
                std::cout << "\nETA Bidomain Solver: Setting up mass matrix";
            }

            {
                using namespace ExpressionAssembly;
                integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(),
                           M_ETFESpacePtr, M_ETFESpacePtr, phi_i * phi_j)
                        >> M_massMatrixPtr->block (0, 0);
            }

            M_massMatrixPtr->globalAssemble();
        }
    }
}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupMassMatrix (
    vector_Type& disp)
{
    if (M_commPtr->MyPID() == 0)
    {
        std::cout << "\nETA Bidomain Solver: Setting up mass matrix with coupling with mechanics ";
    }

    *M_massMatrixPtr *= 0.0;
    ETFESpaceVectorialPtr_Type spaceVectorial ( new ETFESpaceVectorial_Type ( M_localMeshPtr, & (M_feSpacePtr -> refFE() ), M_commPtr ) );

    {
        using namespace ExpressionAssembly;
        BOOST_AUTO_TPL (I, value (M_identity) );
        BOOST_AUTO_TPL (Grad_u, grad (spaceVectorial, disp) );
        BOOST_AUTO_TPL (F, (Grad_u + I) );
        BOOST_AUTO_TPL (J, det (F) );
        integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                   M_ETFESpacePtr,  J * phi_i * phi_j)
                >> M_massMatrixPtr->block (0, 0);
    }

    M_massMatrixPtr->globalAssemble();
}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupLumpedMassMatrix()
{

    M_lumpedMassMatrix = true;
    if (M_displacementPtr)
    {
        setupLumpedMassMatrix (*M_displacementPtr);
    }
    else
    {
        *M_massMatrixPtr *= 0.0;
        if (M_localMeshPtr->comm()->MyPID() == 0)
        {
            std::cout << "\nETA Bidomain Solver: Setting up lumped mass matrix";
        }
        {
            using namespace ExpressionAssembly;

            integrate (elements (M_localMeshPtr), quadRuleTetra4ptNodal,
                       M_ETFESpacePtr, M_ETFESpacePtr, phi_i * phi_j)
                    >> M_massMatrixPtr->block (0, 0);

        }
        M_massMatrixPtr->globalAssemble();

    }
}



template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupLumpedMassMatrix (
    vector_Type& disp)
{
    /*
        *M_massMatrixPtr *= 0.0;
        boost::shared_ptr<FESpace<mesh_Type, MapEpetra> > vectorialSpace(
                new FESpace<mesh_Type, MapEpetra>(M_localMeshPtr, M_elementsOrder,
                        3, M_commPtr));

        vectorPtr_Type dUdx(new vector_Type(M_displacementPtr->map()));
        vectorPtr_Type dUdy(new vector_Type(M_displacementPtr->map()));
        vectorPtr_Type dUdz(new vector_Type(M_displacementPtr->map()));

        *dUdx = GradientRecovery::ZZGradient(vectorialSpace, *M_displacementPtr, 0);
        *dUdy = GradientRecovery::ZZGradient(vectorialSpace, *M_displacementPtr, 1);
        *dUdz = GradientRecovery::ZZGradient(vectorialSpace, *M_displacementPtr, 2);

        vectorPtr_Type J(new vector_Type(M_potentialPtr->map()));
        int n = J->epetraVector().MyLength();
        int i(0);
        int j(0);
        int k(0);
        for (int p(0); p < n; p++) {
            i = dUdx->blockMap().GID(p);
            j = dUdx->blockMap().GID(p + n);
            k = dUdx->blockMap().GID(p + 2 * n);

            Real F11 = 1.0 + (*dUdx)[i];
            Real F12 =       (*dUdy)[i];
            Real F13 =       (*dUdz)[i];
            Real F21 =       (*dUdx)[j];
            Real F22 = 1.0 + (*dUdy)[j];
            Real F23 =       (*dUdz)[j];
            Real F31 =       (*dUdx)[k];
            Real F32 =       (*dUdy)[k];
            Real F33 = 1.0 + (*dUdz)[k];

            (*J)[i] = F11 * ( F22 * F33 - F32 * F23 )
                    - F12 * ( F21 * F33 - F31 * F23 )
                    + F13 * ( F21 * F32 - F31 * F22 );
        }

        if (M_localMeshPtr->comm()->MyPID() == 0) {
            std::cout << "\nETA Bidomain Solver: Setting up lumped mass matrix coupling with mechanics";
        }
    //CHECK IT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
    /*  {
            using namespace ExpressionAssembly;

            integrate(elements(M_localMeshPtr), quadRuleTetra4ptNodal,
                    M_ETFESpacePtr, M_ETFESpacePtr,
                    value(M_ETFESpacePtr, *J) * phi_i * phi_j) >> M_massMatrixPtr;

        }*/

    //Create one Mass Matrix, because it is the same in all the blocks up to an sign
    /*  boost::shared_ptr<blockMatrix_Type> ETcomponentMassMatrix
                (new blockMatrix_Type ( M_ETFESpacePtr->map() ) );
        *ETcomponentMassMatrix *= 0.0;

        {
            using namespace ExpressionAssembly;

            integrate ( elements (M_ETFESpacePtr->mesh() ),
                           quadRuleTetra4ptNodal,
                        ETFESpacePtr,
                           ETFESpacePtr,
                            value(M_ETFESpacePtr, *J) * phi_i * phi_j
                           )
                        >> ETcomponentMassMatrix->block (0, 0);
        }
        ETcomponentMassMatrix->globalAssemble();


        // ---------------------------------------------------------------
        // We copy the blocks
        // ---------------------------------------------------------------
        //result matrix
        // M   -M
        // -M     M
        MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMassMatrix->block (0, 0), *M_massMatrixPtr->block (0, 0) );
        MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMassMatrix->block (0, 0), *M_massMatrixPtr->block (1, 1) );
        *ETsystemMatrixIII *= -1.0;
        MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMassMatrix->block (0, 0), *M_massMatrixPtr->block (0, 1) );
        MatrixEpetraStructuredUtility::copyBlock (*ETcomponentMassMatrix->block (0, 0), *M_massMatrixPtr->block (1, 0) );


        M_massMatrixPtr->globalAssemble();
    */
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupStiffnessMatrix()
{

    if (M_displacementPtr)
    {
        setupStiffnessMatrix (M_displacementPtr);
    }
    else
    {
        setupStiffnessMatrix (M_diffusionTensorIntra, M_diffusionTensorExtra);
    }


}



template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupStiffnessMatrix (
    VectorSmall<3> diffusionIntra, VectorSmall<3> diffusionExtra)
{

    if (M_localMeshPtr->comm()->MyPID() == 0)
    {
        std::cout
                << "\nETA Bidomain Solver: Setting up stiffness matrix (only fiber field)";
    }

    *M_stiffnessMatrixPtr *= 0.0;
    Real sigmal = diffusionIntra[0];
    Real sigmat = diffusionIntra[1];

    ETFESpaceVectorialPtr_Type spaceVectorial (
        new ETFESpaceVectorial_Type (M_localMeshPtr, & (M_feSpacePtr -> refFE() ), M_commPtr) );
    blockMatrixPtr_Type ETcomponentStiffnessMatrix
    (new blockMatrix_Type ( M_ETFESpacePtr->map() ) );
    *ETcomponentStiffnessMatrix *= 0.0;

    //Intracellular conduction
    {
        using namespace ExpressionAssembly;

        BOOST_AUTO_TPL (I, value (M_identity) );
        BOOST_AUTO_TPL (f0, value (spaceVectorial, *M_fiberPtr) );
        BOOST_AUTO_TPL (D,
                        value (sigmat) * I
                        + (value (sigmal) - value (sigmat) )
                        * outerProduct (f0, f0) );

        integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                   M_ETFESpacePtr,
                   dot ( D * grad (phi_i), grad (phi_j) ) )
                >> ETcomponentStiffnessMatrix->block (0, 0);

    }
    ETcomponentStiffnessMatrix->globalAssemble();

    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentStiffnessMatrix->block (0, 0), *M_stiffnessMatrixPtr->block (0, 1) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentStiffnessMatrix->block (0, 0), *M_stiffnessMatrixPtr->block (1, 0) );
    MatrixEpetraStructuredUtility::copyBlock (*ETcomponentStiffnessMatrix->block (0, 0), *M_stiffnessMatrixPtr->block (0, 0) );

    // Extracellular conduction
    sigmal += diffusionExtra[0];
    sigmat += diffusionExtra[1];

    {
        using namespace ExpressionAssembly;

        BOOST_AUTO_TPL (I, value (M_identity) );
        BOOST_AUTO_TPL (f0, value (spaceVectorial, *M_fiberPtr) );
        BOOST_AUTO_TPL (D,
                        (value (sigmat) * I
                         + (value (sigmal) - value (sigmat) )
                         * outerProduct (f0, f0) ) );

        integrate (elements (M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                   M_ETFESpacePtr,
                   dot ( D * grad (phi_i), grad (phi_j) ) )
                >> M_stiffnessMatrixPtr->block (1, 1);

    }
    M_stiffnessMatrixPtr->globalAssemble();
}





template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupStiffnessMatrix (
    vectorPtr_Type disp)
{
    /*
      //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
      //At least the ExtraCellular is wrong. It whould be minus.
        if (M_localMeshPtr->comm()->MyPID() == 0) {
            std::cout
                    << "\nETA Bidomain Solver: Setting up stiffness matrix  coupling with mechanics";
        }

        *M_stiffnessMatrixPtr *= 0.0;

    //Intra Cellular
        Real sigmal = M_diffusionTensorIntra[0];
        Real sigmat = M_diffusionTensorIntra[1];

        ETFESpaceVectorialPtr_Type spaceVectorial(
                new ETFESpaceVectorial_Type(M_localMeshPtr, &(M_feSpacePtr -> refFE()), M_commPtr));

        {
            using namespace ExpressionAssembly;

            BOOST_AUTO_TPL(I, value(M_identity));
            BOOST_AUTO_TPL(Grad_u, grad(spaceVectorial, *disp));
            BOOST_AUTO_TPL(F, (Grad_u + I));
            BOOST_AUTO_TPL(FmT, minusT(F));
            BOOST_AUTO_TPL(Fm1, transpose(FmT));
            BOOST_AUTO_TPL(J, det(F));
            BOOST_AUTO_TPL(Jm23, pow(J, -2. / 3));
            BOOST_AUTO_TPL(f0, value(spaceVectorial, *M_fiberPtr));
            BOOST_AUTO_TPL(D,
                    value(sigmat) * I
                            + (value(sigmal) - value(sigmat))
                                    * outerProduct(f0, f0));

            integrate(elements(M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                    M_ETFESpacePtr,
                    dot(J * Fm1 * D * FmT * grad(phi_i), grad(phi_j)))
                    >> M_stiffnessMatrixPtr->block (0, 0);

        }
    //Extra cellular
        sigmal = M_diffusionTensorExtra[0];
        sigmat = M_diffusionTensorExtra[1];
        {
            using namespace ExpressionAssembly;

            BOOST_AUTO_TPL(I, value(M_identity));
            BOOST_AUTO_TPL(Grad_u, grad(spaceVectorial, *disp));
            BOOST_AUTO_TPL(F, (Grad_u + I));
            BOOST_AUTO_TPL(FmT, minusT(F));
            BOOST_AUTO_TPL(Fm1, transpose(FmT));
            BOOST_AUTO_TPL(J, det(F));
            BOOST_AUTO_TPL(Jm23, pow(J, -2. / 3));
            BOOST_AUTO_TPL(f0, value(spaceVectorial, *M_fiberPtr));
            BOOST_AUTO_TPL(D,
                    value(sigmat) * I
                            + (value(sigmal) - value(sigmat))
                                    * outerProduct(f0, f0));

            integrate(elements(M_localMeshPtr), M_feSpacePtr->qr(), M_ETFESpacePtr,
                    M_ETFESpacePtr,
                    dot(J * Fm1 * D * FmT * grad(phi_i), grad(phi_j)))
                    >> M_stiffnessMatrixPtr->block (0, 0);

        }

        M_stiffnessMatrixPtr->globalAssemble();
    */
}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupGlobalMatrix()
{
    (*M_globalMatrixPtr) *= 0;
    (*M_globalMatrixPtr) = (*M_stiffnessMatrixPtr);
    (*M_globalMatrixPtr) *= 1.0 / M_surfaceVolumeRatio;

    // This line should be right, as long as the membraneCpacitance is intra and extra cellular the same
    (*M_globalMatrixPtr) += ( (*M_massMatrixPtr) * ( M_ionicModelPtr -> membraneCapacitance() / M_timeStep) );
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupLinearSolver (
    GetPot dataFile)
{

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot (dataFile, "prec");
    precPtr.reset (precRawPtr);
    Teuchos::RCP < Teuchos::ParameterList > solverParamList = Teuchos::rcp (
                                                                  new Teuchos::ParameterList);

    std::string xmlpath = dataFile ("electrophysiology/bidomain_xml_path",
                                    "./");
    std::string xmlfile = dataFile ("electrophysiology/bidomain_xml_file",
                                    "BidomainSolverParamList.xml");

    solverParamList = Teuchos::getParametersFromXmlFile (xmlpath + xmlfile);

    M_linearSolverPtr->setCommunicator (M_commPtr);
    M_linearSolverPtr->setParameters (*solverParamList);
    M_linearSolverPtr->setPreconditioner (precPtr);
    M_linearSolverPtr->setOperator (M_globalMatrixPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupLinearSolver (
    GetPot dataFile, list_Type list)
{

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot (dataFile, "prec");
    precPtr.reset (precRawPtr);

    M_linearSolverPtr->setCommunicator (M_commPtr);
    M_linearSolverPtr->setParameters (list);
    M_linearSolverPtr->setPreconditioner (precPtr);
    M_linearSolverPtr->setOperator (M_globalMatrixPtr);
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupGlobalSolution (
    short int ionicSize)
{
    //This is the solution vector for the IonicModel
    M_globalSolution.push_back (M_potentialTransPtr);
    //M_globalSolution.push_back(M_potentialExtraPtr);
    for (int k = 1; k < ionicSize; ++k)
    {
        M_globalSolution.push_back (
            * (new vectorPtr_Type (new VectorEpetra (M_ETFESpacePtr->map() ) ) ) );
    }
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupGlobalRhs (
    short int ionicSize)
{
    M_globalRhs.push_back (M_ionicModelPtr->appliedCurrentPtr() );
    for (int k = 1; k < ionicSize; ++k)
    {
        M_globalRhs.push_back (
            * (new vectorPtr_Type (new VectorEpetra (M_ETFESpacePtr->map() ) ) ) );
    }
}

/****************Solver ********************************************/

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::solveOneReactionStepFE (
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
    //copying th transmembrane potential back
    M_potentialGlobalPtr->subset ( (*M_globalSolution.at (0) ), (M_globalSolution.at (0)->map() ), (UInt) 0, (UInt) 0);
}


template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::solveOneReactionStepRL (
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
    M_potentialGlobalPtr->subset ( (*M_globalRhs.at (0) ), (M_globalRhs.at (0)->map() ), (UInt) 0, (UInt) 0);
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::solveOneDiffusionStepBE()
{
    if (M_displacementPtr)
    {
        M_linearSolverPtr->setOperator (M_globalMatrixPtr);
    }
    M_linearSolverPtr->setRightHandSide (M_rhsPtrUnique);
    M_linearSolverPtr->solve (M_potentialGlobalPtr);

    M_potentialTransPtr->subset ( (const VectorEpetra&) (*M_potentialGlobalPtr->block (0)->vectorPtr() ), M_potentialTransPtr->map(), M_potentialGlobalPtr->block (0)->firstIndex(), static_cast<UInt> (0) );
    M_potentialExtraPtr->subset ( (const VectorEpetra&) (*M_potentialGlobalPtr->block (1)->vectorPtr() ), M_potentialExtraPtr->map(), M_potentialGlobalPtr->block (1)->firstIndex(), static_cast<UInt> (0) );
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::solveOneSplittingStep()
{
    solveOneReactionStepFE();
    (*M_rhsPtrUnique) *= 0.0;
    updateRhs();
    solveOneDiffusionStepBE();
}



/************** EXPORTER */    //////////////
template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupPotentialExporter (
    IOFile_Type& exporter)
{
    exporter.setMeshProcId (M_localMeshPtr, M_commPtr->MyPID() );
    exporter.setPrefix ("Potential");
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "PotentialTransmembrane",
                          M_feSpacePtr, M_globalSolution.at (0), UInt (0) );
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "PotentialExtracellular",
                          M_feSpacePtr, M_potentialExtraPtr, UInt (0) );

}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupPotentialExporter (
    IOFile_Type& exporter, std::string fileName)
{
    exporter.setMeshProcId (M_localMeshPtr, M_commPtr->MyPID() );
    exporter.setPrefix (fileName);
    exporter.exportPID (M_localMeshPtr, M_commPtr);
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "PotentialTransmembrane",
                          M_feSpacePtr, M_globalSolution.at (0), UInt (0) );
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "PotentialExtracellular",
                          M_feSpacePtr, M_potentialExtraPtr, UInt (0) );
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setupExporter (
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
    exporter.addVariable (ExporterData<mesh_Type>::ScalarField, "PotentialExtracellular",
                          M_feSpacePtr, M_potentialExtraPtr, UInt (0) );
}

/********   INITIALIZITION FOR CONSTRUCTOR ****/    //////
template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init()
{
    M_linearSolverPtr.reset (new LinearSolver() );
    M_globalSolution = * (new vectorOfPtr_Type() );
    M_globalRhs = * (new vectorOfPtr_Type() );
    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init (
    ionicModelPtr_Type model)
{
    init();
    M_commPtr.reset (new Epetra_MpiComm (MPI_COMM_WORLD) );
    M_localMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_ionicModelPtr = model;

}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init (commPtr_Type comm)
{
    init();
    M_commPtr = comm;
    M_localMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init (meshPtr_Type meshPtr)
{
    init();
    //TODO change the meshPtr to pass the fullMeshPtr
    M_localMeshPtr = meshPtr;
    M_fullMeshPtr.reset (new mesh_Type (M_commPtr) );
    M_commPtr = meshPtr->comm();
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init (commPtr_Type comm,
                                                       ionicModelPtr_Type model)
{
    init (comm);
    M_ionicModelPtr = model;
}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::init (meshPtr_Type meshPtr,
                                                       ionicModelPtr_Type model)
{
    init (meshPtr);
    M_ionicModelPtr = model;
}

/********* parameter initialization */    ////////
template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setParameters()
{

    // Suggested values from Niederer et al. 2011

    M_surfaceVolumeRatio = 1400.0;

    M_diffusionTensorIntra[0] = 1.7;
    M_diffusionTensorIntra[1] = 0.19;
    M_diffusionTensorIntra[2] = 0.19;

    M_diffusionTensorExtra[0] = 6.2;
    M_diffusionTensorExtra[1] = 2.4;
    M_diffusionTensorExtra[2] = 2.4;

    M_initialTime = 0.0;
    M_endTime = 50.0;
    M_timeStep = 0.01;
    M_elementsOrder = "P1";
    M_lumpedMassMatrix = false;

}

template<typename Mesh, typename IonicModel>
void ElectroETABidomainSolver<Mesh, IonicModel>::setParameters (
    list_Type list)
{

    // Suggested values from Niederer et al. 2011

    M_surfaceVolumeRatio = list.get ("surfaceVolumeRatio", 1400.0);

    M_diffusionTensorIntra[0] = list.get ("longitudinalDiffusionIntra", 1.7);
    M_diffusionTensorIntra[1] = list.get ("transversalDiffusionIntra", 0.19);
    M_diffusionTensorIntra[2] = M_diffusionTensorIntra[1];

    M_diffusionTensorExtra[0] = list.get ("longitudinalDiffusionExtra", 6.2);
    M_diffusionTensorExtra[1] = list.get ("transversalDiffusionExtra", 2.4);
    M_diffusionTensorExtra[2] = M_diffusionTensorExtra[1];

    M_initialTime = list.get ("initialTime", 0.0);
    M_endTime = list.get ("endTime", 50.0);
    M_timeStep = list.get ("timeStep", 0.01);
    M_elementsOrder = list.get ("elementsOrder", "P1");
    M_lumpedMassMatrix = list.get ("LumpedMass", false);

}

} // namespace LifeV

#endif //_BIDOMAINSOLVER_H_






