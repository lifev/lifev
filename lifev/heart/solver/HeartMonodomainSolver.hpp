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
  @brief Class for solving the Monodomain equations in electrophysiology.

  @date 11-2007
  @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>

  @contributors J.Castelneau (INRIA), Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @mantainer Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
  @last update 11-2010
 */

#ifndef _MONODOMAINSOLVER_H_
#define _MONODOMAINSOLVER_H_

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/fem/SobolevNorms.hpp>
#include <lifev/core/fem/GeometricMap.hpp>
#include <lifev/heart/solver/HeartMonodomainData.hpp>
#include <lifev/core/util/LifeChrono.hpp>
#include <boost/shared_ptr.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/heart/solver/HeartStiffnessFibers.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

namespace LifeV
{

//! monodomainSolver - Class featuring the usual solver for monodomain equations

template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >
class HeartMonodomainSolver
{

public:

    //! @name Type definitions
    //@{

    typedef HeartMonodomainData data_type;
    typedef Real ( *Function ) ( const Real& t,
                                 const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const ID& id );
    typedef boost::function < Real ( Real const&, Real const&, Real const&,
                                     Real const&, ID const& ) > source_Type;

    typedef Mesh mesh_Type;
    typedef BCHandler                               bcHandlerRaw_Type;
    typedef boost::shared_ptr<bcHandlerRaw_Type>    bcHandler_Type;

    typedef typename SolverType::matrix_type        matrix_Type;
    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;
    typedef typename SolverType::vector_type        vector_Type;

    typedef typename SolverType::prec_raw_type      preconditionerRaw_Type;
    typedef typename SolverType::prec_type          preconditioner_Type;


    //@}



    //! @name Constructors & Destructor
    //@{

    //! Constructor
    /*!
     * @param dataType
     * @param potential FE space
     * @param bcHandler boundary conditions for potential
     * @param Epetra communicator
     */

    HeartMonodomainSolver ( const data_type& dataType,
                            FESpace<Mesh, MapEpetra>& uFESpace,
                            BCHandler& bcHandler,
                            boost::shared_ptr<Epetra_Comm>& comm);

    //! Destructor
    virtual ~HeartMonodomainSolver() {}

    //@}

    //! @name Methods
    //@{

    //! Updates sources, bc treatment and solve the monodomain system
    virtual void PDEiterate ( bcHandlerRaw_Type& bch );

    //! Sets up the system solver
    virtual void setup        ( const GetPot& dataFile );

    //! Builds time independent parts of PDE system
    virtual void buildSystem();

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem (Real alpha, vector_Type&  sourceVec);

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem ( vector_Type& sourceVec );

    //! Initialize
    void initialize ( const source_Type& );
    void initialize ( const Function&  );
    void initialize ( const vector_Type& );

    //! Returns the local solution vector
    const vector_Type& solutionTransmembranePotential() const
    {
        return M_solutionTransmembranePotential;
    }

    const vector_Type& fiberVector() const
    {
        return M_fiberVector;
    }

    //! Returns the local residual vector
    const vector_Type& residual() const
    {
        return M_residual;
    }

    //! Returns u FE space
    FESpace<Mesh, MapEpetra>& potentialFESpace()
    {
        return M_uFESpace;
    }

    //! Setting of the boundary conditions
    void setBC ( BCHandler& BCh_u )
    {
        M_BChandlerElectric = &BCh_u;
        M_setBC = true;
    }

    //! Postprocessing
    void postProcessing (bool _writeMesh = false);

    void resetPreconditioner()
    {
        M_resetPreconditioner = true;
    }

    //! Return maps
    Epetra_Map const& getRepeatedMapEpetra() const
    {
        return *M_localMap.map (Repeated);
    }

    Epetra_Map const& getRepeatedMapEpetraVec() const
    {
        return *M_localMapVector.map (Repeated);
    }

    MapEpetra const& getMap() const
    {
        return M_localMap;
    }

    void recomputeMatrix (bool const recomp)
    {
        M_recomputeMatrix = recomp;
    }

    matrix_Type& massMatrix()
    {
        return *M_massMatrix;
    }

    //@}
protected:

    //! Solves PDE system
    void solveSystem (  matrixPtr_Type matrFull, vector_Type&   rhsFull );

    //! Apply BC
    void applyBoundaryConditions ( matrix_Type& matrix, vector_Type& rhs, bcHandlerRaw_Type& BCh );

    //! Data
    const data_type&               M_data;

    //! u FE space
    FESpace<Mesh, MapEpetra>&      M_uFESpace;

    //! MPI communicator
    //Epetra_Comm*                   M_comm;
    const boost::shared_ptr<Epetra_Comm> M_comm;
    Int                            M_me;

    //! Monodomain BC
    BCHandler*                     M_BChandlerElectric;
    bool                           M_setBC;

    //! Map
    MapEpetra                      M_localMap;
    MapEpetra                      M_localMapVector;

    //! mass matrix
    matrixPtr_Type                 M_massMatrix;

    //! Stiff matrix: D*stiff
    matrixPtr_Type                 M_stiffnessMatrix;

    matrixPtr_Type                 M_matrNoBC;

    //! Right hand side for the PDE
    vector_Type                    M_rhsNoBC;

    //! Global solution _u
    vector_Type                    M_solutionTransmembranePotential;

    //! Local fibers vector
    vector_Type                    M_fiberVector;

    //! residual
    vector_Type                    M_residual;

    //! Solver
    SolverType                     M_linearSolver;

    preconditioner_Type            M_preconditioner;

    Real                         M_diagonalize;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //! Boolean that indicates if the matrix is updated for the current iteration
    bool                           M_updated;

    //! Boolean that indicates if the precond has to be recomputed
    bool                           M_reusePreconditioner;
    bool                           M_resetPreconditioner;

    //! Integer storing the max number of solver iteration with prec recomputing
    Int                            M_maxIteration;

    //! Boolean that indicates if the matrix has to be recomputed
    bool                           M_recomputeMatrix;

private:

    //! Elementary matrices
    MatrixElemental                        M_stiffnessElementaryMatrix;
    MatrixElemental                        M_massElementaryMatrix;
    Real                           massCoefficient;
    UInt dim_u() const
    {
        return M_uFESpace.dim();
    }

}; // class MonodomainSolver



//
// IMPLEMENTATION
//

//! Constructors
template<typename Mesh, typename SolverType>
HeartMonodomainSolver<Mesh, SolverType>::
HeartMonodomainSolver ( const data_type&          dataType,
                        FESpace<Mesh, MapEpetra>& uFESpace,
                        BCHandler&                BCh_u,
                        boost::shared_ptr<Epetra_Comm>&              comm ) :
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_comm                   ( comm ),
    M_me                     ( M_comm->MyPID() ),
    M_BChandlerElectric      ( &BCh_u ),
    M_setBC                  ( true ),
    M_localMap               ( M_uFESpace.map() ),
    M_localMapVector         (M_localMap + M_localMap + M_localMap),
    M_massMatrix             ( ),
    M_stiffnessMatrix        ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_solutionTransmembranePotential      ( M_localMap ),
    M_fiberVector                  ( M_localMapVector, Repeated ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_preconditioner         ( ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePreconditioner    ( true ),
    M_resetPreconditioner    ( true ),
    M_maxIteration           ( -1 ),
    M_recomputeMatrix        ( false ),
    M_stiffnessElementaryMatrix ( M_uFESpace.fe().nbFEDof(), 1, 1 ),
    M_massElementaryMatrix   ( M_uFESpace.fe().nbFEDof(), 1, 1 )
{

    if (M_data.hasFibers() )
    {
        std::stringstream MyPID;
        ifstream fibers (M_data.fibersFile().c_str() );

        std::cout << "fiber_file: " <<  M_data.fibersFile().c_str() << std::endl;
        UInt NumGlobalElements = M_localMapVector.map (Repeated)->NumGlobalElements();
        std::vector<Real> fiber_global_vector (NumGlobalElements);

        for ( UInt i = 0; i < NumGlobalElements; ++i)
        {
            fibers >> fiber_global_vector[i];
        }
        UInt NumMyElements = M_localMapVector.map (Repeated)->NumMyElements();
        for (UInt j = 0; j < NumMyElements; ++j)
        {
            UInt ig = M_localMapVector.map (Repeated)->MyGlobalElements() [j];
            M_fiberVector[ig] = fiber_global_vector[ig];
        }
        std::cout << std::endl;
        fiber_global_vector.clear();
    }
}


template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::setup ( const GetPot& dataFile )
{

    M_diagonalize = dataFile ( "electric/space_discretization/diagonalize",  1. );

    M_reusePreconditioner   = dataFile ( "electric/prec/reuse", true);

    M_linearSolver.setCommunicator (M_comm);

    M_linearSolver.setDataFromGetPot ( dataFile, "electric/solver" );

    M_maxIteration = dataFile ( "electric/solver/max_iter", -1);

    std::string precType = dataFile ( "electric/prec/prectype", "Ifpack");

    M_preconditioner.reset ( PRECFactory::instance().createObject ( precType ) );
    ASSERT (M_preconditioner.get() != 0, "monodomainSolver : Preconditioner not set");

    M_preconditioner->setDataFromGetPot ( dataFile, "electric/prec" );
}


template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::buildSystem()
{
    M_massMatrix.reset  ( new matrix_Type (M_localMap) );
    M_stiffnessMatrix.reset ( new matrix_Type (M_localMap) );

    if (M_verbose)
    {
        std::cout << "  f-  Computing constant matrices ...        ";
    }

    LifeChrono chrono;

    LifeChrono chronoDer;
    LifeChrono chronoStiff;
    LifeChrono chronoMass;


    LifeChrono chronoStiffAssemble;
    LifeChrono chronoMassAssemble;
    LifeChrono chronoZero;

    M_comm->Barrier();

    chrono.start();

    //! Elementary computation and matrix assembling
    //! Loop on elements

    for ( UInt iVol = 0; iVol < M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoZero.start();

        M_stiffnessElementaryMatrix.zero();
        M_massElementaryMatrix.zero();

        chronoZero.stop();

        chronoStiff.start();
        switch (M_data.heartDiffusionFactor() )
        {
            case 0:
                M_uFESpace.fe().updateFirstDeriv ( M_uFESpace.mesh()->volumeList ( iVol ) );
                if (M_data.hasFibers() )
                {
                    stiff ( M_data.longitudinalConductivity(),
                            M_data.transversalConductivity(),
                            M_fiberVector,
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                else
                {
                    AssemblyElemental::stiff ( M_data.diffusivity(), M_stiffnessElementaryMatrix,  M_uFESpace.fe(), 0, 0 );
                }
                break;

            case 1:
                M_uFESpace.fe().updateFirstDerivQuadPt ( M_uFESpace.mesh()->volumeList ( iVol ) );
                if (M_data.hasFibers() )
                {
                    stiff ( M_data.reducedConductivitySphere(),
                            M_data.longitudinalConductivity(),
                            M_data.transversalConductivity(),
                            M_fiberVector,
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(), 0, 0);
                }
                else
                {
                    stiff ( M_data.reducedConductivitySphere(),
                            M_data.diffusivity(), M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                break;
            case 2:
                M_uFESpace.fe().updateFirstDerivQuadPt ( M_uFESpace.mesh()->volumeList ( iVol ) );
                if (M_data.hasFibers() )
                {
                    stiff ( M_data.reducedConductivityCylinder(),
                            M_data.longitudinalConductivity(),
                            M_data.transversalConductivity(),
                            M_fiberVector,
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                else
                {
                    stiff ( M_data.reducedConductivityCylinder(),
                            M_data.diffusivity(),
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                break;
            case 3:
                M_uFESpace.fe().updateFirstDerivQuadPt ( M_uFESpace.mesh()->volumeList ( iVol ) );
                if (M_data.hasFibers() )
                {
                    stiff ( M_data.reducedConductivityBox(),
                            M_data.longitudinalConductivity(),
                            M_data.transversalConductivity(),
                            M_fiberVector,
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                else
                {
                    stiff ( M_data.reducedConductivityBox(),
                            M_data.diffusivity(),
                            M_stiffnessElementaryMatrix,
                            M_uFESpace.fe(),
                            M_uFESpace.dof(),
                            0,
                            0);
                }
                break;
        }
        chronoStiff.stop();

        chronoMass.start();
        AssemblyElemental::mass ( 1., M_massElementaryMatrix, M_uFESpace.fe(), 0, 0 );
        chronoMass.stop();


        chronoStiffAssemble.start();
        assembleMatrix ( *M_stiffnessMatrix,
                         M_stiffnessElementaryMatrix,
                         M_uFESpace.fe(),
                         M_uFESpace.fe(),
                         M_uFESpace.dof(),
                         M_uFESpace.dof(),
                         0, 0,
                         0, 0);
        chronoStiffAssemble.stop();


        chronoMassAssemble.start();
        assembleMatrix ( *M_massMatrix,
                         M_massElementaryMatrix,
                         M_uFESpace.fe(),
                         M_uFESpace.fe(),
                         M_uFESpace.dof(),
                         M_uFESpace.dof(),
                         0, 0,
                         0, 0);
        chronoMassAssemble.stop();

    }

    massCoefficient = M_data.volumeSurfaceRatio() * M_data.membraneCapacitance() / M_data.timeStep();

    M_comm->Barrier();

    chrono.stop();
    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    }

    if (M_verbose)
    {
        std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    }
    chrono.start();

    M_stiffnessMatrix->globalAssemble();
    M_massMatrix->globalAssemble();

    M_comm->Barrier();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_stiffnessMatrix->meanNumEntries() ) );

    //! Computing 1/dt * M + K

    *M_matrNoBC += *M_stiffnessMatrix;

    *M_matrNoBC += *M_massMatrix * massCoefficient;

    M_matrNoBC->globalAssemble();

    chrono.stop();
    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s." << std::endl;
    }

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diffCumul() << " s.\n"
                  << " Zero           " << chronoZero.diffCumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
                  << " Mass           " << chronoMass.diffCumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
                  << std::endl;

}

template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::
initialize ( const source_Type& u0 )
{

    vector_Type u (M_uFESpace.map() );

    M_uFESpace.interpolate (u0, u, 0.);

    initialize (u);
}



template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::
initialize ( const Function& u0 )
{

    vector_Type u (M_uFESpace.map() );
    M_uFESpace.interpolate (u0, u, 0.);

    initialize (u);
}


template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::
initialize ( const vector_Type& u0 )
{

    M_solutionTransmembranePotential = u0;

}


template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::
updatePDESystem (Real         alpha,
                 vector_Type& sourceVec )
{

    LifeChrono chrono;

    if (M_verbose)
        std::cout << "  f-  Updating mass term and right hand side... "
                  << std::flush;

    chrono.start();

    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();

    chrono.stop();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
    }

    M_updated = false;

    if (M_recomputeMatrix)
    {
        buildSystem();
    }

    if (M_verbose)
        std::cout << "  f-  Copying the matrices ...                 "
                  << std::flush;

    chrono.start();

    M_matrNoBC.reset (new matrix_Type (M_localMap, M_stiffnessMatrix->meanNumEntries() ) );

    *M_matrNoBC += *M_stiffnessMatrix;

    *M_matrNoBC += *M_massMatrix * alpha;


    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                                 << std::flush;



    M_updated = true;
    M_matrNoBC->globalAssemble();

}

template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::
updatePDESystem (vector_Type& sourceVec )
{

    LifeChrono chrono;

    if (M_verbose)
        std::cout << "  f-  Updating right hand side... "
                  << std::flush;

    chrono.start();

    M_rhsNoBC = sourceVec;
    M_rhsNoBC.globalAssemble();

    chrono.stop();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;
    }

}



template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::PDEiterate ( bcHandlerRaw_Type& bch )
{

    LifeChrono chrono;
    chrono.start();

    matrixPtr_Type matrFull ( new matrix_Type (*M_matrNoBC) );
    vector_Type    rhsFull = M_rhsNoBC;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying boundary conditions ...         "
                                 << std::flush;

    chrono.start();

    applyBoundaryConditions ( *matrFull, rhsFull, bch);

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose)
    {
        std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    }

    //! Solving the system
    solveSystem ( matrFull, rhsFull );

    //  M_residual  = M_rhsNoBC;
    //  M_residual -= *M_matrNoBC*M_solutionTransmembranePotential;

} // iterate()



template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::solveSystem ( matrixPtr_Type  matrFull,
                                                            vector_Type&    rhsFull )
{
    LifeChrono chrono;

    if (M_verbose)
    {
        std::cout << "  f-  Setting up the solver ...                ";
    }

    chrono.start();
    M_linearSolver.setMatrix (*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    if ( !M_reusePreconditioner || M_resetPreconditioner )
    {
        chrono.start();

        if (M_verbose)
        {
            std::cout << "  f-  Computing the precond ...                ";
        }

        M_preconditioner->buildPreconditioner (matrFull);

        Real condest = M_preconditioner->condest();

        M_linearSolver.setPreconditioner (M_preconditioner);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPreconditioner = false;
    }
    else
    {
        if (M_verbose)
        {
            std::cout << "f-  Reusing  precond ...                \n" <<  std::flush;
        }
    }


    chrono.start();

    if (M_verbose)
    {
        std::cout << "f-  Solving system ...                                ";
    }

    Int numIter = M_linearSolver.solve (M_solutionTransmembranePotential, rhsFull);

    chrono.stop();

    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }


    if (numIter > M_maxIteration)
    {
        M_resetPreconditioner = true;
    }

    M_comm->Barrier();

}



template<typename Mesh, typename SolverType>
void HeartMonodomainSolver<Mesh, SolverType>::applyBoundaryConditions ( matrix_Type&        matrix,
                                                                        vector_Type&        rhs,
                                                                        bcHandlerRaw_Type&  BCh )
{

    // BC manage for the PDE
    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate ( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    vector_Type rhsFull (M_rhsNoBC, Repeated, Zero);

    bcManage ( matrix, rhs, *M_uFESpace.mesh(), M_uFESpace.dof(),
               BCh, M_uFESpace.feBd(), 1., M_data.time() );

    rhs = rhsFull;
    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize ( 1 * dim_u(),
                             M_diagonalize,
                             rhs,
                             0.);
    }

} // applyBoundaryCondition


} // namespace LifeV


#endif //_MONODOMAINSOLVER_H_
