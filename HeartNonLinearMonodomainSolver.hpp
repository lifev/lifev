//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with LifeV. If not, see <http://www.gnu.org/licences/>.

************************************************************************
*/
//@HEADER
/*!
 * @file
 * @brief Solver class for the monodomain equations including semilinear diffusion
 * @author Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 * @date 11-2009
 * @mantainer Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 */

#ifndef _HEARTNONLINEARMONODOMAINSOLVER_H_
#define _HEARTNONLINEARMONODOMAINSOLVER_H_

#include <fstream>
#include <boost/shared_ptr.hpp>

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/SobolevNorms.hpp>
#include <life/lifefem/GeometricMap.hpp>
#include <life/lifesolver/HeartMonodomainData.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifesolver/HeartStiffnessFibers.hpp>
#include <life/lifefem/bdf_template.hpp>

namespace LifeV
{

//!class HeartNonLinearMonodomainSolver - implements a nonlinear monodomain solver.

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class HeartNonLinearMonodomainSolver
{

public:

   //! @name Type definitions
    //@{

    typedef HeartMonodomainData data_type;

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

    typedef typename SolverType::matrix_type      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>        matrixPtr_Type;
    typedef typename SolverType::vector_type      vector_Type;

    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;


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

    HeartNonLinearMonodomainSolver( const data_type&          dataType,
                         FESpace<Mesh, EpetraMap>& uFESpace,
                         BCHandler&                bcHandler,
                         boost::shared_ptr<Epetra_Comm> comm );


    //! virtual destructor
    virtual ~HeartNonLinearMonodomainSolver();

    //@}

    //! @name Methods
    //@{

    //! Updates sources, bc treatment and solve the monodomain system
    virtual void PDEiterate( bchandler_raw_type& bch );

    //! Sets up the system solver
    virtual void setUp        ( const GetPot& dataFile );

    //! Builds time independent parts of PDE system
    virtual void buildSystem();

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem(Real       alpha,
                              vector_Type& sourceVec
                              );

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem( vector_Type& sourceVec );

    //! Initialize
    void initialize( const source_type& );
    void initialize( const Function&  );
    void initialize( const vector_Type& );

    //! Returns the local solution vector
    const vector_Type& solution_u() const {return M_sol_u;}

     const vector_Type& fiber_vector() const {return M_fiber_vector;}

    //! Returns the local displacements vector
    vector_Type& disp()        { return M_disp; }

    //! Returns the local residual vector
    const vector_Type& residual() const {return M_residual;}

    //! Moves the mesh
    void moveMesh(vector_Type const &dep);



    //! Returns u FE space
    FESpace<Mesh, EpetraMap>& potentialFESpace() {return M_uFESpace;}

    //! Boundary Conditions
    const bool BCset() const {return M_setBC;}

    //! Sets Monodomain BCs
    void setBC(BCHandler &BCh_u)
        {
            M_BCh_electric = &BCh_u; M_setBC = true;
        }

    //! Postprocessing
    void postProcess(bool _writeMesh = false);

    void resetPrec() {M_resetPrec = true;}

    //! Return maps
    Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getMap(Repeated); }

    Epetra_Map const& getRepeatedEpetraMapVec() const { return *M_localMapVec.getMap(Repeated); }

    EpetraMap const& getMap() const { return M_localMap; }

    void recomputeMatrix(bool const recomp){M_recomputeMatrix = recomp;}

    matrix_Type& matrMass()
        {
            return *M_matrMass;
        }
    //@}

protected:

	//! Solves PDE system
    void solveSystem            (  matrixPtr_Type matrFull,
                                   vector_Type&   rhsFull );

    //! Apply BC
    void applyBoundaryConditions(  matrix_Type&        matrix,
                                   vector_Type&        rhs,
                                   bchandler_raw_type& BCh);

    //! Data
    const data_type&               M_data;

    //! u FE space
    FESpace<Mesh, EpetraMap>&      M_uFESpace;

    //! MPI communicator
    boost::shared_ptr<Epetra_Comm> M_comm;
    Int                            M_me;

    //! Monodomain BC
    BCHandler*                     M_BCh_electric;
    bool                           M_setBC;

    //! Map
    EpetraMap                      M_localMap;
    EpetraMap                      M_localMapVec;

    //! mass matrix
    matrixPtr_Type                 M_matrMass;

    //! Stiff matrix: D*stiff
    matrixPtr_Type                 M_matrStiff;

    matrixPtr_Type                 M_matrNoBC;

    //! Right hand side for the PDE
    vector_Type                    M_rhsNoBC;

    //! Global solution _u
    vector_Type                    M_sol_u;

    //! Displacement
    vector_Type                    M_disp;

    //! Local fibers vector
    vector_Type                    M_fiber_vector;

    //! residual
    vector_Type                    M_residual;

    //! Solver
    SolverType                     M_linearSolver;

    prec_type                      M_prec;

    Real                         M_diagonalize;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //! Boolean that indicates if the matrix is updated for the current iteration
    bool                           M_updated;

    //! Boolean that indicates if the precond has to be recomputed
    bool                           M_reusePrec;
    bool                           M_resetPrec;

    //! Integer storing the max number of solver iteration with prec recomputing
    Int                            M_maxIterSolver;

    //! Boolean that indicates if the matrix has to be recomputed
    bool                           M_recomputeMatrix;



private:

    //! Elementary matrices
    ElemMat                        M_elmatStiff;
    ElemMat                        M_elmatMass;
    Real 						   massCoeff;
    UInt dim_u() const           { return M_uFESpace.dim(); }

}; // class HeartNonLinearMonodomainSolver



// ************************************************
// IMPLEMENTATION
// ************************************************

//! Constructors
template<typename Mesh, typename SolverType>
HeartNonLinearMonodomainSolver<Mesh, SolverType>::
HeartNonLinearMonodomainSolver( const data_type&          dataType,
                     FESpace<Mesh, EpetraMap>& uFESpace,
                     BCHandler&                BCh_u,
                     boost::shared_ptr<Epetra_Comm> comm):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_BCh_electric           ( &BCh_u ),
    M_setBC                  ( true ),
    M_comm                   ( comm ),
    M_me                     ( M_comm->MyPID() ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_localMap               ( M_uFESpace.map() ),
    M_localMapVec              (M_localMap+M_localMap+M_localMap),
    M_matrMass               ( ),
    M_matrStiff	             ( ),
    M_matrNoBC               ( ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, 1, 1 ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, 1, 1 ),
    M_rhsNoBC                ( M_localMap ),
    M_sol_u                  ( M_localMap ),
    M_disp                   ( M_localMapVec, Repeated ),
    M_fiber_vector           ( M_localMapVec, Repeated ),
    M_residual               ( M_localMap ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false )
{

 if (M_data.has_fibers() )
	    {
	    	std::stringstream MyPID;
	        ifstream fibers(M_data.fibers_file().c_str());

	        std::cout << "fiber_file: " <<  M_data.fibers_file().c_str() << std::endl;
	        UInt NumGlobalElements= M_localMapVec.getMap(Repeated)->NumGlobalElements();
	        std::vector<Real> fiber_global_vector(NumGlobalElements);

	        for( UInt i=0; i< NumGlobalElements; ++i)
	    		fibers>>fiber_global_vector[i];
	    	 UInt NumMyElements = M_localMapVec.getMap(Repeated)->NumMyElements();
	    	for(UInt j=0; j< NumMyElements; ++j)
	    	{
	    		UInt ig= M_localMapVec.getMap(Repeated)->MyGlobalElements()[j];
	    		M_fiber_vector[ig]= fiber_global_vector[ig-1];
	    		}
	    	std::cout << std::endl;
	    	fiber_global_vector.clear();
	    }

}

template<typename Mesh, typename SolverType>
HeartNonLinearMonodomainSolver<Mesh, SolverType>::
~HeartNonLinearMonodomainSolver()
{
}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{

    M_diagonalize = dataFile( "electric/space_discretization/diagonalize",  1 );

    M_reusePrec   = dataFile( "electric/prec/reuse", true);

    M_linearSolver.setDataFromGetPot( dataFile, "electric/solver" );

    M_maxIterSolver = dataFile( "electric/solver/max_iter", -1);

    std::string precType = dataFile( "electric/prec/prectype", "Ifpack");

    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "NonlinearMonodomain : Preconditioner not set");

    M_prec->setDataFromGetPot( dataFile, "electric/prec" );
}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::buildSystem()
{
    M_matrMass.reset  ( new matrix_Type(M_localMap) );
    M_matrStiff.reset( new matrix_Type(M_localMap) );

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;


    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoZero;

    vector_Type M_sol_uRep(M_sol_u, Repeated );

    M_comm->Barrier();

    chrono.start();

    //! Elementary computation and matrix assembling
    //! Loop on elements

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoDer.start();
        chronoDer.stop();
        chronoZero.start();

        M_elmatStiff.zero();
        M_elmatMass.zero();

        chronoZero.stop();
        chronoStiff.start();

        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
        if (M_data.has_fibers() )
        {
            stiffNL( M_sol_uRep,
                     M_data.sigmal(),
                     M_data.sigmat(),
                     M_fiber_vector,
                     M_elmatStiff,
                     M_uFESpace.fe(),
                     M_uFESpace.dof(),
                     0,
                     0 );
        }
        else
        {
            stiffNL( M_sol_uRep,
                     M_data.D(),
                     M_elmatStiff,
                     M_uFESpace.fe(),
                     M_uFESpace.dof(),
                     0,
                     0 );
            // Only used for the linear case:
            //stiff( M_data.D(), M_elmatStiff,  M_uFESpace.fe(), 0, 0 );
        }
        chronoStiff.stop();
        chronoMass.start();
        mass( 1., M_elmatMass, M_uFESpace.fe(), 0, 0 );
        chronoMass.stop();


        chronoStiffAssemble.start();
        assembleMatrix( *M_matrStiff,
                        M_elmatStiff,
                        M_uFESpace.fe(),
                        M_uFESpace.fe(),
                        M_uFESpace.dof(),
                        M_uFESpace.dof(),
                        0, 0,
                        0, 0);
        chronoStiffAssemble.stop();


        chronoMassAssemble.start();
        assembleMatrix( *M_matrMass,
                        M_elmatMass,
                        M_uFESpace.fe(),
                        M_uFESpace.fe(),
                        M_uFESpace.dof(),
                        M_uFESpace.dof(),
                        0, 0,
                        0, 0);
        chronoMassAssemble.stop();

    }

    massCoeff = M_data.Chi()*M_data.Cm() / M_data.timeStep();

    M_comm->Barrier();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    chrono.start();

    M_matrStiff->GlobalAssemble();
    M_matrMass->GlobalAssemble();

    M_comm->Barrier();

    M_matrNoBC.reset(new matrix_Type(M_localMap, M_matrStiff->getMeanNumEntries() ));

    //! Computing chi*Cm/dt * M + K

    *M_matrNoBC += *M_matrStiff;

    *M_matrNoBC += *M_matrMass*massCoeff;

    M_matrNoBC->GlobalAssemble();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;

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
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::
initialize( const source_type& u0 )
{

    vector_Type u(M_uFESpace.map());

    M_uFESpace.interpolate(u0, u, 0.);

    initialize(u);
}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::
initialize( const Function& u0 )
{

     vector_Type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, 0.);

     initialize(u);
}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::
initialize( const vector_Type& u0 )
{

    M_sol_u = u0;

}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::updatePDESystem(Real alpha,
                                                            vector_Type& sourceVec)
{
    Chrono chrono;

    if (M_verbose)
            std::cout << "  f-  Updating mass term and right hand side... "
                  << std::flush;

    chrono.start();

    M_rhsNoBC = sourceVec;
    M_rhsNoBC.GlobalAssemble();

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;

    M_updated = false;

    if (M_recomputeMatrix)
        buildSystem();

    if (M_verbose)
          std::cout << "  f-  Copying the matrices ...                 "
                    << std::flush;

    chrono.start();

    M_matrNoBC.reset(new matrix_Type(M_localMap, M_matrStiff->getMeanNumEntries() ));

    *M_matrNoBC += *M_matrStiff;

    *M_matrNoBC += *M_matrMass*alpha;

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                             << std::flush;

    M_updated = true;
    M_matrNoBC->GlobalAssemble();

}

template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::
updatePDESystem(vector_Type& sourceVec )
{

    Chrono chrono;

    if (M_verbose)
            std::cout << "  f-  Updating right hand side... "
                  << std::flush;

    chrono.start();

    M_rhsNoBC = sourceVec;
    M_rhsNoBC.GlobalAssemble();

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"  << std::flush;

}



template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::PDEiterate( bchandler_raw_type& bch )
{

    Chrono chrono;
    chrono.start();

    matrixPtr_Type matrFull( new matrix_Type(*M_matrNoBC) );
    vector_Type    rhsFull = M_rhsNoBC;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Updating boundary conditions ...         "
              << std::flush;

    chrono.start();

    applyBoundaryConditions( *matrFull, rhsFull, bch);

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    //! Solving the system
    solveSystem( matrFull, rhsFull );

    M_residual  = M_rhsNoBC;
    M_residual -= *M_matrNoBC*M_sol_u;

} // iterate()



template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::solveSystem( matrixPtr_Type  matrFull,
                                                         vector_Type&    rhsFull )
{
    Chrono chrono;
    for ( Int i = 0 ; i < rhsFull.getEpetraVector().MyLength() ; i++ )
    {
        Int ig=rhsFull.BlockMap().MyGlobalElements()[i];
    }
    if (M_verbose)
        std::cout << "  f-  Setting up the solver ...                ";

    chrono.start();
    M_linearSolver.setMatrix(*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    if ( !M_reusePrec || M_resetPrec )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f-  Computing the precond ...                ";

        M_prec->buildPreconditioner(matrFull);

        Real condest = M_prec->Condest();

        M_linearSolver.setPreconditioner(M_prec);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
          	std::cout << "         Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPrec = false;
    }
    else
    {
        if (M_verbose)
            std::cout << "  f-  Reusing  precond ...                \n" <<  std::flush;
    }


    chrono.start();

    if (M_verbose)
        std::cout << "  f-  Solving system ...                                ";

    Int numIter = M_linearSolver.solve(M_sol_u, rhsFull);
    chrono.stop();

    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }


    if (numIter > M_maxIterSolver)
    {
        M_resetPrec = true;
    }
    M_comm->Barrier();
}


template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::moveMesh(vector_Type const &dep)
{
    std::cout <<"  Moving the mesh ... "<< std::endl;
    M_uFESpace.mesh()->moveMesh(dep, this->M_uFESpace.dof().numTotalDof());
    std::cout << " done.\n"<< std::endl;
    //M_uFESpace->recomputeMatrix(true);
}




template<typename Mesh, typename SolverType>
void HeartNonLinearMonodomainSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_Type&        matrix,
                                                                     vector_Type&        rhs,
                                                                     bchandler_raw_type& BCh )
{
    if ( !BCh.bcUpdateDone() )
    {
        BCh.bcUpdate( *M_uFESpace.mesh(),
                      M_uFESpace.feBd(),
                      M_uFESpace.dof() );
    }

    vector_Type rhsFull(M_rhsNoBC,Repeated,Zero);

    bcManage( matrix,
              rhs,
              *M_uFESpace.mesh(),
              M_uFESpace.dof(),
              BCh,
              M_uFESpace.feBd(),
              1.,
              M_data.time() );

    rhs=rhsFull;

    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize( 1*dim_u(),
                            M_diagonalize,
                            rhs,
                            0.);
    }

} // applyBoundaryCondition


} // namespace LifeV


#endif //_HEARTNONLINEARMONODOMAINSOLVER_H_
