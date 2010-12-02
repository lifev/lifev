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
 * @brief File containing a solver class for the Monodomain equations in electrophysiology
 * @author Lucia Mirabella <lucia.mirabella@mail.polimi.it> and Mauro Perego <mauro.perego@polimi.it>
 * @date 11-2007
 * @contributors Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 * @mantainer Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
 * @last update 11-2010
 */

#ifndef _MONODOMAINSOLVER_H_
#define _MONODOMAINSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataMonodomain.hpp>
#include <boost/shared_ptr.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/stiffnessFibers.hpp>
#include <life/lifefem/bdf_template.hpp>

namespace LifeV
{

//! monodomainSolver - Class featuring the usual solver for monodomain equations

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class MonodomainSolver
{

public:

    //! @name Type definitions
    //@{

    typedef DataMonodomain data_type;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;
    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

    typedef typename SolverType::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;
    typedef typename SolverType::vector_type      vector_type;

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

    MonodomainSolver( const data_type& dataType,
                      FESpace<Mesh, EpetraMap>& uFESpace,
                      BCHandler& bcHandler,
                      Epetra_Comm& comm);

    //! Destructor
    virtual ~MonodomainSolver();

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
                              vector_type& sourceVec
                              );

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem( vector_type& sourceVec );

    //! Initialize
    void initialize( const source_type& );
    void initialize( const Function&  );
    void initialize( const vector_type& );

    //! Returns the local solution vector
    const vector_type& solution_u() const {return M_sol_u;}

    const vector_type& fiber_vector() const {return M_fiber_vector;}

    //! Returns the local residual vector
    const vector_type& residual() const {return M_residual;}

    //! Returns u FE space
    FESpace<Mesh, EpetraMap>& potentialFESpace() {return M_uFESpace;}

    //! Setting of the boundary conditions
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

    matrix_type& matrMass()
        {
            return *M_matrMass;
        }

    //@}
protected:

	//! Solves PDE system
    void solveSystem            (  matrix_ptrtype matrFull,
                                   vector_type&   rhsFull );

    //! Apply BC
    void applyBoundaryConditions(  matrix_type&        matrix,
                                   vector_type&        rhs,
                                   bchandler_raw_type& BCh);

    //! Data
    const data_type&               M_data;

    //! u FE space
    FESpace<Mesh, EpetraMap>&      M_uFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    Int                            M_me;

    //! Monodomain BC
    BCHandler*                     M_BCh_electric;
    bool                           M_setBC;

    //! Map
    EpetraMap                      M_localMap;
    EpetraMap                      M_localMapVec;

    //! mass matrix
    matrix_ptrtype                 M_matrMass;

    //! Stiff matrix: D*stiff
    matrix_ptrtype                 M_matrStiff;

    matrix_ptrtype                 M_matrNoBC;

    //! Right hand side for the PDE
    vector_type                    M_rhsNoBC;

    //! Global solution _u
    vector_type                    M_sol_u;

    //! Local fibers vector
    vector_type                    M_fiber_vector;

    //! residual
    vector_type                    M_residual;

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

}; // class MonodomainSolver



//
// IMPLEMENTATION
//

//! Constructors
template<typename Mesh, typename SolverType>
MonodomainSolver<Mesh, SolverType>::
MonodomainSolver( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& uFESpace,
       BCHandler&                BCh_u,
       Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_BCh_electric           ( &BCh_u ),
    M_setBC                  ( true ),
    M_localMap               ( M_uFESpace.map() ),
    M_localMapVec              (M_localMap+M_localMap+M_localMap),
    M_matrMass               ( ),
    M_matrStiff	             ( ),
    M_matrNoBC               ( ),
    M_rhsNoBC                ( M_localMap ),
    M_sol_u                  ( M_localMap ),
    M_fiber_vector           ( M_localMapVec, Repeated ),
    M_residual               ( M_localMap ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec              ( true ),
    M_resetPrec              ( true ),
    M_maxIterSolver          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, 1, 1 ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, 1, 1 )
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
};



template<typename Mesh, typename SolverType>
MonodomainSolver<Mesh, SolverType>::
~MonodomainSolver()
{
}


template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{

    M_diagonalize = dataFile( "electric/space_discretization/diagonalize",  1. );

    M_reusePrec   = dataFile( "electric/prec/reuse", true);

    M_linearSolver.setDataFromGetPot( dataFile, "electric/solver" );

    M_maxIterSolver = dataFile( "electric/solver/max_iter", -1);

    std::string precType = dataFile( "electric/prec/prectype", "Ifpack");

    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "monodomainSolver : Preconditioner not set");

    M_prec->setDataFromGetPot( dataFile, "electric/prec" );
}


template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::buildSystem()
{
    M_matrMass.reset  ( new matrix_type(M_localMap) );
    M_matrStiff.reset( new matrix_type(M_localMap) );

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;


    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoZero;

    M_comm->Barrier();

    chrono.start();

    //! Elementary computation and matrix assembling
    //! Loop on elements

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoZero.start();

        M_elmatStiff.zero();
        M_elmatMass.zero();

        chronoZero.stop();

        chronoStiff.start();
        switch(M_data.heart_diff_fct() )
        {
        	case 0:
        		M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
        		if (M_data.has_fibers() )
        		{
        			stiff( M_data.sigmal(), M_data.sigmat(), M_fiber_vector, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.dof(), 0, 0);
        		}
        		else
        		{
        			stiff( M_data.D(), M_elmatStiff,  M_uFESpace.fe(), 0, 0 );
        		}
	        break;

        	case 1:
        		M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
        	    if (M_data.has_fibers() )
        	    {
        	        stiff( M_data.red_sigma_sphere, M_data.sigmal(), M_data.sigmat(), M_fiber_vector, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.dof(), 0, 0);
        	    }
        	    else
        	    {
        	        stiff( M_data.red_sigma_sphere, M_data.D(), M_elmatStiff,  M_uFESpace.fe(), M_uFESpace.dof(), 0, 0 );
        	    }
        		break;
        	case 2:
        	    M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
        	    if (M_data.has_fibers() )
        	    {
        	        stiff( M_data.red_sigma_cyl, M_data.sigmal(), M_data.sigmat(), M_fiber_vector, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.dof(), 0, 0);
        	    }
        	    else
        	    {
        	        stiff( M_data.red_sigma_cyl, M_data.D(), M_elmatStiff,  M_uFESpace.fe(), M_uFESpace.dof(), 0, 0 );
        	    }
        	    break;
        	case 3:
        	    M_uFESpace.fe().updateFirstDerivQuadPt( M_uFESpace.mesh()->volumeList( iVol ) );
        	    if (M_data.has_fibers() )
        	    {
        	        stiff( M_data.red_sigma_box, M_data.sigmal(), M_data.sigmat(), M_fiber_vector, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.dof(), 0, 0);
        	    }
        	    else
        	    {
        	        stiff( M_data.red_sigma_box, M_data.D(), M_elmatStiff,  M_uFESpace.fe(), M_uFESpace.dof(), 0, 0 );
        	    }
        	    break;
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

    massCoeff = M_data.Chi() * M_data.Cm() / M_data.getTimeStep();

    M_comm->Barrier();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    chrono.start();

    M_matrStiff->GlobalAssemble();
    M_matrMass->GlobalAssemble();

    M_comm->Barrier();

    M_matrNoBC.reset(new matrix_type(M_localMap, M_matrStiff->getMeanNumEntries() ));

    //! Computing 1/dt * M + K

    *M_matrNoBC += *M_matrStiff;

    *M_matrNoBC += *M_matrMass*massCoeff;

    M_matrNoBC->GlobalAssemble();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s." << std::endl;

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diff_cumul() << " s.\n"
                  << " Zero           " << chronoZero.diff_cumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diff_cumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diff_cumul() << " s.\n"
                  << " Mass           " << chronoMass.diff_cumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diff_cumul() << " s.\n"
                  << std::endl;

}

template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::
initialize( const source_type& u0 )
{

    vector_type u(M_uFESpace.map());

    M_uFESpace.interpolate(u0, u, 0.);

    initialize(u);
}



template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::
initialize( const Function& u0 )
{

     vector_type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, 0.);

     initialize(u);
}


template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::
initialize( const vector_type& u0 )
{

    M_sol_u = u0;

}


template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::
updatePDESystem(Real       alpha,
             vector_type& sourceVec
             )
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

    M_matrNoBC.reset(new matrix_type(M_localMap, M_matrStiff->getMeanNumEntries() ));

    *M_matrNoBC += *M_matrStiff;

    *M_matrNoBC += *M_matrMass*alpha;


    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                             << std::flush;



    M_updated = true;
    M_matrNoBC->GlobalAssemble();

}

template<typename Mesh, typename SolverType>
void MonodomainSolver<Mesh, SolverType>::
updatePDESystem(vector_type& sourceVec )
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
void MonodomainSolver<Mesh, SolverType>::PDEiterate( bchandler_raw_type& bch )
{

    Chrono chrono;


    chrono.start();


    matrix_ptrtype matrFull( new matrix_type(*M_matrNoBC) );
    vector_type    rhsFull = M_rhsNoBC;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying boundary conditions ...         "
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
void MonodomainSolver<Mesh, SolverType>::solveSystem( matrix_ptrtype  matrFull,
                                           vector_type&    rhsFull )
{
    Chrono chrono;

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
          	std::cout << "Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        M_resetPrec = false;
    }
    else
    {
        if (M_verbose)
            std::cout << "f-  Reusing  precond ...                \n" <<  std::flush;
    }


    chrono.start();

    if (M_verbose)
        std::cout << "f-  Solving system ...                                ";

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
void MonodomainSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                       vector_type&        rhs,
                                                       bchandler_raw_type& BCh )
{

    // BC manage for the PDE
    if ( !BCh.bdUpdateDone() )
    {
        BCh.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    bcManage( matrix, rhs, *M_uFESpace.mesh(), M_uFESpace.dof(),
              BCh, M_uFESpace.feBd(), 1.,M_data.getTime() );


    if ( BCh.hasOnlyEssential() && M_diagonalize )
    {
        matrix.diagonalize( 1*dim_u(),
                            M_diagonalize,
                            rhs,
                            0.);
    }

} // applyBoundaryCondition


} // namespace LifeV


#endif //_MONODOMAINSOLVER_H_
