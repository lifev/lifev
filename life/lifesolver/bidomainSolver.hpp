/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
      Date: 2003-06-09

 Copyright (C) 2004 EPFL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
  \file bidomainSolver.hpp
  \author L. Mirabella M. Perego
  \date 11/2007

  \brief This file contains a Oseen equation solver class
*/
#ifndef _BIDOMAINSOLVER_H_
#define _BIODOMAINSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataBidomain.hpp>
#include <boost/shared_ptr.hpp>
#include <life/lifefem/FESpace.hpp>
#include "testsuite/test_heart/stiffness_fibers.hpp"

namespace LifeV
{
/*!
  \class BidomainSolver

  This class implements a bidomain solver.
*/

const UInt nbComp = 2;

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class BidomainSolver
{

public:

    typedef DataBidomain<Mesh> data_type;

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


    //! Constructor
    /*!
      \param dataType
      \param potential FE space
      \param bcHandler boundary conditions for potential
      \param Epetra communicator
    */
    BidomainSolver( const data_type&          dataType,
    	   FESpace<Mesh, EpetraMap>& FESpace,
           FESpace<Mesh, EpetraMap>& uFESpace,
           BCHandler&                bcHandler,
           Epetra_Comm&              comm );

    /*!
      \param dataType
      \param potential FE space
      \param Epetra communicator
    */
    BidomainSolver( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& FESpace,
           FESpace<Mesh, EpetraMap>& uFESpace,
           Epetra_Comm&              comm );

    //! virtual destructor
    virtual ~BidomainSolver();

    //! Updates sources, bc treatment and solve the bidomain system
    virtual void PDEiterate( bchandler_raw_type& bch );

    //! Sets up the system solver
    virtual void setUp        ( const GetPot& dataFile );

    //! Builds time independent parts of PDE system
    virtual void buildSystem();

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem(double       alpha,
                              vector_type& sourceVec
                              );

    //! Updates time dependent parts of PDE system
    virtual void updatePDESystem( vector_type& sourceVec );

    //! Initialize
    void initialize( const Function&, const Function&  );
    void initialize( const vector_type&, const vector_type& );

    //! Returns the local solution vector
    const vector_type& solution_u() const {return M_sol_u;}
    const vector_type& solution_uiue() const {return M_sol_uiue;}

    const vector_type& fiber_vector() const {return M_fiber_vector;}

    //! Returns the local residual vector
    const vector_type& residual() const {return M_residual;}

    //! Returns u FE space
    FESpace<Mesh, EpetraMap>& potentialFESpace() {return M_uFESpace;}

    //! Boundary Conditions
    const bool BCset() const {return M_setBC;}

    //! Sets Bidomain BCs
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
    FESpace<Mesh, EpetraMap>&      M_FESpace;
    FESpace<Mesh, EpetraMap>&      M_uFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    int                            M_me;

    //! Bidomain BC
    BCHandler*                     M_BCh_electric;
    bool                           M_setBC;

    //! Map
    EpetraMap                      M_localMap;
    EpetraMap                      M_localMap_u;
    EpetraMap                      M_localMapVec;

    //! mass matrix
    matrix_ptrtype                 M_matrMass;

    //! Stiff matrix: D*stiff
    matrix_ptrtype                 M_matrStiff;

    matrix_ptrtype                 M_matrNoBC;

    //! Right hand side for the PDE
    vector_type                    M_rhsNoBC;

    //! Global solution _u
    vector_type                    M_sol_uiue;

    vector_type                    M_sol_u;

    //! Local fibers vector
    vector_type                    M_fiber_vector;

    //! residual
    vector_type                    M_residual;

    //! Solver
    SolverType                     M_linearSolver;

    prec_type                      M_prec;

    double                         M_diagonalize;

    //! Boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //! Boolean that indicates if the matrix is updated for the current iteration
    bool                           M_updated;

    //! Boolean that indicates if the precond has to be recomputed
    bool                           M_reusePrec;
    bool                           M_resetPrec;

    //! Integer storing the max number of solver iteration with prec recomputing
    int                            M_maxIterSolver;

    //! Boolean that indicates if the matrix has to be recomputed
    bool                           M_recomputeMatrix;

private:

    //! Elementary matrices
    ElemMat                        M_elmatStiff;
    ElemMat                        M_elmatMass;
    Real 						   massCoeff;
    UInt dim_u() const           { return M_uFESpace.dim(); }

}; // class BidomainSolver



//
// IMPLEMENTATION
//

//! Constructors
template<typename Mesh, typename SolverType>
BidomainSolver<Mesh, SolverType>::
BidomainSolver( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& FESpace,
       FESpace<Mesh, EpetraMap>& uFESpace,
       BCHandler&                BCh_u,
       Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_FESpace               ( FESpace ),
    M_uFESpace               ( uFESpace ),
    M_BCh_electric           ( &BCh_u ),
    M_setBC                  ( true ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_linearSolver           ( ),
    M_prec                   ( ),
    M_localMap               ( M_FESpace.map()),
    M_localMap_u             ( M_uFESpace.map()),
    M_localMapVec            (M_localMap_u+M_localMap_u+M_localMap_u),
    M_matrMass               ( ),
    M_matrStiff	             ( ),
    M_matrNoBC               ( ),
    M_elmatStiff             ( M_FESpace.fe().nbNode, 2, 2 ),
    M_elmatMass              ( M_FESpace.fe().nbNode, 2, 2 ),
    M_rhsNoBC                ( M_localMap ),
    M_sol_uiue                ( M_localMap ),
    M_sol_u                  ( M_localMap_u ),
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
BidomainSolver<Mesh, SolverType>::
BidomainSolver( const data_type&          dataType,
       FESpace<Mesh, EpetraMap>& FESpace,
       FESpace<Mesh, EpetraMap>& uFESpace,
       Epetra_Comm&              comm ):
    M_data               	( dataType ),
    M_FESpace              ( FESpace ),
    M_uFESpace              ( uFESpace ),
    M_setBC                 ( false ),
    M_comm                  ( &comm ),
    M_me                    ( M_comm->MyPID() ),
    M_linearSolver          ( ),
    M_prec                  ( new prec_raw_type() ),
    M_localMap              ( M_FESpace.map()),
    M_localMap_u            ( M_uFESpace.map() ),
    M_matrMass              ( ),
    M_matrStiff             ( ),
    M_matrNoBC              ( ),
    M_elmatStiff            ( M_FESpace.fe().nbNode, 2, 2 ),
    M_elmatMass             ( M_FESpace.fe().nbNode, 2, 2 ),
    M_rhsNoBC               ( M_localMap ),
    M_residual              ( M_localMap ),
    M_sol_uiue               ( M_localMap ),
    M_sol_u                 ( M_localMap_u ),
    M_verbose               ( M_me == 0),
    M_updated               ( false ),
    M_reusePrec             ( true ),
    M_resetPrec             ( true ),
    M_maxIterSolver         ( -1 ),
    M_recomputeMatrix       ( false )
{
}



template<typename Mesh, typename SolverType>
BidomainSolver<Mesh, SolverType>::
~BidomainSolver()
{
}


template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::setUp( const GetPot& dataFile )
{

    M_diagonalize = dataFile( "electric/discretization/diagonalize",  1. );

    M_reusePrec   = dataFile( "electric/prec/reuse", true);

    M_linearSolver.setDataFromGetPot( dataFile, "electric/solver" );

    M_maxIterSolver = dataFile( "electric/solver/max_iter", -1);


    std::string precType = dataFile( "electric/prec/prectype", "Ifpack");

    M_prec.reset( PRECFactory::instance().createObject( precType ) );
    ASSERT(M_prec.get() != 0, "bidomainSolver : Preconditioner not set");

    M_prec->setDataFromGetPot( dataFile, "electric/prec" );

}


template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::buildSystem()
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

    for ( UInt iVol = 1; iVol <= M_FESpace.mesh()->numVolumes(); iVol++ )
    {




        chronoZero.start();

        M_elmatStiff.zero();
        M_elmatMass.zero();

        chronoZero.stop();

        chronoStiff.start();
        switch(M_data.heart_diff_fct() )
        {
        	case 0:
        		chronoDer.start();
        		M_FESpace.fe().updateFirstDeriv( M_FESpace.mesh()->volumeList( iVol ) );
        		chronoDer.stop();
        		if (M_data.has_fibers() )
        		{
        			stiff( M_data.sigmal_i(), M_data.sigmat_i(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 0, 0);
        			stiff( M_data.sigmal_e(), M_data.sigmat_e(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 1, 1);
        		}
        		else
        		{
        			stiff( M_data.D_i(), M_elmatStiff,  M_FESpace.fe(), 0, 0 );
        			stiff( M_data.D_e(), M_elmatStiff,  M_FESpace.fe(), 1, 1 );
	        	}
	        break;
        	case 1:
        		chronoDer.start();
        		M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( iVol ) );
        		chronoDer.stop();
        		if (M_data.has_fibers() )
        	    {
        	    	stiff( M_data.red_sigma_sphere,  M_data.sigmal_i(), M_data.sigmat_i(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 0, 0, 0);
        	    	stiff( M_data.red_sigma_sphere, M_data.sigmal_e(), M_data.sigmat_e(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1);
        	    }
        	    else
        	    {
        	    	stiff( M_data.red_sigma_sphere, M_data.D_i(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 0, 0, 0);
        	        stiff( M_data.red_sigma_sphere, M_data.D_e(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1);
        		}
        	break;
        	case 2:
        	    chronoDer.start();
        	    M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( iVol ) );
        	    chronoDer.stop();
        	    if (M_data.has_fibers() )
        	    {
        	       	stiff( M_data.red_sigma_cyl,  M_data.sigmal_i(), M_data.sigmat_i(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 0, 0, 0);
        	      	stiff( M_data.red_sigma_cyl, M_data.sigmal_e(), M_data.sigmat_e(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1);
        	    }
        	    else
        	    {
        	       	stiff( M_data.red_sigma_cyl, M_data.D_i(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 0, 0 , 0);
        	        stiff( M_data.red_sigma_cyl, M_data.D_e(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1);
        	    }
        	break;
        	case 3:
        	    chronoDer.start();
        	    M_FESpace.fe().updateFirstDerivQuadPt( M_FESpace.mesh()->volumeList( iVol ) );
        	    chronoDer.stop();
        	    if (M_data.has_fibers() )
        	    {
        	      	stiff( M_data.red_sigma_box, M_data.sigmal_i(), M_data.sigmat_i(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 0, 0, 0);
        	       	stiff( M_data.red_sigma_box, M_data.sigmal_e(), M_data.sigmat_e(), M_fiber_vector, M_elmatStiff, M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1);
        	    }
        	    else
        	    {
        	       	stiff( M_data.red_sigma_box, M_data.D_i(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 0, 0, 0 );
        	        stiff( M_data.red_sigma_box, M_data.D_e(), M_elmatStiff,  M_FESpace.fe(), M_FESpace.dof(), 1, 1, 1 );
        	    }
        	    break;
        }
        chronoStiff.stop();

        chronoMass.start();
        mass( 1., M_elmatMass, M_FESpace.fe(), 0, 0 );
        mass( -1., M_elmatMass, M_FESpace.fe(), 0, 1 );
        mass( -1., M_elmatMass, M_FESpace.fe(), 1, 0 );
        mass( 1., M_elmatMass, M_FESpace.fe(), 1, 1 );
        chronoMass.stop();


        chronoStiffAssemble.start();
        for ( UInt iComp = 0; iComp < nbComp; iComp++ ){

        	assembleMatrix( *M_matrStiff,
        	                        M_elmatStiff,
        	                        M_FESpace.fe(),
        	                        M_FESpace.fe(),
        	                        M_FESpace.dof(),
        	                        M_FESpace.dof(),
        	                        iComp, iComp,
        	                        iComp*dim_u(), iComp*dim_u());
        }
        chronoStiffAssemble.stop();


        chronoMassAssemble.start();

        for ( UInt iComp = 0; iComp < nbComp; iComp++ ){
        	for ( UInt jComp = 0; jComp < nbComp; jComp++ ){
        		assembleMatrix( *M_matrMass,
        		                        M_elmatMass,
        		                        M_FESpace.fe(),
        		                        M_FESpace.fe(),
        		                        M_FESpace.dof(),
        		                        M_FESpace.dof(),
        		                        iComp, jComp,
        		                        iComp*dim_u(), jComp*dim_u());
        	}
         }
         chronoMassAssemble.stop();

    }

    massCoeff = 1.0/ M_data.getTimeStep();

    M_comm->Barrier();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    if (M_verbose) std::cout << "  f-  Finalizing the matrices     ...        " << std::flush;
    chrono.start();

    M_matrStiff->GlobalAssemble();
    M_matrMass->GlobalAssemble();

    M_comm->Barrier();

    M_matrNoBC.reset(new matrix_type(M_localMap, M_matrStiff->getMeanNumEntries() ));

    //! Computing 1.0/dt * M + K

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
void BidomainSolver<Mesh, SolverType>::
initialize( const Function& ui0, const Function& ue0 )
{

     vector_type ui(M_uFESpace.map());
     vector_type ue(M_uFESpace.map());

     M_uFESpace.interpolate(ui0, ui, 0.);
     M_uFESpace.interpolate(ue0, ue, 0.);

     initialize(ui, ue);
}


template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::
initialize( const vector_type& ui0, const vector_type& ue0 )
{
   M_sol_uiue = ui0;
   M_sol_uiue.add(ue0, M_uFESpace.dof().numTotalDof());
    for ( int i = 0 ; i < M_sol_u.getEpetraVector().MyLength() ; i++ )
   	{
   		int ig=M_sol_u.BlockMap().MyGlobalElements()[i];
        M_sol_u[ig] = M_sol_uiue[ig] - M_sol_uiue[ig+dim_u()]; // BASEINDEX + 1
               }
}


template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::
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
void BidomainSolver<Mesh, SolverType>::
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
void BidomainSolver<Mesh, SolverType>::PDEiterate( bchandler_raw_type& bch )
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
    M_residual -= *M_matrNoBC*M_sol_uiue;

} // iterate()



template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::solveSystem( matrix_ptrtype  matrFull,
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
        std::cout << "Before condest"<<std::endl<<std::flush;
        double condest = M_prec->Condest();
        std::cout << "After condest"<<std::endl<<std::flush;

        M_linearSolver.setPreconditioner(M_prec);
        std::cout << "After setPreconditioner"<<std::endl<<std::flush;

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

    int numIter = M_linearSolver.solve(M_sol_uiue, rhsFull);

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



    for ( UInt i = 0 ; i < M_sol_u.getEpetraVector().MyLength() ; i++ )
       	{
       		UInt ig=M_sol_u.BlockMap().MyGlobalElements()[i];
            M_sol_u[ig] = M_sol_uiue[ig] - M_sol_uiue[ig+dim_u()]; // BASEINDEX + 1
                   }



    M_comm->Barrier();

}

//occhio: to check!

template<typename Mesh, typename SolverType>
void BidomainSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                       vector_type&        rhs,
                                                       bchandler_raw_type& BCh )
{

    // BC manage for the PDE
    if ( !BCh.bdUpdateDone() )
    {
        BCh.bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );
    }

  //  vector_type rhsFull(M_rhsNoBC,Repeated, Zero);


//    rhsFull.Import(M_rhsNoBC, Zero); // ignoring non-local entries, Otherwise they are summed up lately

    bcManage( matrix, rhs, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.,
              M_data.getTime() );

  //  rhs = rhsFull;


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
