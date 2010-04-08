/**
 * This file is part of the LifeV library
 * Copyright (C) 2009 INRIA, EPFL, Politecnico Milano
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*! 
 * \file ChorinTemamRK.hpp
 * \desc Nearly explicit Chorin-Temam method with: 
 *       order 2 Runge-Kutta like scheme for time discretization and 
 *       explicit interior penalty for convective stabilization.
 */

#ifndef _CHORIN_TEMAM_RK_H_
#define _CHORIN_TEMAM_RK_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>

#include <boost/shared_ptr.hpp>

#include <life/lifefem/FESpace.hpp>

#include <ipStabilization.hpp>
#include <elemOperCT.hpp>

namespace LifeV
{

template< typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >
class ChorinTemamRK

{

public:

    typedef DataNavierStokes<Mesh> data_type;

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

    typedef enum {STEP_1, STEP_2}                 step_type;

    //! Constructor 
    /*!
      \param dataType
      \param velocity FE space
      \param pressure FE space
      \param bcHu boundary conditions for the velocity
      \param bcHp boundary conditions for the pressure
    */
    ChorinTemamRK( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           BCHandler&                bcHu,
	   BCHandler& 		     bcHp,
           Epetra_Comm&              comm );

    //! virtual destructor
    virtual ~ChorinTemamRK();

    // Euler explicit and Crank-Nicholson two RK steps time loop methods
    
    // {
    //! update u system convection and ip stab terms, 
    //  apply u bc's, solve u system, and (if step 2) compute p system rhs
    virtual void iterate_u (bchandler_raw_type& bch_u, step_type step);

    //! apply p bc's, solve p system
    virtual void iterate_p (bchandler_raw_type& bch_p);

    //! update rhs for u system i.e. Chorin-Temam coupling and mass term
    void time_advance (const Real& time, step_type step);
    // }

    //! setup physical and numerical parameters 
    virtual void setUp        ( const GetPot& dataFile );

    void initialize( const Function&, const Function&  );
    void initialize( const vector_type&, const vector_type& );

    //! returns the local solution vector
    const vector_type& solution_u() const {return M_sol_u;}
    const vector_type& solution_p() const {return M_sol_p;}

    //! returns the local residual vector
    const vector_type& residual_u() const {return M_residual_u;}
    const vector_type& residual_p() const {return M_residual_p;}

    FESpace<Mesh, EpetraMap>& velFESpace()   {return M_uFESpace;}
    FESpace<Mesh, EpetraMap>& pressFESpace() {return M_pFESpace;}


    //! Boundary Conditions
    const bool BCset() const {return M_setBC;}
    //! set the fluid BCs
    void setBC(BCHandler &BCh_u, BCHandler &BCh_p)
        {
            M_BCh_fluid_u = &BCh_u; 
	    M_BCh_fluid_p = &BCh_p;
	    M_setBC = true;
        }

    //! set the source term functor
    void setSourceTerm( source_type __s )
        {
            M_source = __s;
        }

    //! get the source term functor
    source_type sourceTerm() const
        {
            return M_source;
        }


    //! Postprocessing
    void postProcess(bool _writeMesh = false);

    void resetPrec_u() {M_resetPrec_u = true;}
    void resetPrec_p() {M_resetPrec_p = true;}

    EpetraMap const& getMap_u() const { return M_localMap_u; }
    EpetraMap const& getMap_p() const { return M_localMap_p; }

    const Epetra_Comm& comm() const {return *M_comm;}

    void recomputeMatrix(bool const recomp){M_recomputeMatrix = recomp;}

    matrix_type& matrNoBC_u()
        {
            return *M_matrNoBC_u;
        }
    matrix_type& matrNoBC_p()
	{
	    return *M_matrNoBC_p;
	} 
    matrix_type& matrMass()
        {
            return *M_matrMass;
        }


protected:

    UInt dim_u() const           { return M_uFESpace.dim(); }
    UInt dim_p() const           { return M_pFESpace.dim(); }

    // build u and p systems constant matrices
    void buildSystem_u_p	();
    // update convective and stab terms for u system
    // rationale : all in one method as the convector betaVec is the same for convection 
    // and stabilization in Euler or Crank-Nicholson step.
    void updateSYS_u (vector_type& betaVec, vector_type& u_rhs, step_type step);
    // compute u like system CT rhs coupling
    void computeCTRHS_u (vector_type& p_vec, vector_type& u_rhs);
    // compute p like system CT rhs coupling 
    void computeCTRHS_p (vector_type& u_vec, vector_type& p_rhs, const Real dt);    

    void solveSystem_u            ( matrix_ptrtype matrFull_u,
                                    vector_type&   rhsFull_u,
				    vector_type&   solFull_u );

    void solveSystem_p 		  ( matrix_ptrtype matrFull_p,
				    vector_type&   rhsFull_p,
				    vector_type&   solFull_p );

    void applyBC_u (matrix_type& matrix_u,
                    vector_type& rhs_u,
                    bchandler_raw_type& BCh_u);
    void applyBC_p (matrix_type& matrix_p,
                    vector_type& rhs_p,
                    bchandler_raw_type& BCh_p );
 
    void echo(std::string message);

    //! removes mean of component comp of vector x
    void removeMean( vector_type& x, UInt comp = 1 );

    //! data for NS solvers
    const data_type&               M_data;

    // FE spaces
    FESpace<Mesh, EpetraMap>&      M_uFESpace;
    FESpace<Mesh, EpetraMap>&      M_pFESpace;

    //! MPI communicator
    Epetra_Comm*                   M_comm;
    int                            M_me;

    //! fluid BC
    BCHandler*                     M_BCh_fluid_u;
    BCHandler*			   M_BCh_fluid_p;
    bool                           M_setBC;

    EpetraMap                      M_localMap_u;
    EpetraMap			   M_localMap_p;

    //! mass matrix
    matrix_ptrtype                 M_matrMass;

    //! Stokes matrix: mu*stiff
    matrix_ptrtype                 M_matrStokes;

    //! pressure matrix 
    matrix_ptrtype		   M_matrPress;

    //! matrix without boundary conditions
    matrix_ptrtype                 M_matrNoBC_u;
    matrix_ptrtype		   M_matrNoBC_p;
    
    //! source term for NS
    source_type                    M_source;

    //! Right hand side for the velocity
    vector_type                    M_rhsNoBC_u;
    vector_type			   M_rhsNoBC_p;

    //! Solutions u and p, and intermediates
    vector_type                    M_sol_u;
    vector_type			   M_sol_p;
    vector_type                    M_sol_u_prev;
    vector_type                    M_sol_u_aux;

    //! Residuals
    vector_type                    M_residual_u;
    vector_type			   M_residual_p;

    SolverType                     M_linearSolver_u;
    SolverType			   M_linearSolver_p;

    boost::shared_ptr<EpetraPreconditioner>             M_prec_u;
    boost::shared_ptr<EpetraPreconditioner>		M_prec_p;

    bool                           M_steady;

    //! IP Stabilization
    bool                           M_stab;
    IPStabilization<Mesh, Dof>     M_ipStab;
    Real                           M_gammaBeta;

    // boolean indicating variational pressure term (velocity system RHS)
    bool                           M_hasVariationalPressure;

    const Function*                M_betaFct;

    bool                           M_divBetaUv;

    double                         M_diagonalize_u;
    double 			   M_diagonalize_p;

    UInt                           M_count;

    //! boolean that indicates if output is sent to cout
    bool                           M_verbose;

    //! boolean that indicates if the matrix is updated for the current iteration
    bool                           M_updated;

    //! boolean that indicates if the precond has to be recomputed
    bool                           M_reusePrec_u;
    int                            M_maxIterForReuse_u;
    bool                           M_resetPrec_u;
    bool 			   M_reusePrec_p;
    int				   M_maxIterForReuse_p;
    bool 			   M_resetPrec_p;

    //! integer storing the max number of solver iteration with prec recomputing
    int                            M_maxIterSolver_u;
    int 			   M_maxIterSolver_p;

    //!
    bool                           M_recomputeMatrix;

private:

    //! Elementary matrices and vectors
    ElemMat                        M_elmatStiff;      // velocity Stokes
    ElemMat                        M_elmatMass;       // velocity mass
    ElemMat			   M_elmatPress;      // pressure laplacian
    ElemVec                        M_elvec_u;         // Elementary velocity vector
    ElemVec                        M_elvec_u_aux;
    ElemVec			   M_elvec_p;         // Elementary pressure vector
    ElemVec			   M_elemrhs_u;
    ElemVec			   M_elemrhs_p;		 

    //! Time related data
    Real			   M_time;	       // current time
    bool			   M_firstTimeStep; 
}; // class ChorinTemamRK


/****************************************************************************
 *
 * Implementation
 *
 ****************************************************************************/


template<typename Mesh, typename SolverType>
ChorinTemamRK<Mesh, SolverType>::
ChorinTemamRK( const data_type&   dataType,
	FESpace<Mesh, EpetraMap>& uFESpace,
	FESpace<Mesh, EpetraMap>& pFESpace,
	BCHandler&                BCh_u,
	BCHandler&	          BCh_p,
	Epetra_Comm&              comm ):
    M_data                   ( dataType ),
    M_uFESpace               ( uFESpace ),
    M_pFESpace               ( pFESpace ),
    M_BCh_fluid_u            ( &BCh_u ),
    M_BCh_fluid_p            ( &BCh_p ),
    M_setBC                  ( true ),
    M_comm                   ( &comm ),
    M_me                     ( M_comm->MyPID() ),
    M_linearSolver_u         ( ),
    M_linearSolver_p	     ( ),
    M_prec_u                 ( ),
    M_prec_p		     ( ),
    M_localMap_u             ( M_uFESpace.map() ),
    M_localMap_p	     ( M_pFESpace.map() ),
    M_matrMass               ( ),
    M_matrStokes             ( ),
    M_matrNoBC_u             ( ),
    M_matrNoBC_p	     ( ),
    M_elmatStiff             ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatMass              ( M_uFESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatPress             ( M_pFESpace.fe().nbNode, 1, 1 ),
    M_elvec_u                ( M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_u_aux            ( M_uFESpace.fe().nbNode, nDimensions ),
    M_elvec_p		     ( M_pFESpace.fe().nbNode, 1 ),
    M_elemrhs_u		     ( M_uFESpace.fe().nbNode, nDimensions ),
    M_elemrhs_p              ( M_pFESpace.fe().nbNode, 1 ),
    M_rhsNoBC_u              ( M_localMap_u ),
    M_rhsNoBC_p		     ( M_localMap_p ),
    M_sol_u                  ( M_localMap_u ),
    M_sol_u_prev             ( M_localMap_u ),
    M_sol_u_aux		     ( M_localMap_u ),
    M_sol_p                  ( M_localMap_p ),
    M_residual_u             ( M_localMap_u ),
    M_residual_p 	     ( M_localMap_p ), 
    M_stab                   ( false ),
    M_ipStab                 ( M_uFESpace.mesh(), 
                               M_uFESpace.dof(),
                               M_uFESpace.refFE(),
			       M_uFESpace.feBd(),
                               M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity()
                             ),
    M_betaFct                ( 0 ),
    M_count                  ( 0 ),
    M_verbose                ( M_me == 0),
    M_updated                ( false ),
    M_reusePrec_u              ( true ),
    M_reusePrec_p	       ( true ),
    M_maxIterForReuse_u        ( -1 ),
    M_maxIterForReuse_p        ( -1 ),
    M_resetPrec_u            ( true ),
    M_resetPrec_p	     ( true ),
    M_maxIterSolver_u          ( -1 ),
    M_maxIterSolver_p          ( -1 ),
    M_recomputeMatrix        ( false ),
    M_time		     (0),
    M_firstTimeStep	     (1),
    M_hasVariationalPressure ( false )
{
     M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE()); 
}

template<typename Mesh, typename SolverType>
ChorinTemamRK<Mesh, SolverType>::
~ChorinTemamRK()
{
}


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_steady      = dataFile( "fluid/miscellaneous/steady",        1  );
    M_gammaBeta   = dataFile ( "fluid/ipstab/gammaBeta", 0. );
//    M_gammaDiv    = dataFile ( "fluid/ipstab/gammaDiv", 0.);
    M_divBetaUv   = dataFile( "fluid/discretization/div_beta_u_v", 0  );
    M_diagonalize_u = dataFile( "fluid/discretization/diagonalizeVel",  1. );
    M_diagonalize_p = dataFile( "fluid/discretization/diagonalizePress",  1. );

    M_hasVariationalPressure = dataFile( "fluid/discretization/variationalPn", false);
    if (M_verbose) {
      if (M_hasVariationalPressure)
        std::cout << "\n  f- Chorin-Temam coupling : variational pressure" << std::endl;
      else
        std::cout << "\n  f- Chorin-Temam coupling : standard" << std::endl;
    }

    M_linearSolver_u.setDataFromGetPot( dataFile, "fluid/solver" );
    M_linearSolver_p.setDataFromGetPot( dataFile, "fluid/solver" );

    // Fill stabilization arguments
    M_ipStab.setGammaBeta(M_gammaBeta);
    M_ipStab.setGammaDiv(0.);	// see if we will need this later
    M_ipStab.setGammaPress(0.);
    if (M_verbose)
      std::cout << "  f- Stabilization convective parm value : " << 
        M_gammaBeta << std::endl;
     

    M_maxIterSolver_u   = dataFile( "fluid/solver/max_iter", -1);
    M_reusePrec_u       = dataFile( "fluid/prec/reuse", true);
    M_maxIterForReuse_u = dataFile( "fluid/prec/max_iter_reuse", M_maxIterSolver_u*8/10);
    M_maxIterSolver_p   = dataFile( "fluid/solver/max_iter", -1);
    M_reusePrec_p       = dataFile( "fluid/prec/reuse", true);
    M_maxIterForReuse_p = dataFile( "fluid/prec/max_iter_reuse", M_maxIterSolver_p*8/10);

    std::string precType = dataFile( "fluid/prec/prectype", "Ifpack" );

    M_prec_u.reset( PRECFactory::instance().createObject(precType) );
    M_prec_p.reset( PRECFactory::instance().createObject(precType) );  

    M_prec_u->setDataFromGetPot( dataFile, "fluid/prec" );
    M_prec_p->setDataFromGetPot( dataFile, "fluid/prec" );
} // setUp

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::buildSystem_u_p()
{

    // We compute u and p systems _constant_ matrices at the same time
    // one loop on geom element w/ possibly different fe for u and p

    M_matrMass.reset  ( new matrix_type(M_localMap_u) );
    M_matrStokes.reset( new matrix_type(M_localMap_u) );
    M_matrPress.reset ( new matrix_type(M_localMap_p) );

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    // See if we are steady
    if (M_verbose) {
	if (M_steady)
		std::cout << "  f-  Steady state ...." << std::endl;
        else
		std::cout << "  f-  Unsteady state ...." << std::endl;
    }

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoStiff;
    Chrono chronoMass;
    Chrono chronoPress;

    Chrono chronoStiffAssemble;
    Chrono chronoMassAssemble;
    Chrono chronoPressAssemble; 
    Chrono chronoStab;
    Chrono chronoZero;

    // Number of velocity components
    UInt nbCompU = nDimensions;

    // Elementary computation and matrix assembling
    // Loop on elements

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
    //UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    chrono.start();

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        chronoDer.start();
        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
        M_pFESpace.fe().updateFirstDeriv( M_pFESpace.mesh()->volumeList( iVol ) );
        chronoDer.stop();

        chronoZero.start();
        M_elmatStiff.zero();
        M_elmatMass.zero();
        M_elmatPress.zero();
        chronoZero.stop();

        // velocity stiffness strain elemental matrix
        chronoStiff.start();
        stiff_strain( 2.0*M_data.viscosity(), M_elmatStiff, M_uFESpace.fe() );
        //stiff( M_data.viscosity(), M_elmatStiff,  M_uFESpace.fe(), 0, 0, nDimensions );
        chronoStiff.stop();

        // velocity mass elemental matrix (1/dt term is taken into account in the matrix)
        if ( !M_steady )
        {
            chronoMass.start();
            mass( M_data.density()/M_data.dataTime()->getTimeStep(), M_elmatMass, M_uFESpace.fe(), 0, 0, nDimensions );
            chronoMass.stop();
        }

	// pressure elemental matrix
        chronoPress.start();
        stiff( 1., M_elmatPress, M_pFESpace.fe() );
        chronoPress.stop();        

        // velocity matrix assembling w/ loop on velocity 3 components 
	for ( UInt iComp = 0; iComp < nbCompU; iComp++ )
        {

                chronoStiffAssemble.start();
                assembleMatrix( *M_matrStokes,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
                chronoStiffAssemble.stop();


            if ( !M_steady )
            {
                chronoMassAssemble.start();
                assembleMatrix( *M_matrMass,
                                M_elmatMass,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof);
		chronoMassAssemble.stop();
            }
	
	}		

	// pressure matrix assembling
		chronoPressAssemble.start();
		assembleMatrix( *M_matrPress, 
				M_elmatPress, 
				M_pFESpace.fe(),
				M_pFESpace.fe(),
				M_pFESpace.dof(),
				M_pFESpace.dof(),
				0, 0, 0, 0);
		chronoPressAssemble.stop();

    }

    // be sure assembling has completed
    M_comm->Barrier();

    chrono.stop();
    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;


    if (M_verbose) std::cout << "  f-  Finalizing the matrices ... " << std::flush;
    chrono.start();

    M_matrStokes->GlobalAssemble();
    M_matrMass->GlobalAssemble();
    M_matrPress->GlobalAssemble();

    chrono.stop();
    if (M_verbose) std::cout << " done in " << chrono.diff() << " s." << std::endl;

} // buildSystem_u_p


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::
initialize( const Function& u0, const Function& p0 )
{
     vector_type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, 0.0);

     vector_type p(M_pFESpace.map());
     M_pFESpace.interpolate(p0, p, 0.0);

     initialize(u, p);
}


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::
initialize( const vector_type& u0, const vector_type& p0 )
{
    M_sol_u = u0;
    M_sol_p = p0;
    M_sol_u_prev *= 0.;
    M_sol_u_aux *= 0.;
}

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::updateSYS_u(vector_type& betaVec,
                                                   vector_type& u_rhs,
						   step_type step)
{
	
    Real dt = M_data.dataTime()->getTimeStep();

    // first update matrices for velocity system

    UInt nbCompU = nDimensions;
    UInt velTotalDof = M_uFESpace.dof().numTotalDof();
    Chrono chrono;

    double normInf;
    betaVec.NormInf(&normInf);

    if (M_verbose)
        std::cout << "  f-  Copying the matrices ..." << std::endl;

    M_matrNoBC_u.reset(new matrix_type(M_localMap_u, M_matrStokes->getMeanNumEntries() ));
    M_matrNoBC_p.reset(new matrix_type(M_localMap_p, M_matrPress->getMeanNumEntries() ));
    *M_matrNoBC_u += *M_matrStokes;
    *M_matrNoBC_p += *M_matrPress;

    if (step == STEP_1)
        *M_matrNoBC_u += *M_matrMass;
    if (step == STEP_2) {
        // divide by 2 the stokes matrix to account for Crank-Nicholson step
        *M_matrNoBC_u *= 0.5;
        *M_matrNoBC_u += *M_matrMass;
    }

    if (normInf == 0.) 
	if (M_verbose)
	  std::cout << "  f-  Convective velocity is zero  " << std::endl;

// SDEB
//    normInf = 0.;
// EDEB

    // update (explicit) convective terms and
    // viscous term when in the Crank-Nicholson step

    if (normInf != 0.)
    {

        if (M_verbose)
            std::cout << "  f-  Updating the convective terms ...        "
                      << std::flush;

        // vector with repeated nodes over the processors
	vector_type betaVecRep( betaVec, Repeated );
	vector_type velVecRep (M_sol_u_prev, Repeated);
	vector_type _u_rhs (u_rhs, Repeated, Zero);
	_u_rhs *= 0.;

        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
            M_elemrhs_u.zero();

            UInt eleID = M_uFESpace.fe().currentLocalId();
            
	    // get elemental vectors on current element
	    for ( UInt iNode = 0 ; iNode < ( UInt ) M_uFESpace.fe().nbNode ; iNode++ )
            {
                UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                {
                    UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*dim_u();
                    M_elvec_u_aux.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = betaVecRep[ig]; // BASEINDEX + 1
		    M_elvec_u.vec()[iloc + iComp * M_uFESpace.fe().nbNode] = velVecRep[ig];
                }
            }

            // compute elemental convective rhs term on current element
	    source_u1gradu2(M_data.density(), M_elvec_u, M_elvec_u_aux, 
	                    M_elemrhs_u, M_uFESpace.fe());

            // assemble convective term into u system rhs 
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                assembleVector(_u_rhs, M_elemrhs_u, M_uFESpace.fe(), 
		                M_uFESpace.dof(), iComp, iComp*velTotalDof);
            }

	    // compute viscous correction when in the Crank-Nicholson step
	    if (step == STEP_2)
	    {
	        M_elemrhs_u.zero();
		source_stiff_strain(- M_data.viscosity(), M_elvec_u, M_elemrhs_u, 
		                      M_uFESpace.fe());
	        for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                    assembleVector(_u_rhs, M_elemrhs_u, M_uFESpace.fe(), 
		                    M_uFESpace.dof(), iComp, iComp*velTotalDof);
	    }

        } // volume loop

        M_comm->Barrier();

	u_rhs += _u_rhs;

        chrono.stop();
        if (M_verbose) std::cout << "done in " << chrono.diff() << " s."
                                 << std::endl;

    // Handles  stabilization terms
   
    if (M_verbose)
        std::cout << "  f-  Adding IP stabilization (explicit) terms" << std::endl;

    M_ipStab.apply_expl(u_rhs, betaVecRep);

    } // if (normInf != 0.)


    M_updated = true;


} // updateSYS_u


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::
computeCTRHS_u(vector_type& p_vec, vector_type& u_rhs)
{

    Chrono chronoCTrhs;

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();

    if (M_verbose)
	std::cout << "  f-  Updating velocity system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start(); 

    // vector with repeated nodes over the processors
    vector_type press_sol (p_vec, Repeated);

    vector_type _rhsNoBC_u ( u_rhs, Repeated, Zero);
    _rhsNoBC_u *= 0.;

    // loop on volumes
    for (UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    { 
	M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList (iVol) );
        M_pFESpace.fe().updateFirstDeriv( M_pFESpace.mesh()->volumeList (iVol) );
	
	// zero for elemental rhs in velocity/pressure
	M_elemrhs_u.zero();
	UInt eleID_p = M_pFESpace.fe().currentLocalId();
	M_elvec_p.zero();
	
	// get local pressure vector 
	for ( UInt iNode = 0; iNode < ( UInt ) M_pFESpace.fe().nbNode; iNode++ )
	{ 
	    UInt iloc = M_pFESpace.fe().patternFirst( iNode );
	    UInt ig = M_pFESpace.dof().localToGlobal( eleID_p, iloc+1 );
	    M_elvec_p.vec()[ iloc ] = press_sol( ig );
	} 

        // case where variational pressure applies
	// get local ( p^n | div v ) term (wrt velocity system) 
        if (M_hasVariationalPressure) {
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp)
            {
	        source_pdivv( 1.0, M_elvec_p, M_elemrhs_u,
		                M_pFESpace.fe(), M_uFESpace.fe(), iComp);
		assembleVector( _rhsNoBC_u, M_elemrhs_u, M_uFESpace.fe(),
		                    M_uFESpace.dof(), iComp, iComp*velTotalDof);
	    }
	} else {
        // otherwise standard Chorin-Temam coupling
	// get local ( grad p^n | v ) rhs term (wrt velocity system)
	    for ( UInt iComp = 0; iComp < nDimensions; ++iComp) 
	    {
	        source_gradpv( -1.0, M_elvec_p, M_elemrhs_u, 
				    M_pFESpace.fe(), M_uFESpace.fe(), iComp);

	        assembleVector( _rhsNoBC_u, M_elemrhs_u, M_uFESpace.fe(), 
				    M_uFESpace.dof(), iComp, iComp*velTotalDof); 
	    } 
        }	
    }

    // be sure assembling has completed
    M_comm->Barrier();
 
    u_rhs += _rhsNoBC_u;

    // note: communication, i.e. GlobalAssembling, will be done after.
    chronoCTrhs.stop();

    if (M_verbose)
	std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_u

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::computeCTRHS_p(vector_type& u_vec, 
                                                     vector_type& p_rhs, 
						     const Real dt)
{
    Chrono chronoCTrhs;

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();

    if (M_verbose)
	std::cout << "  f-  Updating pressure/projector system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start(); 

    vector_type vel_sol( u_vec, Repeated );

    vector_type _rhsNoBC_p ( p_rhs, Repeated, Zero);
    _rhsNoBC_p *= 0.;
    
    // loop on volumes
    for (UInt iVol = 1; iVol <= M_pFESpace.mesh()->numVolumes(); iVol++ )
    { 
        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList (iVol) );
        M_pFESpace.fe().updateFirstDeriv( M_pFESpace.mesh()->volumeList (iVol) );
	
	// zero for elemental rhs in velocity/pressure
        M_elemrhs_p.zero();
	UInt eleID_u = M_uFESpace.fe().currentLocalId();
	M_elvec_u.zero();
	
	// get local velocity vector 
   	for ( UInt iNode = 0; iNode < ( UInt ) M_uFESpace.fe().nbNode; iNode++ )
	{
	    UInt iloc = M_uFESpace.fe().patternFirst( iNode );
	    
	    for ( UInt iComp = 0; iComp < nDimensions; ++iComp)
	    { 
		  UInt ig = M_uFESpace.dof().localToGlobal( eleID_u, iloc+1 ) + 
					iComp*dim_u();
		  M_elvec_u.vec()[ iloc+iComp*M_uFESpace.fe().nbNode  ] = vel_sol( ig );
	    }
	}
	
	// get local -\rho/dt ( div \tilde u^{n+1} | q ) rhs term 
	// (wrt pressure system)
	
	source_divuq( - M_data.density() / dt, M_elvec_u, 
                       M_elemrhs_p, M_uFESpace.fe(), M_pFESpace.fe() );
	assembleVector ( _rhsNoBC_p, M_elemrhs_p, M_pFESpace.fe(), 
			M_pFESpace.dof(),0);

    }

    // be sure assembling has completed
    M_comm->Barrier();
 
    p_rhs += _rhsNoBC_p;

    // note: all communication should be done later  

    chronoCTrhs.stop();

    if (M_verbose)
	std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_p


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::iterate_u(bchandler_raw_type& bch_u, step_type step)
{

    Chrono chrono;

    // update convective term
    if (step == STEP_1)
        updateSYS_u( M_sol_u_prev, M_rhsNoBC_u, step);
    if (step == STEP_2) {
        vector_type _sol_u_tmp (M_sol_u_prev, Repeated);
        _sol_u_tmp *= 0.0;
        _sol_u_tmp += M_sol_u_prev;
        _sol_u_tmp += M_sol_u_aux;
        _sol_u_tmp *= 0.5;
	updateSYS_u( _sol_u_tmp, M_rhsNoBC_u, step);
    }

    // matrix and vector assembling communication
    if (M_verbose)
        std::cout << "  f-  Finalizing the velocity matrix and vectors ...    ";

    chrono.start();

    M_matrNoBC_u->GlobalAssemble();
    M_rhsNoBC_u.GlobalAssemble();

    matrix_ptrtype matrFull_u( new matrix_type( M_localMap_u, M_matrNoBC_u->getMeanNumEntries()));

    *matrFull_u += *M_matrNoBC_u;

    vector_type    rhsFull_u = M_rhsNoBC_u;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying velocity boundary conditions ...         "
              << std::flush;

    chrono.start();

    applyBC_u( *matrFull_u, rhsFull_u, bch_u);

    matrFull_u->GlobalAssemble();

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;
    

    // solve velocity system and, if second step, update pressure system rhs with newly 
    // computed velocity (Chorin-Temam coupling)

    if (step == STEP_1)
    {
        solveSystem_u(matrFull_u, rhsFull_u, M_sol_u_aux);
	M_residual_u  = M_rhsNoBC_u;
	M_residual_u -= *M_matrNoBC_u * M_sol_u_aux;
    }

    if (step == STEP_2)
    {
        solveSystem_u(matrFull_u, rhsFull_u, M_sol_u);
	M_residual_u  = M_rhsNoBC_u;
	M_residual_u -= *M_matrNoBC_u * M_sol_u;
	computeCTRHS_p(M_sol_u, M_rhsNoBC_p, M_data.dataTime()->getTimeStep());
    }

} // iterate_u()

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::iterate_p(bchandler_raw_type& bch_p)
{

    Chrono chrono;

    // matrix and vector assembling communication

    if (M_verbose)
        {
            std::cout << "  f-  Finalizing the pressure matrix and vectors ...    ";
        }

    chrono.start();

    M_matrNoBC_p->GlobalAssemble();
    M_rhsNoBC_p.GlobalAssemble();

    matrix_ptrtype matrFull_p( new matrix_type( M_localMap_p, M_matrNoBC_p->getMeanNumEntries()));
    
    *matrFull_p += *M_matrNoBC_p;

    vector_type    rhsFull_p = M_rhsNoBC_p;

    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // boundary conditions update
    M_comm->Barrier();
    if (M_verbose) std::cout << "  f-  Applying pressure boundary conditions ...         "
              << std::flush;

    chrono.start();

    applyBC_p(*matrFull_p, rhsFull_p, bch_p);

    matrFull_p->GlobalAssemble();

    chrono.stop();

    M_comm->Barrier();

    if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n" << std::flush;

    solveSystem_p(matrFull_p, rhsFull_p, M_sol_p);
    M_residual_p  = M_rhsNoBC_p;
    M_residual_p -= *M_matrNoBC_p * M_sol_p;
 
} // iterate_p()

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::time_advance(Real const& time, step_type step)
{

    // what to do at first time step
    if (M_firstTimeStep) 
    {
	// note: must call initialize() outside the solver for proper init
        buildSystem_u_p();
        resetPrec_u();				
        resetPrec_p();
	M_firstTimeStep = 0;
    }
    
    M_rhsNoBC_u *= 0.;
    M_rhsNoBC_p *= 0.;
    
    // note: the two steps have same pressure/mass force update as it is not 
    //       pure RK2.

// SDEB
//    M_sol_p *= 0.;
// EDEB
    // compute pressure force
    computeCTRHS_u(M_sol_p, M_rhsNoBC_u);
    
    
    // add mass term (here time step is also dt: no pure RK)
    M_rhsNoBC_u += *M_matrMass * M_sol_u;

    // update previous velocity
    if (step == STEP_2)
        M_sol_u_prev = M_sol_u;

    // note: velocity convection updating is done in iterate_u method
 
} // time_advance


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::solveSystem_u( matrix_ptrtype  matrFull,
                                                     vector_type&    rhsFull,
						     vector_type&    solFull)
{
    Chrono chrono;

    if (M_verbose)
        std::cout << "  f-  Setting up the velocity solver ...       ";

    chrono.start();
    // assert (M_matrFull.get() != 0);
    M_linearSolver_u.setMatrix(*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // overlapping schwarz preconditioner

    if ( !M_reusePrec_u || M_resetPrec_u || !M_prec_u->set() )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f-  Computing the velocity precond ...       ";

        M_prec_u->buildPreconditioner(matrFull);

        double condest = M_prec_u->Condest();

        M_linearSolver_u.setPreconditioner(M_prec_u);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-       Estimated condition number = " << condest << "\n" <<  std::flush;
        }


    }
    else
    {
        if (M_verbose)
            std::cout << "  f-  Reusing  velocity precond ...       \n" <<  std::flush;
    }


    chrono.start();

    if (M_verbose)
        std::cout << "  f-  Solving velocity system ...                       ";

    int numIter = M_linearSolver_u.solve(solFull, rhsFull);

    if (numIter >= M_maxIterSolver_u)
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f- Iterative solver failed, recomputing the velocity precond ...       ";

        M_prec_u->buildPreconditioner(matrFull);

        double condest = M_prec_u->Condest();

        M_linearSolver_u.setPreconditioner(M_prec_u);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-      Estimated condition number = " << condest << "\n" <<  std::flush;
        }

        numIter = M_linearSolver_u.solve(solFull, rhsFull);

        if (numIter >= M_maxIterSolver_u) {
            if (M_verbose) { 
              std::cout << "  f- ERROR: Iterative velocity solver failed again.\n";
              std::cout << "  fx-  Exiting ..." << std::endl;
            }
            exit(1);
        }

    }

    M_resetPrec_u = (numIter >= M_maxIterForReuse_u);

    chrono.stop();
    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }

    M_comm->Barrier();

} // solvesystem_u

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::solveSystem_p( matrix_ptrtype  matrFull,
                                                     vector_type&    rhsFull, 
						     vector_type&    solFull)
{
    Chrono chrono;

    if (M_verbose)
        std::cout << "  f-  Setting up the pressure solver ...       ";

    chrono.start();
    // assert (M_matrFull.get() != 0);
    M_linearSolver_p.setMatrix(*matrFull);
    chrono.stop();

    if (M_verbose)
        std::cout << "done in " << chrono.diff() << " s.\n"
                  << std::flush;

    // overlapping schwarz preconditioner

    if ( !M_reusePrec_p || M_resetPrec_p || !M_prec_p->set() )
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f-  Computing the pressure precond ...       ";

        M_prec_p->buildPreconditioner(matrFull);

        double condest = M_prec_p->Condest();

        M_linearSolver_p.setPreconditioner(M_prec_p);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-       Estimated condition number = " << condest << "\n" <<
std::flush;
        }


    }
    else
    {

        if (M_verbose)
            std::cout << "  f-  Reusing  pressure precond ...       \n" <<  std::flush;
    }


    chrono.start();

    if (M_verbose)
        std::cout << "  f-  Solving pressure system ...                       ";

    int numIter = M_linearSolver_p.solve(solFull, rhsFull);

    if (numIter >= M_maxIterSolver_p)
    {
        chrono.start();

        if (M_verbose)
            std::cout << "  f- Iterative pressure solver failed, recomputing the precond ...";

        M_prec_p->buildPreconditioner(matrFull);

        double condest = M_prec_p->Condest();

        M_linearSolver_p.setPreconditioner(M_prec_p);

        chrono.stop();
        if (M_verbose)
        {
            std::cout << "done in " << chrono.diff() << " s.\n";
            std::cout << "  f-      Estimated condition number = " << condest << "\n" << std::flush;
        }

        numIter = M_linearSolver_p.solve(solFull, rhsFull);

        if (numIter >= M_maxIterSolver_p) {
            if (M_verbose) {
              std::cout << "  f- ERROR: Iterative pressure solver failed again.\n";
              std::cout << "  fx-  Exiting ..." << std::endl;
            }
            exit(1);
        }

    }

    M_resetPrec_p = (numIter >= M_maxIterForReuse_p);

    chrono.stop();
    if (M_verbose)
    {
        std::cout << "\ndone in " << chrono.diff()
                  << " s. ( " << numIter << "  iterations. ) \n"
                  << std::flush;
    }

    M_comm->Barrier();

} // solveSystem_p


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::removeMean( vector_type& x, UInt comp )
{
    Real sum1 = 0.;
    Real sum0 = 0.;

    for ( UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    {
        M_pFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
        sum1 += elem_integral( x, M_pFESpace.fe(), M_pFESpace.dof(), comp );
        sum0 += M_pFESpace.fe().measure();
    }
    Real mean = sum1/sum0;

    std::cout << "   Mean pressure value : " << mean << std::endl;
    for ( UInt iNode = (comp - 1)*dim_u(); iNode < comp*dim_u(); ++iNode )
    {
        x[ iNode + 1 ] -= mean; // BASEINDEX + 1
    }

} // removeMean()


template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::applyBC_u( matrix_type&        matrix_u,
                                               vector_type&        rhs_u,
                                               bchandler_raw_type& BCh_u)
{
    if ( !BCh_u.bdUpdateDone() )
    {
        BCh_u.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    vector_type rhsFull_u (rhs_u, Repeated, Zero);

    bcManage( matrix_u, rhsFull_u, *M_uFESpace.mesh(), M_uFESpace.dof(), 
		BCh_u, M_uFESpace.feBd(), 1., M_data.dataTime()->getTime() );

    rhs_u = rhsFull_u;

    if ( BCh_u.hasOnlyEssential() && M_diagonalize_u )
    {
        matrix_u.diagonalize( nDimensions*dim_u(),
                            M_diagonalize_u,
                            rhs_u,
                            0.);
    }
} // applyBC_u

template<typename Mesh, typename SolverType>
void ChorinTemamRK<Mesh, SolverType>::applyBC_p( matrix_type&	   matrix_p,
					       vector_type&	   rhs_p,
					       bchandler_raw_type& BCh_p )
{
    if ( !BCh_p.bdUpdateDone() )
    {
	BCh_p.bdUpdate( *M_pFESpace.mesh(), M_pFESpace.feBd(), M_pFESpace.dof() );
    }

    vector_type rhsFull_p (rhs_p, Repeated, Zero);

    bcManage( matrix_p, rhsFull_p, *M_pFESpace.mesh(), M_pFESpace.dof(), 
		BCh_p, M_pFESpace.feBd(), 1., M_data.dataTime()->getTime() );

    rhs_p = rhsFull_p;

    if ( BCh_p.hasOnlyEssential() && M_diagonalize_p )
    {
	matrix_p.diagonalize( dim_p(),
			    M_diagonalize_p,
			    rhs_p,
			    0.);
    }

} // applyBC_p

// Postprocessing
template <typename Mesh, typename SolverType>
void
ChorinTemamRK<Mesh, SolverType>::postProcess(bool _writeMesh)
{
	// maybe we will add something useful here
}



} // namespace LifeV


#endif //_CHORIN_TEMAM_RK_H_
