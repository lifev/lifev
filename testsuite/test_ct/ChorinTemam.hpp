/*
 * \file ChorinTemam.hpp
 * \desc Chorim-Temam method w/, at the moment, first order Euler backward
 */

#ifndef _CHORIN_TEMAM_H_
#define _CHORIN_TEMAM_H_

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
//#include <life/lifefem/sobolevNorms.hpp>
#include <life/lifefem/geoMap.hpp>
//#include <life/lifesolver/nsipterms.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>

#include <boost/shared_ptr.hpp>

#include "life/lifefem/FESpace.hpp"
//#include <life/lifefem/bdf_template.hpp>

namespace LifeV
{

template< typename Mesh,
          typename SolverType = LifeV::Epetra::SolverTrilinos >
class ChorinTemam

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

    //! Constructor
    /*!
      \param dataType
      \param localMap
      \param velocity FE space
      \param pressure FE space
      \param bcHu boundary conditions for the velocity
      \param bcHp boundary conditions for the pressure
    */
    ChorinTemam( const data_type&          dataType,
           FESpace<Mesh, EpetraMap>& uFESpace,
           FESpace<Mesh, EpetraMap>& pFESpace,
           BCHandler&                bcHu,
	   BCHandler& 		     bcHp,
           Epetra_Comm&              comm );

    //! virtual destructor
    virtual ~ChorinTemam();

    //! update u matrix, apply u bc's, solve u system, compute p system rhs
    virtual void iterate_u (bchandler_raw_type& bch_u);

    //! apply p bc's, solve p system
    virtual void iterate_p (bchandler_raw_type& bch_p);

    //! update rhs for u system (i.e. < \nabla p^n | v>) and add M*u^n
    void time_advance (const Real& time);

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

    //TODO: we do not need these methods here
    //! reduce the local solution solution in global vectors
    void reduceSolution( Vector& u,
                         Vector& p );

    void reduceResidual( Vector& res);


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
    // update convective term for u system
    void updateSystem_u (vector_type& betaVec);
    // fill the p system matrix with pressure matrix (!)
    void updateSystem_p ();
    // compute u system CT coupling term - < \nabla p^n  v >
    void computeCTRHS_u (vector_type& p_sol);
    // compute p system CT coupling term -\rho / dt < \nabla \cdot u^{n+1} | q >
    void computeCTRHS_p (vector_type& u_sol);    

    void solveSystem_u            ( matrix_ptrtype matrFull_u,
                                    vector_type&   rhsFull_u );

    void solveSystem_p 		  ( matrix_ptrtype matrFull_p,
				    vector_type&   rhsFull_p );

    void applyBC_u (matrix_type& matrix_u,
                    vector_type& rhs_u,
                    bchandler_raw_type& BCh_u);
    void applyBC_p (matrix_type& matrix_p,
                    vector_type& rhs_p,
                    bchandler_raw_type& BCh_p );
    
    void echo(std::string message);

    //! removes mean of component comp of vector x
    void removeMean( vector_type& x, UInt comp = 1 );

    //private members

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
    matrix_ptrtype                   M_matrMass;

    //! Stokes matrix: mu*stiff
    matrix_ptrtype                   M_matrStokes;

    //! pressure matrix 
    matrix_ptrtype		    M_matrPress;

    //! matrix without boundary conditions

    matrix_ptrtype                 M_matrNoBC_u;
    matrix_ptrtype		   M_matrNoBC_p;
    
    //! source term for NS
    source_type                    M_source;


    //! Right hand side for the velocity
    vector_type                    M_rhsNoBC_u;
    vector_type			   M_rhsNoBC_p;

    //! Solutions _u and _p
    vector_type                    M_sol_u;
    vector_type			   M_sol_p;

    // solution u at previous time
    vector_type			   M_sol_u_prev;

    //! residual

    vector_type                    M_residual_u;
    vector_type			   M_residual_p;

    SolverType                     M_linearSolver_u;
    SolverType			   M_linearSolver_p;

    prec_type                      M_prec_u;
    prec_type			   M_prec_p;

    bool                           M_steady;

// TODO: sd stabilization for convective effects in u system

    //! Stabilization
    bool                           M_stab;

//    details::
//    IPStabilization<Mesh, Dof>     M_ipStab; 
//    Real                           M_gammaBeta;
//    Real                           M_gammaDiv;
//    Real                           M_gammaPress;

// end TODO

    const Function*                M_betaFct;

    bool                           M_divBetaUv;

    //
    double                         M_diagonalize_u;
    double 			   M_diagonalize_p;

    UInt                           M_count;

    //! boolean that indicates if output is sent to cout

    bool                           M_verbose;

    //! boolean that indicates if the matrix is updated for the current iteration

    bool                           M_updated;

    //! boolean that indicates if le precond has to be recomputed

    bool                           M_reusePrec_u;
    int                            M_maxIterForReuse_u;
    bool                           M_resetPrec_u;
    bool 			   M_reusePrec_p;
    int				   M_maxIterForReuse_p;
    bool 			   M_resetPrec_p;

    //! interger storing the max number of solver iteration with prec recomputing

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
    ElemVec			   M_elvec_p;         // Elementary pressure vector
    ElemVec			   M_elemrhs_u;
    ElemVec			   M_elemrhs_p;		 

    Real			   M_time;	       // current time
    bool			   M_firstTimeStep; 
//    BdfT<vector_type>		   M_bdf_u;	       // velocity bdf 
//    BdfT<vector_type>		   M_bdf_p;              // pressure bdf
}; // class ChorinTemam



/* -------------------------------------------------------------- */
/* IMPLEMENTATION						  */
/* -------------------------------------------------------------- */


template<typename Mesh, typename SolverType>
ChorinTemam<Mesh, SolverType>::
ChorinTemam( const data_type&          dataType,
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
    M_prec_u                 ( new prec_raw_type() ),
    M_prec_p		     ( new prec_raw_type() ),
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
    M_elvec_p		     ( M_pFESpace.fe().nbNode, 1 ),
    M_elemrhs_u		     ( M_uFESpace.fe().nbNode, nDimensions ),
    M_elemrhs_p              ( M_pFESpace.fe().nbNode, 1 ),
    M_rhsNoBC_u              ( M_localMap_u ),
    M_rhsNoBC_p		     ( M_localMap_p ),
    M_sol_u                  ( M_localMap_u ),
    M_sol_p		     ( M_localMap_p ),
    M_residual_u             ( M_localMap_u ),
    M_residual_p 	     ( M_localMap_p ), 
    M_stab                   ( false ),
/*    M_ipStab                 ( M_uFESpace.mesh(),
                               M_uFESpace.dof(), M_uFESpace.refFE(),
                               M_uFESpace.feBd(), M_uFESpace.qr(),
                               0., 0., 0.,
                               M_data.viscosity() ),
*/ // mg: tbm
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
//    M_bdf_u                  (1), 	  // first order time
//    M_bdf_p                  (1)        // first order time
    M_sol_u_prev             ( M_localMap_u )
    
{
    M_stab = (&M_uFESpace.refFE() == &M_pFESpace.refFE());
}

template<typename Mesh, typename SolverType>
ChorinTemam<Mesh, SolverType>::
~ChorinTemam()
{

}


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::setUp( const GetPot& dataFile )
{
    M_steady      = dataFile( "fluid/miscellaneous/steady",        1  );
//    M_gammaBeta   = dataFile( "fluid/ipstab/gammaBeta",            0. );
//    M_gammaDiv    = dataFile( "fluid/ipstab/gammaDiv",             0. );
//    M_gammaPress  = dataFile( "fluid/ipstab/gammaPress",           0. );
    M_divBetaUv   = dataFile( "fluid/discretization/div_beta_u_v", 0  );
    M_diagonalize_u = dataFile( "fluid/discretization/diagonalizeVel",  1. );
    M_diagonalize_p = dataFile( "fluid/discretization/diagonalizePress",  1. );

// TODO modify the following two lines for handling different vel / press solvers ... 
    M_linearSolver_u.setDataFromGetPot( dataFile, "fluid/solver" );
    M_linearSolver_p.setDataFromGetPot( dataFile, "fluid/solver" );

//    M_ipStab.setGammaBeta (M_gammaBeta);
//    M_ipStab.setGammaDiv  (M_gammaDiv);
//    M_ipStab.setGammaPress(M_gammaPress);

// TODO modify the following lines to handle different solvers for press and velocity ... 
    M_maxIterSolver_u   = dataFile( "fluid/solver/max_iter", -1);
    M_reusePrec_u       = dataFile( "fluid/prec/reuse", true);
    M_maxIterForReuse_u = dataFile( "fluid/prec/max_iter_reuse", M_maxIterSolver_u*8/10);
    M_maxIterSolver_p   = dataFile( "fluid/solver/max_iter", -1);
    M_reusePrec_p       = dataFile( "fluid/prec/reuse", true);
    M_maxIterForReuse_p = dataFile( "fluid/prec/max_iter_reuse", M_maxIterSolver_p*8/10);

    M_prec_u->setDataFromGetPot( dataFile, "fluid/prec" );
    M_prec_p->setDataFromGetPot( dataFile, "fluid/prec" );
}

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::buildSystem_u_p()
{

    // NB: we compute u and p systems constant matrices at the same time
    // one loop on geom element w/ possibly different fe for u and p

    M_matrMass.reset  ( new matrix_type(M_localMap_u) );
    M_matrStokes.reset( new matrix_type(M_localMap_u) );
    M_matrPress.reset ( new matrix_type(M_localMap_p) );

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

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
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

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
            mass( M_data.density()/M_data.timestep(), M_elmatMass, M_uFESpace.fe(), 0, 0, nDimensions );
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

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diff_cumul() << " s.\n"
                  << " Stab           " << chronoStab.diff_cumul() << " s.\n"
                  << " Zero           " << chronoZero.diff_cumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diff_cumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diff_cumul() << " s.\n"
                  << " Mass           " << chronoMass.diff_cumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diff_cumul() << " s.\n"
                  << std::endl;

}


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::
initialize( const Function& u0, const Function& p0 )
{
     vector_type u(M_uFESpace.map());
     M_uFESpace.interpolate(u0, u, 0.0);

     vector_type p(M_pFESpace.map());
     M_pFESpace.interpolate(p0, p, 0.0);

     initialize(u, p);
}


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::
initialize( const vector_type& u0, const vector_type& p0 )
{

    M_sol_u = u0;
    M_sol_p = p0;

}

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::updateSystem_u(vector_type& betaVec)
{
	
    Real alpha = 1./M_data.timestep();    
 
    /* copy matrices at the right place */     

    if (M_recomputeMatrix)
        buildSystem_u_p();

    if (M_verbose)
          std::cout << "  f-  Copying the matrices ...                 "
                    << std::flush;

    M_matrNoBC_u.reset(new matrix_type(M_localMap_u, M_matrStokes->getMeanNumEntries() ));
    M_matrNoBC_p.reset(new matrix_type(M_localMap_p, M_matrPress->getMeanNumEntries() ));
    *M_matrNoBC_u += *M_matrStokes;
    *M_matrNoBC_p += *M_matrPress;

    if (alpha != 0. )
    {
        *M_matrNoBC_u += *M_matrMass;
    }

    /* Update convective term for velocity system */

    UInt nbCompU = nDimensions;
    UInt velTotalDof = M_uFESpace.dof().numTotalDof();
    Chrono chrono;

    double normInf;
    betaVec.NormInf(&normInf);

    normInf = 0.; 	// remove this after debug
			// just to prevent convective term in the absence of stabilization

    if (normInf != 0.)
    {

        if (M_verbose)
            std::cout << "  f-  Updating the convective terms ...        "
                      << std::flush;

        // vector with repeated nodes over the processors
	vector_type betaVecRep( betaVec, Repeated );

        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
	    // reuse elmatStiff as space for storing elemtal convective matrix
            M_elmatStiff.zero();

            UInt eleID = M_uFESpace.fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < ( UInt ) M_uFESpace.fe().nbNode ; iNode++ )
            {
                UInt  iloc = M_uFESpace.fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
                {
                    UInt ig = M_uFESpace.dof().localToGlobal( eleID, iloc + 1 ) + iComp*dim_u();
                    M_elvec_u.vec()[ iloc + iComp*M_uFESpace.fe().nbNode ] = M_data.density() * betaVecRep[ig]; // BASEINDEX + 1
                }
            }

            // Stabilising term: div u^n u v
            if ( M_divBetaUv )
            {
                mass_divw( 0.5*M_data.density(), M_elvec_u, M_elmatStiff, M_uFESpace.fe(), 0, 0, nbCompU );
            }

            // loop on components
            for ( UInt iComp = 0; iComp < nbCompU; ++iComp )
            {
                // compute local convective term and assembling
                grad( 0, M_elvec_u, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                grad( 1, M_elvec_u, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );
                grad( 2, M_elvec_u, M_elmatStiff, M_uFESpace.fe(), M_uFESpace.fe(), iComp, iComp );

                assembleMatrix( *M_matrNoBC_u,
                                M_elmatStiff,
                                M_uFESpace.fe(),
                                M_uFESpace.fe(),
                                M_uFESpace.dof(),
                                M_uFESpace.dof(),
                                iComp, iComp,
                                iComp*velTotalDof, iComp*velTotalDof
                                );

            }

        }


        chrono.stop();
        if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                                 << std::flush;

// TODO Modify what follows for handling sd stab
/* 
        if (M_stab)
        {  // if there is a problem of indeces here, you should use the repeated version of beta
            //M_ipStab.applyPressure( *M_matrFull, betaVecRep, M_verbose );
            M_ipStab.applyVelocity( *M_matrNoBC, betaVecRep, M_verbose );
        }
*/
    }

    M_updated = true;

    if (alpha != 0.)
    {
        M_matrNoBC_u->GlobalAssemble();
    }

} // updateSystem_u


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::
computeCTRHS_u(vector_type& p_sol)
{

    Chrono chronoCTrhs;

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();

//    M_updated = false;
    
    if (M_verbose)
	std::cout << "  f-  Updating velocity system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start(); 

    // vector with repeated nodes over the processors
    vector_type press_sol (p_sol, Repeated);

//    vector_type _rhsNoBC_u (M_rhsNoBC_u, Repeated);

    M_rhsNoBC_u *= 0.0;
//    _rhsNoBC_u *= 0.0;			// try this for initialization
    vector_type _rhsNoBC_u ( M_rhsNoBC_u, Repeated, Zero);

    // loop on volumes
    for (UInt iVol = 1; iVol <= M_uFESpace.mesh()->numVolumes(); iVol++ )
    { 
	M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList (iVol) );
        M_pFESpace.fe().updateFirstDeriv( M_pFESpace.mesh()->volumeList (iVol) );
	
	// raz for elemental rhs in velocity/pressure
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

	// get local ( grad p^n | v ) rhs term (wrt velocity system)
	for ( UInt iComp = 0; iComp < nDimensions; ++iComp) 
	{
	    source_gradpv( -1.0, M_elvec_p, M_elemrhs_u, 
				M_pFESpace.fe(), M_uFESpace.fe(), iComp);

	    assembleVector( _rhsNoBC_u, M_elemrhs_u, M_uFESpace.fe(), 
				M_uFESpace.dof(), iComp, iComp*velTotalDof); 

	} 
	
    }

    // be sure assembling has completed
    M_comm->Barrier();
 
    M_rhsNoBC_u = _rhsNoBC_u; 	// we mimic OseenShapeDerivative here, no += 

    //M_rhsNoBC_u.GlobalAssemble();	// these 2 lines have no effect, and after all
					// comm should be done after 

    chronoCTrhs.stop();

    if (M_verbose)
	std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_u

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::computeCTRHS_p(vector_type& u_sol)
{
    Chrono chronoCTrhs;

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();
//    UInt pressTotalDof = M_pFESpace.dof().numTotalDof();

    if (M_verbose)
	std::cout << "  f-  Updating pressure system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start(); 

    vector_type vel_sol( u_sol, Repeated );

//    vector_type _rhsNoBC_p (M_rhsNoBC_p, Repeated);

    M_rhsNoBC_p *= 0.0;			// for test pupose
//    _rhsNoBC_p *= 0.0;			// try this for initialization
    vector_type _rhsNoBC_p ( M_rhsNoBC_p, Repeated, Zero); 
    
    // loop on volumes
    for (UInt iVol = 1; iVol <= M_pFESpace.mesh()->numVolumes(); iVol++ )
    { 
        M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList (iVol) );
        M_pFESpace.fe().updateFirstDeriv( M_pFESpace.mesh()->volumeList (iVol) );
	
	// raz for elemental rhs in velocity/pressure
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
	
	source_divuq( - M_data.density() / M_data.timestep(), M_elvec_u, M_elemrhs_p,
				M_uFESpace.fe(), M_pFESpace.fe() );
	assembleVector ( _rhsNoBC_p, M_elemrhs_p, M_pFESpace.fe(), 
			M_pFESpace.dof(),0);

    }

    // be sure assembling has completed
    M_comm->Barrier();
 
    M_rhsNoBC_p = _rhsNoBC_p;		// taken from Oseen shape derivative

    //M_rhsNoBC_p.GlobalAssemble();	// Z: all communication should be done later  

    chronoCTrhs.stop();

    if (M_verbose)
	std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_p

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::iterate_u( bchandler_raw_type& bch_u)
{

    Chrono chrono;

    // update convective terms
//    vector_type beta(M_sol_u, Repeated);
//    updateSystem_u(beta);
    updateSystem_u(M_sol_u_prev);		// no pb for M_un init here since
						// time_advance is called before

    // matrix and vector assembling communication

    if (M_verbose)
        {
            std::cout << "  f-  Finalizing the velocity matrix and vectors ...    ";
        }

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
    // solving the velocity system

    solveSystem_u( matrFull_u, rhsFull_u );

    M_residual_u  = M_rhsNoBC_u;
    M_residual_u -= *M_matrNoBC_u*M_sol_u;

    // update RHS pressure term with new computed velocity 
    computeCTRHS_p(M_sol_u);

} // iterate_u()

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::iterate_p(bchandler_raw_type& bch_p )
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

    solveSystem_p( matrFull_p, rhsFull_p );

    M_residual_p  = M_rhsNoBC_p;
    M_residual_p -= *M_matrNoBC_p*M_sol_p; 

} // iterate_p()

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::time_advance(Real const& time)
{

    Real dt = M_data.timestep();

    // what to do at first time step
    if (M_firstTimeStep) 
    {
	//XXX: must call initialize in the test_xxx.cpp for proper init
        buildSystem_u_p();			// build contant matrices
        resetPrec_u();				
        resetPrec_p();
	M_firstTimeStep = 0;
    }

    // update velocity system rhs
    
    // first update chorin-temam coupling term 
    computeCTRHS_u(M_sol_p);
    
    // then add the mass term rho/dt * u^n
    if (M_verbose) 
	std::cout << "  f-  Adding velocity mass term on rhs ..." << std::endl;
    
    M_rhsNoBC_u += *M_matrMass * M_sol_u;

    // update old velocity
    M_sol_u_prev = M_sol_u;
  
} // time_advance



template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::solveSystem_u( matrix_ptrtype  matrFull,
                                           vector_type&    rhsFull )
{
    Chrono chrono;

    if (M_verbose)
        std::cout << "  f-  Setting up the velocity solver ...       ";

    chrono.start();
//    assert (M_matrFull.get() != 0);
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

    int numIter = M_linearSolver_u.solve(M_sol_u, rhsFull);

    if (numIter > M_maxIterSolver_u)
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

        numIter = M_linearSolver_u.solve(M_sol_u, rhsFull);

        if (numIter > M_maxIterSolver_u && M_verbose)
            std::cout << "  f- ERROR: Iterative velocity solver failed again.\n";

    }

    M_resetPrec_u = (numIter > M_maxIterForReuse_u);

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
void ChorinTemam<Mesh, SolverType>::solveSystem_p( matrix_ptrtype  matrFull,
                                           vector_type&    rhsFull )
{
    Chrono chrono;

    if (M_verbose)
        std::cout << "  f-  Setting up the pressure solver ...       ";

    chrono.start();
//    assert (M_matrFull.get() != 0);
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

    int numIter = M_linearSolver_p.solve(M_sol_p, rhsFull);

    if (numIter > M_maxIterSolver_p)
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

        numIter = M_linearSolver_p.solve(M_sol_p, rhsFull);

        if (numIter > M_maxIterSolver_p && M_verbose)
            std::cout << "  f- ERROR: Iterative pressure solver failed again.\n";

    }

    M_resetPrec_p = (numIter > M_maxIterForReuse_p);

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
void ChorinTemam<Mesh, SolverType>::reduceSolution( Vector& u,
                                              Vector& p )
{
    // XXX: we keep this function but we donot need it, to be del
    
    vector_type vel(M_sol_u, 0);
    vector_type press(M_sol_p, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
        {
            u[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

        for ( UInt iDof = 0; iDof<dim_p(); ++iDof )
        {
//            p[ iDof ] = vel[ iDof + nDimensions*dim_u() + 1 ]; // BASEINDEX + 1
	    p[ iDof ] = press[ iDof + 1 ];
        }
    }

}

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::reduceResidual( Vector& res )
{
    // same rk as in function reduceSolution above ...
    vector_type vel(M_residual_u, 0);

    if (M_verbose)
    {
        for ( UInt iDof = 0; iDof < nDimensions*dim_u(); ++iDof )
        {
            res[ iDof ] = vel[ iDof + 1 ]; // BASEINDEX + 1
        }

    }
}


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::removeMean( vector_type& x, UInt comp )
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
void ChorinTemam<Mesh, SolverType>::applyBC_u( matrix_type&        matrix_u,
                                               vector_type&        rhs_u,
                                               bchandler_raw_type& BCh_u)
{
    if ( !BCh_u.bdUpdateDone() )
    {
        BCh_u.bdUpdate( *M_uFESpace.mesh(), M_uFESpace.feBd(), M_uFESpace.dof() );
    }

    vector_type rhsFull_u (rhs_u, Repeated, Zero);

    bcManage( matrix_u, rhsFull_u, *M_uFESpace.mesh(), M_uFESpace.dof(), 
		BCh_u, M_uFESpace.feBd(), 1., M_data.time() );

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
void ChorinTemam<Mesh, SolverType>::applyBC_p( matrix_type&	   matrix_p,
					       vector_type&	   rhs_p,
					       bchandler_raw_type& BCh_p )
{
    if ( !BCh_p.bdUpdateDone() )
    {
	BCh_p.bdUpdate( *M_pFESpace.mesh(), M_pFESpace.feBd(), M_pFESpace.dof() );
    }

    vector_type rhsFull_p (rhs_p, Repeated, Zero);

    bcManage( matrix_p, rhsFull_p, *M_pFESpace.mesh(), M_pFESpace.dof(), 
		BCh_p, M_pFESpace.feBd(), 1., M_data.time() ); 

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
ChorinTemam<Mesh, SolverType>::postProcess(bool _writeMesh)
{
    std::ostringstream index;
    std::ostringstream indexMe;
    std::string name;
    std::string me;

    M_count++;

    indexMe << M_me;

    switch ( indexMe.str().size() )
    {
        case 1:
            me = "00" + indexMe.str();
            break;
        case 2:
            me = "0" + indexMe.str();
            break;
        case 3:
            me = indexMe.str();
            break;
    }


//     if (_writeMesh || (M_count / M_data.verbose() == 0) )
//         writeMesh  ("partedMesh." + me + ".mesh", M_pFESpace.mesh() );

//    vector_type velAndPressure(M_sol,*M_localMap.getRepeatedEpetra_Map()); //mg
//    vector_type res(M_residual,*M_localMap.getRepeatedEpetra_Map()); //mg
    vector_type vel( M_sol_u, Repeated );
    vector_type press( M_sol_p, Repeated );
    vector_type res (M_residual_u, Repeated ); 

//         if ( fmod( float( M_count ), float( M_data.verbose() ) ) == 0.0 )
//         {
    if (M_me == 0)
        std::cout << "  F-  Post-processing " << std::flush;

    index << std::setfill('0') << std::setw(3);
    index << ( M_count / M_data.verbose() );
    name = index.str();

    PhysVectUnknown<Vector> u(nDimensions*dim_u());
    ScalUnknown<Vector>     p(dim_p());

    reduceSolution(u, p);

    if (M_me == 0)
    {
        // postprocess data file for medit
//         wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.giveVec(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.giveVec() + this->dim_u(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.giveVec()+2*this->dim_u(),
//                                M_data.mesh()->numGlobalVertices() );
//         wr_medit_ascii_scalar( "press." + name + ".bb", p.giveVec(),
//                                p.size() );

        double dt = M_data.timestep();


       writeMesh("vel_x." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_x." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("vel_y." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_y." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("vel_z." + me + "." + name + ".mesh", *M_uFESpace.mesh());
       writeMesh("resf_z." + me +  "." + name + ".mesh", *M_uFESpace.mesh());

       writeMesh("press." + me +  "." + name + ".mesh", *M_uFESpace.mesh());


        meditSolutionWriter( "vel_x." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), vel, M_uFESpace.dof().numTotalDof()*0);
        meditSolutionWriter( "vel_y." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), vel, M_uFESpace.dof().numTotalDof()*1);
        meditSolutionWriter( "vel_z." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), vel, M_uFESpace.dof().numTotalDof()*2);

//         wr_medit_ascii2("vel_x." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("vel_y." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("vel_z." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );


//         wr_medit_ascii2("resf_x." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("resf_y." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );
//         wr_medit_ascii2("resf_z." + me + "." + name + ".mesh",
//                         *M_uFESpace.mesh(), disp, M_data.factor() );

        meditSolutionWriter( "resf_x." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*0);
        meditSolutionWriter( "resf_y." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*1);
        meditSolutionWriter( "resf_z." + me + "." + name + ".bb",
                             *M_uFESpace.mesh(), res, M_uFESpace.dof().numTotalDof()*2);


        meditSolutionWriter( "press." + me + "." + name + ".bb",
                             *M_pFESpace.mesh(), press, M_uFESpace.dof().numTotalDof()*3);


//    cout << ( "should do ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//             " press." + name + "." + me + ".mesh" ).data()  ;

//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " press." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_x." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_y." + name + ".mesh" ).data() );
//         system( ( "ln -s -f " + M_data.meshDir() + M_data.meshFile() +
//                   " vel_z." + name + ".mesh" ).data() );

    }
//    }



}



} // namespace LifeV


#endif //_CHORIN_TEMAM_H_
