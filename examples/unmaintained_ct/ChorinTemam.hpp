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
 * \file ChorinTemam.hpp
 * \desc Chorin-Temam methods including pressure-correction schemes.
 * \see An Overview of projection methods for incompressible flows, Guermond & al,
 *      Inter. J. Numer. Methods Eng.
 */

#ifndef _CHORIN_TEMAM_H_
#define _CHORIN_TEMAM_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
//#include <life/lifefem/values.hpp>
//#include <life/lifearray/pattern.hpp>
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
#include <life/lifefem/bdf_template.hpp>

#include <sdStabilization.hpp>
#include <ipStabilization.hpp>
#include <elemOperCT.hpp>

namespace LifeV
{

template< typename Mesh,
typename SolverType = LifeV::SolverTrilinos >
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
    typedef enum {SD_STAB, IP_STAB_EXPL, IP_STAB_IMPL}   stab_type;

    //! Constructor
    /*!
      \param dataType
      \param velocity FE space
      \param pressure FE space
      \param bdf order for the velocity ( >=1 and <=3 as is for bdf_template)
      \param bdf order for the pressure ( >=0 and <=3 as 0 is standard Chorin-Temam case)
      \param bcHu boundary conditions for the velocity
      \param bcHp boundary conditions for the pressure
    */
    ChorinTemam( const data_type&          dataType,
                 FESpace<Mesh, EpetraMap>& uFESpace,
                 FESpace<Mesh, EpetraMap>& pFESpace,
                 int 			     uBdfOrder,
                 int			     pBdfOrder,
                 BCHandler&                bcHu,
                 BCHandler& 		     bcHp,
                 Epetra_Comm&              comm );

    //! virtual destructor
    virtual ~ChorinTemam();

    //! update u matrix, apply u bc's, solve u system, compute p system rhs
    virtual void iterate_u (bchandler_raw_type& bch_u);

    //! apply p bc's, solve p system
    virtual void iterate_p (bchandler_raw_type& bch_p);

    //! update rhs for u system i.e. Chorin-Temam coupling and mass term
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

    //! Resistance boundary conditions
    /*!
     * \param resFlag : boundary reference to which apply resistance
     * \param resVal  : resistance values
     * \param explRES : 1 for explicit resistance, 0 for implicit
     */
    void setRES(const std::vector<Real>& resVal, const std::vector<EntityFlag>& resFlag,
                const int explRES);

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

    void recomputeMatrix(bool const recomp) {M_recomputeMatrix = recomp;}

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
    // compute u system CT coupling term - < \nabla p^n  v >
    void computeCTRHS_u (vector_type& p_sol);
    // compute p system CT coupling term -\rho / dt < \nabla \cdot u^{n+1} | q >
    void computeCTRHS_p (vector_type& u_sol);
    // compute explicit resistance (rhs)
    void computeRES_expl (Real resVal, EntityFlag resFlag,
                          vector_type& beta, vector_type& rhs);
    // compute implicit resistance (lhs)
    void computeRES_impl (Real resVal, EntityFlag resFlag,
                          matrix_type& matrix);

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
    // apply all implicit resistance BC (rhs)
    void applyRES_impl (matrix_type& matrix_u);
    // apply all explicit resistance BC (lhs)
    void applyRES_expl (vector_type& rhs_u );

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

    //! Solutions _u and _p, and projector field
    vector_type                    M_sol_u;
    vector_type			   M_sol_p;
    vector_type			   M_sol_p_phi;

    //! Residuals
    vector_type                    M_residual_u;
    vector_type			   M_residual_p;

    SolverType                     M_linearSolver_u;
    SolverType			   M_linearSolver_p;

    boost::shared_ptr<EpetraPreconditioner>             M_prec_u;
    boost::shared_ptr<EpetraPreconditioner>		M_prec_p;

    bool                           M_steady;

    //! Stabilization
    bool                           M_stab;
    stab_type                      M_stabType;
    //! SD Stabilization
    boost::shared_ptr< SDStabilization<Mesh, Dof> >     M_sdStab;
    //! IP Stabilization
    boost::shared_ptr< IPStabilization<Mesh, Dof> >     M_ipStab;
    //! SD and/or IP stabilization parameters
    Real                           M_gammaBeta;
    Real                           M_gammaDiv;

    // boolean indicating resistance bc's
    bool                           M_hasRES;
    // boolean indicating explicitness or implicitness in resistance bc application
    bool			   M_hasRESexpl;
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
    ElemVec			   M_elvec_p;         // Elementary pressure vector
    ElemVec			   M_elemrhs_u;
    ElemVec			   M_elemrhs_p;
    ElemVec			   M_belvec_u;	      // useful for resistance expl
    ElemVec			   M_bele_flow_u;     // elemental flow vector
    ElemMat			   M_belmat_u;        // useful for resistance impl

    //! Time related data
    Real			   M_time;	       // current time
    bool			   M_firstTimeStep;
    int 			   M_uBdfOrder;	       // bdf order for velocity
    int 			   M_pBdfOrder;        // bdf order for pressure
    BdfT<vector_type>		   M_bdf_u;	       // velocity bdf
    BdfT<vector_type>		   M_bdf_p;            // pressure bdf
    BdfT<vector_type>              M_bdf_p_phi;        // projector field bdf

    //! Resistance related data
    std::vector<EntityFlag>	   M_resFlag;	       // resistance flags
    std::vector<Real>		   M_resVal;           // resistance values
}; // class ChorinTemam


/****************************************************************************
 *
 * Implementation
 *
 ****************************************************************************/


template<typename Mesh, typename SolverType>
ChorinTemam<Mesh, SolverType>::
ChorinTemam( const data_type&          dataType,
             FESpace<Mesh, EpetraMap>& uFESpace,
             FESpace<Mesh, EpetraMap>& pFESpace,
             int uBdfOrder,
             int pBdfOrder,
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
        M_elvec_p		     ( M_pFESpace.fe().nbNode, 1 ),
        M_elemrhs_u		     ( M_uFESpace.fe().nbNode, nDimensions ),
        M_elemrhs_p              ( M_pFESpace.fe().nbNode, 1 ),
        M_belvec_u               ( M_uFESpace.feBd().nbNode, nDimensions ),
        M_bele_flow_u            ( M_uFESpace.feBd().nbNode, nDimensions ),
        M_belmat_u               ( M_uFESpace.feBd().nbNode, nDimensions, nDimensions ),
        M_rhsNoBC_u              ( M_localMap_u ),
        M_rhsNoBC_p		     ( M_localMap_p ),
        M_sol_u                  ( M_localMap_u ),
        M_sol_p		     ( M_localMap_p ),
        M_sol_p_phi              ( M_localMap_p ),
        M_residual_u             ( M_localMap_u ),
        M_residual_p 	     ( M_localMap_p ),
        M_stab                   ( false ),
        M_stabType               ( ),
        /*
            M_sdStab                 ( M_uFESpace.mesh(),
                                       M_uFESpace.dof(),
                                       M_uFESpace.refFE(),
                                       M_uFESpace.qr(),
                                       0., 0.,
                                       M_data.viscosity()
                                     ),
        */
        M_sdStab                 ( ),
        M_ipStab                 ( ),
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
        M_uBdfOrder		     ( uBdfOrder ),
        M_pBdfOrder		     ( pBdfOrder ),
        M_bdf_u                  ( M_uBdfOrder ),
        M_bdf_p                  ( M_pBdfOrder ? M_pBdfOrder : 1),
        M_bdf_p_phi              ( M_uBdfOrder ),
        M_hasVariationalPressure ( false ),
        M_hasRES                 ( false ),
        M_hasRESexpl             ( false ),
        M_resVal                 ( 0 ),
        M_resFlag                ( 0 )
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
//    M_gammaBeta   = dataFile ( "fluid/sdstab/gammaBeta", 0. );
//    M_gammaDiv    = dataFile ( "fluid/sdstab/gammaDiv", 0.);
    M_divBetaUv   = dataFile( "fluid/space_discretization/div_beta_u_v", 0  );
    M_diagonalize_u = dataFile( "fluid/space_discretization/diagonalizeVel",  1. );
    M_diagonalize_p = dataFile( "fluid/space_discretization/diagonalizePress",  1. );

    M_hasVariationalPressure = dataFile( "fluid/space_discretization/variationalPn", false);
    if (M_verbose)
    {
        if (M_hasVariationalPressure)
            std::cout << "\n  f- Chorin-Temam coupling : variational pressure" << std::endl;
        else
            std::cout << "\n  f- Chorin-Temam coupling : standard" << std::endl;
    }

    M_linearSolver_u.setDataFromGetPot( dataFile, "fluid/solver" );
    M_linearSolver_p.setDataFromGetPot( dataFile, "fluid/solver" );

    // Fill stabilization arguments
    std::string _stabType = dataFile ( "fluid/stabilization/stabtype", "sd_stab" );
    if (_stabType == "sd_stab")
        M_stabType = SD_STAB;
    else if (_stabType == "ip_stab_expl")
        M_stabType = IP_STAB_EXPL;
    else if (_stabType == "ip_stab_impl")
        M_stabType = IP_STAB_IMPL;
    else
    {
        std::cout << "  -ERROR : Stabilization type undefined -> Exiting ..." << std::endl;
        exit(1);
    }

    if (M_stabType == SD_STAB)
    {
        M_sdStab.reset (new SDStabilization<Mesh, Dof>( M_uFESpace.mesh(),
                                                        M_uFESpace.dof(),
                                                        M_uFESpace.refFE(),
                                                        M_uFESpace.qr(),
                                                        0., 0., M_data.viscosity()) );
        M_gammaBeta = dataFile( "fluid/sdstab/gammaBeta", 0.);
        M_gammaDiv  = dataFile( "fluid/sdstab/gammaDiv", 0.);
        M_sdStab->setGammaBeta(M_gammaBeta);
        M_sdStab->setGammaDiv(M_gammaDiv);
        if (M_verbose)
            std::cout << "  -f SD stabilization is set" << std::endl;
    }
    if (M_stabType == IP_STAB_EXPL || M_stabType == IP_STAB_IMPL)
    {
        M_ipStab.reset (new IPStabilization<Mesh, Dof>( M_uFESpace.mesh(),
                                                        M_uFESpace.dof(),
                                                        M_uFESpace.refFE(),
                                                        M_uFESpace.feBd(),
                                                        M_uFESpace.qr(),
                                                        0., 0., 0., M_data.viscosity()) );
        M_gammaBeta = dataFile( "fluid/ipstab/gammaBeta", 0.);
        M_gammaDiv  = dataFile( "fluid/ipstab/gammaDiv", 0.);
        M_ipStab->setGammaBeta(M_gammaBeta);
        M_ipStab->setGammaDiv(M_gammaDiv);
        M_ipStab->setGammaPress(0.);
        if (M_verbose)
            std::cout << "  f- IP Stabilization is set" << std::endl;
    }


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
void ChorinTemam<Mesh, SolverType>::buildSystem_u_p()
{

    // We compute u and p systems _constant_ matrices at the same time
    // one loop on geom element w/ possibly different fe for u and p

    M_matrMass.reset  ( new matrix_type(M_localMap_u) );
    M_matrStokes.reset( new matrix_type(M_localMap_u) );
    M_matrPress.reset ( new matrix_type(M_localMap_p) );

    if (M_verbose) std::cout << "  f-  Computing constant matrices ...        ";

    // See if we are steady
    if (M_verbose)
    {
        if (M_steady)
            std::cout << "  fd-  Steady state ...." << std::endl;
        else
            std::cout << "  fd-  Unsteady state ...." << std::endl;
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

    if (false)
        std::cout << "partial times:  \n"
                  << " Der            " << chronoDer.diffCumul() << " s.\n"
                  << " Stab           " << chronoStab.diffCumul() << " s.\n"
                  << " Zero           " << chronoZero.diffCumul() << " s.\n"
                  << " Stiff          " << chronoStiff.diffCumul() << " s.\n"
                  << " Stiff Assemble " << chronoStiffAssemble.diffCumul() << " s.\n"
                  << " Mass           " << chronoMass.diffCumul() << " s.\n"
                  << " Mass Assemble  " << chronoMassAssemble.diffCumul() << " s.\n"
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
    M_sol_p_phi *= 0;

    if (M_verbose)
        std::cout << "  f-  Initializing velocity bdf object ..." << std::endl;
    M_bdf_u.initialize_unk(M_sol_u);

    if (M_verbose)
        std::cout << "  f-  Initializing projector field bdf object ..." << std::endl;
    M_bdf_p_phi.initialize_unk(M_sol_p_phi);

    if (M_verbose)
    {
        if (M_pBdfOrder > 0)
            std::cout << "  f-  Initializing pressure bdf object ..." << std::endl;
        else
            std::cout << "  f-  Standard Chorin-Temam projection ..." << std::endl;
    }
    if (M_pBdfOrder > 0)
        M_bdf_p.initialize_unk(M_sol_p);

}

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::updateSystem_u(vector_type& betaVec)
{

    Real alpha = M_bdf_u.coeff_der(0) / M_data.dataTime()->getTimeStep();

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

    // Update convective term for velocity system

    UInt nbCompU = nDimensions;
    UInt velTotalDof = M_uFESpace.dof().numTotalDof();
    Chrono chrono;

    double normInf;
    betaVec.NormInf(&normInf);

    if (normInf == 0.)
        std::cout << "\n  f-  Convective velocity is zero on process : " <<
                  M_me << std::endl;

    if (normInf != 0.)
    {

        if (M_verbose)
            std::cout << "\n  f-  Updating the convective terms ...        "
                      << std::flush;

        // vector with repeated nodes over the processors
        vector_type betaVecRep( betaVec, Repeated );

        chrono.start();

        for ( UInt iVol = 1; iVol<= M_uFESpace.mesh()->numVolumes(); ++iVol )
        {

            M_uFESpace.fe().updateFirstDeriv( M_uFESpace.mesh()->volumeList( iVol ) );
            // reuse elmatStiff as space for storing elemental convective matrix
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

        } // volume for loop


        chrono.stop();
        if (M_verbose) std::cout << "done in " << chrono.diff() << " s.\n"
                                     << std::flush;

        // Handles  stabilization terms
        if (M_verbose)
            std::cout << "  f-  Adding stabilization terms ..." << std::endl;


        if (M_stabType == SD_STAB)
        {
            std::cout << "  f-  Adding SD stabilization (implicit) terms" << std::endl;
            // We must multiply betaVecRep by density before this step
            betaVecRep *= M_data.density();
            M_sdStab->applyCT(M_data.dataTime()->getTimeStep(), *M_matrNoBC_u, betaVecRep);
        }
        else if (M_stabType == IP_STAB_EXPL)
        {
            std::cout << "  f-  Adding IP stabilization (explicit) terms" << std::endl;
            M_ipStab->apply_expl(M_rhsNoBC_u, betaVecRep);
        }
        else if (M_stabType == IP_STAB_IMPL)
        {
            // just to compare with previous code between implicit and explicit IP stab
            std::cout << "  f-  Adding IP stabilization (implicit) terms" << std::endl;
            M_ipStab->apply(*M_matrNoBC_u, betaVecRep, M_verbose);
        }

    } // if (normInf != 0.)

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

    if (M_verbose)
        std::cout << "  f-  Updating velocity system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start();

    // vector with repeated nodes over the processors
    vector_type press_sol (p_sol, Repeated);

    M_rhsNoBC_u *= 0.0;
    vector_type _rhsNoBC_u ( M_rhsNoBC_u, Repeated, Zero);

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
        if (M_hasVariationalPressure)
        {
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp)
            {
                source_pdivv( 1.0, M_elvec_p, M_elemrhs_u,
                              M_pFESpace.fe(), M_uFESpace.fe(), iComp);
                assembleVector( _rhsNoBC_u, M_elemrhs_u, M_uFESpace.fe(),
                                M_uFESpace.dof(), iComp, iComp*velTotalDof);
            }
        }
        else
        {
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

    M_rhsNoBC_u = _rhsNoBC_u;

    // note: communication, i.e. GlobalAssembling, will be done after.
    chronoCTrhs.stop();

    if (M_verbose)
        std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_u

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::computeCTRHS_p(vector_type& u_sol)
{
    Chrono chronoCTrhs;
    Real bdf_coeff = M_bdf_u.coeff_der(0);

    UInt velTotalDof   = M_uFESpace.dof().numTotalDof();

    if (M_verbose)
        std::cout << "  f-  Updating pressure/projector system Chorin-Temam couplings ...      " << std::flush;

    chronoCTrhs.start();

    vector_type vel_sol( u_sol, Repeated );

    M_rhsNoBC_p *= 0.0;			// for test pupose
    vector_type _rhsNoBC_p ( M_rhsNoBC_p, Repeated, Zero);

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

        source_divuq( - bdf_coeff * M_data.density() / M_data.dataTime()->getTimeStep(), M_elvec_u,
                      M_elemrhs_p, M_uFESpace.fe(), M_pFESpace.fe() );
        assembleVector ( _rhsNoBC_p, M_elemrhs_p, M_pFESpace.fe(),
                         M_pFESpace.dof(),0);

    }

    // be sure assembling has completed
    M_comm->Barrier();

    M_rhsNoBC_p = _rhsNoBC_p;

    // note: all communication should be done later

    chronoCTrhs.stop();

    if (M_verbose)
        std::cout << "done in " << chronoCTrhs.diff() << " s." << std::endl;
} // computeCTRHS_p

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::setRES(const std::vector<Real>& resVal, const std::vector<EntityFlag>& resFlag, const int explRES)
{
    if (resVal.size() != resFlag.size())
    {
        if (M_verbose)
        {
            std::cout << " -ERROR : wrong specification of resistance bc's" << std::endl;
            std::cout << " ... Exiting" << std::endl;
        }
        exit(1);
    }

    if (M_verbose)
        std::cout << "  f-  Setting resistance(s) into fluid solver ..." << std::endl;

    M_resVal = resVal;
    M_resFlag = resFlag;
    M_hasRES = true;

    if (explRES == 1)
    {
        if (M_verbose)
        {
            std::cout << "       Resistance coupling is explicit." << std::endl;
        }
        M_hasRESexpl = true;
    }
    else if (explRES == 0)
    {
        if (M_verbose)
        {
            std::cout << "       Resistance coupling is implicit." << std::endl;
        }
        M_hasRESexpl = false;
    }
    else
    {
        if (M_verbose)
        {
            std::cout << "  -ERROR : wrong specification of resistance bc's" << std::endl;
            std::cout << " ... Exiting" << std::endl;
        }
        exit(1);
    }

} // setRES

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::computeRES_expl(Real resVal, EntityFlag resFlag, vector_type& betaVec, vector_type& rhs)
{
    BCBase& bcB = M_BCh_fluid_u->GetBCWithFlag(resFlag);
    vector_type betaVecRep (betaVec, Repeated);
    Real uLocalFlow = 0.0, uFlow = 0.0;
    UInt nNode = M_uFESpace.feBd().nbNode;
    // auxilliary map
    EpetraMap fe_map (M_uFESpace.refFE(), *(M_uFESpace.mesh()), *M_comm);
    // auxilliary (repeated) vector
    vector_type uFlux(rhs, Repeated);
    uFlux *= 0.0;

    // loop on (usually boundary) faces w/ resFlag as reference
    for (UInt i=1; i <=  bcB.list_size(); i++)
    {
        // pointer to the ith identifier in the bc_base list
        const IdentifierNatural* ptr_i = static_cast<const IdentifierNatural*> (bcB(i));
        // index of corresponding face
        UInt iFace = ptr_i->id();

        // update current face
        M_uFESpace.feBd().updateMeasNormalQuadPt(M_uFESpace.mesh()->boundaryFace(iFace));
        //UInt faceID = M_uFESpace.feBd().currentLocalId();

        // get local velocity vector
        M_belvec_u.zero();
        for (UInt iNode=0; iNode < nNode; ++iNode)
        {
            //UInt iloc = M_uFESpace.feBd().patternFirst(iNode);
            int iloc = iNode;
            for (UInt iComp=0; iComp < nDimensions; ++iComp)
            {
                UInt ilocg = ptr_i->bdLocalToGlobal(iloc+1);
                UInt ig = fe_map.getMap(Repeated)->MyGlobalElements()[ilocg];
                ig += iComp*dim_u();
                M_belvec_u.vec()[iloc + iComp*nNode] = betaVecRep[ig];
            }
        }

        // get local flow through current face and compute local FE flow vector
        M_bele_flow_u.zero();
        for (int iNode=0; iNode < nNode; ++iNode)
        {
            int iloc = iNode;
            for (int iComp=0; iComp < nDimensions; ++iComp)
            {
                for (int iQuad=0; iQuad < M_uFESpace.feBd().nbQuadPt; ++iQuad)
                {
                    M_bele_flow_u.vec() [iloc + iComp*nNode] +=
                        M_uFESpace.feBd().phi(iNode, iQuad) *
                        M_uFESpace.feBd().normal(iComp, iQuad) *
                        M_uFESpace.feBd().weightMeas(iQuad);
                    uLocalFlow += M_belvec_u.vec() [iloc + iComp*nNode] *
                                  M_bele_flow_u.vec() [iloc + iComp*nNode];
                }
            }
        }

        // assemble into global rhs vector
        // (no call to assembleVector here since we are on a boundary FE)
        for (UInt iNode = 0; iNode < nNode; ++iNode)
        {
            int iloc = iNode;
            for (UInt iComp=0; iComp < nDimensions; ++iComp)
            {
                UInt ilocg = ptr_i->bdLocalToGlobal(iloc+1); //+ iComp * dim_u();
                UInt ig =  fe_map.getMap(Repeated)->MyGlobalElements()[ilocg] + iComp * dim_u();
                uFlux[ ig ] += M_bele_flow_u.vec() [iloc + iComp*nNode];
            }
        }
    }
    // compute global flow
    fe_map.Comm().SumAll(&uLocalFlow, &uFlow, 1);
    // compute global flux vector
    uFlux.GlobalAssemble();
    // compute actual weak resistance term. Check sign.
    uFlux *= (-resVal*uFlow);
    // update rhs coupling
    rhs += uFlux;
} //computeRES_expl

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::computeRES_impl(Real resVal, EntityFlag resFlag, matrix_type& matrix)
{
    BCBase& bcB = M_BCh_fluid_u->GetBCWithFlag(resFlag);
    UInt nNode = M_uFESpace.feBd().nbNode;
    // auxilliary map
    EpetraMap fe_map (M_uFESpace.refFE(), *(M_uFESpace.mesh()), *M_comm);
    // FE flow vector
    vector_type _flow_u (M_localMap_u);
    _flow_u *= 0.0;

    std::cout << "      RES: computing FE flow vector ..." << std::endl;

    // loop on (usually boundary) faces w/ resFlag as reference
    for (UInt i=1; i <= bcB.list_size(); ++i)
    {
        // pointer to the ith identifier in the bc_base list
        const IdentifierNatural* ptr_i = static_cast<const IdentifierNatural*> (bcB(i));
        // index of corresponding face
        UInt iFace = ptr_i->id();

        // update current face
        M_uFESpace.feBd().updateMeasNormalQuadPt(M_uFESpace.mesh()->boundaryFace(iFace));

        // get local FE flow vector
        M_bele_flow_u.zero();
        for (int iNode = 0; iNode < nNode; ++iNode)
        {
            //UInt iloc = M_uFESpace.feBd().patternFirst(iNode);
            UInt iloc = iNode;
            for (int iComp = 0; iComp < nDimensions; ++iComp)
            {
                for (int iQuad=0; iQuad < M_uFESpace.feBd().nbQuadPt; ++iQuad)
                {
                    M_bele_flow_u.vec() [iloc + iComp*nNode] +=
                        M_uFESpace.feBd().phi(iNode, iQuad) *
                        M_uFESpace.feBd().normal(iComp, iQuad) *
                        M_uFESpace.feBd().weightMeas(iQuad);
                }
            }
        }

        // assemble into _flow_u
        // (no call to assembleVector since we are on a boundary FE)
        for (UInt iNode = 0; iNode < nNode; ++iNode)
        {
            //UInt iloc = M_uFESpace.feBd().patternFirst(iNode);
            UInt iloc = iNode;
            for (UInt iComp=0; iComp < nDimensions; ++iComp)
            {
                UInt ilocg = ptr_i->bdLocalToGlobal(iloc+1);
                UInt ig = fe_map.getMap(Repeated)->MyGlobalElements()[ilocg];
                ig += iComp * dim_u();
                _flow_u[ ig ] += M_bele_flow_u.vec() [iloc + iComp*nNode];
            }
        }
    }

    M_comm->Barrier();
    // compute global FE flow vector
    _flow_u.GlobalAssemble();

    std::cout << "      RES: computing FE flow matrix" << std::endl;
    std::cout << "      RES: __not implemented__ " << std::endl;

    // __todo__ :
    // now modify velocity matrix to add implicit resistance term: this will
    // enlarge matrix profile (non local coupling between boundary dofs
    // flagged as resistance).
    // __note_for_later__ : do not know how to modify matrix/vector product
    // with epetra things either.

    M_comm->Barrier();

} // computeRES_impl

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::applyRES_expl(vector_type& vector)
{
    int nRes = M_resFlag.size();
    vector_type _u_extrap( M_localMap_u );
    _u_extrap *= 0.;
    _u_extrap = M_bdf_u.extrap();

    for (int iRes=0; iRes < nRes; ++iRes)
    {
        computeRES_expl(M_resVal[iRes], M_resFlag[iRes], _u_extrap, vector);
    }
} // applyRES_expl


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::applyRES_impl(matrix_type& matrix)
{
    int nRes = M_resFlag.size();
    for (int iRes=0; iRes < nRes; ++iRes)
    {
        computeRES_impl(M_resVal[iRes], M_resFlag[iRes], matrix);
    }
} // applyRES_impl

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::iterate_u( bchandler_raw_type& bch_u)
{

    Chrono chrono;

    // update convective term with extrapolated velocity
    // little bs as usual w/ epetra things
    vector_type _u_extrap( M_localMap_u );
    _u_extrap *= 0.;
    _u_extrap = M_bdf_u.extrap();
    updateSystem_u( _u_extrap );

    // set resistance boundary conditions if any
    if (M_hasRES && M_hasRESexpl)
    {
        if (M_verbose)
        {
            std::cout << "  f-  Applying resistance bc's (explicit)..." << std::endl;
        }
        applyRES_expl ( M_rhsNoBC_u );
    }

    if (M_hasRES && !M_hasRESexpl)
    {
        if (M_verbose)
        {
            std::cout << "  f-  Applying resistance bc's (implicit)..." << std::endl;
        }
        applyRES_impl (*M_matrNoBC_u);
    }

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

    // update velocity bdf components
    M_bdf_u.shift_right(M_sol_u);

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
    M_residual_p -= *M_matrNoBC_p*M_sol_p_phi;

    // Compute physical pressure multiplier M_sol_p = M_sol_p_phi + M_sol_p_star
    // with M_sol_p_star = extrapolation of M_sol_p at current time.
    // So we have M_sol_p = M_sol_p_phi for pBdfOrder = 0 (std Chorin-Temam projection)
    M_sol_p *= 0;
    M_sol_p += M_sol_p_phi;
    if (M_pBdfOrder > 0)
        M_sol_p += M_bdf_p.extrap();

    // update pressure and projector bdf
    M_bdf_p_phi.shift_right(M_sol_p_phi);
    if (M_pBdfOrder > 0)
        M_bdf_p.shift_right(M_sol_p);

} // iterate_p()

template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::time_advance(Real const& time)
{

    Real dt = M_data.dataTime()->getTimeStep();

    // what to do at first time step
    if (M_firstTimeStep)
    {
        // note: must call initialize() outside the solver for proper init
        buildSystem_u_p();
        resetPrec_u();
        resetPrec_p();
        M_firstTimeStep = 0;
    }

    // update velocity system rhs

    // compute pressure correction to velocity \tilde u equation
    vector_type _p_corr( M_localMap_p );
    _p_corr *= 0.;
    if (M_pBdfOrder > 0)
        _p_corr = M_bdf_p.extrap();
    _p_corr += M_bdf_p_phi.time_der(1);
    _p_corr *= 1./M_bdf_p_phi.coeff_der(0);

    // compute Chorin-Temam coupling term
    computeCTRHS_u(_p_corr);

    // add the mass term rho/dt * u^n
    if (M_verbose)
        std::cout << "  f-  Adding velocity mass term on rhs ..." << std::endl;

    // note: this is 1 and not M_data.getTimeStep() as M_matrMass already contains
    // the timestep term to avoid recomputation (constant matrix)
    M_rhsNoBC_u += *M_matrMass * M_bdf_u.time_der( 1 );

    // note: velocity updating is done in iterate_u method

} // time_advance


template<typename Mesh, typename SolverType>
void ChorinTemam<Mesh, SolverType>::solveSystem_u( matrix_ptrtype  matrFull,
                                                   vector_type&    rhsFull )
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

    int numIter = M_linearSolver_u.solve(M_sol_u, rhsFull);

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

        numIter = M_linearSolver_u.solve(M_sol_u, rhsFull);

        if (numIter >= M_maxIterSolver_u)
        {
            if (M_verbose)
            {
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
void ChorinTemam<Mesh, SolverType>::solveSystem_p( matrix_ptrtype  matrFull,
                                                   vector_type&    rhsFull )
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

    int numIter = M_linearSolver_p.solve(M_sol_p_phi, rhsFull);

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

        numIter = M_linearSolver_p.solve(M_sol_p_phi, rhsFull);

        if (numIter >= M_maxIterSolver_p)
        {
            if (M_verbose)
            {
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
ChorinTemam<Mesh, SolverType>::postProcess(bool _writeMesh)
{
    // maybe we will add something useful here
}



} // namespace LifeV


#endif //_CHORIN_TEMAM_H_
