/* -*- mode: c++ -*-

 This file is part of the LifeV library

 Author(s): Miguel A. Fernandez <miguel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
            Tiziano Passerini <tiziano@mathcs.emory.edu>

 Copyright (C) 2010 EPFL, Emory University

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
  \file AdvectionDiffusionReactionSolver.hpp
  \date 08/2010

  \brief This file contains a solver class for
         Advection-Diffusion-Reaction problems
 */
#ifndef _ADRSOLVER_H_
#define _ADRSOLVER_H_

#include <life/lifecore/chrono.hpp>

#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bdfVariableStep.hpp>

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataADR.hpp>

#include <list>
#include <sstream>

#include <boost/scoped_ptr.hpp>

namespace LifeV
{
/*!
  \typedef What constant matrices should be constructed?
 */
enum ADRConstantMatrices
{
    ADR_NO_CONST_MATRICES, //!< No constant matrices
    ADR_ONLY_MASS,         //!< No diffusion terms
    ADR_ONLY_STIFF,        //!< No reaction terms and steady problem
    ADR_MASS_STIFF         //!< Diffusion and (reaction or unsteady)
};


//! default function for source and advection terms
Real ADRfZero( const Real&, const Real&, const Real&, const Real&, const ID& )
{
    return 0.;
}


/*!
  \class ADRSolver

  This class implements a solver for Advection-Diffusion-Reaction problems.

  The reference PDE problem has the form

  \f$ \frac{ \partial u }{ \partial t } -
   - \mu \triangle u + \beta \cdot \nabla u + \sigma u =
   f \f$,

  where \f$\mu\$ is the diffusion coefficient, \f$\beta\$ is the advection field
  and \f$\sigma\$ is the reaction coefficient. The function \f$\f\f$ is also referred
  to as the source term.

  The parameters \f$\mu\f$ and \f$\sigma\f$ are assumed to be constant in space and time,
  while the vector field \f$\beta\f$ can vary in time and space.

  The advection field and the source term can be specified by the user as functions of type
  ADRSolver::function_type or as vectors of type ADRSolver::vector_type. In the former case,
  the function is evaluated on the quadrature nodes when assemblying respectively the (elementary)
  matrix associated to the advection terms and the (elementary) vector associated to the source
  term in the weak formulation of the problem.

  In the latter case, the vector must contain respectively the right hand side of the system
  (both the source term already assembled and the time advancing terms, if this applies);
  and the advection field, as the vector of coefficients of the expansion wrt the basis of a
  given finite element space. If necessary, i. e. if the finite element space of the advection
  field does not coincide with the finite element space of the solution, the solver will compute
  the interpolant function belonging to the finite element space of the unknown field.
 */
template< typename Mesh,
typename SolverType = LifeV::SolverTrilinos >
class ADRSolver
{
public:

    //! @name Type definitions
    //@{
    typedef Mesh                                  mesh_type;

    typedef typename SolverType::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptr_type;
    typedef typename SolverType::vector_type      vector_type;

    typedef FESpace<Mesh, EpetraMap>              fespace_type;
    typedef fespace_type*                         fespace_rawptr_type;
    typedef boost::shared_ptr<fespace_type>       fespace_ptr_type;

    typedef DataADR                               data_type;
    typedef data_type*                            data_rawptr_type;
    typedef boost::shared_ptr<data_type>          data_ptr_type;

    typedef BdfVS<vector_type>                    timeIntegrator_type;
    typedef boost::shared_ptr<timeIntegrator_type> timeIntegrator_ptr_type;

    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> function_type;

    typedef BCHandler                             bchandler_type;
    typedef boost::shared_ptr<bchandler_type>     bchandler_ptr_type;

    typedef boost::scoped_ptr<ElemMat>            elmat_ptr_type;
    typedef boost::scoped_ptr<ElemVec>            elvec_ptr_type;
    //@}


    //! @name Constructors & Destructor
    //@{

    //! Default constructor
    ADRSolver();

    //! Copy constructor
    /*!
      \param adrsolver the solver to be copied

      not working yet (requires copy constructors for members)
     */
    // ADRSolver( const ADRSolver& adrsolver );

    //! destructor
    ~ADRSolver();
    //@}


    //! @name Operators
    //@{
    //! Operator=
    /*!
     * @param adrsolver an object of type ADRSolver
     *
     * BEWARE this won't work until we have copy constructors
     * for ElemMat & ElemVec!
     */
    // ADRSolver& operator=( const ADRSolver& adrsolver );

    //@}


    //! @name Methods
    //@{

    //! Sets the user-defined parameters
    /*!
     * @param dataFile the GetPot container
     * @section prefix for accessing GetPot container
     */
    void setupLinearSolver( const GetPot& dataFile,
                            const std::string& section = "adr" );

    //! Initializes the solver data structures
    virtual void setup();

    //! Assembles the constant matrices
    virtual void computeConstantMatrices( ADRConstantMatrices whichMatrices = ADR_MASS_STIFF );

    //! Sets the initial condition
    /*!
     * @param the initial value of the unknown field
     */
    void initialize( const function_type& );
    void initialize( const vector_type& );

    //! Fills the right hand side of the linear system
    /*!
     * @param rhsVec the term on the right hand side of the system
     */
    void updateRHS( vector_type const& rhsVec );

    //! Fills the right hand side of the linear system (use the internally stored sourceTerm)
    /*!
     * This method is alternative to updateRHS( rhsVec ). It assembles the source term
     * using its functional form, and manages the time advance if this applies.
     *
     * @param timeIntegratorPtr a pointer to the manager for the time integration
     */
    void updateRHS( timeIntegrator_ptr_type timeIntegratorPtr = timeIntegrator_ptr_type() );

    //! Assembles the time-varying matrices
    /*!
     * @param timeIntegratorPtr a pointer to the manager for the time integration
     */
    virtual void updateMatrix( timeIntegrator_ptr_type timeIntegratorPtr = timeIntegrator_ptr_type() );

    //! Solves the linear system
    /*!
     * @param bch the data structure holding the boundary conditions
     */
    void iterate( bchandler_type& bch );

    //@}


    //! @name Set methods
    //@{
    //! Sets the source term
    /*!
     * @param the source term in functional form
     */
    void setSourceTerm( const function_type& sourceTerm ) { M_sourceTerm = sourceTerm; }

    //! Sets the advection field
    /*!
     * @param the advection term in functional form
     */
    void setAdvectionField( const function_type& advectionField ) { M_advectionField = advectionField; }

    //! Sets the advection field
    /*!
     * @param the advection term in vector form
     */
    void setAdvectionField( const vector_type& advectionVector = vector_type(),
                            fespace_ptr_type advectionFieldFESpacePtr = fespace_ptr_type() );

    //! set the pointer to the data object
    /*!
     * @param dataPtr a pointer to the data object
     */
    void setDataPtr( data_ptr_type dataPtr ) { M_dataPtr = dataPtr; }

    //! set the pointer to the FE space descriptor for the unknown
    /*!
     * @param uFESpacePtr a pointer to the FE space descriptor for the unknown
     */
    void setUFESpacePtr( fespace_ptr_type uFESpacePtr ) { M_uFESpacePtr = uFESpacePtr; }

    //! set the pointer to the FE space descriptor for the advection field
    /*!
     * @param betaFESpacePtr a pointer to the FE space descriptor for the advection field
     */
    // void setBetaFESpacePtr( fespace_ptr_type betaFESpacePtr ) { M_advectionFieldFESpacePtr = betaFESpacePtr; }

    //! set the bool for recomputing the matrix
    /*!
     * @param recomp true if the matrix should be recomputed
     */
    // void setRecomputeMatrix( bool const recomp ){ M_recomputeConstantMatrix = recomp; }

    //@}


    //! @name Get methods
    //@{
    //! get the displayer
    // const Displayer& displayer() const { return M_displayer; }

    //! returns a pointer to the data object
    data_ptr_type dataPtr() const { return M_dataPtr; }

    //! returns a pointer to the FE space descriptor for the unknown
    // fespace_ptr_type uFESpacePtr() const { return M_uFESpacePtr; }

    //! returns a pointer to the FE space descriptor for the advection field
    // fespace_ptr_type betaFESpacePtr() const { return M_advectionFieldFESpacePtr; }

    //! get the pointers to the matrices
    // matrix_ptr_type matrMassPtr() const { return M_matrMassPtr; }
    // matrix_ptr_type matrNoBCPtr() const { return M_matrNoBCPtr; }

    //! returns the local solution vector
    const vector_type& solution() const { return M_sol; }

    //! returns the local residual vector
    const vector_type& residual() const { return M_residual; }

    //! get the source term functor
    // const source_type& sourceTerm() const { return M_sourceTerm; }

    //@}


protected:

    void initializeMatrix();

    void addDiffusionTerms();

    void addReactionTerms();

    void addAdvectionTerms();

    void addStabilizationTerms();

    void addMatrixTimeTerm( timeIntegrator_ptr_type timeIntegratorPtr );

    void applyBoundaryConditions(   matrix_type&        matrix,
                                    vector_type&        rhs,
                                    bchandler_type&     BCh );

    void assembleSourceTerm();

    void addRHSTimeTerm( timeIntegrator_ptr_type timeIntegratorPtr );

    //! Manager for output communications
    Displayer                       M_displayer;

    //! Map for building matrices and vectors
    // EpetraMap                       M_localMap;

    //! data for ADR solver
    data_ptr_type                   M_dataPtr;

    // FE space for the unknown
    fespace_ptr_type                M_uFESpacePtr;

    // FE space for the advection field
    fespace_ptr_type                M_advectionFieldFESpacePtr;

    //! Elementary matrices and vectors
    elmat_ptr_type                  M_elmatMassPtr;       // mass
    elmat_ptr_type                  M_elmatStiffPtr;      // stiffness
    elmat_ptr_type                  M_elmatAdvPtr;        // advection
    elvec_ptr_type                  M_elvecBetaPtr;        // Elementary right hand side
    elmat_ptr_type                  M_elmatStabPtr;
    elvec_ptr_type                  M_elvecRhsPtr;        // Elementary right hand side

    //! mass matrix
    matrix_ptr_type                 M_matrMassPtr;

    //! Stiffness Matrix
    matrix_ptr_type                 M_matrStiffPtr;

    //! Advection Matrix
    matrix_ptr_type                 M_matrAdvPtr;

    //! matrix without boundary conditions
    matrix_ptr_type                 M_matrNoBCPtr;

    //! stabilization matrix
    matrix_ptr_type                 M_matrStabPtr;

    //! Right hand side for the system
    vector_type                     M_rhsNoBC;

    //! Global solution
    vector_type                     M_sol;

    //! Residual
    vector_type                     M_residual;

    //! Source term
    function_type                   M_sourceTerm;

    //! Advection field
    function_type                   M_advectionField;
    vector_type                     M_advectionVector;

    //! The linear solver
    SolverType                      M_linearSolver;

    ADRConstantMatrices             M_whichMatrices;

    bool                            M_recomputeConstantMatrix;

}; // class ADRSolver



//
// IMPLEMENTATION
//

template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
ADRSolver() :
M_displayer              (),
M_dataPtr                ( new data_type() ),
M_whichMatrices          ( ADR_NO_CONST_MATRICES ),
M_recomputeConstantMatrix( false )
{
    this->M_sourceTerm = boost::bind(&ADRfZero, _1, _2, _3, _4, _5);
    this->M_advectionField = boost::bind(&ADRfZero, _1, _2, _3, _4, _5);
}


/* BEWARE this won't work until we have copy constructors
 * for ElemMat & ElemVec!
template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
ADRSolver( const ADRSolver& adrsolver ):
M_displayer              ( adrsolver.M_displayer ),
M_dataPtr                ( adrsolver.M_dataPtr ),
M_uFESpacePtr            ( adrsolver.M_uFESpacePtr ),
M_advectionFieldFESpacePtr         ( adrsolver.M_advectionFieldFESpacePtr ),
M_rhsNoBC                ( adrsolver.M_rhsNoBC ),
M_sol                    ( adrsolver.M_sol ),
M_residual               ( adrsolver.M_residual ),
M_linearSolver           ( adrsolver.M_linearSolver ),
M_steady                 ( adrsolver.M_steady ),
M_stab                   ( adrsolver.M_stab ),
M_dataPtr->stabilizationCoefficient()              ( adrsolver.M_dataPtr->stabilizationCoefficient() ),
M_recomputeMatrix        ( adrsolver.M_recomputeMatrix )
{
    if( adrsolver.*M_elmatMassPtr )
        M_elmatMassPtr.reset( new ElemMat( *adrsolver.*M_elmatMassPtr ) );
    if( adrsolver.*M_elmatStiffPtr )
        M_elmatStiffPtr.reset( new ElemMat( *adrsolver.*M_elmatStiffPtr ) );
    if( adrsolver.*M_elmatAdvPtr )
        M_elmatAdvPtr.reset( new ElemMat( *adrsolver.*M_elmatAdvPtr ) );
    if( adrsolver.M_elvecBeta )
        M_elvecBetaPtr.reset( new ElemMat( *adrsolver.M_elvecBeta ) );
    if( adrsolver.*M_elmatStabPtr )
        M_elmatStabPtr.reset( new ElemMat( *adrsolver.*M_elmatStabPtr ) );

    if( adrsolver.M_matrMassPtr )
            M_matrMassPtr.reset( new matrix_type( *adrsolver.M_matrMassPtr ) );
    if( adrsolver.M_matrStiffPtr )
            M_matrStiffPtr.reset( new matrix_type( *adrsolver.M_matrStiffPtr ) );
    if( adrsolver.M_matrAdvPtr )
            M_matrAdvPtr.reset( new matrix_type( *adrsolver.M_matrAdvPtr ) );
    if( adrsolver.M_matrNoBCPtr )
            M_matrNoBCPtr.reset( new matrix_type( *adrsolver.M_matrNoBCPtr ) );
    if( adrsolver.M_matrStabPtr )
            M_matrStabPtr.reset( new matrix_type( *adrsolver.M_matrStabPtr ) );
}
 */


template<typename Mesh, typename SolverType>
ADRSolver<Mesh, SolverType>::
~ADRSolver()
{
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::setupLinearSolver( const GetPot& dataFile, const std::string& section )
{
    // We want a slash dividing the data file section from the variable name but only
    // when not provided by the user or when not looking in the root of the data file
    std::string corrected_section( section );
    if( ( ! section.empty() ) && ( section[section.length()-1] != '/' ) )
        corrected_section = section + '/';

    M_linearSolver.setDataFromGetPot( dataFile, (corrected_section+"solver").data() );
    M_linearSolver.setUpPrec( dataFile, (corrected_section+"prec").data() );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::setup()
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );

    M_advectionFieldFESpacePtr.reset( new fespace_type( M_uFESpacePtr->mesh(),
                                                        M_dataPtr->solFEType(),
                                                        nDimensions,
                                                        M_uFESpacePtr->map().CommPtr() ) );

    // if ( M_dataPtr->reactionCoefficient() || ( !M_dataPtr->steady() ) )
    M_elmatMassPtr.reset( new ElemMat( M_uFESpacePtr->fe().nbFEDof(),
                                       M_uFESpacePtr->fieldDim(),
                                       M_uFESpacePtr->fieldDim() ) );

    // if ( M_dataPtr->diffusionCoefficient() )
    M_elmatStiffPtr.reset( new ElemMat( M_uFESpacePtr->fe().nbFEDof(),
                                        M_uFESpacePtr->fieldDim(),
                                        M_uFESpacePtr->fieldDim() ) );

    // TODO how to filter out these, when not needed?
    M_elmatAdvPtr.reset( new ElemMat( M_uFESpacePtr->fe().nbFEDof(),
                                      M_uFESpacePtr->fieldDim(),
                                      M_uFESpacePtr->fieldDim() ) );

    M_elvecBetaPtr.reset( new ElemVec( M_uFESpacePtr->fe().nbFEDof(),
                                       nDimensions ) );

    // if( M_dataPtr->stabilizationMethod() != ADR_NO_STABILIZATION )
    M_elmatStabPtr.reset( new ElemMat( M_uFESpacePtr->fe().nbFEDof(),
                                       M_uFESpacePtr->fieldDim(),
                                       M_uFESpacePtr->fieldDim() ) );

    M_elvecRhsPtr.reset( new ElemVec( M_uFESpacePtr->fe().nbFEDof(),
                                      M_uFESpacePtr->fieldDim() ) );

    M_rhsNoBC.setMap( this->M_uFESpacePtr->map() );
    M_rhsNoBC = EpetraVector( M_uFESpacePtr->map() );
    M_sol.setMap( this->M_uFESpacePtr->map() );
    M_sol = EpetraVector( M_uFESpacePtr->map() );
    M_residual.setMap( this->M_uFESpacePtr->map() );
    M_residual = EpetraVector( M_uFESpacePtr->map() );

    M_advectionVector.setMap( this->M_advectionFieldFESpacePtr->map() );
    // M_advectionVector = EpetraVector( M_uFESpacePtr->map() );

    M_displayer.setCommunicator( M_uFESpacePtr->map().CommPtr() );
    M_linearSolver.setCommunicator( M_uFESpacePtr->map().CommPtr() );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::computeConstantMatrices( ADRConstantMatrices whichMatrices )
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );

    // initialize the mass matrix
    if ( whichMatrices!=ADR_ONLY_STIFF )
    {
        ASSERT( M_elmatMassPtr.get(), "\nThe elementary matrix for the mass assembly is undefined!\n" );
        M_matrMassPtr.reset( new matrix_type( M_uFESpacePtr->map() ) );
    }
    // initialize the stiffness matrix
    if ( whichMatrices!=ADR_ONLY_MASS )
    {
        ASSERT( M_elmatStiffPtr.get(), "\nThe elementary matrix for the stiff assembly is undefined!\n" );
        M_matrStiffPtr.reset( new matrix_type( M_uFESpacePtr->map() ) );
    }

    M_displayer.leaderPrint("  adr-  Computing constant matrices ...");

    Chrono chrono;

    Chrono chronoDer;
    Chrono chronoZero;

    Chrono chronoMass;
    Chrono chronoStiff;

    Chrono chronoMassAssemble;
    Chrono chronoStiffAssemble;

    chrono.start();

    if ( whichMatrices!=ADR_NO_CONST_MATRICES )
    {
        // Elementary computation and matrix assembling
        // Loop on elements
        for ( UInt iEl = 1; iEl <= M_uFESpacePtr->mesh()->numElements(); iEl++ )
        {
            // update current FE info
            chronoDer.start();
            M_uFESpacePtr->fe().updateFirstDeriv( M_uFESpacePtr->mesh()->element( iEl ) );
            chronoDer.stop();

            // stiffness matrix
            if ( whichMatrices!=ADR_ONLY_MASS )
            {
                // clear data structures
                chronoZero.start();
                M_elmatStiffPtr->zero();
                chronoZero.stop();

                // elementary matrix
                chronoStiff.start();
                stiff( 1., *M_elmatStiffPtr,  M_uFESpacePtr->fe(), 0, 0, M_uFESpacePtr->fieldDim() );
                chronoStiff.stop();

                // assembly
                for ( UInt iComp = 0; iComp < M_uFESpacePtr->fieldDim(); iComp++ )
                {
                    chronoStiffAssemble.start();
                    assembleMatrix( *M_matrStiffPtr,
                                    *M_elmatStiffPtr,
                                    M_uFESpacePtr->fe(),
                                    M_uFESpacePtr->fe(),
                                    M_uFESpacePtr->dof(),
                                    M_uFESpacePtr->dof(),
                                    iComp, iComp,
                                    iComp*M_uFESpacePtr->dim(), iComp*M_uFESpacePtr->dim() );
                    chronoStiffAssemble.stop();
                }
            } // if (build the mass matrix)

            // mass matrix
            if ( whichMatrices!=ADR_ONLY_STIFF )
            {
                // clear data structures
                chronoZero.start();
                M_elmatMassPtr->zero();
                chronoZero.stop();

                // elementary matrix
                chronoMass.start();
                mass( 1., *M_elmatMassPtr, M_uFESpacePtr->fe(), 0, 0, M_uFESpacePtr->fieldDim());
                chronoMass.stop();

                // assembly
                for ( UInt iComp = 0; iComp < M_uFESpacePtr->fieldDim(); iComp++ )
                {
                    chronoMassAssemble.start();
                    assembleMatrix( *M_matrMassPtr,
                                    *M_elmatMassPtr,
                                    M_uFESpacePtr->fe(),
                                    M_uFESpacePtr->fe(),
                                    M_uFESpacePtr->dof(),
                                    M_uFESpacePtr->dof(),
                                    iComp, iComp,
                                    iComp*M_uFESpacePtr->dim(), iComp*M_uFESpacePtr->dim() );
                    chronoMassAssemble.stop();
                }

            } // if (build the stiff matrix)

        } // loop over mesh elements

    } // if (build const matrices)

    // synchronize the processes
    M_uFESpacePtr->map().Comm().Barrier();

    chrono.stop();
    M_displayer.leaderPrintMax( "\n\t... done in " , chrono.diff());


    M_displayer.leaderPrint( "  adr-  Finalizing the constant matrices ...");
    chrono.start();

    if ( whichMatrices!=ADR_ONLY_STIFF )
        M_matrMassPtr->GlobalAssemble();
    if ( whichMatrices!=ADR_ONLY_MASS )
        M_matrStiffPtr->GlobalAssemble();

    //     M_matrStiffPtr->spy("stiff");
    //     M_matrMassPtr->spy("mass");

    chrono.stop();
    M_displayer.leaderPrintMax("\n\t... done in " , chrono.diff() );;

    // record here the result of the matrices creation method
    M_whichMatrices = whichMatrices;

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
void ADRSolver<Mesh, SolverType>::
initialize( const function_type& u0 )
{
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );
    vector_type u(M_uFESpacePtr->map());
    M_uFESpacePtr->interpolate(u0, u, M_dataPtr->dataTimePtr()->getInitialTime());
    initialize( u );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
initialize( const vector_type& u0 )
{
    M_sol = u0;
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addMatrixTimeTerm( timeIntegrator_ptr_type timeIntegratorPtr )
{
    ASSERT( M_matrMassPtr.get(), "\nThe mass matrix is undefined!\n" );

    // update the matrix for time dependent problems
    M_displayer.leaderPrint( "  adr-  Adding the solution time derivative to the matrix ... ");
    *M_matrNoBCPtr += *M_matrMassPtr * timeIntegratorPtr->coeff_der_dt(0);
    M_displayer.leaderPrint( "done\n");
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
updateMatrix( timeIntegrator_ptr_type timeIntegratorPtr )
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );

    Chrono chrono;
    chrono.start();

    // recompute constant matrices, if needed
    if (M_recomputeConstantMatrix)
        computeConstantMatrices(M_whichMatrices);

    M_displayer.leaderPrint( "  adr-  Updating the full matrix ...\n");

    initializeMatrix();

    if( ( M_dataPtr->diffusionCoefficient() ) &&
            ( ( M_whichMatrices != ADR_ONLY_MASS ) || ( M_whichMatrices != ADR_NO_CONST_MATRICES ) ) )
    {
        M_displayer.leaderPrint( "\t");
        addDiffusionTerms();
    }

    if( ( M_dataPtr->reactionCoefficient() ) &&
            ( ( M_whichMatrices != ADR_ONLY_STIFF ) || ( M_whichMatrices != ADR_NO_CONST_MATRICES ) ) )
    {
        M_displayer.leaderPrint( "\t");
        addReactionTerms();
    }

    // managing the convective term
        M_displayer.leaderPrint( "\t");
        addAdvectionTerms();

    if( ( M_dataPtr->stabilizationMethod() != ADR_NO_STABILIZATION ) &&
            M_advectionVector.size() )
    {
        M_displayer.leaderPrint( "\t");
        addStabilizationTerms();
    }

    if( timeIntegratorPtr.get() )
    {
        M_displayer.leaderPrint( "\t");
        addMatrixTimeTerm( timeIntegratorPtr );
    }

    M_displayer.leaderPrintMax( "\t... done in " , chrono.diff() );

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addRHSTimeTerm( timeIntegrator_ptr_type timeIntegratorPtr )
{
    ASSERT( M_matrMassPtr.get(), "\nThe mass matrix is undefined!\n" );

    M_displayer.leaderPrint( "  adr-  Adding the solution time derivative to rhs ... ");
    // update the RHS for time dependent problems
    M_rhsNoBC += *M_matrMassPtr * timeIntegratorPtr->time_der_dt();
    M_displayer.leaderPrint( "done\n");
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
updateRHS( timeIntegrator_ptr_type timeIntegratorPtr )
{
    Chrono chrono;

    M_displayer.leaderPrint("  adr-  Updating the source term on right hand side ...\n");

    chrono.start();
    // Right hand side for the system
    M_displayer.leaderPrint( "\t" );
    assembleSourceTerm();
    //updateRHS(sourceVec);

    if( timeIntegratorPtr.get() )
    {
        M_displayer.leaderPrint( "\t" );
        addRHSTimeTerm( timeIntegratorPtr );
    }

    chrono.stop();
    M_displayer.leaderPrintMax("\t... done in ", chrono.diff());
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
updateRHS( vector_type const& rhsVec )
{
    ASSERT( M_matrMassPtr.get(), "\nThe mass matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint("  adr-  Updating the source term on right hand side ...\n");

    chrono.start();
    // Right hand side for the system
    M_displayer.leaderPrint( "\t" );

    M_rhsNoBC = *M_matrMassPtr * rhsVec;

    chrono.stop();
    M_displayer.leaderPrintMax("\t... done in ", chrono.diff());
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::iterate( bchandler_type& bch )
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_matrNoBCPtr.get(), "\nThe full matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint("  adr-  Finalizing the full matrix and rhs ... ");
    chrono.start();

    M_matrNoBCPtr->GlobalAssemble();
    if( M_dataPtr->stabilizationMethod() != ADR_NO_STABILIZATION )
    {
        ASSERT( M_matrStabPtr.get(), "\nThe stabilization matrix is undefined!\n" );
        M_matrStabPtr->GlobalAssemble();
    }
    M_rhsNoBC.GlobalAssemble();

    chrono.stop();
    M_displayer.leaderPrintMax("\n\t... done in ", chrono.diff() );

    M_matrNoBCPtr->spy("matrNoBC");

    // we need an empty matrix over which to copy the various terms
    M_displayer.leaderPrint("  adr-  Creating a new full matrix and rhs for BC prescription ...");
    chrono.start();

    matrix_ptr_type matrFull( new matrix_type( M_uFESpacePtr->map(), M_matrNoBCPtr->getMeanNumEntries()));

    *matrFull += *M_matrNoBCPtr;
    if( M_dataPtr->stabilizationMethod() != ADR_NO_STABILIZATION )
        *matrFull += *M_matrStabPtr;

    // the vector to be passed to BCManage is required to have a Unique map
    vector_type rhsFull(M_rhsNoBC, Unique);

    chrono.stop();
    M_displayer.leaderPrintMax("\n\t... done in ", chrono.diff() );

    // boundary conditions update
    M_displayer.leaderPrint("  adr-  Applying boundary conditions ...         ");
    chrono.start();

    M_displayer.leaderPrint("\n\t");
    applyBoundaryConditions( *matrFull, rhsFull, bch );

    M_displayer.leaderPrintMax("\t... done in " , chrono.diff());

    matrFull->spy("matrBC");

    M_displayer.leaderPrint("  adr-  Finalizing the new full matrix and rhs ...");
    chrono.start();

    matrFull->GlobalAssemble();
    rhsFull.GlobalAssemble();

    chrono.stop();
    M_displayer.leaderPrintMax("\n\t... done in " , chrono.diff());

    M_displayer.leaderPrint("  adr-  Invoking the linear solver ...\n");
    // solving the system
    M_linearSolver.setMatrix(*matrFull);

    // M_linearSolver.setReusePreconditioner( M_reusePrec );
    /*int numIter =*/ M_linearSolver.solveSystem( rhsFull, M_sol, matrFull );

    M_residual  = M_rhsNoBC;
    M_residual -= *M_matrNoBCPtr*M_sol;

} // iterate()


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_type&        matrix,
                                                           vector_type&        rhs,
                                                           bchandler_type& BCh )
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );

    // BC manage
    Chrono chrono;

    if ( !BCh.bdUpdateDone() )
    {
        M_displayer.leaderPrint( "  adr-  Updating the BC ... ");
        chrono.start();
        BCh.bdUpdate( *M_uFESpacePtr->mesh(), M_uFESpacePtr->feBd(), M_uFESpacePtr->dof() );
        chrono.stop();
        M_displayer.leaderPrintMax( "done in " , chrono.diff() );
        M_displayer.leaderPrint( "\t");
    }

    //    vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    // vector_type rhsFull(rhs, Unique); // ignoring non-local entries, Otherwise they are summed up lately


    M_displayer.leaderPrint( "  adr-  Managing the BC ... ");
    chrono.start();
    bcManage( matrix, rhs, *M_uFESpacePtr->mesh(), M_uFESpacePtr->dof(), BCh, M_uFESpacePtr->feBd(), 1.,
              M_dataPtr->dataTimePtr()->getTime() );
    chrono.stop();
    M_displayer.leaderPrintMax( "done in " , chrono.diff() );


    // rhs = rhsFull;


} // applyBoundaryCondition



template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
assembleSourceTerm()
{
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );
    ASSERT( M_elvecRhsPtr.get(), "\nThe elementary vector for the rhs is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint( "  adr-  Assemblying the source term on the rhs ... ");

    chrono.start();

    // Number of solution components
    UInt nc_u = this->M_uFESpacePtr->fieldDim();

    // Right hand side for the problem at time
    M_rhsNoBC *= 0.;
    // std::cout << "\nM_rhsNoBC.size = " << M_rhsNoBC.size() << std::endl;

    // loop on volumes: assembling source term
    for ( UInt iEl = 1; iEl <= M_uFESpacePtr->mesh()->numElements(); iEl++ )
    {
        this->M_uFESpacePtr->fe().updateJacQuadPt( this->M_uFESpacePtr->mesh()->element( iEl ) );

        for ( UInt ic = 0; ic < nc_u; ++ic )
        {
            M_elvecRhsPtr->zero();
            compute_vec( M_sourceTerm, *M_elvecRhsPtr, this->M_uFESpacePtr->fe(),
                         M_dataPtr->dataTimePtr()->getTime(), ic ); // compute local vector
            //            assemb_vec( M_rhsNoBC, *M_elvecRhsPtr, this->M_uFESpacePtr->fe(), this->M_uFESpacePtr->dof(), ic ); // assemble local vector into global one
            assembleVector( M_rhsNoBC,
                            this->M_uFESpacePtr->fe().currentLocalId(),
                            *M_elvecRhsPtr,
                            this->M_uFESpacePtr->refFE().nbDof(),
                            this->M_uFESpacePtr->dof(), ic );
        }
    }

    chrono.stop();
    M_displayer.leaderPrintMax( "done in " , chrono.diff() );

}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
initializeMatrix()
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the unknown is undefined!\n" );

    // initialize the system matrices
    if (M_matrNoBCPtr.get())
        M_matrNoBCPtr.reset(new matrix_type(M_uFESpacePtr->map(), M_matrNoBCPtr->getMeanNumEntries() ));
    else
//        M_matrNoBCPtr.reset(new matrix_type(M_uFESpacePtr->map(), M_matrStiffPtr->getMeanNumEntries() ));
        M_matrNoBCPtr.reset(new matrix_type(M_uFESpacePtr->map()));

    if( M_dataPtr->stabilizationMethod() != ADR_NO_STABILIZATION )
        M_matrStabPtr.reset( new matrix_type(M_uFESpacePtr->map()) );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addDiffusionTerms()
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_matrNoBCPtr.get(), "\nThe full matrix is undefined!\n" );
    ASSERT( M_matrStiffPtr.get(), "\nThe stiff matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint( "  adr-  Adding diffusion terms ... ");
    chrono.start();
    *M_matrNoBCPtr += *M_matrStiffPtr*M_dataPtr->diffusionCoefficient();
    chrono.stop();
    M_displayer.leaderPrintMax( "done in ", chrono.diff() );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addReactionTerms()
{
    ASSERT( M_dataPtr.get(), "\nThe data descriptor for ADR solver is undefined!\n" );
    ASSERT( M_matrNoBCPtr.get(), "\nThe full matrix is undefined!\n" );
    ASSERT( M_matrMassPtr.get(), "\nThe mass matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint( "  adr-  Adding reaction terms ... ");
    chrono.start();
    *M_matrNoBCPtr += *M_matrMassPtr*M_dataPtr->reactionCoefficient();
    chrono.stop();
    M_displayer.leaderPrintMax( "done in ", chrono.diff() );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addAdvectionTerms()
{
    ASSERT( M_advectionFieldFESpacePtr.get(), "\nThe FE space for the advection field is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the solution is undefined!\n" );
    ASSERT( M_elmatAdvPtr.get(), "\nThe elementary matrix for the advection terms is undefined!\n" );
    ASSERT( M_matrNoBCPtr.get(), "\nThe full matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint("  adr-  Updating the convective terms ... ");

    // vector with repeated nodes over the processors
    vector_type advectionVecRep( this->M_advectionFieldFESpacePtr->map(), Repeated );
    if( M_advectionVector.size() )
    {
        advectionVecRep = M_advectionVector;
    }

    chrono.start();

    for ( UInt iEl = 1; iEl <= M_uFESpacePtr->mesh()->numElements(); ++iEl )
    {

        M_elmatAdvPtr->zero();

        if( M_advectionVector.size() )
        {
            ASSERT( M_elvecBetaPtr.get(), "\nThe elementary vector for the advection terms is undefined!\n" );

            M_uFESpacePtr->fe().updateFirstDeriv( M_uFESpacePtr->mesh()->element( iEl ) ); //as updateFirstDer

            UInt eleID = M_uFESpacePtr->fe().currentLocalId();
            // Non linear term, Semi-implicit approach
            // M_elvec contains the velocity values in the nodes
            for ( UInt iNode = 0 ; iNode < M_uFESpacePtr->fe().nbFEDof() ; iNode++ )
            {
                // UInt  iloc = betaFESpacePtr->fe().patternFirst( iNode );
                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = M_uFESpacePtr->dof().localToGlobal( eleID, iNode + 1 ) + iComp*nDimensions;
                    M_elvecBetaPtr->vec()[ iNode + iComp*M_uFESpacePtr->fe().nbFEDof() ] = advectionVecRep[ig]; // BASEINDEX + 1
                }
            }

            // compute local convective terms
            advection( 1., *M_elvecBetaPtr, *M_elmatAdvPtr, M_uFESpacePtr->fe(), 0, 0, M_uFESpacePtr->fieldDim() );
        }
        else
        {
            // Update geometrical info on the current FE
            this->M_uFESpacePtr->fe().updateFirstDerivQuadPt( this->M_uFESpacePtr->mesh()->element( iEl ) );

            // compute local convective terms
            advection( 1., M_advectionField, *M_elmatAdvPtr, this->M_uFESpacePtr->fe(), 0, 0,
                       M_uFESpacePtr->fieldDim(), M_dataPtr->dataTimePtr()->getTime() );
        }

        // assembly: loop on components
        for ( UInt iComp = 0; iComp < M_uFESpacePtr->fieldDim(); ++iComp )
        {
            assembleMatrix( *M_matrNoBCPtr,
                            *M_elmatAdvPtr,
                            M_uFESpacePtr->fe(),
                            M_uFESpacePtr->fe(),
                            M_uFESpacePtr->dof(),
                            M_uFESpacePtr->dof(),
                            iComp, iComp,
                            iComp*M_uFESpacePtr->dim(), iComp*M_uFESpacePtr->dim() );
        }

    }
    chrono.stop();
    M_displayer.leaderPrintMax( "done in " , chrono.diff() );
}


template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
addStabilizationTerms()
{
    ASSERT( M_advectionFieldFESpacePtr.get(), "\nThe FE space for the advection field is undefined!\n" );
    ASSERT( M_uFESpacePtr.get(), "\nThe FE space for the solution is undefined!\n" );
    ASSERT( M_elmatStabPtr.get(), "\nThe elementary matrix for the stabilization terms is undefined!\n" );
    ASSERT( M_matrNoBCPtr.get(), "\nThe full matrix is undefined!\n" );

    Chrono chrono;

    M_displayer.leaderPrint("  adr-  Updating the stabilization terms ...        ");

    // vector with repeated nodes over the processors
    vector_type advectionVecRep( this->M_advectionFieldFESpacePtr->map(), Repeated );
    if( M_advectionVector.size() )
    {
        advectionVecRep = M_advectionVector;
    }

    chrono.start();

    for ( UInt iEl = 1; iEl <= M_uFESpacePtr->mesh()->numElements(); ++iEl )
    {

        UInt eleID = M_uFESpacePtr->fe().currentLocalId();

        if( !M_advectionVector.size() )
            this->M_uFESpacePtr->fe().update( M_uFESpacePtr->mesh()->element( iEl ),
                                              UPDATE_ONLY_CELL_NODES );

        // Streamline Diffusion ( only stab for now )
        if( M_dataPtr->stabilizationMethod() == ADR_SD_STABILIZATION )
        {
            M_elmatStabPtr->zero();

            Real advectionFieldValue = 0.;
            Real VLoc_infty = 0.;
            Real VLoc_mean  = 0.;
            Real VLoc_c     = 0.;

            for ( UInt ih_c = 0 ; ih_c < this->M_uFESpacePtr->fe().nbFEDof() ; ih_c++ )
            {
                UInt iloc = this->M_uFESpacePtr->fe().patternFirst( ih_c );

                for ( UInt iComp = 0; iComp < nDimensions; ++iComp)
                {
                    if( M_advectionVector.size() )
                    {
                        UInt ig = M_uFESpacePtr->dof().localToGlobal( eleID, iloc + 1 ) + iComp*M_uFESpacePtr->dim();
                        advectionFieldValue = advectionVecRep[ ig ];
                    }
                    else
                    {
                        Real x( this->M_uFESpacePtr->fe().cellNode( ih_c, 0 ) );
                        Real y( this->M_uFESpacePtr->fe().cellNode( ih_c, 1 ) );
                        Real z( this->M_uFESpacePtr->fe().cellNode( ih_c, 2 ) );
                        advectionFieldValue = M_advectionField( M_dataPtr->dataTimePtr()->getTime(), x, y, z, iComp );
                    }
                    M_elvecBetaPtr->vec()[ iloc + iComp * this->M_uFESpacePtr->fe().nbFEDof() ] = advectionFieldValue;
                    VLoc_c += advectionFieldValue * advectionFieldValue;
                }

                VLoc_c     = sqrt( VLoc_c );
                VLoc_mean += VLoc_c;

                if ( VLoc_c > VLoc_infty )
                    VLoc_infty = VLoc_c;
            }

            VLoc_mean = VLoc_mean / this->M_uFESpacePtr->fe().nbFEDof();

            Real coef_stab = 0;
            coef_stab=M_dataPtr->stabilizationCoefficient()*this->M_uFESpacePtr->fe().diameter()*VLoc_infty; // Alessandro - method

            stiff_sd( coef_stab / ( VLoc_mean*VLoc_mean ), *M_elvecBetaPtr, *M_elmatStabPtr, this->M_uFESpacePtr->fe(), this->M_uFESpacePtr->fe() );

            assembleMatrix( *M_matrStabPtr,
                            *M_elmatStabPtr,
                            M_uFESpacePtr->fe(),
                            M_uFESpacePtr->fe(),
                            M_uFESpacePtr->dof(),
                            M_uFESpacePtr->dof(),
                            0, 0, 0, 0 );
        }


        /*

         //TODO: check both the 2D and 3D implementation of ip-stabilization (on the Laplacian test doesn't work properly)
        if( M_dataPtr->stabilizationMethod() != ADR_IP_STABILIZATION )
        {
            //        M_displayer.leaderPrint("   adr- IP stab");
            //            if ( M_resetStab )
            // {
            const UInt nDof = M_advectionFieldFESpacePtr->dof().numTotalDof();

            CurrentFE fe1(M_uFESpacePtr->refFE(),
                          getGeoMap(*M_uFESpacePtr->mesh()),
                          M_uFESpacePtr->qr());
            CurrentFE fe2(M_uFESpacePtr->refFE(),
                          getGeoMap(*M_uFESpacePtr->mesh()),
                          M_uFESpacePtr->qr());
            CurrentFE fe3(M_advectionFieldFESpacePtr->refFE(),
                          getGeoMap(*M_advectionFieldFESpacePtr->mesh()),
                          M_advectionFieldFESpacePtr->qr());


#ifdef TWODIM
            typedef ID ( *ETOP )( ID const localFace, ID const point );
            ETOP  eToP;
            switch( M_uFESpacePtr->fe().refFE.type )
            {
            case FE_P1_2D:
                eToP = LinearTriangle::eToP;
                break;
            case FE_P2_2D:
                eToP = QuadraticTriangle::eToP;
                break;
            case FE_Q1_2D:
                eToP = LinearQuad::eToP;
                break;
            case FE_Q2_2D:
                eToP = QuadraticQuad::eToP;
                break;
            default:
                eToP=0;
                ERROR_MSG( "This refFE is not allowed with IP stabilization" );
                break;
            }
            for ( UInt iEdge = M_uFESpacePtr->mesh()->numBEdges() + 1; iEdge <= M_uFESpacePtr->mesh()->numEdges();
                    ++iEdge )
            {
                const UInt iElAd1 = M_uFESpacePtr->mesh()->edge( iEdge ).ad_first();
                const UInt iElAd2 = M_uFESpacePtr->mesh()->edge( iEdge ).ad_second();

                if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
                {
                    continue;
                }

         *M_elmatStabPtr.zero();

                M_advectionFieldFESpacePtr->feBd().updateMeas( M_advectionFieldFESpacePtr->mesh()->edge( iEdge ) );
                const Real hK2  = std::pow(M_advectionFieldFESpacePtr->feBd().measure(), 2.);

                M_advectionFieldFESpacePtr->feBd().updateMeasNormal( M_advectionFieldFESpacePtr->mesh()->edge( iEdge ) );
                KNM<Real>& normal = M_advectionFieldFESpacePtr->feBd().normal;

                fe1.updateFirstDeriv( M_uFESpacePtr->mesh()->element( iElAd1 ) );
                fe2.updateFirstDeriv( M_uFESpacePtr->mesh()->element( iElAd2 ) );

                ElemVec beta(M_advectionFieldFESpacePtr->feBd().nbFEDof(), nDimensions);

                // first, get the local trace of the velocity into beta
                // local id of the face in its adjacent element

                UInt iEdEl = M_advectionFieldFESpacePtr->mesh()->edge( iEdge ).pos_first();
                for ( int iNode = 0; iNode < M_advectionFieldFESpacePtr->feBd().nbFEDof(); ++iNode )
                {
                    UInt iloc = eToP( iEdEl, iNode+1 );
                    for ( int iCoor = 0; iCoor < fe1.nbCoor(); ++iCoor )
                    {
                        UInt ig = M_advectionFieldFESpacePtr->dof().localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                        if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
                            beta.vec()[ iCoor*M_advectionFieldFESpacePtr->feBd().nbFEDof() + iNode ] = betaVecRep( ig + 1); // BASEINDEX + 1
                    }
                }
            }
#elif defined THREEDIM
            typedef ID ( *FTOP )( ID const localFace, ID const point );
            FTOP  fToP;
            switch( M_uFESpacePtr->fe().refFE().type() )
            {
            case FE_P1_3D:
            case FE_P1bubble_3D:
                fToP = LinearTetra::fToP;
                break;
            case FE_P2_3D:
                fToP = QuadraticTetra::fToP;
                break;
            case FE_Q1_3D:
                fToP = LinearHexa::fToP;
                break;
            case FE_Q2_3D:
                fToP = QuadraticHexa::fToP;
                break;
            default:
                fToP = 0;
                ERROR_MSG( "This refFE is not allowed with IP stabilisation" );
                break;
            }


            for ( UInt iFace = M_uFESpacePtr->mesh()->numBFaces() + 1; iFace <= M_uFESpacePtr->mesh()->numFaces();
                    ++iFace )
            {

                const UInt iElAd1 = M_uFESpacePtr->mesh()->face( iFace ).ad_first();
                const UInt iElAd2 = M_uFESpacePtr->mesh()->face( iFace ).ad_second();

                if ( iElAd1 == iElAd2 || iElAd1 == 0 || iElAd2 == 0)
                {
                    continue;
                }

                M_elmatStabPtr->zero();

                M_advectionFieldFESpacePtr->feBd().updateMeas( M_advectionFieldFESpacePtr->mesh()->face( iFace ) );
                const Real hK2  = std::pow(M_advectionFieldFESpacePtr->feBd().measure(), 2.);

                M_advectionFieldFESpacePtr->feBd().updateMeasNormal( M_advectionFieldFESpacePtr->mesh()->face( iFace ) );
                KNM<Real>& normal = M_advectionFieldFESpacePtr->feBd().normal;

                fe1.updateFirstDeriv( M_uFESpacePtr->mesh()->element( iElAd1 ) );
                fe2.updateFirstDeriv( M_uFESpacePtr->mesh()->element( iElAd2 ) );

                Real bn   = 0;
                Real bmax = 0;
         */
        // Old version, removed by SQ
        /*
                    ElemVec beta( M_advectionFieldFESpacePtr->feBd().nbFEDof(), nDimensions);

                    // first, get the local trace of the velocity into beta
                    // local id of the face in its adjacent element

                    UInt iFaEl = M_advectionFieldFESpacePtr->mesh()->face( iFace ).pos_first();
                    for ( int iNode = 0; iNode < M_advectionFieldFESpacePtr->feBd().nbFEDof(); ++iNode )
                    {
                        UInt iloc = fToP( iFaEl, iNode+1 );
                        for ( int iCoor = 0; iCoor < fe1.nbCoor(); ++iCoor )
                        {
                            UInt ig = M_advectionFieldFESpacePtr->dof().localToGlobal( iElAd1, iloc + 1 ) - 1 +iCoor*nDof;
                            if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
                  {
                Real value( betaVecRep( ig + 1) );
                                beta.vec()[ iCoor*M_advectionFieldFESpacePtr->feBd().nbFEDof() + iNode ] = value; // BASEINDEX + 1
                  };
                        }
            }*/
        /*
                // New version, added by SQ : beta on the domain, not the boundary!
                // See elemOper.cpp for justification of this usage
                ElemVec beta( M_advectionFieldFESpacePtr->fe().nbFEDof(), nDimensions);

                for ( UInt iNode = 0; iNode < M_advectionFieldFESpacePtr->fe().nbFEDof(); ++iNode )
                {
                    UInt  iloc = M_advectionFieldFESpacePtr->fe().patternFirst( iNode );
                    for ( UInt iCoor = 0; iCoor < fe1.nbCoor(); ++iCoor )
                    {
                        UInt ig = M_advectionFieldFESpacePtr->dof().localToGlobal( iElAd1, iloc + 1 ) + iCoor*nDof;
                        beta.vec()[ iloc + iCoor*M_advectionFieldFESpacePtr->fe().nbFEDof() ] = betaVecRep[ig]; // BASEINDEX + 1

                    }
                }

                // second, calculate its max norm
                for ( int l = 0; l < int( M_advectionFieldFESpacePtr->fe().nbCoor()*M_advectionFieldFESpacePtr->fe().nbFEDof() ); ++l ) // SQ: feBd->fe
                {
                    if ( bmax < fabs( beta.vec()[ l ] ) )
                        bmax = fabs( beta.vec()[ l ] );
                }

                UInt iFaEl = M_advectionFieldFESpacePtr->mesh()->face( iFace ).pos_first();
                for ( int iNode = 0; iNode < M_advectionFieldFESpacePtr->feBd().nbNode; ++iNode )
                {
                    UInt iloc = fToP( iFaEl, iNode + 1 );
                    for ( UInt iCoor = 0; iCoor < nDimensions; ++iCoor )
                    {
                        UInt ig = M_advectionFieldFESpacePtr->dof().localToGlobal( iElAd1, iloc + 1 ) - 1 + iCoor*nDof;
                        if (betaVecRep.BlockMap().LID(ig + 1) >= 0)
                            bn += normal(iNode, (int)iCoor)*betaVecRep( ig + 1 );
                    }
                }

                Real coeffBeta = hK2*M_dataPtr->stabilizationCoefficient()*abs(bn);


#endif
                //                    Real coeffBeta = M_dataPtr->stabilizationCoefficient();
                //                     std::cout << coeffBeta << std::endl;

                ipstab_bagrad( coeffBeta,
         *M_elmatStabPtr,
                               fe1,
                               fe2,
                               fe3,
                               beta,
                               M_uFESpacePtr->feBd(),
                               0, 0);

                assembleMatrix( *M_matrStabPtr,
         *M_elmatStabPtr,
                                fe1,
                                fe2,
                                M_uFESpacePtr->dof(),
                                M_uFESpacePtr->dof(),
                                0, 0, 0, 0 );

            }
        }

         */
        chrono.stop();
        M_displayer.leaderPrintMax( "done in " , chrono.diff() );
    }

}

template<typename Mesh, typename SolverType>
void ADRSolver<Mesh, SolverType>::
setAdvectionField( const vector_type& betaVec,
                   fespace_ptr_type advectionFieldFESpacePtr )
{
    ASSERT( M_advectionFieldFESpacePtr.get(), "\nThe FE space for the advection field is undefined!\n" );

    // if needed, interpolate over the solution FE space
    vector_type betaVecInterpolated( this->M_advectionFieldFESpacePtr->map(), Unique );

    // interpolate over the solution FE
    if( advectionFieldFESpacePtr )
    {
        betaVecInterpolated = M_advectionFieldFESpacePtr->FeToFeInterpolate( *advectionFieldFESpacePtr, betaVec );
        M_advectionVector = betaVecInterpolated;
    }
    else
        M_advectionVector = betaVec;
}


} // namespace LifeV


#endif //_ADR_H_
