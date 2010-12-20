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
 *  @file
 *  @brief This file contains solver for St. Venant-Kirchhof materials.
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @version 1.1
 *  @date 01-03-2010
 *  @author Gilles Fourestey <gilles.fourestey@cscs.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  A more detailed description of the file (if necessary)
 */

#ifndef _NLVENANTKIRCHHOFSOLVER_H_
#define _NLVENANTKIRCHHOFSOLVER_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/AssemblyElemental.hpp>
//#include <life/lifefem/values.hpp>
//#include <life/lifearray/pattern.hpp>
#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifecore/chrono.hpp>

#include <life/lifealg/nonLinRichardson.hpp>
#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>

#include <life/lifecore/displayer.hpp>

namespace LifeV
{

#define nonlinear

/*!
  \class NonLinearVenantKirchhofSolver
  \brief
  This class solves the linear elastodynamics equations for a  St. Venant-Kirchoff material.

*/
template <typename Mesh, typename SolverType = LifeV::SolverTrilinos >
class NonLinearVenantKirchhofSolver : public VenantKirchhofSolver<Mesh, SolverType>
{
public:

    //! @name Type definition
    //@{

    typedef VenantKirchhofSolver<Mesh, SolverType> super;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_Type;

    typedef typename super::bchandler_Type                  bchandler_Type;

    typedef SolverType                                      solver_Type;

    typedef typename solver_Type::matrix_type               matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                  matrixPtr_Type;
    typedef typename solver_Type::vector_type               vector_Type;
    typedef boost::shared_ptr<vector_Type>                  vectorPtr_Type;


    typedef typename SolverType::prec_raw_type              precRaw_Type;
    typedef typename SolverType::prec_type                  prec_Type;

    typedef DataElasticStructure                            data_Type;

    //@}


    //! @name Constructor
    //@{

    NonLinearVenantKirchhofSolver();

    //@}

    //! @name Methods
    //@{

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param comm the comunicator parameter
      \param monolithicMap the EpetraMap
      \param offset the offset parameter
    */
    void setup( boost::shared_ptr<data_Type> data,
                const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
                boost::shared_ptr<Epetra_Comm>&     comm,
                const boost::shared_ptr<const EpetraMap>&   monolithicMap,
                UInt       offset=0
		);
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param comm
    */
    void setup(boost::shared_ptr<data_Type>        data,
	       const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
	       boost::shared_ptr<Epetra_Comm>&     comm
	       );

    //! Computes the Linear part of the Jacobian and the Mass Matrix
    /*!

      \param matrix the matrix containing the mass and the linear part of the stiffness matrix
      \param factor scaling factor
     */
    void buildSystem(matrixPtr_Type matrix, Real const& factor =1.);

    //! Computes the Linear part of the Jacobian and the Mass Matrix
    void buildSystem( );

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 space function describing the initial velocity
      \param a0 space function describing the initial acceleration
    */
    void initialize   ( const Function& d0, const Function& w0, const Function& a0 );

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param w0 space function describing the initial velocity
    */
    void initialize   ( vectorPtr_Type d0,  vectorPtr_Type w0 = vectorPtr_Type());

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 empty vector
      \param a0 empty vector
    */
    void initializeVel( const vector_Type& w0);

    //! Updates the system at the end of each time step when the matrix is passed from outside
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateSystem( matrixPtr_Type& stiff );

    //! Updates the system at the end of each time step when the matrix is passed from outside
    void updateSystem( );

    //! Update the right  hand side  for time advancing given a source term
    /*!
      \param source volumic source
      \param time present time
    */
    void updateSystem( source_Type const& source, Real t );

    //! Updates the nonlinear terms in the matrix at the end of each Newton iteration
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateNonlinearMatrix( matrixPtr_Type& stiff );

    //! Updates the nonlinear terms on the stiffness matrix at the end of each time step
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateNonlinearTerms( matrixPtr_Type& stiff );

    //! Solve the non-linear system

    /*! They compute the solution solving the non linear system

      \param bch BCHandler object with the applied boundary conditions

    */
    void iterate( bchandler_Type& bch );

    //!It computes the velocity vector at the n-th time step.
    void updateVel();

    //! computes the global matrix for the residual computation (calls updateNonlinearMatrix)
    /*!
      \param stiff stiffness matrix provided from outside
      \param solution vector
      \param factor scaling factor
    */
    void computeMatrix( matrixPtr_Type& stiff, const vector_Type& sol, Real const& factor );

    //! computes the global matrix for the residual computation (calls updateNonlinearMatrix)
    /*!
      \param solution vector
      \param factor scaling factor
    */
    void computeMatrix( const vector_Type& sol, Real const& factor );

    //! evaluates residual for newton interations
    /*!
      \param res residal vector that is update every time the method is called
      \param sol solution vector from which the residual is computed
      \param iter iteration of the nonLinearRichardson method
    */
    void evalResidual( vector_Type &res, const vector_Type& sol, Int iter);

    //! Update Jacobian at each nonLinearRichardson iteration
    /*!
      \param sol the current solution at the k-th iteration of Newton method
      \param jac the Jacobian matrix that must be updated
    */
    void updateJacobian( vector_Type& sol, matrixPtr_Type& jacobian  );

    //! It is called by the NonLinear Richarson method. It calls the solveJacobian method
    /*!

      \param step the current step at the k-th iteration
      \param res the residual at the k-th iteration of the nonlinear Richardson method
      \param linear_rel_tol send for the relative tolerance to the linear solver is therefore eta. eta is determined
             by the modified Eisenstat-Walker formula
     */
    void solveJac( vector_Type&       step,
                   const vector_Type& res,
                   Real&            linear_rel_tol);

    //! solves the tangent problem with custom BC

    //!  It solves the tangent problem with custom BC. It solves the linear system J*step=-Res at each k-th iteration. The first three arguments are the same passed to sovleJac
    /*!
      \param BCd BCHandler object
    */
    void solveJacobian( vector_Type&       step,
                        const vector_Type& res,
                        Real&            linear_rel_tol,
                        bchandler_Type&    BCd);

    //@}

    //! @name Get Methods
    //@{

    //! Get the acceleration
    vector_Type& getAcceleration()  { return *M_acc; }

    //! Get the solidMatrix
    void getSolidMatrix( matrixPtr_Type& matrix);

    //@}

private:

    //!It applies the boundary conditions before solving the linear system J*step=-Res. It is called in solveJacobian
    void applyBoundaryConditions(matrix_Type &matrix,
                                 vector_Type &rhs,
                                 bchandler_Type& BCh,
                                 UInt         offset=0);

    //! time scheme coefficients
    // _theta and _zeta are the coefficient of the time discretization with Newmark scheme
    Real                            M_zeta;
    Real                            M_theta;

    Real                            M_pressure;

    //! acceleration
    vectorPtr_Type                    M_acc;


    //! right  hand  side acceleration
    vectorPtr_Type                    M_rhsA;

};

//============================================
// Implementations
//=============================================

//=============================================
// Constructor
//=============================================

template <typename Mesh, typename SolverType>
NonLinearVenantKirchhofSolver<Mesh, SolverType>::
NonLinearVenantKirchhofSolver( ) :
        super                        ( ),
        M_zeta                       ( 0.75 ),
        M_theta                      ( 0.7 ),
        M_acc			 ( ),
        M_rhsA			 ( )
{
    //    M_BCh->setOffset(M_offset);
}

//=============================================
// Methods
//=============================================


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
setup(
      boost::shared_ptr<data_Type>        data,
      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
      boost::shared_ptr<Epetra_Comm>&     comm,
      const boost::shared_ptr<const EpetraMap>&  monolithicMap,
      UInt                                offset
      )
{
    super::setup(data, dFESpace, comm, monolithicMap, offset);

    M_acc.reset(new vector_Type(*this->M_localMap));
    M_rhsA.reset(new vector_Type(*this->M_localMap));

}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
setup(
      boost::shared_ptr<data_Type>        data,
      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
      boost::shared_ptr<Epetra_Comm>&     comm
      )
{
    super::setup(data, dFESpace, comm);
}

template <typename Mesh, typename SolverType>
void
NonLinearVenantKirchhofSolver<Mesh, SolverType>::
buildSystem(matrixPtr_Type massStiff, Real const & factor)
{
    UInt totalDof = this->M_FESpace->dof().numTotalDof();

    this->M_Displayer->leaderPrint( "NonLin S-  Building the system             ... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    UInt nc = nDimensions;

    //inverse of dt:
    Real PTemp =  this->M_data->getDataTime()->timeStep();
    Real dti2  = 2.0 / ( PTemp * PTemp );

    // Elementary computation and matrix assembling
    // Loop on elements
    // These Matrices are assembled according to the Newmark Scheme!
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        Int marker    = this->M_FESpace->mesh()->volumeList( i ).marker();

        Real lambda = this->M_data->getLambda(marker);
        Real mu     = this->M_data->getMu    (marker);

        this->M_elmatK->zero();
        this->M_elmatM->zero();


        // stiffness
        stiff_strain (        2.0*mu, *this->M_elmatK, this->M_FESpace->fe() );

        stiff_div   ( lambda, *this->M_elmatK, this->M_FESpace->fe() );

        this->M_elmatC->mat() = this->M_elmatK->mat();

        // mass
        mass( dti2 * this->M_data->getRho(), *this->M_elmatM, this->M_FESpace->fe(), 0, 0, nDimensions );

        this->M_elmatC->mat() += this->M_elmatM->mat();

        // assembling
        for ( UInt ic = 0; ic < nc; ic++ )
        {
            for ( UInt jc = 0; jc < nc; jc++ )
            {
                assembleMatrix( *this->M_linearStiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->fe(), this->M_FESpace->dof(), this->M_FESpace->dof(),  ic,  jc,  this->M_offset +ic*totalDof, this->M_offset + jc*totalDof );

                assembleMatrix( *massStiff  , *this->M_elmatC, this->M_FESpace->fe(), this->M_FESpace->fe(), this->M_FESpace->dof(), this->M_FESpace->dof(),  ic,jc, this->M_offset +ic*totalDof, this->M_offset + jc*totalDof );
            }

            //mass
            assembleMatrix( *this->M_mass, *this->M_elmatM, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, ic, this->M_offset +  ic*totalDof, this->M_offset +  ic*totalDof);
        }
    }

    this->M_linearStiff->globalAssemble();
    massStiff->globalAssemble();
    //*M_massStiff *= factor; //Used in monolithic
    this->M_mass->globalAssemble();

    chrono.stop();

    this->M_Displayer->leaderPrintMax( " done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
buildSystem( )
{
    super::buildSystem( );
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
initialize( const Function& d0, const Function& w0, const Function& a0 )
{

    this->M_FESpace->interpolate(d0, *this->M_disp, 0.0);
    this->M_FESpace->interpolate(w0, *this->M_vel , 0.0);
    this->M_FESpace->interpolate(a0, *this->M_vel , 0.0);
}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem( matrixPtr_Type& stiff )
{
    this->M_Displayer->leaderPrint(" NonLin S-  Updating mass term on right hand side... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    //    UInt nc = nDimensions;
    Real coef;

    // start of the non linear part

#ifdef nonlinear


    stiff.reset(new matrix_Type(*this->M_localMap));

    matrixPtr_Type tmp( new matrix_Type(*this->M_localMap, 1) );
    updateNonlinearTerms(tmp);
    tmp->globalAssemble();
    *stiff += *tmp;
    *stiff += *this->M_linearStiff;

    stiff->globalAssemble();

    //_rhsContributionSecondDerivativeithoutBC -= _K * this->_d;

#endif

    // end of the nonlinear part
    //Computation of the right hand sides

    Real DeltaT    = this->M_data->getDataTime()->timeStep();
    vector_Type z = *this->M_disp;

     z            +=  DeltaT*(*this->M_vel);

    coef= (1.0-this->M_zeta);

    *this->M_rhsNoBC  = *this->M_mass * z;

    *this->M_rhsNoBC -= (*stiff) * coef * (*this->M_disp);

    // acceleration rhs
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);

    // velocity rhs
    *this->M_rhsW = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);

    std::cout << std::endl;

    std::cout << "rhsContributionSecondDerivativeithoutBC norm = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "rhs_w norm        = " << this->M_rhsW->norm2() << std::endl;
    std::cout << "    w norm        = " << this->M_vel->norm2() << std::endl;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());


}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem(  )
{
    updateSystem(this->M_stiff);
}



template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem(  source_Type const& source, Real t )
{

    this->M_Displayer->leaderPrint(" NonLin S-  Updating mass term on right hand side... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    //    UInt nc = nDimensions;

    Real coef;

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    *this->M_rhsNoBC *= 0.0;

    // start of the non linear part

#ifdef nonlinear
    ElemVec dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type disp(*this->M_disp);

    vector_Type dRep(disp, Repeated);

    this->M_stiff.reset(new matrix_Type(*this->M_localMap));

    matrixPtr_Type tmp(new matrix_Type(*this->M_localMap, 1));
    updateNonlinearTerms( tmp );
    //updateNonlinearMatrix( tmp );
    tmp->GlobalAssemble();
    *this->M_stiff += *tmp;

    *this->M_stiff += *this->M_linearStiff;

    this->M_stiff->GlobalAssemble();

    //_rhsContributionSecondDerivativeithoutBC -= _K * this->_d;


#endif
    // End of the nonlinear part

    ElemVec M_elvec(this->M_FESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

// loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec( source, M_elvec, this->M_FESpace->fe(),  this->M_data->getDataTime()->getTime(), ic ); // compute local vector
            assembleVector( *this->M_rhsNoBC, M_elvec, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, ic*this->M_FESpace->dim() ); // assemble local vector into global one
        }
    }

    this->M_rhsNoBC->GlobalAssemble();

    //Computation of the right hand sides

    Real DeltaT    = this->M_data->getDataTime()->getTimeStep();
    vector_Type z  = *this->M_disp;

    z             +=  DeltaT * (*this->M_vel);

    coef= ( 1.0-M_zeta );

    *this->M_rhsNoBC += *this->M_mass * z;

    *this->M_rhsNoBC -= (*this->M_stiff) * coef * (*this->M_disp);

    // acceleration rhs
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);

    // velocity rhs
    *this->M_rhsContributionSecondDerivative = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);

    std::cout << std::endl;

    std::cout << "rhsContributionSecondDerivativeithoutBC norm = " << this->M_rhsNoBC->Norm2() << std::endl;
    std::cout << "rhs_w norm       = " << this->M_rhsContributionSecondDerivative->Norm2() << std::endl;
    std::cout << "    w norm       = " << this->M_vel->Norm2() << std::endl;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());


}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateNonlinearMatrix( matrixPtr_Type& stiff )
{
    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    ElemVec dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type disp(*this->M_disp);

    vector_Type dRep(disp, Repeated);

    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim();
                dk_loc[ iloc + iComp*this->M_FESpace->fe().nbFEDof() ] = dRep[ig]; // BASEINDEX + 1
            }
        }


        this->M_elmatK->zero();

        //non-linear terms of the stiffness matrix
	// the coefficient M_zeta is passed to the matrix in order to have the right
	// assembling of the matrix. Another solution could be export the moltiplication
	// of the nonlinear part by M_zeta outside. In this case, the moltiplication
        // must be done after the GlobalAssemble() and NOT BEFORE!!

        // 3) 1/2 * \lambda * ( \tr { [\grad d^k]^T \grad d }, \div v  )
        stiff_derdiv( this->M_zeta * 0.5 * this->M_data->getLambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //4)  \mu * ( [\grad d^k]^T \grad d : \grad v  )
        stiff_dergradbis( this->M_zeta* this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //  5):\lambda * (div u_k) \grad d : \grad v
        stiff_divgrad( this->M_zeta * this->M_data->getLambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 6)  \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
        stiff_gradgrad( this->M_zeta *  0.5 * this->M_data->getLambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
        stiff_dergrad_gradbis( this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
        stiff_dergrad_gradbis_Tr( this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 8) // \mu * (  \grad d^k [\grad d^k]^T \grad d : \grad v  )
        stiff_gradgradTr_gradbis( this->M_zeta * this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        UInt totalDof   = this->M_FESpace->dof().numTotalDof();

        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            // stiff is the nonlinear matrix of the bilinear form in the weak formulation
            for ( UInt jc = 0; jc < nDimensions; jc++ )
                assembleMatrix( *stiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);

        }
    }
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateNonlinearTerms( matrixPtr_Type& stiff )
{

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    ElemVec dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );
    vector_Type disp(*this->M_disp);

    vector_Type dRep(disp, Repeated);

    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );


        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim() + this->M_offset;
                dk_loc[ iloc + iComp*this->M_FESpace->fe().nbFEDof() ] = dRep[ig]; // BASEINDEX + 1
            }
        }

        this->M_elmatK->zero();

        // non-linear terms of the stiffness matrix

        // 3) 1/2 * \lambda  ( \tr { [\grad d^k]^T \grad d }, \div v  )
        stiff_derdiv( 0.5 * this->M_data->getLambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //4)  \mu *( [\grad d^k]^T \grad d : \grad v  )
        stiff_dergradbis( this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //  5): \lambda * (div u_k) \grad d : \grad v
        stiff_divgrad( this->M_data->getLambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 6) 1/2  * \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
        stiff_gradgrad( 0.5 * this->M_data->getLambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
        stiff_dergrad_gradbis( this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
        // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
        stiff_dergrad_gradbis_Tr( this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 8) // \mu *  (  \grad d^k [\grad d^k]^T \grad d : \grad v  )
        stiff_gradgradTr_gradbis( this->M_data->getMu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            // stiff is the nonlinear matrix of the bilinear form
            for ( UInt jc = 0; jc < nDimensions; jc++ )

                assembleMatrix( *stiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);
        }

    }
}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
iterate( bchandler_Type& bch )
{
    Chrono chrono;

    // matrix and vector assembling communication
    this->M_Displayer->leaderPrint("  NonLin S-  Solving the system ... \n");

    this->M_BCh = bch;

    Real abstol  = 1.e-5;
    Real reltol  = 1.e-5;
    UInt   maxiter = 200;
    Real etamax  = 0;
    Int linesearch = 0;

    Real time = this->M_data->time();

    Int status = 0;

    status = nonLinRichardson( *this->M_disp, *this, abstol, reltol, maxiter, etamax, linesearch, this->M_out_res, this->M_data->getDataTime()->time() );

    if ( status == 1 )
    {
        std::ostringstream ex;
        ex << "VenantKirchhofSolver::iterate() Inners nonLinearRichardson iterations failed to converge\n";
        throw std::logic_error( ex.str() );
    }
    else // if status == 0 NonLinearrRichardson converges
    {
        std::cout << std::endl;

        std::cout <<" Number of inner iterations       : " << maxiter <<  std::endl;

        std::cout <<" We are at the time step          : "  << this->M_data->getDataTime()->time() << std::endl;

        this->M_out_iter << time << " " << maxiter << std::endl;
    }

    updateVel();

    std::cout << "iterate: d norm       = " << this->M_disp->norm2() << std::endl;
    std::cout << "iterate: w norm       = " << this->M_vel->norm2() << std::endl;
    std::cout << "iterate: a norm       = " << M_acc->norm2() << std::endl;

    *this->M_residual_d = *this->M_mass*(*this->M_disp);
    *this->M_residual_d -= *this->M_rhsNoBC;

    //This part should be checked with Non Linear Structure!!
    //    *M_residual_d  = *M_massStiff*(M_disp);
    //    *M_residual_d -= *M_rhsNoBC;

} // iterate()


template <typename Mesh, typename SolverType> // for monolithic
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
updateVel()
{
    Real DeltaT = this->M_data->getDataTime()->timeStep();

    *M_acc = (2.0 /( M_zeta * pow(DeltaT,2) ))  * (*this->M_disp)  - *M_rhsA;
    *this->M_vel = *this->M_rhsW + M_theta * DeltaT * (*M_acc) ;
}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::computeMatrix( matrixPtr_Type& stiff, const vector_Type& sol,  Real const& factor)
{
    this->M_Displayer->leaderPrint( "    NonLin S- Computing residual ... \t\t\t");

    Chrono chrono;
    chrono.start();

    // Matrices initialization

    stiff.reset(new matrix_Type(*this->M_localMap));
    *stiff += *this->M_linearStiff;

    Real coef;
    coef=M_zeta;
    *stiff *= coef;


    ElemVec dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type dRep(sol, Repeated);


#ifdef nonlinear

    matrixPtr_Type tmp(new matrix_Type(*this->M_localMap, 1));
    updateNonlinearMatrix( tmp );
    tmp->globalAssemble();
    *tmp *= factor;
    *stiff += *tmp;

#endif

    *stiff += *this->M_mass;
    stiff->globalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );

}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::computeMatrix( const vector_Type& sol,  Real const& factor)
{
    computeMatrix(this->M_stiff, sol, factor);
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::evalResidual( vector_Type &res, const vector_Type& sol, Int /*iter*/)
{
    //this->M_stiff.reset(new matrix_Type(this->M_localMap));
    computeMatrix(this->M_stiff, sol, 1.);

    this->M_Displayer->leaderPrint("    NonLin S- Updating the boundary conditions ... \t");
    Chrono chrono;
    chrono.start();
    if ( !this->M_BCh->bcUpdateDone() )
        this->M_BCh->bcUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    bcManageMatrix( *this->M_stiff, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), 1.0 );

    vector_Type rhsFull(*this->M_rhsNoBC, Unique); // ignoring non-local entries, Otherwise they are summed up lately

    bcManageVector( rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(),  this->M_data->getDataTime()->time(), 1.0 );

    *this->M_rhs = rhsFull;

    res  = *this->M_stiff*sol;
    res -= *this->M_rhs;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}



template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateJacobian( vector_Type & sol, matrixPtr_Type& jacobian  )
{
    this->M_Displayer->leaderPrint("  NonLin S-  Solid: Updating JACOBIAN... ");

    Chrono chrono;
    chrono.start();


    // copy of the linear part
    UInt totalDof = this->M_FESpace->dof().numTotalDof();

    jacobian.reset(new matrix_Type(*this->M_localMap));

    *jacobian += *this->M_linearStiff;

    Real coef;
    coef = this->M_zeta;
    *jacobian *= coef;


    UInt ig;

    ElemVec dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type dRep(sol, Repeated);

    // Number of displacement components
    UInt nc = nDimensions;

#ifdef nonlinear
    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        this->M_elmatK->zero();

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim() + this->M_offset;
                dk_loc[iloc + iComp*this->M_FESpace->fe().nbFEDof()] = dRep[ig]; // BASEINDEX + 1
            }
        }

       //non-linear terms of the stiffness matrix
	// the coefficient M_zeta is passed to the matrix in order to have the right
	// assembling of the matrix. Another solution could be export the moltiplication
	// of the nonlinear part by M_zeta outside. In this case, the moltiplication
        // must be done after the GlobalAssemble() and NOT BEFORE!!

        //  3):  \lambda * ( \tr { [\grad d^k]^T \grad \delta d }, \div v  )
        stiff_derdiv( this->M_zeta * this->M_data->getLambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  4):  \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )
        stiff_dergrad( this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // the sum of these terms is the Jacobian of the divgrad term
        // 5):  \lambda * ( (\div u_k) \grad \delta u : \grad v  )
        stiff_divgrad(  this->M_zeta * this->M_data->getLambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  \lambda * ( (\div u) \grad u_k : \grad v  )
        stiff_divgrad_2(  this->M_zeta * this->M_data->getLambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // the sum of these terms is the Jacobian of the gradgrad term
        // 6): 1/2 * \lambda * ( \grad u_k : \grad  u_k) *( \grad \delta u : \grad v  )
        stiff_gradgrad(  this->M_zeta * 0.5 * this->M_data->getLambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //\lambda * ( \grad u_k : \grad \delta u) *( \grad u_k : \grad v  )
        stiff_gradgrad_2(  this->M_zeta * this->M_data->getLambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // the sum of these terms is he jacobian of the stiff_dergrad_gradbis term
        // 7A) : \mu *  ( \grad u^k \grad \delta u : \grad v  )
        stiff_dergrad_gradbis(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  \mu *  ( \grad \delta u \grad u^k : \grad v  )
        stiff_dergrad_gradbis_2(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  the sum of these terms is he jacobian of the stiff_dergrad_gradbis_Tr term
        // 7B) :  \mu *  ( \grad u^k [\grad \delta u]^T : \grad v  )
        stiff_dergrad_gradbis_Tr(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // \mu *  ( \grad \delta u [\grad u^k]^T : \grad v  )
        stiff_dergrad_gradbis_Tr_2(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //   the sum of these terms is he jacobian of the stiff_gradgradTr_gradbis term
        // 8) :   \mu * (  \grad d^k [\grad d^k]^T \grad \delta d : \grad v  )
        stiff_gradgradTr_gradbis(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  \mu * (  \grad d^k [\grad \delta d]^T \grad d^k : \grad v  )
        stiff_gradgradTr_gradbis_2(  this->M_zeta * this->M_data->getMu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  \mu * (  \grad \delta u [\grad u^k]^T \grad u^k : \grad v  )
        stiff_gradgradTr_gradbis_3(  this->M_zeta * this->M_data->getMu() , dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // assembling
        for ( UInt ic = 0; ic < nc; ++ic )
            for ( UInt jc = 0; jc < nc; jc++ )
                assembleMatrix( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof  );
    }

#endif

    *jacobian += *this->M_mass;

    jacobian->globalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("   ... done in ", chrono.diff() );
}


//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
solveJac( vector_Type &step, const vector_Type& res, Real& linear_rel_tol)
{
    solveJacobian(step,  res, linear_rel_tol, this->M_BCh);
}


//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
solveJacobian( vector_Type&           step,
               const vector_Type&     res,
               Real&                /*linear_rel_tol*/,
               bchandler_Type&        BCh)
{
    Chrono chrono;

    updateJacobian( *this->M_disp, this->M_jacobian );

    this->M_Displayer->leaderPrint("\tS'-  Solving the linear system ... \n");

    //*this->M_f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;

    this->M_Displayer->leaderPrint("\tS'-  Applying boundary conditions      ... ");

    this->M_rhsNoBC->globalAssemble();
    this->M_rhsW->globalAssemble();

    vector_Type rhsFull (res);


//    bcManageMatrix( *this->M_jacobian, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), tgv );
    applyBoundaryConditions( *this->M_jacobian, rhsFull, BCh);


    this->M_Displayer->leaderPrintMax( "done in ", chrono.diff() );

    this->M_Displayer->leaderPrint("\tS'-  Solving system                    ... \n");
    chrono.start();

    this->M_linearSolver->setMatrix(*this->M_jacobian);

    Int numIter = this->M_linearSolver->solveSystem( rhsFull, step, this->M_jacobian );

    chrono.stop();

    *this->M_residual_d= *this->M_mass*step;


////////////////////////////////////////////////////////////////////////////////////////
//    *this->M_residual_d = *this->M_massStiff*step; // - *this->M_rhsNoBC;
}


template<typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
applyBoundaryConditions(matrix_Type&        matrix,
                        vector_Type&        rhs,
                        bchandler_Type&     BCh,
                        UInt                offset)
{

    // BC manage for the velocity
    if (offset)
        BCh->setOffset(offset);
    if ( !BCh->bcUpdateDone() )
        BCh->bcUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    // vector_Type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    //    vector_Type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts


    //In the original versione it was not commented, modified by Paolo Tricerri
//  bcManage( matrix, rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), BCh, this->M_FESpace->feBd(), 1., this->M_data->time() );


    bcManageMatrix( matrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *BCh, this->M_FESpace->feBd(), 1.,  this->M_data->time() );

    // matrix should be GlobalAssembled by  bcManage

    //The Boundary conditions should not be applied to rhs
    // rhs = rhsFull;

} // applyBoundaryCondition

template<typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
getSolidMatrix( matrixPtr_Type& matrix)
{
    //updateSystem(/*matrix*/);
    matrix.reset(new matrix_Type(*this->M_localMap));
    //*this->M_stiff *= this->M_data->dataTime()->getTimeStep() * this->M_rescaleFactor;
    *matrix  += *this->M_stiff;
    matrix->GlobalAssemble();
    matrix *= this->M_data->getDataTime()->getTimeStep() * this->M_rescaleFactor;

}

}

#endif
