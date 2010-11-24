/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
/*!
  \file NonLinearVenantKirchhofSolver.hpp

  \author Gilles Fourestey   <Gille's email>
  \author Mariarita Deluca   <Mariarita's email>
  \date   \03/2010
  \version 1.0

////////////////////////////Original/////////////////////////////////////////
  \author M.A. Fernandez
  \date 6/2003
  \version 1.0
/////////////////////////////////////////////////////////////////////////////

  \brief
  This file contains solvers for St. Venant-Kirchhof materials.

*/
#ifndef _NLVENANTKIRCHHOFSOLVER_H_
#define _NLVENANTKIRCHHOFSOLVER_H_

#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <life/lifefem/elemOper.hpp>
//#include <life/lifefem/values.hpp>
//#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifecore/chrono.hpp>

//#include <life/lifealg/dataNewton.hpp>
#include <life/lifealg/nonLinRichardson.hpp>
//#include <life/lifealg/newton.hpp>
#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>

#include <life/lifecore/displayer.hpp>


#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

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

    typedef VenantKirchhofSolver<Mesh, SolverType> super;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    typedef typename super::bchandler_type                  bchandler_type;

    typedef SolverType                             solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;


    typedef typename SolverType::prec_raw_type    prec_raw_type;
    typedef typename SolverType::prec_type        prec_type;

    typedef DataElasticStructure                   data_type;



    NonLinearVenantKirchhofSolver();

    void setup( boost::shared_ptr<data_type> data,
                const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
                boost::shared_ptr<Epetra_Comm>&     comm,
                const boost::shared_ptr<const EpetraMap>&   monolithicMap,
                UInt       offset=0
                );
void setup(
           boost::shared_ptr<data_type>        data,
           const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
           boost::shared_ptr<Epetra_Comm>&     comm
           );


    //! Update the right  hand side  for time advancing given a source term
    /*!
      \param source volumic source
      \param time present time
    */
    void updateSystem( source_type const& source, Real t );

    //! Updates the system at the end of each time step when the matrix is passed from outside
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateSystem( matrix_ptrtype& stiff );

    //! Updates the system at the end of each time step when the matrix is passed from outside
    void updateSystem( );

    //! Updates the nonlinear terms on the stiffness matrix at the end of each time step
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateNonlinearTerms( matrix_ptrtype& stiff );

    //! Updates the nonlinear terms in the matrix at the end of each Newton iteration
    /*!
      \param stiff stiffness matrix provided from outside
     */
    void updateNonlinearMatrix( matrix_ptrtype& stiff );

    //! computes the global matrix for the residual computation (calls updateNonlinearMatrix)
    /*!
      \param stiff stiffness matrix provided from outside
      \param solution vector
      \param factor scaling factor
    */
    void computeMatrix( matrix_ptrtype& stiff, const vector_type& sol, Real const& factor );

    //! computes the global matrix for the residual computation (calls updateNonlinearMatrix)
    /*!
      \param solution vector
      \param factor scaling factor
    */
    void computeMatrix( const vector_type& sol, Real const& factor );

    //! Computes the Linear part of the Jacobian and the Mass Matrix
    /*!

      \param matrix the matrix containing the mass and the linear part of the stiffness matrix
      \param factor scaling factor
     */
    void buildSystem(matrix_ptrtype matrix, Real const& factor =1.);

    //! Computes the Linear part of the Jacobian and the Mass Matrix
    void buildSystem( );

    //! Solve the non-linear system

    /*! They compute the solution solving the non linear system

      \param bch BCHandler object with the applied boundary conditions

    */
    void iterate( bchandler_type& bch );

    //! Update Jacobian at each Newton iteration
    /*!
      \param sol the current solution at the k-th iteration of Newton method
    */
    void updateJacobian( vector_type& sol, matrix_ptrtype& jac  );

    //! It is called by the NonLinear Richarson method. It calls the solveJacobian method
    /*!

      \param step the current step at the k-th iteration
      \param res the residual at the k-th iteration of the nonlinear Richardson method
      \param linear_rel_tol send for the relative tolerance to the linear solver is therefore eta. eta is determined
             by the modified Eisenstat-Walker formula
     */
     void solveJac( vector_type&       step,
                    const vector_type& res,
                    Real&            linear_rel_tol);

    //! solves the tangent problem with custom BC

    //!  It solves the tangent problem with custom BC. It solves the linear system J*step=-Res at each k-th iteration. The first three arguments are the same passed to sovleJac
    /*!
      \param BCd BCHandler object
    */
    void solveJacobian( vector_type&       step,
                        const vector_type& res,
                        Real&            linear_rel_tol,
                        bchandler_type&    BCd);

//! evaluates residual for newton interations
    void evalResidual( vector_type &res, const vector_type& sol, int iter);

    //! returns the acceleration
    vector_type& acc()         { return *M_acc; }

    //!Initializers of the initial displacement, velocity and acceleration
    void initialize   ( const Function& d0, const Function& w0, const Function& a0 );
    void initialize   ( vector_ptrtype d0,  vector_ptrtype w0 = vector_ptrtype());
    void initializeVel( const vector_type& w0);

    //!It computes the velocity vector at the n-th time step.
    void updateVel();

    //!Getter of the Offset parameter. It is taken into account when the boundary conditions are applied and the matrices are assembled.
    void getSolidMatrix( matrix_ptrtype& matrix);

private:


    //! time scheme coefficients
    // _theta and _zeta are the coefficient of the time discretization with Newmark scheme
    Real                            M_zeta;
    Real                            M_theta;

    Real                            M_pressure;


    //! linearized velocity
    vector_ptrtype                    M_acc;


    //! right  hand  side acceleration
    vector_ptrtype                    M_rhsA;

    //
    //! methods
    //

    //!It applies the boundary conditions before solving the linear system J*step=-Res. It is called in solveJacobian
    void applyBoundaryConditions(matrix_type &matrix,
                                 vector_type &rhs,
                                 bchandler_type& BCh,
                                 UInt         offset=0);

};


//
//                                         IMPLEMENTATION
//
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

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
setup(
      boost::shared_ptr<data_type>        data,
      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
      boost::shared_ptr<Epetra_Comm>&     comm,
      const boost::shared_ptr<const EpetraMap>&  monolithicMap,
      UInt                                offset
      )
{
    super::setup(data, dFESpace, comm, monolithicMap, offset);

    M_acc.reset(new vector_type(*this->M_localMap));
    M_rhsA.reset(new vector_type(*this->M_localMap));

}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
setup(
      boost::shared_ptr<data_type>        data,
      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
      boost::shared_ptr<Epetra_Comm>&     comm
      )
{
    super::setup(data, dFESpace, comm);
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
buildSystem( )
{
    super::buildSystem( );
}


template <typename Mesh, typename SolverType>
void
NonLinearVenantKirchhofSolver<Mesh, SolverType>::
buildSystem(matrix_ptrtype massStiff, Real const & factor)
{
    UInt totalDof = this->M_FESpace->dof().numTotalDof();

    this->M_Displayer->leaderPrint( "NonLin S-  Building the system             ... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    UInt nc = nDimensions;

    //inverse of dt:
    Real PTemp =  this->M_data->dataTime()->getTimeStep();
    Real dti2  = 2.0 / ( PTemp * PTemp );
//    Real dti2 = 2.0 / ( this->M_data.getTimeStep() * this->M_data.getTimeStep() );

    // Elementary computation and matrix assembling
    // Loop on elements
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        int marker    = this->M_FESpace->mesh()->volumeList( i ).marker();

        Real lambda = this->M_data->lambda(marker);
        Real mu     = this->M_data->mu    (marker);

        this->M_elmatK->zero();
        this->M_elmatM->zero();


        // stiffness
        stiff_strain (        2.0*mu, *this->M_elmatK, this->M_FESpace->fe() );

        stiff_div   ( lambda, *this->M_elmatK, this->M_FESpace->fe() );

        this->M_elmatC->mat() = this->M_elmatK->mat();

        // mass
        mass( dti2 * this->M_data->rho(), *this->M_elmatM, this->M_FESpace->fe(), 0, 0, nDimensions );

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

    this->M_linearStiff->GlobalAssemble();
    massStiff->GlobalAssemble();
    //*M_massStiff *= factor;
    this->M_mass->GlobalAssemble();

    chrono.stop();

    this->M_Displayer->leaderPrintMax( " done in ", chrono.diff() );
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
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem( matrix_ptrtype& stiff )
{
    this->M_Displayer->leaderPrint(" NonLin S-  Updating mass term on right hand side... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    //    UInt nc = nDimensions;


    Real coef;



    chrono.stop();

    // non linear part

#ifdef nonlinear


    stiff.reset(new matrix_type(*this->M_localMap));

    matrix_ptrtype tmp( new matrix_type(*this->M_localMap, 1) );
    updateNonlinearTerms(tmp);
    tmp->GlobalAssemble();
    *stiff += *tmp;
    *stiff += *this->M_linearStiff;

    stiff->GlobalAssemble();

    //_rhsWithoutBC -= _K * this->_d;


#endif

    // end
    //Computation of the right hand sides

    Real DeltaT    = this->M_data->dataTime()->getTimeStep();
    vector_type _z = *this->M_disp;

    _z            +=  DeltaT*(*this->M_vel);
    //    _z            += this->M_data->getTimeStep()*this->M_vel;

    coef= (1.0-this->M_zeta);

    *this->M_rhsNoBC  = *this->M_mass*_z;

    *this->M_rhsNoBC -= (*stiff)*coef*(*this->M_disp);


    //*this->M_rhsNoBC -= *this->M_linearStiff*(this->M_disp);

    // acceleration rhs
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * _z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);
    //    M_rhsA = (2.0 / ( M_zeta * pow(M_data->getTimeStep(),2) )) * _z + ((1.0 - M_zeta ) / ( M_zeta )) * M_acc;

    // velocity rhs
    *this->M_rhsW = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);
    //    M_rhsW = M_vel + ( 1 - M_theta  ) * M_data->getTimeStep() *  M_acc;

    /*    M_rhsW  = coef*(M_disp);
          M_rhsW += M_vel;
          coef =  1.0/(M_data->getTimeStep() * M_theta);
          M_rhsW  = coef*(M_disp);

          coef =  ( 1.0 / M_theta  - 1.0  );
          M_rhsW += coef*(M_vel);*/


    std::cout << std::endl;

    std::cout << "rhsWithoutBC norm = " << this->M_rhsNoBC->Norm2() << std::endl;
    std::cout << "rhs_w norm       = " << this->M_rhsW->Norm2() << std::endl;
    std::cout << "    w norm       = " << this->M_vel->Norm2() << std::endl;

    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());


}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem(  )
{
    updateSystem(this->M_stiff);
}



template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateNonlinearMatrix( matrix_ptrtype& stiff )
{
    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    ElemVec dk_loc( this->M_FESpace->fe().nbNode, nDimensions );
    vector_type disp(*this->M_disp);
    vector_type dRep(disp, Repeated);

    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbNode ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim();
                dk_loc[ iloc + iComp*this->M_FESpace->fe().nbNode ] = dRep[ig]; // BASEINDEX + 1
            }
        }


        this->M_elmatK->zero();
        // stiffness for non-linear terms


        // 3) 1/2 * \lambda * ( \tr { [\grad d^k]^T \grad d }, \div v  )

        stiff_derdiv( this->M_zeta * 0.5 * this->M_data->lambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_derdiv( 0.5 * this->M_data.lambda(), dk_loc, this->M_elmatK,  this->M_FESpace->fe() );


        // i risultati degli integrali sono memorizzati
        // di volta in volta in elmatK, fino a costruire la matrice K
        // della forma bilineare

        //4)  \mu * ( [\grad d^k]^T \grad d : \grad v  )
        //stiff_dergradbis( this->M_data.mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5

        stiff_dergradbis( this->M_zeta* this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_dergradbis( this->M_data.mu() , dk_loc, this->M_elmatK,  this->M_FESpace->fe() );

        // *******************aggiunti Rita emory june 2008

        //  5):\lambda * (div u_k) \grad d : \grad v
        //  stiff_divgrad( this-> _lambda * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5

        stiff_divgrad( this->M_zeta * this->M_data->lambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_divgrad( this->M_data->lambda(), dk_loc, this->M_elmatK,  this->M_FESpace->fe() );


        // 6)  \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
        // stiff_gradgrad( this-> _lambda * 0.25, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5

        stiff_gradgrad( this->M_zeta *  0.5 * this->M_data->lambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_gradgrad( 0.5 * this->M_data->lambda() , dk_loc, this->M_elmatK,  this->M_FESpace->fe() );


        // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
        // stiff_dergrad_gradbis( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() );//zeta=0.5

        stiff_dergrad_gradbis( this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_dergrad_gradbis( this->M_data->mu(), dk_loc, this->M_elmatK,  this->M_FESpace->fe() );

        // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
        // stiff_dergrad_gradbis_Tr( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() );//zeta=0.5

        stiff_dergrad_gradbis_Tr( this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_dergrad_gradbis_Tr( this->M_data->mu(), dk_loc, this->M_elmatK,  this->M_FESpace->fe() );


        // 8) // \mu * (  \grad d^k [\grad d^k]^T \grad d : \grad v  )

        stiff_gradgradTr_gradbis( this->M_zeta * this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
//        stiff_gradgradTr_gradbis(  this->M_data->mu() , dk_loc, this->M_elmatK,  this->M_FESpace->fe() );


        //std::cout << "---------------------Esce dalla parte nonlineare del TIMEADVANCE----------------------------"<< std::endl;

        UInt totalDof   = this->M_FESpace->dof().numTotalDof();

        //--------------------------------------------------------------------------------
        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            // stiff la matrice non lineare della forma bilineare
            for ( UInt jc = 0; jc < nDimensions; jc++ )
                assembleMatrix( *stiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);

        }
    }
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateNonlinearTerms( matrix_ptrtype& stiff )
{

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    ElemVec dk_loc( this->M_FESpace->fe().nbNode, nDimensions );
    vector_type disp(*this->M_disp);
    //    disp *= this->M_data->dataTime()->getTimeStep();
    vector_type dRep(disp, Repeated);

    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); i++ )
    {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );


        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbNode ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim() + this->M_offset;
                dk_loc[ iloc + iComp*this->M_FESpace->fe().nbNode ] = dRep[ig]; // BASEINDEX + 1
            }
        }

        this->M_elmatK->zero();

        // stiffness for non-linear terms

        // 3) 1/2 * \lambda  ( \tr { [\grad d^k]^T \grad d }, \div v  )
        // stiff_derdiv( this->M_data->lambda() * 0.25, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5
        stiff_derdiv( 0.5 * this->M_data->lambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
        // i risultati degli integrali sono memorizzati
        // di volta in volta in elmatK, fino a costruire la matrice K
        // della forma bilineare

        //4)  \mu *( [\grad d^k]^T \grad d : \grad v  )
        //stiff_dergradbis( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5
        stiff_dergradbis( this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // *******************aggiunti Rita emory june 2008

        //  5): \lambda * (div u_k) \grad d : \grad v
        // stiff_divgrad( this-> _lambda * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5
        stiff_divgrad( this->M_data->lambda(), dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 6) 1/2  * \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
        // stiff_gradgrad( this-> _lambda * 0.25, dk_loc, _elmatK,  this->M_FESpace->fe() ); //zeta=0.5
        stiff_gradgrad( 0.5 * this->M_data->lambda() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
        // stiff_dergrad_gradbis( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() );//zeta=0.5
        stiff_dergrad_gradbis( this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
        // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
        // stiff_dergrad_gradbis_Tr( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() );//zeta=0.5
        stiff_dergrad_gradbis_Tr( this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 8) // \mu *  (  \grad d^k [\grad d^k]^T \grad d : \grad v  )
        // stiff_gradgradTr_gradbis( this->M_data->mu() * 0.5, dk_loc, _elmatK,  this->M_FESpace->fe() );   // commentiamo la costruzione di K per il caso statico
        stiff_gradgradTr_gradbis( this->M_data->mu() , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //std::cout << "---------------------Esce dalla parte nonlineare del TIMEADVANCE----------------------------"<< std::endl;

        //--------------------------------------------------------------------------------
        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            // K è la matrice non lineare della forma bilineare
            for ( UInt jc = 0; jc < nDimensions; jc++ )
                //assemb_mat( _K, this->M_elmatK,  this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc );
                assembleMatrix( *stiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);
        }

    }
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateSystem(  source_type const& source, Real t )
{

    this->M_Displayer->leaderPrint(" NonLin S-  Updating mass term on right hand side... ");

    Chrono chrono;
    chrono.start();

    // Number of displacement components
    //    UInt nc = nDimensions;


    Real coef;

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    *this->M_rhsNoBC *= 0.0;


    chrono.stop();

    // non linear part

#ifdef nonlinear
    ElemVec dk_loc( this->M_FESpace->fe().nbNode, nDimensions );

    vector_type disp(*this->M_disp);
    //    disp *= this->M_data->dataTime()->getTimeStep();
    vector_type dRep(disp, Repeated);


    this->M_stiff.reset(new matrix_type(*this->M_localMap));

    matrix_ptrtype tmp(new matrix_type(*this->M_localMap, 1));
    updateNonlinearTerms( tmp );
    //updateNonlinearMatrix( tmp );
    tmp->GlobalAssemble();
    *this->M_stiff += *tmp;

    *this->M_stiff += *this->M_linearStiff;

    this->M_stiff->GlobalAssemble();

    //_rhsWithoutBC -= _K * this->_d;


#endif

    ElemVec M_elvec(this->M_FESpace->fe().nbNode, nDimensions);
    UInt nc = nDimensions;

// loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
      {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        M_elvec.zero();

        int nothing;
        for ( UInt ic = 0; ic < nc; ++ic )
        {
          compute_vec( source, M_elvec, this->M_FESpace->fe(),  this->M_data->dataTime()->getTime(), ic ); // compute local vector
          assembleVector( *this->M_rhsNoBC, M_elvec, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, ic*this->M_FESpace->dim() ); // assemble local vector into global one
        }
      }

    this->M_rhsNoBC->GlobalAssemble();
    // end
    //Computation of the right hand sides

    Real DeltaT    = this->M_data->dataTime()->getTimeStep();
    vector_type _z = *this->M_disp;

    _z            +=  DeltaT*(*this->M_vel);
    //    _z            += this->M_data->getTimeStep()*this->M_vel;

    coef= (1.0-M_zeta);

    *this->M_rhsNoBC += *this->M_mass*_z;

    *this->M_rhsNoBC -= (*this->M_stiff)*coef*(*this->M_disp);

    //*this->M_rhsNoBC -= *this->M_linearStiff*(this->M_disp);

    // acceleration rhs
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * _z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);
    //    M_rhsA = (2.0 / ( M_zeta * pow(M_data->getTimeStep(),2) )) * _z + ((1.0 - M_zeta ) / ( M_zeta )) * M_acc;

    // velocity rhs
    *this->M_rhsW = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);
    //    M_rhsW = M_vel + ( 1 - M_theta  ) * M_data->getTimeStep() *  M_acc;

    /*    M_rhsW  = coef*(M_disp);
          M_rhsW += M_vel;
          coef =  1.0/(M_data->getTimeStep() * M_theta);
          M_rhsW  = coef*(M_disp);

          coef =  ( 1.0 / M_theta  - 1.0  );
          M_rhsW += coef*(M_vel);*/


    std::cout << std::endl;

    std::cout << "rhsWithoutBC norm = " << this->M_rhsNoBC->Norm2() << std::endl;
    std::cout << "rhs_w norm       = " << this->M_rhsW->Norm2() << std::endl;
    std::cout << "    w norm       = " << this->M_vel->Norm2() << std::endl;

    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());


}



template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
iterate( bchandler_type& bch )
{
    Chrono chrono;

    // matrix and vector assembling communication
    this->M_Displayer->leaderPrint("  NonLin S-  Solving the system ... \n");

    this->M_BCh = bch;

    Real abstol  = 1.e-5;
    Real reltol  = 1.e-5;
    UInt   maxiter = 200;
    Real etamax  = 0;
    int linesearch = 0;


    Real time = this->M_data->getTime();

    int status = 0;

    status = nonLinRichardson( *this->M_disp, *this, abstol, reltol, maxiter, etamax, linesearch, this->M_out_res, this->M_data->dataTime()->getTime() );



    if ( status == 1 )
    {
        std::ostringstream __ex;
        __ex << "VenantKirchhofSolver::iterate() Inners newton iterations failed to converge\n";
        throw std::logic_error( __ex.str() );
    }
    else // if status == 0 vuol dire che Newton converge
    {
    	std::cout << std::endl;

        std::cout <<" Number of inner iterations       : " << maxiter <<  std::endl;

        std::cout <<" We are at the time step          : "  << this->M_data->dataTime()->getTime() << std::endl;
        this->M_out_iter << time << " " << maxiter << std::endl;
    }

    updateVel();
    //   this->M_acc = (2.0 /( this->M_zeta * pow(this->M_data->getTimeStep(),2) ))  * this->M_disp  - this->M_rhsA;
    //   this->M_vel = this->M_rhsW + this->M_theta * this->M_data->getTimeStep() * this->M_acc ;

    std::cout << "iterate: d norm       = " << this->M_disp->Norm2() << std::endl;
    std::cout << "iterate: w norm       = " << this->M_vel->Norm2() << std::endl;
    std::cout << "iterate: a norm       = " << M_acc->Norm2() << std::endl;
    ///////////////////////////////////////////////////////////////////////////////
    //    M_vel  = ( 1.0 / M_data->getTimeStep() )*M_theta*(M_disp);
    //    M_vel -= M_rhsW;


    *this->M_residual_d = *this->M_mass*(*this->M_disp);
    *this->M_residual_d -= *this->M_rhsNoBC;

    //    *M_residual_d  = *M_massStiff*(M_disp);
    //    *M_residual_d -= *M_rhsNoBC;

} // iterate()


template <typename Mesh, typename SolverType> // for monolithic
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
updateVel()
{
    Real DeltaT = this->M_data->dataTime()->getTimeStep();
    *M_acc = (2.0 /( M_zeta * pow(DeltaT,2) ))  * (*this->M_disp)  - *M_rhsA;
    *this->M_vel = *this->M_rhsW + M_theta * DeltaT * (*M_acc) ;
}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::computeMatrix( matrix_ptrtype& stiff, const vector_type& sol,  Real const& factor)
{
    this->M_Displayer->leaderPrint( "    NonLin S- Computing residual ... \t\t\t");
    /////////////////PaoloTricerri/////////////////////
    std::cout << "The norm of the solution is:" << sol.Norm2() << std::endl;
    ///////////////////////////////////////////////////

    Chrono chrono;
    chrono.start();

    // Matrices initialization

    stiff.reset(new matrix_type(*this->M_localMap));
    *stiff += *this->M_linearStiff;

    Real coef;
    coef=M_zeta;
    *stiff *= coef;


    ElemVec dk_loc( this->M_FESpace->fe().nbNode, nDimensions );

    vector_type dRep(sol, Repeated);


#ifdef nonlinear

    matrix_ptrtype tmp(new matrix_type(*this->M_localMap, 1));
    updateNonlinearMatrix( tmp );
    tmp->GlobalAssemble();
    *tmp *= factor;
    *stiff += *tmp;

#endif

    *stiff += *this->M_mass;
    stiff->GlobalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );

}


template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::computeMatrix( const vector_type& sol,  Real const& factor)
{
    computeMatrix(this->M_stiff, sol, factor);
}

template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::evalResidual( vector_type &res, const vector_type& sol, int /*iter*/)
{
    //this->M_stiff.reset(new matrix_type(this->M_localMap));
    computeMatrix(this->M_stiff, sol, 1.);

    //    int nothing;
    //    std::cin >> nothing;
    this->M_Displayer->leaderPrint("    NonLin S- Updating the boundary conditions ... \t");
    Chrono chrono;
    chrono.start();
    if ( !this->M_BCh->bdUpdateDone() )
        this->M_BCh->bdUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    bcManageMatrix( *this->M_stiff, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), 1.0 );

    vector_type rhsFull(*this->M_rhsNoBC, Unique); // ignoring non-local entries, Otherwise they are summed up lately

    bcManageVector( rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(),  this->M_data->dataTime()->getTime(), 1.0 );
    *this->M_rhs = rhsFull;

    res  = *this->M_stiff*sol;
    res -= *this->M_rhs;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}



template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::updateJacobian( vector_type & sol, matrix_ptrtype& jacobian  )
{
    this->M_Displayer->leaderPrint("  NonLin S-  Solid: Updating JACOBIAN... ");

    Chrono chrono;
    chrono.start();


    // copy of the linear part
    UInt totalDof = this->M_FESpace->dof().numTotalDof();

    jacobian.reset(new matrix_type(*this->M_localMap));

    *jacobian += *this->M_linearStiff;

    Real coef;
    coef = this->M_zeta;
    *jacobian *= coef;


    UInt ig;

    ElemVec dk_loc( this->M_FESpace->fe().nbNode, nDimensions );

    vector_type dRep(sol, Repeated);

    // Number of displacement components
    UInt nc = nDimensions;

#ifdef nonlinear
        // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        this->M_elmatK->zero();

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbNode ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobal( eleID, iloc + 1 ) + iComp*this->dim() + this->M_offset;
                dk_loc[iloc + iComp*this->M_FESpace->fe().nbNode] = dRep[ig]; // BASEINDEX + 1
            }
        }

        //*****************JACOBIAN*******************************************************

        // Costruzione della matrice _K
        // con il paramentro generico dello schema di Newmark zeta
        // stiffness for non-linear terms

        //  3):  \lambda * ( \tr { [\grad d^k]^T \grad \delta d }, \div v  )

        //stiff_derdiv( this->M_data->lambda(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_derdiv( this->M_zeta * this->M_data->lambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //  4):  \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )

        //stiff_dergrad( this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_dergrad( this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // -------------------------- aggiunti Rita emory 2008---------------------------------------------------------------

        // la somma di questi due termini rappresenta lo jacobiano del termine divgrad
        // 5):  \lambda * ( (\div u_k) \grad \delta u : \grad v  )

        //stiff_divgrad(this->M_data->lambda(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_divgrad(  this->M_zeta * this->M_data->lambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        //  \lambda * ( (\div u) \grad u_k : \grad v  )

        //stiff_divgrad_2(this->M_data->lambda(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_divgrad_2(  this->M_zeta * this->M_data->lambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        // la somma di questi due termini rappresenta lo jacobiano del termine gradgrad
        // 6): 1/2 * \lambda * ( \grad u_k : \grad  u_k) *( \grad \delta u : \grad v  )

        //stiff_gradgrad(0.5 * this->M_data->lambda(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_gradgrad(  this->M_zeta * 0.5 * this->M_data->lambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        //\lambda * ( \grad u_k : \grad \delta u) *( \grad u_k : \grad v  )

        //stiff_gradgrad_2(this->M_data->lambda(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_gradgrad_2(  this->M_zeta * this->M_data->lambda(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //   la somma di questi due termini rappresenta lo jacobiano del termine stiff_dergrad_gradbis
        // 7A) : \mu *  ( \grad u^k \grad \delta u : \grad v  )

        //stiff_dergrad_gradbis( this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_dergrad_gradbis(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        //  \mu *  ( \grad \delta u \grad u^k : \grad v  )

        //stiff_dergrad_gradbis_2(this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_dergrad_gradbis_2(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //   la somma di questi due termini rappresenta lo jacobiano del termine stiff_dergrad_gradbis_Tr
        // 7B) :  \mu *  ( \grad u^k [\grad \delta u]^T : \grad v  )

        //stiff_dergrad_gradbis_Tr(this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_dergrad_gradbis_Tr(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        // \mu *  ( \grad \delta u [\grad u^k]^T : \grad v  )

        //stiff_dergrad_gradbis_Tr_2(this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_dergrad_gradbis_Tr_2(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //   la somma di questi tre termini rappresenta lo jacobiano del termine stiff_gradgradTr_gradbis
        // 8) :   \mu * (  \grad d^k [\grad d^k]^T \grad \delta d : \grad v  )

        //stiff_gradgradTr_gradbis(this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_gradgradTr_gradbis(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        //  \mu * (  \grad d^k [\grad \delta d]^T \grad d^k : \grad v  )

        //siff_gradgradTr_gradbis_2(this->M_data->mu(), dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_gradgradTr_gradbis_2(  this->M_zeta * this->M_data->mu(), dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
        //  \mu * (  \grad \delta u [\grad u^k]^T \grad u^k : \grad v  )
        //stff_gradgradTr_gradbis_3(this->M_data->mu() , dk_loc, this->M_elmatK, this->M_FESpace->fe() );
        stiff_gradgradTr_gradbis_3(  this->M_zeta * this->M_data->mu() , dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

        //----------------------------------------------------------------------------------------------------------

        // assembling
        for ( UInt ic = 0; ic < nc; ++ic )
            for ( UInt jc = 0; jc < nc; jc++ )
                assembleMatrix( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof  );
    }

#endif

    *jacobian += *this->M_mass;

    jacobian->GlobalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("   ... done in ", chrono.diff() );
}


//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
solveJac( vector_type &step, const vector_type& res, Real& linear_rel_tol)
{
    solveJacobian(step,  res, linear_rel_tol, this->M_BCh);
}


//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
solveJacobian( vector_type&           step,
               const vector_type&     res,
               Real&                /*linear_rel_tol*/,
               bchandler_type&        BCh)
{
    Chrono chrono;

    updateJacobian( *this->M_disp, this->M_jacobian );

    this->M_Displayer->leaderPrint("\tS'-  Solving the linear system ... \n");

    //*this->M_f = res;

    // for BC treatment (done at each time-step)
    Real tgv = 1.0;

    this->M_Displayer->leaderPrint("\tS'-  Applying boundary conditions      ... ");

    this->M_rhsNoBC->GlobalAssemble();
    this->M_rhsW->GlobalAssemble();

    vector_type rhsFull (res);


//    bcManageMatrix( *this->M_jacobian, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), tgv );
    applyBoundaryConditions( *this->M_jacobian, rhsFull, BCh);


    this->M_Displayer->leaderPrintMax( "done in ", chrono.diff() );

    this->M_Displayer->leaderPrint("\tS'-  Solving system                    ... \n");
    chrono.start();

    this->M_linearSolver->setMatrix(*this->M_jacobian);

    int numIter = this->M_linearSolver->solveSystem( rhsFull, step, this->M_jacobian );



    //step *= -1.;
    //    this->M_linearSolver.solve( step , _f);
    chrono.stop();

    *this->M_residual_d= *this->M_mass*step;


////////////////////////////////////////////////////////////////////////////////////////
 //    *this->M_residual_d = *this->M_massStiff*step; // - *this->M_rhsNoBC;
}


template<typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
applyBoundaryConditions(matrix_type&        matrix,
                        vector_type&        rhs,
                        bchandler_type&     BCh,
                        UInt                offset)
{

    // BC manage for the velocity
    if(offset)
        BCh->setOffset(offset);
    if ( !BCh->bdUpdateDone() )
        BCh->bdUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    // vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    //    vector_type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts


    //In the original versione it was not commented, modified by Paolo Tricerri
//  bcManage( matrix, rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), BCh, this->M_FESpace->feBd(), 1., this->M_data->getTime() );


    bcManageMatrix( matrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *BCh, this->M_FESpace->feBd(), 1.,  this->M_data->getTime() );

    // matrix should be GlobalAssembled by  bcManage

   //The Boundary conditions should not be applied to rhs
   // rhs = rhsFull;

} // applyBoundaryCondition

template<typename Mesh, typename SolverType>
void NonLinearVenantKirchhofSolver<Mesh, SolverType>::
getSolidMatrix( matrix_ptrtype& matrix)
{
    //updateSystem(/*matrix*/);
    matrix.reset(new matrix_type(*this->M_localMap));
    //*this->M_stiff *= this->M_data->dataTime()->getTimeStep() * this->M_rescaleFactor;
    *matrix  += *this->M_stiff;
    matrix->GlobalAssemble();
    matrix *= this->M_data->dataTime()->getTimeStep() * this->M_rescaleFactor;
}

}

#endif
