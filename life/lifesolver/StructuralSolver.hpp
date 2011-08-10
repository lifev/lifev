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
 *  @brief This file contains solvers for different materials (St. Venant-Kirchhoff materials right now )
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _STRUCTURALSOLVER_H_
#define _STRUCTURALSOLVER_H_ 1

//#include<boost/scoped_ptr.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifecore/LifeChrono.hpp>
#include <life/lifecore/Displayer.hpp>

#include <life/lifearray/MatrixElemental.hpp>
#include <life/lifearray/VectorElemental.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>

#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/BCManage.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifealg/SolverAztecOO.hpp>
#include <life/lifealg/NonLinearRichardson.hpp>

#include <life/lifesolver/VenantKirchhoffElasticData.hpp>
#include <life/lifesolver/StructuralMaterial.hpp>


namespace LifeV
{

/*!
  \class StructuralSolver
  \brief
  This class solves the linear elastodynamics equations for different kinds of materials
  (St. Venant-Kirchoff materials right now)

*/
template <typename Mesh,
     	  typename SolverType = LifeV::SolverAztecOO >

class StructuralSolver
{
public:

  //!@name Type definitions
  //@{
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_Type;

    typedef StructuralMaterial<Mesh>               material_Type;
    typedef boost::shared_ptr<material_Type>       materialPtr_Type;

    typedef BCHandler                              bchandlerRaw_Type;
    typedef boost::shared_ptr<bchandlerRaw_Type>   bchandler_Type;

    typedef SolverType                             solver_Type;

    typedef typename solver_Type::matrix_type      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;
    typedef typename solver_Type::vector_type      vector_Type;
    typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

    typedef typename SolverType::prec_raw_type     precRaw_Type;
    typedef typename SolverType::prec_type         prec_Type;

    typedef VenantKirchhoffElasticData             data_Type;

  //@}


  //! @name Constructor & Destructor
  //@{

  StructuralSolver();

  virtual ~StructuralSolver() {};

   //@}

  //!@name Methods
  //@{

  //! Setup the created object of the class Venantkirchhof
  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param BCh boundary conditions for the displacement
    \param comm the Epetra Comunicator
  */
  void setup( boost::shared_ptr<data_Type> data,
	      const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   FESpace,
	      bchandler_Type&       BCh,
	      boost::shared_ptr<Epetra_Comm>&     comm
	      );

  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param comm the Epetra Comunicator
  */
  void setup( boost::shared_ptr<data_Type> data,
	      const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   FESpace,
	      boost::shared_ptr<Epetra_Comm>&     comm
	      );

  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param comm the comunicator parameter
    \param monolithicMap the MapEpetra
    \param offset the offset parameter
  */
  void setup( boost::shared_ptr<data_Type> data,
		      const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   dFESpace,
		      boost::shared_ptr<Epetra_Comm>&     comm,
		      const boost::shared_ptr<const MapEpetra>&       monolithicMap,
		      UInt       offset=0
		      //boost::shared_ptr<FESpace<Mesh, MapEpetra> >   uFESpace=0
		      );


  //! Updates the system at the end of each time step
  void updateSystem( void );

  //! Updates the system at the end of each time step when the matrix is passed from outside
  /*!
    \param stiff stiffness matrix provided from outside
  */
  void updateSystem(matrixPtr_Type& stiff);

  //! Updates the system at the end of each time step given a source term
  /*!
    \param source volumic source
    \param time present time
  */
  void updateSystem( source_Type const& source );


  //! Updates the system at the end of each time step given a source term
  /*!
    \param source volumic source
    \param time present time
  */
  void updateSourceTerm( source_Type const& source );

  //! Updates the rhs at the start of each time step
  /*!
  \param rhs: solid  right hand side
  !*/
  void updateRightHandSide(const vector_Type& rightHandSide) { *M_rhsNoBC = rightHandSide;};

  //! Comuptes the right hand side in the updateSystem methods
  void computeRightHandSide( void );

    //! Compute the mass matrix and it calls the method to build the linear part of the stiffness matrix of the material class
    void buildSystem( const Real& coefficient );

    //void buildSystem(matrix_Type & bigMatrixStokes, const Real& timeAdvanceCoefficient, const Real& factor); // used for monolithic

  //! Compute the mass matrix and the linear part of the stiffness matrix
  /*!
    \param matrix the matrix containing the mass matrix and the linear part of he stiffness matrix
    \param rescale factor for FSI problems
  */
  void computeMassMatrix( const Real& factor=1.);

  //! Solve the non-linear system
  /*!
    \param bch BCHander object containing the boundary conditions
  */
  void iterate( bchandler_Type& bch );

  //! Solve the linearized problem. Used in FSI segregated in ExactJacobian
  /*!
    \param bch BCHander object containing the boundary conditions
  */
  void iterateLin( bchandler_Type& bch );


  //! Output
  /*!
    \param c output file
  */
  void showMe( std::ostream& c = std::cout ) const;

  //! Update the Jacobian Matrix at each iteration of the nonLinearRichardson method
  /*!
    \param sol the current solution at each iteration of the nonLinearRichardson method
    \param jac the Jacobian matrix that must be updated
  */
  void updateJacobian( vector_Type& solution, matrixPtr_Type& jacobian );

  //! Solves the tangent problem for newton iterations
  /*!
    \param step the vector containing the solution of the sistem J*step=-Res
    \param res the vector conteining the residual
    \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
           eta is determined by the modified Eisenstat-Walker formula
  */
  void solveJac( vector_Type&       step,
		 const vector_Type& residual,
		 Real&            linear_rel_tol) ;
  //    void solveJac( const Vector& res, Real& linear_rel_tol, Vector &step);

  //! Solves the tangent problem with custom BC
  /*!
    \param step the vector containing the solution of the sistem J*step=-Res
    \param res the vector conteining the residual
    \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
           eta is determined by the modified Eisenstat-Walker formula
    \param BCd BCHandler object containing the boundary condition
  */
  void solveJacobian( vector_Type&       step,
		      const vector_Type& residual,
		      Real&            linear_rel_tol,
		      bchandler_Type&    BCd ) ;


  //! Evaluates residual for newton interations
  /*!
    \param res residal vector that is update every time the method is called
    \param sol solution vector from which the residual is computed
    \param iter iteration of the nonLinearRichardson method
  */
  void evalResidual( vector_Type &residual, const vector_Type& solution, Int iter);

  //! Evaluates residual of the displacement for FSI problems
  /*!
    \param sol, the current displacement of he sturcture
  */
  void evalResidualDisplacement( const vector_Type& solution );

  //! Evaluates residual of the displacement in the Linearized problem of ExactJcobian. FSI problems
  /*!
    \param sol, the current displacement of he sturcture
  */
  void evalResidualDisplacementLin( const vector_Type& solution );

  void evalConstraintTensor();

  //! Sets the initial displacement, velocity, acceleration
  /*!
    \param d0 space function describing the initial displacement
    \param w0 space function describing the initial velocity
    \param a0 space function describing the initial acceleration
  */
  void initialize( const Function& d0, const Function& w0, const Function& a0 );

  //! Sets the initial velocity
  /*!
    \param w0 space function describing the initial velocity
  */
  void initializeVel( const vector_Type& w0);

  //! Sets the initial displacement, velocity, acceleration
  /*!
    \param d0 space function describing the initial displacement
    \param w0 empty vector
    \param a0 empty vector
  */
  void initialize( vectorPtr_Type d0,  vectorPtr_Type w0 = vectorPtr_Type(),  vectorPtr_Type a0 = vectorPtr_Type() );

  //! Computes the velocity and acceleration vector at the n-th time step
    //void updateVelAndAcceleration();

  //! Reduce the complete solution to the solution on the pocessor with rank 0
  /*!
    \param disp displacement solution
    \param vel velocity solution
  */
  void reduceSolution( Vector& displacement, Vector& velocity );

  //! Multiply the mass matrix and the linear stiffness matrix by the rescaleFactor
    //  void rescaleMatrices(); // used for monolithic

  /**
     in the linear case the solid matrix is constant, thus it does not need to be recomputed.
  */

  //! Update (in the case of nonlinear material) the solid matrix
  /*!
    \param stiff stiffness matrix
    \param sol the current solution
    \param factor the rescaleFactor
  */
  void computeMatrix( matrixPtr_Type& stiff, const vector_Type& sol, Real const& factor );

  //void updateMatrix(matrix_Type & bigMatrixStokes);// used for monolithic
  //void updateCoupling(matrix_Type couplingMatrix);// used for monolithic

  //@}

  //! @name Set Methods
  //@{

  //!Setters
  //! Set the BCHandler object
  void setBC(bchandler_Type& BCd)   {M_BCh = BCd;}

  //! Set the source object
  void setSourceTerm( source_Type const& s ) { M_source = s; }

  //! Set the preconditioner
  void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }

  //! Set the displacement
  virtual void setDisp(const vector_Type& disp) {*M_disp = disp;} // used for monolithic

  //! Set the recur parameter
  void setRecur(UInt recur) {M_recur = recur;}

  //! Set the data fields with the Getpot data file
  void setDataFromGetPot( const GetPot& dataFile );

  //@}


  //! @name Get Methods
  //@{

  //! Getters
  //! Get the Epetramap
  MapEpetra   const& map()       const { return *M_localMap; }

  //! Get the Displayer object
    Displayer   const& displayer() const { return *M_Displayer; }

    boost::shared_ptr<const Displayer>   const& displayerPtr() const { return M_Displayer; }

  //! Get the matrix containing the mass mtrix and the linear part of the stiffness matrix
  //matrixPtr_Type const MassStiff() const {return M_massStiff; }

  //! Get the mass matrix
  matrixPtr_Type const Mass() const {return M_mass; }

  //! Get the FESpace object
  FESpace<Mesh, MapEpetra>& dFESpace() {return M_FESpace;}

  //! Get the bCHandler object
  bchandler_Type const & BChandler() const {return M_BCh;}

  //! Get the residual
  vector_Type& residual()             {return *M_residual_d;}

  //! Get the source term
  source_Type const& sourceTerm() const { return M_source; }

  //! Get the displacement
  vector_Type& displacement()        { return *M_disp; }

  //! Get the velocity
  //vector_Type& velocity()         { return *M_vel; }

  //! Get the velocity
  //vector_Type& acceleration()         { return *M_acc; }

  //! Get the right hand sde without BC
  vectorPtr_Type& rhsWithoutBC() { return M_rhsNoBC; }

  //! Get the comunicator object
  boost::shared_ptr<Epetra_Comm> const& getComunicator() const {return M_Displayer->comm();}

  //! Get the rescaleFactor
  Real rescaleFactor() {return M_rescaleFactor;}

  /*! Get the offset parameter. It is taken into account when the boundary conditions
    are applied and the matrices are assembled.
  */
  const UInt& offset() const { return M_offset; }

  /*! Get the offset parameter. It is taken into account when the boundary conditions
    are applied and the matrices are assembled.
  */
    const materialPtr_Type& material() const { return M_material; }

  /**
     Do nothing in the linear case: the matrix remains constant. Otherwise substitute the matrix with an updated one
  */
  //! Get the Solid Matrix
  void solidMatrix( matrixPtr_Type& /*matrix*/ )
  {
  }

  // Physic constant
  //! Get the thickness
  Real thickness() const { return M_data->thickness(); }

  //! Get the Young modulus
  Real young( UInt material )            const { return M_data->young( material ); }

  //! Get the Poisson coefficient
  Real poisson( UInt material )          const { return M_data->poisson( material ); }

  //! Get the density
  Real rho()       const { return M_data->rho(); }

  //@}

protected:

  //! Apply boundary condition
  /*!
    \param matrix the matrix of the system
    \param rhs the right hand side of the system
    \param BCh BCHandler object
    \param offset the offset parameter
  */
  void applyBoundaryConditions(matrix_Type &matrix,
			       vector_Type &rhs,
			       bchandler_Type& BCh,
			       UInt         offset=0);

  void applyBoundaryConditionsLin(matrix_Type &matrix,
			       vector_Type &rhs,
			       bchandler_Type& BCh,
			       UInt         offset=0);


  UInt dim() const { return M_FESpace->dim(); }

  //!Protected Members

  boost::shared_ptr<data_Type>         M_data;

  boost::shared_ptr<FESpace<Mesh, MapEpetra> >      M_FESpace;

  boost::shared_ptr<const Displayer>         M_Displayer;

  Int                                  M_me;

  //! data for solving tangent problem with aztec
  boost::shared_ptr<solver_Type>       M_linearSolver;

  //! Elementary matrices and vectors
  boost::shared_ptr<MatrixElemental>           M_elmatM; // mass

  //! linearized velocity
  vectorPtr_Type                       M_disp;
  //vectorPtr_Type                       M_vel;
  //vectorPtr_Type                       M_acc;

  //! right  hand  side displacement
  vectorPtr_Type                       M_rhs;

  //! right  hand  side velocity
  //  vectorPtr_Type                       M_rhsW;

  //! right  hand  side velocity
  //vectorPtr_Type                       M_rhsA;

  //! right  hand  side
  vectorPtr_Type                       M_rhsNoBC;

  //! right  hand  side
  //boost::shared_ptr<vector_Type>       M_f;

  //! residual
  boost::shared_ptr<vector_Type>       M_residual_d;

  //! Components of the Constraint Tensor
  vectorPtr_Type                       M_sxx;
  vectorPtr_Type                       M_syy;
  vectorPtr_Type                       M_szz;

  //! files for lists of iterations and residuals per timestep
  std::ofstream                        M_out_iter;
  std::ofstream                        M_out_res;

  //! BCHandler object
  bchandler_Type                       M_BCh;

  //! Map Epetra
  boost::shared_ptr<const MapEpetra>   M_localMap;

  //! Matrix M: mass
  matrixPtr_Type                       M_mass;

  //! Matrix mass: M * xi /(dt*dt)
  matrixPtr_Type                       M_massTimeAdvanceCoefficient;

  //! Matrix Temp: Temporary matrix to compute residuals, store jacobian
  matrixPtr_Type                       M_tempMatrix;
  //! Jacobian Matrix: Matrix to store the jacobian of the newton method
  matrixPtr_Type                       M_jacobian;

  //! level of recursion for Aztec (has a sens with FSI coupling)
  UInt                                 M_recur;

  source_Type                          M_source;

  UInt                                 M_offset;
  Real                                 M_rescaleFactor;
//  Real                                 M_zeta;
//  Real                                 M_theta;

  //! Material class
  materialPtr_Type                     M_material;

};

//====================================
// Constructor
//=====================================

template <typename Mesh, typename SolverType>
StructuralSolver<Mesh, SolverType>::StructuralSolver( ):
  M_data                       ( ),
  M_FESpace                    ( ),
  M_Displayer                  ( ),
  M_me                         ( 0 ),
  M_linearSolver               ( ),
  M_elmatM                     ( ),
  M_disp                       ( ),
 // M_vel                        ( ),
 // M_acc                        ( ),
 // M_rhs                        ( /*new vector_Type(M_localMap)*/),//useful
 // M_rhsW                       ( ),
 // M_rhsA                       ( ),
  M_rhsNoBC                    ( ),
  M_residual_d                 ( ),
  M_sxx                        (/*M_localMap*/),//useless
  M_syy                        (/*M_localMap*/),//useless
  M_szz                        (/*M_localMap*/),//useless
  M_out_iter                   ( "out_iter_solid" ),
  M_out_res                    ( "out_res_solid" ),
  M_BCh                        (),
  M_localMap                   ( ),
  M_mass                       ( ),
  M_massTimeAdvanceCoefficient ( ),
  M_tempMatrix                 ( ),
  M_jacobian                   ( ),
  M_recur                      ( ),
  M_source                     ( ),
  M_offset                     ( 0 ),
  M_rescaleFactor              ( 1. ),
  //M_zeta                       ( 0.75 ),
 // M_theta                      ( 0.7 ),
  M_material                   ( )
{
  std::cout << "I am in the constructor for the solver" << std::endl;
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>          data,
					      const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
					      bchandler_Type&                BCh,
					      boost::shared_ptr<Epetra_Comm>&              comm)
{
  setup(data, dFESpace, comm);
  M_BCh = BCh;
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>        data,
					  const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
					  boost::shared_ptr<Epetra_Comm>&     comm)
{
  setup( data, dFESpace, comm, dFESpace->mapPtr(), (UInt)0 );

  M_rhs.reset                        ( new vector_Type(*M_localMap));
  M_residual_d.reset                 ( new vector_Type(*M_localMap));
  M_sxx.reset                        ( new vector_Type(*M_localMap) );
  M_syy.reset                        ( new vector_Type(*M_localMap) );
  M_szz.reset                        ( new vector_Type(*M_localMap) );
  M_linearSolver.reset               ( new SolverType( comm ) );
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>        data,
                                          const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                          boost::shared_ptr<Epetra_Comm>&     comm,
                                          const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                          UInt                                offset)
{
  M_data                            = data;
  M_FESpace                         = dFESpace;
  M_Displayer.reset                 (new Displayer(comm));
  M_me                              = comm->MyPID();
  M_elmatM.reset                    ( new MatrixElemental( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
  M_localMap                        = monolithicMap;
  M_disp.reset                      (new vector_Type(*M_localMap));
  //M_vel.reset                       (new vector_Type(*M_localMap));
  // M_acc.reset                       (new vector_Type(*M_localMap));
  // M_rhsW.reset                      ( new vector_Type(*M_localMap) );
  // M_rhsA.reset                      ( new vector_Type(*M_localMap) );
  // M_rhsNoBC.reset                   ( new vector_Type(*M_localMap) );

  M_mass.reset                      (new matrix_Type(*M_localMap));
  M_massTimeAdvanceCoefficient.reset(new matrix_Type(*M_localMap));
  M_tempMatrix.reset                (new matrix_Type(*M_localMap));
  M_jacobian.reset                  (new matrix_Type(*M_localMap));
  M_offset                          = offset;

  //M_theta                           = 2.0 * M_data->dataTime()->theta();
  //M_zeta                            = M_data->dataTime()->gamma();

  M_material.reset( material_Type::StructureMaterialFactory::instance().createObject( M_data->getSolidType()));
  M_material->setup(dFESpace,M_localMap,M_offset);
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( void )
{
  updateSystem(M_tempMatrix);
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( matrixPtr_Type& stiff )
{
    this->M_Displayer->leaderPrint(" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix
    M_material->computeMatrix(*this->M_disp, M_rescaleFactor, this->M_data, this->M_Displayer);


    //stiff.reset(new matrix_Type(*this->M_localMap));
    *stiff += *this->M_material->stiff();
    stiff->globalAssemble();

  /*
  //Matteo: timeAdvance method updates the update Right Hand Side;
   *this->M_rhsNoBC *= 0.0;

    computeRightHandSide();

  //  Real DeltaT    = this->M_data->dataTime()->timeStep();
    vector_Type z = *this->M_disp;

    z            +=  DeltaT*(*this->M_vel);

   //  Real coef;
    coef= (1.0 - this->M_zeta);

    *this->M_rhsNoBC  = *this->M_mass * z;
    *this->M_rhsNoBC -= (*stiff) * coef * (*this->M_disp);

    // acceleration rhs
    *this->M_rhsA = (2.0 / ( this->M_zeta * pow(DeltaT,2) )) * z + ((1.0 - this->M_zeta ) / ( this->M_zeta )) * (*M_acc);

  // velocity rhs

   *this->M_rhsW = *this->M_vel + ( 1 - this->M_theta  ) * DeltaT *  (*M_acc);

   std::cout << std::endl;

    std::cout << "rhsNoBC norm    = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "rhs_w   norm    = " << this->M_rhsW->norm2() << std::endl;
    std::cout << "    w   norm    = " << this->M_vel->norm2() << std::endl;
*/
    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());

}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSourceTerm( source_Type const& source )
{
   vector_Type rhs(vector_Type(*M_localMap));

    VectorElemental M_elvec(this->M_FESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
      {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
	  {
            compute_vec( source, M_elvec, this->M_FESpace->fe(),  this->M_data->dataTime()->getTime(), ic ); // compute local vector
            assembleVector( *rhs, M_elvec, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, ic*this->M_FESpace->getDim() ); // assemble local vector into global one
	  }
      }
   M_rhsNoBC +=rhs;
}

// Matteo this method isn't necessary it is replaced to updateRHS + updateSourceTerm;
/*
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( source_Type const& source )
{
   this->M_Displayer->leaderPrint(" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix

    M_material->computeNewMatrix(*this->M_disp, M_rescaleFactor, this->M_data, this->M_Displayer);

    M_tempMatrix.reset(new matrix_Type(*this->M_localMap));
    *M_tempMatrix += *this->M_material->stiff();
    M_tempMatrix->globalAssemble();

    *this->M_rhsNoBC *= 0.0;

    VectorElemental M_elvec(this->M_FESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

     // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->M_FESpace->mesh()->numVolumes(); ++i )
      {

        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
	  {
            compute_vec( source, M_elvec, this->M_FESpace->fe(),  this->M_data->dataTime()->getTime(), ic ); // compute local vector
            assembleVector( *this->M_rhsNoBC, M_elvec, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, ic*this->M_FESpace->getDim() ); // assemble local vector into global one
	  }
      }

   this->M_rhsNoBC->GlobalAssemble();

   computeRightHandSide();

}
*/
/*
//Matteo  this method isn't necessary. The method updateRHSContribution of timeAdvance method  does it;
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::computeRightHandSide( void )
{

    Real DeltaT    = this->M_data->dataTime()->timeStep();
    vector_Type z  = *this->M_disp;

    z             +=  DeltaT * (*this->M_vel);

    Real coef;
    coef= ( 1.0-M_zeta );

    *this->M_rhsNoBC += *this->M_mass * z;

    *this->M_rhsNoBC -= (*this->M_tempMatrix) * coef * (*this->M_disp);

    // acceleration rhs
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);

    // velocity rhs
    *this->M_rhsW = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);

    std::cout << std::endl;

    std::cout << "rhsNoBC norm    = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "rhs_w   norm    = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "    w   norm    = " << this->M_vel->norm2() << std::endl;
}

*/


template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::buildSystem( const Real& coefficient )
{
  M_Displayer->leaderPrint("  S-  Computing constant matrices ...          ");
  LifeChrono chrono;
  chrono.start();

  computeMassMatrix( coefficient );

  // Matteo  compute \xi_0 *Mass* \eta/(dt*dt)
  //*M_massTimeAdvanceCoefficient = *M_mass;
  //*M_massTimeAdvanceCoefficient *= /*timeAdvanceCoefficient*/1 / ( M_data->dataTime()->timeStep() * M_data->dataTime()->timeStep() );

  M_material->computeLinearStiffMatrix(this->M_data);

  chrono.stop();
  M_Displayer->leaderPrintMax( "done in ", chrono.diff() );
}


// template <typename Mesh, typename SolverType>
// void  StructuralSolver<Mesh, SolverType>::buildSystem(matrix_Type & bigMatrixStokes, const Real& timeAdvanceCoefficient, const Real& factor)
// {}
// ;


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::computeMassMatrix( const Real& factor)
{
  UInt totalDof = M_FESpace->dof().numTotalDof();

  // Number of displacement components
  UInt nc = nDimensions;

  //inverse of dt square:

  // Elementary computation and matrix assembling
  // Loop on elements
  for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); i++ )
    {

      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );

      M_elmatM->zero();

      // mass
      mass( factor * M_data->rho(), *M_elmatM, M_FESpace->fe(), 0, 0, nDimensions );

      // assembling
      for ( UInt ic = 0; ic < nc; ic++ )
      {
          //mass
          assembleMatrix( *M_mass, *M_elmatM, M_FESpace->fe(), M_FESpace->dof(),ic,ic, M_offset +  ic*totalDof, M_offset +  ic*totalDof);
      }
    }

  //getComunicator()->Barrier();

  M_mass->globalAssemble();

  //*massStiff *= factor; //M_data.dataTime()->timeStep() * M_rescaleFactor;
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::iterate( bchandler_Type& bch )
{
    LifeChrono chrono;

    // matrix and vector assembling communication
    this->M_Displayer->leaderPrint("  S-  Solving the system ... \n");

    this->M_BCh = bch;

    Real abstol  = 1.e-5;
    Real reltol  = 1.e-5;
    UInt maxiter = 50;
    Real etamax  = 0;
    Int NonLinearLineSearch = 0;

    Real time = this->M_data->time();

    Int status = 0;

    status = NonLinearRichardson( *this->M_disp, *this, abstol, reltol, maxiter, etamax, NonLinearLineSearch, this->M_out_res, this->M_data->dataTime()->time() );

    if ( status == 1 )
    {
        std::ostringstream ex;
        ex << "StructuralSolver::iterate() Inners nonLinearRichardson iterations failed to converge\n";
        throw std::logic_error( ex.str() );
    }
    else // if status == 0 NonLinearrRichardson converges
    {
        std::cout << std::endl;

        std::cout <<" Number of inner iterations       : " << maxiter <<  std::endl;

        std::cout <<" We are at the time step          : "  << this->M_data->dataTime()->time() << std::endl;

        this->M_out_iter << time << " " << maxiter << std::endl;
    }

   // updateVelAndAcceleration();

    std::cout << "iterate: d norm       = " << this->M_disp->norm2() << std::endl;
    //std::cout << "iterate: w norm       = " << this->M_vel->norm2() << std::endl;
    //std::cout << "iterate: a norm       = " << this->M_acc->norm2() << std::endl;

    //These two lines mut be checked fo FSI. With the linear solver, they have a totally
    //different expression. For structural problems it is not used.
    evalResidualDisplacement(*M_disp);

    //*this->M_residual_d = *this->M_mass*(*this->M_disp);
    //*this->M_residual_d -= *this->M_rhsNoBC;


}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::iterateLin( bchandler_Type& bch )
{
  LifeChrono chrono;

  matrixPtr_Type matrFull( new matrix_Type( *M_localMap, M_tempMatrix->meanNumEntries()));
\
  // matrix and vector assembling communication
  this->M_Displayer->leaderPrint("  S-  Solving the system in iteratLin... \n");

  // First Approximation: The Jacobian of P is equal to its linear part.
  //this->M_tempMatrix.reset(new matrix_Type(*this->M_localMap));
  //*matrFull += *this->M_material->linearStiff(); //it returns just the linear part
  //*matrFull *= M_zeta;
  //*matrFull += *this->M_mass; // Global Assemble is done inside BCManageMatrix
  ///End First Approximantion

  // Use of the complete Jacobian
  *matrFull += *this->M_jacobian;
 // Matteo:
   // *matrFull *= M_zeta;
  *matrFull += *this->M_massTimeAdvanceCoefficient; // Global Assemble is done inside BCManageMatrix

  this->M_Displayer->leaderPrint("\tS'-  Solving the linear system in iterateLin... \n");

  // for BC treatment (done at each time-step)

  this->M_Displayer->leaderPrint("\tS'-  Applying boundary conditions      ... ");

  vector_Type rhsFull (M_rhsNoBC->map());

  applyBoundaryConditionsLin( *matrFull, rhsFull, bch);


  this->M_Displayer->leaderPrintMax( "done in ", chrono.diff() );

  this->M_Displayer->leaderPrint("\tS'-  Solving system                    ... \n");
  chrono.start();

  this->M_linearSolver->setMatrix(*matrFull);

  this->M_linearSolver->solveSystem( rhsFull, *M_disp, matrFull );

  chrono.stop();

  //This line must be checked for FSI. In VenantKirchhoffSolver.hpp it has a
  //totally different expression.For structural problems it is not used
  evalResidualDisplacementLin(*M_disp);

}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::showMe( std::ostream& c  ) const
{
  c << "\n*** StructuralSolver::showMe method" << std::endl;
  c << "****** Data of the Material************" << std::endl;
  c << "Thickness:    " << M_data->thickness() << std::endl;
  c << "Density:      " << M_data->rho() << std::endl;
  c << "Young:        " << M_data->young(1) << std::endl;
  c << "Poisson:      " << M_data->poisson(1) << std::endl;
//  c << "Theta:        " << M_theta << std::endl;
//  c << "Zeta:         " << M_zeta << std::endl;
  c << "***************************************" << std::endl;
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::computeMatrix( matrixPtr_Type& stiff, const vector_Type& sol,  Real const& factor)
{
    this->M_Displayer->leaderPrint( " Computing residual ... \t\t\t");

    LifeChrono chrono;
    chrono.start();

    //It is right to do globalAssemble() inside the M_material class
    M_material->computeMatrix( sol, 1., this->M_data, M_Displayer);

    stiff.reset(new matrix_Type(*this->M_localMap));
    *stiff +=*this->M_material->stiff();
    // Matteo
    //*stiff *= M_zeta;
    *stiff += *this->M_massTimeAdvanceCoefficient;
    stiff->globalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );

}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidual( vector_Type &residual, const vector_Type& solution, Int /*iter*/)
{

    computeMatrix(this->M_tempMatrix, solution, 1.);

    this->M_Displayer->leaderPrint("    S- Updating the boundary conditions ... \t");
    LifeChrono chrono;
    chrono.start();
    if ( !this->M_BCh->bcUpdateDone() )
        this->M_BCh->bcUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    bcManageMatrix( *this->M_tempMatrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), 1.0 );

    vector_Type rhsFull(*this->M_rhsNoBC, Unique); // ignoring non-local entries, Otherwise they are summed up lately

    bcManageVector( rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(),  this->M_data->dataTime()->time(), 1.0 );

    *this->M_rhs = rhsFull;

    residual  = *this->M_tempMatrix*solution;
    residual -= *this->M_rhs;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidualDisplacement( const vector_Type& solution )
{

    computeMatrix(this->M_tempMatrix, solution, 1.);

    this->M_Displayer->leaderPrint("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    *this->M_residual_d  = *this->M_tempMatrix*solution;
    *this->M_residual_d -= *this->M_rhsNoBC;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidualDisplacementLin( const vector_Type& solution )
{

  //This is consisten with the previous first approximation in iterateLin
  this->M_tempMatrix.reset (new matrix_Type(*this->M_localMap));
  *this->M_tempMatrix += *this->M_material->linearStiff();
  // Matteo
  //*this->M_tempMatrix *= M_zeta;
  *this->M_tempMatrix += *this->M_massTimeAdvanceCoefficient;
  this->M_tempMatrix->globalAssemble();

  this->M_Displayer->leaderPrint("    S- Computing the residual displacement for the structure..... \t");
  LifeChrono chrono;
  chrono.start();

  //This definition of residual_d is similar to the one of iterateLin in VenantKirchhoffSolver
  *this->M_residual_d  = *this->M_tempMatrix*solution;

  chrono.stop();
  this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}



template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalConstraintTensor()
{
  vector_Type count(*M_localMap);

  *M_sxx *= 0.;
  *M_syy *= 0.;
  *M_szz *= 0.;

  for ( UInt ielem = 0; ielem < M_FESpace->mesh()->numVolumes(); ielem++ )
    {
      //UInt elem = M_FESpace->mesh()->volumeList( ielem ).id();
      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( ielem ) );

      //int    marker = M_FESpace->mesh()->volumeList( ielem ).marker();
      Real s      = 0;
      Real volume = M_FESpace->fe().detJac(0);

      for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
	{
        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	    {
	      Int i    = M_FESpace->fe().patternFirst(k);
	      Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	      s+= (2*M_data->mu(1) + M_data->lambda(1))*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s+= M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s+= M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 2 , ig )*
		(*M_disp)[idof + 2*M_FESpace->dim()];

	      count[idof]++;
	    }
	}

      for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	{
	  Int i    = M_FESpace->fe().patternFirst(k);
	  Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	  (*M_sxx)[idof] += s/M_FESpace->fe().detJac(0);
	}

      s = 0;

      for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
	{
        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	    {
	      Int i    = M_FESpace->fe().patternFirst(k);
	      Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	      s += M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += (2*M_data->mu(1) + M_data->lambda(1))*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 2 , ig )*
		(*M_disp)[idof + 2*M_FESpace->dim()];
	      //                         M_sxx[idof] += s;

	      //                M_syy[idof] += s/volume;
	    }
	}

      for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	{
	  Int i    = M_FESpace->fe().patternFirst(k);
	  Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	  (*M_syy)[idof] += s/volume;
	}


      s = 0;

      for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
	{
        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	    {
	      Int i    = M_FESpace->fe().patternFirst(k);
	      Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	      s += M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += M_data->lambda(1)*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += (2*M_data->mu(1) + M_data->lambda(1))*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 2 , ig )*
		(*M_disp)[idof + 2*M_FESpace->dim()];

	      //                         M_sxx[idof] += s;
	    }
	}

      for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
	{
	  Int i    = M_FESpace->fe().patternFirst(k);
	  Int idof = M_FESpace->dof().localToGlobal(M_FESpace->fe().currentLocalId(), i + 1);

	  (*M_szz)[idof] += s/M_FESpace->fe().detJac(0);
	}

    }

  for (UInt ii = 1; ii <= M_FESpace->dim(); ++ii)
    {
      (*M_sxx)[ii] /= count[ii];
      (*M_syy)[ii] /= count[ii];
      (*M_szz)[ii] /= count[ii];
    }
}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::initialize( vectorPtr_Type disp, vectorPtr_Type vel, vectorPtr_Type /*acc*/)
{
  *M_disp = *disp;
//  if (vel.get())
//    initializeVel(*vel);
}
/*
template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::initializeVel( const vector_Type& vel)
{
  *M_vel = vel;
}
*/

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::initialize( const Function& d0, const Function& w0, const Function& a0 )
{
  this->M_FESpace->interpolate(d0, *M_disp, 0.0);
  //this->M_FESpace->interpolate(w0, *M_vel , 0.0);
  // this->M_FESpace->interpolate(a0, *M_acc , 0.0);
}

/*
//Matteo this method isn't necessary timeAdvance compute the accelerate and velocity
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateVelAndAcceleration()
{
    Real DeltaT = this->M_data->dataTime()->timeStep();

    *this->M_acc = (2.0 /( this->M_zeta * pow(DeltaT,2) ))  * (*this->M_disp)  - *this->M_rhsA;
    *this->M_vel = *this->M_rhsW + this->M_theta * DeltaT * (*M_acc) ;
}
*/

template<typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::reduceSolution( Vector& displacement, Vector& velocity )
{
  vector_Type disp(*M_disp, 0);
  //vector_Type vel(*M_vel , 0);

  if ( getComunicator()->MyPID() == 0 )
    {
      for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
	{
	  disp[ iDof ] = displacement[ iDof + 1 ];
	  //vel [ iDof ] = velocity    [ iDof + 1 ];
	}
    }
}


// template <typename Mesh, typename SolverType>
// void
// StructuralSolver<Mesh, SolverType>::rescaleMatrices()
// {
//     *M_mass *=(M_data->dataTime()->timeStep()*M_rescaleFactor);
//     //*M_linearStiff *= (M_data->dataTime()->timeStep()*M_rescaleFactor);
// }

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::setDataFromGetPot( const GetPot& dataFile )
{
  M_linearSolver->setDataFromGetPot( dataFile, "solid/solver" );
  M_linearSolver->setupPreconditioner(dataFile, "solid/prec");

  UInt marker = M_FESpace->mesh()->volumeList( 1 ).marker();
  if (!M_data->young(marker))
    M_data->setYoung(dataFile( "solid/physics/young", 0. ), marker);
  if (!M_data->poisson(marker))
    M_data->setPoisson(dataFile( "solid/physics/poisson", 0. ), marker);
}


//Method UpdateJacobian
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateJacobian( vector_Type & sol, matrixPtr_Type& jacobian  )
{
    this->M_Displayer->leaderPrint("  S-  Solid: Updating JACOBIAN... ");

    LifeChrono chrono;
    chrono.start();

    M_material->updateJacobianMatrix(sol, this->M_data, this->M_Displayer);

    M_jacobian.reset(new matrix_Type(*this->M_localMap));
    *M_jacobian += *this->M_material->stiff();
    //This is necessary since the matrix has to be multiplied by a constant in iterateLin
    M_jacobian->globalAssemble();

    jacobian.reset(new matrix_Type(*this->M_localMap));
    *jacobian += *this->M_material->stiff();
     // Matteo
    //   *jacobian *= M_zeta;
    *jacobian += *this->M_massTimeAdvanceCoefficient;
    jacobian->globalAssemble();

    chrono.stop();
    this->M_Displayer->leaderPrintMax("   ... done in ", chrono.diff() );
}

//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::
solveJac( vector_Type& step, const vector_Type& res, Real& linear_rel_tol)
{
    solveJacobian(step,  res, linear_rel_tol, this->M_BCh);
}


//Method SolveJacobian
template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::
solveJacobian( vector_Type&           step,
               const vector_Type&     res,
               Real&                /*linear_rel_tol*/,
               bchandler_Type&        BCh)
{
    LifeChrono chrono;

    updateJacobian( *this->M_disp, this->M_tempMatrix );

    this->M_Displayer->leaderPrint("\tS'-  Solving the linear system ... \n");

    // for BC treatment (done at each time-step)

    this->M_Displayer->leaderPrint("\tS'-  Applying boundary conditions      ... ");

    this->M_rhsNoBC->globalAssemble();
   // this->M_rhsW->globalAssemble();

    vector_Type rhsFull (res);


//    bcManageMatrix( *this->M_tempMatrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), tgv );
    applyBoundaryConditions( *this->M_tempMatrix, rhsFull, BCh);


    this->M_Displayer->leaderPrintMax( "done in ", chrono.diff() );

    this->M_Displayer->leaderPrint("\tS'-  Solving system                    ... \n");
    chrono.start();

    this->M_linearSolver->setMatrix(*this->M_tempMatrix);

    this->M_linearSolver->solveSystem( rhsFull, step, this->M_tempMatrix );

    chrono.stop();

    //This line must be checked for FSI. In VenantKirchhoffSolver.hpp it has a
    //totally different expression.For structural problems it is not used
    *this->M_residual_d= *this->M_massTimeAdvanceCoefficient*step;

}


template<typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_Type&        matrix,
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
  vector_Type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

  //bcManage( matrix, rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *BCh, M_FESpace->feBd(), 1., M_data->dataTime()->time() );
  bcManageMatrix( matrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *BCh, this->M_FESpace->feBd(), 1., this->M_data->dataTime()->time() );

  // matrix should be GlobalAssembled by  bcManage

  rhs = rhsFull;

}

template<typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::applyBoundaryConditionsLin( matrix_Type&        matrix,
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
  vector_Type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

  bcManage( matrix, rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *BCh, M_FESpace->feBd(), 1., M_data->dataTime()->time() );

  // matrix should be GlobalAssembled by  bcManage

  rhs = rhsFull;

}



}
#endif
