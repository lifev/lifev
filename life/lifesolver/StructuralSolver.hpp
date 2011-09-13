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
 *  @brief This file contains solvers for different materials. WARNING!!!!This is the most important issue related with this class. At the moment, the BC are applied on the matrix and on rhsNoBc for VK models but for NH and EXP they are applied on the residual directly. This does not work for nonhomogeneus Dirichlet conditions!!
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _STRUCTURALSOLVER_H_
#define _STRUCTURALSOLVER_H_ 1

#include<boost/scoped_ptr.hpp>

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
  void updateSystem( matrixPtr_Type& mat_stiff, vectorPtr_Type& vec_stiff );

  //! Updates the system at the end of each time step given a source term
  /*!
    \param source volumic source
    \param time present time
  */
  void updateSystem( source_Type const& source );

  //! Comuptes the right hand side in the updateSystem methods
  void computeRHSNoBC( void );

  //! Compute the mass matrix and it calls the method to build the linear part of the stiffness matrix of the material class
  void buildSystem( void );

  //void buildSystem(matrix_Type & bigMatrixStokes); // used for monolithic

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
		 Real&            linear_rel_tol ) ;
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
  void updateVelAndAcceleration();

  //! Reduce the complete solution to the solution on the pocessor with rank 0
  /*!
    \param disp displacement solution
    \param vel velocity solution
  */
  void reduceSolution( Vector& displacement, Vector& velocity );

  //! Multiply the mass matrix and the linear stiffness matrix by the rescaleFactor
  void rescaleMatrices(); // used for monolithic

  /**
     in the linear case the solid matrix is constant, thus it does not need to be recomputed.
  */

  //! Update (in the case of nonlinear material) the solid matrix
  /*!
    \param stiff stiffness matrix
    \param sol the current solution
    \param factor the rescaleFactor
  */
  void computeMatrix( const vector_Type& sol, Real const& factor );

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
  vector_Type& velocity()         { return *M_vel; }

  //! Get the velocity
  vector_Type& acceleration()         { return *M_acc; }

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

  boost::scoped_ptr<Displayer>         M_Displayer;

  Int                                  M_me;

  //! data for solving tangent problem with aztec
  boost::shared_ptr<solver_Type>       M_linearSolver;

  //! Elementary matrices and vectors
  boost::shared_ptr<MatrixElemental>   M_elmatM;

  //! linearized velocity
  vectorPtr_Type                       M_disp;
  vectorPtr_Type                       M_vel;
  vectorPtr_Type                       M_acc;

  //! right  hand  side displacement
  vectorPtr_Type                       M_rhs;

  //! right  hand  side velocity
  vectorPtr_Type                       M_rhsW;

  //! right  hand  side velocity
  vectorPtr_Type                       M_rhsA;

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

  //! Matrix Temp: Temporary matrix to compute residuals or rhs
  matrixPtr_Type                       M_tempMatrix;

  //! Matrix Temp: Temporary auxiliary matrix to compute residuals or rhs
  matrixPtr_Type                       M_tempMatrixWithoutZeta;

  //! Jacobian Matrix: Matrix to store the jacobian of the newton method
  matrixPtr_Type                       M_jacobian;

  //! Stiffness vector for NH and Exp. It is used to compute residuals or rhs 
  vectorPtr_Type                       M_tempVect;

  //! Stiffness auxiliary vector for NH and Exp. It is used to compute residuals or rhs 
  vectorPtr_Type                       M_tempVectWithoutZeta;

  //! level of recursion for Aztec (has a sens with FSI coupling)
  UInt                                 M_recur;

  source_Type                          M_source;

  UInt                                 M_offset;
  Real                                 M_rescaleFactor;

  Real                                 M_zeta;
  Real                                 M_theta;

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
  M_vel                        ( ),
  M_acc                        ( ),
  M_rhs                        ( /*new vector_Type(M_localMap)*/),//useful
  M_rhsW                       ( ),
  M_rhsA                       ( ),
  M_rhsNoBC                    ( ),
  M_residual_d                 ( ),
  M_sxx                        (/*M_localMap*/),//useless
  M_syy                        (/*M_localMap*/),//useless
  M_szz                        (/*M_localMap*/),//useless
  M_out_iter                   ( "out_iter_solid" ),
  M_out_res                    ( "out_res_solid" ),
  M_BCh                        ( ),
  M_localMap                   ( ),
  M_mass                       ( ),
  M_tempMatrix                 ( ),
  M_tempMatrixWithoutZeta      ( ),
  M_jacobian                   ( ),
  M_tempVect                   ( ),
  M_tempVectWithoutZeta        ( ),
  M_recur                      ( ),
  M_source                     ( ),
  M_offset                     ( 0 ),
  M_rescaleFactor              ( 1. ),
  M_zeta                       ( 0.75 ),
  M_theta                      ( 0.7 ),
  M_material                   ( )
{
  std::cout << "I am in the constructor for the solver" << std::endl;
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>          data,
					      const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
					      bchandler_Type&                	BCh,
					      boost::shared_ptr<Epetra_Comm>&   comm)
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
  M_disp.reset                      ( new vector_Type(*M_localMap) );
  M_vel.reset                       ( new vector_Type(*M_localMap) );
  M_acc.reset                       ( new vector_Type(*M_localMap) );
  M_rhsW.reset                      ( new vector_Type(*M_localMap) );
  M_rhsA.reset                      ( new vector_Type(*M_localMap) );
  M_rhsNoBC.reset                   ( new vector_Type(*M_localMap) );
  M_mass.reset                      ( new matrix_Type(*M_localMap) );
  M_tempMatrix.reset                ( new matrix_Type(*M_localMap) );
  M_tempMatrixWithoutZeta.reset     ( new matrix_Type(*M_localMap) );
  M_jacobian.reset                  ( new matrix_Type(*M_localMap) );

  //Vector of Stiffness for NH and Exp
  //This vector stores the stiffness vector both in
  //updateSystem and in evalresidual. That's why it is called tempVect
  M_tempVect.reset                  (new vector_Type(*M_localMap));
  M_tempVectWithoutZeta.reset       (new vector_Type(*M_localMap));
  M_offset                          = offset;

  //M_theta                           = 2.0 * M_data->dataTime()->theta();
  //M_zeta                            = M_data->dataTime()->gamma();

  M_material.reset( material_Type::StructureMaterialFactory::instance().createObject( M_data->solidType() ) );
  M_material->setup( dFESpace,M_localMap,M_offset );
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( void )
{
  updateSystem(M_tempMatrix, M_tempVect);
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( matrixPtr_Type& mat_stiff, vectorPtr_Type& vec_stiff )
{
    this->M_Displayer->leaderPrint(" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix
    M_material->computeStiffness(*this->M_disp, M_rescaleFactor, this->M_data, this->M_Displayer);

    if ( this->M_data->solidType() == "linearVenantKirchhoff" || this->M_data->solidType() == "nonlinearVenantKirchhoff" )
      {
	mat_stiff.reset(new matrix_Type(*this->M_localMap));
	*mat_stiff += *this->M_material->stiffMatrix();
	mat_stiff->globalAssemble();
      }
    else
      {
	vec_stiff.reset(new vector_Type(*this->M_localMap));
	*vec_stiff += *this->M_material->stiffVector();
	vec_stiff->globalAssemble();
      }
    
    *this->M_rhsNoBC *= 0.0;

    computeRHSNoBC();

    std::cout << std::endl;

    std::cout << "rhsNoBC norm    = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "rhs_w   norm    = " << this->M_rhsW->norm2() << std::endl;
    std::cout << "    w   norm    = " << this->M_vel->norm2() << std::endl;

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff());

}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateSystem( source_Type const& source )
{
    this->M_Displayer->leaderPrint(" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix

    M_material->computeStiffness(*this->M_disp, M_rescaleFactor, this->M_data, this->M_Displayer);

    M_tempMatrix.reset(new matrix_Type(*this->M_localMap));
    *M_tempMatrix += *this->M_material->stiffMatrix();
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

    computeRHSNoBC();

}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::computeRHSNoBC( void )
{
  
    Real coef;
    coef= ( 1.0 - M_zeta );

    Real DeltaT    = this->M_data->dataTime()->timeStep();
    vector_Type z  = *this->M_disp;

    z             +=  DeltaT * (*this->M_vel);

    *this->M_rhsNoBC += *this->M_mass * z;

    //Are the multiplication optimized in Trilinos? Right I do as if
    //they are not!
    if ( this->M_data->solidType() == "linearVenantKirchhoff" || this->M_data->solidType() == "nonlinearVenantKirchhoff" )
	{
      	*this->M_rhsNoBC -= (*this->M_tempMatrix) * coef * (*this->M_disp);
	std::cout<<"\nMAT_stiff in UPSYSTEM";
	std::cout<<"\nrhsNoBC = "<<this->M_rhsNoBC->normInf()<<std::endl;
	}
    else
	{
      	*this->M_rhsNoBC -= (*this->M_tempVect)*coef;
	std::cout<<"\nVEC_stiff in UPSYSTEM";
	std::cout<<"\nrhsNoBC = "<<this->M_rhsNoBC->normInf()<<std::endl;
	}

    //! Acceleration right-hand side
    *M_rhsA = (2.0 / ( M_zeta * pow(DeltaT,2) )) * z + ((1.0 - M_zeta ) / ( M_zeta )) * (*M_acc);
    std::cout<<"\nrhsA = "<<this->M_rhsA->normInf()<<std::endl;

    //! Velocity right-hand side
    *this->M_rhsW = *this->M_vel + ( 1 - M_theta  ) * DeltaT *  (*M_acc);
    std::cout<<"\nrhsW = "<<this->M_rhsW->normInf()<<std::endl;

    std::cout << std::endl;

    std::cout << "rhsNoBC norm    = " << this->M_rhsNoBC->norm2() << std::endl;
    std::cout << "rhsA    norm    = " << this->M_rhsA->norm2()    << std::endl;
    std::cout << "rhsW    norm    = " << this->M_rhsW->norm2()    << std::endl;
    std::cout << "   W    norm    = " << this->M_vel->norm2() 	  << std::endl;

    std::cout << "\n   ZETA         = " << this->M_zeta  	   << std::endl;
    std::cout << "   THETA          = " << this->M_theta 	   << std::endl;
    std::cout << "   YOUNG          = " << this->M_data->young()   << std::endl;
    std::cout << "   POISSON        = " << this->M_data->poisson() << std::endl;
    std::cout << "   DENSITY        = " << this->M_data->rho()     << std::endl;
    std::cout << "   BULK MODULUS   = " << this->M_data->bulk()    << std::endl;
    std::cout << "   ALPHA          = " << this->M_data->alpha()   << std::endl;
    std::cout << "   GAMMA          = " << this->M_data->gamma()   << std::endl;
}



template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::buildSystem( void )
{
  M_Displayer->leaderPrint("  S-  Computing constant matrices ...          ");
  LifeChrono chrono;
  chrono.start();

  computeMassMatrix();
  M_material->computeLinearStiff(this->M_data);

  chrono.stop();
  M_Displayer->leaderPrintMax( "done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::computeMassMatrix( const Real& /*factor*/)
{
  UInt totalDof = M_FESpace->dof().numTotalDof();

  // Number of displacement components
  UInt nc = nDimensions;

  //inverse of dt square:
  Real dti2 = 2.0 / ( M_data->dataTime()->timeStep() * M_data->dataTime()->timeStep() );

  // Elementary computation and matrix assembling
  // Loop on elements
  for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); i++ )
    {

      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );

      M_elmatM->zero();

      // mass
      mass( dti2 * M_data->rho(), *M_elmatM, M_FESpace->fe(), 0, 0, nDimensions );

      // assembling
      for ( UInt ic = 0; ic < nc; ic++ )
      {
          //mass
          assembleMatrix( *M_mass, *M_elmatM, M_FESpace->fe(), M_FESpace->dof(),ic,ic, M_offset +  ic*totalDof, M_offset +  ic*totalDof);
      }
    }

  getComunicator()->Barrier();

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

    Real time = this->M_data->dataTime()->time();

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

    updateVelAndAcceleration();

    std::cout << "iterate: d norm       = " << this->M_disp->norm2() << std::endl;
    std::cout << "iterate: w norm       = " << this->M_vel->norm2() << std::endl;
    std::cout << "iterate: a norm       = " << this->M_acc->norm2() << std::endl;

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
  //*matrFull *= M_zeta;
  *matrFull += *this->M_mass; // Global Assemble is done inside BCManageMatrix
 
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
  c << "Young:        " << M_data->young() << std::endl;
  c << "Poisson:      " << M_data->poisson() << std::endl;
  c << "Theta:        " << M_theta << std::endl;
  c << "Zeta:         " << M_zeta << std::endl;
  c << "***************************************" << std::endl;
}

template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::computeMatrix( const vector_Type& sol,  Real const& /*factor*/)
{
    this->M_Displayer->leaderPrint( " Computing residual ... \t\t\t");

    LifeChrono chrono;
    chrono.start();

    //! It is right to do globalAssemble() inside the M_material class
    M_material->computeStiffness( sol, 1., this->M_data, this->M_Displayer);

    if ( this->M_data->solidType() == "linearVenantKirchhoff" || this->M_data->solidType() == "nonlinearVenantKirchhoff" )
      {
	M_tempMatrix.reset(new matrix_Type(*this->M_localMap));
	*M_tempMatrix +=*this->M_material->stiffMatrix()*M_zeta;
	*M_tempMatrix += *this->M_mass;
	M_tempMatrix->globalAssemble();
      }
    else
      {
	M_tempVect.reset(new vector_Type(*this->M_localMap));
	*M_tempVect = *this->M_material->stiffVector()*M_zeta;
	std::cout<< "\nVEC_stiff pre global = "<<M_tempVect->normInf()<<std::endl;

	M_tempVect->globalAssemble();
	std::cout<< "\nVEC_stiff post global = "<<M_tempVect->normInf()<<std::endl;
      }

    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidual( vector_Type &residual, const vector_Type& solution, Int /*iter*/)
{

    //This method call the M_material computeStiffness 
    computeMatrix(solution, 1.);

    this->M_Displayer->leaderPrint("    S- Updating the boundary conditions ... \t");
    LifeChrono chrono;
    
    //This is the most important issue related with this class.
    //At the moment, the BC are applied on the matrix and on rhsNoBc for VK models
    //but for NH and EXP they are applied on the residual directly. This does not work for
    //nonhomogeneus Dirichlet conditions!!

    if ( !this->M_BCh->bcUpdateDone() )
      this->M_BCh->bcUpdate( *this->M_FESpace->mesh(), this->M_FESpace->feBd(), this->M_FESpace->dof() );

    if ( this->M_data->solidType() == "linearVenantKirchhoff" || this->M_data->solidType() == "nonlinearVenantKirchhoff" )
      {
	chrono.start();

	bcManageMatrix( *this->M_tempMatrix, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), 1.0 );

	vector_Type rhsFull(*this->M_rhsNoBC, Unique); // ignoring non-local entries, Otherwise they are summed up lately

	bcManageVector( rhsFull, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(),  this->M_data->dataTime()->time(), 1.0 );
	
	*this->M_rhs = rhsFull;
    
	residual  = *this->M_tempMatrix * solution;
	residual -= *this->M_rhs;

	chrono.stop();
	this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
      }
    else //NH and Exp
      {
    	chrono.start();
	
    	residual = *this->M_mass * solution;
std::cout<<"\n RES1 NOBC = "<<residual.normInf()<<std::endl;
    	residual += *this->M_tempVect;
std::cout<<"\n RES2 NOBC = "<<residual.normInf()<<std::endl;
    	residual -= *this->M_rhsNoBC;
std::cout<<"\n RES3 NOBC = "<<residual.normInf()<<std::endl;

    	//! Apply the boundary conditions
    	bcManageVector( residual, *this->M_FESpace->mesh(), this->M_FESpace->dof(), *this->M_BCh, this->M_FESpace->feBd(), this->M_data->dataTime()->time(), 1.0);
std::cout<<"\n RES4 BC = "<<residual.normInf()<<std::endl;

    	chrono.stop();
    	this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );	
      }
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidualDisplacement( const vector_Type& solution )
{

    computeMatrix(solution, 1.);

    this->M_Displayer->leaderPrint("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    if ( this->M_data->solidType() == "linearVenantKirchhoff" || this->M_data->solidType() == "nonlinearVenantKirchhoff" )
      {    
	*this->M_residual_d  = *this->M_tempMatrix * solution;
	*this->M_residual_d -= *this->M_rhsNoBC;
      }
    else
      {
	*this->M_residual_d  = *this->M_tempVect;
	*this->M_residual_d -= *this->M_rhsNoBC;
      }
    chrono.stop();
    this->M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}


template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::evalResidualDisplacementLin( const vector_Type& solution )
{

  //This is consisten with the previous first approximation in iterateLin
  this->M_tempMatrix.reset (new matrix_Type(*this->M_localMap));
  *this->M_tempMatrix += *this->M_material->jacobian();
  //*this->M_tempMatrix *= M_zeta;
  *this->M_tempMatrix += *this->M_mass;
  this->M_tempMatrix->globalAssemble();

  this->M_Displayer->leaderPrint("    S- Computing the residual displacement for the structure..... \t");
  LifeChrono chrono;
  chrono.start();

  //This definition of residual_d is similar to the one of iterateLin in VenantKirchhoffSolver
  *this->M_residual_d  = *this->M_tempMatrix * solution;

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

	      s+= (2*M_data->mu() + M_data->lambda())*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s+= M_data->lambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s+= M_data->lambda()*
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

	      s += M_data->lambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += (2*M_data->mu() + M_data->lambda())*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += M_data->lambda()*
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

	      s += M_data->lambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += M_data->lambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += (2*M_data->mu() + M_data->lambda())*
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
  if (vel.get())
    initializeVel(*vel);
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::initializeVel( const vector_Type& vel)
{
  *M_vel = vel;
}

template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::initialize( const Function& d0, const Function& w0, const Function& a0 )
{
  this->M_FESpace->interpolate(d0, *M_disp, 0.0);
  this->M_FESpace->interpolate(w0, *M_vel , 0.0);
  this->M_FESpace->interpolate(a0, *M_acc , 0.0);
}


template <typename Mesh, typename SolverType>
void StructuralSolver<Mesh, SolverType>::updateVelAndAcceleration()
{
    Real DeltaT = this->M_data->dataTime()->timeStep();

    *this->M_acc = (2.0 /( this->M_zeta * pow(DeltaT,2) ))  * (*this->M_disp)  - *this->M_rhsA;
    *this->M_vel = *this->M_rhsW + this->M_theta * DeltaT * (*M_acc) ;
}


template<typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::reduceSolution( Vector& displacement, Vector& velocity )
{
  vector_Type disp(*M_disp, 0);
  vector_Type vel(*M_vel , 0);

  if ( getComunicator()->MyPID() == 0 )
    {
      for ( UInt iDof = 0; iDof < nDimensions*dim(); ++iDof )
	{
	  disp[ iDof ] = displacement[ iDof + 1 ];
	  vel [ iDof ] = velocity    [ iDof + 1 ];
	}
    }
}

/*
template <typename Mesh, typename SolverType>
void
StructuralSolver<Mesh, SolverType>::rescaleMatrices()
{

  *M_mass *=(M_data->dataTime()->timeStep()*M_rescaleFactor);
  *M_linearStiff *= (M_data->dataTime()->timeStep()*M_rescaleFactor);
}
*/
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

    std::cout << "\n   ZETA         = " << this->M_zeta  << std::endl;
    std::cout << "   THETA          = " << this->M_theta << std::endl;
    std::cout << "   YOUNG          = " << this->M_data->young() << std::endl;
    std::cout << "   POISSON        = " << this->M_data->poisson() << std::endl;
    std::cout << "   DENSITY        = " << this->M_data->rho() << std::endl;
    std::cout << "   BULK MODULUS   = " << this->M_data->bulk() << std::endl;
    std::cout << "   ALPHA          = " << this->M_data->alpha() << std::endl;
    std::cout << "   GAMMA          = " << this->M_data->gamma() << std::endl;

    M_material->updateJacobianMatrix(sol, this->M_data, this->M_Displayer);

    M_jacobian.reset(new matrix_Type(*this->M_localMap));
    *M_jacobian += *this->M_material->jacobian();
    //This is necessary since the matrix has to be multiplied by a constant in iterateLin
    M_jacobian->globalAssemble(); 

    jacobian.reset(new matrix_Type(*this->M_localMap));
    M_tempMatrixWithoutZeta.reset(new matrix_Type(*this->M_localMap));
    *M_tempMatrixWithoutZeta += *this->M_material->jacobian();
    M_tempMatrixWithoutZeta->globalAssemble();
    *M_tempMatrixWithoutZeta *= M_zeta;

    *jacobian += *M_tempMatrixWithoutZeta;

    *jacobian += *this->M_mass;
    jacobian->globalAssemble();

//std::string stringaP="M_jacobianSSPaolo";
//jacobian->spy(stringaP);

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
    this->M_rhsW->globalAssemble();

    vector_Type rhsFull (res);

//std::string stringaM_temp="M_tempMatrix";
//this->M_tempMatrix->spy(stringaM_temp);
    applyBoundaryConditions( *this->M_tempMatrix, rhsFull, BCh);
//std::string stringaPBC="M_jacobianSSPaoloBC";
//this->M_tempMatrix->spy(stringaPBC);

    this->M_Displayer->leaderPrintMax( "done in ", chrono.diff() );

    this->M_Displayer->leaderPrint("\tS'-  Solving system                    ... \n");
    chrono.start();

    this->M_linearSolver->setMatrix(*this->M_tempMatrix);

    this->M_linearSolver->solveSystem( rhsFull, step, this->M_tempMatrix );

    chrono.stop();

    //This line must be checked for FSI. In VenantKirchhoffSolver.hpp it has a 
    //totally different expression.For structural problems it is not used
    *this->M_residual_d = *this->M_mass*step;
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
