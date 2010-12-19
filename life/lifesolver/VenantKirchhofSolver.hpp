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
 *  @brief This file contains solvers for St. Venant-Kirchhof materials (linear for the moment)
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @version 1.1
 *  @date 01-09-2010
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _VENANTKIRCHHOFSOLVER_H_
#define _VENANTKIRCHHOFSOLVER_H_

#include<boost/scoped_ptr.hpp>

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

#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifecore/chrono.hpp>

#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataElasticStructure.hpp>
#include <life/lifecore/displayer.hpp>

namespace LifeV
{


/*!
  \class VenantKirchhofSolver
  \brief
  This class solves the linear elastodynamics equations for a (only linear right now)
  St. Venant-Kirchoff material

*/
template <typename Mesh,
	  typename SolverType = LifeV::SolverTrilinos >

class VenantKirchhofSolver
{
public:

  //!@name Type definitions
  //@{
  typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
  typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_Type;

  typedef BCHandler                              bchandlerRaw_Type;
  typedef boost::shared_ptr<bchandlerRaw_Type>   bchandler_Type;

  typedef SolverType                             solver_Type;

  typedef typename solver_Type::matrix_type      matrix_Type;
  typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;
  typedef typename solver_Type::vector_type      vector_Type;
  typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

  typedef typename SolverType::prec_raw_type     precRaw_Type;
  typedef typename SolverType::prec_type         prec_Type;

  typedef DataElasticStructure                   data_Type;

  typedef singleton<factory<VenantKirchhofSolver,  std::string> >  StructureSolverFactory;

  //@}


  //! @name Constructor
  //@{

  VenantKirchhofSolver();

  //@}

  //!@name Methods
  //@{

  //! Setup the created object of the class Venantkirchhof
  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param BCh boundary conditions for the displacement
    \param comm
  */
  void setup( boost::shared_ptr<data_Type> data,
	      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   FESpace,
	      bchandler_Type&       BCh,
	      boost::shared_ptr<Epetra_Comm>&     comm
	      );

  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param comm
  */
  void setup( boost::shared_ptr<data_Type> data,
	      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   FESpace,
	      boost::shared_ptr<Epetra_Comm>&     comm
	      );

  /*!
    \param data_file GetPot data file
    \param refFE reference FE for the displacement
    \param comm the comunicator parameter
    \param monolithicMap the EpetraMap
    \param offset the offset parameter
  */
  virtual void setup( boost::shared_ptr<data_Type> data,
		      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >&   dFESpace,
		      boost::shared_ptr<Epetra_Comm>&     comm,
		      const boost::shared_ptr<const EpetraMap>&       monolithicMap,
		      UInt       offset=0
		      //boost::shared_ptr<FESpace<Mesh, EpetraMap> >   uFESpace=0
		      );


  //! Updates the system at the end of each time step given a source term
  /*!
    \param source volumic source
    \param time present time
  */
  //void updateSystem( source_Type const& source );

  //! Updates the system at the end of each time step
  void updateSystem();

  //! Updates the system at the end of each time step when the matrix is passed from outside
  /*!
    \param stiff stiffness matrix provided from outside
  */
  virtual void updateSystem(matrixPtr_Type& stiff);

  //! Compute the mass matrix and the linear part of the stiffness matrix
  void buildSystem();

  //void buildSystem(matrix_Type & bigMatrixStokes); // used for monolithic

  //! Compute the mass matrix and the linear part of the stiffness matrix
  /*!
    \param matrix the matrix containing the mass matrix and the linear part of he stiffness matrix
    \param factor rescale matrices factor used in monolithic
  */
  virtual void buildSystem(matrixPtr_Type matrix, const Real& factor=1.);

  //! Solve the non-linear system
  /*! They compute the solution solving the non linear system
    \param sol vector containing the solution at the current time
  */
  void iterate( vector_Type& solution );

  //! Solve the non-linear system
  /*!
    \param bch BCHander object containing the boundary conditions
  */
  virtual void iterate( bchandler_Type& bch );

  //! Solve the non-linear system
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
  virtual void updateJacobian( vector_Type& solution, matrixPtr_Type& jacobian )=0;

  //! Solves the tangent problem for newton iterations
  /*!
    \param step the vector containing the solution of the sistem J*step=-Res
    \param res the vector conteining the residual
    \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
           eta is determined by the modified Eisenstat-Walker formula
  */
  virtual void solveJac( vector_Type&       step,
			 const vector_Type& residual,
			 Real&            linear_rel_tol)=0 ;
  //    void solveJac( const Vector& res, Real& linear_rel_tol, Vector &step);

  //! Solves the tangent problem with custom BC
  /*!
    \param step the vector containing the solution of the sistem J*step=-Res
    \param res the vector conteining the residual
    \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
           eta is determined by the modified Eisenstat-Walker formula
    \param BCd BCHandler object containing the boundary condition
  */
  virtual void solveJacobian( vector_Type&       step,
			      const vector_Type& residual,
			      Real&            linear_rel_tol,
			      bchandler_Type&    BCd )=0 ;


  //! Evaluates residual for newton interations
  /*!
    \param res residal vector that is update every time the method is called
    \param sol solution vector from which the residual is computed
    \param iter iteration of the nonLinearRichardson method
  */
  void evalResidual( vector_Type &residual, const vector_Type& solution, Int iter);

  void evalConstraintTensor();

  //! Sets the initial displacement, velocity, acceleration
  /*!
    \param d0 space function describing the initial displacement
    \param w0 space function describing the initial velocity
    \param a0 space function describing the initial acceleration
  */
  virtual void initialize( const Function& d0, const Function& w0, const Function& a0 = Function() );

  //! Sets the initial displacement, velocity, acceleration
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

  //! Computes the velocity vector at the n-th time step
  virtual void updateVel();

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
  void computeMatrix( matrixPtr_Type& stiff, const vector_Type& sol, Real const& factor )
  {
  }
  //! Update (in the case of nonlinear material) the solid matrix
  /*!
    \param sol the current solution
    \param factor the rescaleFactor
  */
  void computeMatrix( const vector_Type& sol, Real const& factor )
  {
  }

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
  EpetraMap   const& getMap()       const { return *M_localMap; }

  //! Get the Displayer object
  Displayer   const& getDisplayer() const { return *M_Displayer; }

  //! Get the matrix containing the mass mtrix and the linear part of the stiffness matrix
  matrixPtr_Type const getMassStiff() const {return M_massStiff; }

  //! Get the mass matrix
  matrixPtr_Type const getMass() const {return M_mass; }

  //! Get the linear part of the stiffness matrix
  matrixPtr_Type const getLinearStiff() const {return M_linearStiff; }

  //! BCHandler getter and setter
  //    LIFEV_DEPRECATED BCHandler const & BC_solid() const {return BCh_solid();}

  //! Get the FESpace object
  FESpace<Mesh, EpetraMap>& getDFESpace() {return M_FESpace;}

  //! Get the bCHandler object
  bchandler_Type const & getBChandler() const {return M_BCh;}

  //! Get the residual
  vector_Type& getResidual()             {return *M_residual_d;}

  //! Get the source term
  source_Type const& getSourceTerm() const { return M_source; }

  //! Get the displacement
  vector_Type& getDisplacement()        { return *M_disp; }

  //! Get the velocity
  vector_Type& getVelocity()         { return *M_vel; }

  //! Get the right hand sde without BC

  vectorPtr_Type& getRhsWithoutBC() { return M_rhsNoBC; }

  //const Dof& dDof() const { return M_FESpace.dof(); }

  //const Mesh& mesh() const { return M_FESpace.mesh(); }

  //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

  //! Get the comunicator object
  boost::shared_ptr<Epetra_Comm> const& getComunicator() const {return M_Displayer->comm();}

  //! Get the rescaleFactor
  Real getRescaleFactor() {return M_rescaleFactor;}

  /*! Get the offset parameter. It is taken into account when the boundary conditions
    are applied and the matrices are assembled.
  */
  const UInt& getOffset() const { return M_offset; }

  /**
     Do nothing in the linear case: the matrix remains constant. Otherwise substitute the matrix with an updated one
  */
  //! Get the Solid Matrix
  void getSolidMatrix( matrixPtr_Type& matrix)
  {
  }

  // Physic constant
  //! Get the thickness
  const Real& getThickness() const { return M_data->getThickness(); }

  //! Get the density
  //  const Real& density()   const { return M_data->rho(); }

  //! Get the Young modulus
  const Real& getYoung()     const { return M_data->getYoung(); }

  //! Get the Poisson coefficient
  const Real& getPoisson()   const { return M_data->getPoisson(); }

  //! Get the density
  const Real& getRho()       const { return M_data->getRho(); }

  //@}

  //Deprecated Interfaces

  //! Get the density //Deprecated
  const Real& __attribute__ ((__deprecated__)) rho() const
  { return getRho(); }

  //! Get the Poisson coefficient  //Deprecated
  const Real& __attribute__ ((__deprecated__)) poisson() const
  { return getPoisson(); }

  //! Get the Young modulus //Deprecated
  const Real& __attribute__ ((__deprecated__)) young() const
  { return getYoung(); }

  //! Get the thickness //Deprecated
  const Real& __attribute__ ((__deprecated__)) thickness() const
  { return getThickness(); }

  /*! Get the offset parameter. It is taken into account when the boundary conditions
    are applied and the matrices are assembled.
  */ //Deprecated

    const UInt& __attribute__ ((__deprecated__)) offset() const
  { return getOffset(); }

  //! Get the comunicator object //Deprecated
  boost::shared_ptr<Epetra_Comm> const& __attribute__ ((__deprecated__)) comm() const
  {return getComunicator();}
 
protected:

  //! Apply boundary condition
  /*!
    \param matrix the matrix of the system
    \param rhs the right hand side of the system
    \param BCh BCHandler object
    \param offset the offset parameter
  */
  virtual void applyBoundaryConditions(matrix_Type &matrix,
				       vector_Type &rhs,
				       bchandler_Type& BCh,
				       UInt         offset=0);

  UInt getDim() const { return M_FESpace->dim(); }

  //Deprecated
  UInt __attribute__ ((__deprecated__))dim() const
  { return getDim(); }


  //!Protected Members

  boost::shared_ptr<data_Type>   M_data;

  boost::shared_ptr<FESpace<Mesh, EpetraMap> >      M_FESpace;

  boost::scoped_ptr<Displayer>   M_Displayer;

  Int                            M_me;

  //! data for solving tangent problem with aztec
  boost::shared_ptr<solver_Type>                    M_linearSolver;

  //! Elementary matrices and vectors
  boost::shared_ptr<ElemMat>                        M_elmatK; // stiffnes
  boost::shared_ptr<ElemMat>                        M_elmatM; // mass
  boost::shared_ptr<ElemMat>                        M_elmatC; // mass + stiffness
  //    ElemVec                        M_elvec;  // Elementary right hand side
  //    ElemVec                        M_dk_loc; // Local displacement

  //! linearized velocity

  vectorPtr_Type                    M_disp;
  vectorPtr_Type                    M_vel;

  //! right  hand  side displacement
  vectorPtr_Type                    M_rhs;

  //! right  hand  side velocity

  vectorPtr_Type                    M_rhsW;

  //! right  hand  side
  vectorPtr_Type                    M_rhsNoBC;

  //! right  hand  side
  boost::shared_ptr<vector_Type>                    M_f;

  //! residual
  boost::shared_ptr<vector_Type>                    M_residual_d;

  //    vector_Type*                   M_sxx;

  vectorPtr_Type                    M_sxx;
  vectorPtr_Type                    M_syy;
  vectorPtr_Type                    M_szz;

  //! files for lists of iterations and residuals per timestep
  std::ofstream                  M_out_iter;
  std::ofstream                  M_out_res;

  bchandler_Type   M_BCh;

  boost::shared_ptr<const EpetraMap>                      M_localMap;

  //! Matrix M: mass
  matrixPtr_Type                    M_mass;

  //! Matrix C: mass + linear stiffness
  matrixPtr_Type                    M_massStiff;

  //! Matrix Knl: stiffness non-linear
  matrixPtr_Type                    M_stiff;

  //! Matrix Kl: stiffness linear
  matrixPtr_Type                    M_linearStiff;

  //! Matrix J: jacobian
  matrixPtr_Type                    M_jacobian;

  //! level of recursion for Aztec (has a sens with FSI coupling)
  UInt                              M_recur;

  source_Type                    M_source;

  Int                            M_count;

  UInt                           M_offset;
  Real                           M_rescaleFactor;

  matrixPtr_Type                  M_matrFull;

};

//====================================
// Constructor
//=====================================

template <typename Mesh, typename SolverType>
VenantKirchhofSolver<Mesh, SolverType>::VenantKirchhofSolver( ):
  M_data                       ( ),
  M_FESpace                    ( ),
  M_Displayer                  ( ),
  M_me                         ( 0 ),
  M_linearSolver               ( ),
  M_elmatK                     ( ),
  M_elmatM                     ( ),
  M_elmatC                     ( ),
  M_disp                       ( ),
  M_vel                        ( ),
  M_rhs                        ( /*new vector_Type(M_localMap)*/),//useful
  M_rhsW                       ( ),
  M_rhsNoBC                    ( ),
  M_f                          ( ),//useless
  M_residual_d                 ( ),//useless
  M_sxx                        (/*M_localMap*/),//useless
  M_syy                        (/*M_localMap*/),//useless
  M_szz                        (/*M_localMap*/),//useless
  M_out_iter                   ( "out_iter_solid" ),
  M_out_res                    ( "out_res_solid" ),
  M_BCh                        (),
  M_localMap                   ( ),
  M_mass                       ( ),
  M_massStiff                  ( /*new matrix_Type(monolithicMap) */),//constructed outside
  M_stiff                      ( /*new matrix_Type(M_localMap)*/ ),
  M_linearStiff                ( ),
  M_jacobian                   (/*M_localMap*/),
  M_recur                       (),
  M_source                     (),
  M_count                      ( 0 ),//useless
  M_offset                     ( 0 ),
  M_rescaleFactor              ( 1. ),
  M_matrFull                   ( )//useless
{
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>          data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      bchandler_Type&                BCh,
					      boost::shared_ptr<Epetra_Comm>&              comm)
{
  setup(data, dFESpace, comm);
  M_BCh = BCh;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>        data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      boost::shared_ptr<Epetra_Comm>&     comm)
{
  setup( data, dFESpace, comm, dFESpace->mapPtr(), (UInt)0 );
  M_stiff.reset                      ( new matrix_Type(*M_localMap) );
  M_massStiff.reset                  ( new matrix_Type(*M_localMap) );
  M_jacobian.reset                   ( new matrix_Type(*M_localMap) );
  M_rhs.reset                        ( new vector_Type(*M_localMap));
  M_f.reset                          ( new vector_Type(*M_localMap));
  M_residual_d.reset                 ( new vector_Type(*M_localMap));
  M_sxx.reset                        ( new vector_Type(*M_localMap) );
  M_syy.reset                        ( new vector_Type(*M_localMap) );
  M_szz.reset                        ( new vector_Type(*M_localMap) );
  M_linearSolver.reset               ( new SolverType( comm ) );

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setup(boost::shared_ptr<data_Type>        data,
					      const boost::shared_ptr< FESpace<Mesh, EpetraMap> >& dFESpace,
					      boost::shared_ptr<Epetra_Comm>&     comm,
					      const boost::shared_ptr<const EpetraMap>&  monolithicMap,
					      UInt                                offset)
{
  M_data                            = data;
  M_FESpace                         = dFESpace;
  M_Displayer.reset                 (new Displayer(comm));
  M_me                              = comm->MyPID();
  M_elmatK.reset                    ( new ElemMat( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
  M_elmatM.reset                    ( new ElemMat( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
  M_elmatC.reset                    ( new ElemMat( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
  M_localMap                        = monolithicMap;

  M_disp.reset                      (new vector_Type(*M_localMap));
  M_vel.reset                       (new vector_Type(*M_localMap));
  M_rhsW.reset                      ( new vector_Type(*M_localMap) );
  M_rhsNoBC.reset                   ( new vector_Type(*M_localMap) );
  M_mass.reset                      (new matrix_Type(*M_localMap));
  M_linearStiff.reset               (new matrix_Type(*M_localMap));
  M_offset                          = offset;
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::updateSystem(  )
{
  updateSystem(M_stiff);
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::updateSystem( matrixPtr_Type& /*stiff*/ )
{
  M_Displayer->leaderPrint("  S-  Updating mass term on right hand side... ");

  Chrono chrono;
  chrono.start();

  //M_stiff = M_linearStiff;

  // Number of displacement components
  //    UInt nc = nDimensions;

  Real coef;
  coef = (Real) M_data->getDataTime()->timeStep();


  vector_Type z = *M_disp;
  z             += coef*(*M_vel);

  std::cout<< "M_disp in solid" << M_disp->norm2()<<std::endl;

  *M_rhsNoBC  = *M_mass * z;

  std::cout<< "rhsNoBC in solid 1" << M_rhsNoBC->norm2()<<std::endl;

  *M_rhsNoBC -= *M_linearStiff * (*M_disp);

  std::cout<< "rhsNoBC in solid 2" << M_rhsNoBC->norm2()<<std::endl;

  coef = 2.0/M_data->getDataTime()->timeStep();
  *M_rhsW  = coef * (*M_disp);
  *M_rhsW += *M_vel;

  chrono.stop();

  M_Displayer->leaderPrintMax("done in ", chrono.diff());
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::buildSystem()
{
  M_Displayer->leaderPrint("  S-  Computing constant matrices ...          ");
  Chrono chrono;
  chrono.start();

  M_massStiff.reset(new matrix_Type(*M_localMap));
  buildSystem(M_massStiff);
  M_massStiff->globalAssemble();

  chrono.stop();
  M_Displayer->leaderPrintMax( "done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::buildSystem(matrixPtr_Type massStiff, const Real& factor)
{
  UInt totalDof = M_FESpace->dof().numTotalDof();

  // Number of displacement components
  UInt nc = nDimensions;

  //inverse of dt square:
  Real dti2 = 2.0 / ( M_data->getDataTime()->timeStep() * M_data->getDataTime()->timeStep() );

  // Elementary computation and matrix assembling
  // Loop on elements
  for ( UInt i = 1; i <= M_FESpace->mesh()->numVolumes(); i++ )
    {

      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );

      M_elmatK->zero();
      M_elmatM->zero();


      // stiffness
      UInt marker = M_FESpace->mesh()->volumeList( i ).marker();
      stiff_strain(     M_data->getMu(marker), *M_elmatK, M_FESpace->fe() );
      stiff_div   ( 0.5*M_data->getLambda(marker), *M_elmatK, M_FESpace->fe() );

      M_elmatC->mat() = M_elmatK->mat();

      // mass
      mass( dti2 * M_data->getRho(), *M_elmatM, M_FESpace->fe(), 0, 0, nDimensions );

      M_elmatC->mat() += M_elmatM->mat();

      // assembling
      for ( UInt ic = 0; ic < nc; ic++ )
	{
	  for ( UInt jc = 0; jc < nc; jc++ )
	    {
	      assembleMatrix( *M_linearStiff, *M_elmatK, M_FESpace->fe(), M_FESpace->fe(), M_FESpace->dof(), M_FESpace->dof(),  ic,  jc,  M_offset +ic*totalDof, M_offset + jc*totalDof );

	      assembleMatrix( *massStiff  , *M_elmatC, M_FESpace->fe(), M_FESpace->fe(), M_FESpace->dof(), M_FESpace->dof(),  ic,jc, M_offset +ic*totalDof,M_offset + jc*totalDof );
	    }

	  //mass
	  assembleMatrix( *M_mass, *M_elmatM, M_FESpace->fe(), M_FESpace->dof(),ic,ic, M_offset +  ic*totalDof, M_offset +  ic*totalDof);
	}
    }

  getComunicator()->Barrier();

  M_linearStiff->globalAssemble();
  massStiff->globalAssemble();
  M_mass->globalAssemble();
  *massStiff *= factor; //M_data.dataTime()->timeStep() * M_rescaleFactor;
}

template <typename Mesh, typename SolverType>
void VenantKirchhofSolver<Mesh, SolverType>::iterate(vector_Type& solution)
{
  Int status;

  vector_Type sol(*M_localMap);

  *M_vel  = ( 2.0 / M_data->getDataTime()->timeStep() ) * (*M_disp);
  *M_vel -= *M_rhsW;


  //M_Displayer->leaderPrint("sol norm = ", norm(this->sol));

  *M_residual_d  = M_massStiff*sol;
  *M_residual_d -= *M_rhsNoBC;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::iterate( bchandler_Type& bch )
{
  // matrix and vector assembling communication
  matrixPtr_Type matrFull( new matrix_Type( *M_localMap, M_massStiff->meanNumEntries()));

  *matrFull += *M_massStiff;

  M_rhsNoBC->globalAssemble();
  M_rhsW->globalAssemble();

  vector_Type rhsFull (*M_rhsNoBC);

  // boundary conditions update
  //M_comm->Barrier();

  M_Displayer->leaderPrint("  S-  Applying boundary conditions ...         ");
  Chrono chrono;
  chrono.start();

  applyBoundaryConditions( *matrFull, rhsFull, bch);

  chrono.stop();
  M_Displayer->leaderPrintMax("done in " , chrono.diff());


  // solving the system
  M_linearSolver->setMatrix(*matrFull);

  M_linearSolver->solveSystem( rhsFull, *M_disp, matrFull);

  // computing the velocity vector and the residual

  *M_vel  = ( 2.0 / M_data->getDataTime()->timeStep() ) * (*M_disp);
  *M_vel -= *M_rhsW;

  *M_residual_d =  *M_massStiff * (*M_disp);
  *M_residual_d -= *M_rhsNoBC;

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::iterateLin( bchandler_Type& bch )
{
  matrixPtr_Type matrFull( new matrix_Type( *M_localMap, M_massStiff->meanNumEntries()));

  *matrFull += *M_massStiff;

  M_rhsNoBC->globalAssemble();
  M_rhsW->globalAssemble();

  vector_Type rhsFull (M_rhsNoBC->map());

  M_Displayer->leaderPrint(" LS-  Applying boundary conditions ...         ");
  Chrono chrono;
  chrono.start();

  // boundary conditions update
  applyBoundaryConditions( *matrFull, rhsFull, bch);

  chrono.stop();
  M_Displayer->leaderPrintMax("done in " , chrono.diff());

  //M_Displayer->leaderPrint("rhs_dz norm = ", rhsFull.NormInf() );

  M_linearSolver->setMatrix(*matrFull);
  Int numIter = M_linearSolver->solveSystem(  rhsFull, *M_disp, matrFull );

  //M_Displayer->leaderPrintMax("dz norm     = " , M_disp.NormInf() );

  numIter = abs(numIter);

  *M_residual_d =  *M_massStiff*(*M_disp);
  //    M_residual_d -= M_rhsNoBC;

}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::showMe( std::ostream& c  ) const
{
  c << "\n*** VenantKirchhof::showMe method" << std::endl;
  c << "****** Data of the Material************" << std::endl;
  c << "Thickness:   " << M_data->getThickness();
  c << "Density:   " << M_data->getRho();
  c << "Young:   " << M_data->getYoung();
  c << "Poisson:   " << M_data->getPoisson();
  c << "***************************************" << std::endl;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::evalResidual( vector_Type &residual, const vector_Type& solution, Int /*iter*/)
{
  M_Displayer->leaderPrint("  S-  Computing residual ...                   ");
  Chrono chrono;
  chrono.start();

  // Matrices initialization
  M_stiff = M_massStiff;

  M_Displayer->leaderPrint("updating the boundary conditions ... ");
  if ( !M_BCh->bcUpdateDone() )
    M_BCh->bcUpdate( M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );

  bcManageMatrix( M_stiff, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), 1.0 );

  *M_rhs = *M_rhsNoBC;

  bcManageVector( *M_rhs, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), M_data->getDataTime()->time(), 1.0 );
  residual  = M_stiff * solution;
  //    res -= M_rhs;

  chrono.stop();
  M_Displayer->leaderPrintMax("done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::evalConstraintTensor()
{
  vector_Type count(*M_localMap);

  *M_sxx *= 0.;
  *M_syy *= 0.;
  *M_szz *= 0.;

  for ( UInt ielem = 1; ielem <= M_FESpace->mesh()->numVolumes(); ielem++ )
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

	      s+= (2*M_data->getMu() + M_data->getLambda())*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s+= M_data->getLambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s+= M_data->getLambda()*
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

	      s += M_data->getLambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += (2*M_data->getMu() + M_data->getLambda())*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += M_data->getLambda()*
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

	      s += M_data->getLambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 0 , ig )*
		(*M_disp)[idof + 0*M_FESpace->dim()];

	      s += M_data->getLambda()*
		M_FESpace->fe().weightDet( ig )*
		M_FESpace->fe().phiDer( k, 1 , ig )*
		(*M_disp)[idof + 1*M_FESpace->dim()];

	      s += (2*M_data->getMu() + M_data->getLambda())*
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
VenantKirchhofSolver<Mesh, SolverType>::initialize( vectorPtr_Type disp, vectorPtr_Type vel, vectorPtr_Type /*acc*/)
{
  *M_disp = *disp;
  if (vel.get())
    initializeVel(*vel);
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::initializeVel( const vector_Type& vel)
{
  *M_vel = vel;
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::initialize( const Function& d0, const Function& w0, const Function& /*a0*/ )
{
  M_FESpace->interpolate(d0, *M_disp, 0.0);
  M_FESpace->interpolate(w0, *M_vel , 0.0);
}


template <typename Mesh, typename SolverType> // for monolithic
void
VenantKirchhofSolver<Mesh, SolverType>::updateVel()
{
  *M_vel  = ( 2.0 / M_data->getDataTime()->timeStep() ) * (*M_disp);
  *M_vel -= *M_rhsW;
}

template<typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::reduceSolution( Vector& displacement, Vector& velocity )
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

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::rescaleMatrices()
{
  *M_mass *=(M_data->getDataTime()->timeStep()*M_rescaleFactor);
  *M_linearStiff *= (M_data->getDataTime()->timeStep()*M_rescaleFactor);
}

template <typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::setDataFromGetPot( const GetPot& dataFile )
{
  M_linearSolver->setDataFromGetPot( dataFile, "solid/solver" );
  M_linearSolver->setUpPrec(dataFile, "solid/prec");

  UInt marker = M_FESpace->mesh()->volumeList( 1 ).marker();
  if (!M_data->getYoung(marker))
    M_data->setYoung(dataFile( "solid/physics/young", 0. ), marker);
  if (!M_data->getPoisson(marker))
    M_data->setPoisson(dataFile( "solid/physics/poisson", 0. ), marker);
}


template<typename Mesh, typename SolverType>
void
VenantKirchhofSolver<Mesh, SolverType>::applyBoundaryConditions( matrix_Type&        matrix,
								 vector_Type&        rhs,
								 bchandler_Type&     BCh,
								 UInt                offset)
{
  // BC manage for the velocity
  if (offset)
    BCh->setOffset(offset);
  if ( !BCh->bcUpdateDone() )
    BCh->bcUpdate( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );

  // vector_Type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
  vector_Type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

  bcManage( matrix, rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *BCh, M_FESpace->feBd(), 1.,
              M_data->getDataTime()->time() );

  // matrix should be GlobalAssembled by  bcManage

  rhs = rhsFull;

}


}

#endif
