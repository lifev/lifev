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
  \file SecondOrderSolver.h
  \author M.A. Fernandez
  \date 6/2003
  \version 1.0

  \brief
  This file contains solvers for St. Venant-Kirchhof materials (linear for the moment)

*/
#ifndef _SecondOrderSolver_H_
#define _SecondOrderSolver_H_

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

#include <life/lifealg/dataNewton.hpp>
#include <life/lifealg/newton.hpp>
#include <life/lifealg/SolverTrilinos.hpp>

#include <life/lifesolver/dataSecondOrder.hpp>
#include <life/lifecore/displayer.hpp>



#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

namespace LifeV
{


/*!
  \class SecondOrderSolver
  \brief
This class solves the following waves equation:

 M \frac{d^2 u}{d t^2} + D(u,\frac{d u}{d t} ) \frac{\partial u}{\partial t} + A(u, frac{d u}{dt}) = f 

where M is mass matrix, A  stiffness matrix and D is  damping matrix given by alpha M + beta A.

 If A and D depend on u and du we  linearize the problem using suitable extrapolations,
see timeAdvance.pdf notes.
*/

template <typename Mesh,
          typename SolverType = LifeV::SolverTrilinos >

class SecondOrderSolver
{
public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;


    typedef SolverType                    solver_type;

    typedef typename solver_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>         matrix_ptrtype;
    typedef typename solver_type::vector_type      vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;


    typedef typename SolverType::prec_raw_type     prec_raw_type;
    typedef typename SolverType::prec_type         prec_type;

    typedef DataSecondOrder                                      data_type;

    //! Constructors
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
      \param BCh boundary conditions for the displacement
    */

    SecondOrderSolver( const data_type& data,
                          FESpace<Mesh, EpetraMap>&   FESpace,
                          BCHandler&       BCh,
                          boost::shared_ptr<Epetra_Comm>&     comm,
                          UInt             offset=0);

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param Qr volumic quadrature rule
      \param bdQr surface quadrature rule
    */

    SecondOrderSolver( const data_type& data,
                          FESpace<Mesh, EpetraMap>&   FESpace,
                          boost::shared_ptr<Epetra_Comm>&     comm,
                          UInt       offset=0);

     //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */

    void updateSystem();

  //updateSystem with time depend term 
 
  void updateSystem( double coefficient,
		     const vector_type& sourceVec
		     );
  
   void buildSystem();

    //! Solve the non-linear system

    void iterate( vector_type& sol );
    void iterate( bchandler_raw_type& bch );
 
     //! getters
    EpetraMap   const& getMap()       const { return M_localMap; }
    Displayer   const& getDisplayer() const { return M_Displayer; }
    matrix_ptrtype const matrMassStiff() const { return M_massStiff; }
   
    matrix_type&  matrMass()         { return *M_mass; }
    matrix_type&  matrDamping()      { return *M_damping; }
    matrix_type&  matrLinearStiff()  { return *M_linearStiff; }

    //! BCHandler getter and setter
    FESpace<Mesh, EpetraMap>& uFESpace(){return M_FESpace;}
    BCHandler const & BChandler() const {return M_BCh;}
    //! residual getter
    vector_type& residual()             {return *M_residual_d;}
    // end of getters

    //!setters
    void setBC(BCHandler& BCd)   {M_BCh = &BCd;}
    void setSourceTerm( source_type const& __s ) { M_source = __s; }
    void resetPrec() {M_resetPrec = true;}
   
    //! evaluates residual for newton interations
    void evalResidual( vector_type &res, const vector_type& sol, int iter);

    source_type const& sourceTerm() const { return M_source; }

    vector_type& u()        { return M_u; }
 
    vector_ptrtype& rhsWithoutBC() { return M_rhsNoBC; }

    void setUp( const GetPot& dataFile );

    void initialize( const Function& u0);
 
  
   // Matteo
    void solver0( vector_type& sol0,  bchandler_raw_type& bch, vector_type  rhs0);


    void postProcess();


  //Matteo: updateRHS
  virtual void updateRHS(vector_type const& rhs)
  {
    *M_rhsNoBC = rhs;
    M_rhsNoBC->GlobalAssemble();
  }

   //Epetra_Map const& getRepeatedEpetraMap() const { return *M_localMap.getRepeatedEpetra_Map(); }

    boost::shared_ptr<Epetra_Comm> const& comm()         const {return M_Displayer.comm();}

    const UInt& offset() const { return M_offset; }

private:

    const data_type&               M_data;

    FESpace<Mesh, EpetraMap>&      M_FESpace;

    Displayer                      M_Displayer;

    int                            M_me;

    BCHandler*                     M_BCh;

    EpetraMap                      M_localMap;


    //! Matrix M: mass
    matrix_ptrtype                    M_mass;

    //! Matrix Kl: stiffness linear
    matrix_ptrtype                    M_linearStiff;

    //! Matrix Knl: stiffness non-linear
    matrix_ptrtype                    M_stiff;

    //! Matrix C: mass + linear stiffness
    matrix_ptrtype                    M_massStiff;

  //! Matrix D: alpha* mass + beta * linear stiffness
    matrix_ptrtype                    M_damping;
  
    //! Elementary matrices and vectors
    ElemMat                        M_elmatK; // stiffnes
    ElemMat                        M_elmatM; // mass
    ElemMat                        M_elmatC; // mass + stiffness
    ElemMat                        M_elmatD; // mass + stiffness

    //! unknowns vector
    vector_type                    M_u;
   
    //! right  hand  side displacement
    vector_ptrtype                    M_rhs;

      //! right  hand  side
    vector_ptrtype                    M_rhsNoBC;

    //! right  hand  side
    boost::shared_ptr<vector_type>                    M_f;

    //! residual
    boost::shared_ptr<vector_type>                    M_residual_d;

    //! files for lists of iterations and residuals per timestep
    std::ofstream                  M_out_iter;
    std::ofstream                  M_out_res;

    source_type                    M_source;

    //! data for solving tangent problem with aztec
    boost::shared_ptr<solver_type>                    M_linearSolver;

    bool                         M_reusePrec;
    int                            M_maxIterForReuse;
    bool                         M_resetPrec;

    int                            M_maxIterSolver;
    int                            M_count;

    UInt                           M_offset;
    Real                          M_rescaleFactor;

    //
    //! methods
    //

    matrix_ptrtype                  M_matrFull;

    void applyBoundaryConditions(matrix_type &matrix,
                                 vector_type &rhs,
                                 bchandler_raw_type& BCh,
                                 UInt         offset=0);


    UInt dim() const { return M_FESpace.dim(); }

};


//
//                                         IMPLEMENTATION
//
template <typename Mesh, typename SolverType>
SecondOrderSolver<Mesh, SolverType>::
SecondOrderSolver( const data_type&          data,
                      FESpace<Mesh, EpetraMap>& FESpace,
                      BCHandler&                BCh,
                      boost::shared_ptr<Epetra_Comm>&              comm,
                     UInt                      offset
                      ) :
    M_data                       ( data ),
    M_FESpace                    ( FESpace ),
    M_BCh                        ( &BCh ),
    M_Displayer                  ( comm ),
    M_me                         ( comm->MyPID() ),
    M_linearSolver               ( new SolverType( comm ) ),
    M_localMap                   ( M_FESpace.map() ),
    M_mass                       ( new matrix_type(M_localMap) ),
    M_linearStiff                ( new matrix_type(M_localMap) ),
    M_stiff                      ( new matrix_type(M_localMap) ),
    M_massStiff                  ( new matrix_type(M_localMap) ),
    M_damping                    ( new matrix_type(M_localMap) ),
    M_elmatK                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatM                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatC                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatD                     ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_u                          ( M_localMap ),
    M_rhs                        ( new vector_type(M_localMap) ),
    M_rhsNoBC                    ( new vector_type(M_localMap) ),
    M_f                          ( new vector_type(M_localMap)),
    M_residual_d                 (  new vector_type(M_localMap)),
    M_out_iter                   ( "out_iter_problem" ),
    M_out_res                    ( "out_res_problem" ),
    M_reusePrec                  ( true ),
    M_maxIterForReuse            ( -1 ),
    M_resetPrec                  ( true ),
    M_maxIterSolver              ( -1 ),
    M_offset                     (offset),

    M_matrFull                   ( )
{
    //    M_BCh->setOffset(M_offset);
}

template <typename Mesh, typename SolverType>
SecondOrderSolver<Mesh, SolverType>::
SecondOrderSolver( const data_type& data,
                      FESpace<Mesh, EpetraMap>&   FESpace,
                     boost::shared_ptr<Epetra_Comm>&     comm,
                      UInt             /*offset*/
                      ) :
    M_data                        ( data ),
    M_FESpace                 ( FESpace ),
    M_Displayer               ( comm ),
    M_me                          ( comm->MyPID() ),
    M_localMap                ( M_FESpace.map() ),
    M_mass                       ( new matrix_type(M_localMap) ),
    M_linearStiff               ( new matrix_type(M_localMap) ),
    M_stiff                         ( new matrix_type(M_localMap) ),
    M_massStiff               ( new matrix_type(M_localMap) ),
    M_damping                ( new matrix_type(M_localMap) ),
    M_elmatK                   ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatM                  ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatC                   ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_elmatD                   ( M_FESpace.fe().nbNode, nDimensions, nDimensions ),
    M_u                              ( M_localMap),
    M_rhs                           ( new vector_type(M_localMap) ),
    M_rhsNoBC                 ( new vector_type(M_localMap) ),
    M_f                               ( new vector_type(M_localMap)),
    M_residual_d              ( new vector_type(M_localMap)),
    M_out_iter                   ( "out_iter_problem" ),
    M_out_res                    ( "out_res_problem" ),
    M_linearSolver               ( new SolverType( comm ) ),
    //M_prec                       ( ),//new prec_raw_type() ),
    M_reusePrec                  ( true ),
    M_maxIterForReuse            ( -1 ),
    M_resetPrec                  ( true ),
    M_maxIterSolver              ( -1 ),
    M_count                      ( 0 ),
    M_offset                     ( 0 ),
    M_rescaleFactor            (1.),
    M_matrFull                   ( )
{
    //    M_BCh->setOffset(0);
}


template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
setUp( const GetPot& dataFile )
{
    M_Displayer.leaderPrint("\n P-   unknowns: ",  M_FESpace.dof().numTotalDof() );
    M_Displayer.leaderPrint(" P-  Computing mass and linear strain matrices ... \n");

    M_linearSolver->setDataFromGetPot( dataFile, "problem/solver" );

    M_reusePrec     = dataFile( "problem/prec/reuse", true);
    M_maxIterSolver = dataFile( "problem/solver/max_iter", -1);
    M_maxIterForReuse = dataFile( "problem/solver/max_iter_reuse", M_maxIterSolver*8/10);

    M_linearSolver->setUpPrec(dataFile, "problem/prec");

}



template <typename Mesh, typename SolverType>
void
SecondOrderSolver<Mesh, SolverType>::
buildSystem()
{
    UInt totalDof = M_FESpace.dof().numTotalDof();

    M_Displayer.leaderPrint( "P-  Building the system             ... ");

    Chrono chrono;
    chrono.start();
 
    // Number of displacement components
    UInt nc = nDimensions;

     // Elementary computation and matrix assembling
    // Loop on elements

    for ( UInt iVol = 1; iVol <= M_FESpace.mesh()->numVolumes(); iVol++ )
      {
        M_FESpace.fe().updateFirstDeriv( M_FESpace.mesh()->element( iVol ) );

	int marker    = M_FESpace.mesh()->volumeList( iVol ).marker();
	
	//	chrono.start();
        M_elmatK.zero();
        M_elmatM.zero();

	// building stiffness matrix
 	stiff_strain( 2.0, M_elmatK, M_FESpace.fe() );

        // building mass matrix
	mass( 1., M_elmatM, M_FESpace.fe(), 0, 0);
      
        assembleMatrix( *M_linearStiff,
                        M_elmatK,
                        M_FESpace.fe(),
                        M_FESpace.fe(),
                        M_FESpace.dof(),
                        M_FESpace.dof(),
                        0, 0, 0, 0);

	assembleMatrix( *M_mass,
			M_elmatM,
			M_FESpace.fe(),
			M_FESpace.fe(),
			M_FESpace.dof(),
			M_FESpace.dof(),
			0, 0, 0, 0);
      }

    if(M_data.isDamping())
      {
	M_Displayer.leaderPrint( "P-  Building the system   Damping matrix          ... ");
	
	for ( UInt iVol = 1; iVol <= M_FESpace.mesh()->numVolumes(); iVol++ )
	  {
	    M_FESpace.fe().updateFirstDeriv( M_FESpace.mesh()->element( iVol ) );
	    
	    int marker    = M_FESpace.mesh()->volumeList( iVol ).marker();
	    
	    double alpha = M_data.alpha(marker);
	    double beta  = M_data.beta(marker);
	    //building damping matrix
	    
	    M_elmatD.zero();
 
	    stiff_strain( 2.0 * beta, M_elmatD, M_FESpace.fe() );
	    mass( alpha, M_elmatD, M_FESpace.fe(), 0, 0);
	    
	    assembleMatrix( *M_damping,
			    M_elmatD,
			    M_FESpace.fe(),
			    M_FESpace.fe(),
			    M_FESpace.dof(),
			    M_FESpace.dof(),
			    0, 0, 0, 0);
	  }
      }
    
    //M_comm->Barrier();

    //Global Assemble :
    M_linearStiff->GlobalAssemble();
    M_mass->GlobalAssemble();
    M_damping->GlobalAssemble();

    M_Displayer.leaderPrintMax( " done in ", chrono.diff() );
}

template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
solver0( vector_type& sol0,  bchandler_raw_type& bch, vector_type  rhs0 )
{
  vector_type rhsFull(rhs0);
  
  prec_type prec;
  
  std::string precType ="Ifpack" ;
  
  prec.reset( PRECFactory::instance().createObject(precType) );
  
  prec->buildPreconditioner(M_mass);
  
  Real condest = prec->Condest();
  
  M_linearSolver->setPreconditioner(prec);
  
  M_linearSolver->setMatrix(*M_mass);
  
  int numIter =  M_linearSolver->solve(sol0, rhs0);

}


template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
initialize( const Function& u0)
{
  M_FESpace.interpolate(u0, M_u, 0.0);
}

template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
updateSystem(const      double  coefficient,
	     const vector_type& sourceVec   )
{
  M_Displayer.leaderPrint("  P-  Updating mass term on right hand side... "); 
  
  Chrono chrono;
  chrono.start();
  
  updateRHS(sourceVec);
  
  *M_stiff *=0;
  *M_massStiff *=0;
   
  *M_massStiff += *M_linearStiff;
  
  if (coefficient!= 0. )
    *M_massStiff += *M_mass * coefficient;
   
  M_massStiff->GlobalAssemble();
  chrono.stop();
 
  M_Displayer.leaderPrintMax("done in ", chrono.diff()); 
}


template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
iterate( bchandler_raw_type& bch )
{
  Chrono chrono;
  
  // matrix and vector assembling communication
  M_Displayer.leaderPrint("  P-  Solving the system ... \n");
  
  chrono.start();
  
    matrix_ptrtype matrFull( new matrix_type( M_localMap, M_massStiff->getMeanNumEntries()));
   *matrFull += *M_massStiff;
  
   if (M_data.isDamping())
      *matrFull += *M_damping;
    
    M_rhsNoBC->GlobalAssemble();
    

    vector_type rhsFull (*M_rhsNoBC);
    
    // boundary conditions update
    M_Displayer.leaderPrint("  P-  Applying boundary conditions ...         ");

    chrono.start();
    applyBoundaryConditions( *matrFull, rhsFull, bch);
    
    chrono.stop();
    
    M_Displayer.leaderPrintMax("done in " , chrono.diff());
    
    M_linearSolver->precReset();
    M_resetPrec = false;
    
    // solving the system
    M_linearSolver->setMatrix(*matrFull);
    int numIter = M_linearSolver->solveSystem( rhsFull, M_u, matrFull);

    numIter = abs(numIter);

    if(numIter >= M_maxIterForReuse || numIter >= M_maxIterSolver)
    {
        resetPrec();
    }

    *M_residual_d =  *M_massStiff * (M_u);
    *M_residual_d -= *M_rhsNoBC;

} // iterate()


template <typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
evalResidual( vector_type &res, const vector_type& sol, int /*iter*/)
{
    M_Displayer.leaderPrint( "    P-    Computing residual... \n");
    Chrono chrono;
    chrono.start();

    // Matrices initialization
    M_stiff = M_massStiff;

    M_Displayer.leaderPrint("updating the boundary conditions ... ");
    if ( !M_BCh->bdUpdateDone() )
        M_BCh->bdUpdate( M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    bcManageMatrix( M_stiff, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), 1.0 );

    *M_rhs = *M_rhsNoBC;

    bcManageVector( *M_rhs, *M_FESpace.mesh(), M_FESpace.dof(), *M_BCh, M_FESpace.feBd(), M_data.getTime(), 1.0 );

    res  = M_stiff * sol;

    chrono.stop();
    M_Displayer.leaderPrintMax("done in ", chrono.diff() );
}

template<typename Mesh, typename SolverType>
void SecondOrderSolver<Mesh, SolverType>::
applyBoundaryConditions(matrix_type&        matrix,
                        vector_type&        rhs,
                        bchandler_raw_type& BCh,
                        UInt                offset)
{

    // BC manage for the velocity
    if(offset)
        BCh.setOffset(offset);
    if ( !BCh.bdUpdateDone() )
        BCh.bdUpdate( *M_FESpace.mesh(), M_FESpace.feBd(), M_FESpace.dof() );

    // vector_type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_type rhsFull(rhs, Unique);  // bcManages now manages the also repeated parts

    bcManage( matrix, rhsFull, *M_FESpace.mesh(), M_FESpace.dof(), BCh, M_FESpace.feBd(), 1.,
              M_data.getTime() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

} // applyBoundaryCondition


}
#endif
