/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s):  A. Fumagalli  <alessio.fumagalli@mail.polimi.it>
       Date: 2010-03-24

  Copyright (C) 2010 EPFL, Politecnico di Milano

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file darcy.cpp
   \author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   \date
 */

/*!
  Simple 3D non-linear Darcy test with Dirichlet, Neumann and Robin Boundary condition.
  <br>
  Solve the problem in dual-mixed form
  \f[
  \left\{
  \begin{array}{l l l }
  \Lambda^{-1}(p) \sigma + \nabla p = 0 & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  \nabla \cdot \sigma - f = 0        & \mathrm{in} & \Omega\,,  \vspace{0.2cm} \\
  p = g_D                            & \mathrm{on} & \Gamma_D\,,\vspace{0.2cm} \\
  \sigma \cdot n + h p = g_R         & \mathrm{on} & \Gamma_R\,, \vspace{0.2cm} \\
  \sigma \cdot n = g_n               & \mathrm{on} & \Gamma_N\,.
  \end{array}
  \right\.
  \f]
where \f$ \Omega \f$ is the unit cube with
  \f[
  \begin{array}{l}
  \Gamma_R = \left\{ y = 1 \right\}\,, \vspace{0.2cm} \\
  \Gamma_N = \left\{ x = 1 \right\} \cup \left\{ x = 0 \right\}\,, \vspace{0.2cm} \\
  \Gamma_D = \partial [0,1]^3 \setminus \left( \Gamma_R \cup \Gamma_D \right)\,,
  \end{array}
  \f]
and the data are
  \f[
  \left\{
  \begin{array}{l}
  f(x,y,z) =  -(2(x^2y^2+6x+5z)(2xy^2+6)(2xy^2+6)+2((x^2y^2+6x+5z)^2+2)y^2+8x^2+2x^2)\,, \vspace{0.2cm} \ \
  g_D(x,y,z) = x^2y^2 + 6x + 5z\,, \vspace{0.2cm} \\
  h(x,y,z) = 1\,, \vspace{0.2cm} \\
  g_R(x,y,z) =  -(2xy^2+6+2x^2y) +  x^2y^2 + 6x + 5z\,, \vspace{0.2cm} \\
  g_N(x,y,z) = \pm (((x^2y^2+6x+5z)^2+2)(2xy^2+6)+2x^2y )\,, \vspace{0.2cm} \\
  K(x,y,z;p) = \left\[
  \begin{array}{c c c}
  p^2 + 2 & 1 & 0 \\
  1       & 1 & 0 \\
  0       & 0 & 2
  \end{array}
  \right\]
  \end{array}
  \right\.
  \f]
The analytical solutions are
  \f[
  p(x,y,z) = x^2y^2 + 6x + 5z\,, \vspace{0.2cm} \\
  \sigma(x,y,z) = \left(
  \begin{array}{l}
  ((x^2y^2+6x+5z)^2+2)(2xy^2+6)+2x^2y \\
  2xy^2+6+2x^2y \\
  10
  \end{array}
  \right)\,.
  \f]
*/

// ===================================================
//! Includes
// ===================================================

#include "darcy.hpp"

// ===================================================
//! Namespaces & define
// ===================================================

using namespace LifeV;

enum BCNAME
{
    BACK   = 1,
    FRONT  = 2,
    LEFT   = 3,
    RIGHT  = 4,
    BOTTOM = 5,
    TOP    = 6
};

// ===================================================
//! User functions
// ===================================================

// Analytical solution
Real analyticalSolution( const Real& /*t*/,
                         const Real& x,
                         const Real& y,
                         const Real& z,
                         const ID& /*ic*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Analytical flux
Real analyticalFlux( const Real& /*t*/,
                     const Real& x,
                     const Real& y,
                     const Real& z,
                     const ID& icomp)
{
    switch(icomp)
    {
    case 1: // \frac{\partial }{\partial x}
        return -1.*(((x*x*y*y + 6.*x + 5.*z)*(x*x*y*y + 6.*x + 5.*z) + 1) * (2.*x*y*y + 6.) + 2.*x*x*y);
    case 2: // \frac{\partial }{\partial y}
        return -1.*(2.*x*y*y + 6. + 2.*x*x*y);
    case 3: // \frac{\partial }{\partial z}
        return -10.;
    default:
        return 0.;
    }
}

// Non-linear inverse of permeability matrix.
/* In this case the permeability matrix is
K = [p^2+2 1   0
     1     1   0
     0     0   2]
*/
Matrix nonLinearInversePermeability( const Real& primalOld,
                                     const Real& /*t*/,
                                     const Real& x,
                                     const Real& y,
                                     const Real& z )
{
    Matrix inversePermeabilityMatrix( static_cast<UInt>(3), static_cast<UInt>(3) );

    // First row
    Real Entry00 = 1./(primalOld*primalOld+1);
    Real Entry01 = -1./(primalOld*primalOld+1);
    Real Entry02 = 0.;

    // Second row
    Real Entry11 = (primalOld*primalOld+2)/(primalOld*primalOld+1);
    Real Entry12 = 0.;

    // Third row
    Real Entry22 = 1./2.;

    // Fill in of the inversePermeabilityMatrix
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(0) ) = Entry00;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(1) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(0), static_cast<UInt>(2) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(0) ) = Entry01;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(1) ) = Entry11;
    inversePermeabilityMatrix( static_cast<UInt>(1), static_cast<UInt>(2) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(0) ) = Entry02;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(1) ) = Entry12;
    inversePermeabilityMatrix( static_cast<UInt>(2), static_cast<UInt>(2) ) = Entry22;

    return inversePermeabilityMatrix;
}

// Zero iteration primal variable
Real zeroItarationPrimal( const Real& /*t*/,
                          const Real& x,
                          const Real& y,
                          const Real& /*z*/,
                          const ID& /*ic*/)
{
    return 0;
}

// Boundary condition of Dirichlet
Real dirichlet( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&   /*icomp*/)
{
    return x*x*y*y + 6.*x + 5.*z;
}

// Boundary condition of Neumann
Real neumann1( const Real& /*t*/,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
	switch(icomp){
  		case 1:   //! Dx
            return -1.*( ((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*(2*x*y*y+6)+2*x*x*y );
		break;
		case 2:   //! Dy
            return 0.;
	 	break;
  		case 3:   //! Dz
    		return 0.;
    	break;
  	}
	return 0.;
}

// Boundary condition of Neumann
Real neumann2( const Real& /*t*/,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   icomp)
{
	switch(icomp){
  		case 1:   //! Dx
            return ((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*(2*x*y*y+6)+2*x*x*y;
		break;
		case 2:   //! Dy
            return 0.;
	 	break;
  		case 3:   //! Dz
    		return 0.;
    	break;
  	}
	return 0.;
}

// Boundary condition of Robin
Real mixte( const Real& /*t*/,
            const Real& x,
            const Real& y,
            const Real& z,
            const ID&   /*icomp*/)
{
    return -1.*(2*x*y*y+6+2*x*x*y) +  x*x*y*y + 6.*x + 5.*z;
}

// Source term
Real source_in( const Real& /*t*/,
                const Real& x,
                const Real& y,
                const Real& z,
                const ID&  /*icomp*/)
{
    return -1.*(2*(x*x*y*y+6*x+5*z)*(2*x*y*y+6)*(2*x*y*y+6)+2*((x*x*y*y+6*x+5*z)*(x*x*y*y+6*x+5*z)+2)*y*y+8*x*y+2*x*x);
}

// Standard functions
Real UOne( const Real& /* t */,
           const Real& /* x */,
           const Real& /* y */,
           const Real& /* z */,
           const ID&   /* icomp */)
{
    return 1.;
}

Real UZero( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID&   /* icomp */)
{
    return 0.;
}

// ===================================================
//! Private Members
// ===================================================

struct darcy::Private
{
    Private() {}

    // Policy for functions
    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )>
                                     fct_type;
    // Policy for matrices
    typedef boost::function<Matrix ( const Real&, const Real&, const Real&,
                                     const Real&, const Real& )>
                                     nonLinearMatrix_type;

    std::string    data_file_name;
    std::string    discretization_section;

    boost::shared_ptr<Epetra_Comm>   comm;

    // Function Types

    fct_type getUOne()
    {
    	fct_type f;
    	f = boost::bind( &UOne, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getUZero()
    {
    	fct_type f;
    	f = boost::bind( &UZero, _1, _2, _3, _4, _5 );
    	return f;
    }

    fct_type getAnalyticalSolution()
    {
        fct_type f;
        f = boost::bind( &analyticalSolution, _1, _2, _3, _4, _5 );
        return f;
    }

    fct_type getAnalyticalFlux()
    {
        fct_type f;
        f = boost::bind( &analyticalFlux, _1, _2, _3, _4, _5 );
        return f;
    }

    nonLinearMatrix_type getNonLinearInversePermeability()
    {
        nonLinearMatrix_type m;
        m = boost::bind( &nonLinearInversePermeability, _1, _2, _3, _4, _5 );
        return m;
    }

};

// ===================================================
//! Constructors
// ===================================================

darcy::darcy( int argc,
              char** argv,
              LifeV::AboutData const& /*ad*/,
              LifeV::po::options_description const& /*od*/ ): Members( new Private )
{
    GetPot command_line(argc, argv);
    const string data_file_name = command_line.follow("data", 2, "-f", "--file");
    GetPot dataFile( data_file_name );

    Members->data_file_name = data_file_name;
    Members->discretization_section = "darcy";

	#ifdef EPETRA_MPI
		std::cout << "Epetra Initialization" << std::endl;
		Members->comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );
        int ntasks;
        MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	#else
		Members->comm.reset( new Epetra_SerialComm() );
	#endif
}

// ===================================================
//! Methods
// ===================================================

Real
darcy::run()
{
    typedef RegionMesh3D<LinearTetra>                       RegionMesh;
    typedef SolverTrilinos                                  solver_type;
    typedef DarcySolverNonLinear< RegionMesh, solver_type > ds;
	typedef ds::vector_type                                 vector_type;
    typedef boost::shared_ptr<vector_type>                  vector_ptrtype;

    Chrono chronoTotal;
    Chrono chronoReadAndPartitionMesh;
    Chrono chronoBoundaryCondition;
    Chrono chronoFiniteElementSpace;
    Chrono chronoProblem;
    Chrono chronoProcess;
    Chrono chronoError;

    // Start chronoTotal for measure the total time for the computation
    chronoTotal.start();

    // Reading from data file
    GetPot dataFile( Members->data_file_name.c_str() );

    // Create the leader process, i.e. the process with MyPID equal to zero
    bool isLeader = ( Members->comm->MyPID() == 0 );

    //
    // The Darcy Solver
    //

    if ( isLeader )
        std::cout << "The Darcy non-linear solver" << std::endl << std::flush;

    // Start chronoReadAndPartitionMesh for measure the total time for the creation of the local meshes
    chronoReadAndPartitionMesh.start();

    // Create the data file
	DataDarcy<RegionMesh> dataDarcy;

    // Set up the data
    dataDarcy.setup( dataFile );

    // Create the mesh file handler
    DataMesh dataMesh;

    // Set up the mesh file
    dataMesh.setup( dataFile,  Members->discretization_section + "/space_discretization");

    // Create the mesh
    boost::shared_ptr<RegionMesh> fullMeshPtr( new RegionMesh );

    // Set up the mesh
    readMesh( *fullMeshPtr, dataMesh );

    // Partition the mesh using ParMetis
    partitionMesh< RegionMesh >  meshPart( fullMeshPtr, Members->comm );

    // Stop chronoReadAndPartitionMesh
    chronoReadAndPartitionMesh.stop();

    // The leader process print chronoReadAndPartitionMesh
    if( isLeader )
        std::cout << "Time for read and partition the mesh " <<
                     chronoReadAndPartitionMesh.diff() << std::endl << std::flush;

    // Create the boundary conditions

    // Start chronoBoundaryCondition for measure the total time for create the boundary conditions
    chronoBoundaryCondition.start();

    BCFunctionBase dirichletBDfun, neumannBDfun1, neumannBDfun2;
   	BCFunctionMixte mixteBDfun;

    dirichletBDfun.setFunction( dirichlet );
	neumannBDfun1.setFunction( neumann1 );
	neumannBDfun2.setFunction( neumann2 );
	// dp/dn = first_parameter + second_parameter * p
	mixteBDfun.setFunctions_Mixte( mixte, Members->getUOne() );

	BCHandler bcDarcy( 6 );

    bcDarcy.addBC( "Top",     TOP,     Natural,    Full,    neumannBDfun1, 1 );
    bcDarcy.addBC( "Bottom",  BOTTOM,  Mixte,      Scalar,  mixteBDfun      );
    //bcDarcy.addBC(   "Top",    TOP,    Essential,  Scalar,  dirichletBDfun  );
    //bcDarcy.addBC("Bottom", BOTTOM,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy.addBC(  "Left",   LEFT,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy.addBC( "Right",  RIGHT,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy.addBC( "Front",  FRONT,    Essential,  Scalar,  dirichletBDfun  );
    //bcDarcy.addBC(  "Back",   BACK,    Essential,  Scalar,  dirichletBDfun  );
    bcDarcy.addBC( "Back",    BACK,    Natural,    Full,    neumannBDfun2, 1 );

    // Stop chronoBoundaryCondition
    chronoBoundaryCondition.stop();

    // The leader process print chronoBoundaryCondition
    if ( isLeader )
    {
        std::cout << "Time for create the boundary conditions handler " <<
                     chronoBoundaryCondition.diff() << std::endl << std::flush;

    }

    // Create the solution spaces

    // Start chronoFiniteElementSpace for measure the total time for create the finite element spaces
    chronoFiniteElementSpace.start();

    // Primal solution parameters
    const RefFE*    refFE_primal ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR_primal    ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr_primal  ( static_cast<QuadRule*>(NULL) );

    refFE_primal = &feTetraP0;
    qR_primal    = &quadRuleTetra15pt;
    bdQr_primal  = &quadRuleTria4pt;

	// Dual solution parameters
	const RefFE*    refFE_dual   ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR_dual      ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr_dual    ( static_cast<QuadRule*>(NULL) );

    refFE_dual = &feTetraRT0;
    qR_dual    = &quadRuleTetra15pt;
    bdQr_dual  = &quadRuleTria4pt;

    // Interpolate of dual solution parameters
	const RefFE*    refFE_dualInterpolate ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR_dualInterpolate    ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr_dualInterpolate  ( static_cast<QuadRule*>(NULL) );

    refFE_dualInterpolate = &feTetraP0;
    qR_dualInterpolate    = &quadRuleTetra15pt;
    bdQr_dualInterpolate  = &quadRuleTria4pt;

    // Hybrid solution parameters
    // hybrid
    const RefFE*    refFE_hybrid ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR_hybrid    ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr_hybrid  ( static_cast<QuadRule*>(NULL) );

    refFE_hybrid = &feTetraRT0Hyb;
    qR_hybrid    = &quadRuleTetra15pt;
    bdQr_hybrid  = &quadRuleTria4pt;

    // dual dot outward unit normal
    const RefFE*    refFE_VdotN ( static_cast<RefFE*>(NULL) );
    const QuadRule* qR_VdotN    ( static_cast<QuadRule*>(NULL) );
    const QuadRule* bdQr_VdotN  ( static_cast<QuadRule*>(NULL) );

    refFE_VdotN = &feTetraRT0VdotNHyb;
    qR_VdotN    = &quadRuleTetra15pt;
    bdQr_VdotN  = &quadRuleTria4pt;


    // Finite element space of the primal variable
    FESpace< RegionMesh, EpetraMap > p_FESpace( meshPart,
												*refFE_primal,
                                                *qR_primal,
                                                *bdQr_primal,
                                                1,
                                                Members->comm );

    // Finite element space of the dual variable
    FESpace< RegionMesh, EpetraMap > u_FESpace( meshPart,
                                                *refFE_dual,
                                                *qR_dual,
                                                *bdQr_dual,
                                                1,
                                                Members->comm );

    // Finite element space of the interpolation of dual variable
    FESpace< RegionMesh, EpetraMap > uInterpolate_FESpace( meshPart,
                                                           *refFE_dualInterpolate,
                                                           *qR_dualInterpolate,
                                                           *bdQr_dualInterpolate,
                                                           3,
                                                           Members->comm );

    // Vector for the interpolated dual solution
    vector_ptrtype dualInterpolated( new vector_type ( uInterpolate_FESpace.map(), Repeated ) );

    // Finite element space of the hybrid variable
    FESpace< RegionMesh, EpetraMap > hybrid_FESpace( meshPart,
                                                    *refFE_hybrid,
                                                    *qR_hybrid,
                                                    *bdQr_hybrid,
                                                    1,
                                                    Members->comm );

    // Finite element space of the  outward unit normal variable
    FESpace< RegionMesh, EpetraMap > VdotN_FESpace( meshPart,
                                                    *refFE_VdotN,
                                                    *qR_VdotN,
                                                    *bdQr_VdotN,
                                                    1,
                                                    Members->comm );

    // Stop chronoFiniteElementSpace
    chronoFiniteElementSpace.stop();

    // The leader process print chronoFiniteElementSpace
    if ( isLeader )
        std::cout << "Time for create the finite element spaces " <<
                      chronoFiniteElementSpace.diff() << std::endl << std::flush;

    // Start chronoProblem for measure the total time for create the problem
    chronoProblem.start();

    // Instantiation of the DarcySolverNonLinear class
    ds darcySolver ( dataDarcy,
                     p_FESpace,
                     u_FESpace,
                     hybrid_FESpace,
                     VdotN_FESpace,
                     Members->comm );

    // Stop chronoProblem
    chronoProblem.stop();

    // The leader process print chronoProblem
    darcySolver.getDisplayer().leaderPrint( "Time for create the problem ",
                                            chronoProblem.diff(), "\n" );

    // Process the problem

    // Start chronoProcess for measure the total time for the simulation
    chronoProcess.start();

    // Setup phase for the linear solver
    darcySolver.setup();

    // Set the initial primal variable
    darcySolver.setZeroItarationPrimal( zeroItarationPrimal );

	// Set the source term
    darcySolver.setSourceTerm( source_in );

    // Set the non-linear inverse of the permeability
    darcySolver.setNonLinearInversePermeability( Members->getNonLinearInversePermeability() );

    // Set the boudary conditions
    darcySolver.setBC( bcDarcy );

    // Set the exporter for the solution
    boost::shared_ptr< Exporter< RegionMesh > > exporter;

    // Shared pointer used in the exporter for the primal solution
    vector_ptrtype primalExporter;

    // Shared pointer used in the exporter for the dual solution
    vector_ptrtype dualExporter;

    // Type of the exporter
    std::string const exporterType =  dataFile( "exporter/type", "ensight");

    // Choose the exporter
#ifdef HAVE_HDF5
    if ( exporterType.compare("hdf5") == 0 )
    {
        exporter.reset( new Hdf5exporter< RegionMesh > ( dataFile, "PressureVelocity" ) );

        // Set directory where to save the solution
        exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

        exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
    }
    else
#endif
    {
        if ( exporterType.compare("none") == 0 )
        {
            exporter.reset( new NoExport< RegionMesh > ( dataFile, "PressureVelocity" ) );

            // Set directory where to save the solution
            exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
        }
        else
        {
            exporter.reset( new Ensight< RegionMesh > ( dataFile, "PressureVelocity" ) );

            // Set directory where to save the solution
            exporter->setDirectory( dataFile( "exporter/folder", "./" ) );

            exporter->setMeshProcId( meshPart.mesh(), Members->comm->MyPID() );
        }
    }

    // Set the exporter primal pointer
    primalExporter.reset( new vector_type ( *darcySolver.primalSolution(),
                                            exporter->mapType() ) );

    // Add the primal variable to the exporter
    exporter->addVariable( ExporterData::Scalar,
                           "Pressure",
                           primalExporter,
                           static_cast<UInt>( 0 ),
                           static_cast<UInt>( p_FESpace.dof().numTotalDof() ),
                           static_cast<UInt>( 0 ),
                           ExporterData::Cell );

    // Set the exporter dual pointer
    dualExporter.reset( new vector_type ( *dualInterpolated, exporter->mapType() ) );

    // Add the variable to the exporter
    exporter->addVariable( ExporterData::Vector,
                           "Velocity",
                           dualExporter,
                           static_cast<UInt>( 0 ),
                           static_cast<UInt>( uInterpolate_FESpace.dof().numTotalDof() ),
                           static_cast<UInt>( 0 ),
                           ExporterData::Cell );

    // Compute the total number of unknowns
    EpetraMap fullDarcyMap( darcySolver.getMap() );
    if ( isLeader )
    	std::cout << "Number of unknowns : " << hybrid_FESpace.map().getMap(Unique)->NumGlobalElements() << std::endl << std::flush;

    // Solve the problem

    // Start the fixed point simulation
    darcySolver.fixedPointScheme();

    // Save the solution

    // Copy the primal solution to the exporter
    *primalExporter = *darcySolver.primalSolution();

    // Interpolate the dual vector field spammed as Raviart-Thomas into a P0 vector field
    *dualInterpolated = uInterpolate_FESpace.FeToFeInterpolate( u_FESpace,
                                                                *darcySolver.dualSolution() );

    // Copy the dual interpolated solution to the exporter
    *dualExporter = *dualInterpolated;

    // Save the solution into the exporter
    exporter->postProcess( static_cast<Real>(0) );

    // Stop chronoProcess
    chronoProcess.stop();

    // The leader process print chronoProcess
    darcySolver.getDisplayer().leaderPrint( "Time for process ",
                                            chronoProcess.diff(), "\n" );

    // Start chronoError for measure the total time for computing the errors
    chronoError.start();

    // The leader process print chronoProcess
    darcySolver.getDisplayer().leaderPrint( "Time for process ",
                                            chronoProcess.diff(), "\n" );

    // Compute the errors

    // Start chronoError for measure the total time for computing the errors.
    chronoError.start();

    // Compute the error L2 norms
    Real primalL2Norm(0), exactPrimalL2Norm(0), primalL2Error(0), primalL2RelativeError(0);
    Real dualL2Norm(0), exactDualL2Norm(0), dualL2Error(0), dualL2RelativeError(0);

    // Norms and errors for the pressure
    darcySolver.getDisplayer().leaderPrint( "\nPRESSURE ERROR\n" );

    // Compute the L2 norm for the primal solution
    primalL2Norm = p_FESpace.L2Norm( *darcySolver.primalSolution() );

    // Display the L2 norm for the primal solution
    darcySolver.getDisplayer().leaderPrint( " L2 norm of primal unknown:            ",
                                            primalL2Norm, "\n" );

    // Compute the L2 norm for the analytical primal
    exactPrimalL2Norm = p_FESpace.L2NormFunction( Members->getAnalyticalSolution(),
                                                  dataDarcy.dataTime()->getEndTime() );

    // Display the L2 norm for the analytical primal
    darcySolver.getDisplayer().leaderPrint( " L2 norm of primal exact:              ",
                                            exactPrimalL2Norm, "\n" );

    // Compute the L2 error for the primal solution
    primalL2Error = p_FESpace.L2ErrorWeighted( Members->getAnalyticalSolution(),
                                               *darcySolver.primalSolution(),
                                               Members->getUOne(),
                                               dataDarcy.dataTime()->getEndTime() );

    // Display the L2 error for the primal solution
    darcySolver.getDisplayer().leaderPrint( " L2 error of primal unknown:           ",
                                            primalL2Error, "\n" );

    // Compute the L2 realative error for the primal solution
    primalL2RelativeError = primalL2Error / exactPrimalL2Norm;

    // Display the L2 relative error for the primal solution
    darcySolver.getDisplayer().leaderPrint( " L2 relative error of primal unknown:  ",
                                            primalL2RelativeError, "\n" );

    // Norms and errors for the interpolated dual
    darcySolver.getDisplayer().leaderPrint( "\n\nINTERPOLATED DARCY VELOCITY ERROR\n" );

    // Compute the L2 norm for the interpolated dual solution
    dualL2Norm = uInterpolate_FESpace.L2Norm( *dualInterpolated );

    // Display the L2 norm for the interpolated dual solution
    darcySolver.getDisplayer().leaderPrint( " L2 norm of dual unknown:              ",
                                            dualL2Norm, "\n" );

    // Compute the L2 norm for the analytical dual
    exactDualL2Norm = uInterpolate_FESpace.L2NormFunction( Members->getAnalyticalFlux(),
                                                           dataDarcy.dataTime()->getEndTime() );

    // Display the L2 norm for the analytical dual
    darcySolver.getDisplayer().leaderPrint( " L2 norm of dual exact:                ",
                                            exactDualL2Norm, "\n" );

    // Compute the L2 error for the dual solution
    dualL2Error = uInterpolate_FESpace.L2Error( Members->getAnalyticalFlux(),
                                                *dualInterpolated,
                                                dataDarcy.dataTime()->getEndTime(),
                                                NULL );

    // Display the L2 error for the dual solution
    darcySolver.getDisplayer().leaderPrint( " L2 error of dual unknown:             ",
                                            dualL2Error, "\n" );

    // Compute the L2 relative error for the dual solution
    dualL2RelativeError = dualL2Error / exactDualL2Norm;

    // Display the L2 relative error for the dual solution
    darcySolver.getDisplayer().leaderPrint( " L2 relative error of Dual unknown:    ",
                                            dualL2RelativeError, "\n" );

    // Stop chronoError
    chronoError.stop();

    // The leader process print chronoError
    darcySolver.getDisplayer().leaderPrint( "Time for compute errors ",
                                            chronoError.diff(), "\n" );

    // Stop chronoTotal
    chronoTotal.stop();

    // The leader process print chronoTotal
    darcySolver.getDisplayer().leaderPrint( "Total time for the computation ",
                                            chronoTotal.diff(), "\n" );

    // Return the error, needed for the succes/failure of the test
    return primalL2Error;

}
