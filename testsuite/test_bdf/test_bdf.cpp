/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s):  Umberto Villa <uvilla@emory.edu>
 Date: 2010-04-14

 Copyright (C) 2009 EPFL

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
 \file test_bdf.cpp
 \author U. Villa <uvilla@emory.edu>
 \date 2010-04-14
 */

// ===================================================
//! Includes
// ===================================================
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#include <mpi.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfVariableStep.hpp>
#include <life/lifefilters/ensight.hpp>

#include "ud_functions.hpp"
#include "test_bdf.hpp"

#define POSTPROCESS 1

// ===================================================
//! Namespaces & Define
// ===================================================
using namespace LifeV;

#ifdef TWODIM

const int TOP = 3;
const int BOTTOM = 1;
const int LEFT = 4;
const int RIGHT = 2;

typedef RegionMesh2D<LinearTriangle> RegionMesh;

const std::string discretization_section="discretization2D";

#elif defined THREEDIM

const int TOP = 6;
const int BOTTOM = 5;
const int LEFT = 3;
const int RIGHT = 4;
const int FRONT = 2;
const int BACK = 1;

typedef RegionMesh3D<LinearTetra> RegionMesh;

const std::string discretization_section = "discretization3D";

#endif

// ===================================================
//! Private Members
// ===================================================
struct test_bdf::Private {
	Private() {
	}

	typedef boost::function<Real(Real const&, Real const&, Real const&,
			Real const&, ID const&)> fct_type;

	std::string data_file_name;
	boost::shared_ptr<Epetra_Comm> comm;

};

// ===================================================
//! Constructors
// ===================================================
test_bdf::test_bdf(int argc, char** argv, LifeV::AboutData const& /*ad*/,
		LifeV::po::options_description const& /*od*/) :
	Members(new Private) {
	GetPot command_line(argc, argv);
	const string data_file_name =
			command_line.follow("data", 2, "-f", "--file");
	GetPot dataFile(data_file_name);

	Members->data_file_name = data_file_name;

#ifdef EPETRA_MPI
	std::cout << "Epetra Initialization" << std::endl;
	Members->comm.reset(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
	Members->comm.reset(new Epetra_SerialComm() );
#endif
}

// ===================================================
//! Methods
// ===================================================
void test_bdf::run() {

	//Useful typedef
	typedef SolverTrilinos solver_type;
	typedef EpetraVector vector_type;
	typedef boost::shared_ptr<vector_type> vector_ptrtype;

	// Reading from data file
	GetPot dataFile(Members->data_file_name.c_str());
	// Select Pid(0) to write on std output
	bool verbose = (Members->comm->MyPID() == 0);

	if (verbose)
		std::cout << "The BDF Solver" << std::flush;

	//the forcing term
	SourceFct sf;
	// the boundary conditions
	BCFunctionBase g_Ess(AnalyticalSol::u);

#ifdef TWODIM

	BCHandler bc( 4, BCHandler::HINT_BC_ONLY_ESSENTIAL );

	bc.addBC( "Top", TOP, Essential, Full, g_Ess, 1 );
	bc.addBC( "Bottom", BOTTOM, Essential, Full, g_Ess, 1 );
	bc.addBC( "Left", LEFT, Essential, Full, g_Ess, 1 );
	bc.addBC( "Right", RIGHT, Essential, Full, g_Ess, 1 );

#elif defined THREEDIM

	BCHandler bc(6, BCHandler::HINT_BC_ONLY_ESSENTIAL);

	bc.addBC("Top", TOP, Essential, Full, g_Ess, 1);
	bc.addBC("Bottom", BOTTOM, Essential, Full, g_Ess, 1);
	bc.addBC("Left", LEFT, Essential, Full, g_Ess, 1);
	bc.addBC("Right", RIGHT, Essential, Full, g_Ess, 1);
	bc.addBC("Front", FRONT, Essential, Full, g_Ess, 1);
	bc.addBC("Back", BACK, Essential, Full, g_Ess, 1);
#endif

	//=============================================================================
	//Mesh stuff
	Members->comm->Barrier();
	DataMesh dataMesh(dataFile, ("bdf/" + discretization_section).c_str());
	boost::shared_ptr<RegionMesh> meshPtr( new RegionMesh() );
	readMesh(*meshPtr,dataMesh);
	partitionMesh<RegionMesh> meshPart(meshPtr, Members->comm);

	//=============================================================================
	//finite element space of the solution
	FESpace<RegionMesh, EpetraMap> FeSpace(meshPart,
			dataFile(("bdf/"+discretization_section).c_str(), "P2"), 1, Members->comm);

	if (verbose)
		std::cout << "  Number of unknowns : "
				<< FeSpace.map().getMap(Unique)->NumGlobalElements()
				<< std::endl;

	bc.bdUpdate(*(FeSpace.mesh()), FeSpace.feBd(), FeSpace.dof());

	//=============================================================================
	//Fe Matrices and vectors
	ElemMat elmat(FeSpace.fe().nbNode, 1, 1); //local matrix
	EpetraMatrix<double> matM(FeSpace.map()); //mass matrix
	boost::shared_ptr<EpetraMatrix<double> > matA_ptr(
			new EpetraMatrix<double> (FeSpace.map())); //stiff matrix
	EpetraVector u(FeSpace.map(), Unique); // solution vector
	EpetraVector f(FeSpace.map(), Unique); // forcing term vector

	Chrono chrono;
	//Assembling Matrix M
	Members->comm->Barrier();
	chrono.start();
	for (UInt iVol = 1; iVol <= FeSpace.mesh()->numElements(); iVol++)
	{
		FeSpace.fe().updateJac(FeSpace.mesh()->element(iVol));
		elmat.zero();
		mass(1., elmat, FeSpace.fe(), 0, 0);
		assembleMatrix(matM, elmat, FeSpace.fe(), FeSpace.fe(), FeSpace.dof(),
				FeSpace.dof(), 0, 0, 0, 0);
	}
	matM.GlobalAssemble();
	Members->comm->Barrier();
	chrono.stop();
	if (verbose)
		std::cout << "\n \n -- Mass matrix assembling time = " << chrono.diff() << std::endl
				<< std::endl;

	// ==========================================
	// Definition of the time integration stuff
	// ==========================================
	Real Tfin = dataFile("bdf/endtime", 10.0);
	Real delta_t = dataFile("bdf/timestep", 0.5);
	Real t0 = 1.;
	UInt ord_bdf = dataFile("bdf/order", 3);
	BdfVS<EpetraVector> bdf(ord_bdf);

	//Initialization
	bdf.initialize_unk(AnalyticalSol::u, u, FeSpace, t0, delta_t);
	if(verbose) bdf.showMe();
	Members->comm->Barrier();

	//===================================================
	// post processing setup
#if POSTPROCESS
	boost::shared_ptr<Exporter<RegionMesh> > exporter;
	exporter.reset(new Ensight<RegionMesh> (dataFile, meshPart.mesh(), "u",
			Members->comm->MyPID()));

	boost::shared_ptr<EpetraVector> u_display_ptr(new EpetraVector(
			FeSpace.map(), Repeated));
	exporter->addVariable(ExporterData::Scalar, "u", u_display_ptr,
                          UInt(0),
                          UInt(FeSpace.dof().numTotalDof()));
	*u_display_ptr = u;
	exporter->postProcess(0);
#endif

	//===================================================
	//Definition of the linear solver
	SolverTrilinos az_A(Members->comm);
	az_A.setDataFromGetPot(dataFile, "bdf/solver");
	az_A.setUpPrec(dataFile, "bdf/prec");

	//===================================================
	// TIME LOOP
	//===================================================
	for (Real t = t0 + delta_t; t <= Tfin; t += delta_t)
	{
		Members->comm->Barrier();
		if (verbose)
			cout << "Now we are at time " << t << endl;

		matA_ptr.reset(new EpetraMatrix<double> (FeSpace.map()));
//		matA_ptr->clear();	/*Reset the EpetraMatrix pointer until we discuss how to clear a EpetraMatrix*/

		chrono.start();
		//Assemble A
		Real coeff = bdf.coeff_der(0) / delta_t;
		Real visc = nu(t);
		Real s = sigma(t);
		for (UInt i = 1; i <= FeSpace.mesh()->numElements(); i++) {
			FeSpace.fe().updateFirstDerivQuadPt(FeSpace.mesh()->element(i));
			elmat.zero();
			mass(coeff + s, elmat, FeSpace.fe());
			stiff(visc, elmat, FeSpace.fe());
			assembleMatrix(*matA_ptr, elmat, FeSpace.fe(), FeSpace.fe(),
					FeSpace.dof(), FeSpace.dof(), 0, 0, 0, 0);
		}

		chrono.stop();
		if(verbose)
			cout << "A has been constructed in "<< chrono.diff() << "s." << endl;

		// Handling of the right hand side
		FeSpace.interpolate(sf, f, t);
		Members->comm->Barrier();
		f += bdf.time_der_dt();
		f = matM * f;

		// Treatment of the Boundary conditions
		if(verbose) cout << "*** BC Management: " << endl;
		Real tgv = 1.;
		chrono.start();
		bcManage(*matA_ptr, f, *FeSpace.mesh(), FeSpace.dof(), bc, FeSpace.feBd(), tgv, t);
		matA_ptr->GlobalAssemble();
		chrono.stop();
		if(verbose) cout << chrono.diff() << "s." << endl;

		Members->comm->Barrier();
		chrono.start();
		az_A.setMatrix(*matA_ptr);
		az_A.setReusePreconditioner(false);

		Members->comm->Barrier();
		az_A.solveSystem(f, u, matA_ptr);
		chrono.stop();

		bdf.shift_right(u);

		if (verbose)
			cout << "*** Solution computed in " << chrono.diff() << "s."
					<< endl;

		Members->comm->Barrier();
		// Error L2 and H1 Norms
		AnalyticalSol uExact;
		vector_type uComputed(u, Repeated);

		Real L2_Error, L2_RelError;

		L2_Error = FeSpace.L2Error(AnalyticalSol::u, uComputed, t, &L2_RelError);

		if (verbose)
			std::cout << "Error Norm L2: " << L2_Error
					<< "\nRelative Error Norm L2: " << L2_RelError << std::endl;

#if POSTPROCESS
		//store solution at time t
		*u_display_ptr = u;
		exporter->postProcess(t);
#endif
	}

}
