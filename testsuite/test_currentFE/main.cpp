//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Test for the CurrentFE update

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 16-04-2010

    Just a fast test to see if the update procedure is
    working correctly in the CurrentFE class.
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
    #include <mpi.h>
    #include <Epetra_MpiComm.h>
#else
    #include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <life/lifefem/currentFE.hpp>


using namespace LifeV;

int
main( int argc, char** argv )
{
    //MPI communicator initialization
    boost::shared_ptr<Epetra_Comm> comm;

#ifdef HAVE_MPI
    std::cout << "MPI Initialization" << std::endl;
    MPI_Init( &argc, &argv );
#endif

    //MPI Preprocessing
#ifdef EPETRA_MPI

    int nprocs;
    int rank;

    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    if ( rank == 0 )
    {
        std::cout << "MPI processes: " << nprocs << std::endl;
        std::cout << "MPI Epetra Initialization ... " << std::endl;
    }
    comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

    comm->Barrier();

#else

    std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
    comm.reset( new Epetra_SerialComm() );

#endif

    
// Test the current FE...
	
	Real test_tolerance(1e-10);	

// We define a fictious tetra
	std::vector<Real> P0 (3); P0[0]=1; P0[1]=0; P0[2]=0;
	std::vector<Real> P1 (3); P1[0]=1.1; P1[1]=0; P1[2]=0;
	std::vector<Real> P2 (3); P2[0]=1; P2[1]=0.2; P2[2]=0;
	std::vector<Real> P3 (3); P3[0]=1; P3[1]=0; P3[2]=0.15;
	
	std::vector< std::vector< Real > > Tetra1(4);
	Tetra1[0]=P0; Tetra1[1]=P1; Tetra1[2]=P2; Tetra1[3]=P3;

// 1. Part: flags & update

	std::cout << " Checking updates ... " << std::endl;
	
	CurrentFE test_CFE(feTetraP2,geoLinearTetra,quadRuleTetra4pt);

	Real notInterestedIn(0.0);

	test_CFE.update(Tetra1,UPDATE_QUAD_NODES);
	notInterestedIn = test_CFE.quadNode(0,0);

	test_CFE.update(Tetra1,UPDATE_PHI);
	notInterestedIn = test_CFE.phi(1,1);

	test_CFE.update(Tetra1,UPDATE_DPHI);
	notInterestedIn = test_CFE.dphi(2,0,3);

	test_CFE.update(Tetra1,UPDATE_D2PHI);
	notInterestedIn = test_CFE.d2phi(1,0,1,2);

	test_CFE.update(Tetra1,UPDATE_WDET);
	notInterestedIn = test_CFE.wDetJacobian(0);

// 2. Part: check values
	
	// Check the partition of unity

	std::cout << " Checking partition of unity ... " << std::endl;
	
	test_CFE.update(Tetra1,UPDATE_PHI | UPDATE_DPHI);
	Real sum_phi;
	Real sum_dphi;
	for (UInt i(0); i<test_CFE.nbQuadPt(); ++i)
	{
		sum_phi=0.0;
		sum_dphi=0.0;

		for (UInt n(0); n<test_CFE.nbFEDof(); ++n)
		{
			sum_phi += test_CFE.phi(n,i);
			sum_dphi += test_CFE.dphi(n,0,i);
		}

		if (std::abs(sum_phi-1) > test_tolerance)
		{		
			std::cerr << " Sum of the basis functions : " << sum_phi << std::endl;	
			return EXIT_FAILURE;
		};
		if (std::abs(sum_dphi) > test_tolerance)
		{	
			std::cerr << " Sum of the derivatives of the basis functions : " << sum_dphi << std::endl;	
			return EXIT_FAILURE;
		};
	}

	// Check the measures

	std::cout << " Checking measure ... " << std::endl;

	test_CFE.update(Tetra1,UPDATE_WDET );

	Real meas(test_CFE.measure());
	if (std::abs(meas-0.0005) > test_tolerance)
	{
		std::cerr << " Measure : " << meas << std::endl;	
		return EXIT_FAILURE;
	}


// 3. Part: change of quadrature

	std::cout << " Checking setQuadRule ... " << std::endl;
	test_CFE.setQuadRule(quadRuleTetra5pt);

	test_CFE.update(Tetra1,UPDATE_QUAD_NODES);
	notInterestedIn = test_CFE.quadNode(3,0);

	test_CFE.update(Tetra1,UPDATE_PHI);
	notInterestedIn = test_CFE.phi(1,4);

	test_CFE.update(Tetra1,UPDATE_DPHI);
	notInterestedIn = test_CFE.dphi(2,0,4);

	test_CFE.update(Tetra1,UPDATE_D2PHI);
	notInterestedIn = test_CFE.d2phi(1,0,1,4);

	test_CFE.update(Tetra1,UPDATE_WDET);
	notInterestedIn = test_CFE.wDetJacobian(4);


    // ----- End of test calls -----

#ifdef HAVE_MPI
    std::cout << "MPI Finalization" << std::endl;
    MPI_Finalize();
#endif

    // If everything runs correctly we return EXIT_SUCCESS flag
    return EXIT_SUCCESS;
}
