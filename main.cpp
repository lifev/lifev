//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 * @file
 * @brief ContainerOfVectors Test
 *
 * @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 * @date 29-09-2009
 *
 *  This is a test to verify that the ContainerOfVectors works correctly.
 */

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
	#include <mpi.h>
#endif

#include <iomanip>
#include <string>

#include <life/lifecore/application.hpp>
#include <life/lifecore/life.hpp>

#include <life/lifealg/EpetraMap.hpp>
#include <life/lifearray/EpetraVector.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/shared_ptr.hpp>

#include "TestFunction.hpp"

LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "ContainerOfVectors Test" ,
                            "A ContainerOfVectors tester" ,
                            "1.0",
                            "Some practical examples about how the ContainerOfVectors is working!",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2009 EPFL");

    about.addAuthor("Cristiano Malossi", "Developer", "cristiano.malossi@epfl.ch", "");
    return about;
}

using namespace LifeV;

int
main( int argc, char** argv )
{
	#ifdef HAVE_MPI
		std::cout << "MPI Initialization" << std::endl;
		MPI_Init( &argc, &argv );
	#endif

	// ===================================================
	// TEST: LIFEV::EPETRAVECTOR
	// ===================================================

	//Setup main communicator
	boost::shared_ptr<Epetra_Comm>	comm;

	int	nprocs;
	int rank;

	//MPI Preprocessing
	#ifdef EPETRA_MPI

		MPI_Comm_size( MPI_COMM_WORLD, &nprocs );
		MPI_Comm_rank( MPI_COMM_WORLD, &rank );

		if ( rank == 0 )
		{	std::cout << "MPI processes: " << nprocs << std::endl;

			std::cout << "MPI Epetra Initialization ... " << std::endl;
		}
		comm.reset( new Epetra_MpiComm( MPI_COMM_WORLD ) );

		comm->Barrier();

	#else

		std::cout << "MPI SERIAL Epetra Initialization ... " << std::endl;
		comm.reset( new Epetra_SerialComm() );
	#endif

    std::cout << std::setprecision(15) << std::endl;

	//BASE VECTORS
	typedef EpetraVector                        Vector;
	typedef boost::shared_ptr<Vector>           Vector_ptr;

	int MyGlobalIElementsA[3],MyGlobalIElementsB[2];
	MyGlobalIElementsA[0] = 0;  MyGlobalIElementsA[1] = 1;  MyGlobalIElementsA[2] = 2;
	MyGlobalIElementsB[0] = 0;  MyGlobalIElementsB[1] = 1;

	EpetraMap mapA( 3, 3, &MyGlobalIElementsA[0], 0, comm);
	EpetraMap mapB( 2, 2, &MyGlobalIElementsB[0], 0, comm);

	Vector_ptr A1, B1, A2, B2, A3, B3, A4, B4;
	A1.reset( new Vector( mapA, Unique ) );
	B1.reset( new Vector( mapB, Unique ) );
	A2.reset( new Vector( mapA, Unique ) );
	B2.reset( new Vector( mapB, Unique ) );
	A3.reset( new Vector( mapA, Unique ) );
	B3.reset( new Vector( mapB, Unique ) );
    A4.reset( new Vector( mapA, Unique ) );
    B4.reset( new Vector( mapB, Unique ) );

	//Fill the vectors
	for ( UInt i(0) ; i < static_cast<UInt> (A1->size()) ; ++i )
		(*A1)[i] = i;
	for ( UInt i(0) ; i < static_cast<UInt> (B1->size()) ; ++i )
		(*B1)[i] = A1->size() + i;

	for ( UInt i(0) ; i < static_cast<UInt> (A2->size()) ; ++i )
		(*A2)[i] = i*10;
	for ( UInt i(0) ; i < static_cast<UInt> (B2->size()) ; ++i )
		(*B2)[i] = (A2->size() + i)*10;

    for ( UInt i(0) ; i < static_cast<UInt> (A3->size()) ; ++i )
        (*A3)[i] = i*10;
    for ( UInt i(0) ; i < static_cast<UInt> (B3->size()) ; ++i )
        (*B3)[i] = (A3->size() + i)*100;

    for ( UInt i(0) ; i < static_cast<UInt> (A4->size()) ; ++i )
        (*A4)[i] = i*10;
    for ( UInt i(0) ; i < static_cast<UInt> (B4->size()) ; ++i )
        (*B4)[i] = (A4->size() + i)*1000;

	UInt result=TestFunction( A1, B1, A2, B2, A3, B3, A4, B4 );



	// ===================================================
	// TEST: BOOST::NUMERIC::UBLAS::VECTOR
	// ===================================================

	typedef boost::numeric::ublas::vector<Real> VectorBOOST;
	typedef boost::shared_ptr<VectorBOOST>      VectorBOOST_ptr;

	// BASE Vectors
	VectorBOOST_ptr A1_BOOST, B1_BOOST, A2_BOOST, B2_BOOST;

	A1_BOOST.reset( new VectorBOOST(2) );
	(*A1_BOOST)[0] = 0;
	(*A1_BOOST)[1] = 1;

	B1_BOOST.reset( new VectorBOOST(3) );
	(*B1_BOOST)[0] = 2;
	(*B1_BOOST)[1] = 3;
	(*B1_BOOST)[2] = 4;

	A2_BOOST.reset( new VectorBOOST(2) );
	(*A2_BOOST)[0] = 0;
	(*A2_BOOST)[1] = 10;

	B2_BOOST.reset( new VectorBOOST(3) );
	(*B2_BOOST)[0] = 20;
	(*B2_BOOST)[1] = 30;
	(*B2_BOOST)[2] = 40;

	// NOTE:
	// To perform a test with BOOST we first have to implement a derived class
	// from ContainerOfVectors in which re-implement operator= and operator*

	//TestFunction( A1_BOOST, B1_BOOST, A2_BOOST, B2_BOOST );

	#ifdef HAVE_MPI
		std::cout << std::endl << "MPI Finalization" << std::endl;
		MPI_Finalize();
	#endif

	return( result );
}
