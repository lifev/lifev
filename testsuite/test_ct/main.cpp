/*!
 * \file main.cpp
 */

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include <life/lifecore/life.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include <ct.hpp>
#include <mpi.h>


namespace LifeV
{
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}
}

using namespace LifeV;



int
main( int argc, char** argv )
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  if ( Comm.MyPID() == 0 )
      cout << "  a-  Using MPI" << endl;
#else
  Epetra_SerialComm Comm;
  cout << "  a-  Using serial version" << endl;
#endif

    CT ct( argc, argv );
    ct.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  return( EXIT_SUCCESS );
}
