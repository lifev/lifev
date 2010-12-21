/*!
 * \file main.cpp
 */

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

//Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <life/lifecore/Life.hpp>

#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>

#include <ctRK.hpp>


using namespace LifeV;

namespace LifeV
{
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}
}


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

    CTRK ctRK( argc, argv );
    ctRK.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return( EXIT_SUCCESS );
}
