// Main programm for coupled calculation of fluid-dynamics
// and mass transport in the arterial lumen


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

#include <life/lifecore/life.hpp>

#include <life/lifealg/PreconditionerIfpack.hpp>
#include <life/lifealg/PreconditionerML.hpp>


#include "masstransport.hpp"


using namespace LifeV;

namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
}



int main(int argc, char** argv)
{

// MPI
#ifdef HAVE_MPI
    MPI_Init(&argc, &argv);
    Epetra_MpiComm Comm(MPI_COMM_WORLD);
    if ( Comm.MyPID() == 0 )
        cout << "% using MPI" << endl;
#else
    Epetra_SerialComm Comm;
    cout << "% using serial Version" << endl;
#endif



// Main

//**************** cylinder
//    MPI_Init(&argc,&argv);

    MassTransport mt( argc, argv );
    mt.run();




// MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
//
    return( 0 );
}
