/*
 * test_genAlpha.cpp
 *
 *  Created on: Jul 27, 2010
 *      Author: uvilla
 */

#include <Epetra_ConfigDefs.h>
#include <Epetra_Comm.h>

#include <lifev/core/LifeV.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/util/LifeChrono.hpp>


using namespace LifeV;

int run (GetPot& dataFile, bool verbose, boost::shared_ptr<Epetra_Comm>& comm);

// Do not edit
int main (int argc, char** argv)
{
    using namespace LifeV;
#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    std::cout << "MPI Initialization\n";
#endif

#ifdef EPETRA_MPI
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> comm (new Epetra_SerialComm() );
#endif
    bool verbose = comm->MyPID() == 0;

    std::string dataFileName;
    GetPot command_line (argc, argv);
    dataFileName = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (dataFileName);

    int check;
    check = run (dataFile, verbose, comm);

    comm.reset();

#ifdef HAVE_MPI
    MPI_Finalize();
    std::cout << "MPI Finalization \n";
#endif

    return check;
}
