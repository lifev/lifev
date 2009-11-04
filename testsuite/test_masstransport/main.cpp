// Main programm for coupled calculation of fluid-dynamics
// and mass transport in the arterial lumen


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>


#include "masstransport.hpp"
#include <mpi.h>


LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "Mass_Transport" ,
                            "Mass_Transport" ,
                            "0.1",
                            "Mass Transport Test Case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2005 EPFL");

    about.addAuthor("Gilles Fourestey", "developer", "gilles.fourestey@imag.fr", "");
    return about;

}


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

    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");
    MassTransport mt( argc, argv, makeAbout(), desc );
    mt.run();




// MPI
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
//
    return( 0 );
}
