/**
 *  \file main.cpp
 *  NB : we follow the "test_cylinder" (ie mixt NS FEM test) structure :
 *      - short main
 *      - all problem defined in the xxx.{cpp,hpp,class} : elements, bdf, update params,
 *        bcs, ...
 */

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <boost/program_options.hpp>

#include <life/lifecore/life.hpp>
#include <life/lifecore/application.hpp>

#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>

#include <ct.hpp>
#include "mpi.h"


LifeV::AboutData
makeAbout()
{
    LifeV::AboutData about( "test_ct" ,
                            "test_chorintemam" ,
                            "0.1",
                            "3D test case",
                            LifeV::AboutData::License_GPL,
                            "Copyright (c) 2008 INRIA/EPFL");
    return about;

}

using namespace LifeV;



namespace
{
static bool regIF = (PRECFactory::instance().registerProduct( "Ifpack", &createIfpack ));
static bool regML = (PRECFactory::instance().registerProduct( "ML", &createML ));
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


    LifeV::po::options_description desc("Specific options");
    desc.add_options()("file,f", LifeV::po::value<std::string>()->default_value( "data" ), "data file name");

    CT ct( argc, argv, makeAbout(), desc );
    ct.run();

#ifdef HAVE_MPI
    MPI_Finalize();
#endif
  return( EXIT_SUCCESS );
}
