#include <iostream>
#include "GetPot.hpp"
#include "darcySolver.hpp"
#include "chrono.hpp"

using namespace std;

/*
  
  Darcy soler using Mixed Hybrid finite element

  
  usage: darcy              : read the data file "data" and run the solver
         darcy -f otherdata : read the data file "otherdata" and run the solver
	 darcy -h           : read the data file, print help and exit
	 darcy -i           : read the data file, print the read values and exit

*/
int main(int argc, char** argv)
{
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);
  if( command_line.search(2, "-i","--info") ) {
    data_file.print();
    exit(0);
  }
  Chrono chrono;
  //
  cout << "*** Initialisation --->" << endl;
  chrono.start();
  DarcySolver pb(data_file);
  if(command_line.search(2,".hpp","--help")){
    cout << endl << endl;
    cout <<"usage: darcy              : read the data file 'data' \n";
    cout <<"       darcy -f otherdata : read the data file 'otherdata' \n";
    cout <<"	   darcy -h           : help and exit\n";
    cout <<"       darcy -i           : : read the data file, print the read values and exit\n";
    cout << endl;
    cout << "Help for the data file:\n";
    pb.dataAztecHelp();
    pb.dataAztecShowMe();
    pb.dataDarcyHelp();
    pb.dataDarcyShowMe();
    exit(0);
  }
  //
  chrono.stop();
  cout << "<--- Initialisation done in " << chrono.diff() << "s." << endl; 
  
  if(pb.verbose)
    cout << "*** Compute the matrix --->" << endl;
  chrono.start();
  pb.computeHybridMatrix();
  chrono.stop();
  if(pb.verbose)
    cout << "<--- matrix computation done in "<< chrono.diff() << "s." << endl;
  pb.applyBC();
  //
  if(pb.verbose) cout << "*** Resolution of the hybrid system --->\n";
  chrono.start();
  pb.solveDarcy();
  chrono.stop();
  if(pb.verbose)cout << "<--- Linear system solved in " << chrono.diff()
        << "s." << endl << endl;
  
  if(pb.verbose) cout << "*** Compute pressure and flux --->" << endl;
  chrono.start();
  pb.computePresFlux();
  chrono.stop();
  if(pb.verbose)
    cout << "<---  done in " << chrono.diff() << "s." << endl << endl;
  if(pb.verbose) cout << "*** Postproc --->" << endl;
  chrono.start();
  pb.postProcessPressureQ0();
  pb.postProcessPressureQ1();
  pb.postProcessVelocityQ1();
  chrono.stop();
  if(pb.verbose)
    cout << "<---  done in " << chrono.diff() << "s." << endl << endl;
  //
  
}


