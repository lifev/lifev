#include "life.hpp"
#include "NavierStokesAleSolverPC.new.hpp"
#include "regionMesh3D_ALE.hpp"
#include "ud_functions.hpp"
#include "bcVector.hpp"

#define PVM
//#undef PVM

#ifdef PVM
#include <pvm3.h>
#endif

int fsiStatus=-1;

#ifdef PVM
int masterId;
#endif

int main(int argc, char** argv)
{

  using namespace LifeV;
  using namespace std;


#ifdef PVM
  cout << "Je suis au debut du fluid" << endl << flush;
#endif

  // Reading from data file
  //
  GetPot command_line(argc,argv);
  const char* data_file_name = command_line.follow("data", 2, "-f","--file");
  GetPot data_file(data_file_name);

  // Number of boundary conditions for the mesh motion and velocity
  //

#ifdef PVM
  BCHandler BCh_u;
#else
  BCHandler BCh_u;
#endif
  BCHandler BCh_mesh;



  // NavierStokes ALE solver
  NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > nsAle(data_file, feTetraP1bubble,
								 feTetraP1,quadRuleTetra64pt,
								 quadRuleTria3pt, quadRuleTetra64pt,
								 quadRuleTria3pt, BCh_u,BCh_mesh);

  nsAle.showMe();


  // partie initialisation pour pvm
#ifdef PVM
  masterId = pvm_parent();
  if(masterId==PvmNoParent){
    cout << "I am a poor lonesome job\n";
    exit(1);
  } else
    if(masterId==0){
      cout << "PVM is not ok\n";
      exit(1);
    }
  cout <<"%%%%%%%%%%% \n";
  cout << " I am the slave " << pvm_mytid()<< endl;
  cout << " My master is " << masterId << endl;
  cout << " Avant rv " << endl << flush;
  nsAle.rvFromMasterFSI(140);
#endif


  // Boundary conditions for the fluid velocity u
  //
  cout << "Main - BC" << endl << flush;

  BCFunctionBase bcf(fZero); // BCFunction object holding the null function



#ifdef PVM
  BCFunctionBase in_flow(g);
#endif


  // BCVector object holding the velocity of the fluid domain
  BCVector u_wall(nsAle.wInterpolated(),nsAle.uDof().numTotalDof());

  BCh_u.addBC("Wall",    1,  Essential, Full, u_wall,  3);
  BCh_u.addBC("Wall",   20, Essential, Full, u_wall,   3);
  BCh_u.addBC("Wall",   10, Essential, Full, u_wall,   3);
  BCh_u.addBC("OutFlow", 3,  Natural,   Full, bcf, 3);

#ifdef PVM
  BCh_u.addBC("InFlow", 2,  Natural,   Full, in_flow, 3);
#else
  BCh_u.addBC("InFlow", 2,  Natural,   Full, bcf, 3);
#endif




  // Boundary conditions for the mesh displacement (harmonic extension)

  UInt dim_mesh = nsAle.dofMesh().numTotalDof(); // Number of Degree of Freedom DOF for the displacement of the fluid mesh (number of points

  Vector analytic_disp(3*dim_mesh); // The vector holding the displacement of the fluid domain
  analytic_disp = ZeroVector( analytic_disp.size() ); // initalization to zero


  BCVector dispbd(analytic_disp,dim_mesh); // The BCVector object holding the displacement of the fluid domain d^s
  vector<ID> compDir(1); // array containg the number of the component corresponding to the axial direction
  compDir[0]=3;
  BCh_mesh.addBC("Wall", 1, Essential, Full, dispbd, 3);  //  displacement given on the lateral points
  BCh_mesh.addBC("Wall", 20, Essential, Full, dispbd, 3);  //   "       "
  BCh_mesh.addBC("Wall", 10, Essential, Full, dispbd, 3);
  BCh_mesh.addBC("Axis", 2, Essential, Component, bcf, compDir); // Fixing axial ("3th-component) displacement
                                                                 // on the upper and bottom faces
  BCh_mesh.addBC("Axis", 3, Essential, Component, bcf, compDir); //     "                               "


  // Initialization
  //
  cout << "Main - Initialisation" << endl << flush;
  //
  Real dt = nsAle.timestep();  // The time step
  Real T  = nsAle.endtime();   // The end time
  nsAle.initialize(u0); // Initial condition

  Real x,y,z;

  //case of f-s interaction
  int max_iter; max_iter=100;
  int iter=0,istep=1;
  Real time=dt; // pas de temps a regler

do{

 cout << "Main : boucle -do-" << endl << flush;

#ifdef PVM
     nsAle.rvFromMasterFSI(110);
#else
    fsiStatus=1;
#endif
    if(fsiStatus == -1) { cout << "Normal fluid exit\n" << flush;  exit(1); }
    if(fsiStatus == 1) { // nouveau pas de temps
#ifdef PVM
      if ( iter >= 1 )
	nsAle.postProcess();
#endif
      nsAle.timeAdvance(f,time);
      cout << "PVM's timeAdvance ok!" << endl << flush;
      // Updating the given imposed displacement at this time step
#ifdef PVM
      nsAle.dep_fluid_interf(analytic_disp);
#else
      for (UInt i=0; i < dim_mesh; ++i) {
	if ( nsAle.mesh().point(i+1).marker() == 1 ||
	     nsAle.mesh().point(i+1).marker() == 20 ||
	     nsAle.mesh().point(i+1).marker() == 10   ) { // only on lateral points
	  x=nsAle.mesh().point(i+1).x();
	  y=nsAle.mesh().point(i+1).y();
	  z=nsAle.mesh().point(i+1).z();
	  for (UInt j=0; j<3; ++j)
	    analytic_disp(i + j*dim_mesh) = bdDisp(time, x, y, z, j+1);
	}
      }
#endif
      cout << "analytic_disp ok!" << endl << flush;
      nsAle.updateMesh(time); // compute the harmonic extension d^f and updating the mesh points
      nsAle.iterate(time);

#ifdef PVM
      nsAle.sdToMasterFSI();
#else
      nsAle.postProcess();
#endif

      time+=dt; //pas de temps a regler
      iter++;
    }
#ifdef PVM
    if(fsiStatus == 0) {  // on reste au meme pas de temps (iterations internes)
      // depl maillage au pas n-1 ??? a-t-on vraiment besoin ??? NON
      // Updating the given imposed displacement at this time step
      nsAle.dep_fluid_interf(analytic_disp);
      cout << "analytic_disp for the sub-iteration ok!" << endl << flush;
      nsAle.updateMesh(); // compute the harmonic extension d^f and updating the mesh points
      nsAle.iterate();
      nsAle.sdToMasterFSI();
      istep++;
    }
#endif
}
#ifdef PVM
    while ( 1 ); // c'est le maitre qui decide en renvoyant fsiStatus = -1
#else
    while( time<= T  && iter <= max_iter );
#endif

    return 0;
}

