#include "lifeV.hpp"
#include "NavierStokesSolverPC.hpp"

#include <vector>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <string>

namespace LifeV
{
class Coord{
public:
  Coord(){init(0);};
  Coord &operator=(const Coord & A){init(A.coord);return *this;};
  float &operator[](int i){return coord[i];};

private:

  void init(const float* a){
    for(int i=0;i<3;++i)
      coord[i]=(a != 0) ? a[i]:-100;};

  float coord[3];
};


template <typename RegionMesh3D>
bool outensight7Mesh3D(RegionMesh3D & mesh, PhysVectUnknown<Vector> & u, ScalUnknown<Vector> & p, Real & time){

   vector<Coord> grid; // coordinates of Grid nodes
  grid.resize(mesh.numVertices()+mesh.numVolumes());
//  grid.resize(mesh.numVertices()+mesh.numEdges());


  for(ID i=0; i < mesh.numVertices(); i++){
     grid[i][0]=(float)mesh.point(i+1).x();
     grid[i][1]=(float)mesh.point(i+1).y();
     grid[i][2]=(float)mesh.point(i+1).z();}
//     cout << grid[i][0] << ", " << grid[i][1]  << ", " << grid[i][2] << endl; }

// 6 additional mesh points for Tetra P2 (mittle points of edges)
//  typename RegionMesh3D::EdgeType * pe=0;

//  ID i1,i2;

//  for(ID i=1; i <= mesh.numEdges(); i++){
//     pe=& mesh.edge(i);
//     i1=(pe->point(1)).id();
//     i2=(pe->point(2)).id();
//     grid[mesh.numVertices()+i-1][0]=(float)(mesh.point(i1).x()+mesh.point(i2).x())*0.5;
//     grid[mesh.numVertices()+i-1][1]=(float)(mesh.point(i1).y()+mesh.point(i2).y())*0.5;
//     grid[mesh.numVertices()+i-1][2]=(float)(mesh.point(i1).z()+mesh.point(i2).z())*0.5;}
//     cout << grid[i][0] << ", " << grid[i][1]  << ", " << grid[i][2] << endl;}


  char buffer[80];
  vector <int> idnode, idelem;

  fstream File3("test.geo",ios::out | ios::binary);

  strcpy(buffer,"C Binary");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"test cube - inria mesh");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"LinearTetra - elements");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"node id given");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"element id given");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"coordinates");
  File3.write((char *)&buffer,sizeof(buffer));
  int i=mesh.numVertices();
//  int i=mesh.numVertices()+mesh.numVolumes();
//  int i=mesh.numVertices() + mesh.numEdges();
  File3.write((char *)&i,sizeof(int));
  for (ID i = 1; i <= mesh.numVertices(); i++){  // load node id-vector for nodes
//  for (ID i = 1; i <= mesh.numVertices()+mesh.numVolumes(); i++){  // load node id-vector for nodes
//  for (ID i = 1; i <= mesh.numVertices()+mesh.numEdges(); i++){  // load node id-vector for nodes
    idnode.push_back(i);
  }
  File3.write((char *)&idnode.front(),idnode.size()*sizeof(int));
  File3.write((char *)&grid.front(),3*mesh.numVertices()*sizeof(float));
//  File3.write((char *)&grid.front(),3*(mesh.numVertices()+mesh.numVolumes())*sizeof(float));
//  File3.write((char *)&grid.front(),3*(mesh.numVertices()+mesh.numEdges())*sizeof(float));
  strcpy(buffer,"part 1");
  File3.write((char *)&buffer,sizeof(buffer));
  strcpy(buffer,"description line");
  File3.write((char *)&buffer,sizeof(buffer));
  if(mesh.numLocalVertices() == 4)
     strcpy(buffer,"tetra4");
//     strcpy(buffer,"tetra10");
  else
     strcpy(buffer,"hexa8");
  File3.write((char *)&buffer,sizeof(buffer));
  int e=mesh.storedVolumes();
//  int e=4*mesh.storedVolumes();
  File3.write((char *)&e,sizeof(int));
  for (int i = 1; i <= e; i++){  // load node id-vector for elements
    idelem.push_back(i);
  }
  File3.write((char *)&idelem.front(),idelem.size()*sizeof(int));

// ------------------   Output P1 Elements (Tetra or Hexahedra)   -------------------------

  for( ID k = 0; k < mesh.storedVolumes(); k++){
      for (ID j = 0; j < mesh.numLocalVertices(); j++){
	File3.write((char *)&mesh.volume(k+1).point(j+1).id(),sizeof(ID));
      }
  }

// ------------------     Output P2 Elements (Tetra)     -----------------------------------

// for( ID k = 0; k < mesh.storedVolumes(); k++){
//      for (ID j = 0; j < mesh.numLocalVertices(); j++){
//	File3.write((char *)&mesh.volume(k+1).point(j+1).id(),sizeof(ID));
//      }
//      for (ID j = 0; j < mesh.numLocalEdges(); j++){
//	 ID i = mesh.localEdgeId(k+1,j+1)+mesh.numVertices();
//	File3.write((char *)&i,sizeof(ID));
//      }
//  }


  File3.close();

  cout << "output in ensight7 format" << endl;
  cout << "geometry file is test.geo" << endl;

// read pressure von acsii-file ./Post/presQ1.bb and convert it into ensight7 binary format

  vector<float> pressure; // presure values at Grid nodes
  vector<Coord> velocity; // velocity values at Grid nodes

  ostringstream index;
  string name, vname;
  Coord nodevel;

    index << (time*100);

    switch( index.str().size() ) {
    case 1:
      name = "pressure.res00"+index.str();
      vname = "velocity.res00"+index.str();
      break;
    case 2:
      name = "pressure.res0"+index.str();
      vname = "velocity.res0"+index.str();
      break;
    case 3:
      name = "pressure.res"+index.str();
      vname = "velocity.res"+index.str();
      break;
    }

    for(ID i=0; i< p.size(); i++){
       pressure.push_back((float)(p(i)));}

    for(ID i=0; i< mesh.numVertices(); i++){
       //   for(ID i=0; i< mesh.numVertices()+mesh.numEdges(); i++){
       nodevel[0]=(float)u(i);
       nodevel[1]=(float)u(i+(u.size()/3));
       nodevel[2]=(float)u(i+2*(u.size()/3));
       velocity.push_back(nodevel);}

  std::fstream File4(name.c_str(),ios::out | ios::binary);

  strcpy(buffer,"concentration distribution timestep ");
  File4.write((char*)&buffer,sizeof(buffer));

  File4.write((char*)&pressure.front(),p.size()*sizeof(float));

  File4.close();

  std::fstream File5(vname.c_str(),ios::out | ios::binary);

  strcpy(buffer,"velocity field  timestep 1");
  File5.write((char*)&buffer,sizeof(buffer));

  File5.write((char*)&velocity.front(),3*velocity.size()*sizeof(float));

  File5.close();

  cout << "result files are velocity.res*** and pressure.res*** " << endl;

  return true;

}

}
