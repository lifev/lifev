
#include "markers.hpp"
#include "readMesh3D.hpp"
#include "regionMesh3D.hpp"
#include "bareItems.hpp"

int main()
  {
    
    string dirname="../data/mesh/mesh++/";
    string fname=dirname+"cube_47785.m++";
    string outfile="testBuilders.dat";
    ofstream ofile(outfile.c_str());
    if (ofile.fail()) {cerr<<" Error: Cannot creat output file"<<endl; abort();}
  
    //RegionMesh3D<LinearTetra> aMesh;

    //template class RegionMesh3D<QuadraticTetra>; /gnu c+ +does not support it yet
    
    RegionMesh3D<QuadraticTetra> aMesh;
    
    //    aMesh.test3DBuilder();
    //    aMesh.readMppFile(mystream, id, m);
    ID m=1;
    readMppFile(aMesh,fname, m);

    cout <<" **********************************************************"<<endl;
    
    cout<< "  ****************** CHECKING MESH WITH INTERNAL CHECKER"<<endl;
    
    aMesh.check(0,true,true);
    Switch sw;

    
    cout <<" **********************************************************"<<endl;
    
    cout<< "  ****************** CLEANING edges and faces "<<endl;
    
    aMesh.edgeList.clear();
    aMesh.faceList.clear();
    
    checkMesh3D(aMesh, sw,
		       true, true,
		cerr, cout, ofile);
    

    aMesh.showMe();

    
    cout<< " Now building local Edges/faces Stuff"<<endl<<endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    cout<< " Now cleaning local Edges/faces Stuff"<<endl<<endl;
    aMesh.cleanElementEdges();
    aMesh.cleanElementFaces();
    aMesh.showMe();
    UInt bedges_found, iedges_found;
    cout <<" **********************************************************"<<endl;
    
    cout<< "  ****************** BUILDING ALL EDGES"<<endl;
    buildEdges(aMesh, ofile, cerr, bedges_found,
	       iedges_found, true, true,true);
    
    cout <<" **********************************************************"<<endl;
    
    cout<< "  ****************** BUILDING ALL FACES"<<endl;
    UInt bfaces_found, ifaces_found;
    buildFaces(aMesh, ofile, cerr, bfaces_found,
	       ifaces_found, true,true,true);
    
    aMesh.showMe();

    cout<< " Now building again local Edges/faces Stuff"<<endl<<endl;
    aMesh.updateElementEdges();
    aMesh.updateElementFaces();
    aMesh.showMe();
    checkMesh3D(aMesh, sw,
		       true, true,
		cerr, cout, ofile);
    ofile.close();

    ///THIS PART IS ONLY TO VERIFY IF THESE ROUTINES COMPILE PROPERLY
      cerr<<"Fixing bpoints"<<endl;
      
      fixBPoints(aMesh, ofile,cerr, true);
      cerr<<"Fixing edge markers"<<endl;
      setBEdgesMarker(aMesh, ofile,cerr, true);
      cerr<<"Fixing faces marker"<<endl;
      setBFacesMarker(aMesh, ofile,cerr, true);
      cerr<<"Fixing points marker"<<endl;
      setBPointsMarker(aMesh, ofile,cerr, true);
      cerr<<endl;
      vector<Real> disp(3*aMesh.numPoints());
      aMesh.moveMesh(disp);
      
      //    SimpleVect<bool>elsign;
    //aMesh.volume(2).swapPoints(1,3);
    //UInt lp=GeoND<QuadraticTetra>::numLocalPoints;
    //    UInt lp=GeoND<QuadraticTetra>::numLocalPoints;
    //for(UInt i = 1; i<=aMesh.numVolumes(); i++){
      //      for(UInt j=1; j<=lp;j++){
      //cout << "Point " << j << ":"
      //  
      //     << "("<< aMesh.volumeList(i).point(j).id() << ") ,"
      //     << aMesh.volumeList(i).point(j).coordinate(1) <<","
      //     << aMesh.volumeList(i).point(j).coordinate(2) <<","
      //     << aMesh.volumeList(i).point(j).coordinate(3) << endl;}}
    //    UInt i1,i2,i3,i4;
    
    //for (UInt i=1; i<=aMesh.numVolumes(); ++i){
      //  cout << "Element # "<< i << ": "<<aMesh.volume(i).point(1).id()<<","
      //   <<aMesh.volume(i).point(2).id()
    //   <<","<<aMesh.volume(i).point(3).id()<<","
    //   <<aMesh.volume(i).point(4).id()<<endl;
    //for(UInt k=1;k<=aMesh.numLocalEdges();++k){
    //i1=aMesh.localEdgeId(i,k);
    //cout << i1<< "-";
    //}
    //cout<<endl;
    //}
    
//for (UInt i=1; i<=aMesh.numVolumes(); ++i){
//    cout << "Element # "<< i << ": "<<aMesh.volume(i).point(1).id()<<","
//   <<aMesh.volume(i).point(2).id()
//   <<","<<aMesh.volume(i).point(3).id()<<","
//   <<aMesh.volume(i).point(4).id()<<endl;
//    for(UInt k=1;k<=aMesh.numLocalFaces();++k){
//i1=aMesh.localFaceId(i,k);
    //cout << i1<< "-";
    //}
    //cout<<endl;
    //}

    // TESTING FIELDS
    //Dof<FE_P1_Tetra_Base,LinearTetra> dof;
    //Dof<FE_P2_Tetra_Base,LinearTetra> dof;
      //FiniteEle<FE_P2_Tetra_4pt> fe;
    ///Dof<FE_P2_Tetra_Base,QuadraticTetra> dof;
    //Dof<FE_P2_Tetra,QuadraticTetra> dof();
    //dof.update(aMesh);
    // dof.showMe();
    //    for (UInt i=1; i<=aMesh.numVolumes(); ++i){
    //cout << "Element # "<< i << ": "<<aMesh.volume(i).point(1).id()<<","
    //   <<aMesh.volume(i).point(2).id()
    //   <<","<<aMesh.volume(i).point(3).id()<<","
    //   <<aMesh.volume(i).point(4).id()<<endl;
    //for(UInt k=1;k<=dof.numLocalDof();++k){
    //i1=dof.localToGlobal(i,k);
    //cout << i1<< "-";
    //}
    //cout<<endl;
    //}
    /*
    RegionMesh3D<LinearTetra> aMesh;
    aMesh.test3DBuilder();
    aMesh.showMe();
    */
  }
 
