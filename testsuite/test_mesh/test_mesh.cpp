/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#include <life/lifecore/GetPot.hpp>

#include <life/lifemesh/markers.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifemesh/bareItems.hpp>


int main()
{
    using namespace LifeV;
    using namespace std;

    GetPot datafile( "data" );
    string dirname=datafile( "mesh_dir","." );//"../data/mesh/mesh++/";
    string fname=dirname+datafile( "mesh_file","cube_47785.m++" );//dirname+"cube_47785.m++";
    string outfile="testBuilders.dat";
    ofstream ofile(outfile.c_str());
    if (ofile.fail()) {cerr<<" Error: Cannot creat output file"<<endl; abort();}

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
    cout <<" **********************************************************"<<endl;

    cout<< "  ****************** BUILDING ALL EDGES"<<endl;
    UInt bedges_found, iedges_found;
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
}

