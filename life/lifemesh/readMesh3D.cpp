/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano
  
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  
  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
  
  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "readMesh3D.hpp"
# include <stdlib.h>
# include <stdio.h>

// here I should cut and paste the end of the readMesh3D.h
// but I have troubles when linking (see the comments in readMesh3D.h)
//----------------------------------------------------------------------
// 
// Problem: the functions that follows should be in a readMesh3D.cc
// nevertheless, if I do that, I can compile the library, but I
// obtain several errors during the linking with a main.cc.
// JFG, 07/09/2002
//======================================================================

bool
readMppFileHead(ifstream & mystream, UInt & numVertices, UInt & numBVertices, UInt & numBFaces, UInt & numBEdges, UInt & numVolumes)
{
  unsigned done=0;
  string line;
  Real x,y,z;
  int ity,ity_id;
  UInt p1,p2,p3;
  UInt i,ibc;
  //streampos start=mystream.tellg();
  
  while (next_good_line(mystream,line).good()){
    if (line.find("odes") != string::npos){
      string node_s=line.substr(line.find_last_of(":")+1);
      numVertices=atoi(node_s);
      done++;
      numBVertices=0;
      for(i=0;i<numVertices;i++) {
	mystream >> x >> y >> z >> ity >> ibc;
	if(ity !=3 ){
#ifndef OLDMPPFILE
	  mystream>>ibc;
#endif
	  numBVertices++;
	}
      }
    }
    
    
    if (line.find("iangular") != string::npos){
      string node_s=line.substr(line.find_last_of(":")+1);
      numBFaces=atoi(node_s);
      done++;
      for(i=0;i<numBFaces;i++){
#ifdef OLDMPPFILE
	mystream>>p1 >> p2 >> p3 >> ity >> ibc;
#else
	mystream>>p1 >> p2 >> p3 >> ity >> ity_id>> ibc;
#endif
      }
    }
    
    if (line.find("oundary") != string::npos) {
      string node_s=line.substr(line.find_last_of(":")+1);
      numBEdges=atoi(node_s);
      for(i=0;i<numBEdges;i++){
#ifdef OLDMPPFILE
	mystream>>p1>>p2>>ity>>ibc;
#else
	mystream>>p1>>p2>>ity>> ity_id>>ibc;
#endif
      }
      done++;
    }
    if (line.find("etrahedral") != string::npos){
      string node_s=line.substr(line.find_last_of(":")+1);
      numVolumes=atoi(node_s);
      done++;
    }
    }
  return done==4 ;
}



// ****************      INRIA mesh  readers   **********************************



int nextIntINRIAMeshField(string const & line, istream & mystream)
{
  // first control if line has something. If so use atoi (the version from util_string.h) to extract
  // the integer. Otherwise get if from the input stream
  for (std::string::const_iterator is=line.begin(); is!= line.end(); ++is)
    if(*is != ' ') return atoi(line);
  int dummy;
  mystream >> dummy;
  return dummy;
}

//! Reads all basic info from INRIA MESH file
//! so as to be able to properly dimension all arrays
bool
readINRIAMeshFileHead(ifstream & mystream, UInt & numVertices, UInt & numBVertices, UInt & numBFaces, UInt & numBEdges, UInt & numVolumes, ReferenceShapes & shape)
{
  unsigned done=0;
  string line;
  Real x,y,z;
  UInt p1,p2,p3,p4,p5,p6,p7,p8;
  UInt i,ibc;
  
  int idummy;
  shape=NONE;
  //streampos start=mystream.tellg();
  
  while (next_good_line(mystream,line).good()){

    if (line.find("MeshVersionFormatted") != string:: npos){
      idummy=nextIntINRIAMeshField(line.substr(line.find_last_of("d")+1),mystream);
      ASSERT_PRE0(idummy == 1, "I can read only formatted INRIA Mesh files, sorry");
    }

    if (line.find("Dimension") != string:: npos){
      idummy=nextIntINRIAMeshField(line.substr(line.find_last_of("n")+1),mystream);
      ASSERT_PRE0(idummy == 3, "I can read only 3D INRIA Mesh files, sorry");
    }

    // I assume that internal vertices have their Ref value set to 0 (not clear from medit manual)
    if (line.find("Vertices") != string::npos){
      numVertices=nextIntINRIAMeshField(line.substr(line.find_last_of("s")+1),mystream);
      done++;
      numBVertices=0;
      for(i=0;i<numVertices;i++) {
	mystream >> x >> y >> z >> ibc;
	if(ibc !=0 ) numBVertices++;
      }
    }
    
    // I am assuming we are storing only boundary faces    
    if (line.find("Triangles") != string::npos){
      ASSERT_PRE0(shape != HEXA," Cannot have triangular faces in an HEXA INRIA  MESH");
      shape=TETRA;
      numBFaces=nextIntINRIAMeshField(line.substr(line.find_last_of("s")+1),mystream);
      done++;
      for(i=0;i<numBFaces;i++){
	mystream>>p1 >> p2 >> p3 >> ibc;
      }
    }
    
    if (line.find("Quadrilaterals") != string::npos){
      ASSERT_PRE0(shape != TETRA," Cannot have quad faces in an TETRA INRIA MESH");
      shape=HEXA;
      numBFaces=nextIntINRIAMeshField(line.substr(line.find_last_of("s")+1),mystream);
      done++;
      for(i=0;i<numBFaces;i++){
	mystream>>p1 >> p2 >> p3 >> p4>> ibc;
      }
    }
    // To cope with a mistake int INRIA Mesh files    
    if (line.find("Tetrahedra") != string::npos){
      ASSERT_PRE0(shape != HEXA," Cannot have tetras  in a HEXA INRIA MESH");
      shape=TETRA;
      numVolumes=nextIntINRIAMeshField(line.substr(line.find_last_of("a")+1),mystream);
      done++;
      for(i=0;i<numVolumes;i++){
	mystream>>p1 >> p2 >> p3 >> p4>>ibc;
      }
    }
    
    if (line.find("Hexahedra") != string::npos){
      ASSERT_PRE0(shape != TETRA," Cannot have Hexahedra in a TETRA INRIA MESH");
      shape=HEXA;
      numVolumes=nextIntINRIAMeshField(line.substr(line.find_last_of("a")+1),mystream);
      done++;
      for(i=0;i<numVolumes;i++){
	mystream>>p1 >> p2 >> p3 >> p4>> p5>>p6>>p7>>p8>>ibc;
      }
    }
    // I assume we are storing only boundary edges
    if (line.find("Edges") != string::npos){
      numBEdges=nextIntINRIAMeshField(line.substr(line.find_last_of("a")+1),mystream);
      done++;
      for(i=0;i<numBEdges;i++){
	mystream>>p1 >> p2 >>ibc;
      }
    }
  }
  return true ;
}

