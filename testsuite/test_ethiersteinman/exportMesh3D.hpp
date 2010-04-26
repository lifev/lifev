//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains methods to export mesh into different format.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 26 Apr 2010

    This file contains the following export formats:
    - .mesh (default)
    - .m, Matlab/Octave script that generates a structure mesh which contains points, edges, faces, elements and flags.
 */

#ifndef EXPORTMESH3D_H
#define EXPORTMESH3D_H 1

#include <life/lifecore/life.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <fstream>
#include <string>

namespace LifeV {

//! @name Methods
//@{

#define MESH_FORMAT 1
#define MATLAB_FORMAT 2

//! Export the mesh in the specified format.
/*!
  @param mesh mesh to export
  @param format export format


*/
template <typename GeoShape, typename MC>
void exportMesh3D( RegionMesh3D<GeoShape,MC>& mesh, const std::string& fileName, const UInt& format=MESH_FORMAT ){
    std::string outName(fileName);
    if((format&MATLAB_FORMAT)!=0)
    {
        // Export for Matlab
        outName.append(".m");
        std::ofstream ofile (outName.c_str());
        if (ofile.is_open())
        {
            ofile << "#Output of the mesh structure.\n";

            // Writing the points
            ofile << "mesh.p=[\n";
            for(UInt i(1);i<=mesh.storedPoints();++i){
                ofile << mesh.point(i).x() << ",";
                ofile << mesh.point(i).y() << ",";
                ofile << mesh.point(i).z() << ";\n";
            }
            ofile << "];\n";

            ofile << "mesh.pflag=[\n";
            for(UInt i(1);i<=mesh.storedPoints();++i){
                ofile << mesh.point(i).marker() << ";\n";
            }
            ofile << "];\n";

            // Writing the edges
            ofile << "mesh.e=[\n";
            for(UInt i(1);i<=mesh.storedEdges();++i){
                ofile << mesh.edge(i).point(1).id() << ",";
                ofile << mesh.edge(i).point(2).id() << ";\n";
            }
            ofile << "];\n";

            ofile << "mesh.eflag=[\n";
            for(UInt i(1);i<=mesh.storedEdges();++i){
                ofile << mesh.edge(i).marker() << ";\n";
            }
            ofile << "];\n";

            // Writing the faces
            ofile << "mesh.f=[\n";
            for(UInt i(1);i<=mesh.storedFaces();++i){
                ofile << mesh.face(i).point(1).id() << ",";
                ofile << mesh.face(i).point(2).id() << ",";
                ofile << mesh.face(i).point(3).id() << ";\n";
            }
            ofile << "];\n";

            ofile << "mesh.fflag=[\n";
            for(UInt i(1);i<=mesh.storedFaces();++i){
            ofile << mesh.face(i).marker() << ";\n";
            }
            ofile << "];\n";

            // Writing the volumes
            ofile << "mesh.v=[\n";
            for(UInt i(1);i<=mesh.storedVolumes();++i){
                ofile << mesh.volume(i).point(1).id() << ",";
                ofile << mesh.volume(i).point(2).id() << ",";
                ofile << mesh.volume(i).point(3).id() << ",";
                ofile << mesh.volume(i).point(4).id() << ";\n";
            }
            ofile << "];\n";

            ofile << "mesh.vflag=[\n";
            for(UInt i(1);i<=mesh.storedVolumes();++i){
                ofile << mesh.volume(i).marker() << ";\n";
            }
            ofile << "];\n";
            ofile.close();
        }
        else std::cout << "Unable to open file" << std::endl;
    }else if((format&MESH_FORMAT)!=0){
        //Export to a mesh file
        outName.append(".mesh");
        std::ofstream ofile (outName.c_str());
        if (ofile.is_open())
        {
            ofile<< "MeshVersionFormatted 1\n\n";
            ofile<< "Dimension 3\n\n";

            // Writing the points
            ofile<< "Vertices\n";
            ofile<< mesh.storedPoints() << std::endl;
            for(UInt i(1);i<=mesh.storedPoints();++i){
                ofile << mesh.point(i).x() << " ";
                ofile << mesh.point(i).y() << " ";
                ofile << mesh.point(i).z() << " ";
                ofile << mesh.point(i).marker() << std::endl;
            }

            // Writing the volumes
            ofile<< "\nTetrahedra\n";
            ofile<< mesh.storedVolumes() << std::endl;
            for(UInt i(1);i<=mesh.storedVolumes();++i){
                ofile << mesh.volume(i).point(1).id() << " ";
                ofile << mesh.volume(i).point(2).id() << " ";
                ofile << mesh.volume(i).point(3).id() << " ";
                ofile << mesh.volume(i).point(4).id() << " ";
                ofile << mesh.volume(i).marker() << std::endl;
            }
            ofile.close();
        }
        else std::cout << "Unable to open file" << std::endl;
    }
}


//@}



} // Namespace LifeV

#endif /* EXPORTMESH3D_H */
