/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#ifndef _MEDIT_WRTRS_H
#define _MEDIT_WRTRS_H
#include <fstream>

#include <life/lifecore/life.hpp>
#include <life/lifemesh/basisElSh.hpp>

namespace LifeV
{

/*!
  \file medit_wrtrs.h
  \author J.-F. Gerbeau & V. Martin
  \date 04/11/2002 - 27/02/2003

  \note Miguel 11/2002: writer for mesh

  \brief Write a mesh (Inria medit format ("toto.mesh")),
  and the solution (Inria medit format ("OUTtoto.bb").

  Also read a solution given in the Inria medit format.
  ("INtoto.bb")
*/

/*!
  A simple medit scalar writer.

  \param fname : name of the solution file.
  \warning medit assumes that, if the mesh file is name.mesh,
  the solution file is name.bb

  \param U : the array containing the solution

  \param Usize : the size of U

  \param type : 1 -> solution given per element
  2 -> solution given per vertex (default)

  \par Example of call:
  ScalUnknown<Vector> pressure(nbDof);
  wr_medit_ascii_scalar("meshname.bb",pressure.giveVec(),pressure.size());
*/

    void wr_medit_ascii_scalar( std::string fname, Real const* U, int Usize, int type = 2 );

/*!
  A simple medit vector writer.

  \param fname : name of the solution file.
  \warning medit assumes that, if the mesh file is name.mesh,
  the solution file is name.bb

  \param U : the array containing the solution

  \param Usize : the size of U

  \param type : 1 -> solution given per element
  2 -> solution given per vertex (default)

  \par Example of call:
  PhysVectUnknown<Vector> velocity(nbDof);
  wr_medit_ascii_vector("meshname.bb",velocity.giveVec(),velocity.size());
*/
    void wr_medit_ascii_vector( std::string fname, Real const* U, int Usize, int type = 2 );

/*!
  A simple medit scalar reader.
  It does NOT change the size of the given vector!

  \param fname : name of the solution file.
  \warning medit assumes that, if the mesh file is name.mesh,
  the solution file is name.bb

  \param U : the array containing the solution.
  The array must be correctly initialized previously.

  \param Usize : the size of U (correctly set.)

  \param type : 1 -> solution given per element
  2 -> solution given per vertex (default)

  \par Example of call:
  ScalUnknown<Vector> pressure(nbDof);
  wr_medit_ascii_scalar("meshname.bb",pressure.giveVec(),pressure.size());
*/

    void rd_medit_ascii_scalar( std::string fname, Real * U, const UInt& Usize, UInt& type );

/*!
  A simple medit vector reader.
  It does NOT change the size of the given vector!

  \param fname : name of the solution file.
  \warning medit assumes that, if the mesh file is name.mesh,
  the solution file is name.bb

  \param U : the array containing the solution
  The array must be correctly initialized previously.

  \param Usize : the size of U (correctly set.)

  \param type : 1 -> solution given per element
  2 -> solution given per vertex (default)

  \par Example of call:
  PhysVectUnknown<Vector> velocity(nbDof);
  wr_medit_ascii_vector("meshname.bb",velocity.giveVec(),velocity.size());
*/
    void rd_medit_ascii_vector( std::string fname, Real * U, const UInt& Usize, UInt& type );


/*!
  A simple medit mesh writer.

  \param fname : name of the mesh file.
  \param mesh : the Mesh object
*/

     template <typename Mesh>
     void writeMesh  (std::string fname, Mesh& mesh )
     {

         std::ofstream ofile( fname.c_str() );

         ASSERT( ofile, "Error: Output file cannot be open" );

         ofile << "MeshVersionFormatted 1\n";
         ofile << "Dimension\n3\n";
         ofile << "\n";
         ofile << "Vertices\n";

         UInt nV = mesh.numVertices();
         ofile << nV << "\n";

         int ig;

         std::map<int, int> localToGlobalNode = mesh.localToGlobalNode();
         std::map<int, int>::iterator im;


         for ( UInt i = 1; i <= nV; ++i )
         {
             ofile << mesh.pointList( i ).x() << " "
                   << mesh.pointList( i ).y() << " "
                   << mesh.pointList( i ).z() << " "
                   << mesh.pointList( i ).marker() << "\n";
         }
         ofile << "\n";

//         std::map<int, int>::iterator im;
//    UInt iface;

         typedef typename Mesh::FaceShape FaceShape;

         switch ( FaceShape::Shape )
         {
             case QUAD:
                 ofile << "Quadrilaterals\n";
                 break;
             case TRIANGLE:
                 ofile << "Triangles\n";
                 break;
             default:
                 ERROR_MSG( "BdShape not implement in MEDIT writer" );
         }

         UInt nBdF = mesh. numBFaces();
         ofile << nBdF << "\n";

         UInt nVpF = FaceShape::numVertices;


         for ( ID k = 1; k <= nBdF; ++k )
         {
             for ( ID i = 1; i <= nVpF; ++i )
             {
//            iface = mesh.boundaryFace( k ).point( i ).id();
                 ofile << mesh.boundaryFace( k ).point( i ).id()
                       << " ";
             }
             ofile << mesh.boundaryFace( k ).marker() << "\n";
         }
         ofile << "\n";

         typedef typename Mesh::VolumeShape ElementShape;

         switch ( ElementShape::Shape )
         {
             case HEXA:
                 ofile << "Hexaedra\n";
                 break;
             case TETRA:
                 ofile << "Tetrahedra\n";
                 break;
             default:
                 ERROR_MSG( "Shape not implement in MEDIT writer" );
         }

         UInt nE = mesh.numVolumes();
         ofile << nE << "\n";

         UInt nVpE = ElementShape::numVertices;

         UInt ielem;

         for ( ID k = 1; k <= nE; ++k )
         {
             for ( ID i = 1; i <= nVpE; ++i )
             {
                 ielem =  mesh.volume( k ).point( i ).localId();

                 im = localToGlobalNode.find(ielem);
                 ig = im->second;
//                  std::cout << ig << " <-> " << ielem << std::endl;

//                 ielem =  mesh.volume( k ).point( i ).localId();
                 ofile << ielem
                       << " ";
             }
             ofile << mesh.volume( k ).marker() << "\n";
         }
        ofile.close();
     }


/*!
  A simple medit mesh solutionb writer.

  \param fname : name of the mesh file.
  \param mesh  : the Mesh object
  \param sol   : solution,
  \param iStart: offset in the solution array.

*/




    template<typename Mesh, typename VectorType>
    void meditSolutionWriter( std::string        fname,
                              Mesh&              mesh,
                              const VectorType&  sol,
                              const int          iStart)
    {
        std::ofstream ofile( fname.c_str() );
        ASSERT( ofile, "Error: Output file cannot be opened." );

        int size = sol.size();
        int type = 2;

        ofile << nDimensions << " 1 " << size << " " << type << "\n";

        ofile.precision(12);

        UInt nVert       = mesh.numVertices();

        std::map<int, int> localToGlobalNode = mesh.localToGlobalNode();
        std::map<int, int>::iterator im;

        int ig;

        for(int i = 1; i <= (int)nVert; ++i)
        {
            int id = mesh.pointList( i ).id();

            im = localToGlobalNode.find(id);
            ig = im->second;

            ofile << sol[id + iStart] << "\n";
        }

        ofile.close();
    }







    template <typename Mesh>
    void wr_medit_ascii( std::string fname, const Mesh& mesh )
    {

        std::ofstream ofile( fname.c_str() );

        ASSERT( ofile, "Error: Output file cannot be open" );

        ofile << "MeshVersionFormatted 1\n";
        ofile << "Dimension 3\n";
        ofile << "\n";
        ofile << "Vertices\n";

        UInt nV = mesh.numVertices();
        ofile << nV << "\n";

        for ( UInt i = 1; i <= nV; ++i )
            ofile << mesh.pointList( i ).x() << " "
                  << mesh.pointList( i ).y() << " "
                  << mesh.pointList( i ).z() << " "
                  << mesh.pointList( i ).marker() << "\n";
        ofile << "\n";

        typedef typename Mesh::FaceShape FaceShape;

        switch ( FaceShape::Shape )
        {
            case QUAD:
                ofile << "Quadrilaterals\n";
                break;
            case TRIANGLE:
                ofile << "Triangles\n";
                break;
            default:
                ERROR_MSG( "BdShape not implement in MEDIT writer" );
        }

        UInt nBdF = mesh. numBFaces();
        ofile << nBdF << "\n";

        UInt nVpF = FaceShape::numVertices;


        for ( ID k = 1; k <= nBdF; ++k )
        {
            for ( ID i = 1; i <= nVpF; ++i )
                ofile << mesh.boundaryFace( k ).point( i ).id() << " ";
            ofile << mesh.boundaryFace( k ).marker() << "\n";
        }
        ofile << "\n";

        typedef typename Mesh::VolumeShape ElementShape;

        switch ( ElementShape::Shape )
        {
            case HEXA:
                ofile << "Hexaedra\n";
                break;
            case TETRA:
                ofile << "Tetrahedra\n";
                break;
            default:
                ERROR_MSG( "Shape not implement in MEDIT writer" );
        }

        UInt nE = mesh.numVolumes();
        ofile << nE << "\n";

        UInt nVpE = ElementShape::numVertices;

        for ( ID k = 1; k <= nE; ++k )
        {
            for ( ID i = 1; i <= nVpE; ++i )
                ofile << mesh.volume( k ).point( i ).id() << " ";
            ofile << mesh.volume( k ).marker() << "\n";
        }
        ofile.close();
    }


/*!
  A medit mesh writer.

  \param fname : name of the mesh file.
  \param mesh : the Mesh object
*/
    template <typename Mesh, typename Vector>
    void wr_medit_ascii( std::string fname, const Mesh& mesh, const Vector& disp, const Real& factor )
    {

        std::cout << "writing the mesh" << std::endl;
        std::ofstream ofile( fname.c_str() );

        ASSERT( ofile, "Error: Output file cannot be open" );

        ofile << "MeshVersionFormatted 1\n";
        ofile << "Dimension 3\n";
        ofile << "\n";
        ofile << "Vertices\n";

        UInt nV = mesh.numVertices();
        ofile << nV << "\n";

        for ( UInt i = 1; i <= nV; ++i )
        {
            int id = mesh.pointList(i).id();
            ofile << ( mesh.pointList( i ).x() - disp[ id  + 0*nV ] ) + factor * disp[ id + 0*nV ] << " "
                  << ( mesh.pointList( i ).y() - disp[ id  + 1*nV ] ) + factor * disp[ id + 1*nV ] << " "
                  << ( mesh.pointList( i ).z() - disp[ id  + 2*nV ] ) + factor * disp[ id + 2*nV ] << " "
                  << mesh.pointList( i ).marker() << "\n";
        }
        ofile << "\n";

        typedef typename Mesh::FaceShape FaceShape;

        switch ( FaceShape::Shape )
        {
            case QUAD:
                ofile << "Quadrilaterals\n";
                break;
            case TRIANGLE:
                ofile << "Triangles\n";
                break;
            default:
                ERROR_MSG( "BdShape not implement in MEDIT writer" );
        }

        UInt nBdF = mesh. numBFaces();
        ofile << nBdF << "\n";

        UInt nVpF = FaceShape::numVertices;


        for ( ID k = 1; k <= nBdF; ++k )
        {
            for ( ID i = 1; i <= nVpF; ++i )
                ofile << mesh.boundaryFace( k ).point( i ).id() << " ";
            ofile << mesh.boundaryFace( k ).marker() << "\n";
        }
        ofile << "\n";

        typedef typename Mesh::VolumeShape ElementShape;

        switch ( ElementShape::Shape )
        {
            case HEXA:
                ofile << "Hexaedra\n";
                break;
            case TETRA:
                ofile << "Tetrahedra\n";
                break;
            default:
                ERROR_MSG( "Shape not implement in MEDIT writer" );
        }

        UInt nE = mesh.numVolumes();
        ofile << nE << "\n";

        UInt nVpE = ElementShape::numVertices;

        for ( ID k = 1; k <= nE; ++k )
        {
            for ( ID i = 1; i <= nVpE; ++i )
                ofile << mesh.volume( k ).point( i ).id() << " ";
            ofile << mesh.volume( k ).marker() << "\n";
        }
        ofile.close();
    }


/*!
  A medit mesh writer.

  \param fname : name of the mesh file.
  \param mesh : the Mesh object
*/
    template <typename Mesh, typename Vector>
    void wr_medit_ascii2( std::string fname, const Mesh& mesh, const Vector& disp, const Real& factor )
    {

        std::ofstream ofile( fname.c_str() );

        ASSERT( ofile, "Error: Output file cannot be open" );

        ofile << "MeshVersionFormatted 1\n";
        ofile << "Dimension 3\n";
        ofile << "\n";
        ofile << "Vertices\n";

        UInt nVlocal  = mesh.numVertices();
        UInt nVglobal = mesh.numGlobalVertices();

        ofile << nVlocal << "\n";

        for ( UInt i = 1; i <= nVlocal; ++i )
        {
            int id = mesh.pointList(i).id() + 1;


            ofile << mesh.pointList( i ).x() + factor * disp[ id - 1 + 0 * nVglobal ] << " "
                  << mesh.pointList( i ).y() + factor * disp[ id - 1 + 1 * nVglobal ] << " "
                  << mesh.pointList( i ).z() + factor * disp[ id - 1 + 2 * nVglobal ] << " "
                  << mesh.pointList( i ).marker() << "\n";
        }

        ofile << "\n";

        typedef typename Mesh::FaceShape FaceShape;

        switch ( FaceShape::Shape )
        {
            case QUAD:
                ofile << "Quadrilaterals\n";
                break;
            case TRIANGLE:
                ofile << "Triangles\n";
                break;
            default:
                ERROR_MSG( "BdShape not implement in MEDIT writer" );
        }

        UInt nBdF = mesh.numBFaces();
        ofile << nBdF << "\n";

        UInt nVpF = FaceShape::numVertices;


        for ( ID k = 1; k <= nBdF; ++k )
        {
            for ( ID i = 1; i <= nVpF; ++i )
                ofile << mesh.boundaryFace( k ).point( i ).localId() << " ";
            ofile << mesh.boundaryFace( k ).marker() << "\n";
        }
        ofile << "\n";

        typedef typename Mesh::VolumeShape ElementShape;

        switch ( ElementShape::Shape )
        {
            case HEXA:
                ofile << "Hexaedra\n";
                break;
            case TETRA:
                ofile << "Tetrahedra\n";
                break;
            default:
                ERROR_MSG( "Shape not implement in MEDIT writer" );
        }

        UInt nE = mesh.numVolumes();
        ofile << nE << "\n";

        UInt nVpE = ElementShape::numVertices;

        for ( ID k = 1; k <= nE; ++k )
        {
            for ( ID i = 1; i <= nVpE; ++i )
                ofile << mesh.volume( k ).point( i ).localId() << " ";
            ofile << mesh.volume( k ).marker() << "\n";
        }
        ofile.close();
    }
}

#endif
