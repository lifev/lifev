//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @file
 *  @brief File containing a class for handling spatial discretization.
 *
 *  \author M.A. Fernandez
 *  \date 01/2003
 *  \version 1.0
 *
 *  \version 1.3
 *  \date 06/2009
 *  \author Cristiano Malossi<cristiano.malossi@epfl.ch>
 */

#ifndef _DATAMESH_H_
#define _DATAMESH_H_

#include <string>
#include <ostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

#ifdef TWODIM
    #include <life/lifemesh/regionMesh2D.hpp>
    #include <life/lifefilters/readMesh2D.hpp>
#elif defined THREEDIM
    #include <life/lifemesh/regionMesh3D.hpp>
    #include <life/lifefilters/readMesh3D.hpp>
#endif

namespace LifeV {

//! DataTime - Class for handling temporal discretization.
/*!
 *  @author M.A. Fernandez, Cristiano Malossi
 *
 *  The class is a container for time information.
 */
template <typename Mesh>
class DataMesh
{
public:

    // Typedef
    typedef Mesh                                Mesh_Type;
    typedef boost::shared_ptr<Mesh_Type>        Mesh_ptrType;


    /** @name Constructors, Destructor
     */
    //@{

    //! Empty Constructor
    DataMesh();

    //! Constructor
    /*!
      \param section the section in the data file
    */
    DataMesh( const GetPot& dataFile, const std::string& section = "space_discretization" );

    DataMesh( const DataMesh& dataMesh );

    //! Virtual destructor
    virtual ~DataMesh() {};

    //@}


    /** @name Methods
     */
    //@{

    void setup( const GetPot& dataFile, const std::string& section );

    void readMesh( );

    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    /** @name Set Methods
     */
    //@{

    void  setMesh( const Mesh_ptrType& mesh ) { M_mesh = mesh; }

    //@}


    /** @name Get Methods
     */
    //@{

    const Mesh_ptrType   mesh()      const { return M_mesh; }

    const std::string    meshDir()   const { return M_mesh_dir; }

    const std::string    meshFile()  const { return M_mesh_file; }

    //@}

private:

    Mesh_ptrType    M_mesh;         // the mesh

    std::string     M_mesh_dir;     // mesh dir
    std::string     M_mesh_file;    // mesh file
    std::string     M_mesh_type;    // mesh type

    bool            M_verbose;
};



// ===================================================
// Constructors & Destructor
// ===================================================
template <typename Mesh>
DataMesh<Mesh>::DataMesh( ):
    M_mesh      ( new Mesh ),
    M_mesh_dir  ( "./" ),
    M_mesh_file ( "mesh.mesh" ),
    M_mesh_type ( ".mesh" ),
    M_verbose   ( false )
{}

template <typename Mesh>
DataMesh<Mesh>::
DataMesh( const GetPot& dataFile, const std::string& section ):
    M_mesh      ( new Mesh ),
    M_mesh_dir  ( dataFile( ( section + "/mesh_dir"  ).data(), "./" ) ),
    M_mesh_file ( dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" ) ),
    M_mesh_type ( dataFile( ( section + "/mesh_type" ).data(), ".mesh" ) ),
    M_verbose   ( dataFile( ( section + "/verbose" ).data(), false ) )
{
    readMesh();
}

template <typename Mesh>
DataMesh<Mesh>::DataMesh( const DataMesh& dataMesh ):
    M_mesh        ( dataMesh.M_mesh ),
    M_mesh_dir    ( dataMesh.M_mesh_dir ),
    M_mesh_file   ( dataMesh.M_mesh_file ),
    M_mesh_type   ( dataMesh.M_mesh_type ),
    M_verbose     ( dataMesh.M_verbose )
{
    M_mesh->updateElementEdges();
    M_mesh->updateElementFaces();
}



// ===================================================
// Methods
// ===================================================
template <typename Mesh>
void
DataMesh<Mesh>::setup( const GetPot& dataFile, const std::string& section )
{
    M_mesh_dir  = dataFile( ( section + "/mesh_dir"  ).data(), "./" );
    M_mesh_file = dataFile( ( section + "/mesh_file" ).data(), "mesh.mesh" );
    M_mesh_type = dataFile( ( section + "/mesh_type" ).data(), ".mesh" );
    M_verbose   = dataFile( ( section + "/verbose" ).data(), 0 );

    readMesh();
}

template <typename Mesh>
void
DataMesh<Mesh>::readMesh( )
{
    if ( M_verbose )
        std::cout << "\nBuilding mesh ... ";

#ifdef TWODIM

    if ( M_mesh_type == ".msh" )
        readFreeFemFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, verbatim );
    else
        ERROR_MSG( "Sorry, this mesh file can not be loaded" );

    //Update Edges
    M_mesh->updateElementEdges(true);

#elif defined( THREEDIM )

    if ( M_mesh_type == ".mesh" )
        readINRIAMeshFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, M_verbose );
    else if ( M_mesh_type == ".m++" )
        readMppFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, M_verbose );
    else if ( M_mesh_type == ".msh" )
        readGmshFile( *M_mesh, M_mesh_dir + M_mesh_file, 1 );
    else if ( M_mesh_type == ".vol" )
        readNetgenMesh( *M_mesh, M_mesh_dir + M_mesh_file, 1, M_verbose );
    else
        ERROR_MSG( "Sorry, this mesh file can not be loaded" );

    //Update Edges & Faces
    M_mesh->updateElementEdges( true, M_verbose );
    M_mesh->updateElementFaces( true, M_verbose );

#endif

    if ( M_verbose )
        std::cout << "ok.\n" << std::endl;
}

template <typename Mesh>
void DataMesh<Mesh>::showMe( std::ostream& output ) const
{
    output << "mesh_dir   = " << M_mesh_dir << std::endl;
    output << "mesh_file  = " << M_mesh_file << std::endl;
    output << "mesh_type  = " << M_mesh_type << std::endl;
}

} // namespace LifeV

#endif
