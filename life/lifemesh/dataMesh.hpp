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
#else // THREEDIM
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#endif

namespace LifeV
{

//! DataTime - Class for handling temporal discretization.
/*!
 *  @author M.A. Fernandez, Cristiano Malossi
 *
 *  The class is a container for time information.
 */
class DataMesh
{
public:

    //! @name Constructors, Destructor
    //@{

    //! Empty Constructor
    DataMesh();

    //! Constructor
    /*!
      @param section the section in the data file
    */
    DataMesh( const GetPot& dataFile, const std::string& section = "space_discretization" );

    DataMesh( const DataMesh& dataMesh );

    //! Virtual destructor
    virtual ~DataMesh() {};

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section file section
     */
    void setup( const GetPot& dataFile, const std::string& section );

    //! Display the values
    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    const std::string&   meshDir()   const { return M_mesh_dir; }
    const std::string&   meshFile()  const { return M_mesh_file; }
    const std::string&   meshType()  const { return M_mesh_type; }

    const bool&          verbose()   const { return M_verbose; }

    //@}

private:

    std::string     M_mesh_dir;     // mesh dir
    std::string     M_mesh_file;    // mesh file
    std::string     M_mesh_type;    // mesh type

    bool            M_verbose;
};

template <typename Mesh>
void readMesh( Mesh& mesh, const DataMesh& data )
{
    if ( data.verbose() )
        std::cout << "\nBuilding mesh ... ";

#ifdef TWODIM

    if ( data.meshType() == ".msh" )
        readFreeFemFile( mesh, data.meshDir() + M_meshFile(), 1, verbatim );
    else
        ERROR_MSG( "Sorry, this mesh file can not be loaded" );

    //Update Edges
    M_mesh->updateElementEdges(true);

#else // THREEDIM

    bool updateEdgesAndFaces(true);

    if ( data.meshType() == ".mesh" )
        readINRIAMeshFile( mesh, data.meshDir() + data.meshFile(), 1, data.verbose() );
    else if ( data.meshType() == ".m++" )
        readMppFile( mesh, data.meshDir() + data.meshFile(), 1, data.verbose() );
    else if ( data.meshType() == ".msh" )
        readGmshFile( mesh, data.meshDir() + data.meshFile(), 1 );
    else if ( data.meshType() == ".vol" )
        readNetgenMesh( mesh, data.meshDir() + data.meshFile(), 1, data.verbose() );
    // else if ( data.meshType() == "structured" )
    // {
    //     // Reading the parameters
    //     std::string section("fluid/space_discretization");
    //     UInt m_x( dataFile( ( section + "/mesh_mx"  ).data(), 3   ) );
    //     UInt m_y( dataFile( ( section + "/mesh_my"  ).data(), 3   ) );
    //     UInt m_z( dataFile( ( section + "/mesh_mz"  ).data(), 3   ) );
    //     Real l_x( dataFile( ( section + "/mesh_lx"  ).data(), 1.0 ) );
    //     Real l_y( dataFile( ( section + "/mesh_ly"  ).data(), 1.0 ) );
    //     Real l_z( dataFile( ( section + "/mesh_lz"  ).data(), 1.0 ) );
    //     Real t_x( dataFile( ( section + "/mesh_tx"  ).data(), 0.0 ) );
    //     Real t_y( dataFile( ( section + "/mesh_ty"  ).data(), 0.0 ) );
    //     Real t_z( dataFile( ( section + "/mesh_tz"  ).data(), 0.0 ) );
    //     regularMesh3D( *M_mesh, 1, m_x, m_y, m_z, true, l_x, l_y, l_z, t_x, t_y, t_z);
    // }
    // else if ( data.meshType() == "insideCode" )
    // {
    //     // Do nothing
    //     updateEdgesAndFaces = false;
    //}
    else
        ERROR_MSG( "Sorry, this mesh file can not be loaded" );

    //Update Edges & Faces
    if (updateEdgesAndFaces)
    {
        mesh.updateElementEdges( true, data.verbose() );
        mesh.updateElementFaces( true, data.verbose() );
    }

#endif

    if ( data.verbose() )
        std::cout << "ok.\n" << std::endl;
}

} // namespace LifeV

#endif
