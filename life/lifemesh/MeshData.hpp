//@HEADER
/*
 *******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

 *******************************************************************************
 */
//@HEADER

/*!
    @file
    @brief   File containing a class for handling spatial discretization.

    @author M.A. Fernandez
    @date 01-2003
    @author Cristiano Malossi<cristiano.malossi@epfl.ch>
    @date 06-2009
    @author Gilles Fourestey
    @date 06-2010
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com> Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Lucia Mirabella <lucia.mirabell@gmail.com>

 */

#ifndef MESHDATA_H
#define MESHDATA_H

#include <string>
#include <ostream>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/LifeV.hpp>

#ifdef TWODIM
#include <life/lifemesh/RegionMesh2D.hpp>
#include <life/lifefilters/ImporterMesh2D.hpp>
#else // THREEDIM
#include <life/lifemesh/RegionMesh3D.hpp>
#include <life/lifefilters/ImporterMesh3D.hpp>
#include <life/lifemesh/RegionMesh3DStructured.hpp>
#endif

namespace LifeV
{

//! MeshData - class for handling spatial discretization.
/*!
	@author M.A. Fernandez
	@author Cristiano Malossi

	The class is a container for mesh information.
 */

class MeshData
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    MeshData();

    //! Constructor
    /*!
      @param dataFile data file
      @param section the section in the data file
     */
    MeshData( const GetPot& dataFile, const std::string& section = "space_discretization" );

    //! Copy constructor
    /*!
     */
    MeshData( const MeshData& meshData );

    //! Virtual destructor
    virtual ~MeshData() {};

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the members
    /*!
     @param dataFile data file
     @param section file section
     */
    void setup( const GetPot& dataFile, const std::string& section );

    //! Display the values
    virtual void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    void setMeshDir  ( const std::string& dir )   { M_meshDir  = dir; }
    void setMeshFile ( const std::string& file )  { M_meshFile = file; }
    void setMeshType ( const std::string& type )  { M_meshType = type; }
    void setMOrder   ( const std::string& order ) { M_order    = order; }
    void setVerbose  ( const bool& isVerbose )    { M_verbose  = isVerbose; }

    //@}


    //! @name Get Methods
    //@{

    const std::string&   meshDir()   const { return M_meshDir; }
    const std::string&   meshFile()  const { return M_meshFile; }
    const std::string&   meshType()  const { return M_meshType; }
    const std::string&   mOrder()    const { return M_order; }
    const bool&          verbose()   const { return M_verbose; }

    //@}

private:

    std::string     M_meshDir;     //!< mesh directory
    std::string     M_meshFile;    //!< mesh file
    std::string     M_meshType;    //!< mesh type
    std::string     M_order;       //!< mesh type

    bool            M_verbose;		//!< verbose output?
};

template <typename Mesh>
void readMesh( Mesh& mesh, const MeshData& data )
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
        std::cout << "mesh read.\n" << std::endl;
}

} // namespace LifeV

#endif  /* MESHDATA_H */
