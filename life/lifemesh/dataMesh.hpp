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
    @contributor Lucia Mirabella <lucia.mirabell@gmail.com>
    @maintainer Lucia Mirabella <lucia.mirabell@gmail.com>

 */

#ifndef DATAMESH_H
#define DATAMESH_H

#include <string>
#include <ostream>
#include <life/lifefilters/GetPot.hpp>
#include <life/lifecore/life.hpp>

#ifdef TWODIM
#include <life/lifemesh/regionMesh2D.hpp>
#include <life/lifefilters/ImporterMesh2D.hpp>
#else // THREEDIM
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/ImporterMesh3D.hpp>
#include <life/lifemesh/structuredMesh3D.hpp>
#endif

namespace LifeV
{

//! DataMesh - class for handling spatial discretization.
/*!
	@author M.A. Fernandez
	@author Cristiano Malossi

	The class is a container for mesh information.
 */

class DataMesh
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    DataMesh();

    //! Constructor
    /*!
      @param dataFile data file
      @param section the section in the data file
     */
    DataMesh( const GetPot& dataFile, const std::string& section = "space_discretization" );

    //! Copy constructor
    /*!
     */
    DataMesh( const DataMesh& dataMesh );

    //! Virtual destructor
    virtual ~DataMesh() {};

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

    //@}


    //! @name Get Methods
    //@{

    const std::string&   meshDir()   const { return M_meshDir; }
    const std::string&   meshFile()  const { return M_meshFile; }
    const std::string&   meshType()  const { return M_meshType; }

    const bool&          verbose()   const { return M_verbose; }

    //@}

private:

    std::string     M_meshDir;     //!< mesh directory
    std::string     M_meshFile;    //!< mesh file
    std::string     M_meshType;    //!< mesh type

    bool            M_verbose;		//!< verbose output?
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

#endif
