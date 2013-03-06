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
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with LifeV. If not, see <http://www.gnu.org/licenses/>.
*******************************************************************************
*/
//@HEADER
/*!
 * @file
 * @brief Import mesh data formats into LifeV mesh data structure.
 *
 *
 * @date 06-11-2004
 *
 * @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
 *
 *
 * @contributor Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 *
 * @mantainer Alessio Fumagalli <alessio.fumagalli@mail.polimi.it>
 *
 */

#ifndef _IMPORTER_H
#define _IMPORTER_H 1

#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

/*! @enum MeshFormat
  List of reading mesh format.
*/
enum MeshFormat
{
    MESHPP, /*!< Meshpp type mesh */
    INRIA,  /*!< INRIA type mesh */
    GMSH,   /*!< Gmsh type mesh */
    NETGEN, /*!< NetGen type mesh */
    FREEFEM /*!< FreeFem type mesh */
};

namespace detail
{

// Import function for 3D mesh
template<typename Elt>
void
import3D ( std::string const&  fileName,
           MeshFormat const&   format,
           RegionMesh<Elt>&  mesh,
           markerID_Type     regionFlag )
{
    // Select the right mesh format
    switch ( format )
    {
        case MESHPP:
            readMppFile ( mesh, fileName, regionFlag );
            break;

        case INRIA:
            readINRIAMeshFile ( mesh, fileName, regionFlag );
            break;

        case GMSH:
            readGmshFile ( mesh, fileName, regionFlag );
            break;

        case NETGEN:
            readNetgenMesh ( mesh, fileName, regionFlag );
            break;
        default:
        {
            std::ostringstream ostr;
            ostr << "Unsupported 2D file format";
            throw std::invalid_argument ( ostr.str() );
        }
    }
} // import

//Import function for 2D mesh
template<typename Elt>
void
import2D ( std::string const& fileName,
           MeshFormat const&  format,
           RegionMesh<Elt>& mesh,
           markerID_Type    regionFlag )
{
    // Select the right mesh format, only Gmsh allowed
    switch ( format )
    {
        case FREEFEM:
            readFreeFemFile ( mesh, fileName, regionFlag );
            break;
        default:
        {
            std::ostringstream ostr;
            ostr << "Unsupported 2D file format";
            throw std::invalid_argument ( ostr.str() );
        }
    }
} // import

} // Namespace detail

//! Importer General interface for read different types of mesh.
/*!
  @author Christophe Prud'homme <christophe.prudhomme@epfl.ch>

  Import different type of mesh data formats into Life mesh data structure.

*/
class Importer
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty constructor, use GMSH as default mesh format
    Importer() :
        M_fileName ( ),
        M_format   ( GMSH )
    {}

    //! Constructor with name and format
    /*!
      @param filename mesh filename to import
      @param format format of the file
    */
    Importer ( std::string const& fileName, MeshFormat const&  format ) :
        M_fileName ( fileName ),
        M_format   ( format )
    {}

    //! Copy constructor
    /*!
      @param import Importer object to be copied
    */
    Importer ( const Importer& importer ) :
        M_fileName ( importer.M_fileName ),
        M_format   ( importer.M_format )
    {}

    //@}

    //! @name Operators
    //@{

    //! Assign opertor overloading
    /*!
      @param import Importer object to be copied
    */
    Importer& operator= ( const Importer& importer );

    //@}

    //! @name Methods
    //@{

    //! Import mesh with tetrahedras
    /*!
      @param mesh mesh data structure to fill in
      @param regionFlag marker for the region to load
    */
    void import ( RegionMesh<LinearTetra>& mesh, markerID_Type regionFlag );


    //! Import mesh with linear hexahedras
    /*!
      @param mesh mesh data structure to fill in
      @param regionFlag marker for the region to load
    */
    void import ( RegionMesh<LinearHexa>& mesh, markerID_Type regionFlag );

    //! Import mesh with linear triangles
    /*!
      @param mesh mesh data structure to fill in
      @param regionFlag marker for the region to load
    */
    void import ( RegionMesh<LinearTriangle>& mesh, markerID_Type regionFlag );


    //! Import mesh with linear quadrangles
    /*!
      @param mesh mesh data structure to fill in
      @param regionFlag marker for the region to load
    */
    void import ( RegionMesh<LinearQuad>& mesh, markerID_Type regionFlag );

    //! Print attributes of the class
    /*!
      @param output Stream to put the output
    */
    void showMe ( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Set the file name
    /*!
      @param fileName of the mesh file
    */
    inline void setFileName ( std::string const& fileName )
    {
        M_fileName = fileName;
    }

    //! Set the format of the mesh file
    /*!
      @param format format of the mesh file
    */
    inline void setFormat ( MeshFormat const& format )
    {
        M_format = format;
    }

    //@}

private:

    //! Name of the file to import
    std::string M_fileName;

    //! Format of the file to import
    MeshFormat M_format;
};

} // Namespace LifeV

#endif /* _IMPORTER_H */
