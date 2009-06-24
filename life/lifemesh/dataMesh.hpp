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
/*!
  \file dataMesh.hpp

  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \version 1.3
  \date 06/2009
  \author Cristiano Malossi<cristiano.malossi@epfl.ch>

  \brief File containing a class for handling spatial discretization.
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

namespace LifeV
{

/*!
 * \class DataMesh
 * \brief DataMeshBase class which holds data concerning spatial discretization.
 *
 * @author M.A. Fernandez, Cristiano Malossi
 * @see
 */
template <typename Mesh>
class DataMesh
{
public:

	// Typedef
    typedef Mesh								mesh_raw_type;
    typedef boost::shared_ptr<mesh_raw_type> 	mesh_ptr_type;



    /** @name Constructors, Destructor
     */
	//@{

    //! Constructor, Destructor
    /*!
      \param section the section in the data file
    */
    DataMesh( const GetPot& dfile, const std::string& section = "discretization" );

    DataMesh( const DataMesh& dataMesh );

    //! Virtual destructor
    virtual ~DataMesh() {};

    //@}



    /** @name Get Functions
     */
	//@{

    const mesh_ptr_type mesh()		const { return M_mesh; }

    const std::string	meshDir()	const { return M_mesh_dir; }

    const std::string	meshFile()	const { return M_mesh_file; }

    //@}



    /** @name Set Functions
     */
	//@{

    void  setMesh( const mesh_ptr_type& mesh ) { M_mesh = mesh; }

    //@}



    /** @name General Functions
     */
	//@{

    virtual void showMe( std::ostream& c = std::cout ) const;

    //@}

private:

    mesh_ptr_type	M_mesh;			// the mesh

    std::string 	M_mesh_dir;		// mesh dir
    std::string 	M_mesh_file;	// mesh file
    std::string 	M_mesh_type;	// mesh type

};





//
// IMPLEMENTATION
//

// Constructor
template <typename Mesh>
DataMesh<Mesh>::
DataMesh( const GetPot& dfile, const std::string& section ):
    M_mesh		( new Mesh ),
    M_mesh_dir  ( dfile( ( section + "/mesh_dir"  ).data(), "./" ) ),
    M_mesh_file ( dfile( ( section + "/mesh_file" ).data(), "mesh.mesh" ) ),
    M_mesh_type ( dfile( ( section + "/mesh_type" ).data(), ".mesh" ) )
{
    bool verbose = dfile( ( section + "/verbose" ).data(), 0 );

	#ifdef TWODIM
		if ( M_mesh_type == ".msh" )
			readFreeFemFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, verbose );
		else
			ERROR_MSG( "Sorry, this mesh file can not be loaded" );

		//Update Edges
		M_mesh->updateElementEdges(true);



	#elif defined(THREEDIM)
		if ( M_mesh_type == ".mesh" )
			readINRIAMeshFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, verbose );
		else if ( M_mesh_type == ".m++" )
			readMppFile( *M_mesh, M_mesh_dir + M_mesh_file, 1, verbose );
		else if ( M_mesh_type == ".msh" )
			readGmshFile( *M_mesh, M_mesh_dir + M_mesh_file, 1 );
		else if ( M_mesh_type == ".vol" )
			readNetgenMesh( *M_mesh, M_mesh_dir + M_mesh_file, 1, verbose );
		else
			ERROR_MSG( "Sorry, this mesh file can not be loaded" );

		//Update Edges & Faces
		M_mesh->updateElementEdges( true );
		M_mesh->updateElementFaces( true );
	#endif
}



template <typename Mesh>
DataMesh<Mesh>::DataMesh( const DataMesh& dataMesh ):
	M_mesh        (dataMesh.M_mesh),
    M_mesh_dir    (dataMesh.M_mesh_dir),
    M_mesh_file   (dataMesh.M_mesh_file),
    M_mesh_type   (dataMesh.M_mesh_type)
{
	M_mesh->updateElementEdges();
	M_mesh->updateElementFaces();
}


//Output
template <typename Mesh>
void DataMesh<Mesh>::showMe( std::ostream& c ) const
{
    c << "mesh_dir   = " << M_mesh_dir << std::endl;
    c << "mesh_file  = " << M_mesh_file << std::endl;
    c << "mesh_type  = " << M_mesh_type << std::endl;
}

} // namespace LifeV

#endif
