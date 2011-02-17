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


#include <life/lifefilters/ImporterMesh3D.hpp>
#include <life/lifefilters/ImporterMesh2D.hpp>
#include <life/lifefilters/Importer.hpp>


namespace LifeV
{

// ===================================================
// Operators
// ===================================================

// Assign opertor overloading
Importer& Importer::operator= ( const Importer& importer )
{
    // Avoid auto-copy
    if ( this != &importer )
    {
        M_fileName = importer.M_fileName;
        M_format   = importer.M_format;
    }

    return *this;
}

// ===================================================
// Methods
// ===================================================

// Import mesh with tetrahedras
void
Importer::import( RegionMesh3D<LinearTetra> & mesh,
                  entityFlag_Type             regionFlag )
{
    detail::import( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear hexahedras
void
Importer::import( RegionMesh3D<LinearHexa> & mesh,
                  entityFlag_Type            regionFlag )
{
    detail::import( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear triangles
void
Importer::import( RegionMesh2D<LinearTriangle> & mesh,
                  entityFlag_Type                regionFlag )
{
    detail::import( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear quadrangles
void
Importer::import( RegionMesh2D<LinearQuad> & mesh,
                  entityFlag_Type            regionFlag )
{
    detail::import( M_fileName, M_format, mesh, regionFlag );
} // import

// Print attributes of the class
void
Importer::showMe( std::ostream& output ) const
{
    output << "Class importer" << std::endl
           << "File Name   " << M_fileName << std::endl
           << "File format " << M_format << std::endl
           << std::flush;
} // showMe

} // Namespace LifeV
