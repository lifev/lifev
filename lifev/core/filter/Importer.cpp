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


#include <lifev/core/filter/ImporterMesh3D.hpp>
#include <lifev/core/filter/ImporterMesh2D.hpp>
#include <lifev/core/filter/Importer.hpp>


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
Importer::import ( RegionMesh<LinearTetra>& mesh,
                   markerID_Type             regionFlag )
{
    detail::import3D ( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear hexahedras
void
Importer::import ( RegionMesh<LinearHexa>& mesh,
                   markerID_Type            regionFlag )
{
    detail::import3D ( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear triangles
void
Importer::import ( RegionMesh<LinearTriangle>& mesh,
                   markerID_Type                regionFlag )
{
    detail::import2D ( M_fileName, M_format, mesh, regionFlag );
} // import

// Import mesh with linear quadrangles
void
Importer::import ( RegionMesh<LinearQuad>&, markerID_Type )
{
    ERROR_MSG ("Importer:No importers available for this type of mesh");
} // import

// Print attributes of the class
void
Importer::showMe ( std::ostream& output ) const
{
    output << "Class importer" << std::endl
           << "File Name   " << M_fileName << std::endl
           << "File format " << M_format << std::endl
           << std::flush;
} // showMe

} // Namespace LifeV
