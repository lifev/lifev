/********************************************************************************
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
********************************************************************************/

/**
 * @file   ParserGmsh.hpp
 * @brief  Gmsh files (.gmsh, .msh) reader.
 * @author Simone Pezzuto <simone.pezzuto@mail.polimi.it>
 * @date   10/2012
**/

#ifndef PARSER_GMSH_HPP__
#define PARSER_GMSH_HPP__

#include <lifev/core/LifeV.hpp>
#include <lifev/core/mesh/BareMesh.hpp>

#include <fstream>
#include <deque>

namespace LifeV
{

namespace MeshIO
{

namespace
{
// Gmsh doesn't specify anything on int type
typedef int gmsh_int_t;
typedef double gmsh_float_t;

typedef struct
{
    gmsh_int_t id;
    short type;
    std::vector<gmsh_int_t> tags;
    std::vector<gmsh_int_t> nodes;
} gmsh_elm_t;

static const LifeV::UInt elm_nodes_num[] =
{
    2,  3,  4,  4,  8, 6,  5,  3,   6,  9, 10, 27, 18, 14,
    1,  8, 20, 15, 13, 9, 10, 12,  15, 15, 21,  4,  5,  6,
    20, 35, 56,  0,  0, 0,  0,  0,   0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0, 0,  0,  0,   0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0, 0,  0,  0,   0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0, 0,  0,  0,   0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0, 0,  0, 64, 125
};

namespace utils
{
namespace detail
{
// int to type (see std::integral_constant)
template <int id>
struct shape_id
{
    static const int value = id;
};

// type container
template <typename T>
struct shape_type
{
    typedef T type;
};
}

// <shape, int> map
// the friend function is dummy, but it's useful
// because if shape is associated to more than one
// int (or viceversa), an error is raised
template <typename shape, int id>
struct register_shape
        : detail::shape_id<id>, detail::shape_type<shape>
{
private:
    friend detail::shape_type<shape> registered (detail::shape_id<id>)
    {
        return detail::shape_type<shape>();
    }
};

// this is the register
template <typename geoshape>
struct id_of;

// shortcut macro
#       define GMSH_REGISTER_SHAPE(shape, id) \
        template <> struct id_of<shape> : register_shape<shape, id> {}

// gmsh shapes (same for legacy format)
GMSH_REGISTER_SHAPE (LifeV::nullShape, 0);
GMSH_REGISTER_SHAPE (LifeV::LinearLine, 1);
GMSH_REGISTER_SHAPE (LifeV::LinearTriangle, 2);
GMSH_REGISTER_SHAPE (LifeV::LinearQuad, 3);
GMSH_REGISTER_SHAPE (LifeV::LinearTetra, 4);
GMSH_REGISTER_SHAPE (LifeV::LinearHexa, 5);
GMSH_REGISTER_SHAPE (LifeV::QuadraticLine, 8);
GMSH_REGISTER_SHAPE (LifeV::QuadraticTriangle, 9);
GMSH_REGISTER_SHAPE (LifeV::QuadraticQuad, 10);
GMSH_REGISTER_SHAPE (LifeV::QuadraticTetra, 11);
GMSH_REGISTER_SHAPE (LifeV::QuadraticHexa, 12);
GMSH_REGISTER_SHAPE (LifeV::GeoPoint, 15);

// Remove the macro
#       undef GMSH_REGISTER_SHAPE

// LifeV shapes
template <typename GeoShape>
struct adm_shapes
{
    typedef GeoShape elem_type;
    typedef typename  elem_type::GeoBShape facet_type;
    typedef typename facet_type::GeoBShape ridge_type;
    typedef typename ridge_type::GeoBShape  peak_type;
    static const int dim = GeoShape::S_nDimensions;
    enum
    {
        elem_id  = id_of<elem_type>::value,
        facet_id = id_of<facet_type>::value,
        ridge_id = id_of<ridge_type>::value,
        peak_id  = id_of<peak_type>::value
    };
};
} // anonymous
} // utils

/**
 * @fn     ReadGmshFile
 * @brief  Reads a .msh file (ASCII, binary or legacy).
 *
 * @param  filename, name of the file
 * @param  baremesh, a baremesh object to be filled
 * @param  regionFlag, the id of the mesh (default 0)
 *
 * @return true if everything is ok
**/
template <typename GeoShape>
bool ReadGmshFile (const std::string& filename,
                   LifeV::BareMesh<GeoShape>& baremesh,
                   LifeV::ID regionFlag = 0,
                   bool verbose = false)
{

    // Check if the file is valid
    // ==========================
    std::ifstream ifile (filename.c_str(), std::ios::in | std::ios::binary);
    if (ifile.fail() )
    {
        // Error: file not found
        std::cerr << "[ERROR:GmshIO] File '"
                  << filename << "' not found." << std::endl;
        return false;
    }

    // Clean baremesh object
    baremesh.clear();

    // ============
    // begin.HEADER
    // ============
    bool is_binary = false, is_legacy = false;
    bool is_differ_endian = false;
    {
        bool status_ok = false;
        std::string line;
        std::getline (ifile, line);
        if (line == "$MeshFormat")
        {
            // Parse the next line
            std::getline (ifile, line);
            std::stringstream iss (line);

            float version;
            int doublesize;
            iss >> version >> is_binary >> doublesize;
            status_ok = iss.eof();

            // Endianness
            if (is_binary)
            {
                // Check file endianness
                union
                {
                    gmsh_int_t i;
                    char c[sizeof (gmsh_int_t)];
                } one;
                ifile.read (one.c, sizeof (gmsh_int_t) );
                ifile.ignore (1);
                is_differ_endian = (one.i != 1);
                // Check if system is little-endian
                if (is_differ_endian)
                {
                    std::cerr << "[ERROR:GmshIO] Different endianness of system"
                              << " and file not supported" << std::endl;
                    return false;
                }
            }

            if (verbose)
                std::clog << "[INFO:GmshIO] Found gmsh header with "
                          << (is_binary ? "binary" : "ASCII")
                          << " data." << std::endl;
        }
        else if (line == "$NOD")
        {
            // Legacy format
            is_legacy = true;
            status_ok = true;

            if (verbose)
                std::clog << "[INFO:GmshIO] Found LEGACY gmsh header"
                          << std::endl;
        }

        // Status should be fine if everything has been found
        if (!status_ok)
        {
            std::cerr << "[ERROR:GmshIO] Error during parsing header"
                      << std::endl;
            return false;
        }
    }

    // ===========
    // begin.NODES
    // ===========
    {
        // Back to the initial position
        ifile.seekg (0, std::ios::beg);

        // Look for the node section
        std::string node_tag = (is_legacy) ? "$NOD" : "$Nodes";

        // Seeking to the correct position
        std::string line;
        while (std::getline (ifile, line) && line != node_tag) {};
        if (ifile.eof() )
        {
            std::cerr << "[ERROR:GmshIO] Nodes section not found" << std::endl;
            return false;
        }

        // Next line is the number of nodes
        std::getline (ifile, line);
        std::size_t num_nodes = atoi (line.c_str() );
        if (!num_nodes)
        {
            std::cerr << "[ERROR:GmshIO] Number of nodes not correct" << std::endl;
            return false;
        }

        // Reshape the points datastructure: 3 x nnodes
        baremesh.points.reshape (3, num_nodes);
        baremesh.pointMarkers.resize (num_nodes);

        std::stringstream iss;

        // Loading the nodes (same for legacy)
        for (std::size_t i = 0; i < num_nodes; ++i)
        {
            gmsh_int_t id;
            gmsh_float_t p[3];

            if (is_binary)
            {
                // Binary case
                ifile.read (reinterpret_cast<char*> (&id), sizeof (gmsh_int_t) );
                ifile.read (reinterpret_cast<char*> (p), 3 * sizeof (gmsh_float_t) );
            }
            else
            {
                // ASCII case (also if legacy)
                std::getline (ifile, line);
                iss.clear();
                iss.str (line.c_str() );
                // Parsing
                iss >> id >> p[0] >> p[1] >> p[2];
            }

            // Adds the point
            baremesh.points (0, id - 1) = p[0];
            baremesh.points (1, id - 1) = p[1];
            baremesh.points (2, id - 1) = p[2];

        }
        if (is_binary)
        {
            ifile.ignore (1);
        }

        // Check if next line is $EndNodes
        std::getline (ifile, line);
        if (line != "$EndNodes" && line != "$ENDNOD")
        {
            std::cerr << "[ERROR:GmshIO] Something wrong with Nodes section"
                      << std::endl;
            return false;
        }

        if (verbose)
            std::clog << "[INFO:GmshIO] Found " << num_nodes << " nodes."
                      << std::endl;
    }

    // ==============
    // begin.ELEMENTS
    // ==============
    {
        // A gmsh elements can be anything, i.e. volume, edge, face, point
        // with an arbitrary number of tags.
        std::string elm_tag = (is_legacy) ? "$ELM" : "$Elements";

        // Seeking to the correct position
        std::string line;
        while (std::getline (ifile, line) && line != elm_tag) {};
        if (ifile.eof() )
        {
            std::cerr << "[ERROR:GmshIO] Elements section not found" << std::endl;
            return false;
        }
        // Next line is the number of elements
        std::getline (ifile, line);
        std::size_t num_elms = atoi (line.c_str() );

        // Don't know how many elements of each type there are, so we cannot
        // resize the data structures yet.
        baremesh.nDimensions = utils::adm_shapes<GeoShape>::dim;

        // Admissible shapes ids inherited from GeoShape
        const short s_ids[4] =
        {
            utils::adm_shapes<GeoShape>::elem_id,
            utils::adm_shapes<GeoShape>::facet_id,
            utils::adm_shapes<GeoShape>::ridge_id,
            utils::adm_shapes<GeoShape>::peak_id
        };

        // Stores the (gmsh) elements nodes and the type
        std::deque<gmsh_elm_t> elems[4];


        // Loading the elements
        bool status_ok = false;
        if (is_binary)
        {
            // Binary format
            gmsh_int_t elm_count = 0;
            do
            {
                // header = { elm_type, num_elems, num_tags }
                gmsh_int_t header[3];
                ifile.read (reinterpret_cast<char*> (header), 3 * sizeof (gmsh_int_t) );
                // The type of the element, from 0 (elem) to 3 (peak)
                short s_type = static_cast<short> (std::find (s_ids, s_ids + 4, header[0]) - s_ids);
                // Is the shape admissible?
                if (s_type == 4)
                {
                    // Found a shape that cannot be stored
                    std::cerr << "[ERROR:GmshIO] Found elements "
                              << " of type " << header[0]
                              << " , which is not an allowed shape." << std::endl;
                    status_ok = false;
                    break;
                }
                // Nodes of element of this type
                gmsh_int_t nnodes = elm_nodes_num[header[0] - 1];
                // Next we have the elements { tags, nodes, id } * num_elems
                gmsh_int_t total_size = header[1] * (1 + header[2] + nnodes);
                std::vector<gmsh_int_t> e (total_size);
                // Read at once
                ifile.read (reinterpret_cast<char*> (&e[0]), total_size * sizeof (gmsh_int_t) );
                // Increment counter
                elm_count += header[1];
                // Fill the queue
                for (gmsh_int_t k = 0; k < header[1]; ++k)
                {
                    gmsh_elm_t tmp;
                    tmp.type = header[0];
                    // Fill the tags
                    tmp.tags.reserve (header[2]);
                    std::copy (e.begin(), e.begin() + header[2] + 1, std::back_inserter (tmp.tags) );
                    // Fill the nodes
                    tmp.nodes.reserve (nnodes);
                    std::copy (e.begin() + header[2] + 1, e.end(), std::back_inserter (tmp.nodes) );
                    // Save
                    elems[s_type].push_back (tmp);
                }
                status_ok = true;
            }
            while (elm_count < static_cast<gmsh_int_t> (num_elms) );
            // Skip return
            ifile.ignore (1);
        }
        else
        {
            // ASCII format
            typedef std::istream_iterator<gmsh_int_t> iss_it;
            std::stringstream iss;
            for (std::size_t i = 0; i < num_elms; ++i)
            {
                gmsh_int_t id, type, num_tags;

                // Gathering id, type and number of tags
                std::getline (ifile, line);
                iss.clear();
                iss.str (line);
                iss >> id >> type >> num_tags;
                if (is_legacy)
                {
                    num_tags = 2;
                }

                // The type of the element, from 0 (elem) to 3 (peak)
                short s_type = static_cast<short> (std::find (s_ids, s_ids + 4, type) - s_ids);
                // Is the shape admissible?
                if (s_type == 4)
                {
                    // Found a shape that cannot be stored
                    std::cerr << "[ERROR:GmshIO] Element " << id
                              << ", of type " << type
                              << ", has a not allowed shape." << std::endl;
                    status_ok = false;
                    break;
                }

                // Parsing nodes and tags
                gmsh_elm_t tmp;
                tmp.type = type;
                // Fill the tags
                tmp.tags.resize (num_tags);
                for (int k = 0; k < num_tags; ++k)
                {
                    iss >> tmp.tags[k];
                }
                // Fill the values
                int nnodes = elm_nodes_num[type - 1];
                tmp.nodes.reserve (nnodes);
                std::copy (iss_it (iss), iss_it(), std::back_inserter (tmp.nodes) );
                // Saving
                elems[s_type].push_back (tmp);
                status_ok = iss.eof();
            }
        }
        // Check if next line is $EndElements
        std::getline (ifile, line);
        if (!status_ok || (line != "$EndElements" && line != "$ENDELM") )
        {
            std::cerr << "[ERROR:GmshIO] Something wrong with Elements section." << std::endl;
            return false;
        }

        // Now we can update the mesh
        if (verbose)
            std::clog << "[INFO:GmshIO] Highest dimension should be "
                      << baremesh.nDimensions << "." << std::endl;
        if (!elems[0].size() )
        {
            std::cerr << "[ERROR:GmshIO] No "
                      << baremesh.nDimensions << "d elements found!" << std::endl;
            return false;
        }

        //  BareMesh data structures
        LifeV::ArraySimple<LifeV::UInt>* bare_elm_ptr[3] =
        {
            &baremesh.elements,
            &baremesh.facets,
            &baremesh.ridges
        };
        std::vector<LifeV::ID>* bare_mrk_ptr[4] =
        {
            &baremesh.elementMarkers,
            &baremesh.facetMarkers,
            &baremesh.ridgeMarkers,
            &baremesh.pointMarkers
        };

        for (LifeV::UInt s = 0; s < baremesh.nDimensions + 1; ++s)
        {
            // (n-s)-dimensional element
            LifeV::UInt num_s_elems = elems[s].size();
            LifeV::UInt num_s_nodes = elm_nodes_num[s_ids[s] - 1];
            // Allocating memory inside baremesh
            if (s < 3)
            {
                bare_elm_ptr[s]->reshape (num_s_nodes, num_s_elems);
            }
            bare_mrk_ptr[s]->resize (num_s_elems);

            if (verbose)
                std::clog << "[INFO:GmshIO] Found " << num_s_elems
                          << " marked " << (baremesh.nDimensions - s)
                          << "d elements." << std::endl;

            // Fill the values
            for (LifeV::UInt i = 0; i < num_s_elems; ++i)
            {
                // Marker (only the first one)
                bare_mrk_ptr[s]->at (i) = elems[s][i].tags[1];
                // Nodes
                if (s < 3)
                {
                    for (LifeV::UInt k = 0; k < num_s_nodes; ++k)
                    {
                        (*bare_elm_ptr[s]) (k, i) = elems[s][i].nodes[k] - 1;
                    }
                }
            }
        }

        // In general we have no information about boundary points and facets
        // So we assume that all are boundary facets, and let the mesh checker
        // repair it if necessary
        baremesh.numBoundaryPoints = 0;
        baremesh.numBoundaryFacets = baremesh.facets.numberOfColumns();

    }

    // ============
    // begin.OTHERS
    // ============
    baremesh.regionMarkerID = regionFlag;

    // Everything seems fine
    return true;
}


} // GmshIO

} // LifeV

#endif // PARSER_GMSH_HPP__
