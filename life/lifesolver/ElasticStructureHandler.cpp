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
 *  @file
 *  @brief This file contains an abstract class for Elastic Structures.
 *
 *  @version 1.0
 *  @date 01-06-2003
 *  @author Miguel Angel Fernandez
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#include <life/lifesolver/ElasticStructureHandler.hpp>

//==================================================
// IMPLEMENTATION
//=================================================



namespace LifeV
{

//=========================================
// Constructor
//=========================================
template <typename Mesh>
ElasticStructureHandler<Mesh>::
ElasticStructureHandler( const DataElasticStructure<Mesh>& data,
                         const RefFE&                      refFE,
                         const QuadRule&                   Qr,
                         const QuadRule&                   bdQr,
                         BCHandler&                        BCh ):
        DataElasticStructure<Mesh>( data ),
        M_refFE                    ( refFE ),
        M_dof                      ( this->mesh(), M_refFE ),
        M_dim                      ( M_dof.numTotalDof() ),
        M_Qr                       ( Qr ),
        M_bdQr                     ( bdQr ),
        M_fe                       ( M_refFE, getGeoMap( this->mesh() ), M_Qr ),
        M_feBd                     ( M_refFE.boundaryFE(), getGeoMap( this->mesh() ).boundaryMap(), M_bdQr ),
        M_d                        ( M_dim ),
        M_dRhs                     ( M_dim ),
        M_w                        ( M_dim ),
        M_time                     ( 0 ),
        M_count                    ( 0 ),
        M_BCh_solid               ( &BCh )
{}

template <typename Mesh>
ElasticStructureHandler<Mesh>::
ElasticStructureHandler( const DataElasticStructure<Mesh>& data,
                         const RefFE&                      refFE,
                         const QuadRule&                   Qr,
                         const QuadRule&                   bdQr):
        DataElasticStructure<Mesh>( data ),
        M_refFE                    ( refFE ),
        M_dof                      ( this->mesh(), M_refFE ),
        M_dim                      ( M_dof.numTotalDof() ),
        M_Qr                       ( Qr ),
        M_bdQr                     ( bdQr ),
        M_fe                       ( M_refFE, getGeoMap( this->mesh() ), M_Qr ),
        M_feBd                     ( M_refFE.boundaryFE(), getGeoMap( this->mesh() ).boundaryMap(), M_bdQr ),
        M_d                        ( M_dim ),
        M_dRhs                     ( M_dim ),
        M_w                        ( M_dim ),
//    _BCh( new BCHandler(0)),
        M_time                     ( 0 ),
        M_count                    ( 0 ),
        M_BCh_solid               ( 0 )
{}


// Returns the displacement vector
template <typename Mesh>
PhysVectUnknown<Vector> &
ElasticStructureHandler<Mesh>::disp()
{
    return M_d;
}


// Returns the velocity vector
template <typename Mesh>
PhysVectUnknown<Vector> &
ElasticStructureHandler<Mesh>::w()
{
    return M_w;
}


// Postprocessing
template <typename Mesh>
void
ElasticStructureHandler<Mesh>::postProcess()
{
    std::ostringstream index;
    std::string name, namedef;

    ++M_count;

    if ( fmod( float( M_count ), float( this->_verbose ) ) == 0.0 )
    {
        std::cout << "  S-  Post-processing \n";
        index << ( M_count / this->_verbose );

        switch ( index.str().size() )
        {
        case 1:
            name = "00" + index.str();
            break;
        case 2:
            name = "0" + index.str();
            break;
        case 3:
            name = index.str();
            break;
        }

        namedef = "defor." + name + ".mesh";
        wr_medit_ascii_scalar( "dep_x." + name + ".bb", M_d.giveVec(), this->mesh().numVertices() );
        wr_medit_ascii_scalar( "dep_y." + name + ".bb", M_d.giveVec() + M_dim, this->mesh().numVertices() );
        wr_medit_ascii_scalar( "dep_z." + name + ".bb", M_d.giveVec() + 2 * M_dim, this->mesh().numVertices() );
        wr_medit_ascii2( namedef, this->mesh(), M_d, this->_factor );

        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),this->mesh().numVertices(),_dim_u);

        system( ( "ln -sf " + namedef + " dep_x." + name + ".mesh" ).data() );
        system( ( "ln -sf " + namedef + " dep_y." + name + ".mesh" ).data() );
        system( ( "ln -sf " + namedef + " dep_z." + name + ".mesh" ).data() );

        // system(("ln -s "+this->mesh()_file+" veloc."+name+".mesh").data());

        wr_medit_ascii_scalar( "veld_x." + name + ".bb", M_w.giveVec(), this->mesh().numVertices() );
        wr_medit_ascii_scalar( "veld_y." + name + ".bb", M_w.giveVec() + M_dim, this->mesh().numVertices() );
        wr_medit_ascii_scalar( "veld_z." + name + ".bb", M_w.giveVec() + 2 * M_dim, this->mesh().numVertices() );

        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),this->mesh().numVertices(),_dim_u);

        system( ( "ln -sf " + namedef + " veld_x." + name + ".mesh" ).data() );
        system( ( "ln -sf " + namedef + " veld_y." + name + ".mesh" ).data() );
        system( ( "ln -sf " + namedef + " veld_z." + name + ".mesh" ).data() );

        // system(("ln -s "+this->mesh()_file+" veloc."+name+".mesh").data());
    }
}

// Sets the initial condition
template <typename Mesh>
void
ElasticStructureHandler<Mesh>::initialize( const Function& d0, const Function& w0 )
{


    // Initialize velocity

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = M_refFE.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = M_refFE.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = M_refFE.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = M_refFE.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    ID nbComp = M_d.nbcomp(); // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= this->mesh().numVolumes(); ++iElem )
    {

        M_fe.updateJac( this->mesh().volume( iElem ) );

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on element vertices
            for ( ID iVe = 1; iVe <= nElemV; ++iVe )
            {

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVe - 1 ) * nDofpV + l; // Local dof in this element

                    // Nodal coordinates
                    M_fe.coorMap( x, y, z, M_fe.refFE.xi( lDof - 1 ), M_fe.refFE.eta( lDof - 1 ), M_fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        M_w( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on element edges
            for ( ID iEd = 1; iEd <= nElemE; ++iEd )
            {

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    M_fe.coorMap( x, y, z, M_fe.refFE.xi( lDof - 1 ), M_fe.refFE.eta( lDof - 1 ), M_fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        M_w( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }

        // Face based Dof
        if ( nDofpF )
        {

            // loop on element faces
            for ( ID iFa = 1; iFa <= nElemF; ++iFa )
            {

                // Loop on number of Dof per face
                for ( ID l = 1; l <= nDofpF; ++l )
                {

                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l; // Local dof in the adjacent Element

                    // Nodal coordinates
                    M_fe.coorMap( x, y, z, M_fe.refFE.xi( lDof - 1 ), M_fe.refFE.eta( lDof - 1 ), M_fe.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {
                        M_d( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                        M_w( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
                    }
                }
            }
        }
        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            M_fe.coorMap( x, y, z, M_fe.refFE.xi( lDof - 1 ), M_fe.refFE.eta( lDof - 1 ), M_fe.refFE.zeta( lDof - 1 ) );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {
                M_d( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = d0( 0.0, x, y, z, icmp + 1 );
                M_w( icmp * M_dim + M_dof.localToGlobal( iElem, lDof ) - 1 ) = w0( 0.0, x, y, z, icmp + 1 );
            }
        }
    }
}


// Sets the initial condition
template <typename Mesh>
void
ElasticStructureHandler<Mesh>::initialize( const std::string& depName,
                                           const std::string& velName,
                                           double             startT)
{
    std::cout << "  S- restarting at time = " << startT << std::endl;

    M_count = (int) (startT/this->timestep() - 0.5);

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= this->mesh().numVolumes(); ++iElem )
    {
        M_fe.updateJac( this->mesh().volume( iElem ) );
    }
    readUnknown(depName, _d);
    readUnknown(velName, _w);
}


template <typename Mesh>
void
ElasticStructureHandler<Mesh>::readUnknown( const std::string       &name,
                                            PhysVectUnknown<Vector> &unknown)
{
    std::string sdummy;
    std::string ext;
    UInt nsx, nsy, nsz;
    int ndim;

    int nDof = M_dim;

    std::string filenamex = name;
    ext = "_x.bb";
    filenamex.insert(filenamex.end(), ext.begin(), ext.end());

    std::ifstream filex(filenamex.c_str(), std::ios::in);

    if (!filex)
    {
        std::cout << "Reading file " << filenamex
                  << " impossible" << std::endl;
        exit(1);
    }

    filex >> ndim;
    filex >> sdummy;
    filex >> nsx;

    if (nsx != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit(1);
    }

    filex >> sdummy;

    for (UInt ix = 0; ix < nsx; ++ix)
    {
        filex >> unknown[ix + 0*nDof];
//            unknown[ix + 0*nDof] /= factor;
    }

    filex.close();

    std::string filenamey = name;
    ext = "_y.bb";
    filenamey.insert(filenamey.end(), ext.begin(), ext.end());

    std::ifstream filey(filenamey.c_str(), std::ios::in);

    if (!filey)
    {
        std::cout << "Reading file " << filenamey
                  << " impossible" << std::endl;
        exit(1);
    }

    filey >> ndim;
    filey >> sdummy;
    filey >> nsy;

    if (nsy != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit(1);
    }

    filey >> sdummy;

    for (UInt iy = 0; iy < nsy; ++iy)
    {
        filey >> unknown[iy + 1*nDof];
//        unknown[iy + 1*nDof] /= factor;
    }

    filey.close();

    std::string filenamez = name;
    ext = "_z.bb";
    filenamez.insert(filenamez.end(), ext.begin(), ext.end());

//     std::cout << "Reading INRIA solid file   (" << filenamez << ")"
//               << ":" << std::endl;

    std::ifstream filez(filenamez.c_str(), std::ios::in);

    if (!filez)
    {
        std::cout << "Reading mesh file " << filenamez
                  << " impossible" << std::endl;
        exit(1);
    }

    filez >> ndim;
    filez >> sdummy;
    filez >> nsz;

    if (nsz != M_dim)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        exit(1);
    }

    filez >> sdummy;

    for (UInt iz = 0; iz < nsz; ++iz)
    {
        filez >> unknown[iz + 2*nDof];
//        unknown[iz + 2*nDof] /= factor;
    }

    filez.close();
}
}

