/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
  \file NavierStokesHandler.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \brief This file contains an abstract class for NavierStokes solvers.

*/

#ifndef _NAVIERSTOKESHANDLER_H_
#define _NAVIERSTOKESHANDLER_H_

#include <sstream>


#include "lifeV.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "geoMap.hpp"
#include "dataNavierStokes.hpp"
#include "dataAztec.hpp"
#include "medit_wrtrs.hpp"
#include "gmv_wrtrs.hpp"
#include "bcCond.hpp"
#include "bdfNS.hpp"
#include "post_proc.hpp"
#include "openDX_wrtrs.hpp"
#include <cmath>
#include <sstream>
#include <ext/slist>
#include "SimpleVect.hpp"
#include <utility>
using std::pair;

namespace LifeV
{

/*!
  \class NavierStokesHandler

  Abstract class which defines the general structure of a NavierStokes solver.
  For each new NavierStokes solver  we have to implement the corresponding timeAdvance and an iterate methods

*/

template <typename Mesh>
class NavierStokesHandler:
            public DataNavierStokes<Mesh>
{

public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p volumic quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_u boundary conditions for the velocity
      \param ord_bdf order of the bdf time advancing scheme and incremental pressure approach (default: Backward Euler)
    */
    NavierStokesHandler( const GetPot& data_file,
                         const RefFE& refFE_u,
                         const RefFE& refFE_p,
                         const QuadRule& Qr_u,
                         const QuadRule& bdQr_u,
                         const QuadRule& Qr_p,
                         const QuadRule& bdQr_p,
                         BCHandler& BCh_u );

    //! Sets initial condition for the velocity (here the initial time is 0.0)
    void initialize( const Function& u0 );

    //! Sets initial condition for the velocity and the pressure (incremental approach): the initial time is t0, the time step dt
    void initialize( const Function& u0, const Function& p0, Real t0, Real dt );

    //! Sets initial condition for the velocity and the pressure from file
    void initialize( const std::string & vname );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    virtual void timeAdvance( const Function source, const Real& time ) = 0;

    //! Update convective term, bc treatment and solve the linearized ns system
    virtual void iterate( const Real& time ) = 0;

    //! Returns the velocity vector
    PhysVectUnknown<Vector>& u();

    //! Returns the pressure
    ScalUnknown<Vector>& p();

    //! Returns the velocity Dof
    const Dof& uDof() const;

    //! Returns the pressure Dof
    const Dof& pDof() const;

    //! Postprocessing
    void postProcess();


    //! OpenDX writers
    void dx_write_sol( std::string const& file_sol, std::string const& fe_type_vel, std::string const& fe_type_pre );

    //! Returns the BDF Time Advancing stuff
    const BdfNS& bdf() const;

    //! Returns the  Post Processing  stuff
    PostProc<Mesh>& post_proc();

    //! Initialization of Post Processing structures
    void post_proc_set_area();
    void post_proc_set_normal();
    void post_proc_set_phi();

    //! Computes the flux on a given part of the boundary
    //! \param flag the mesh flag identifying the part of the mesh where the flux is computed
    Real flux( const EntityFlag& flag );

    //! calculate L2 pressure error for given exact pressure function
    //! takes into account a possible offset by a constant
    Real pErrorL2( const Function& pexact, Real time );

    //! calculate L2 velocity error for given exact velocity function
    Real uErrorL2( const Function& uexact, Real time );

    //! Do nothing destructor
    virtual ~NavierStokesHandler()
    {}


protected:

    //! Reference FE for the velocity
    const RefFE& _refFE_u;

    //! Reference FE for the pressure
    const RefFE& _refFE_p;

    //! The Dof object associated with the velocity
    Dof _dof_u;

    //! The Dof object associated with the pressure
    Dof _dof_p;

    //! The number of total velocity dofs
    UInt _dim_u;

    //! The number of total pressure dofs
    UInt _dim_p;

    //! Quadrature rule for velocity volumic elementary computations
    const QuadRule& _Qr_u;

    //! Quadrature rule for velocity surface elementary computations
    const QuadRule& _bdQr_u;

    //! Quadrature rule for pressure volumic elementary computations
    const QuadRule& _Qr_p;

    //! Quadrature rule for pressure surface elementary computations
    const QuadRule& _bdQr_p;

    //! Current FE for the velocity u
    CurrentFE _fe_u;
    CurrentBdFE _feBd_u;

    //! Current FE for the pressure p
    CurrentFE _fe_p;

    //! The velocity
    PhysVectUnknown<Vector> _u;

    //! The pressure
    ScalUnknown<Vector> _p;

    //! The BC handler
    BCHandler& _BCh_u;

    // ! The BDF Time Advance Method + Incremental Pressure
    BdfNS _bdf;


    //***** Prova di Agosto 2003
    PostProc<Mesh> _ns_post_proc;

    //! Aux. var. for PostProc
    UInt _count;
};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
NavierStokesHandler<Mesh>::
NavierStokesHandler( const GetPot& data_file, const RefFE& refFE_u,
                     const RefFE& refFE_p, const QuadRule& Qr_u, const QuadRule& bdQr_u,
                     const QuadRule& Qr_p, const QuadRule& bdQr_p, BCHandler& BCh_u ) :
    DataNavierStokes<Mesh>( data_file ),
    _refFE_u( refFE_u ),
    _refFE_p( refFE_p ),
    _dof_u( this->_mesh, _refFE_u ),
    _dof_p( this->_mesh, _refFE_p ),
    _dim_u( _dof_u.numTotalDof() ),
    _dim_p( _dof_p.numTotalDof() ),
    _Qr_u( Qr_u ),
    _bdQr_u( bdQr_u ),
    _Qr_p( Qr_p ),
    _bdQr_p( bdQr_p ),
    _fe_u( _refFE_u, getGeoMap( this->_mesh ), _Qr_u ),
    _feBd_u( _refFE_u.boundaryFE(), getGeoMap( this->_mesh ).boundaryMap(), _bdQr_u ),
    _fe_p( _refFE_p, getGeoMap( this->_mesh ), _Qr_p ),
    _u( _dim_u ),
    _p( _dim_p ),
    _BCh_u( BCh_u ),
    _bdf( this->_order_bdf ), _ns_post_proc( this->_mesh, _feBd_u, _dof_u, NDIM ),  // /******************
    _count( 0 )
{}


// Returns the velocity vector
template <typename Mesh>
PhysVectUnknown<Vector>&
NavierStokesHandler<Mesh>::u()
{
    return _u;
}

// Returns the pressure
template <typename Mesh>
ScalUnknown<Vector>&
NavierStokesHandler<Mesh>::p()
{
    return _p;
}

// Returns the velocity Dof
template <typename Mesh>
const Dof&
NavierStokesHandler<Mesh>::uDof() const
{
    return _dof_u;
}

// Returns the pressure Dof
template <typename Mesh>
const Dof&
NavierStokesHandler<Mesh>::pDof() const
{
    return _dof_p;
}

// Returns the BDF Time Advancing stuff
template <typename Mesh>
const BdfNS&
NavierStokesHandler<Mesh>::bdf() const
{
    return _bdf;
}

// Returns the Post Processing structure
template <typename Mesh>
PostProc<Mesh>&
NavierStokesHandler<Mesh>::post_proc()
{
    return _ns_post_proc;
}

// Set up of post processing structures
template <typename Mesh>
void
NavierStokesHandler<Mesh>::post_proc_set_area()
{
    _ns_post_proc.set_area( _feBd_u, _dof_u );
}

template <typename Mesh>
void
NavierStokesHandler<Mesh>::post_proc_set_normal()
{
    _ns_post_proc.set_normal( _feBd_u, _dof_u );
}

template <typename Mesh>
void
NavierStokesHandler<Mesh>::post_proc_set_phi()
{
    _ns_post_proc.set_phi( _feBd_u, _dof_u );
}

// Postprocessing
template <typename Mesh>
void
NavierStokesHandler<Mesh>::postProcess()
{
    std::ostringstream index;
    std::string name;

    ++_count;

    if ( fmod( float( _count ), float( this->_verbose ) ) == 0.0 )
    {
        std::cout << "  o-  Post-processing \n";
        index << ( _count / this->_verbose );

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

        // postprocess data file for GMV
        wr_gmv_ascii( "test." + name + ".inp", this->_mesh, _dim_u, _u.giveVec(), _p.giveVec() );

        // postprocess data file for medit
        wr_medit_ascii_scalar( "press." + name + ".bb", _p.giveVec(), _p.size() );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", _u.giveVec(), this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", _u.giveVec() + _dim_u, this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", _u.giveVec() + 2 * _dim_u, this->_mesh.numVertices() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file + " press." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file + " vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file + " vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file + " vel_z." + name + ".mesh" ).data() );
    }
}


// Writing (DX)
// ! Write the solution in DX format
template <typename Mesh>
void
NavierStokesHandler<Mesh>::dx_write_sol( std::string const& file_sol, std::string const& fe_type_vel, std::string const& fe_type_pre )
{

    std::string file_vel = file_sol + "_vel.dx";
    std::string file_pre = file_sol + "_pre.dx";

    wr_opendx_header( file_pre, this->_mesh, _dof_p, _fe_p, fe_type_pre );
    wr_opendx_scalar( file_pre, "pression", _p );

    wr_opendx_header( file_vel, this->_mesh, _dof_u, _fe_u, fe_type_vel );
    wr_opendx_vector( file_vel, "velocity", _u, _u.nbcomp() );

}

// Set the initial condition
// ! Initialize when only initial conditions on the velocity are given
template <typename Mesh>
void
NavierStokesHandler<Mesh>::initialize( const Function& u0 )
{

    // Initialize pressure
    _p = 0.0;

    // ********** initialize in the pressure BDF structure
    _bdf.bdf_p().initialize_unk( _p );

    // Initialize velocity

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = _refFE_u.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = _refFE_u.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = _refFE_u.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = _refFE_u.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element

    ID nbComp = _u.nbcomp(); // Number of components of the mesh velocity

    Real x, y, z;

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= this->_mesh.numVolumes(); ++iElem )
    {

        _fe_u.updateJac( this->_mesh.volume( iElem ) );

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
                    _fe_u.coorMap( x, y, z, _fe_u.refFE.xi( lDof - 1 ), _fe_u.refFE.eta( lDof - 1 ), _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        _u( icmp * _dim_u + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = u0( 0.0, x, y, z, icmp + 1 );
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
                    _fe_u.coorMap( x, y, z, _fe_u.refFE.xi( lDof - 1 ), _fe_u.refFE.eta( lDof - 1 ), _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        _u( icmp * _dim_u + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = u0( 0.0, x, y, z, icmp + 1 );
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
                    _fe_u.coorMap( x, y, z, _fe_u.refFE.xi( lDof - 1 ), _fe_u.refFE.eta( lDof - 1 ), _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        _u( icmp * _dim_u + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = u0( 0.0, x, y, z, icmp + 1 );
                }
            }
        }

        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            lDof = nDofElemF + nDofElemE + nDofElemV + l; // Local dof in the Element

            // Nodal coordinates
            _fe_u.coorMap( x, y, z, _fe_u.refFE.xi( lDof - 1 ), _fe_u.refFE.eta( lDof - 1 ), _fe_u.refFE.zeta( lDof - 1 ) );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                _u( icmp * _dim_u + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = u0( 0.0, x, y, z, icmp + 1 );
        }
    }
    //****** Initialize in the BDF structure
    _bdf.bdf_u().initialize_unk( _u );

    _bdf.bdf_u().showMe();
    _bdf.bdf_p().showMe();

}

// ! Initialize when  initial conditions both on the velocity and the pressure (for incremental schemes) are given
// ! Useful for test cases when the analytical solution is known
template <typename Mesh>
void
NavierStokesHandler<Mesh>::initialize( const Function& u0, const Function& p0, Real t0, Real dt )
{

    ID nbComp = _u.nbcomp(); // Number of components of the velocity


    _bdf.bdf_u().initialize_unk( u0, this->_mesh, _refFE_u, _fe_u, _dof_u, t0, dt, nbComp );
    _u = *( _bdf.bdf_u().unk().begin() ); // initialize _u with the first element in bdf_u.unk (=last value)

    _bdf.bdf_p().initialize_unk( p0, this->_mesh, _refFE_p, _fe_p, _dof_p, t0, dt, 1 );
    _p = *( _bdf.bdf_p().unk().begin() ); // initialize _u with the first element in bdf_u.unk (=last value)


    _bdf.bdf_u().showMe();
    _bdf.bdf_p().showMe();

}

// ! Initialize when initial values for the velocity and the pressure are read from file (M. Prosi)
template <typename Mesh>
void
NavierStokesHandler<Mesh>::initialize( const std::string & vname )
{

    std::fstream resfile( vname.c_str(), std::ios::in | std::ios::binary );
    if ( resfile.fail() )
    {
        std::cerr << " Error in initialization: File not found or locked" << std::endl;
        abort();
    }
    resfile.read( ( char* ) & _u( 1 ), _u.size() * sizeof( double ) );
    resfile.read( ( char* ) & _p( 1 ), _p.size() * sizeof( double ) );
    resfile.close();

    _bdf.bdf_u().initialize_unk( _u );
    _bdf.bdf_p().initialize_unk( _p );

    _bdf.bdf_u().showMe();
    _bdf.bdf_p().showMe();

}


//! Computes the flux on a given part of the boundary
template <typename Mesh>
Real
NavierStokesHandler<Mesh>::flux( const EntityFlag& flag )
{

    typedef typename Mesh::VolumeShape GeoShape;
    typedef typename Mesh::FaceShape GeoBShape;

    // Some useful local variables, to save some typing
    UInt nDofpV = _refFE_u.nbDofPerVertex; // number of Dof per vertices
    UInt nDofpE = _refFE_u.nbDofPerEdge;   // number of Dof per edges
    UInt nDofpF = _refFE_u.nbDofPerFace;   // number of Dof per faces

    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges

    UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
    UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

    UInt nDofF = nDofFV + nDofFE + nDofpF; // number of total Dof on a face

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element

    UInt bdnF = this->_mesh.numBFaces();    // number of faces on boundary

    Real flux = 0.0;
    __gnu_cxx::slist<std::pair<ID, SimpleVect<ID> > > faces;
    ID ibF;
    UInt iElAd, iVeEl, iFaEl, iEdEl;
    ID lDof, gDof, numTotalDof = _dof_u.numTotalDof();

    EntityFlag marker;
    typedef __gnu_cxx::slist<pair<ID, SimpleVect<ID> > >::iterator Iterator;

    //
    // Loop on boundary faces: List of boundary faces
    // with marker = flag
    //
    for ( ID i = 1 ; i <= bdnF; ++i )
    {
        marker = this->_mesh.boundaryFace( i ).marker();
        if ( marker == flag )
        {
            faces.push_front( make_pair( i, SimpleVect<ID>( nDofF ) ) );
        }
    }

    //
    // Loop on faces: building the local to global vector
    // for these boundary faces
    //
    for ( Iterator j = faces.begin(); j != faces.end(); ++j )
    {

        ibF = j->first;

        iElAd = this->_mesh.boundaryFace( ibF ).ad_first();  // id of the element adjacent to the face
        iFaEl = this->_mesh.boundaryFace( ibF ).pos_first(); // local id of the face in its adjacent element

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on face vertices
            for ( ID iVeFa = 1; iVeFa <= nFaceV; ++iVeFa )
            {

                iVeEl = GeoShape::fToP( iFaEl, iVeFa ); // local vertex number (in element)

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    lDof = ( iVeFa - 1 ) * nDofpV + l ; // local Dof j-esimo grado di liberta' su una faccia
                    gDof = _dof_u.localToGlobal( iElAd, ( iVeEl - 1 ) * nDofpV + l ); // global Dof
                    j->second( lDof ) = gDof; // local to global on this face
                }
            }
        }

        // Edge based Dof
        if ( nDofpE )
        {

            // loop on face edges
            for ( ID iEdFa = 1; iEdFa <= nFaceE; ++iEdFa )
            {

                iEdEl = GeoShape::fToE( iFaEl, iEdFa ).first; // local edge number (in element)
                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {

                    lDof = nDofFV + ( iEdFa - 1 ) * nDofpE + l ; // local Dof sono messi dopo gli lDof dei vertici
                    gDof = _dof_u.localToGlobal( iElAd, nDofElemV + ( iEdEl - 1 ) * nDofpE + l ); // global Dof
                    j->second( lDof ) = gDof; // local to global on this face
                }
            }
        }
        // Face based Dof
        if ( nDofpF )
        {

            // Loop on number of Dof per face
            for ( ID l = 1; l <= nDofpF; ++l )
            {
                lDof = nDofFE + nDofFV + l; // local Dof sono messi dopo gli lDof dei vertici e dopo quelli degli spigoli
                gDof = _dof_u.localToGlobal( iElAd, nDofElemE + nDofElemV + ( iFaEl - 1 ) * nDofpF + l ); // global Dof
                j->second( lDof ) = gDof; // local to global on this face
            }
        }
    }

    // Number of velocity components
    UInt nc_u = _u.nbcomp();

    // Nodal values of the velocity in the current face
    std::vector<Real> u_local( nc_u * nDofF );

    // Loop on faces
    for ( Iterator j = faces.begin(); j != faces.end(); ++j )
    {

        // Extracting lodal values of the velocity in the current face
        for ( UInt ic = 0; ic < nc_u; ++ic )
        {
            for ( ID l = 1; l <= nDofF; ++l )
            {
                gDof = j->second( l );
                u_local[ ic * nDofF + l - 1 ] = _u( ic * numTotalDof + gDof - 1 );
            }
        }

        // Updating quadrature data on the current face
        _feBd_u.updateMeasNormalQuadPt( this->_mesh.boundaryFace( j->first ) );

        // Quadrautre formula
        // Loop on quadrature points
        for ( int iq = 0; iq < _feBd_u.nbQuadPt; ++iq )
        {

            // Dot product
            // Loop on components
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {

                // Interpolation
                // Loop on local dof
                for ( ID l = 1; l <= nDofF; ++l )
                    flux += _feBd_u.weightMeas( iq ) * u_local[ ic * nDofF + l - 1 ] * _feBd_u.phi( int( l - 1 ), iq ) * _feBd_u.normal( int( ic ), iq );
            }
        }
    }

    return flux;
}

template<typename Mesh>
Real NavierStokesHandler<Mesh>::pErrorL2( const Function& pexact, Real time )
{
    Real sum2 = 0.;
    Real sum1 = 0.;
    Real sum0 = 0.;
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_p.updateFirstDeriv( _mesh.volumeList( iVol ) );
        sum2 += elem_L2_diff_2( _p, pexact, _fe_p, _dof_p, time, 1 );
        sum1 += elem_integral_diff( _p, pexact, _fe_p, _dof_p, time, 1 );
        sum0 += _fe_p.measure();
    }
    return sqrt( sum2 - sum1*sum1/sum0 );
}

template<typename Mesh>
Real NavierStokesHandler<Mesh>::uErrorL2( const Function& uexact, Real time )
{
    Real normU = 0.;
    UInt nbCompU = _u.nbcomp();
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_u.updateFirstDeriv( _mesh.volumeList( iVol ) );
        normU += elem_L2_diff_2( _u, uexact, _fe_u, _dof_u, time,
                                 int( nbCompU ) );
    }
    return sqrt( normU );
}

}
#endif
