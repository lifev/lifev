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

#include <boost/function.hpp>


#include <life/lifecore/life.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifealg/dataAztec.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifefilters/gmv_wrtrs.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bdfNS.hpp>
#include <life/lifefem/postProc.hpp>
#include <life/lifefilters/openDX_wrtrs.hpp>
#include <cmath>
#include <sstream>
#include <ext/slist>
#include <life/lifearray/SimpleVect.hpp>
#include <utility>
using std::pair;

namespace LifeV
{

/*!
  \class NavierStokesHandler

  Abstract class which defines the general structure of a NavierStokes solver.
  For each new NavierStokes solver  we have to implement the corresponding
  timeAdvance and an iterate methods

*/

template <typename Mesh, typename DataType = DataNavierStokes<Mesh> >
class NavierStokesHandler:
            public DataType
{

public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;

    //! type used for flux computations (see FacesOnSections)
    typedef __gnu_cxx::slist< std::pair< ID, SimpleVect< ID > > > face_dof_type;

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

    //! Sets initial condition for the velocity and the pressure
    //! (incremental approach): the initial time is t0, the time step dt
    void initialize( const Function& u0, const Function& p0, Real t0, Real dt);

    //! Sets initial condition for the velocity and the pressure from file
    void initialize( const std::string & vname );

    //! Sets initial condition for the velocity and the pressure
    //! from medit file
    void initialize( const std::string& velName,
                     const std::string& pressName);
    //! set the source term functor
    void setSourceTerm( source_type __s )
        {
            _M_source = __s;
        }

    //! get the source term functor
    source_type sourceTerm() const
        {
            return _M_source;
        }

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    virtual void timeAdvance( source_type const& source, Real const& time ) =0;

    //! Update convective term, bc treatment and solve the linearized ns system
    virtual void iterate( const Real& time ) = 0;

    //! returns the mesh
    mesh_type& mesh() { return _mesh;}

    //! returns the BCHandler
    BCHandler& bcHandler()
        {
            return _BCh_u;
        }
    //! Returns the velocity vector
    PhysVectUnknown<Vector>& u();

    //! Returns the pressure
    ScalUnknown<Vector>& p();

    //! Returns the velocity Dof
    const Dof& uDof() const;

    //! Returns the pressure Dof
    const Dof& pDof() const;

    //! returns the FE ref for the velocity

    const RefFE& refFEu() const {return _refFE_u;}

    //! returns the FE ref for the pressure

    const RefFE& refFEp() const {return _refFE_p;}

    //! returns the current FE for the velocity u
    CurrentFE&   fe_u()   {return _fe_u;}
    CurrentBdFE& feBd_u() {return _feBd_u;}

    //! returns the current FE for the pressure p
    CurrentFE& fe_p() {return _fe_p;}

    //! Postprocessing
    void postProcess();


    //! OpenDX writers
    void dx_write_sol( std::string const& file_sol,
                       std::string const& fe_type_vel,
                       std::string const& fe_type_pre );

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
    //! NB: all the connectivity is recomputed at each call of this function
    Real flux( const EntityFlag& flag );

    /*! Computes the area and flux on a given part of the boundary
      \param flag the mesh flag identifying the part of the mesh where the flux is computed
      \param __reffe : reference fe.
      \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face). contains the local (in the face) to global dof.
    */
    void FacesOnFlag( const EntityFlag& flag, const RefFE& __reffe,
                      face_dof_type & __faces_on_flag );


    /*! Get the faces that live on a horizontal section (z=__z_section)
      up to a given tolerance (z = __z_section +/- __tolerance ).

      \param __z_section : position of the planar section
      \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face). contains the local (in the face) to global dof.

      The first is for velocity ("_u") and the second
      for pressure ("_p").
     \param __tolerance_section :
                  all faces between the two planes  z = __z_section - __tolerance_section
                  and z = __z_section + __tolerance_section are kept.
      (we only check the first 3 points of the faces)

      \param __x_frontier : given point on the boundary.
      \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
          close to a given boundary point __x_frontier.

      This function is VERY mesh dependent!! Use it with caution.
    */
    void FacesOnSection( const Real&  __z_section,
                         face_dof_type & __faces_on_section_u,
                         face_dof_type & __faces_on_section_p,
                         const Real& __tolerance_section,
                         const Real&  __x_frontier,
                         std::pair<ID, ID> & __point_on_boundary );

protected:
    /*! Update the local (face dof) to global dof vector

       \param __faces_on_section : list of faces that were found:
       first:  global face number (input)
       second: output vector of size nDofF (nb of dof per face).
       On output it contains the local (in the face) to global dof.
       \param __reffe : reference fe.
   */
    void _updateDofFaces( face_dof_type & __faces_on_section, const RefFE& __reffe );

public:
    /*! compute the area in the case of a cylindric domain (axis centered on the origin),
      given a boundary point (i.e. the radius).
      The area is approximated using the regular polygonal area formula.

      \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
      close to a given boundary point.

      \param nb_edge_polygon : nb of edges of the polygon for the formula
    */
    Real AreaCylindric( const std::pair<ID, ID> & __point_on_boundary,
                        const int & __nb_edge_polygon );

    /*! compute the Area and Flux
      \param __faces_on_section_u : list of (faces <-> veloc dof) that were found:
                                first:  global face number
                                second: vector of size nDofF (nb of dof per face)
                                        contains the local (in the face) to global dof.
      \param modify_sign_normal : if == true, the computed flux is positive if (u \dot z) is positive
    */
    std::pair< Real, Real > AreaAndFlux( const face_dof_type & __faces_on_section_u,
                                         const bool& modify_sign_normal = false );
    /*! compute the Area and Flux
      \param __faces_on_section_p : list of (faces <-> pressure dof) that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face)
      contains the local (in the face) to global dof.
    */
    Real MeanPressure( const face_dof_type & __faces_on_section_p );

    /*! compute the Areas and Fluxes for all sections and write them on a file
      in the plotmtv format
    */
    void PostProcessPressureAreaAndFlux( const Real & __time );

    //! Interpolate a given velocity function nodally onto a velocity vector
    void uInterpolate( const Function& uFct, Vector& uVect, Real time );

    //! calculate L2 pressure error for given exact pressure function
    //! takes into account a possible offset by a constant
    //! \param pexact the exact pressure as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in
    Real pErrorL2( const Function& pexact, Real time, Real* relError=0 );

    //! calculate L2 velocity error for given exact velocity function
    //! \param pexact the exact velocity as a function
    //! \param time the time
    //! \param relError Real* to store the relative error in
    Real uErrorL2( const Function& uexact, Real time, Real* relError=0 );

    //! Do nothing destructor
    virtual ~NavierStokesHandler()
    {}

protected:

    //! source term for NS
    source_type _M_source;

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

    //! The BDF Time Advance Method + Incremental Pressure
    BdfNS _bdf;


    //***** Prova di Agosto 2003
    PostProc<Mesh> _ns_post_proc;

    //! Aux. var. for PostProc
    UInt _count;

protected:
    //---------------
    //! stuff to compute the fluxes at each section for a cylindrical tube (mesh: tube20.mesh)
    //! This is hard coded as it strongly depends on the mesh
    //---------------
    const int M_nb_sections;
    std::vector<Real> M_z_section; //! position of the sections
    std::vector< face_dof_type > M_list_of_faces_on_section_velocity;
    std::vector< face_dof_type > M_list_of_faces_on_section_pressure;
    //! points that are on the external boundary (used to compute the area)
    std::vector< std::pair<ID, ID> > M_list_of_points_on_boundary;

    std::ofstream M_out_areas;
    std::ofstream M_out_areas_polygon;
    std::ofstream M_out_fluxes;
    std::ofstream M_out_pressure;

};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh, typename DataType>
NavierStokesHandler<Mesh, DataType>::
NavierStokesHandler( const GetPot& data_file, const RefFE& refFE_u,
                     const RefFE& refFE_p, const QuadRule& Qr_u,
                     const QuadRule& bdQr_u, const QuadRule& Qr_p,
                     const QuadRule& bdQr_p, BCHandler& BCh_u ) :
    DataType( data_file ),
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
    _feBd_u( _refFE_u.boundaryFE(), getGeoMap( this->_mesh ).boundaryMap(),
             _bdQr_u ),
    _fe_p( _refFE_p, getGeoMap( this->_mesh ), _Qr_p ),
    _u( _dim_u ),
    _p( _dim_p ),
    _BCh_u( BCh_u ),
    _bdf( this->_order_bdf ), _ns_post_proc( this->_mesh, _feBd_u, _dof_u,
                                             NDIM ),  // /******************
    _count( 0 ),
    //! stuff to compute the fluxes at each section (ex. of mesh: tube20.mesh)
    M_nb_sections( NbZSections() ),
    M_z_section( M_nb_sections ),
    M_list_of_faces_on_section_velocity( M_nb_sections ),
    M_list_of_faces_on_section_pressure( M_nb_sections ),
    M_list_of_points_on_boundary( M_nb_sections ),
    M_out_areas("Areas.res"), M_out_areas_polygon("AreasPolygon.res"),
    M_out_fluxes("Fluxes.res"), M_out_pressure("Pressure.res")
{
    if ( this->computeMeanValuesPerSection() == 1 ) {
        //---------------
        //! stuff to compute the fluxes at each section for a cylindrical tube
        //! (ex. of mesh: tube20.mesh)
        //---------------

        if ( ! this->_mesh.hasInternalFaces() )
            ERROR_MSG("The mesh must have all internal faces built up. Check that 'mesh_faces = all' in the data file.");
        if ( M_nb_sections < 2 )
            ERROR_MSG("We can't compute the mean values on less than 2 sections.");
        ASSERT( ZSectionFinal() - ZSectionInit() > 0,
                "We can't compute the mean values on less than 2 sections.");

        for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
            M_z_section[ izs ] = ZSectionInit() + Real(izs) * ( ZSectionFinal() - ZSectionInit() )
                / Real(M_nb_sections-1); // mesh dependent (length of tube=5)
        }

        for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
            //! for velocity
            this->FacesOnSection( M_z_section[izs],
                                  M_list_of_faces_on_section_velocity[izs],
                                  M_list_of_faces_on_section_pressure[izs],
                                  ToleranceSection(),
                                  XSectionFrontier(), M_list_of_points_on_boundary[izs] );
            if ( M_list_of_faces_on_section_velocity[izs].size() == 0 ||
                 M_list_of_faces_on_section_pressure[izs].size() == 0  ) {
                std::cout << "section z=" << M_z_section[izs] << " size="
                          << M_list_of_faces_on_section_velocity[izs].size() << std::endl;
                ERROR_MSG("For this section, no faces found.");
            }
        }
        if ( M_list_of_faces_on_section_velocity.size() == 0 ||
             M_list_of_faces_on_section_pressure.size() == 0 ) {
            ERROR_MSG("No list of faces found.");
        }
        //---------------
        //! end of stuff to compute the fluxes
        //---------------
    }

}


// Returns the velocity vector
template <typename Mesh, typename DataType>
PhysVectUnknown<Vector>&
NavierStokesHandler<Mesh, DataType>::u()
{
    return _u;
}

// Returns the pressure
template <typename Mesh, typename DataType>
ScalUnknown<Vector>&
NavierStokesHandler<Mesh, DataType>::p()
{
    return _p;
}

// Returns the velocity Dof
template <typename Mesh, typename DataType>
const Dof&
NavierStokesHandler<Mesh, DataType>::uDof() const
{
    return _dof_u;
}

// Returns the pressure Dof
template <typename Mesh, typename DataType>
const Dof&
NavierStokesHandler<Mesh, DataType>::pDof() const
{
    return _dof_p;
}

// Returns the BDF Time Advancing stuff
template <typename Mesh, typename DataType>
const BdfNS&
NavierStokesHandler<Mesh, DataType>::bdf() const
{
    return _bdf;
}

// Returns the Post Processing structure
template <typename Mesh, typename DataType>
PostProc<Mesh>&
NavierStokesHandler<Mesh, DataType>::post_proc()
{
    return _ns_post_proc;
}

// Set up of post processing structures
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::post_proc_set_area()
{
    _ns_post_proc.set_area( _feBd_u, _dof_u );
}

template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::post_proc_set_normal()
{
    _ns_post_proc.set_normal( _feBd_u, _dof_u );
}

template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::post_proc_set_phi()
{
    _ns_post_proc.set_phi( _feBd_u, _dof_u );
}

// Postprocessing
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::postProcess()
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
        wr_gmv_ascii( "test." + name + ".inp", this->_mesh, _dim_u,
                      _u.giveVec(), _p.giveVec() );

        // postprocess data file for medit
        wr_medit_ascii_scalar( "press." + name + ".bb", _p.giveVec(),
                               _p.size() );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", _u.giveVec(),
                               this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", _u.giveVec() + _dim_u,
                               this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", _u.giveVec()+2*_dim_u,
                               this->_mesh.numVertices() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file +
                  " press." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file +
                  " vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file +
                  " vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s " + this->_mesh_dir + this->_mesh_file +
                  " vel_z." + name + ".mesh" ).data() );
    }
}


// Writing (DX)
//! Write the solution in DX format
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::dx_write_sol( std::string const& file_sol,
                                         std::string const& fe_type_vel,
                                         std::string const& fe_type_pre )
{

    std::string file_vel = file_sol + "_vel.dx";
    std::string file_pre = file_sol + "_pre.dx";

    wr_opendx_header( file_pre, this->_mesh, _dof_p, _fe_p, fe_type_pre );
    wr_opendx_scalar( file_pre, "pression", _p );

    wr_opendx_header( file_vel, this->_mesh, _dof_u, _fe_u, fe_type_vel );
    wr_opendx_vector( file_vel, "velocity", _u, _u.nbcomp() );

}

// Set the initial condition
//! Initialize when only initial conditions on the velocity are given
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::initialize( const Function& u0 )
{

    // Initialize pressure
    _p = ZeroVector( _p.size() );

    // ********** initialize in the pressure BDF structure
    _bdf.bdf_p().initialize_unk( _p );

    // Initialize velocity
    uInterpolate( u0, _u, 0.0 );

    //****** Initialize in the BDF structure
    _bdf.bdf_u().initialize_unk( _u );

    _bdf.bdf_u().showMe();
    _bdf.bdf_p().showMe();

}

//! Initialize when  initial conditions both on the velocity and the pressure
//! (for incremental schemes) are given
//! Useful for test cases when the analytical solution is known
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::initialize( const Function& u0, const Function& p0,
                                       Real t0, Real dt )
{

    ID nbComp = _u.nbcomp(); // Number of components of the velocity


    _bdf.bdf_u().initialize_unk( u0, this->_mesh, _refFE_u, _fe_u, _dof_u, t0,
                                 dt, nbComp );

    // initialize _u with the first element in bdf_u.unk (=last value)
    _u = *( _bdf.bdf_u().unk().begin() );

    _bdf.bdf_p().initialize_unk( p0, this->_mesh, _refFE_p, _fe_p, _dof_p, t0,
                                 dt, 1 );

    // initialize _p with the first element in bdf_p.unk (=last value)
    _p = *( _bdf.bdf_p().unk().begin() );


    _bdf.bdf_u().showMe();
    _bdf.bdf_p().showMe();

}

//! Initialize when initial values for the velocity and the pressure are read
//! from file (M. Prosi)
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::initialize( const std::string & vname )
{

    std::fstream resfile( vname.c_str(), std::ios::in | std::ios::binary );
    if ( resfile.fail() )
    {
        std::cerr << " Error in initialization: File not found or locked"
                  << std::endl;
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

//! Initialize the fluid with a solution file
//! written in MEDIT format
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::initialize( const std::string& velName,
                                       const std::string& pressName)
{
    std::string sdummy;
    std::string ext;
    int nsx, nsy, nsz, nsp;
    int ndim;

    int nDof = uDof().numTotalDof();

    std::string filenamex = velName;

    ext = "_x.bb";
    filenamex.insert(filenamex.end(), ext.begin(), ext.end());

    std::cout << "Reading INRIA fluid file   (" << filenamex << ")"
              << ":" << std::endl;

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
    filex >> sdummy;

    for (int ix = 0; ix < nsx; ++ix)
        filex >> _u[ix + 0*nDof];

    filex.close();

    std::string filenamey = velName;
    ext = "_y.bb";
    filenamey.insert(filenamey.end(), ext.begin(), ext.end());

    std::cout << "Reading INRIA fluid file   (" << filenamey << ")"
              << ":" << std::endl;

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
    filey >> sdummy;

    for (int iy = 0; iy < nsy; ++iy)
        filey >> _u[iy + 1*nDof];

    filey.close();

    std::string filenamez = velName;

    ext = "_z.bb";
    filenamez.insert(filenamez.end(), ext.begin(), ext.end());

    std::cout << "Reading INRIA fluid file   (" << filenamez << ")"
              << ":" << std::endl;

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
    filez >> sdummy;

    for (int iz = 0; iz < nsz; ++iz)
        filez >> _u[iz + 2*nDof];

    filez.close();

    std::string filenamep = pressName;
    ext = ".bb";
    filenamep.insert(filenamep.end(), ext.begin(), ext.end());
    std::cout << "Reading INRIA fluid file   (" << filenamep << ")"
              << ":" << std::endl;

    std::ifstream filep(filenamep.c_str(), std::ios::in);

    if (!filep)
    {
        std::cout << "Reading file " << filenamep
                  << " impossible" << std::endl;
        exit(1);
    }

    filep >> ndim;
    filep >> sdummy;
    filep >> nsp;
    filep >> sdummy;

    for (int ip = 0; ip < nsp; ++ip)
        filep >> _p[ip];

    filep.close();
}

/*! Computes the flux on a given part of the boundary
    call of function AreaAndFlux( )
    NOT optimal! recompute the connectivity at each time step.
    (To do something faster, call once FacesOnFlag and store the
    connectivity... see ex. with FacesOnSection)
*/
template <typename Mesh, typename DataType>
Real
NavierStokesHandler<Mesh, DataType>::flux( const EntityFlag& flag )
{
    face_dof_type faces_on_flag;
    FacesOnFlag( flag, _refFE_u, faces_on_flag ); //! reconstruct all the connectivity

    std::pair< Real, Real > area_flux = AreaAndFlux( faces_on_flag,
                                                     false ); //! keep the orientation of the normal
    return area_flux.second;
}

/*! Get the faces that are on a part of the boundary which is defined by its flag

   \param flag : __flag of the faces that are sought
   \param __reffe : reference fe.
   \param __faces_on_section : list of faces that were found:
   first:  global face number
   second: vector of size nDofF (nb of dof per face).
   contains the local (in the face) to global dof.
*/
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::FacesOnFlag( const EntityFlag& __flag ,
                                        const RefFE& __reffe,
                                        face_dof_type & __faces_on_flag )
{

    typedef typename Mesh::FaceShape GeoBShape;

    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

    // Some useful local variables, to save some typing
    //! velocity
    // number of total Dof on a face
    UInt nDofF =
        __reffe.nbDofPerVertex * nFaceV
        + __reffe.nbDofPerEdge * nFaceE
        + __reffe.nbDofPerFace;

    UInt bdnF = this->_mesh.numBFaces();    // number of faces on boundary

    EntityFlag marker;
    typedef face_dof_type::iterator FaceDofIterator;

    //
    // Loop on boundary faces: List of boundary faces
    // with marker = __flag
    //
    for ( ID i = 1 ; i <= bdnF; ++i )
    {
        marker = this->_mesh.boundaryFace( i ).marker();
        if ( marker == __flag )
        {
            __faces_on_flag.push_front( make_pair( i, SimpleVect<ID>( nDofF ) ) );
        }
    }

    //! building the local to global vector for these boundary faces
    _updateDofFaces( __faces_on_flag, __reffe );

}

/*! Get the faces that live on a horizontal section (z=__z_section)
    up to a given tolerance (z = __z_section +/- __tolerance ).

    \param __z_section : position of the planar section
    \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face).
      contains the local (in the face) to global dof.

      The first is for velocity ("_u") and the second
      for pressure ("_p").
   \param __tolerance_section :
                  all faces between the two planes  z = __z_section - __tolerance_section
                  and z = __z_section + __tolerance_section are kept.
   (we only check the first 3 points of the faces)

   \param __x_frontier : given point on the boundary.
   \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
          close to a given boundary point __x_frontier.

   This function is VERY mesh dependent!! Use it with caution.
*/
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::FacesOnSection( const Real&  __z_section,
                                           face_dof_type & __faces_on_section_u,
                                           face_dof_type & __faces_on_section_p,
                                           const Real& __tolerance_section ,
                                           const Real&  __x_frontier,
                                           std::pair<ID, ID> & __point_on_boundary )
{
    typedef typename Mesh::FaceShape GeoBShape;

    //! geometrical data for faces
    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges
    UInt nF = this->_mesh.numFaces();    // number of faces (all faces)

    //! dof per face for u
    UInt nDofF_u =
        _refFE_u.nbDofPerVertex * nFaceV
        + _refFE_u.nbDofPerEdge * nFaceE
        + _refFE_u.nbDofPerFace;

    //! dof per face for p
    UInt nDofF_p =
        _refFE_p.nbDofPerVertex * nFaceV
        + _refFE_p.nbDofPerEdge * nFaceE
        + _refFE_p.nbDofPerFace;

    ID  iglobface;

    bool _found_point = false;

    //! search the faces on a given horizontal section
    //! loop on all faces
    for( ID iface = 1 ; iface <= nF ; iface ++ ) {
        iglobface = this->_mesh.faceList(iface).id();

        //! 3 points on the faces
        Real z_VeFa1 = this->_mesh.faceList(iface).point( 1 ).z();
        Real z_VeFa2 = this->_mesh.faceList(iface).point( 2 ).z();
        Real z_VeFa3 = this->_mesh.faceList(iface).point( 3 ).z();

        //! check the tolerance here (depends on the mesh)
        if ( std::fabs( __z_section - z_VeFa1 ) < __tolerance_section &&
             std::fabs( __z_section - z_VeFa2 ) < __tolerance_section &&
             std::fabs( __z_section - z_VeFa3 ) < __tolerance_section
             ) {
            __faces_on_section_u.push_front( make_pair( iglobface, SimpleVect<ID>( nDofF_u ) ) );
            __faces_on_section_p.push_front( make_pair( iglobface, SimpleVect<ID>( nDofF_p ) ) );

            ID iVeFa = 1;
            while ( ! _found_point && iVeFa <= nFaceV ) {
                Real x_VeFa = this->_mesh.faceList(iface).point( iVeFa ).x();
                if ( std::fabs( x_VeFa - __x_frontier ) < __tolerance_section ) {
                    __point_on_boundary = std::pair<ID, ID>( iglobface , iVeFa );
                    _found_point = true;
                }
                iVeFa ++;
            }
        }
    }

    //! building the local to global vector for these boundary faces (veloc and pressure)
    _updateDofFaces( __faces_on_section_u, _refFE_u );
    _updateDofFaces( __faces_on_section_p, _refFE_p );
}

/*! Update the local (face dof) to global dof vector

    \param __faces_on_section : list of faces that were found:
    first:  global face number (input)
    second: output vector of size nDofF (nb of dof per face).
    On output it contains the local (in the face) to global dof.
    \param __reffe : reference fe.
*/
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::_updateDofFaces( face_dof_type & __faces_on_section,
                                            const RefFE& __reffe )
{
    typedef typename Mesh::VolumeShape GeoShape;
    typedef typename Mesh::FaceShape GeoBShape;

    // Some useful local variables, to save some typing
    UInt nDofpV = __reffe.nbDofPerVertex; // number of Dof per vertices
    UInt nDofpE = __reffe.nbDofPerEdge;   // number of Dof per edges
    UInt nDofpF = __reffe.nbDofPerFace;   // number of Dof per faces

    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges

    UInt nDofFV = nDofpV * nFaceV; // number of vertex's Dof on a face
    UInt nDofFE = nDofpE * nFaceE; // number of edge's Dof on a face

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element

    //! search the faces on a given horizontal section
    ID iFace;
    UInt iElAd, iVeEl, iFaEl, iEdEl;
    ID lDof, gDof;

    typedef face_dof_type::iterator FaceDofIterator;


    //
    // Loop on faces: building the local to global vector
    // for these faces
    //
    for ( FaceDofIterator j = __faces_on_section.begin(); j != __faces_on_section.end(); ++j )
    {

        iFace = j->first;

        // id of the element adjacent to the face
        iElAd = this->_mesh.boundaryFace( iFace ).ad_first();

        // local id of the face in its adjacent element
        iFaEl = this->_mesh.boundaryFace( iFace ).pos_first();

        // Vertex based Dof
        if ( nDofpV )
        {

            // loop on face vertices
            for ( ID iVeFa = 1; iVeFa <= nFaceV; ++iVeFa )
            {
                // local vertex number (in element)
                iVeEl = GeoShape::fToP( iFaEl, iVeFa );

                // Loop number of Dof per vertex
                for ( ID l = 1; l <= nDofpV; ++l )
                {
                    // local Dof j-th degree of freedom on a face
                    lDof = ( iVeFa - 1 ) * nDofpV + l ;

                    // global Dof
                    gDof = _dof_u.localToGlobal( iElAd,
                                                 ( iVeEl - 1 ) * nDofpV + l );

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
                // local edge number (in element)
                iEdEl = GeoShape::fToE( iFaEl, iEdFa ).first;

                // Loop number of Dof per edge
                for ( ID l = 1; l <= nDofpE; ++l )
                {
                    // local Dof are put after the lDof of the vertices
                    lDof = nDofFV + ( iEdFa - 1 ) * nDofpE + l ;

                    // global Dof
                    gDof = _dof_u.localToGlobal( iElAd,
                                                 nDofElemV +
                                                 ( iEdEl - 1 ) * nDofpE + l );

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
                // local Dof are put after the lDof of the vertices and edges
                lDof = nDofFE + nDofFV + l;

                // global Dof
                gDof = _dof_u.localToGlobal( iElAd,
                                             nDofElemE + nDofElemV +
                                             ( iFaEl - 1 ) * nDofpF + l );

                j->second( lDof ) = gDof; // local to global on this face
            }
        }
    }
}

/*! compute the area in the case of a cylindric domain (axis centered on the origin),
    given a boundary point (i.e. the radius).
    The area is approximated using the regular polygonal area formula.

   \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
          close to a given boundary point.

   \param nb_edge_polygon : nb of edges of the polygon for the formula

 */
template <typename Mesh, typename DataType>
Real
NavierStokesHandler<Mesh, DataType>::AreaCylindric( const std::pair<ID, ID> & __point_on_boundary,
                                          const int & __nb_edge_polygon )
{
    ID iglobface = __point_on_boundary.first;
    ID iVeFa     = __point_on_boundary.second;

    Real radius = this->_mesh.faceList(iglobface).point( iVeFa ).x();

    return  0.5 * __nb_edge_polygon * radius * radius * std::sin( 2. * LifeV::Pi / __nb_edge_polygon );
}



/*! compute the Area and Flux
  \param __faces_on_section : list of (faces <-> veloc dof) that were found:
  first:  global face number
  second: vector of size nDofF (nb of dof per face)
  contains the local (in the face) to global dof.
  \param modify_sign_normal : if == true, the computed flux is positive if (u \dot z) is positive
*/
template <typename Mesh, typename DataType>
std::pair< Real, Real >
NavierStokesHandler<Mesh, DataType>::AreaAndFlux( const face_dof_type & __faces_on_section_u,
                                        const bool & modify_sign_normal )
{
    typedef typename Mesh::FaceShape GeoBShape;

    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

    // Some useful local variables, to save some typing
    //! velocity
    // number of total Dof on a face
    UInt nDofF =
        _refFE_u.nbDofPerVertex * nFaceV
        + _refFE_u.nbDofPerEdge * nFaceE
        + _refFE_u.nbDofPerFace;

    //! unknows
    Real __area = 0.0;
    Real __flux = 0.0;

    ID gDof, numTotalDof = _dof_u.numTotalDof();

    typedef face_dof_type::const_iterator constFaceDofIterator;

    // Number of velocity components
    UInt nc_u = _u.nbcomp();

    // Nodal values of the velocity in the current face
    std::vector<Real> u_local( nc_u * nDofF );

    // Loop on faces
    for ( constFaceDofIterator j = __faces_on_section_u.begin();
          j != __faces_on_section_u.end(); ++j )
    {

        // Extracting local values of the velocity in the current face
        for ( UInt ic = 0; ic < nc_u; ++ic )
        {
            for ( ID l = 1; l <= nDofF; ++l )
            {
                gDof = j->second( l );
                u_local[ ic * nDofF + l - 1 ] =
                    _u( ic * numTotalDof + gDof - 1 );
            }
        }

        // Updating quadrature data on the current face
        _feBd_u.updateMeasNormalQuadPt( this->_mesh.boundaryFace( j->first ) );

        __area += _feBd_u.measure();


        // Quadrature formula
        // Loop on quadrature points for velocities
        for ( int iq = 0; iq < _feBd_u.nbQuadPt; ++iq )
        {
            //! check the orientation of the normal
            //! the flux will be positive if (u \dot z) is positive
            Real sign_normal_dot_z = 1.0;
            if ( modify_sign_normal && _feBd_u.normal( 2, iq ) < 0 ) { //! component=2 -> z axis
                sign_normal_dot_z = -1.0;
            }

            // Dot product
            // Loop on components
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {

                // Interpolation
                // Loop on local dof
                for ( ID l = 1; l <= nDofF; ++l )
                    __flux += _feBd_u.weightMeas( iq ) *
                        u_local[ ic * nDofF + l - 1 ] *
                        _feBd_u.phi( int( l - 1 ), iq ) *
                        _feBd_u.normal( int( ic ), iq ) *
                        sign_normal_dot_z;
            }
        }
    }

    return  std::pair< Real, Real > ( __area , __flux );
}

/*! compute the Area and Flux
  \param __faces_on_section_p : list of (faces <-> pressure dof) that were found:
  first:  global face number
  second: vector of size nDofF (nb of dof per face)
  contains the local (in the face) to global dof.
*/
template <typename Mesh, typename DataType>
Real
NavierStokesHandler<Mesh, DataType>::MeanPressure( const face_dof_type & __faces_on_section_p )
{
    typedef typename Mesh::FaceShape GeoBShape;

    UInt nFaceV = GeoBShape::numVertices; // Number of face's vertices
    UInt nFaceE = GeoBShape::numEdges;    // Number of face's edges

    // Some useful local variables, to save some typing
    //! pressure
    // number of total Dof on a face
    UInt nDofF =
        _refFE_p.nbDofPerVertex * nFaceV
        + _refFE_p.nbDofPerEdge * nFaceE
        + _refFE_p.nbDofPerFace;


    //! unknows
    Real __pressure = 0.0;

    ID gDof;
    //ID numTotalDof = _dof_p.numTotalDof();

    typedef face_dof_type::const_iterator constFaceDofIterator;

    // Nodal values of the pressure in the current face
    std::vector<Real> p_local( nDofF );

    //! define the boundary fe for the pressure
    CurrentBdFE __ThefeBd_p( _refFE_p.boundaryFE(), getGeoMap( this->_mesh ).boundaryMap(),
                             _bdQr_p );

    // Loop on faces
    for ( constFaceDofIterator j = __faces_on_section_p.begin();
          j != __faces_on_section_p.end(); ++j )
    {

        for ( ID l = 1; l <= nDofF; ++l ){
            gDof = j->second( l );
            p_local[ l - 1 ] = _p( gDof - 1 );
        }

        // Updating quadrature data on the current face
        __ThefeBd_p.updateMeasQuadPt( this->_mesh.boundaryFace( j->first ) );

        // Quadrature formula
        // Loop on quadrature points for velocities
        for ( int iq = 0; iq < __ThefeBd_p.nbQuadPt; ++iq )
        {
            // Interpolation
            // Loop on local dof
            for ( ID l = 1; l <= nDofF; ++l )
                __pressure += //__ThefeBd_p.weightMeas( iq ) *
                    p_local[ l - 1 ] *
                    __ThefeBd_p.phi( int( l - 1 ), iq );
        }
    }

    return  __pressure;
}

/*! compute the mean Pressures, Areas and Fluxes for all sections and write them on a file
    in the plotmtv format
*/
template <typename Mesh, typename DataType>
void NavierStokesHandler<Mesh, DataType>::PostProcessPressureAreaAndFlux( const Real & __time )
{
    if ( computeMeanValuesPerSection() != 1 )
        ERROR_MSG("This function is disabled, if you don't ask to compute the Mean Values. (in data file)");

    M_out_areas  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                 << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Area'\n";
    M_out_fluxes << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                 << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Flux'\n";
    M_out_pressure  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                    << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Pressure'\n";
    M_out_areas_polygon  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                         << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='AreaPolygonal'\n";

    for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
        //! all normals are oriented along z axis -> "true"
        std::pair<Real, Real> AQmid  = this->AreaAndFlux( M_list_of_faces_on_section_velocity[izs], true );
        M_out_areas  << M_z_section[izs] << "\t" << AQmid.first << "\n";
        M_out_fluxes << M_z_section[izs] << "\t" << AQmid.second << "\n";
        M_out_pressure << M_z_section[izs] << "\t" << this->MeanPressure( M_list_of_faces_on_section_pressure[izs] ) << "\n";

        Real area_polygon = this->AreaCylindric( M_list_of_points_on_boundary[izs], NbPolygonEdges() );
        M_out_areas_polygon  << M_z_section[izs] << "\t" << area_polygon << "\n";
    }
    M_out_areas << std::endl;
    M_out_fluxes << std::endl;
    M_out_pressure << std::endl;
    M_out_areas_polygon << std::endl;
}


template <typename Mesh, typename DataType>
void NavierStokesHandler<Mesh, DataType>::uInterpolate( const Function& uFct,
                                                        Vector& uVect, Real time )
{
    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV  = _refFE_u.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE  = _refFE_u.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF  = _refFE_u.nbDofPerFace;   // number of Dof per face
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
                    // Local dof in this element
                    lDof = ( iVe - 1 ) * nDofpV + l;

                    // Nodal coordinates
                    _fe_u.coorMap( x, y, z,
                                   _fe_u.refFE.xi( lDof - 1 ),
                                   _fe_u.refFE.eta( lDof - 1 ),
                                   _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        uVect( icmp * _dim_u +
                               _dof_u.localToGlobal( iElem, lDof ) - 1 ) =
                            uFct( time, x, y, z, icmp + 1 );
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
                    // Local dof in the adjacent Element
                    lDof = nDofElemV + ( iEd - 1 ) * nDofpE + l;

                    // Nodal coordinates
                    _fe_u.coorMap( x, y, z,
                                   _fe_u.refFE.xi( lDof - 1 ),
                                   _fe_u.refFE.eta( lDof - 1 ),
                                   _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        uVect( icmp * _dim_u +
                               _dof_u.localToGlobal( iElem, lDof ) - 1 ) =
                            uFct( time, x, y, z, icmp + 1 );
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

                    // Local dof in the adjacent Element
                    lDof = nDofElemE + nDofElemV + ( iFa - 1 ) * nDofpF + l;

                    // Nodal coordinates
                    _fe_u.coorMap( x, y, z,
                                   _fe_u.refFE.xi( lDof - 1 ),
                                   _fe_u.refFE.eta( lDof - 1 ),
                                   _fe_u.refFE.zeta( lDof - 1 ) );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                        uVect( icmp * _dim_u +
                               _dof_u.localToGlobal( iElem, lDof ) - 1 ) =
                            uFct( time, x, y, z, icmp + 1 );
                }
            }
        }

        // Element based Dof
        // Loop on number of Dof per Element
        for ( ID l = 1; l <= nDofpEl; ++l )
        {
            // Local dof in the Element
            lDof = nDofElemF + nDofElemE + nDofElemV + l;

            // Nodal coordinates
            _fe_u.coorMap( x, y, z,
                           _fe_u.refFE.xi( lDof - 1 ),
                           _fe_u.refFE.eta( lDof - 1 ),
                           _fe_u.refFE.zeta( lDof - 1 ) );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                uVect( icmp * _dim_u +
                       _dof_u.localToGlobal( iElem, lDof ) - 1 ) =
                    uFct( time, x, y, z, icmp + 1 );
        }
    }
}

template <typename Mesh, typename DataType>
Real NavierStokesHandler<Mesh, DataType>::pErrorL2( const Function& pexact,
                                          Real time,
                                          Real* relError )
{
    Real sum2 = 0.;
    Real sum1 = 0.;
    Real sum0 = 0.;
    Real sumExact2 = 0.;
    Real sumExact1 = 0.;
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_p.updateFirstDeriv( _mesh.volumeList( iVol ) );
        sum2 += elem_L2_diff_2( _p, pexact, _fe_p, _dof_p, time, 1 );
        sum1 += elem_integral_diff( _p, pexact, _fe_p, _dof_p, time, 1 );
        sum0 += _fe_p.measure();
        if (relError)
        {
            sumExact2 += elem_L2_2( pexact, _fe_p, time, 1 );
            sumExact1 += elem_integral( pexact, _fe_p, time, 1 );
        }
    }
    Real absError = sqrt( sum2 - sum1*sum1/sum0 );
    if (relError)
    {
        Real normExact = sqrt( sumExact2 - sumExact1*sumExact1/sum0 );
        *relError = absError / normExact;
    }
    return absError;
}

template <typename Mesh, typename DataType>
Real NavierStokesHandler<Mesh, DataType>::uErrorL2( const Function& uexact,
                                          Real time,
                                          Real* relError )
{
    Real normU = 0.;
    UInt nbCompU = _u.nbcomp();
    Real sumExact = 0.;
    for ( UInt iVol = 1; iVol <= _mesh.numVolumes(); iVol++ )
    {
        _fe_u.updateFirstDeriv( _mesh.volumeList( iVol ) );
        normU += elem_L2_diff_2( _u, uexact, _fe_u, _dof_u, time,
                                 int( nbCompU ) );
        if (relError)
        {
            sumExact += elem_L2_2( uexact, _fe_u, time, int( nbCompU ) );
        }
    }
    if (relError)
    {
        *relError = sqrt( normU / sumExact );
    }
    return sqrt( normU );
}

}
#endif
