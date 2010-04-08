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
  \file NavierStokesHandler.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \brief This file contains an abstract class for NavierStokes solvers.

*/

#warning: NavierStokesHandler.hpp NavierStokesSolverPC.hpp  are obselete should be removed from lifev-parallel

#ifndef _NAVIERSTOKESHANDLER_H_
#define _NAVIERSTOKESHANDLER_H_

#include <sstream>

#include <boost/function.hpp>


#include <life/lifecore/life.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifefem/geoMap.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
//#include <life/lifealg/dataAztec.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifefilters/gmv_wrtrs.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bdfNS.hpp>
#include <life/lifefem/postProc.hpp>
#include <life/lifefilters/openDX_wrtrs.hpp>
#include <life/lifearray/pattern.hpp>
#include <cmath>
#include <sstream>
#include <iomanip>
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
class NavierStokesHandler//:
//            public DataType
{

public:

    typedef DataType data_type;

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&,
                                 const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&,
                                   Real const&, ID const& )> source_type;

    typedef Mesh mesh_type;
    typedef BCHandler                             bchandler_raw_type;
    typedef boost::shared_ptr<bchandler_raw_type> bchandler_type;

    //! type used for flux computations (see FacesOnSections)
    typedef __gnu_cxx::slist< std::pair< ID, SimpleVect< ID > > > face_dof_type;
    typedef face_dof_type::const_iterator const_face_dof_iterator_type;

//    typedef typename NSStabilization NSStabilization;
    //! Constructors
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the velocity
      \param refFE_p reference FE for the pressure
      \param Qr_u volumic quadrature rule for the velocity
      \param bdQr_u surface quadrature rule for the velocity
      \param Qr_p volumic quadrature rule for the pressure
      \param bdQr_p surface quadrature rule for the pressure
      \param BCh_fluid boundary conditions for the fluid
      \param ord_bdf order of the bdf time advancing scheme and incremental pressure approach (default: Backward Euler)
    */
    NavierStokesHandler( const GetPot&   data_file,
                         const RefFE&    refFE_u,
                         const RefFE&    refFE_p,
                         const QuadRule& Qr_u,
                         const QuadRule& bdQr_u,
                         const QuadRule& Qr_p,
                         const QuadRule& bdQr_p,
                         BCHandler&      BCh_u );

    NavierStokesHandler( const DataType&        dataNavierStokes,
                         const RefFE&           refFE_u,
                         const RefFE&           refFE_p,
                         const QuadRule&        Qr_u,
                         const QuadRule&        bdQr_u,
                         const QuadRule&        Qr_p,
                         const QuadRule&        bdQr_p,
                         BCHandler&             BCh_u );

    /*! constructor without BCs */
    NavierStokesHandler( const GetPot&   data_file,
                         const RefFE&    refFE_u,
                         const RefFE&    refFE_p,
                         const QuadRule& Qr_u,
                         const QuadRule& bdQr_u,
                         const QuadRule& Qr_p,
                         const QuadRule& bdQr_p);


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

    //! checking if BC are set
    const bool setFluidBC() const {return M_setBC;}
    //! set the fluid BCs
    void setFluidBC(BCHandler &BCh_u){M_BCh_fluid = &BCh_u; M_setBC = true;}
    //! returns the BCHandler
    //BCHandler& BCh_fluid() {return *M_BCh_fluid;}
    //! deprecated
    BCHandler& bcHandler() {return *M_BCh_fluid;}
    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    virtual void timeAdvance( source_type const& source, Real const& time ) =0;

    //! Update convective term, bc treatment and solve the linearized ns system
    virtual void iterate( const Real& time ) = 0;

    // DataStokes accessors
    //! returns the mesh
//@@    mesh_type& mesh() { return this->_mesh;}
    const DataType & dataType() const {return M_dataType;}

    mesh_type& mesh()             {return M_dataType.dataMesh()->mesh();}
    const mesh_type& mesh() const {return M_dataType.dataMesh()->mesh();}

    std::string meshDir() {return M_dataType.dataMesh()->meshDir();}
    std::string meshFile(){return M_dataType.dataMesh()->meshFile();}

    Real        dt()   {return M_dataType.dataTime()->getTimeStep();}
    Real  timestep()   {return M_dataType.dataTime()->getTimeStep();}
    Real  inittime()   {return M_dataType.dataTime()->getInitialTime();}
    Real   endtime()   {return M_dataType.dataTime()->getEndTime();}
    Real   density()   {return M_dataType.density();}
    Real viscosity()   {return M_dataType.viscosity();}

    unsigned int order_bdf() const {return M_dataType.dataTime()->getBDF_order();}

    void showMe()  {return M_dataType.showMe();}
    UInt verbose() {return M_dataType.verbose();}
    Real factor()  {return M_dataType.factor();}

    NSStabilization stabilization(){return M_dataType.stabilization();}

    UInt computeMeanValuesPerSection(){return M_dataType.computeMeanValuesPerSection();}

    UInt NbZSections()                {return M_dataType.NbZSections();}
    Real ToleranceSection()           {return M_dataType.ToleranceSection();}
    Real XSectionFrontier()           {return M_dataType.XSectionFrontier();}
    Real ZSectionInit()               {return M_dataType.ZSectionInit();}
    Real ZSectionFinal()              {return M_dataType.ZSectionFinal();}
    UInt NbPolygonEdges()             {return M_dataType.NbPolygonEdges();}

    //! Returns the velocity vector
    PhysVectUnknown<Vector>& u();

    //! Returns the pressure
    ScalUnknown<Vector>& p();

    //! Returns the velocity Dof
    const Dof& uDof() const;
    Dof& uDof();

    //! Returns the pressure Dof
    const Dof& pDof() const;
    Dof& pDof();

    //! returns the FE ref for the velocity
    const RefFE& refFEu() const {return _refFE_u;}

    //! returns the FE ref for the pressure
    const RefFE& refFEp() const {return _refFE_p;}

    //! returns the current FE for the velocity u
    CurrentFE&   fe_u()   {return _fe_u;}

    //! returns the current boundary FE for the velocity u
    CurrentBdFE& feBd_u() {return _feBd_u;}

    //! returns the current FE for the pressure p
    CurrentFE&   fe_p()   {return _fe_p;}

    //! returns the current boundary FE for the pressure p
    CurrentBdFE& feBd_p() {return _feBd_p;}

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
      \param __dof : dof.
      \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face). contains the local (in the face) to global dof.
    */
    void FacesOnFlag( const EntityFlag& flag, const RefFE& __reffe,
                      const Dof& __dof, face_dof_type & __faces_on_flag );


    /*! Get the faces that live on a horizontal section (z=__z_section)
      up to a given tolerance (z = __z_section +/- __tolerance ).

      \param __planar_section : position of the planar section (plane
      given by its 4 (THREED) components (a,b,c,d) such that
      a x + b y + c z + d = 0
      \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face). contains the local (in the face) to global dof.

      The first is for velocity ("_u") and the second
      for pressure ("_p").

      \param __tolerance_section :
      all faces between the two planes
      a x + b y + c z + d = +/- __tolerance_section
      are kept.
      (we only check the first 3 points of the faces)

      \param __x_frontier : given point on the boundary.
      \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
          close to a given boundary point __x_frontier.

      This function is VERY mesh dependent!! Use it with caution.
    */
    void FacesOnSection( const SimpleVect<Real,0> &  __planar_section,
                         face_dof_type & __faces_on_section_u,
                         face_dof_type & __faces_on_section_p,
                         const Real& __tolerance_section,
                         const Real&  __x_frontier,
                         std::pair<ID, ID> & __point_on_boundary );

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
     \param __modify_sign_normal : if == true, the computed flux is positive if (u \dot __reference_normal) is positive
     \param __reference_normal : vector provided to give a fixed orientation to the normal

    */
    std::pair< Real, Real > AreaAndFlux( const face_dof_type & __faces_on_section_u,
                                         const bool& __modify_sign_normal,
                                         const SimpleVect<Real,0>& __reference_normal);
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



    BasePattern::PatternType patternType();


    //! Do nothing destructor
    virtual ~NavierStokesHandler()
    {}

protected:

    //! data for NS solvers

    data_type     M_dataType;

    //! source term for NS
    source_type   _M_source;

    //! Reference FE for the velocity
    const RefFE&  _refFE_u;

    //! Reference FE for the pressure
    const RefFE&  _refFE_p;

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
    CurrentBdFE _feBd_p;

    //! The velocity
    PhysVectUnknown<Vector> _u;

    //! The pressure
    ScalUnknown<Vector> _p;



    //! The BDF Time Advance Method + Incremental Pressure
    BdfNS _bdf;


    //***** Prova di Agosto 2003
    PostProc<Mesh> _ns_post_proc;

    //! Aux. var. for PostProc
    UInt _count;

    //---------------
    //! stuff to compute the fluxes at each section for a cylindrical tube (mesh: tube20.mesh)
    //! This is hard coded as it strongly depends on the mesh
    //---------------
    const int M_nb_sections;
    std::vector<Real> M_z_section; //! position of the sections
    std::vector< face_dof_type > M_list_of_faces_on_section_velocity;
    std::vector< face_dof_type > M_list_of_faces_on_section_pressure;
    //! points that are on the external boundary (used to compute the area)
    std::vector< std::pair<ID, ID> >  M_list_of_points_on_boundary;

    std::ofstream M_out_areas;
#if COMPUTE_POLYGONAL_AREA
    std::ofstream M_out_areas_polygon;
#endif
    std::ofstream M_out_fluxes;
    std::ofstream M_out_pressure;

    /*! Update the local (face dof) to global dof vector

       \param __faces_on_section : list of faces that were found:
       first:  global face number (input)
       second: output vector of size nDofF (nb of dof per face).
       On output it contains the local (in the face) to global dof.
       \param __reffe : reference fe.
   */
    void _updateDofFaces( face_dof_type & __faces_on_section,
                          const RefFE& __reffe, const Dof& __dof );


private:

//! private members
    bool          M_setBC;

//! private methods
    void initializeMeanValuesPerSection();
    void initializeSectionsBifurc();

    //! The BC handler
    BCHandler    *M_BCh_fluid;

    BasePattern::PatternType M_patternType;


};



//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh, typename DataType>
NavierStokesHandler<Mesh, DataType>::
NavierStokesHandler( const GetPot& data_file, const RefFE& refFE_u,
                     const RefFE& refFE_p, const QuadRule& Qr_u,
                     const QuadRule& bdQr_u, const QuadRule& Qr_p,
                     const QuadRule& bdQr_p, BCHandler& BCh_u ) :
    M_dataType                           ( data_file ),
    _refFE_u                           ( refFE_u ),
    _refFE_p                           ( refFE_p ),
    _dof_u                             ( this->mesh(), _refFE_u ),
    _dof_p                             ( this->mesh(), _refFE_p ),
    _dim_u                             ( _dof_u.numTotalDof() ),
    _dim_p                             ( _dof_p.numTotalDof() ),
    _Qr_u                              ( Qr_u ),
    _bdQr_u                            ( bdQr_u ),
    _Qr_p                              ( Qr_p ),
    _bdQr_p                            ( bdQr_p ),
    _fe_u                              ( _refFE_u,
                                         getGeoMap( this->mesh() ),
                                         _Qr_u ),
    _feBd_u                            ( _refFE_u.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_u ),
    _fe_p                              ( _refFE_p,
                                         getGeoMap( this->mesh() ),
                                         _Qr_p ),
    _feBd_p                            ( _refFE_p.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_p ),
    _u                                 ( _dim_u ),
    _p                                 ( _dim_p ),
    _bdf                               ( this->order_bdf() ),
    _ns_post_proc                      ( this->mesh(), _feBd_u, _dof_u, NDIM ),
    _count                             ( 0 ),
    //! stuff to compute the fluxes at each section (ex. of mesh: tube20.mesh)
    M_nb_sections                      ( this->NbZSections() ),
    M_z_section                        ( M_nb_sections ),
    M_list_of_faces_on_section_velocity( M_nb_sections ),
    M_list_of_faces_on_section_pressure( M_nb_sections ),
    M_list_of_points_on_boundary       ( M_nb_sections ),
    M_out_areas                        ("Areas.res"),
#if COMPUTE_POLYGONAL_AREA
    M_out_areas_polygon                ("AreasPolygon.res"),
#endif
    M_out_fluxes                       ("Fluxes.res"),
    M_out_pressure                     ("Pressure.res"),
    M_BCh_fluid                        ( &BCh_u )
{
    std::cout << "Old fluid constructor ... " << std::endl;
    if ( this->computeMeanValuesPerSection() == 1 )
        initializeMeanValuesPerSection();
    /*
    //! for simple bifurcation mesh
        initializeSectionsBifurc();
    */
}


// Constructors
template <typename Mesh, typename DataType>
NavierStokesHandler<Mesh, DataType>::
NavierStokesHandler( const DataType&  dataNavierStokes,
                     const RefFE&     refFE_u,
                     const RefFE&     refFE_p,
                     const QuadRule&  Qr_u,
                     const QuadRule&  bdQr_u,
                     const QuadRule&  Qr_p,
                     const QuadRule&  bdQr_p,
                     BCHandler& BCh_u ) :
    M_dataType                         ( dataNavierStokes ),
    _refFE_u                           ( refFE_u ),
    _refFE_p                           ( refFE_p ),
    _dof_u                             ( this->mesh(), _refFE_u ),
    _dof_p                             ( this->mesh(), _refFE_p ),
    _dim_u                             ( _dof_u.numTotalDof() ),
    _dim_p                             ( _dof_p.numTotalDof() ),
    _Qr_u                              ( Qr_u ),
    _bdQr_u                            ( bdQr_u ),
    _Qr_p                              ( Qr_p ),
    _bdQr_p                            ( bdQr_p ),
    _fe_u                              ( _refFE_u,
                                         getGeoMap( this->mesh() ),
                                         _Qr_u ),
    _feBd_u                            ( _refFE_u.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_u ),
    _fe_p                              ( _refFE_p,
                                         getGeoMap( this->mesh() ),
                                         _Qr_p ),
    _feBd_p                            ( _refFE_p.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_p ),
    _u                                 ( _dim_u ),
    _p                                 ( _dim_p ),
    _bdf                               ( this->order_bdf() ),
    _ns_post_proc                      ( this->mesh(), _feBd_u, _dof_u, NDIM ),
    _count                             ( 0 ),
    //! stuff to compute the fluxes at each section (ex. of mesh: tube20.mesh)
    M_nb_sections                      ( this->NbZSections() ),
    M_z_section                        ( M_nb_sections ),
    M_list_of_faces_on_section_velocity( M_nb_sections ),
    M_list_of_faces_on_section_pressure( M_nb_sections ),
    M_list_of_points_on_boundary       ( M_nb_sections ),
    M_out_areas                        ("Areas.res"),
#if COMPUTE_POLYGONAL_AREA
    M_out_areas_polygon                ("AreasPolygon.res"),
#endif
    M_out_fluxes                       ("Fluxes.res"),
    M_out_pressure                     ("Pressure.res"),
    M_BCh_fluid                        ( &BCh_u )
{
    std::cout << "New fluid constructor ... " << std::flush << std::endl;
    if ( this->computeMeanValuesPerSection() == 1 )
        initializeMeanValuesPerSection();
    /*
    //! for simple bifurcation mesh
        initializeSectionsBifurc();
    */
}



template <typename Mesh, typename DataType>
NavierStokesHandler<Mesh, DataType>::
NavierStokesHandler( const GetPot&   data_file,
                     const RefFE&    refFE_u,
                     const RefFE&    refFE_p,
                     const QuadRule& Qr_u,
                     const QuadRule& bdQr_u,
                     const QuadRule& Qr_p,
                     const QuadRule& bdQr_p) :
    M_dataType                         ( data_file ),
    _refFE_u                           ( refFE_u ),
    _refFE_p                           ( refFE_p ),
    _dof_u                             ( this->mesh(), _refFE_u ),
    _dof_p                             ( this->mesh(), _refFE_p ),
    _dim_u                             ( _dof_u.numTotalDof() ),
    _dim_p                             ( _dof_p.numTotalDof() ),
    _Qr_u                              ( Qr_u ),
    _bdQr_u                            ( bdQr_u ),
    _Qr_p                              ( Qr_p ),
    _bdQr_p                            ( bdQr_p ),
    _fe_u                              ( _refFE_u,
                                         getGeoMap( this->mesh() ),
                                         _Qr_u ),
    _feBd_u                            ( _refFE_u.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_u ),
    _fe_p                              ( _refFE_p,
                                         getGeoMap( this->mesh() ),
                                         _Qr_p ),
    _feBd_p                            ( _refFE_p.boundaryFE(),
                                         getGeoMap( this->mesh() ).boundaryMap(),
                                         _bdQr_p ),
    _u                                 ( _dim_u ),
    _p                                 ( _dim_p ),
    _bdf                               ( this->order_bdf() ),
    _ns_post_proc                      ( this->mesh(), _feBd_u, _dof_u, NDIM ),
    _count                             ( 0 ),
    //! stuff to compute the fluxes at each section (ex. of mesh: tube20.mesh)
    M_nb_sections                      ( this->NbZSections() ),
    M_z_section                        ( M_nb_sections ),
    M_list_of_faces_on_section_velocity( M_nb_sections ),
    M_list_of_faces_on_section_pressure( M_nb_sections ),
    M_list_of_points_on_boundary       ( M_nb_sections ),
    M_out_areas                        ("Areas.res"),
#if COMPUTE_POLYGONAL_AREA
    M_out_areas_polygon                ("AreasPolygon.res"),
#endif
    M_out_fluxes                       ("Fluxes.res"),
    M_out_pressure                     ("Pressure.res"),
    // fluid Bc conditions
    M_BCh_fluid                        ( 0 )
{
    std::cout << "New fluid constructor (2) ... " << std::endl;
    if ( this->computeMeanValuesPerSection() == 1 )
        initializeMeanValuesPerSection();
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

template <typename Mesh, typename DataType>
Dof&
NavierStokesHandler<Mesh, DataType>::uDof()
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


template <typename Mesh, typename DataType>
Dof&
NavierStokesHandler<Mesh, DataType>::pDof()
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

    if ( fmod( float( _count ), float( this->verbose() ) ) == 0.0 )
    {
        std::cout << "  F-  Post-processing \n";
        index << std::setfill('0') << std::setw(3);
        index << ( _count / this->verbose() );
        name = index.str();

        // postprocess data file for GMV
        wr_gmv_ascii( "test." + name + ".inp", this->mesh(), _dim_u,
                      _u.giveVec(), _p.giveVec() );

        // postprocess data file for medit
        wr_medit_ascii_scalar( "press." + name + ".bb", _p.giveVec(),
                               _p.size() );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", _u.giveVec(),
                               this->mesh().numVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", _u.giveVec() + _dim_u,
                               this->mesh().numVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", _u.giveVec()+2*_dim_u,
                               this->mesh().numVertices() );
        system( ( "ln -s -f " + this->meshDir() + this->meshFile() +
                  " press." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + this->meshDir() + this->meshFile() +
                  " vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + this->meshDir() + this->meshFile() +
                  " vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s -f " + this->meshDir() + this->meshFile() +
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

    wr_opendx_header( file_pre, this->mesh(), _dof_p, _fe_p, fe_type_pre );
    wr_opendx_scalar( file_pre, "pression", _p );

    wr_opendx_header( file_vel, this->mesh(), _dof_u, _fe_u, fe_type_vel );
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


    _bdf.bdf_u().initialize_unk( u0, this->mesh(), _refFE_u, _fe_u, _dof_u, t0,
                                 dt, nbComp );

    // initialize _u with the first element in bdf_u.unk (=last value)
    _u = *( _bdf.bdf_u().unk().begin() );

    _bdf.bdf_p().initialize_unk( p0, this->mesh(), _refFE_p, _fe_p, _dof_p, t0,
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

    resfile.read( ( char* ) & _u( 0 ), _u.size() * sizeof( double ) );
    resfile.read( ( char* ) & _p( 0 ), _p.size() * sizeof( double ) );
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
//    int dotpos;

    int nDof = uDof().numTotalDof();

    std::string filenamex = velName;

    ext = "_x.bb";
//    dotpos = filenamex.find(".");

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
        std::cout << "Reading file " << filenamez
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


    _bdf.bdf_u().initialize_unk( _u );
    _bdf.bdf_p().initialize_unk( _p );


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
    FacesOnFlag( flag, _refFE_u, _dof_u, faces_on_flag ); //! reconstruct all the connectivity

    //! reference normal vector (unused, as we keep the positive outward normal on the boundary)
    SimpleVect<Real,0> unused(3);
    unused(0) = 0.;
    unused(1) = 0.;
    unused(2) = 0.;

    std::pair< Real, Real > area_flux = AreaAndFlux( faces_on_flag,
                                                     false, unused ); //! keep the orientation of the normal
    return area_flux.second;
}

/*! Get the faces that are on a part of the boundary which is defined by its flag

   \param flag : __flag of the faces that are sought
   \param __reffe : reference fe.
   \param __dof : dof.
   \param __faces_on_section : list of faces that were found:
   first:  global face number
   second: vector of size nDofF (nb of dof per face).
   contains the local (in the face) to global dof.
*/
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::FacesOnFlag( const EntityFlag& __flag ,
                                                  const RefFE& __reffe,  const Dof& __dof,
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

    UInt bdnF = this->mesh().numBFaces();    // number of faces on boundary

    EntityFlag marker;
    typedef face_dof_type::iterator FaceDofIterator;

    //
    // Loop on boundary faces: List of boundary faces
    // with marker = __flag
    //
    for ( ID i = 1 ; i <= bdnF; ++i )
    {
        marker = this->mesh().boundaryFace( i ).marker();
        if ( marker == __flag )
        {
            __faces_on_flag.push_front( make_pair( i, SimpleVect<ID>( nDofF ) ) );
        }
    }

    //! building the local to global vector for these boundary faces
    _updateDofFaces( __faces_on_flag, __reffe, __dof );

}

/*! Get the faces that live on a horizontal section (z=__z_section)
    up to a given tolerance (z = __z_section +/- __tolerance ).

    \param __planar_section : position of the planar section (plane
    given by its 4 (THREED) components (a,b,c,d) such that
    a x + b y + c z + d = 0

    \param __faces_on_section : list of faces that were found:
      first:  global face number
      second: vector of size nDofF (nb of dof per face).
      contains the local (in the face) to global dof.

      The first is for velocity ("_u") and the second
      for pressure ("_p").
   \param __tolerance_section :
                  all faces between the two planes
		  a x + b y + c z + d = +/- __tolerance_section
		  are kept.
   (we only check the first 3 points of the mesh faces)

   \param __x_frontier : given point on the boundary.
   \param __point_on_boundary is a pair<global face ID, Vertex ID> determining a point
          close to a given boundary point __x_frontier.

   This function is VERY mesh dependent!! Use it with caution.
*/
template <typename Mesh, typename DataType>
void
NavierStokesHandler<Mesh, DataType>::FacesOnSection( const SimpleVect<Real,0> &  __planar_section,
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
    UInt nF = this->mesh().numFaces();    // number of faces (all faces)

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
    //! default value
    __point_on_boundary = std::pair<ID, ID>( 0 , 0 );

    //! coefficients defining the plane ax + by + cz + d = 0
    Real a_plane = __planar_section(0);
    Real b_plane = __planar_section(1);
    Real c_plane = __planar_section(2);
    Real d_plane = __planar_section(3);

    /*
    int verbose = 0;
    if ( verbose > 1 )
        std::cout << "a " << a_plane << " b " << b_plane << " c "
                  << c_plane << " d " << d_plane << std::endl;
    */

    //! search the faces on a given horizontal section
    //! loop on all faces
    for( ID iface = 1 ; iface <= nF ; iface ++ ) {
        iglobface = this->mesh().faceList(iface).id();

        //! 3 points on the faces
        Real x_VeFa1 = this->mesh().faceList(iface).point( 1 ).x();
        Real y_VeFa1 = this->mesh().faceList(iface).point( 1 ).y();
        Real z_VeFa1 = this->mesh().faceList(iface).point( 1 ).z();

        Real x_VeFa2 = this->mesh().faceList(iface).point( 2 ).x();
        Real y_VeFa2 = this->mesh().faceList(iface).point( 2 ).y();
        Real z_VeFa2 = this->mesh().faceList(iface).point( 2 ).z();

        Real x_VeFa3 = this->mesh().faceList(iface).point( 3 ).x();
        Real y_VeFa3 = this->mesh().faceList(iface).point( 3 ).y();
        Real z_VeFa3 = this->mesh().faceList(iface).point( 3 ).z();

        //! check the tolerance here (depends on the mesh)
        if ( std::fabs( a_plane * x_VeFa1 + b_plane * y_VeFa1 + c_plane * z_VeFa1 + d_plane ) < __tolerance_section &&
             std::fabs( a_plane * x_VeFa2 + b_plane * y_VeFa2 + c_plane * z_VeFa2 + d_plane ) < __tolerance_section &&
             std::fabs( a_plane * x_VeFa3 + b_plane * y_VeFa3 + c_plane * z_VeFa3 + d_plane ) < __tolerance_section
             ) {
            //! GEOMETRY : additional inequality constraint:
            /*
            //! 3<y<3.5 for sections x=+/- 0.5 and
            //! -0.5<x<0.5 for section y=4
            if ( ! ( a_plane == 1. && std::fabs( d_plane ) == 0.5 ) ||
                 ! ( b_plane == 1. && d_plane == - 4. ) ||
                 ( a_plane == 1. && std::fabs( d_plane ) == 0.5 &&
                   y_VeFa1 > 3. - __tolerance_section && y_VeFa1 < 3.5 + __tolerance_section &&
                   y_VeFa2 > 3. - __tolerance_section && y_VeFa2 < 3.5 + __tolerance_section &&
                   y_VeFa3 > 3. - __tolerance_section && y_VeFa3 < 3.5 + __tolerance_section    ) ||
                 ( b_plane == 1. && d_plane == -4. &&
                   std::fabs( x_VeFa1 ) < 0.5 + __tolerance_section &&
                   std::fabs( x_VeFa2 ) < 0.5 + __tolerance_section &&
                   std::fabs( x_VeFa3 ) < 0.5 + __tolerance_section )
                 ) {
            */
            //! plane z=2.5 : only y<0
            if ( ! ( std::fabs( c_plane - 1. ) < __tolerance_section
                     && std::fabs( d_plane + 2.5 ) < __tolerance_section ) ||
                 ( std::fabs( c_plane - 1. ) < __tolerance_section &&
                   std::fabs( d_plane + 2.5 ) < __tolerance_section &&
                   y_VeFa1 < 0. + __tolerance_section &&
                   y_VeFa2 < 0. + __tolerance_section &&
                   y_VeFa3 < 0. + __tolerance_section )
                 ) {
                /*
                if ( verbose > 1 )
                    std::cout << "Found:" << iglobface << " x " << x_VeFa1 << " y " << y_VeFa1
                              << " z " << z_VeFa1
                              << std::endl;
                */

                __faces_on_section_u.push_front( make_pair( iglobface, SimpleVect<ID>( nDofF_u ) ) );
                __faces_on_section_p.push_front( make_pair( iglobface, SimpleVect<ID>( nDofF_p ) ) );

                ID iVeFa = 1;
                while ( ! _found_point && iVeFa <= nFaceV ) {
                    Real x_VeFa = this->mesh().faceList(iface).point( iVeFa ).x();
                    if ( std::fabs( x_VeFa - __x_frontier ) < __tolerance_section ) {
                        __point_on_boundary = std::pair<ID, ID>( iglobface , iVeFa );
                        _found_point = true;
                    }
                    iVeFa ++;
                }
            }
        }
    }

    //! building the local to global vector for these boundary faces (veloc and pressure)
    _updateDofFaces( __faces_on_section_u, _refFE_u, _dof_u );
    _updateDofFaces( __faces_on_section_p, _refFE_p, _dof_p );
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
                                                      const RefFE& __reffe, const Dof& __dof)
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
        iElAd = this->mesh().boundaryFace( iFace ).ad_first();

        // local id of the face in its adjacent element
        iFaEl = this->mesh().boundaryFace( iFace ).pos_first();

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
                    //!!!!!!!!!! BUG : GENERAL DOF here CORRECT this!!!!!!
                    gDof = __dof.localToGlobal( iElAd,
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
                    gDof = __dof.localToGlobal( iElAd,
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
                gDof = __dof.localToGlobal( iElAd,
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

    Real radius = this->mesh().faceList(iglobface).point( iVeFa ).x();

    return  0.5 * __nb_edge_polygon * radius * radius * std::sin( 2. * LifeV::Pi / __nb_edge_polygon );
}



/*! compute the Area and Flux
  \param __faces_on_section : list of (faces <-> veloc dof) that were found:
  first:  global face number
  second: vector of size nDofF (nb of dof per face)
  contains the local (in the face) to global dof.
  \param __modify_sign_normal : if == true, the computed flux is positive if (u \dot __reference_normal) is positive
  \param __reference_normal : vector provided to give a fixed orientation to the normal
*/
template <typename Mesh, typename DataType>
std::pair< Real, Real >
NavierStokesHandler<Mesh, DataType>::AreaAndFlux( const face_dof_type & __faces_on_section_u,
                                                  const bool & __modify_sign_normal,
                                                  const SimpleVect<Real,0>& __reference_normal )
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

    // Number of velocity components
    UInt nc_u = _u.nbcomp();

    // Nodal values of the velocity in the current face
    std::vector<Real> u_local( nc_u * nDofF );

    //    int verbose = 2;

    // Loop on faces
    for ( const_face_dof_iterator_type j = __faces_on_section_u.begin();
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
        _feBd_u.updateMeasNormalQuadPt( this->mesh().boundaryFace( j->first ) );

        __area += _feBd_u.measure();

        /*
        if ( verbose > 1)
            std::cout << "area/flux comp : " << j->first << " area= " << __area << std::endl;
        */

        //! check the orientation of the normal
        //! the flux will be positive if (u \dot positive_normal) is positive
        Real sign_normal = 1.0;
        Int i_quadr = 0; //! check only at the 1rst quad point
        //(assuming that the element is not too distorted)
        Real n_dot_ref = 0.;
        for ( int ic = 0; ic < (int)nc_u; ++ic ) { //! dot product
            n_dot_ref =_feBd_u.normal( ic, i_quadr ) * __reference_normal( ic );
        }
        if ( __modify_sign_normal && n_dot_ref < 0 ) {
            sign_normal = -1.0;
        }

        // Quadrature formula
        // Loop on quadrature points for velocities
        for ( int iq = 0; iq < _feBd_u.nbQuadPt; ++iq )
        {
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
                        sign_normal;
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

    // Nodal values of the pressure in the current face
    std::vector<Real> p_local( nDofF );

    Real area = 0.;

    // Loop on faces
    for ( const_face_dof_iterator_type j = __faces_on_section_p.begin();
          j != __faces_on_section_p.end(); ++j )
    {

        for ( ID l = 1; l <= nDofF; ++l ){
            gDof = j->second( l );
            p_local[ l - 1 ] = _p( gDof - 1 );
        }

        // Updating quadrature data on the current face
        _feBd_p.updateMeasQuadPt( this->mesh().boundaryFace( j->first ) );
        area += _feBd_p.measure();

        // Quadrature formula
        // Loop on quadrature points for velocities
        for ( int iq = 0; iq < _feBd_p.nbQuadPt; ++iq )
        {
            // Interpolation
            // Loop on local dof
            for ( ID l = 1; l <= nDofF; ++l )
                __pressure += _feBd_p.weightMeas( iq ) *
                    p_local[ l - 1 ] *
                    _feBd_p.phi( int( l - 1 ), iq );

        }
    }
    // divide by the (approximated in the case of FSI) area to obtain the mean.
    __pressure = __pressure / area;

    return  __pressure;
}

/*! compute the mean Pressures, Areas and Fluxes for all sections and write them on a file
    in the plotmtv format
*/
template <typename Mesh, typename DataType>
void NavierStokesHandler<Mesh, DataType>::PostProcessPressureAreaAndFlux( const Real & __time )
{
    if ( this->computeMeanValuesPerSection() != 1 )
        ERROR_MSG("This function is disabled, if you don't ask to compute the Mean Values. (in data file)");

    M_out_areas  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                 << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Area'\n";
    M_out_fluxes << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                 << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Flux'\n";
    M_out_pressure  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                    << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='Pressure'\n";
#if COMPUTE_POLYGONAL_AREA
    M_out_areas_polygon  << "$ DATA = CURVE2D\n %% xlabel='z'\n"
                         << "%% toplabel='Section,time=" << __time << "'\n %% ylabel='AreaPolygonal'\n";
#endif

    SimpleVect<Real,0> reference_normal(3);
    reference_normal(0) = 1.;
    reference_normal(1) = 0.;
    reference_normal(2) = 0.; // flux positive if (u dot reference_normal) > 0

    for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
        //! all normals oriented according to the reference_normal -> "true"
        std::pair<Real, Real> AQmid  = this->AreaAndFlux( M_list_of_faces_on_section_velocity[izs], true, reference_normal );
        M_out_areas  << M_z_section[izs] << "\t" << AQmid.first << "\n";
        M_out_fluxes << M_z_section[izs] << "\t" << AQmid.second << "\n";
        M_out_pressure << M_z_section[izs] << "\t" << this->MeanPressure( M_list_of_faces_on_section_pressure[izs] ) << "\n";

#if COMPUTE_POLYGONAL_AREA
        //! remove the polygonal area??
        Real area_polygon = this->AreaCylindric( M_list_of_points_on_boundary[izs], NbPolygonEdges() );
        M_out_areas_polygon  << M_z_section[izs] << "\t" << area_polygon << "\n";
#endif
    }
    M_out_areas << std::endl;
    M_out_fluxes << std::endl;
    M_out_pressure << std::endl;
#if COMPUTE_POLYGONAL_AREA
    M_out_areas_polygon << std::endl;
#endif
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
    for ( ID iElem = 1; iElem <= this->mesh().numVolumes(); ++iElem )
    {

        _fe_u.updateJac( this->mesh().volume( iElem ) );

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
    for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
    {
        _fe_p.updateFirstDeriv( this->mesh().volumeList( iVol ) );
        sum2 += elem_L2_diff_2( this->_p, pexact, this->_fe_p, this->_dof_p, time, 1 );
        sum1 += elem_integral_diff( this->_p, pexact, this->_fe_p, this->_dof_p, time, 1 );
        sum0 += _fe_p.measure();
        if (relError)
        {
            sumExact2 += elem_L2_2( pexact, this->_fe_p, time, 1 );
            sumExact1 += elem_integral( pexact, this->_fe_p, time, 1 );
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
    for ( UInt iVol = 1; iVol <= this->mesh().numVolumes(); iVol++ )
    {
        _fe_u.updateFirstDeriv( this->mesh().volumeList( iVol ) );
        normU += elem_L2_diff_2( _u, uexact, this->_fe_u, this->_dof_u, time,
                                 int( nbCompU ) );
        if (relError)
        {
            sumExact += elem_L2_2( uexact, this->_fe_u, time, int( nbCompU ) );
        }
    }
    if (relError)
    {
        *relError = sqrt( normU / sumExact );
    }
    return sqrt( normU );
}



  template <typename Mesh, typename DataType>
  BasePattern::PatternType NavierStokesHandler<Mesh, DataType>::patternType() {

    BasePattern::PatternType  pt = BasePattern::STANDARD_PATTERN;
    if ( M_dataType.stabilization() == IP_STABILIZATION )
      pt = BasePattern::EDGE_COUPLING_PATTERN;
    return pt;
  }





template <typename Mesh, typename DataType>
void NavierStokesHandler<Mesh, DataType>::initializeMeanValuesPerSection()
{
    //---------------
    //! stuff to compute the fluxes at each section for a cylindrical tube
    //! (ex. of mesh: tube20.mesh)
    //---------------

    if ( ! this->mesh().hasInternalFaces() )
        ERROR_MSG("The mesh must have all internal faces built up.");// Check that 'mesh_faces = all' in the data file.");
    if ( M_nb_sections < 2 )
        ERROR_MSG("We can't compute the mean values on less than 2 sections.");
    ASSERT( this->ZSectionFinal() - this->ZSectionInit() > 0,
            "Error on the z given to compute the sections.");

    for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
        M_z_section[ izs ] = this->ZSectionInit() + Real(izs) * ( this->ZSectionFinal() - this->ZSectionInit() )
            / Real(M_nb_sections-1); // mesh dependent (length of tube=5)
    }

    //! coefficients defining the planar section: ax + by + cz + d = 0
    SimpleVect<Real,0> PlaneCoeff(4);
    //! Assuming a plane parallel to the Oz axis
    PlaneCoeff(0) = 0.; // a
    PlaneCoeff(1) = 0.; // b
    PlaneCoeff(2) = 1.; // c
    PlaneCoeff(3) = 0.; // d

    for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
        PlaneCoeff(3) = - M_z_section[izs];
        //! for velocity
        this->FacesOnSection( PlaneCoeff,
                              M_list_of_faces_on_section_velocity[izs],
                              M_list_of_faces_on_section_pressure[izs],
                              this->ToleranceSection(),
                              this->XSectionFrontier(), M_list_of_points_on_boundary[izs] );
        if ( M_list_of_faces_on_section_velocity[izs].size() == 0 ||
             M_list_of_faces_on_section_pressure[izs].size() == 0  ) {
            std::cout << "section z=" << M_z_section[izs] << " size="
                      << M_list_of_faces_on_section_velocity[izs].size() << std::endl;
            ERROR_MSG("For this section, no faces found.");
        }
#if COMPUTE_POLYGONAL_AREA
        if ( M_list_of_points_on_boundary[izs].first  == 0 &&
             M_list_of_points_on_boundary[izs].second == 0 ) {
            std::cout << "section z=" << M_z_section[izs]
                      << " (point=" <<  XSectionFrontier() << ")" << std::endl;
            ERROR_MSG("For this section, no boundary points found.");
        }
#endif
    }
    if ( M_list_of_faces_on_section_velocity.size() == 0 ||
         M_list_of_faces_on_section_pressure.size() == 0 ) {
        ERROR_MSG("No list of faces found.");
    }
}

template <typename Mesh, typename DataType>
void NavierStokesHandler<Mesh, DataType>::initializeSectionsBifurc()
{
    //---------------
    //! To obtain geometrically the internal faces representing the stent
    //! (ex. of mesh: bif3d_*.mesh)
    //---------------

    if ( ! this->mesh().hasInternalFaces() )
        ERROR_MSG("The mesh must have all internal faces built up.");// Check that 'mesh_faces = all' in the data file."); Check that 'mesh_faces = all' in the data file.");
    if ( M_nb_sections != 3 )
        ERROR_MSG("We a priori know that the stent lives on 3 planar sections.");

    //! coefficients defining the planar section: ax + by + cz + d = 0
    SimpleVect<Real,0> PlaneCoeff(4);

    //! values of "d"
    M_z_section[ 0 ] = 0.5;  // plane x=0.5
    M_z_section[ 1 ] = -0.5; // plane x=-0.5
    M_z_section[ 2 ] = 4;    // plane y=4

    for ( int izs = 0; izs < M_nb_sections ; izs ++  ){
        switch ( izs ) {
        case 0:
            //! Right part of the stent :  x=0.5
        case 1:
            //! Left part of the stent :  x=-0.5
            PlaneCoeff(0) = 1.; // a
            PlaneCoeff(1) = 0.; // b
            PlaneCoeff(2) = 0.; // c
            PlaneCoeff(3) = - M_z_section[ izs ]; // d
            break;
        case 2:
            //! Aneurismal part of the stent :  y=4
            PlaneCoeff(0) = 0.; // a
            PlaneCoeff(1) = 1.; // b
            PlaneCoeff(2) = 0.; // c
            PlaneCoeff(3) = - M_z_section[ izs ]; // d
            break;
        default:
            ERROR_MSG("No such planar section");
        }

        this->FacesOnSection( PlaneCoeff,
                              M_list_of_faces_on_section_velocity[izs],
                              M_list_of_faces_on_section_pressure[izs],
                              this->ToleranceSection(),
                              this->XSectionFrontier(), M_list_of_points_on_boundary[izs] );
        if ( M_list_of_faces_on_section_velocity[izs].size() == 0 ||
             M_list_of_faces_on_section_pressure[izs].size() == 0  ) {
            std::cout << "section z=" << M_z_section[izs] << " size="
                      << M_list_of_faces_on_section_velocity[izs].size() << std::endl;
            ERROR_MSG("For this section, no faces found.");
        }
#if COMPUTE_POLYGONAL_AREA
        if ( M_list_of_points_on_boundary[izs].first  == 0 &&
             M_list_of_points_on_boundary[izs].second == 0 ) {
            std::cout << "section z=" << M_z_section[izs]
                      << " (point=" <<  XSectionFrontier() << ")" << std::endl;
            ERROR_MSG("For this section, no boundary points found.");
        }
#endif
    }

    if ( M_list_of_faces_on_section_velocity.size() == 0 ||
         M_list_of_faces_on_section_pressure.size() == 0 ) {
        ERROR_MSG("No list of faces found.");
    }
}


}
#endif
