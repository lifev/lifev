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
  \file postProc.hpp
  \author A. Veneziani
  \author T. Passerini
  \date 09/2008

  \brief File containing a class with all the methods for a post processing of the solution.

 */
#ifndef _POST_PROC_H
#define _POST_PROC_H

#include <string>
#include <iostream>
#include <sstream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>

namespace LifeV
{
  /*!
  \class PostProc

  \brief class with all the methods for post processing the solution.

  The data structures are based on vector _bdLtoG:
    - the size is equal to the number of boundary elements in the mesh
    - each component is a vector of identifiers, corresponding to the "local"
      numbering of the face dof's
    - in the previous sentence, "local" is intended as "local to this class":
      an additional structure gives back the "traditional" local to global mapping
    - careful: the "local" numbering of faces starts from 1

  Vector _fBdToIn has this usage:
    - the size is equal to the total number of boundary dof's
    - the i-th component is the global numbering of the i-th boundary dof

  How to use _bdLtoG and _fBdToIn together?
    _bdLtoG[ j ] is a vector containing the dof "local" identifiers for boundary face j
    _fBdToIn[ _bdLtoG[ j ][ i ] ] is the global ID associated to the i-th dof on j-th face

  NB: here "global" is intended wrt the global mesh (prior to partitioning, if this applies)
   */
  template <typename Mesh>
  class PostProc
  {
    typedef typename Mesh::VolumeShape GeoShape;
    typedef typename GeoShape::GeoBShape GeoBShape;
    typedef boost::shared_ptr<Mesh> mesh_ptrtype;

  public:
    /*! @name Constructors */
    //@{
    /*!
     \brief Constructor

     In this general case we allow to pass an arbitrary number (nvar) of dof/fe classes

     NOTE: in most parts of the code, only feBd[0] and dof[0] will be considered

     \param mesh the mesh
     \param feBd the finite element prototype for boundary elements
     \param dof class for the description of degrees of freedom
     \param _map the map describing the partition over different processors
     \param nvar number of elements in vectors feBd, dof
     */
    PostProc<Mesh>( mesh_ptrtype mesh,
        std::vector<CurrentBdFE* > feBd,
        std::vector<Dof* > dof,
        const EpetraMap& _map, UInt nvar = 1 );

    /*!
     \brief Constructor for the case in which we have only one feBd and dof
     */
    PostProc<Mesh>( mesh_ptrtype mesh,
        CurrentBdFE* feBd, Dof* dof,
        const EpetraMap& _map );

    /*!
     \brief Constructor for tha case in which we have two feBd's and dof's

     This is the case of NS problem, in which we will want to consider both
     velocity and pressure dof/fe classes
     */
    PostProc<Mesh>( mesh_ptrtype mesh,
        CurrentBdFE* feBdu, Dof* dofu,
        CurrentBdFE* feBdp, Dof* dofp,
        const EpetraMap& _map );

    //@}

    /*! @name Methods */
    //@{
    /*!
     These methods build the data structures describing patches of boundary elements
     with respect to patch area, patch normal and patch test function
     */
    void set_area();
    void set_normal();
    void set_phi();

    /*!
     These methods compute quantities on boundary sections associated to a given flag
     */
    Real area( const EntityFlag& flag );
    /*!
      \tparam V Vector type
      \param sol is intended to be a vector, this method computes the flux of sol
        across section "flag"
     */
    template< typename V >
    Real flux( const V& sol, const EntityFlag& flag, UInt feSpace = 0, UInt nDim = nDimensions);
    /*!
      \tparam V Vector type. Basic policy for type V: operator[] available

      \param sol is intended to be a vector or a scalar (pressure in NS problem),
       this method computes the average value of sol on section "flag"

       \return the averaged vector
     */
    template< typename V >
    Vector average( const V& sol, const EntityFlag& flag, UInt feSpace = 0, UInt nDim = 1);

// NOT READY!
#if 0
    /*!
      This procedure computes the tangential stresses on the boundary of the domain
      starting from an estimate of the viscous stresses on the boundary.

      The basic approach resorts to weak (residual-based) computing of the stresses.

      \param r contains the estimate of the viscous stresses
      \param residual switches two possibilities:
        If (residual == true), vector r contains the residual of the fluid problem,
        i. e. the integral over the domain of the product of the stress and the test functions
          r_i = \int_Omega tau \cdot \phi_i domega

        If (residual == false), vector r represents an estimate of the stress
        in each dof of the _boundary_
          r_i = tau_i
     */
    template<typename V>
    V compute_sstress( const V& r, UInt ncomp, bool residual ) const;

    /*!
      Here we compute the stress over the boundary, _given_ the velocity gradient

      s(n) = nu ( grad(u) + grad(u)^T ) \cdot n
      component-wise:
      s_i = nu ( d u_i / d x_j + d u_j / d x_i ) n_j

      \param grad the velocity gradient
        grad[ idGlobalDof-1 + ncomp*dim + scomp*nDimensions*dim ] =
         = ( d u_{scomp} / d x_ncomp +  d u_{ucomp} / d x_scomp )
         (at node idGlobalDof)
     */
    template<typename V>
    V compute_stress( const V& grad, const UInt& dim, const Real& visc  ) const;
#endif
    //@}

    /*! @name Methods */
    //@{
    /*!
     These methods print to screen information about the class data structures
     */
    void show_bdLtoG();

    void show_area();

    void show_normal();

    void show_phi();
    //@}

    /*! @name Getters */
    //@{
    /*!
     Access to private members
     */
    std::vector<Real> patch_area() const { return _area; }

    std::vector<Real> patch_normal() const { return _normal; }

    std::vector< ID > fBdToIn() const { return _fBdToIn; }

    //! Number of boundary DOF for the mesh at hand
    UInt nBdDof() const { return _nBdDof; }
    //@}


  private:
    void build_vectors();

    UInt _nvar;

    std::vector<UInt> _nBdDof;
    std::vector<UInt> _M_nDofpV, _M_nDofpE, _M_nDofpF;
    UInt _M_nFaceV, _M_nFaceE;
    UInt _M_nElemV, _M_nElemE;
    std::vector<UInt> _M_nDofFV, _M_nDofFE;
    std::vector<UInt> _M_nDofF;
    std::vector<UInt> _M_nDofElemV, _M_nDofElemE;
    UInt _M_bdnF;
    std::vector<UInt> _M_nTotalDof;
    // vector whose ith-component is the area of the patch of the ith-node on the boundary
    std::vector<std::vector<Real> > _area;
    // vector with the components of the average normal vector of each boundary node
    std::vector<std::vector<Real> > _normal;
    // vector with \int_{Patch(i)} \phi_i dx where \phi_i is the basis function.
    std::vector<std::vector<Real> > _phi;

    // store once for all a map, with key={boundary flag}, value={ID list}
    std::map< EntityFlag, std::list<ID> > _faces_with_flag;

    // for each boundary face, it contains the numbering of the dof of the face
    std::vector< std::vector< SimpleVect<ID> > > _bdLtoG;
    // it converts from a local numeration over the boundary faces on the global numeration of the mesh
    std::vector< std::vector< ID > > _fBdToIn;

    // reference to the current boundary FE
    std::vector<CurrentBdFE* > _feBd;
    // reference to a Dof class
    std::vector<Dof* > _dof;
    // pointer to the mesh
    mesh_ptrtype _mesh;
    // pointer to the processor mapping
    boost::shared_ptr<EpetraMap>  _epetraMap;

  };

  //
  // IMPLEMENTATIONS
  //
  template <typename Mesh>
  PostProc<Mesh>::PostProc( mesh_ptrtype mesh,
      std::vector<CurrentBdFE* > feBd,
      std::vector<Dof* > dof,
      const EpetraMap& _map, UInt nvar )
  : _nvar(nvar),
  _nBdDof(_nvar), _M_nDofpV(_nvar), _M_nDofpE(_nvar), _M_nDofpF(_nvar),
  _M_nDofFV(_nvar), _M_nDofFE(_nvar), _M_nDofF(_nvar),
  _M_nDofElemV(_nvar), _M_nDofElemE(_nvar), _M_nTotalDof(_nvar),
  _area(_nvar), _normal(_nvar), _phi(_nvar),
  _bdLtoG(_nvar), _fBdToIn(_nvar), _feBd(feBd), _dof(dof),
  _mesh( mesh ), _epetraMap( new EpetraMap(_map) )
  {
    for(UInt var=0; var<_nvar; ++var){
      _bdLtoG[_nvar].clear(); _fBdToIn[var].clear();
    }
    // Some useful local variables, to save some typing
    _M_nFaceV = GeoBShape::numVertices; // Number of face's vertices
    _M_nFaceE = GeoBShape::numEdges;    // Number of face's edges

    _M_nElemV = GeoShape::numVertices; // Number of element's vertices
    _M_nElemE = GeoShape::numEdges;    // Number of element's edges

    _M_bdnF = _mesh->numBFaces();    // number of faces on boundary

    // Construction of _bdLtoG & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
  }

  template <typename Mesh>
  PostProc<Mesh>::PostProc( mesh_ptrtype mesh,
      CurrentBdFE* feBd, Dof* dof,
      const EpetraMap& _map )
  : _nvar(1),
  _nBdDof(_nvar), _M_nDofpV(_nvar), _M_nDofpE(_nvar), _M_nDofpF(_nvar),
  _M_nDofFV(_nvar), _M_nDofFE(_nvar), _M_nDofF(_nvar),
  _M_nDofElemV(_nvar), _M_nDofElemE(_nvar), _M_nTotalDof(_nvar),
  _area(_nvar), _normal(_nvar), _phi(_nvar),
  _bdLtoG(_nvar), _fBdToIn(_nvar), _feBd(_nvar), _dof(_nvar),
  _mesh( mesh ), _epetraMap( new EpetraMap(_map) )
  {
    _feBd[0]=feBd; _dof[0]=dof;

    // Some useful local variables, to save some typing
    _M_nFaceV = GeoBShape::numVertices; // Number of face's vertices
    _M_nFaceE = GeoBShape::numEdges;    // Number of face's edges

    _M_nElemV = GeoShape::numVertices; // Number of element's vertices
    _M_nElemE = GeoShape::numEdges;    // Number of element's edges

    _M_bdnF = _mesh->numBFaces();    // number of faces on boundary

    // Construction of _bdLtoG & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
  }

  template <typename Mesh>
  PostProc<Mesh>::PostProc( mesh_ptrtype mesh,
      CurrentBdFE* feBdu, Dof* dofu,
      CurrentBdFE* feBdp, Dof* dofp,
      const EpetraMap& _map )
  : _nvar(2),
  _nBdDof(_nvar), _M_nDofpV(_nvar), _M_nDofpE(_nvar), _M_nDofpF(_nvar),
  _M_nDofFV(_nvar), _M_nDofFE(_nvar), _M_nDofF(_nvar),
  _M_nDofElemV(_nvar), _M_nDofElemE(_nvar), _M_nTotalDof(_nvar),
  _area(_nvar), _normal(_nvar), _phi(_nvar),
  _bdLtoG(_nvar), _fBdToIn(_nvar), _feBd(_nvar), _dof(_nvar),
  _mesh( mesh ), _epetraMap( new EpetraMap(_map) )
  {
    _feBd[0] = feBdu; _dof[0] = dofu;
    _feBd[1] = feBdp; _dof[1] = dofp;

    // Some useful local variables, to save some typing
    _M_nFaceV = GeoBShape::numVertices; // Number of face's vertices
    _M_nFaceE = GeoBShape::numEdges;    // Number of face's edges

    _M_nElemV = GeoShape::numVertices; // Number of element's vertices
    _M_nElemE = GeoShape::numEdges;    // Number of element's edges

    _M_bdnF = _mesh->numBFaces();    // number of faces on boundary

    // Construction of _bdLtoG & other data structures
    build_vectors();

    set_area();
    set_normal();
    set_phi();
  }

  // Area of faces with a certain marker
  template<typename Mesh>
  void PostProc<Mesh>::build_vectors()
  {
    SimpleVect<ID> bdltg;

    UInt iElAd, iVeEl, iFaEl, iEdEl;
    ID lDof, gDof, auxDof;
    std::vector<ID> bDof(_nvar);
    std::vector<ID>::iterator vidit;
    EntityFlag _flag;

    for(UInt var=0; var<_nvar; ++var) {

      bDof[var] = 1;

      // Some useful local variables, to save some typing
      _M_nDofpV[var] = _dof[var]->fe.nbDofPerVertex(); // number of Dof per vertices
      _M_nDofpE[var] = _dof[var]->fe.nbDofPerEdge();   // number of Dof per edges
      _M_nDofpF[var] = _dof[var]->fe.nbDofPerFace();   // number of Dof per faces

      _M_nDofFV[var] = _M_nDofpV[var] * _M_nFaceV; // number of vertex's Dof on a face
      _M_nDofFE[var] = _M_nDofpE[var] * _M_nFaceE; // number of edge's Dof on a face

      // number of total Dof on a face
      _M_nDofF[var] = _M_nDofFV[var] + _M_nDofFE[var] + _M_nDofpF[var];

      _M_nDofElemV[var] = _M_nElemV * _M_nDofpV[var]; // number of vertex's Dof on a Element
      _M_nDofElemE[var] = _M_nElemE * _M_nDofpE[var]; // number of edge's Dof on a Element

      _M_nTotalDof[var] = _dof[var]->numTotalDof();
    }

    // ===================================================
    // Loop on boundary faces
    // ===================================================
    for ( ID ibF = 1 ; ibF <= _M_bdnF; ++ibF )
      {

        iElAd = _mesh->bElement( ibF ).ad_first();  // id of the element adjacent to the face
        iFaEl = _mesh->bElement( ibF ).pos_first(); // local id of the face in its adjacent element

        _flag = _mesh->bElement(ibF).marker();
        _faces_with_flag[_flag].push_back( ibF ); // fill the flag-to-faceIdList map

        for(UInt var=0; var<_nvar; ++var) {

          bdltg.clean(); bdltg.resize( _M_nDofF[var] );

          _feBd[var]->updateMeas( _mesh->bElement( ibF ) ); // updating finite element information

          // ===================================================
          // Vertex based Dof
          // ===================================================
          if ( _M_nDofpV[var] )
            {

              // loop on face vertices
              for ( ID iVeFa = 1; iVeFa <= _M_nFaceV; ++iVeFa )
                {
                  iVeEl = GeoShape::fToP( iFaEl, iVeFa ); // local vertex number (in element)

                  // Loop number of Dof per vertex
                  for ( ID l = 1; l <= _M_nDofpV[var]; ++l )
                    {
                      lDof = ( iVeFa - 1 ) * _M_nDofpV[var] + l ; // local Dof
                      gDof = _dof[var]->localToGlobal( iElAd, ( iVeEl - 1 ) * _M_nDofpV[var] + l ); // global Dof
                      vidit = find( _fBdToIn[var].begin(), _fBdToIn[var].end(), gDof );
                      if ( vidit == _fBdToIn[var].end() )
                        { // the gDof has been encountered for the first time
                          bdltg( lDof ) = bDof[var];
                          _fBdToIn[var].push_back( gDof ); // local to boundary global on this face
                          bDof[var]++;
                        }
                      else
                        { // the gDof has been already inserted in the _fBdToIn vector
                          auxDof = ( ID ) ( ( vidit - _fBdToIn[var].begin() ) ) + 1;
                          bdltg( lDof ) = auxDof; // local to boundary global on this face
                        }
                    }
                }
            }
          // ===================================================
          // Edge based Dof
          // ===================================================
          if ( _M_nDofpE[var] )
            {
              // loop on face edges
              for ( ID iEdFa = 1; iEdFa <= _M_nFaceV; ++iEdFa )
                {
                  iEdEl = GeoShape::fToE( iFaEl, iEdFa ).first; // local edge number (in element)
                  // Loop number of Dof per edge
                  for ( ID l = 1; l <= _M_nDofpE[var]; ++l )
                    {
                      lDof = _M_nDofFV[var] + ( iEdFa - 1 ) * _M_nDofpE[var] + l ; // local Dof
                      gDof = _dof[var]->localToGlobal( iElAd, _M_nDofElemV[var] + ( iEdEl - 1 )
                          * _M_nDofpE[var] + l ); // global Dof
                      vidit = find( _fBdToIn[var].begin(), _fBdToIn[var].end(), gDof );
                      if ( vidit == _fBdToIn[var].end() )
                        { // the gDof has been encountered for the first time
                          bdltg( lDof ) = bDof[var];
                          _fBdToIn[var].push_back( gDof ); // local to boundary global on this face
                          bDof[var]++;
                        }
                      else
                        { // the gDof has been already inserted in the _fBdToIn vector
                          auxDof = ( ID ) ( vidit - _fBdToIn[var].begin() ) + 1;
                          bdltg( lDof ) = auxDof; // local to boundary global on this face
                        }
                    }
                }
            }
          // ===================================================
          // Face based Dof
          // ===================================================
          // Loop on number of Dof per face
          for ( ID l = 1; l <= _M_nDofpF[var]; ++l )
            {
              lDof = _M_nDofFE[var] + _M_nDofFV[var] + l; // local Dof
              gDof = _dof[var]->localToGlobal( iElAd, _M_nDofElemE[var] + _M_nDofElemV[var]
                  + ( iFaEl - 1 ) * _M_nDofpF[var] + l ); // global Dof
              vidit = find( _fBdToIn[var].begin(), _fBdToIn[var].end(), gDof );
              if ( vidit == _fBdToIn[var].end() )
                { // the gDof has been encountered for the first time
                  bdltg( lDof ) = bDof[var];
                  _fBdToIn[var].push_back( gDof ); // local to boundary global on this face
                  bDof[var]++;
                }
              else
                { // the gDof has been already inserted in the _fBdToIn vector
                  auxDof = ( ID ) ( vidit - _fBdToIn[var].begin() ) + 1;
                  bdltg( lDof ) = auxDof; // local to boundary global on this face
                }
            }

          _bdLtoG[var].push_back( bdltg );
        }
      }

    for(UInt var=0; var<_nvar; ++var) {
      // each processor holds information on HIS OWN patches
      _nBdDof[var] = _fBdToIn[var].size();

      _area[var].resize( _nBdDof[var] );
      for ( std::vector<Real>::iterator it = _area[var].begin();it<_area[var].end();it++ )
        *it = 0.0;

      _normal[var].resize( _nBdDof[var] * NDIM );
      for ( std::vector<Real>::iterator it = _normal[var].begin();it<_normal[var].end();it++ )
        *it = 0.0;

      _phi[var].resize( _nBdDof[var] );
      for ( std::vector<Real>::iterator it = _phi[var].begin();it<_phi[var].end();it++ )
        *it = 0.0;
    }
  }

  ///////////////////////////////////////////////


  // Area of faces with a certain marker
  template<typename Mesh>
  Real PostProc<Mesh>::area( const EntityFlag& flag )
  {
    // Each processor computes the area across his own flagged faces --> area_scatter
    // At the end I'll reduce the process areas --> area
    Real area_scatter(0.0), area(0.);

    std::list<ID> faces( _faces_with_flag[flag] );
    typedef std::list<ID>::iterator Iterator;

    //
    // Loop on flagged processor faces
    //
    for (Iterator j=faces.begin(); j != faces.end(); ++j) {

      _feBd[0]->updateMeas( _mesh->bElement( *j ) );  // updating finite element information

      area_scatter += _feBd[0]->measure();

    }

    // reducing per-processor information
    _epetraMap->Comm().SumAll( &area_scatter, &area, 1 );

    return area;
  }


  // flux of vector field "sol" through faces with a certain marker
  template<typename Mesh>
  template<typename V>
  Real PostProc<Mesh>::flux( const V& sol, const EntityFlag& flag, UInt feSpace,
      UInt nDim )
  {
    // Each processor computes the flux across his own flagged faces --> flux_scatter
    // At the end I'll reduce the process fluxes --> flux
    Real flux_scatter(0.0), flux(0.);

    // I need the global Dof ID to query the vector
    // idLocalDof is the id of the dof in the data structure of PostProc class
    // idGlobalDof is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt idLocalDof, idGlobalDof;

    // list of flagged faces on current processor
    std::list<ID> faces( _faces_with_flag[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of sol in the current face
    Vector sol_localDof(nDim * _M_nDofF[feSpace]);

    // Loop on faces
    for (Iterator j=faces.begin(); j != faces.end(); ++j) {

      // Updating quadrature data on the current face
      _feBd[feSpace]->updateMeasNormalQuadPt(_mesh->bElement(*j));

      // Quadrature formula
      // Loop on quadrature points
      for(int iq=0; iq< _feBd[feSpace]->nbQuadPt; ++iq) {

        // Dot product
        // Loop on components
        for (UInt ic =0; ic<nDim; ++ic) {

          // Interpolation
          // Loop on local dof
          for (ID idofF=1; idofF<=_M_nDofF[feSpace]; ++idofF) {

            // Extracting nodal values of sol in the current face
            idLocalDof = _bdLtoG[feSpace][ ( UInt ) *j - 1 ][ idofF - 1 ];
            idGlobalDof = _fBdToIn[feSpace][idLocalDof-1]; // this is in the GLOBAL mesh

            sol_localDof[ic*_M_nDofF[feSpace]+idofF-1] = sol[ic*_M_nTotalDof[feSpace]+idGlobalDof];

            flux_scatter += _feBd[feSpace]->weightMeas(iq)
            * sol_localDof[ic*_M_nDofF[feSpace]+idofF-1]
                      * _feBd[feSpace]->phi(int(idofF-1),iq)
                      * _feBd[feSpace]->normal(int(ic),iq);
          }
        }
      }
    }
    // Reducing per-processor values
    _epetraMap->Comm().SumAll( &flux_scatter, &flux, 1 );

    return flux;
  }


  // Average value of sol on faces with a certain marker
  template<typename Mesh>
  template<typename V>
  Vector PostProc<Mesh>::average( const V& sol, const EntityFlag& flag,
      UInt feSpace, UInt nDim )
  {
    // Each processor computes the average value on his own flagged faces --> sol_avg_scatter
    // At the end I'll reduce the process values --> sol_avg
    Vector sol_avg_scatter(nDim), sol_avg(nDim), sol_face(nDim);
    // basic policy for type V: operator[] available
    for( UInt ic=0; ic < nDim; ++ic ) {
      sol_avg_scatter[ic] = 0.; sol_avg[ic] = 0.; sol_face[ic] = 0.;
    }

    // The total area of the considered faces
    Real area_scatter(0.), area;

    // I need the global Dof ID to query the Oseen solution vector
    // idLocalDof is the id of the dof in the data structure of PostProc class
    // idGlobalDof is the corresponding ID in the GLOBAL mesh (prior to partitioning)
    UInt idLocalDof, idGlobalDof;

    // list of flagged faces on current processor
    std::list<ID> faces( _faces_with_flag[flag] );
    typedef std::list<ID>::iterator Iterator;

    // Nodal values of sol in the current face
    std::vector<Real> sol_localDof(_M_nDofF[feSpace]);

    // Loop on faces
    for (Iterator j=faces.begin(); j != faces.end(); ++j) {

      // basic policy for type V: operator[] available
      for( UInt ic=0; ic < nDim; ++ic ) sol_face[ic] = 0.;

      // Updating quadrature data on the current face
      _feBd[feSpace]->updateMeasNormalQuadPt(_mesh->bElement(*j));

      // Loop on components
      for (UInt ic =0; ic<nDim; ++ic) {

        // Quadrature formula
        // Loop on quadrature points
        for(int iq=0; iq< _feBd[feSpace]->nbQuadPt; ++iq) {

          // Interpolation
          // Loop on local dof
          for (ID idofF=1; idofF<=_M_nDofF[feSpace]; ++idofF) {

            // Extracting nodal values of sol in the current face
            idLocalDof = _bdLtoG[feSpace][ ( UInt ) *j - 1 ][ idofF - 1 ];
            idGlobalDof = _fBdToIn[feSpace][idLocalDof-1]; // this is in the GLOBAL mesh

            // basic policy for type V: operator[] available
            sol_localDof[idofF-1] = sol[ic*_M_nTotalDof[feSpace]+idGlobalDof];

            sol_face[ic] += _feBd[feSpace]->weightMeas(iq)
            * sol_localDof[idofF-1] * _feBd[feSpace]->phi(int(idofF-1),iq);
          }
        }
        // Computing the sol integral over the boundary faces
        sol_avg_scatter[ic] += sol_face[ic];
      }

      // Computing the area
      area_scatter += _feBd[feSpace]->measure();
    }

    _epetraMap->Comm().SumAll( &area_scatter, &area, 1 );

    // Reducing per-processor values
    for( UInt ic=0; ic < nDim; ++ic ) {
//      sol_avg_scatter[ic] /= area;
      _epetraMap->Comm().SumAll( &sol_avg_scatter[ic], &sol_avg[ic], 1 );
    }

    return sol_avg / area;
  }


  // Area of patches on the boundary
  template<typename Mesh>
  void PostProc<Mesh>::set_area()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      // area of the mesh face
      Real loc_area;

      // ID of the considered dof in this class' vectors
      ID idLocalDof;

      // ===================================================
      // Loop on boundary faces
      // ===================================================
      for ( ID ibF = 1 ; ibF <= _M_bdnF; ++ibF )
        {

          _feBd[var]->updateMeas( _mesh->bElement( ibF ) );  // updating finite element information

          loc_area = _feBd[var]->measure();
          // Loop on the total Dof per Face
          for ( ID idofF = 1; idofF <= _M_nDofF[var]; ++idofF )
            {
              // Extracting local ID of idofF
              idLocalDof = _bdLtoG[var][ ibF - 1 ][ idofF - 1 ];
              _area[var][idLocalDof-1] += loc_area;
            }
        }
    }
  }

  template <typename Mesh>
  void PostProc<Mesh>::show_area()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      std::cout << "\n***** Post Proc: Area of the patches *****"
      << "\n\tfor variable " << var << std::endl;
      ID count = 1;

      for ( std::vector<Real>::iterator it = _area[var].begin(); it<_area[var].end(); it++ )
        {
          std::cout << "Boundary Dof: " << count
          << ", corresponding to Global Dof: " << _fBdToIn[var][ count - 1 ]
          << " has patch area: " << *it << std::endl;
          count++;
        }
    }
  }

  template<typename Mesh>
  void PostProc<Mesh>::show_bdLtoG()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      int count = 0;
      std::cout << "\n***** Post Proc: Bd Local To Global *****"
      << "\n\tfor variable " << var << std::endl;
      std::cout << _bdLtoG[var].size() << std::endl;
      for ( std::vector<SimpleVect<ID> >::iterator it1 = _bdLtoG[var].begin();
      it1<_bdLtoG[var].end();it1++ )
        {
          count++;
          std::cout << "Bd Face " << count << std::endl;
          for ( SimpleVect<ID>::iterator it2 = it1->begin();it2<it1->end();it2++ )
            {
              std::cout << *it2 << ",";
            }
          std::cout << std::endl;
        }

      std::cout << "***** Post Proc: From Boundary Faces to Global Dof *****" << std::endl;
      std::cout << _fBdToIn[var].size() << std::endl;

      for ( std::vector<ID>::iterator it3 = _fBdToIn[var].begin();
      it3<_fBdToIn[var].end(); it3++ )
        {
          std::cout << "Index :" << it3 - _fBdToIn[var].begin()
          << ", Global Dof: " << *it3 << std::endl;
        }
    }
  }

  /////////////////////////////////////////////////

  ///////////////////////////////////////////////


  // Normal vectors of patches on the boundary
  template<typename Mesh>
  void PostProc<Mesh>::set_normal()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      // for each patch, the average of each component of the normal vector
      Real sum;

      // ID of the considered dof in this class' vectors
      ID idLocalDof;

      // ===================================================
      // Loop on boundary faces
      // ===================================================
      for ( ID ibF = 1 ; ibF <= _M_bdnF; ++ibF )
        {
          // updating finite element information
          _feBd[var]->updateMeasNormal( _mesh->bElement( ibF ) );

          // Loop on the components
          for ( int icomp = 0; icomp < NDIM;icomp++ )
            {
              sum = 0.;
              // Loop on the quadrature points
              for ( int l = 0; l < _feBd[var]->nbQuadPt; ++l )
                {
                  sum += _feBd[var]->normal( icomp, l ) * _feBd[var]->weightMeas( l );
                }
              for ( ID idofF = 1; idofF <= _M_nDofF[var]; ++idofF )
                {
                  // Extracting local ID of idofF
                  idLocalDof = _bdLtoG[var][ ibF - 1 ][ idofF - 1 ];
                  _normal[var][ icomp * _nBdDof[var] + idLocalDof - 1 ] += sum;
                }
            }
        }
      // Normalization of the averaged normals with the patch area
      for ( UInt inorm = 0; inorm < _nBdDof[var]; ++inorm )
        {
          Real loc_area = _area[var][ inorm ];
          for ( int icc = 0;icc < NDIM;icc++ )
            _normal[var][ icc * _nBdDof[var] + inorm ] *= 1. / loc_area;
        }
    }
  }


  template <typename Mesh>
  void PostProc<Mesh>::show_normal()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      std::cout << "\n***** Post Proc: Normal vector on the patches *****"
      << "\n\tfor variable " << var << std::endl;

      ID count = 1;

      for ( std::vector<Real>::iterator it = _area[var].begin();it<_area[var].end();it++ )
        {
          std::cout << "Boundary Dof: " << count
          << ", corresponding to Global Dof: " << _fBdToIn[var][ count - 1 ]
          << " has patch area: " << *it << std::endl;
          std::cout << "and normal components " ;
          for ( int icomp = 0; icomp<NDIM; icomp++ )
            std::cout << _normal[var][ icomp * _nBdDof[var] + count - 1 ] << " ";

          std::cout << std::endl;
          count++;
        }
    }

    std::cout << "End SHOW NORMAL" << std::endl;
  }

  //////////////////////////////////////////////////
  //
  ///////////////////////////////////////////////////

  // Vector with the integral of the shape functions on the patches on the boundary
  template<typename Mesh>
  void PostProc<Mesh>::set_phi()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      // sum contributions from each face of the patch
      Real sum;

      // ID of the considered dof in this class' vectors
      ID idLocalDof;

      // ===================================================
      // Loop on boundary faces
      // ===================================================
      for ( ID ibF = 1 ; ibF <= _M_bdnF; ++ibF )
        {

          _feBd[var]->updateMeas( _mesh->bElement( ibF ) );  // updating finite element information

          for ( ID idofF = 1; idofF <= _M_nDofF[var]; ++idofF )
            {
              sum = 0.0;
              // global dof
              idLocalDof = _bdLtoG[var][ ( UInt ) ibF - 1 ][ idofF - 1 ];

              // Loop on the quadrature points
              for ( int l = 0; l < _feBd[var]->nbQuadPt; ++l )
                {
                  sum += _feBd[var]->phi( ( int ) ( idofF - 1 ), l ) * _feBd[var]->weightMeas( l );
                }
              _phi[var][ idLocalDof - 1 ] += sum;
            }
        }
    }
  }


  template <typename Mesh>
  void PostProc<Mesh>::show_phi()
  {
    for( UInt var=0; var<_nvar; ++var ) {

      std::cout << "\n***** Post Proc: Average phi on the patches *****"
      << "\n\tfor variable " << var << std::endl;

      ID count = 0;

      for ( std::vector<Real>::iterator it = _area[var].begin();it<_area[var].end();it++ )
        {
          std::cout << "Boundary Dof: " << count + 1
          << ", corresponding to Global Dof: " << _fBdToIn[var][ count - 1 ]
          << " has patch area: " << *it << std::endl;
          std::cout << "and average phi  " << _phi[var][ count ] << std::endl ;
          count++;
        }
    }
    std::cout << "End SHOW PHI" << std::endl;
  }

// the following part is definitely not ready yet TP 09/2008
#if 0
  ////////////////////////////////
  ////////////////////////////////
  ////////////////////////////////
  template<typename Mesh>
  template<typename V>
  V PostProc<Mesh>::compute_sstress( const V& r, UInt ncomp, bool residual = true ) const
  {
    ASSERT( ncomp = NDIM, "Error: Shear stress computation possible only for vector unknowns" );

    // prepare the vectors
    V stress( _nBdDof[0] * NDIM ); stress.clear();
    V nstress( _nBdDof[0] * NDIM ); nstress.clear();
    V sstress( _nBdDof[0] * NDIM ); sstress.clear();
    ID count, idGlobalDof;

    // number of DOFs for each component
    UInt dim = r.size() / ncomp;

    // helper structures to avoid "if" statement
    // if residual==true, vector r has ncomp*_nBdDof components
    // if residual==false, vector r has ncomp*_M_nTotalDof components
    Vector coef(2);
    std::vector<ID> index(2);

    // loop on locally stored DOFs
    for ( count = 0; count < _nBdDof[0]; ++count )
      {
        // find global ID for boundary dof
        idGlobalDof = _fBdToIn[0][count]; // this is in the GLOBAL mesh

        coef[0] = 1.; index[0] = count;
        coef[1] = 1./ _phi[0][count]; index[1] = idGlobalDof - 1;

        for ( UInt ind_comp = 0; ind_comp < ncomp; ++ind_comp ) {
          stress[ count + ind_comp * _nBdDof[0] ] = //damned conventions trouble : 0 or 1
            coef[residual] * r[ index[residual] + ind_comp * dim ];
        }
      }

    Real sn = 0.;
    ///// Normal stress
    for ( count = 0; count < _nBdDof[0]; count++ )
      {
        sn = 0.;
        for ( UInt ind_comp = 0;ind_comp < ncomp;ind_comp++ )
          sn += stress[ count + ind_comp * _nBdDof[0] ] * _normal[0][ count + ind_comp * _nBdDof[0] ];

        for ( UInt ind_comp = 0;ind_comp < ncomp;ind_comp++ )
          nstress[ count + ind_comp * _nBdDof[0] ] = sn * _normal[0][ count + ind_comp * _nBdDof[0] ];
      }

    // Shear Stress: this vector lives on the patches (ncomp*_nBdDof components)
    sstress = stress - nstress;

    return sstress;
  } // compute_sstress



  ////////////////////////////////
  ////////////////////////////////
  ///////////////////////////////
  template<typename Mesh>
  template<typename V>
  V PostProc<Mesh>::compute_stress( const V& grad, const UInt& dim,
      const Real& visc ) const
      {
        V stress( _nBdDof[0] * NDIM ); stress.clear();

        ID idGlobalDof;

        // cycle over boundary dof
        for ( ID bound_dof = 0; bound_dof < _nBdDof[0]; ++bound_dof )
          {
            // find global ID for boundary dof
            idGlobalDof = _fBdToIn[0][bound_dof]; // this is in the GLOBAL mesh

            // cycle over stress components
            for( UInt scomp=0; scomp<NDIM; ++scomp )

              // cycle over normal components
              for ( UInt ncomp = 0; ncomp < NDIM; ++ncomp ) {
                stress[ bound_dof + scomp*_nBdDof[0] ] += visc *
                // grad!
                ( grad[ idGlobalDof-1 + ncomp*dim + scomp*nDimensions*dim ] +
                    // transpose grad!
                    grad[ idGlobalDof-1 + scomp*dim + ncomp*nDimensions*dim ] )
                    * _normal[ bound_dof + ncomp*_nBdDof[0] ];
              }
          }


        return stress;
      } // compute_stress
#endif
} // namespace LifeV

#endif /* _POST_PROC_H */

