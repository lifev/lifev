/*
 This file is part of the LifeV library
 Copyright (C) 2004 EPFL, INRIA and Politechnico di Milano

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
  \file convDiffReactHandler.h
  \author M. Prosi
  \date 03/2004
  \version 1.0

  \brief This file contains an abstract class for the Convection-Diffusion-Reaktion equation  solvers.

*/

#ifndef _CONVDIFFREACTHANDLER_H_
#define _CONVDIFFREACTHANDLER_H_


#include "lifeV.hpp"
#include "dataConvDiffReact.hpp"
#include "geoMap.hpp"
#include "dataAztec.hpp"
#include "refFE.hpp"
#include "dof.hpp"
#include "medit_wrtrs.hpp"
#include "bcCond.hpp"
#include "bdf.hpp"
#include "post_proc.hpp"
#include "openDX_wrtrs.hpp"
#include <cmath>
#include <sstream>


namespace LifeV
{
/*!
  \class convDiffReactHandler

  Abstract class which defines the general structure of a Convection-Diffusion-Reaction solver.
  For each new Convection-Diffusion-Reaction solver  we have to implement the corresponding
  timeAdvance and an iterate methods

*/

template <typename Mesh>
class ConvDiffReactHandler:
            public DataConvDiffReact<Mesh>
{

public:

    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE_c reference FE for the concentration
      \param Qr_c volumic quadrature rule for the concentration
      \param bdQr_c surface quadrature rule for the concentration
      \param BCh_c boundary conditions for the concentration
      \param ord_bdf order of the bdf time advancing scheme (default: Backward Euler)
    */
    ConvDiffReactHandler( const GetPot& data_file, const RefFE& refFE_c,
                          const QuadRule& Qr_c, const QuadRule& bdQr_c, BCHandler& BCh_c );

    //! Sets initial condition for the concentration (incremental approach): the initial time is t0, the time step dt
    void initialize( const Function& c0, Real t0, Real dt );

    //! Sets initial condition for the concentration from file
    void initialize( const std::string & vname );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    virtual void timeAdvance( const Function source, const Real& time ) = 0;

    //! Update convective term, bc treatment and solve the linearized cdr system
    virtual void iterate( const Real& time ) = 0;

    //! Returns the concentration vector
    ScalUnknown<Vector>& c();

    //! Returns the velocity on the concentration nodes
    PhysVectUnknown<Vector>& u_c();

    //! Returns the concentration Dof
    const Dof& cDof() const;

    //! Returns the BDF Time Advancing stuff
    const Bdf& bdf() const;

    //! Do nothing destructor
    virtual ~ConvDiffReactHandler()
    {}


protected:

    //! Reference FE for the concentration
    const RefFE& _refFE_c;

    //! The Dof object associated with the concentration
    Dof _dof_c;

    //! The number of total concentration dofs
    UInt _dim_c;

    //! Quadrature rule for concentration volumic elementary computations
    const QuadRule& _Qr_c;

    //! Quadrature rule for concentration surface elementary computations
    const QuadRule& _bdQr_c;

    //! Current FE for the concentration c
    CurrentFE _fe_c;
    CurrentBdFE _feBd_c;

    //! The concenration
    ScalUnknown<Vector> _c;

    //! velocity vector on the concentration nodes
    PhysVectUnknown<Vector> _u_c;


    //! Structure to hold the interpolation values of concentration nodes in the velocity grid
    struct intpolcoord
    {
        Real b[ 4 ];
        ID ele;
    };
    std::vector<intpolcoord> _u_to_c;

    //! The BC handler
    BCHandler& _BCh_c;

    // ! The BDF Time Advance Method
    Bdf _bdf;

};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
ConvDiffReactHandler<Mesh>::
ConvDiffReactHandler( const GetPot& data_file, const RefFE& refFE_c,
                      const QuadRule& Qr_c, const QuadRule& bdQr_c, BCHandler& BCh_c ) :
        DataConvDiffReact<Mesh>( data_file ),
        _refFE_c( refFE_c ),
        _dof_c( _mesh, _refFE_c ),
        _dim_c( _dof_c.numTotalDof() ),
        _Qr_c( Qr_c ),
        _bdQr_c( bdQr_c ),
        _fe_c( _refFE_c, getGeoMap( _mesh ), _Qr_c ),
        _feBd_c( _refFE_c.boundaryFE(), getGeoMap( _mesh ).boundaryMap(), _bdQr_c ),
        _c( _dim_c ),
        _u_c( _dim_c ),
        _BCh_c( BCh_c ),
        _bdf( _order_bdf )
{}


// Returns the concentration
template <typename Mesh>
ScalUnknown<Vector>&
ConvDiffReactHandler<Mesh>::c()
{
    return _c;
}

// Returns the velocity
template <typename Mesh>
PhysVectUnknown<Vector>&
ConvDiffReactHandler<Mesh>::u_c()
{
    return _u_c;
}


// Returns the concentration Dof
template <typename Mesh>
const Dof&
ConvDiffReactHandler<Mesh>::cDof() const
{
    return _dof_c;
}


// Returns the BDF Time Advancing stuff
template <typename Mesh>
const Bdf&
ConvDiffReactHandler<Mesh>::bdf() const
{
    return _bdf;
}

// ! Initialize when  initial conditions concentration
template <typename Mesh>
void
ConvDiffReactHandler<Mesh>::initialize( const Function& c0, Real t0, Real dt )
{

    _bdf.initialize_unk( c0, _mesh, _refFE_c, _fe_c, _dof_c, t0, dt, 1 );
    _c = *( _bdf.unk().begin() );

    _bdf.showMe();

}

// ! Initialize when initial values for the concentration are read from file
template <typename Mesh>
void
ConvDiffReactHandler<Mesh>::initialize( const std::string & vname )
{


    std::fstream Resfile( vname.c_str(), std::ios::in | std::ios::binary );
    if ( Resfile.fail() )
    {
        std::cerr << " Error in initialize: File not found or locked" << std::endl;
        abort();
    }
    Resfile.read( ( char* ) & _c( 1 ), _c.size() * sizeof( double ) );
    Resfile.close();

    _bdf.initialize_unk( _c );

    _bdf.showMe();

}
}
#endif
