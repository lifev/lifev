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
  \file NavierStokesAleHandler.h
  \author M.A. Fernandez
  \date 01/2003
  \version 1.0

  \brief This file contains an abstract class for NavierStokes solvers in moving domains.
  An ALE formulation is used

*/

#ifndef _NAVIERSTOKESALEHANDLER_H_
#define _NAVIERSTOKESALEHANDLER_H_
#include "NavierStokesHandler.hpp"
#include "meshMotion.hpp"

namespace LifeV
{
/*!
  \class NavierStokesAleHandler

  Abstract class which defines the general structure of a NavierStokes solver with moving domain.
  For each new NavierStokes solver  we have to implement the corresponding timeAdvance and an iterate methods.

  \note This class inherits from NavierStokesHandler, new methods concern mesh motion staff

*/
template <typename Mesh>
class NavierStokesAleHandler:
            public NavierStokesHandler<Mesh>,
            public HarmonicExtension
{
public:

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
      \param BCh_mesh boundary conditions for the motion harmonic extension
     */
    NavierStokesAleHandler( const GetPot& data_file,
                            const RefFE& refFE_u,
                            const RefFE& refFE_p,
                            const QuadRule& Qr_u,
                            const QuadRule& bdQr_u,
                            const QuadRule& Qr_p,
                            const QuadRule& bdQr_p,
                            BCHandler& BCh_u,
                            BCHandler& BCh_mesh );

    //! Do nothing destructor
    virtual ~NavierStokesAleHandler()
    {}

    //! Interpolated mesh velocity
    PhysVectUnknown<Vector>& wInterpolated();

    //! Interpolated mesh velocity
    PhysVectUnknown<Vector>& w();

    //! Interpolated mesh velocity
    PhysVectUnknown<Vector>& dwInterpolated();

    //! Updating mesh
    void updateMesh( const Real& time );

    void updateMeshTransp( const Real& time );

    void updateDispVelo();

    //! Postprocessing
    void postProcess();

protected:
    //! The previous extension of the displacement
    PhysVectUnknown<Vector> _dispOld;

    //! The mesh velocity
    PhysVectUnknown<Vector> _w;

    //! The mesh velocity
    PhysVectUnknown<Vector> _wInterp;

    //! The interpolated displacement dervivative (right hand sized linearized)
    PhysVectUnknown<Vector> _dInterp;

    //! The interpolated mesh velocity mesh velocity (right hand sized linearized)
    PhysVectUnknown<Vector> _dwInterp;


    //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
    void _interpolate( const UInt nbcomp, const Vector& w, Vector& wInterp );
};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
NavierStokesAleHandler<Mesh>::
NavierStokesAleHandler( const GetPot& data_file, const RefFE& refFE_u,
                        const RefFE& refFE_p, const QuadRule& Qr_u, const QuadRule& bdQr_u,
                        const QuadRule& Qr_p, const QuadRule& bdQr_p, BCHandler& BCh_u, BCHandler& BCh_mesh ) :
        NavierStokesHandler<Mesh>( data_file,
                                   refFE_u,
                                   refFE_p,
                                   Qr_u,
                                   bdQr_u,
                                   Qr_p,
                                   bdQr_p,
                                   BCh_u ),
        HarmonicExtension( _mesh,
                           1.0,
                           quadRuleTetra4pt,
                           quadRuleTria3pt,
                           BCh_mesh ),
        _dispOld( _dof_mesh.numTotalDof() ),
        _w( _dof_mesh.numTotalDof() ),
        _wInterp( _dim_u ),
        _dInterp( _dim_u ),
        _dwInterp( _dim_u )
{
    _dispOld = ZeroVector( _dispOld.size() );
    _w = ZeroVector( _w.size() );
    _wInterp = ZeroVector( _dim_u );
    _dInterp = ZeroVector( _dim_u );
    _dwInterp = ZeroVector( _dim_u );
}


// Mesh and grid velocity update
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateMesh( const Real& time )
{
    // Updating mesh displacement and velocity
    updateExtension( _mesh, time );

    Real dti = 1.0 / _dt;
    _w = ( _disp - _dispOld ) * dti;
    _interpolate( _w.nbcomp(), _w, _wInterp );
    // Updating mesh points
    _mesh.moveMesh( _disp );
}


// Mesh and grid velocity update
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateMeshTransp( const Real& time )
{

    updateExtensionTransp( _mesh, time );
    Real dti = 1.0 / _dt;
    _w = ( _disp - _dispOld ) * dti;
    _interpMeshVelocity();
}


//  Updating variations for grid velocity and displacement
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateDispVelo()
{
    // Updating mesh displacement and velocity
    updateExtension( _mesh, 0.0, 1 );

    Real dti = 1.0 / _dt;

    std::cout << " max norm dx = " << norm_inf( _disp ) << std::endl;

    _interpolate( _w.nbcomp(), _disp, _dInterp );

    std::cout << " max norm dxInterp = " << norm_inf( _dInterp ) << std::endl;

    _dwInterp = _dInterp * dti;

    std::cout << " max norm dwInterp = " << norm_inf( _dwInterp ) << std::endl;
}


// Postprocessing pressure
template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::postProcess()
{

    std::ostringstream index;
    std::string name;

    ++_count;
    if ( fmod( float( _count ), float( _verbose ) ) == 0.0 )
    {
        std::cout << "  o-  Post-processing \n";
        index << ( _count / _verbose );

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


        wr_medit_ascii_scalar( "press." + name + ".bb", _p.giveVec(), _p.size() );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", _u.giveVec(), _mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", _u.giveVec() + _dim_u, _mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", _u.giveVec() + 2 * _dim_u, _mesh.numVertices() );
        //   wr_medit_ascii_scalar("velw_x."+name+".bb",_w.giveVec(),_mesh.numVertices());
        // wr_medit_ascii_scalar("velw_y."+name+".bb",_w.giveVec() + _mesh.numVertices(),_mesh.numVertices());
        //wr_medit_ascii_scalar("velw_z."+name+".bb",_w.giveVec() + 2*_mesh.numVertices(),_mesh.numVertices());


        //wr_medit_ascii("press."+name+".mesh", _mesh);
        wr_medit_ascii( "press." + name + ".mesh", _mesh, _disp, 12 );
        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(),_mesh.numVertices(),_dim_u);
        system( ( "ln -s press." + name + ".mesh vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s press." + name + ".mesh vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s press." + name + ".mesh vel_z." + name + ".mesh" ).data() );
        //system(("ln -s press."+name+".mesh velw_x."+name+".mesh").data());
        //system(("ln -s press."+name+".mesh velw_y."+name+".mesh").data());
        //system(("ln -s press."+name+".mesh velw_z."+name+".mesh").data());

        // system(("ln -s "+_mesh_file+" veloc."+name+".mesh").data());
    }
}




// The interpolated mesh velocity
template <typename Mesh>
PhysVectUnknown<Vector>& NavierStokesAleHandler<Mesh>::
wInterpolated()
{
    return _wInterp;
}

// The interpolated mesh velocity
template <typename Mesh>
PhysVectUnknown<Vector>& NavierStokesAleHandler<Mesh>::
dwInterpolated()
{
    return _dwInterp;
}


// The mesh velocity
template <typename Mesh>
PhysVectUnknown<Vector>& NavierStokesAleHandler<Mesh>::w()
{
    return _w;
}


// This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
_interpolate( const UInt nbComp, const Vector& w, Vector& wInterp )
{

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = _refFE_u.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = _refFE_u.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = _refFE_u.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = _refFE_u.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElem = _mesh.getRefFE().nbDof; // Number of local dof per element of the _mesh (_mesh.getRefFE().nbDof)

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element


    Real x, y, z;

    KN<Real> wLoc( nDofElem * nbComp );

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= _mesh.numVolumes(); ++iElem )
    {

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            for ( ID idof = 0; idof < nDofElem; ++idof )
                wLoc( icmp * nDofElem + idof ) = w( icmp * _dof_mesh.numTotalDof() + _dof_mesh.localToGlobal( iElem, idof + 1 ) - 1 );

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
                    x = _refFE_u.xi( lDof - 1 );
                    y = _refFE_u.eta( lDof - 1 );
                    z = _refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the element
                            __sum += wLoc( icmp * nDofElem + idof ) * _mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolated mesh velocity
                        wInterp( icmp * _dof_u.numTotalDof() + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
                    x = _refFE_u.xi( lDof - 1 );
                    y = _refFE_u.eta( lDof - 1 );
                    z = _refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )   // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * _mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        wInterp( icmp * _dof_u.numTotalDof() + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
                    x = _refFE_u.xi( lDof - 1 );
                    y = _refFE_u.eta( lDof - 1 );
                    z = _refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * _mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        wInterp( icmp * _dof_u.numTotalDof() + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
            x = _refFE_u.xi( lDof - 1 );
            y = _refFE_u.eta( lDof - 1 );
            z = _refFE_u.zeta( lDof - 1 );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {

                // Interpolating data at the nodal point
                double __sum = 0;
                for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the adjacent element
                    __sum += wLoc( icmp * nDofElem + idof ) * _mesh.getRefFE().phi( idof, x, y, z );

                // Updating interpolating vector
                wInterp( icmp * _dof_u.numTotalDof() + _dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
            }
        }
    }
}
}
#endif
