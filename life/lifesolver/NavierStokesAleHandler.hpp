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
#include <life/lifesolver/NavierStokesHandler.hpp>
#include <life/lifefem/meshMotion.hpp>
#include <life/lifecore/life.hpp>
#include <boost/function.hpp>

namespace LifeV
{
/*!
  \class NavierStokesAleHandler

  Abstract class which defines the general structure of a NavierStokes solver with moving domain.
  For each new NavierStokes solver  we have to implement the corresponding timeAdvance and an iterate methods.

  \note This class inherits from NavierStokesHandler, new methods concern mesh motion staff

*/
template <typename Mesh>
class NavierStokesAleHandler:public NavierStokesHandler<Mesh>
//    ,
//                             public HarmonicExtension
{
public:

    typedef typename NavierStokesHandler<Mesh>::source_type source_type;

    typedef typename NavierStokesHandler<Mesh>::Function Function;

    //! Constructors
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

    NavierStokesAleHandler( const GetPot&   data_file,
                            const RefFE&    refFE_u,
                            const RefFE&    refFE_p,
                            const QuadRule& Qr_u,
                            const QuadRule& bdQr_u,
                            const QuadRule& Qr_p,
                            const QuadRule& bdQr_p);

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
                     const std::string& pressName,
                     const std::string& velwName,
                     double             startT = 0.);
    //! Postprocessing
    void postProcess();

    //! HarmonicExtension getter
    const HarmonicExtension& harmonicExtension() {return M_harmonicExtension;}
    //! HarmonicExtension BC setter
    void setHarmonicExtensionBC( BCHandler &bc_he)
        {M_harmonicExtension.setHarmonicExtensionBC(bc_he);}
private:
    HarmonicExtension    M_harmonicExtension;

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


    //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > this->_mesh.getRefFE().nbNodes)
    void _interpolate( const UInt nbcomp, const Vector& w, Vector& wInterp );

private:

    Real _factor; // amplification factor for deformed mesh
    void readUnknown( const std::string       &name,
                      PhysVectUnknown<Vector> &unknown);

};



//
// IMPLEMENTATION
//


// Constructors
template <typename Mesh>
NavierStokesAleHandler<Mesh>::
NavierStokesAleHandler( const GetPot&   data_file,
                        const RefFE&    refFE_u,
                        const RefFE&    refFE_p,
                        const QuadRule& Qr_u,
                        const QuadRule& bdQr_u,
                        const QuadRule& Qr_p,
                        const QuadRule& bdQr_p,
                        BCHandler&      BCh_u,
                        BCHandler&      BCh_mesh ) :
        NavierStokesHandler<Mesh>( data_file,
                                   refFE_u,
                                   refFE_p,
                                   Qr_u,
                                   bdQr_u,
                                   Qr_p,
                                   bdQr_p,
                                   BCh_u ),
        M_harmonicExtension        ( data_file,
                                     this->_mesh,
                                     1.0,
                                     quadRuleTetra4pt,
                                     quadRuleTria3pt,
                                     BCh_mesh ),
        _dispOld                 ( M_harmonicExtension.dofMesh().numTotalDof() ),
        _w                       ( M_harmonicExtension.dofMesh().numTotalDof() ),
        _wInterp                 ( this->_dim_u ),
        _dInterp                 ( this->_dim_u ),
        _dwInterp                ( this->_dim_u )
{
    _factor   = data_file ( "fluid/miscellaneous/factor", 1.0 );
    _dispOld  = ZeroVector( _dispOld.size() );
    _w        = ZeroVector( _w.size() );
    _wInterp  = ZeroVector( _wInterp.size() );
    _dInterp  = ZeroVector( _dInterp.size() );
    _dwInterp = ZeroVector( _dwInterp.size() );
}


template <typename Mesh>
NavierStokesAleHandler<Mesh>::
NavierStokesAleHandler( const GetPot&   data_file,
                        const RefFE&    refFE_u,
                        const RefFE&    refFE_p,
                        const QuadRule& Qr_u,
                        const QuadRule& bdQr_u,
                        const QuadRule& Qr_p,
                        const QuadRule& bdQr_p ) :
        NavierStokesHandler<Mesh>( data_file,
                                   refFE_u,
                                   refFE_p,
                                   Qr_u,
                                   bdQr_u,
                                   Qr_p,
                                   bdQr_p ),
        M_harmonicExtension        ( data_file,
                                     this->_mesh,
                                     1.0,
                                     quadRuleTetra4pt,
                                     quadRuleTria3pt ),
        _dispOld                 ( M_harmonicExtension.dofMesh().numTotalDof() ),
        _w                       ( M_harmonicExtension.dofMesh().numTotalDof() ),
        _wInterp                 ( this->_dim_u ),
        _dInterp                 ( this->_dim_u ),
        _dwInterp                ( this->_dim_u )
{
    _factor   = data_file ( "fluid/miscellaneous/factor", 1.0 );
    _dispOld  = ZeroVector( _dispOld.size() );
    _w        = ZeroVector( _w.size() );
    _wInterp  = ZeroVector( _wInterp.size() );
    _dInterp  = ZeroVector( _dInterp.size() );
    _dwInterp = ZeroVector( _dwInterp.size() );
}


// Mesh and grid velocity update
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateMesh( const Real& time )
{
    // Updating mesh displacement and velocity
    std::cout << "Updating the harmonic extension mesh ... " << std::flush;
    M_harmonicExtension.updateExtension( this->_mesh, time );

    Real dti = 1.0 / this->_dt;
    _w = ( M_harmonicExtension.getDisplacement() - _dispOld ) * dti;
    _interpolate( _w.nbcomp(), _w, _wInterp );
    // Updating mesh points
    this->_mesh.moveMesh( M_harmonicExtension.getDisplacement() );
    std::cout << "ok." << std::flush << std::endl;
}


// Mesh and grid velocity update
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateMeshTransp( const Real& time )
{

    M_harmonicExtension.updateExtensionTransp( this->_mesh, time );
    Real dti = 1.0 / this->_dt;
    _w = ( M_harmonicExtension.getDisplacement() - _dispOld ) * dti;
    this->_interpMeshVelocity();
}


//  Updating variations for grid velocity and displacement
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
updateDispVelo()
{
    // Updating mesh displacement and velocity
    M_harmonicExtension.updateExtension( this->_mesh, 0.0, 1 );

    Real dti = 1.0 / this->_dt;

    std::cout << " max norm dx = " << norm_inf( M_harmonicExtension.getDisplacement() ) << std::endl;

    _interpolate( _w.nbcomp(), M_harmonicExtension.getDisplacement(), _dInterp );

    std::cout << " max norm dxInterp = " << norm_inf( _dInterp ) << std::endl;

    _dwInterp = _dInterp * dti;

    std::cout << " max norm dwInterp = " << norm_inf( _dwInterp ) << std::endl;
}

template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::initialize( const Function& u0 )
{
    NavierStokesHandler<Mesh>::initialize(u0);
}


template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::initialize( const Function& u0, const Function& p0,
                                          Real t0, Real dt )
{
    NavierStokesHandler<Mesh>::initialize(u0, p0, t0, dt);
}

template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::initialize(  const std::string & vname )
{
     NavierStokesHandler<Mesh>::initialize( vname );
}



//! Initialize the fluid with a solution file
//! written in MEDIT format
template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::initialize( const std::string& velName,
                                          const std::string& pressName,
                                          const std::string& velwName,
                                          double             startT)
{
    std::cout << "  F- restarting at time = " << startT << " " << this->_dt << std::endl;

    this->_count = (UInt) (startT/this->_dt - 0.5);

    NavierStokesHandler<Mesh>::initialize( velName, pressName);

    UInt nnode = this->_mesh.pointList.size();

    std::string filenamep = pressName;
    std::string ext       = ".mesh";

    filenamep.insert(filenamep.end(), ext.begin(), ext.end());

//     std::cout << "    F - Reading INRIA fluid mesh file   (" << filenamep << ")"
//               << ":" << std::endl;

    std::ifstream file(filenamep.c_str(), std::ios::in);

    if (!file)
    {
        std::cout << "    F- initialize: Reading file " << filenamep
                  << " impossible" << std::endl;
        exit(1);
    }

    std::string sdummy;
    UInt        ns;

    file >> sdummy >> sdummy;
    file >> sdummy >> sdummy;
    file >> sdummy;
    file >> ns;

    if (ns != nnode)
    {
        std::cout << "     F- initialize: non-matching grids " << ns << " != " << nnode << std::endl;
        exit(1);
    }

    PhysVectUnknown<Vector> disp(nnode);

    UInt   idummy;
    UInt   inode;

    for (inode = 0; inode < nnode; ++inode)
    {
        file >> disp[inode + 0*nnode]
             >> disp[inode + 1*nnode]
             >> disp[inode + 2*nnode]
             >> idummy;

        disp[inode + 0*nnode] = (disp[inode + 0*nnode] - this->_mesh.pointList(inode + 1).x())/_factor;
        disp[inode + 1*nnode] = (disp[inode + 1*nnode] - this->_mesh.pointList(inode + 1).y())/_factor;
        disp[inode + 2*nnode] = (disp[inode + 2*nnode] - this->_mesh.pointList(inode + 1).z())/_factor;
    }

    M_harmonicExtension.setDisplacement(disp);

    file.close();

    //_mesh.moveMesh(disp);

    readUnknown(velwName, _w);

    //now we must compute the "old" displacement

    _dispOld = disp - _w*this->_dt;

//    this->_mesh.moveMesh(_disp);

}


template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::readUnknown( const std::string       &name,
                                            PhysVectUnknown<Vector> &unknown)
//                                            double                  factor)
{
    std::string sdummy;
    std::string ext;
    int nsx, nsy, nsz;
    int ndim;

    int nDof = this->_mesh.pointList.size();

    std::string filenamex = name;
    ext = "_x.bb";
    filenamex.insert(filenamex.end(), ext.begin(), ext.end());

//     std::cout << "Reading INRIA solid file   (" << filenamex << ")"
//               << ":" << std::endl;

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

    if (nsx != nDof)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamex << std::endl;
        std::cout << nDof << " != " << nsx << std::endl;
        exit(1);
    }

    filex >> sdummy;

    for (int ix = 0; ix < nsx; ++ix)
        {
            filex >> unknown[ix + 0*nDof];
        }

    filex.close();

    std::string filenamey = name;
    ext = "_y.bb";
    filenamey.insert(filenamey.end(), ext.begin(), ext.end());

//     std::cout << "Reading INRIA solid file   (" << filenamey << ")"
//               << ":" << std::endl;

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

    if (nsy != nDof)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamey << std::endl;
        std::cout << nDof << " != " << nsy << std::endl;
        exit(1);
    }

    filey >> sdummy;

    for (int iy = 0; iy < nsy; ++iy)
    {
        filey >> unknown[iy + 1*nDof];
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

    if (nsz != nDof)
    {
        std::cout << "restart: nonmatching dimension in file " << filenamez << std::endl;
        std::cout << nDof << " != " << nsz << std::endl;
        exit(1);
    }

    filez >> sdummy;

    for (int iz = 0; iz < nsz; ++iz)
    {
        filez >> unknown[iz + 2*nDof];
    }

    filez.close();
}



// Postprocessing pressure
template <typename Mesh>
void
NavierStokesAleHandler<Mesh>::postProcess()
{

    std::ostringstream index;
    std::string name;

    ++this->_count;
    if ( fmod( float( this->_count ), float( this->_verbose ) ) == 0.0 )
    {
        std::cout << "  F-  Post-processing \n";
        index << ( this->_count / this->_verbose );

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


        wr_medit_ascii_scalar( "press." + name + ".bb", this->_p.giveVec(), this->_p.size() );

//         wr_medit_ascii_scalar( "disp_x." + name + ".bb", _disp.giveVec(), this->_mesh.numVertices() );
//         wr_medit_ascii_scalar( "disp_y." + name + ".bb", _disp.giveVec() + this->_dim_u, this->_mesh.numVertices() );
//         wr_medit_ascii_scalar( "disp_z." + name + ".bb", _disp.giveVec() + 2 * this->_dim_u, this->_mesh.numVertices() );

        wr_medit_ascii_scalar( "vel_x." + name + ".bb", this->_u.giveVec(), this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", this->_u.giveVec() + this->_dim_u, this->_mesh.numVertices() );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", this->_u.giveVec() + 2 * this->_dim_u, this->_mesh.numVertices() );

        wr_medit_ascii_scalar("velw_x."+name+".bb",_w.giveVec(),this->_mesh.numVertices());
        wr_medit_ascii_scalar("velw_y."+name+".bb",_w.giveVec() + this->_mesh.numVertices(), this->_mesh.numVertices());
        wr_medit_ascii_scalar("velw_z."+name+".bb",_w.giveVec() + 2*this->_mesh.numVertices(), this->_mesh.numVertices());


        //wr_medit_ascii("press."+name+".mesh", this->_mesh);
        wr_medit_ascii( "press." + name + ".mesh", this->_mesh, M_harmonicExtension.getDisplacement(), _factor );
        // wr_medit_ascii_vector("veloc."+name+".bb",_u.giveVec(), this->_mesh.numVertices(),_dim_u);
        system( ( "ln -s press." + name + ".mesh vel_x." + name + ".mesh" ).data() );
        system( ( "ln -s press." + name + ".mesh vel_y." + name + ".mesh" ).data() );
        system( ( "ln -s press." + name + ".mesh vel_z." + name + ".mesh" ).data() );

        system(("ln -s press."+name+".mesh velw_x."+name+".mesh").data());
        system(("ln -s press."+name+".mesh velw_y."+name+".mesh").data());
        system(("ln -s press."+name+".mesh velw_z."+name+".mesh").data());

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


// This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > this->_mesh.getRefFE().nbNodes)
template <typename Mesh>
void NavierStokesAleHandler<Mesh>::
_interpolate( const UInt nbComp, const Vector& w, Vector& wInterp )
{

    typedef typename Mesh::VolumeShape GeoShape; // Element shape

    UInt nDofpV = this->_refFE_u.nbDofPerVertex; // number of Dof per vertex
    UInt nDofpE = this->_refFE_u.nbDofPerEdge;   // number of Dof per edge
    UInt nDofpF = this->_refFE_u.nbDofPerFace;   // number of Dof per face
    UInt nDofpEl = this->_refFE_u.nbDofPerVolume; // number of Dof per Volume

    UInt nElemV = GeoShape::numVertices; // Number of element's vertices
    UInt nElemE = GeoShape::numEdges;    // Number of element's edges
    UInt nElemF = GeoShape::numFaces;    // Number of element's faces

    UInt nDofElem = this->_mesh.getRefFE().nbDof; // Number of local dof per element of the this->_mesh (_mesh.getRefFE().nbDof)

    UInt nDofElemV = nElemV * nDofpV; // number of vertex's Dof on a Element
    UInt nDofElemE = nElemE * nDofpE; // number of edge's Dof on a Element
    UInt nDofElemF = nElemF * nDofpF; // number of face's Dof on a Element


    Real x, y, z;

    KN<Real> wLoc( nDofElem * nbComp );

    ID lDof;

    // Loop on elements of the mesh
    for ( ID iElem = 1; iElem <= this->_mesh.numVolumes(); ++iElem )
    {

        // Updating the local mesh velocity in this mesh elment
        for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            for ( ID idof = 0; idof < nDofElem; ++idof )
                wLoc( icmp * nDofElem + idof ) =
                    w( icmp * M_harmonicExtension.dofMesh().numTotalDof()
                       + M_harmonicExtension.dofMesh().localToGlobal( iElem, idof + 1 ) - 1 );

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
                    x = this->_refFE_u.xi( lDof - 1 );
                    y = this->_refFE_u.eta( lDof - 1 );
                    z = this->_refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the element
                            __sum += wLoc( icmp * nDofElem + idof ) * this->_mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolated mesh velocity
                        wInterp( icmp * this->_dof_u.numTotalDof() + this->_dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
                    x = this->_refFE_u.xi( lDof - 1 );
                    y = this->_refFE_u.eta( lDof - 1 );
                    z = this->_refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )   // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * this->_mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        wInterp( icmp * this->_dof_u.numTotalDof() + this->_dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
                    x = this->_refFE_u.xi( lDof - 1 );
                    y = this->_refFE_u.eta( lDof - 1 );
                    z = this->_refFE_u.zeta( lDof - 1 );

                    // Loop on data vector components
                    for ( UInt icmp = 0; icmp < nbComp; ++icmp )
                    {

                        // Interpolating data at the nodal point
                        double __sum = 0;
                        for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the adjacent element
                            __sum += wLoc( icmp * nDofElem + idof ) * this->_mesh.getRefFE().phi( idof, x, y, z );

                        // Updating interpolating vector
                        wInterp( icmp * this->_dof_u.numTotalDof() + this->_dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
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
            x = this->_refFE_u.xi( lDof - 1 );
            y = this->_refFE_u.eta( lDof - 1 );
            z = this->_refFE_u.zeta( lDof - 1 );

            // Loop on data vector components
            for ( UInt icmp = 0; icmp < nbComp; ++icmp )
            {

                // Interpolating data at the nodal point
                double __sum = 0;
                for ( ID idof = 0; idof < nDofElem; ++idof )  // Loop on local Dof on the adjacent element
                    __sum += wLoc( icmp * nDofElem + idof ) * this->_mesh.getRefFE().phi( idof, x, y, z );

                // Updating interpolating vector
                wInterp( icmp * this->_dof_u.numTotalDof() + this->_dof_u.localToGlobal( iElem, lDof ) - 1 ) = __sum;
            }
        }
    }
}
}
#endif
