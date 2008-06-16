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
  \file convDiffReactSolverPC.h
  \author M. Prosi
  \date 03/2004
  \version 1.0

  \brief This file contains a solver class for the Convection-Diffusion-Reaction equation
*/

#ifndef _CONVDIFFREACTSOLVERPC_H_
#define _CONVDIFFREACTSOLVERPC_H_

#include <life/lifesolver/convDiffReactHandler.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifearray/elemVec.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/assembGeneric.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifealg/algebraic_facto.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifecore/chrono.hpp>
//#include <life/lifealg/dataAztec.hpp>
#include <life/lifefem/bdf.hpp>
#include <life/lifefilters/openDX_wrtrs.hpp>
#include <string>

namespace LifeV
{
/*!
  \class convDiffReactSolverPC

   This class implements a solver for the Convection-Diffusion-Reaction equation

*/
template <typename Mesh>
class ConvDiffReactSolverPC:
            public ConvDiffReactHandler<Mesh>
{

public:

    typedef typename ConvDiffReactHandler<Mesh>::Function Function;
    typedef typename ConvDiffReactHandler<Mesh>::source_type source_type;

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE_c reference FE for the concentration
      \param Qr_c volumic quadrature rule for the concentration
      \param bdQr_c surface quadrature rule for the concentration
      \param BCh_c boundary conditions for the concentration
    */
    ConvDiffReactSolverPC( const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
                           const QuadRule& bdQr_c, BCHandler& BCh_c );

    //! Update the right  hand side  for time advancing
    /*!
      \param source volumic source
      \param time present time
    */
    void timeAdvance( source_type const& source, const Real& time );

    //! Update convective term, bc treatment and solve the linearized ns system
    void iterate( const Real& time );

    //! Projection of the velocity on grid of concentration discretization
    template <typename RegionMesh3D>
    void getvel( RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BCHandler& BCh_u, const Real& time );

    //! Calculate the local coordinates of concentration gridpoints in the
    //! velocity grid (is needed for the Projection)
    template <typename RegionMesh3D>
    void getcoord( RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BCHandler& BCh_u );

    //! Calculate the volume of a tetrahedra given by its corner nodes
    Real calcvol( SimpleVect<GeoElement0D<DefMarkerCommon> > points  );

    //! tests if point (xp, yp, zp) is in the tetrahedra (x[4], y[4], z[4]) and returns
    //! interpolation coefficents (1-b1-b2-b3, b1, b2, b3)
      void test( SimpleVect<GeoElement0D<DefMarkerCommon> > points_tet, SimpleVect<GeoElement0D<DefMarkerCommon> > point, Real & b1, Real & b2, Real & b3 );
    template <typename RegionMesh3D>
    void checkallvolumes( RegionMesh3D & umesh, SimpleVect<GeoElement0D<DefMarkerCommon> > point, Real & b1r, Real & b2r, Real & b3r );

protected:

    // inherited from parent class
    typedef typename ConvDiffReactHandler<Mesh>::intpolcoord intpolcoord;

private:

    //! Pattern of M
    MSRPatt _pattM;

    //! Matrix C:  1/dt*Cmass + D*Cstiff operator + r*Cmass
    MSRMatr<double> _DR;

    //! Matrix C:  1/dt*Cmass + D*Cstiff operator + Convective_transport term + r*Cmass
    MSRMatr<double> _CDR;

    //! Matrix C_u: Cmass
    MSRMatr<double> _M_c;

    //! Elementary matrices and vectors
    ElemMat _elmatC; //Concentration stiffnes
    ElemMat _elmatM_c; //Concentration mass
    ElemMat _elmatMR_c; // Concentration mass and reaction
    ElemVec _elvec; // Elementary right hand side
    ElemVec _elvec_u; // Elementary velocity for convection term

    //! Right  hand  side for the concentration
    ScalUnknown<Vector> _f_c;
    ScalUnknown<Vector> _f_c_noBC;

    DataAztec _dataAztec_o;
};

//
//                                         IMPLEMENTATION
//
template <typename Mesh>
ConvDiffReactSolverPC<Mesh>::
ConvDiffReactSolverPC( const GetPot& data_file, const RefFE& refFE_c, const QuadRule& Qr_c,
                       const QuadRule& bdQr_c, BCHandler& BCh_c ) :
        ConvDiffReactHandler<Mesh>( data_file, refFE_c, Qr_c, bdQr_c, BCh_c ),
        _pattM( this->_dof_c ),
        _DR( _pattM ),
        _CDR( _pattM ),
        _M_c( _pattM ),
        _elmatC( this->_fe_c.nbNode, 1, 1 ),
        _elmatM_c( this->_fe_c.nbNode, 1, 1 ),
        _elmatMR_c( this->_fe_c.nbNode, 1, 1 ),
        _elvec( this->_fe_c.nbNode, 1 ),
        _elvec_u( this->_fe_c.nbNode, nDimensions ),
        _f_c( this->_dim_c ),
        _dataAztec_o( data_file, "masstransport/aztec_o" )
{

    std::cout << std::endl;
    std::cout << "O-  Concentration unknowns: " << this->_dim_c << std::endl;
    std::cout << "O-  Computing mass and stiffness matrices... ";

    Chrono chrono;
    chrono.start();

    // Matrices initialization
    _DR.zeros();
    _CDR.zeros();
    _M_c.zeros();

    //inverse of dt:
    Real dti = 1. / this->timestep();

    // *******************************************************
    // Coefficient of the mass term at time t^{n+1}
    Real first_coeff = this->_bdf.coeff_der( 0 );
    std::cout << std::endl;
    std::cout << "Bdf CDR first coeff " << first_coeff << std::endl;

    this->_bdf.showMe();

    // Elementary computation and matrix assembling

    for ( UInt i = 1; i <= this->mesh().numVolumes(); i++ )
    {          // Loop on elements

        this->_fe_c.updateFirstDerivQuadPt( this->mesh().volumeList( i ) );

        _elmatC.zero();
        _elmatM_c.zero();
        _elmatMR_c.zero();

        stiff( this->_diffusivity, _elmatC, this->_fe_c );

        if ( this->_stationary )
        {
            mass( ( -this->_react ), _elmatMR_c, this->_fe_c );
        }
        else
        {
            mass( ( first_coeff * dti - this->_react ), _elmatMR_c, this->_fe_c );
            mass( dti, _elmatM_c, this->_fe_c );
        }

        _elmatC.mat() += _elmatMR_c.mat();

        // stiffness + mass + reaction term
        assemb_mat( _DR, _elmatC, this->_fe_c, this->_dof_c );

        // mass
        assemb_mat( _M_c, _elmatM_c, this->_fe_c, this->_dof_c );

    }

    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;

}

template <typename Mesh>
void ConvDiffReactSolverPC<Mesh>::
timeAdvance( source_type const& source, const Real& time )
{

    std::cout << "  o-  Updating mass term on right hand side (concentration)... ";

    Chrono chrono;
    chrono.start();

    // Right hand side for the velocity at time
    _f_c = ZeroVector( _f_c.size() );

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= this->mesh().numVolumes(); ++i )
    {
        _elvec.zero();
        this->_fe_c.update( this->mesh().volumeList( i ) );

        compute_vec( source, _elvec, this->_fe_c, time, 0 ); // compute local vector
        assemb_vec( _f_c, _elvec, this->_fe_c, this->_dof_c, 0 ); // assemble local vector into global one
    }

    // *******************************************************
    _f_c += _M_c * this->_bdf.time_der(); //_M_u is the mass matrix divided by the time step
    _f_c_noBC = _f_c;
    chrono.stop();
    std::cout << "done in " << chrono.diff() << " s." << std::endl;
}


template <typename Mesh>
void ConvDiffReactSolverPC<Mesh>::
iterate( const Real& time )
{

    UInt nc_u = this->_u_c.nbcomp();

    Chrono chrono;

    // CDR = DR + convective term (C)
    chrono.start();
    _CDR = _DR;
    _f_c = _f_c_noBC;
    chrono.stop();


    std::cout << "  o-  Diffusion-Reaction matrix was copied in " << chrono.diff() << "s." << std::endl;
    std::cout << "  o-  Updating convective transport... ";

    chrono.start();

    // loop on volumes
    for ( UInt i = 1; i <= this->mesh().numVolumes(); ++i )
    {

        this->_fe_c.updateFirstDeriv( this->mesh().volumeList( i ) ); // as updateFirstDer

        _elmatC.zero();

        UInt eleID = this->_fe_c.currentId();

        // ********** copy global velocity vector to local velocity vector *************
        // ********** assuming velocity is given on concentration mesh *****************

        for ( UInt k = 0 ; k < ( UInt ) this->_fe_c.nbNode ; k++ )
        {
            UInt iloc = this->_fe_c.patternFirst( k );
            for ( UInt ic = 0; ic < nc_u; ++ic )
            {
                UInt ig = this->_dof_c.localToGlobal( eleID, iloc + 1 ) - 1 + ic * this->_dim_c;
                _elvec_u[ iloc + ic * this->_fe_c.nbNode ] = this->_u_c( ig );
            }
        }

        grad( 0, _elvec_u, _elmatC, this->_fe_c, this->_fe_c, this->_fe_c );
        grad( 1, _elvec_u, _elmatC, this->_fe_c, this->_fe_c, this->_fe_c );
        grad( 2, _elvec_u, _elmatC, this->_fe_c, this->_fe_c, this->_fe_c );


        // *************************************** Upwind ******************************

        Real VLoc_infty = 0.;
        Real VLoc_mean = 0.;
        Real VLoc_c = 0.;
        for ( UInt ih_c = 0 ; ih_c < ( UInt ) this->_fe_c.nbNode ; ih_c++ )
        {
            UInt iloc = this->_fe_c.patternFirst( ih_c );
            for ( UInt ic = 0; ic < nc_u;++ic )
            {
                UInt ig = this->_dof_c.localToGlobal( eleID, iloc + 1 ) - 1 + ic * this->_dim_c;
                _elvec_u[ iloc + ic * this->_fe_c.nbNode ] = this->_u_c( ig );
                VLoc_c += this->_u_c( ig ) * this->_u_c( ig );
            }
            VLoc_c = sqrt( VLoc_c );
            VLoc_mean += VLoc_c;
            if ( VLoc_c > VLoc_infty )
                VLoc_infty = VLoc_c;
        }
        VLoc_mean = VLoc_mean / this->_fe_c.nbNode;

        Real coef_stab, Pe_loc;
        //      coef_stab=this->_fe_c.diameter()*VLoc_infty; // Alessandro - method

        Pe_loc = VLoc_infty * this->_fe_c.diameter() / ( 2.0 * this->_diffusivity );

        //      coef_stab=(1.0/tanh(Pe_loc))-(1.0/Pe_loc); // classical approach

        if ( Pe_loc < -3.0 )
            coef_stab = -1.0;
        else
        {
            if ( Pe_loc > 3.0 )
                coef_stab = 1.0;
            else
                coef_stab = Pe_loc / 3.0;
        }

        // ******************************* STREAMLINEUPWIND ****************************
        stiff_sd( coef_stab / ( VLoc_mean * VLoc_mean ), _elvec_u, _elmatC, this->_fe_c, this->_fe_c );

        // ************************* Assembling ****************************************

        assemb_mat( _CDR, _elmatC, this->_fe_c, this->_dof_c );

    }
    chrono.stop();
    std::cout << "done in " << chrono.diff() << "s." << std::endl;

    // for BC treatment (done at each time-step)
    Real tgv = 1.e02;

    std::cout << "  o-  Applying boundary conditions... ";
    chrono.start();
    // BC manage for the concentration
    if ( !this->_BCh_c.bdUpdateDone() )
        this->_BCh_c.bdUpdate( this->mesh(), this->_feBd_c, this->_dof_c );
    bcManage( _CDR, _f_c, this->mesh(), this->_dof_c, this->_BCh_c, this->_feBd_c, tgv, time );
    chrono.stop();

    std::cout << "done in " << chrono.diff() << "s." << std::endl;


    int proc_config_o[ AZ_PROC_SIZE ]; // Processor information:
    //  proc_config[AZ_node] = node name
    //  proc_config[AZ_N_procs] = # of nodes
    int options_o[ AZ_OPTIONS_SIZE ]; // Array used to select solver options.
    double params_o[ AZ_PARAMS_SIZE ];   // User selected solver paramters.
    int *data_org_o;                // Array to specify data layout
    double status_o[ AZ_STATUS_SIZE ];   // Information returned from AZ_solve()
    // indicating success or failure.
    // other declarations for AZTEC
    int *update_o,                   // vector elements updated on this node.
    *external_o;                // vector elements needed by this node.
    int *update_index_o;            // ordering of update[] and external[]
    int *extern_index_o;            // locally on this processor.
    //  int    *bindx;                 // Sparse matrix to be solved is stored
    //  double *val;                   // in these MSR arrays.
    int N_update_o;                 // # of unknowns updated on this node

    AZ_set_proc_config( proc_config_o, AZ_NOT_MPI );

    AZ_read_update( &N_update_o, &update_o, proc_config_o, this->_dim_c, 1, AZ_linear );
    AZ_defaults( options_o, params_o );
    _dataAztec_o.aztecOptionsFromDataFile( options_o, params_o );
    AZ_transform( proc_config_o, &external_o,
                  ( int * ) _pattM.giveRaw_bindx(), _CDR.giveRaw_value(),
                  update_o, &update_index_o,
                  &extern_index_o, &data_org_o, N_update_o, NULL, NULL, NULL, NULL,
                  AZ_MSR_MATRIX );

    chrono.start();
    //    init_options_c(options_o,params_o);

    AZ_solve( this->_c.giveVec(), _f_c.giveVec(), options_o, params_o, NULL,
              ( int * ) _pattM.giveRaw_bindx(), NULL, NULL, NULL,
              _CDR.giveRaw_value(), data_org_o, status_o, proc_config_o );

    chrono.stop();
    std::cout << "*** Solution (Concentration) computed in " << chrono.diff() << "s." << std::endl;
    this->_bdf.shift_right( this->_c );

}



template <typename Mesh>
template <typename RegionMesh3D>
void ConvDiffReactSolverPC<Mesh>::
getvel( RegionMesh3D & umesh, PhysVectUnknown<Vector> & u, BCHandler& BCh_u, const Real& time )
{


    for ( UInt i = 0; i < this->mesh().numVertices(); i++ )
    {

        if ( this->_u_to_c[ i ].ele == 0 )
        {
       // Dirichlet boundary for the velocity -> get velocity for boundary function
            Real xp = this->mesh().point( i + 1 ).x();
            Real yp = this->mesh().point( i + 1 ).y();
            Real zp = this->mesh().point( i + 1 ).z();
            for ( ID jj = 0; jj < 3; ++jj )
                this->_u_c( i + jj * this->_u_c.size() / 3 ) = BCh_u[ ( int ) this->_u_to_c[ i ].b[ 0 ] ] ( time, xp, yp, zp, BCh_u[ ( int ) this->_u_to_c[ i ].b[ 0 ] ].component( jj + 1 ) );

        }
        else
        {
       // Velocity interpolation
            this->_u_c( i ) = this->_u_to_c[ i ].b[ 0 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 1 ).id() - 1 ) +
                        this->_u_to_c[ i ].b[ 1 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 2 ).id() - 1 ) +
                        this->_u_to_c[ i ].b[ 2 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 3 ).id() - 1 ) +
                        this->_u_to_c[ i ].b[ 3 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 4 ).id() - 1 );

            this->_u_c( i + this->_dim_c ) = this->_u_to_c[ i ].b[ 0 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 1 ).id() - 1 + u.size() / 3 ) +
                                 this->_u_to_c[ i ].b[ 1 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 2 ).id() - 1 + u.size() / 3 ) +
                                 this->_u_to_c[ i ].b[ 2 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 3 ).id() - 1 + u.size() / 3 ) +
                                 this->_u_to_c[ i ].b[ 3 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 4 ).id() - 1 + u.size() / 3 );

            this->_u_c( i + 2 * this->_dim_c ) = this->_u_to_c[ i ].b[ 0 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 1 ).id() - 1 + 2 * u.size() / 3 ) +
                                     this->_u_to_c[ i ].b[ 1 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 2 ).id() - 1 + 2 * u.size() / 3 ) +
                                     this->_u_to_c[ i ].b[ 2 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 3 ).id() - 1 + 2 * u.size() / 3 ) +
                                     this->_u_to_c[ i ].b[ 3 ] * u( umesh.volume( this->_u_to_c[ i ].ele ).point( 4 ).id() - 1 + 2 * u.size() / 3 );

        }


    }

}

template <typename Mesh>
template <typename RegionMesh3D>
void ConvDiffReactSolverPC<Mesh>::
getcoord( RegionMesh3D & umesh, PhysVectUnknown<Vector> & /*u*/, BCHandler& BCh_u )
{

    Real b1, b2, b3;
    intpolcoord localcoord;

    Chrono chrono;
    chrono.start();

    this->_u_c = ScalarVector( this->_u_c.size(), -100 );

    UInt vid, adjacent, minadjacent;
    LinearTetra ele;

    SimpleVect<GeoElement0D<DefMarkerCommon> > tetrapoints, tetrapoints_found;
    tetrapoints.resize( 4 );
    tetrapoints_found.resize( 4 );

    SimpleVect<GeoElement3D<LinearTetra> >::iterator iv = umesh.volumeList.begin();

   for ( UInt i = 0; i < this->mesh().numVertices(); i++ )
    {
    tetrapoints( 1 ) = this->mesh().point( i + 1 );

        if ( this->mesh().point( i + 1 ).boundary() )
        {
            UInt l;

            for ( UInt k = 0; k < BCh_u.size(); k++ )
            {
                if ( BCh_u[ k ].flag() == this->mesh().point( i + 1 ).marker() )
                {
                    l = k;
                }
            }

            switch ( BCh_u[ l ].type() )
            {
            case Essential:
                localcoord.b[ 0 ] = ( Real ) l;
                localcoord.b[ 1 ] = ( Real ) l;
                localcoord.b[ 2 ] = ( Real ) l;
                localcoord.b[ 3 ] = ( Real ) l;
                localcoord.ele = 0;
                this->_u_to_c.push_back( localcoord );
                break;
            default:
            for ( Int found = 0;found == 0; )
            {
                vid = iv -> id();
                Real volume, minvolume = 100.0, minminvolume = 100.0;
                UInt jk = UInt( -1 );
                for ( UInt jj = 1; jj <= umesh.numLocalFaces();jj++ )
                {

            tetrapoints( 2 ) = umesh.point( ( iv->point( ele.fToP( jj, 1 ) ) ).id() );
            tetrapoints( 3 ) = umesh.point( ( iv->point( ele.fToP( jj, 2 ) ) ).id() );
            tetrapoints( 4 ) = umesh.point( ( iv->point( ele.fToP( jj, 3 ) ) ).id() );
            volume = calcvol( tetrapoints );

                        if ( umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_first() == vid )
                            {
                                adjacent = umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_second();
                            }
                        else
                            {
                                adjacent = umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_first();
                            }

                       if (volume < minminvolume)
                        {
               minminvolume = volume;
            }
                       if (volume < minvolume)
                        {
                            if (adjacent != 0)
                                {
                                    minvolume = volume;
                                    minadjacent=adjacent;
                                    jk = jj;
                                }
                        }
        }

                if ( minvolume < -1.0e-15 )
                    {
                        iv = iv + minadjacent - vid;
                    }

                else
                    {
                        if(minminvolume < -1.0e-15)
                            {
                                found = 1;
                                checkallvolumes( umesh, tetrapoints, b1, b2, b3 );
                                localcoord.b[ 0 ] = 1.0 - b1 - b2 - b3;
                                localcoord.b[ 1 ] = b1;
                                localcoord.b[ 2 ] = b2;
                                localcoord.b[ 3 ] = b3;
                                localcoord.ele = vid;
                                this->_u_to_c.push_back( localcoord );
                            }
                        else
                            {

                                found = 1;
                tetrapoints_found( 1 ) = umesh.point( ( iv->point( 1 ) ).id() );
                tetrapoints_found( 2 ) = umesh.point( ( iv->point( 2 ) ).id() );
                tetrapoints_found( 3 ) = umesh.point( ( iv->point( 3 ) ).id() );
                tetrapoints_found( 4 ) = umesh.point( ( iv->point( 4 ) ).id() );

                test( tetrapoints_found, tetrapoints, b1, b2, b3 );
                                localcoord.b[ 0 ] = 1.0 - b1 - b2 - b3;
                                localcoord.b[ 1 ] = b1;
                                localcoord.b[ 2 ] = b2;
                                localcoord.b[ 3 ] = b3;
                                localcoord.ele = vid;
                                this->_u_to_c.push_back( localcoord );
                            }
                    }
        }

                   break;
            }
        }
        else
        {
            for ( Int found = 0;found == 0; )
            {
                vid = iv -> id();
                Real volume, minvolume = 100.0, minminvolume = 100.0;
                UInt jk = UInt( -1 );
                for ( UInt jj = 1; jj <= umesh.numLocalFaces();jj++ )
                {

            tetrapoints( 2 ) = umesh.point( ( iv->point( ele.fToP( jj, 1 ) ) ).id() );
            tetrapoints( 3 ) = umesh.point( ( iv->point( ele.fToP( jj, 2 ) ) ).id() );
            tetrapoints( 4 ) = umesh.point( ( iv->point( ele.fToP( jj, 3 ) ) ).id() );
            volume = calcvol( tetrapoints );

                        if ( umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_first() == vid )
                            {
                                adjacent = umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_second();
                            }
                        else
                            {
                                adjacent = umesh.faceList( umesh.localFaceId( vid, jj ) ).ad_first();
                            }

                       if (volume < minminvolume)
                        {
               minminvolume = volume;
            }
                       if (volume < minvolume)
                        {
                            if (adjacent != 0)
                                {
                                    minvolume = volume;
                                    minadjacent=adjacent;
                                    jk = jj;
                                }
                        }


                }

                if ( minvolume < -1.0e-15 )
                    {
                        iv = iv + minadjacent - vid;
                    }

                else
                    {
                        if(minminvolume < -1.0e-15)
                            {
                                found = 1;
                                checkallvolumes( umesh, tetrapoints, b1, b2, b3 );
                                localcoord.b[ 0 ] = 1.0 - b1 - b2 - b3;
                                localcoord.b[ 1 ] = b1;
                                localcoord.b[ 2 ] = b2;
                                localcoord.b[ 3 ] = b3;
                                localcoord.ele = vid;
                                this->_u_to_c.push_back( localcoord );
                            }
                        else
                            {
                                found = 1;
                tetrapoints_found( 1 ) = umesh.point( ( iv->point( 1 ) ).id() );
                tetrapoints_found( 2 ) = umesh.point( ( iv->point( 2 ) ).id() );
                tetrapoints_found( 3 ) = umesh.point( ( iv->point( 3 ) ).id() );
                tetrapoints_found( 4 ) = umesh.point( ( iv->point( 4 ) ).id() );

                test( tetrapoints_found, tetrapoints, b1, b2, b3 );
                                localcoord.b[ 0 ] = 1.0 - b1 - b2 - b3;
                                localcoord.b[ 1 ] = b1;
                                localcoord.b[ 2 ] = b2;
                                localcoord.b[ 3 ] = b3;
                                localcoord.ele = vid;
                                this->_u_to_c.push_back( localcoord );
                            }
                    }
            }

        }
    }

    chrono.stop();
    std::cout << " Calculation of the projection coordinates " << chrono.diff() << "s." << std::endl;

}

template <typename Mesh>
Real ConvDiffReactSolverPC<Mesh>::
calcvol( SimpleVect<GeoElement0D<DefMarkerCommon> > points  )
{
    Real volume = 0.0;
    volume += ( ( points(2).x() - points(3).x() ) * ( points(2).y() - points(4).y() ) * ( points(2).z() - points(1).z() ) );
    volume += ( ( points(2).y() - points(3).y() ) * ( points(2).z() - points(4).z() ) * ( points(2).x() - points(1).x() ) );
    volume += ( ( points(2).z() - points(3).z() ) * ( points(2).x() - points(4).x() ) * ( points(2).y() - points(1).y() ) );
    volume -= ( ( points(2).x() - points(1).x() ) * ( points(2).y() - points(4).y() ) * ( points(2).z() - points(3).z() ) );
    volume -= ( ( points(2).y() - points(1).y() ) * ( points(2).z() - points(4).z() ) * ( points(2).x() - points(3).x() ) );
    volume -= ( ( points(2).z() - points(1).z() ) * ( points(2).x() - points(4).x() ) * ( points(2).y() - points(3).y() ) );
    return volume;
}



template <typename Mesh>
void ConvDiffReactSolverPC<Mesh>::
test( SimpleVect<GeoElement0D<DefMarkerCommon> > points_tet, SimpleVect<GeoElement0D<DefMarkerCommon> > point, Real & b1, Real & b2, Real & b3 )
{

    Real a11, a12, a13, a21, a22, a23, a31, a32, a33, zw;

    a11 = points_tet(2).x() - points_tet(1).x();
    a12 = points_tet(3).x() - points_tet(1).x();
    a13 = points_tet(4).x() - points_tet(1).x();
    a21 = points_tet(2).y() - points_tet(1).y();
    a22 = points_tet(3).y() - points_tet(1).y();
    a23 = points_tet(4).y() - points_tet(1).y();
    a31 = points_tet(2).z() - points_tet(1).z();
    a32 = points_tet(3).z() - points_tet(1).z();
    a33 = points_tet(4).z() - points_tet(1).z();

    b1 = point(1).x() - points_tet(1).x();
    b2 = point(1).y() - points_tet(1).y();
    b3 = point(1).z() - points_tet(1).z();

    //  Solve the equation system (without loop - faster)

    if ( fabs( a11 ) < std::max<Real>( fabs( a21 ), fabs( a31 ) ) )
    {
        if ( fabs( a21 ) > fabs( a31 ) )
        {
            zw = a11;
            a11 = a21;
            a21 = zw;
            zw = a12;
            a12 = a22;
            a22 = zw;
            zw = a13;
            a13 = a23;
            a23 = zw;
            zw = b1;
            b1 = b2;
            b2 = zw;
        }
        else
        {
            zw = a11;
            a11 = a31;
            a31 = zw;
            zw = a12;
            a12 = a32;
            a32 = zw;
            zw = a13;
            a13 = a33;
            a33 = zw;
            zw = b1;
            b1 = b3;
            b3 = zw;
        }
    }

    zw = a21 / a11;
    a22 = a22 - zw * a12;
    a23 = a23 - zw * a13;
    b2 = b2 - zw * b1;
    zw = a31 / a11;
    a32 = a32 - zw * a12;
    a33 = a33 - zw * a13;
    b3 = b3 - zw * b1;

    if ( fabs( a32 ) > fabs( a22 ) )
    {
        zw = a22;
        a22 = a32;
        a32 = zw;
        zw = a23;
        a23 = a33;
        a33 = zw;
        zw = b2;
        b2 = b3;
        b3 = zw;
    }

    zw = a32 / a22;
    a33 = a33 - zw * a23;

    b3 = b3 - zw * b2;
    b3 = b3 / a33;
    b2 = ( b2 - a23 * b3 ) / a22;
    b1 = ( b1 - a12 * b2 - a13 * b3 ) / a11;

 }

template <typename Mesh>
template <typename RegionMesh3D>
void ConvDiffReactSolverPC<Mesh>::
checkallvolumes( RegionMesh3D & umesh, SimpleVect<GeoElement0D<DefMarkerCommon> > point, Real & b1r, Real & b2r, Real & b3r )
{
    Real b1, b2, b3, bmin = -100.0;
    SimpleVect<GeoElement0D<DefMarkerCommon> > tetrapoints;
    tetrapoints.resize(4);

    b1r=100.0;
    b2r=100.0;
    b3r=100.0;

    for ( UInt i = 1; i <= umesh.numVolumes(); i++ )
    {
       tetrapoints( 1 ) = umesh.point( umesh.volume(i).point(1).id() );
       tetrapoints( 2 ) = umesh.point( umesh.volume(i).point(2).id() );
       tetrapoints( 3 ) = umesh.point( umesh.volume(i).point(3).id() );
       tetrapoints( 4 ) = umesh.point( umesh.volume(i).point(4).id() );

       test( tetrapoints, point, b1, b2, b3 );

       if( (b1 >= bmin) && (b2 >= bmin) && (b3 >= bmin) && (1.0-b1-b2-b3 >= bmin))
       {
      b1r = b1;
      b2r = b2;
      b3r = b3;
      bmin = std::min<Real>(b1,std::min<Real>(b2,std::min<Real>(b3,1.0-b1-b2-b3)));
       }

    }
}

}

#endif
