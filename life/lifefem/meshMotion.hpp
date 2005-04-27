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
  \file meshMotion.h
  \brief Classes to hold algorithms for the mesh motion, for instance, involved in a ALE formulation.
  \version 1.0
  \author M.A. Fernandez
  \date 11/2002

  This file contains classes which may be used to compute the extension inside the reference domain of a given
  displacement at a specified interface

*/

#ifndef __MESHMOTION_HH__
#define __MESHMOTION_HH__

#include <life/lifecore/life.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifearray/elemMat.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefem/bcManage.hpp>
#include <life/lifealg/SolverAztec.hpp>



namespace LifeV
{
/*!
  \class HarmonicExtension

  Base class which provides the harmonic extension of a given displacement on a specified part
  of the mesh boundary

  In order to deal with harmonic extensions, we have to provide a mesh (to  be moved), the parameters
  involved in the laplacian discretization: viscosity, quadrature rules and boundary conditions.
  This class contains a   PhysVectUnknown objet wich will hold the extension of the interface displacement.
  The constructor of the class built the global matrix of the discretized laplacian. The extension of the
  displacement is computed by calling the public method update. Finally, this extension can be recovered by
  calling method getDisplacement.

*/

class HarmonicExtension
{
public:

    //! Constructors for an harmonics extensions
    /*!
      \param mesh the mesh of the reference domain to be moved
      \param mu the "viscosity" in the laplacian operator
      \param Qr the quadrature rule for volumic elementary computations
      \param bdQr the quadrature rule for surface elementary computations
      \param mesh_BCh the list of boundary conditions involved in the harmonic extension

      \note The BCHandler objet (bch) holds the displacement imposed on moving boundary
      in the mesh trhough a BCVetor_Interface objet.

    */

    template <typename Mesh>

    HarmonicExtension(  const GetPot& data_file, Mesh& mesh,
                       const Real& diffusion,
                       const QuadRule& Qr,
                       const QuadRule& bdQr,
                       BCHandler&      mesh_BCh );

    template <typename Mesh>
    HarmonicExtension( Mesh&           mesh,
                       const Real&     diffusion,
                       const QuadRule& Qr,
                       const QuadRule& bdQr );

    //! This method updates the extension of the displacement, i.e. it solves the laplacian proglem

    template <typename Mesh>
    void updateExtension( Mesh& mesh, const Real& time = 0.0, const UInt recur = 0 );
    template <typename Mesh>


    void updateExtensionTransp( Mesh& mesh, const Real& time = 0.0 );
    //! This method gives a reference to the computed harmonic extension.
    Vector& getDisplacement();

    //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
    template <typename Mesh>
    void interpMeshVelocity( Mesh& mesh, const RefFE& refFE_u, const Dof& dof_u, Vector& wInterp );

    //! Returns a reference to the corresponding Dof object
    const Dof& dofMesh() const;

    //! checking if BC are set
    const bool setHarmonicExtensionBC() const {return M_setBC;}
    //! set the mesh BCs
    void setHarmonicExtensionBC(BCHandler &BCh_HarmonicExtension);
    //! returns the BCHandler
    BCHandler& BCh_HarmonicExtension() {return *M_BCh_HarmonicExtension;}

protected:

    //! Diffusion coefficient for the laplacian operator
    Real _diffusion;

    //! Quadrature rule for volumic elementary computations
    const QuadRule& _Qr;

    //! Quadrature rule for surface elementary computations
    const QuadRule& _bdQr;


    //! The Dof object associated with the displacement computations
    Dof _dof_mesh;

  //! The pattern of the FE discretized vector laplacian matrix
    MSRPatt _aPatt;

    //! The matrix holding the values
    MSRMatr<double> _a;

    //! Current element
    CurrentFE _fe;

    //! Current element  on the boundary
    CurrentBdFE _feBd;

    //! Elementary matrix : 3 blocks
    ElemMat _elmat;

    //! The actual extension of the displacement
    Vector _disp;

    //! Auxiliary vector holding the second right hand of the system
    PhysVectUnknown<Vector> _f;


    SolverAztec _linearSolver;


private:

    //! BC holding the imposed boundary displacement
    BCHandler    *M_BCh_HarmonicExtension;

    bool          M_setBC;


};


//! Constructor for an harmonics extensions
/*!
  \param mesh the mesh of the reference domain to be moved
  \param diffusion the "viscosity" in the laplacian operator
  \param Qr the quadrature rule for volumic elementary computations
  \param bdQr the quadrature rule for surface elementary computations
  \param bch the list of boundary conditions involved in the harmonic extension

  \note The BCHandler objet (bch) holds the displacement imposed on moving boundary
  in the mesh trhough a BCVetor_Interface objet.

*/
template <typename Mesh>
HarmonicExtension::
HarmonicExtension(  const GetPot& data_file,Mesh& mesh,
                   const Real& diffusion,
                   const QuadRule& Qr,
                   const QuadRule& bdQr,
                   BCHandler& mesh_BCh ) :
        _diffusion( diffusion ),
        _Qr( Qr ),
        _bdQr( bdQr ),
        _dof_mesh( mesh, mesh.getRefFE() ),
        _aPatt( _dof_mesh, nDimensions ),
        _a( _aPatt ),
        _fe( mesh.getRefFE(), mesh.getGeoMap(), _Qr ),
        _feBd( mesh.getRefFE().boundaryFE(), mesh.getGeoMap().boundaryMap(), _bdQr ),
        _elmat( _fe.nbNode, nDimensions, nDimensions ),
        _disp( nDimensions * _dof_mesh.numTotalDof() ),
        _f( _dof_mesh.numTotalDof() ),
        M_BCh_HarmonicExtension ( &mesh_BCh )
 {
   // Loop on elements
   for ( UInt i = 1; i <= mesh.numVolumes(); ++i )
     {
       // Updating derivatives
       _fe.updateFirstDerivQuadPt( mesh.volumeList( i ) );
       _elmat.zero();
       stiff( _diffusion, _elmat, _fe, 0, 0, 3 );
       // Assembling
       for ( UInt j = 0; j < 3; ++j )
	 {
	   assemb_mat( _a, _elmat, _fe, _dof_mesh, j, j );
	 }
     }
   
   // Initializations
   _disp = ZeroVector( _disp.size() );

  _linearSolver.setOptionsFromGetPot( data_file, "mesh_motion/aztec" );
  _linearSolver.setMatrix( _a );


 }


template <typename Mesh>
HarmonicExtension::
HarmonicExtension( Mesh& mesh,
                   const Real& diffusion,
                   const QuadRule& Qr,
                   const QuadRule& bdQr ) :
        _diffusion( diffusion ),
        _Qr       ( Qr ),
        _bdQr     ( bdQr ),
        _dof_mesh ( mesh, mesh.getRefFE() ),
	_aPatt( _dof_mesh, nDimensions ),
        _a        ( _aPatt ),
        _fe       ( mesh.getRefFE(), mesh.getGeoMap(), _Qr ),
        _feBd     ( mesh.getRefFE().boundaryFE(), mesh.getGeoMap().boundaryMap(), _bdQr ),
        _elmat    ( _fe.nbNode, nDimensions, nDimensions ),
        _disp     ( nDimensions * _dof_mesh.numTotalDof() ),
        _f        ( _dof_mesh.numTotalDof() ),
        M_BCh_HarmonicExtension ( 0 )
{

    // Loop on elements
    for ( UInt i = 1; i <= mesh.numVolumes(); ++i )
    {
        // Updating derivatives
        _fe.updateFirstDerivQuadPt( mesh.volumeList( i ) );
        _elmat.zero();
        stiff( _diffusion, _elmat, _fe, 0, 0, 3 );
        // Assembling
        for ( UInt j = 0; j < 3; ++j )
        {
            assemb_mat( _a, _elmat, _fe, _dof_mesh, j, j );
        }
    }

    // Initializations
    _disp = ZeroVector( _disp.size() );
}

// This method updates the extension of the displacement, i.e. it solves the laplacian problem

template <typename Mesh>
void HarmonicExtension::updateExtension( Mesh& mesh, const Real& time, const UInt recur )
{
    if ( !BCh_HarmonicExtension().bdUpdateDone() )
    {
        // BC boundary information update
        BCh_HarmonicExtension().bdUpdate( mesh, _feBd, _dof_mesh );

        // Boundary conditions treatment on the matrix
        bcManageMatrix( _a, mesh, _dof_mesh, BCh_HarmonicExtension(), _feBd, 1.0 );
    }

    // Initializations
    _f = ZeroVector( _f.size() );
    _disp = ZeroVector( _disp.size() );

    // Boundary conditions treatment
    bcManageVector( _f, mesh, _dof_mesh, BCh_HarmonicExtension(), _feBd, time, 1.0 );

    _linearSolver.setRecursionLevel( recur );
  
    _linearSolver.solve( _disp, _f, SolverAztec::SAME_PRECONDITIONER);

}



// This method updates the extension of the displacement, i.e. it solves the laplacian proglem
template <typename Mesh>
void HarmonicExtension::updateExtensionTransp( Mesh& mesh, const Real& time )
{

    // Boundary conditions treatment
    bcManageVector( _disp, mesh, _dof_mesh, BCh_HarmonicExtension(), _feBd, time );
}
}
#endif
