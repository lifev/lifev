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

#include "lifeV.hpp"
#include "dof.hpp"
#include "pattern.hpp"
#include "elemMat.hpp"
#include "elemOper.hpp"
#include "refFE.hpp"
#include "values.hpp"
#include "assemb.hpp"
#include "bc_manage.hpp"


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

class HarmonicExtension {
 public:

  //! Constructor for an harmonics extensions
  /*!
    \param mesh the mesh of the reference domain to be moved
    \param mu the "viscosity" in the laplacian operator
    \param Qr the quadrature rule for volumic elementary computations
    \param bdQr the quadrature rule for surface elementary computations
    \param mesh_BCh the list of boundary conditions involved in the harmonic extension

    \note The BC_Handler objet (bch) holds the displacement imposed on moving boundary
    in the mesh trhough a BCVetor_Interface objet.

  */

  template<typename Mesh>
    HarmonicExtension(Mesh& mesh, const Real& diffusion, const QuadRule& Qr, const QuadRule& bdQr, BC_Handler& mesh_BCh);

  //! This method updates the extension of the displacement, i.e. it solves the laplacian proglem

  template<typename Mesh>
    void updateExtension(Mesh& mesh, const Real& time=0.0, const UInt recur=0 );
    template<typename Mesh>


    void updateExtensionTransp(Mesh& mesh, const Real& time=0.0);
  //! This method gives a reference to the computed harmonic extension.
  Vector& getDisplacement();

  //! This method interpolates the mesh velocity when necessary (refFE_u.nbNodes > _mesh.getRefFE().nbNodes)
  template<typename Mesh>
    void interpMeshVelocity(Mesh& mesh, const RefFE& refFE_u, const Dof& dof_u, Vector& wInterp);

  //! Returns a reference to the corresponding Dof object
  const Dof& dofMesh() const;

 protected:

  //! Diffusion coefficient for the laplacian operator
  Real  _diffusion;

  //! Quadrature rule for volumic elementary computations
  const QuadRule& _Qr;

  //! Quadrature rule for surface elementary computations
  const QuadRule& _bdQr;

  //! BC holding the imposed boundary displacement
  BC_Handler& _mesh_BCh;

  //! The Dof object associated with the displacement computations
  Dof _dof_mesh;

  //! The pattern of the FE discretized scalar laplacian matrix
  MSRPatt _aPattB;

  //! The pattern of the FE discretized vector laplacian matrix
  MixedPattern<nDimensions,nDimensions,MSRPatt> _aPatt;

  //! The matrix holding the values
  MixedMatr<nDimensions,nDimensions,MSRPatt,Real> _a;

  //! Current element
  CurrentFE _fe;

  //! Current element  on the boundary
  CurrentBdFE _feBd;

  //! Elementary matrix : 3 blocks
  ElemMat _elmat;

  //! The actual extension of the displacement
  PhysVectUnknown<Vector> _disp;

  //! Auxiliary vector holding the second right hand of the system
  PhysVectUnknown<Vector> _f;

};


//! Constructor for an harmonics extensions
/*!
  \param mesh the mesh of the reference domain to be moved
  \param diffusion the "viscosity" in the laplacian operator
  \param Qr the quadrature rule for volumic elementary computations
  \param bdQr the quadrature rule for surface elementary computations
  \param bch the list of boundary conditions involved in the harmonic extension

  \note The BC_Handler objet (bch) holds the displacement imposed on moving boundary
  in the mesh trhough a BCVetor_Interface objet.

*/
template<typename Mesh> HarmonicExtension::
HarmonicExtension(Mesh& mesh, const Real& diffusion, const QuadRule& Qr, const QuadRule& bdQr, BC_Handler& mesh_BCh):
     _diffusion(diffusion),
     _Qr(Qr),
     _bdQr(bdQr),
     _mesh_BCh(mesh_BCh),
     _dof_mesh(mesh,mesh.getRefFE()),
     _aPattB(_dof_mesh),
     _aPatt(_aPattB,"diag"),
     _a(_aPatt),
     _fe(mesh.getRefFE(), mesh.getGeoMap(), _Qr),
     _feBd(mesh.getRefFE().boundaryFE(), mesh.getGeoMap().boundaryMap(), _bdQr),
     _elmat(_fe.nbNode, nDimensions,  nDimensions),
     _disp(_dof_mesh.numTotalDof()),
     _f(_dof_mesh.numTotalDof()) {

  // Loop on elements
  for(UInt i = 1; i<=mesh.numVolumes(); ++i){

    // Updating derivatives
    _fe.updateFirstDerivQuadPt(mesh.volumeList(i));
    _elmat.zero();
    stiff(_diffusion, _elmat, _fe, 0, 0, 3);
    // Assembling
    for (UInt j=0; j<3; ++j) {
      assemb_mat(_a, _elmat, _fe, _dof_mesh, j, j);
    }
  }

  // Initializations
  _disp  = 0.0;
}

// This method updates the extension of the displacement, i.e. it solves the laplacian proglem
template<typename Mesh>
void HarmonicExtension::updateExtension(Mesh& mesh, const Real& time, const UInt recur) {


  if ( !_mesh_BCh.bdUpdateDone() ) {
    // BC boundary information update
    _mesh_BCh.bdUpdate(mesh, _feBd, _dof_mesh);

    // Boundary conditions treatment on the matrix
    bc_manage_matrix(_a, mesh, _dof_mesh, _mesh_BCh, _feBd, 1.0);
  }

  // Number to total dof
  UInt dim = _dof_mesh.numTotalDof();

  // Initializations
  _f=0.0;

  // Boundary conditions treatment
  bc_manage_vector(_f, mesh, _dof_mesh, _mesh_BCh, _feBd, time,1.0);

  // AZTEC stuff
  int    proc_config[AZ_PROC_SIZE];// Processor information:
  int    options[AZ_OPTIONS_SIZE]; // Array used to select solver options.
  double params[AZ_PARAMS_SIZE];   // User selected solver paramters.
  int    *data_org;                // Array to specify data layout
  double status[AZ_STATUS_SIZE];   // Information returned from AZ_solve()
                                   // indicating success or failure.
  int    *update,                  // vector elements updated on this node.
         *external;                // vector elements needed by this node.
  int    *update_index;            // ordering of update[] and external[]
  int    *extern_index;            // locally on this processor.
  int    N_update;                 // # of unknowns updated on this node

  AZ_set_proc_config(proc_config, AZ_NOT_MPI );
  AZ_read_update(&N_update, &update, proc_config, dim, 1, AZ_linear);
  AZ_transform(proc_config, &external, (int *)_aPattB.giveRaw_bindx(), _a.giveRaw_value(0,0),
	       update, &update_index, &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);
  AZ_transform(proc_config, &external, (int *)_aPattB.giveRaw_bindx(), _a.giveRaw_value(1,1),
	       update, &update_index, &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);
  AZ_transform(proc_config, &external, (int *)_aPattB.giveRaw_bindx(), _a.giveRaw_value(2,2),
	       update, &update_index, &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
	       AZ_MSR_MATRIX);

  // Recovering AZTEC defaults options and params
  AZ_defaults(options,params);

  // Fixed Aztec options for this linera system
  options[AZ_solver]   = AZ_gmres;
  options[AZ_output]   = AZ_warnings;
  options[AZ_poly_ord] = 5;
  options[AZ_kspace]   = 40;
  options[AZ_precond]  = AZ_dom_decomp;
  params[AZ_tol]       = 1.00e-10;
  params[AZ_drop]      = 1.00e-4;
  params[AZ_ilut_fill] = 5;

  options[AZ_recursion_level] = recur;

  // Solving first component
  AZ_solve(_disp.giveVec(), _f.giveVec(), options, params, NULL,
	   (int *)_aPattB.giveRaw_bindx(), NULL, NULL, NULL,
	   _a.giveRaw_value(0,0), data_org, status, proc_config);

  // Solving second component
  AZ_solve(_disp.giveVec()+dim, _f.giveVec()+dim, options, params, NULL,
	   (int *)_aPattB.giveRaw_bindx(), NULL, NULL, NULL,
	   _a.giveRaw_value(1,1), data_org, status, proc_config);

  // Solving third component
  AZ_solve(_disp.giveVec()+2*dim, _f.giveVec()+2*dim, options, params, NULL,
           (int *)_aPattB.giveRaw_bindx(), NULL, NULL, NULL,
	   _a.giveRaw_value(2,2), data_org, status, proc_config);

}



// This method updates the extension of the displacement, i.e. it solves the laplacian proglem
template<typename Mesh>
void HarmonicExtension::updateExtensionTransp(Mesh& mesh, const Real& time) {

  // Boundary conditions treatment
  bc_manage_vector(_disp, mesh, _dof_mesh, _mesh_BCh, _feBd, time);
}
}
#endif
