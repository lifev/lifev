/* -*- mode: c++ -*-
   This program is part of the LifeV library

  Author(s): Vincent Martin <vincent.martin@mate.polimi.it>
       Date: 2004-10-11

  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef _DARCYHANDLER_H_
#define _DARCYHANDLER_H_
#include <life/lifesolver/dataDarcy.hpp>
#include <life/lifealg/dataAztec.hpp>
#include <life/lifefilters/readMesh3D.hpp>
#include <life/lifefem/refFE.hpp>
#include <life/lifefem/refHdivFE.hpp>
#include <life/lifefem/refHybridFE.hpp>
#include <life/lifefem/elemOper.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/dof.hpp>
#include <life/lifearray/pattern.hpp>
#include <life/lifefem/values.hpp>
#include <life/lifefem/assemb.hpp>
#include <life/lifefilters/vtk_wrtrs.hpp>

namespace LifeV
{
/*
  \brief Basic objects for a Darcy solver using Mixed Hybrid finite elements
  \file darcyHandler.h
  \author J.-F. Gerbeau and V. Martin
  \data 11/2002

  \par 1) The equation:

  K^{-1} u + grad p = f
  div u  = 0

  (K : diffusion)

  Boundary conditions:

  p = p_d         (Essential B.C.)
  grad p . n = g (Natural B.C.)


  \par 2) Approximation:

  u is approximated in RT_k
  p is approximated in Q_k

  Note: we have just tested the case k=0,
  but the general case could be considered

  \par 3) Method of resolution: hybridization

  We relax the continuity of the flux on the faces of the elements, and we
  impose this continuity with lagrange multipliers denoted by TP (which
  correspond to the trace of the pressure on the faces).

  The linear system (sym. pos. def.) on TP is assemble by eliminating
  at the *element level* the other unknows. For efficiency, all the
  manipulations done on the element matrices are performed
  using BLAS and LAPACK.

  Once the TP problem solved, U and P are recovered by manipulation
  at the *element level*, one more time using BLAS and LAPACK.

  Note that the degree of freedom stored in U are the fluxes
  \int_\Sigma u\cdot n through the faces of the element. The result
  of this computation is therefore specially suitable in a finite
  volume framework, or with discontinous finite element.
*/
template <typename Mesh>
class DarcyHandler:
        public DataDarcy<Mesh>,
        public DataAztec
{
public:

    //! Constructor
    /*!
      \param data_file GetPot data file
      \param refFE_u reference FE for the RTk velocity
      \param refFE_p reference FE for the Pk/Qk pressure
      \param refFE_tp reference FE for the "hybrid" trace of pressure
      \param refFE_vdotn reference FE for the velocity on the faces
      \param refFE_pnodal reference FE for the P1/Q1 pressure (used for post process)
      \param Qr_u volumic quadrature rule
      \param bdQr_u surface quadrature rule
    */
    DarcyHandler( const GetPot& data_file, const RefHdivFE& refFE_u,
                  const RefFE& refFE_p, const RefHybridFE& refFE_tp,
                  const RefHybridFE& refFE_vdotn, const RefFE& refFE_pnodal,
                  const QuadRule& qr_u, const QuadRule& bdqr_u );


protected:
    const UInt nbCoor; //!< = 3 in 3D, 2 in 2D
    const GeoMap& geoMap;
    const QuadRule& qr;
    const GeoMap& geoMapBd;
    const QuadRule& qrBd;
    const RefHdivFE& refVFE; //!< finite element for u (RT0 hexa/tetra)
    const RefFE &  refPFE; //!< finite element for p (Q0 or P0)
    const RefHybridFE &  refTPFE; //!< finite element for TP (Q0 or P0 on faces)
    const RefHybridFE &  refVdotNFE; //!< finite element for V dot N on faces
    const RefFE& refBdFE; //!< finite element on the boundary (Q0 or P0)
    const RefFE& refPFEnodal; //!< nodal finite element for p (Q1 or P1) (for post process only)

    CurrentHdivFE vfe;
    CurrentFE pfe;
    CurrentBdFE feBd;
    Dof vdof; //! the degree of freedom for the velocity (RT0 hexa/tetra)
    Dof pdof; //! the degree of freedom for the pressure (Q0 or P0)
    Dof tpdof;//! the degree of freedom for the trace of the pressure (Q0 or P0)
    UInt dimPdof; //! number of pressure dof
    UInt dimVdof; //! number of velocity dof
    UInt dimTPdof;//! number of trace of pressure dof
    UInt numFacesPerVolume; //! number of faces per volume


};

// --------------------------------------------------
// IMPLEMENTATION
// --------------------------------------------------

// constructor
template <typename Mesh>
DarcyHandler<Mesh>::DarcyHandler( const GetPot& data_file, const RefHdivFE& refFE_u,
                                  const RefFE& refFE_p, const RefHybridFE& refFE_tp,
                                  const RefHybridFE& refFE_vdotn, const RefFE& refFE_pnodal,
                                  const QuadRule& qr_u, const QuadRule& bdqr_u):
    DataDarcy<Mesh>(data_file),
    DataAztec(data_file, "darcy/aztec"),
    nbCoor(nDimensions),
    geoMap( getGeoMap( this->_mesh ) ),
    qr( qr_u ),
    geoMapBd( geoMap.boundaryMap() ),
    qrBd( bdqr_u ),
    refVFE( refFE_u ),
    refPFE( refFE_p ),
    refTPFE( refFE_tp ),
    refVdotNFE( refFE_vdotn ),
    refBdFE( refPFE.boundaryFE() ),
    refPFEnodal( refFE_pnodal ),
    vfe( refVFE , geoMap , qr ),
    pfe( refPFE , geoMap , qr ),
    feBd( refBdFE , geoMapBd , qrBd ),
    vdof( refVFE ),
    pdof( refPFE ),
    tpdof( refTPFE ),
    numFacesPerVolume( this->_mesh.volumeList(1).numLocalFaces )
        /* we assume that all element have the same number
           of faces, so we just look at the first one */
{
    //! This function should be already called in DataMesh construction
    // this->_mesh.updateElementFaces(true);
    /*the "true" flag is to build the faceList
      of all faces (and not only the boundary) */

    if(this->verbose>2) this->_mesh.showMe();
    // build dof
    vdof.update(this->_mesh);
    pdof.update(this->_mesh);
    tpdof.update(this->_mesh);

    dimTPdof = tpdof.numTotalDof();
    dimPdof = pdof.numTotalDof();
    dimVdof = vdof.numTotalDof();

    if(verbose>0){
        std::cout << "Number of TP dof : " << dimTPdof << std::endl;
        std::cout << "Number of  P dof : " << dimPdof << std::endl;
        std::cout << "Number of  V dof : " << dimVdof << std::endl;
    }

    /*
    // check the mesh after b.c.
      BE CAREFUL: calling
      this->_mesh.check(true,true);
      after updateElementFaces(true) erase faceElement !!!!
    */
    //
    if(verbose>2) tpdof.showMe();
}

}
#endif

