/* -*- mode: c++ -*-
   This program is part of the LifeV library
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
#include "dataDarcy.hpp"
#include "dataAztec.hpp"
#include "readMesh3D.hpp"
#include "refFE.hpp"
#include "refHdivFE.hpp"
#include "refHybridFE.hpp"
#include "elemOper.hpp"
#include "bcHandler.hpp"
#include "dof.hpp"
#include "pattern.hpp"
#include "values.hpp"
#include "assemb.hpp"
#include "vtk_wrtrs.hpp"
#include "user_fct.hpp"
#include "user_diffusion.hpp"

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
class DarcyHandler:
  public DataDarcy,
  public DataAztec
{
public:
  const UInt nbCoor; //!< = 3 in 3D, 2 in 2D
  const GeoMap& geoMap;
  const QuadRule& qr;
  const GeoMap& geoMapBd;
  const QuadRule& qrBd;
  const RefFE& refBdFE;
  const RefHdivFE& refVFE; //!< finite element for u (RT0)
  const RefFE &  refPFE; //!< finite element for p (Q0)
  const RefHybridFE &  refTPFE; //!< finite element for TP (Q0 on faces)
  const RefHybridFE &  refVdotNFE; //!< finite element for V dot N on faces

  CurrentHdivFE vfe;
  CurrentFE pfe;
  CurrentBdFE feBd;
  Dof vdof; //! the degree of freedom for the velocity (RT0)
  Dof pdof; //! the degree of freedom for the pressure (Q0)
  Dof tpdof;//! the degree of freedom for the trace of the pressure (Q0)
  UInt dimPdof; //! number of pressure dof
  UInt dimVdof; //! number of velocity dof
  UInt dimTPdof;//! number of trace of pressure dof
  UInt numFacesPerVolume; //! number of faces per volume
  RegionMesh3D<LinearHexa> mesh; // the mesh
  int nb_bc;                //!< number of boundary conditions
  BCHandler bc;            //!< boundary conditions handler

  //! boundary conditions functions
  //! on a cube one might use them as follows:
  BCFunctionBase bc_fct1;  //!< low  X
  BCFunctionBase bc_fct2;  //!< high X
  BCFunctionBase bc_fct3;  //!< low  Y
  BCFunctionBase bc_fct4;  //!< high Y
  BCFunctionBase bc_fct5;  //!< low  Z
  BCFunctionBase bc_fct6;  //!< high Z

  BCFunctionMixte bc_fct_rob;  //!< a mixte (or Robin) bc function

public:
  DarcyHandler(const GetPot& data_file);
};
}
#endif

