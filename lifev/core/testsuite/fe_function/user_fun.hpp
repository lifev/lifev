//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/**
   @file
   @author A. Fumagalli <alessio.fumagalli@mail.polimi.it>
   @date 2010-07-29
*/

#ifndef __user_fun_H
#define __user_fun_H 1

#include "fefct.hpp"

namespace dataProblem
{

using namespace LifeV;

// ===================================================
//!                     Typedef
// ===================================================

typedef LinearTetra geoElement_Type;
typedef RegionMesh < geoElement_Type > regionMesh_Type;
typedef boost::shared_ptr < regionMesh_Type > regionMeshPtr_Type;

typedef MapEpetra map_Type;

typedef Exporter < regionMesh_Type > exporter_Type;
typedef boost::shared_ptr < exporter_Type > exporterPtr_Type;

typedef FESpace < regionMesh_Type, map_Type > FESpace_Type;
typedef boost::shared_ptr < FESpace_Type > FESpacePtr_Type;

typedef FEVectorField < regionMesh_Type, map_Type > FEVectorField_Type;
typedef boost::shared_ptr < FEVectorField_Type > FEVectorFieldPtr_Type;

typedef FEScalarField < regionMesh_Type, map_Type > FEScalarField_Type;
typedef boost::shared_ptr < FEScalarField_Type > FEScalarFieldPtr_Type;

typedef FEFunction < regionMesh_Type, map_Type, Vector > FEVectorFct_Type;
typedef boost::shared_ptr < FEVectorFct_Type > FEVectorFctPtr_Type;

typedef FEFunction < regionMesh_Type, map_Type, Real > FEScalarFct_Type;
typedef boost::shared_ptr < FEScalarFct_Type > FEScalarFctPtr_Type;

// ===================================================
//!                       Data
// ===================================================

class MyFun : public FEScalarFct_Type
{
public:
    virtual Real eval ( const UInt& iElem,
                        const FEScalarFct_Type::point_Type& P,
                        const Real& time = 0. ) const;
};

} // Namespace dataProblem

#endif /* __user_fun_H */
