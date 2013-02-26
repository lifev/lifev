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

/*!
    @file
    @brief File for the implementation of the QuadratureRule class.

    @author Jean-Frederic Gerbeau
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01-06-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#include <lifev/core/fem/QuadratureRule.hpp>

namespace LifeV
{

Real QuadratureRule::S_exactnessTol = 1e-8;

// ===================================================
// Constructors & Destructor
// ===================================================

QuadratureRule::QuadratureRule()
    :    M_pt (0), M_shape (NONE), M_name (""), M_nbQuadPt (0), M_degOfExact (0), M_dimension (0)
{}

QuadratureRule::QuadratureRule ( const QuadraturePoint* pt, int /*id*/, std::string name,
                                 ReferenceShapes shape, UInt nbQuadPt, UInt degOfExact ) :
    M_pt ( nbQuadPt),
    M_shape ( shape ), M_name ( name ),
    M_nbQuadPt ( nbQuadPt ), M_degOfExact ( degOfExact )
{
    M_dimension = 3;
    for (UInt i (0); i < nbQuadPt; ++i)
    {
        M_pt[i] = pt[i];
    }
}

QuadratureRule::QuadratureRule ( const QuadratureRule& qr ) :
    M_pt ( qr.M_pt ), M_shape ( qr.M_shape ), M_name ( qr.M_name ),
    M_nbQuadPt ( qr.M_nbQuadPt ), M_degOfExact ( qr.M_degOfExact ), M_dimension (qr.M_dimension)
{
}

QuadratureRule::QuadratureRule ( const QuadratureRule& qr, const UInt dim) :
    M_pt ( qr.M_nbQuadPt ), M_shape ( qr.M_shape ), M_name ( qr.M_name ),
    M_nbQuadPt ( qr.M_nbQuadPt ), M_degOfExact ( qr.M_degOfExact ), M_dimension (dim)
{
    ASSERT (dim >= shapeDimension (M_shape), " Downgrading quadrature rule is forbidden ")

    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        M_pt[i] = QuadraturePoint (qr.M_pt[i], dim);
    }

}

QuadratureRule::QuadratureRule (std::string name, ReferenceShapes shape, UInt dimension, UInt degreeOfExactness, UInt nbQuadPt, ... ) :
    M_pt ( nbQuadPt),
    M_shape ( shape),
    M_name ( name),
    M_nbQuadPt ( nbQuadPt),
    M_degOfExact ( degreeOfExactness),
    M_dimension ( dimension)
{
    ASSERT (dimension >= shapeDimension (shape), " Downgrading quadrature rule is forbidden ");

    va_list quadList;
    va_start (quadList, nbQuadPt);
    for (UInt iterArg (0); iterArg < nbQuadPt; ++iterArg)
    {
        GeoVector nextPoint (dimension);
        for (UInt i (0); i < dimension ; ++i)
        {
            nextPoint[i] = va_arg (quadList, Real);
        };
        M_pt[iterArg] = QuadraturePoint (nextPoint, va_arg (quadList, double) );
    }
    va_end (quadList);
}

QuadratureRule::~QuadratureRule()
{
}

// ===================================================
// Operators
// ===================================================


std::ostream& operator << ( std::ostream& c, const QuadratureRule& qr )
{
    c << " name: " << qr.M_name << std::endl;
    c << " shape:" << ( int ) qr.M_shape << std::endl;
    c << " nbQuadPt: " << qr.M_nbQuadPt << std::endl;
    c << " Points: \n";
    for ( UInt i (0); i < qr.M_nbQuadPt; ++i )
    {
        c << qr.M_pt[ i ] << std::endl;
    }
    return c;
}

// ===================================================
// Methods
// ===================================================


void QuadratureRule::showMe ( std::ostream& output) const
{
    output << " Name  : " << M_name << std::endl;
    output << " Shape : " << M_shape << std::endl;
    output << " Size  : " << M_nbQuadPt << std::endl;
    output << " --- Points --- " << std::endl;
    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        output << M_pt[i] << std::endl;
    }
}

UInt QuadratureRule::checkExactness() const
{
    switch (M_shape)
    {
        case LINE:
            return checkExactnessSegment();
            break;

        case TRIANGLE:
            return checkExactnessTriangle();
            break;

        case TETRA:
            return checkExactnessTetra();
            break;

        default:
            std::cerr << " No check for this shape ..." << std::endl;
            return 0;
    }
}

void QuadratureRule::vtkExport ( const std::string& filename) const
{
    std::ofstream output (filename.c_str() );
    ASSERT (!output.fail(), " Unable to open the file for the export of the quadrature ");

    // Header
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "LifeV : Quadrature export" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << M_nbQuadPt << " float" << std::endl;

    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        output << M_pt[i].coor (0) << " " << M_pt[i].coor (1) << " " << M_pt[i].coor (2) << std::endl;
    };

    output << "VERTICES " << M_nbQuadPt << " " << 2 * M_nbQuadPt << std::endl;

    for (UInt i (0); i < M_nbQuadPt; ++i)
    {
        output << 1 << " " << i << std::endl;
    };

    output.close();
}

// ===================================================
// Set Methods
// ===================================================

void QuadratureRule::setPoints (const std::vector<QuadraturePoint>& pts)
{
    M_pt.clear();
    M_pt.resize (pts.size() );
    for (UInt i (0); i < pts.size(); ++i)
    {
        M_pt[i] = QuadraturePoint (pts[i], M_dimension);
    }

    M_nbQuadPt = pts.size();
}

void QuadratureRule::setPoints (const std::vector<GeoVector>& coordinates, const std::vector<Real>& weights)
{
    ASSERT (coordinates.size() == weights.size(), "Non matching length of the arguments");

    M_pt.clear();
    M_pt.resize (coordinates.size() );
    for (UInt i (0); i < M_pt.size(); ++i)
    {
        M_pt[i] = QuadraturePoint (coordinates[i], weights[i], M_dimension);
    }

    M_nbQuadPt = M_pt.size();
}

void QuadratureRule::setName (const std::string& newName)
{
    M_name = newName;
}

void QuadratureRule::setExactness (const UInt& exactness)
{
    M_degOfExact = exactness;
}

void QuadratureRule::setDimensionShape (const UInt& newDim, const ReferenceShapes& newShape)
{
    ASSERT (newDim >= shapeDimension (newShape), " Impossible shape-dimension combinaison ");
    M_dimension = newDim;
    M_shape = newShape;

    // Change also the dimension of the points if they are already set!
    for (UInt i (0); i < M_pt.size(); ++i)
    {
        M_pt[i] = QuadraturePoint (M_pt[i], newDim);
    }
}

// ===================================================
// Private Methods
// ===================================================

UInt QuadratureRule::checkExactnessTetra() const
{
    // Degre 0: f=1 => exact value : 1/6
    Real partialSum (0.0);

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 6.0) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x+y+2z => exact value: 1/6
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) + quadPointCoor (iterQ, 1) + 2 * quadPointCoor (iterQ, 2) ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 6.0) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x + y*z + y*y => exact value: 1/24
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 2)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 24.0) > S_exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*x + y*y + y*y*z + x => exact value: 5/72
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 2)
                       + quadPointCoor (iterQ, 0)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 5.0 / 72.0) > S_exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*x +y*y*z*z + z*z*z => exact value: 1/72
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 2) * quadPointCoor (iterQ, 2)
                       + quadPointCoor (iterQ, 2) * quadPointCoor (iterQ, 2) * quadPointCoor (iterQ, 2)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 72.0) > S_exactnessTol)
    {
        return 3;
    }

    return 4;

}

UInt QuadratureRule::checkExactnessTriangle() const
{

    // Degre 0: f=1 => exact value : 1/2
    Real partialSum (0.0);

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += weight (iterQ);
    }

    if ( std::fabs (partialSum - 0.5) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x+y => exact value : 1/3
    partialSum = 0.0;

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) + quadPointCoor (iterQ, 1) )
                      * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 3.0) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x+3*x*y => exact value : 5/24
    partialSum = 0.0;

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       + 3 * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1) )
                      * weight (iterQ);
    }

    if ( std::fabs (partialSum - 5.0 / 24.0) > S_exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*y-5*x*y => exact value :-23/120
    partialSum = 0.0;

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1)
                       - 5 * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1) )
                      * weight (iterQ);
    }

    if ( std::fabs (partialSum + 23.0 / 120.0) > S_exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*y - x*y*y + y*y*y*y => exact value : 1/40
    partialSum = 0.0;

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1)
                       - quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1)
                       + quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1) )
                      * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 40.0) > S_exactnessTol)
    {
        return 3;
    }

    // Degre 5: f= x^4*y - x*y^2 + y^5 => exact value : 1/84
    partialSum = 0.0;

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += ( std::pow (quadPointCoor (iterQ, 0), 4.0) * quadPointCoor (iterQ, 1)
                        - quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 1) * quadPointCoor (iterQ, 1)
                        + std::pow (quadPointCoor (iterQ, 1), 5.0)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 84.0) > S_exactnessTol)
    {
        return 4;
    }


    return 5;
}

UInt QuadratureRule::checkExactnessSegment() const
{
    // Degre 0: f=1 => exact value : 1
    Real partialSum (0.0);

    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x => exact value: 1/2
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 2.0) > S_exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x => exact value: 1/3
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 3.0) > S_exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*x => exact value: 1/4
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 4.0) > S_exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*x => exact value: 1/5
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 5.0) > S_exactnessTol)
    {
        return 3;
    }

    // Degre 5: f=x*x*x*x*x => exact value: 1/6
    partialSum = 0.0;
    for (UInt iterQ (0); iterQ < M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       * quadPointCoor (iterQ, 0) * quadPointCoor (iterQ, 0)
                       * quadPointCoor (iterQ, 0)
                      ) * weight (iterQ);
    }

    if ( std::fabs (partialSum - 1.0 / 6.0) > S_exactnessTol)
    {
        return 4;
    }

    return 5;
}



}
