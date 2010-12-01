/*-*- mode: c++ -*-
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

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
#include <life/lifefem/quadRule.hpp>

namespace LifeV
{

Real QuadRule::exactnessTol=1e-8;

QuadRule::QuadRule()
        :    M_pt(0), M_shape(NONE), M_name(""), M_nbQuadPt(0), M_degOfExact(0), M_dimension(0)
{}

QuadRule::QuadRule( const QuadPoint* pt, int /*id*/, std::string name,
                    ReferenceShapes shape, UInt nbQuadPt, UInt degOfExact ) :
        M_pt( nbQuadPt),
        M_shape( shape ), M_name( name ),
        M_nbQuadPt( nbQuadPt ), M_degOfExact( degOfExact )
{
    CONSTRUCTOR( "QuadRule" );
    M_dimension = 3;
    for (UInt i(0); i<nbQuadPt; ++i)
    {
        M_pt[i]=pt[i];
    }
}

QuadRule::QuadRule( const QuadRule& qr ) :
        M_pt( qr.M_pt ), M_shape( qr.M_shape ), M_name( qr.M_name ),
        M_nbQuadPt( qr.M_nbQuadPt ), M_degOfExact( qr.M_degOfExact ), M_dimension(qr.M_dimension)
{
    CONSTRUCTOR( "QuadRule" );
}

QuadRule::QuadRule( const QuadRule& qr, const UInt dim) :
        M_pt( qr.M_nbQuadPt ), M_shape( qr.M_shape ), M_name( qr.M_name ),
        M_nbQuadPt( qr.M_nbQuadPt ), M_degOfExact( qr.M_degOfExact ), M_dimension(dim)
{
    ASSERT(dim >= getReferenceDimension(M_shape)," Downgrading quadrature rule is forbidden ")

    for (UInt i(0); i<M_nbQuadPt; ++i)
    {
        M_pt[i] = QuadPoint(qr.M_pt[i],dim);
    }

    CONSTRUCTOR( "QuadRule" );
}

QuadRule::QuadRule(std::string name, ReferenceShapes shape, UInt dimension, UInt degreeOfExactness, UInt nbQuadPt, ... ) :
        M_pt( nbQuadPt),
        M_shape( shape),
        M_name( name),
        M_nbQuadPt( nbQuadPt),
        M_degOfExact( degreeOfExactness),
        M_dimension( dimension)
{
    ASSERT(dimension >= getReferenceDimension(shape)," Downgrading quadrature rule is forbidden ");

    CONSTRUCTOR( "QuadRule" );
    va_list quadList;
    va_start(quadList,nbQuadPt);
    for (UInt iterArg(0); iterArg<nbQuadPt; ++iterArg)
    {
        GeoVector nextPoint(dimension);
        for (UInt i(0); i<dimension ; ++i)
        {
            nextPoint[i] = va_arg(quadList,Real);
        };
        M_pt[iterArg] = QuadPoint(nextPoint,va_arg(quadList,double));
    }
    va_end(quadList);
}

QuadRule::~QuadRule()
{
    DESTRUCTOR( "QuadRule" );
}

void QuadRule::showMe( std::ostream& output) const
{
    output << " Name  : " << M_name << std::endl;
    output << " Shape : " << M_shape << std::endl;
    output << " Size  : " << M_nbQuadPt << std::endl;
    output << " --- Points --- " << std::endl;
    for (UInt i(0); i<M_nbQuadPt; ++i)
    {
        output << M_pt[i] << std::endl;
    }
}

UInt QuadRule::checkExactness() const
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

UInt QuadRule::checkExactnessTetra() const
{
    // Degre 0: f=1 => exact value : 1/6
    Real partialSum(0.0);

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/6.0) > exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x+y+2z => exact value: 1/6
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0) + quadPointCoor(iterQ,1) + 2*quadPointCoor(iterQ,2) )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/6.0) > exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x + y*z + y*y => exact value: 1/24
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,2)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)
                      )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/24.0) > exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*x + y*y + y*y*z + x => exact value: 5/72
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,2)
                       + quadPointCoor(iterQ,0)
                      )*weight(iterQ);
    }

    if ( std::abs(partialSum - 5.0/72.0) > exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*x +y*y*z*z + z*z*z => exact value: 1/72
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,2)*quadPointCoor(iterQ,2)
                       + quadPointCoor(iterQ,2)*quadPointCoor(iterQ,2)*quadPointCoor(iterQ,2)
                      )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/72.0) > exactnessTol)
    {
        return 3;
    }

    return 4;

}

UInt QuadRule::checkExactnessTriangle() const
{

    // Degre 0: f=1 => exact value : 1/2
    Real partialSum(0.0);

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += weight(iterQ);
    }

    if ( std::abs(partialSum - 0.5) > exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x+y => exact value : 1/3
    partialSum=0.0;

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0) + quadPointCoor(iterQ,1))
                      * weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/3.0) > exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x+3*x*y => exact value : 5/24
    partialSum=0.0;

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       + 3*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1))
                      * weight(iterQ);
    }

    if ( std::abs(partialSum - 5.0/24.0) > exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*y-5*x*y => exact value :-23/120
    partialSum=0.0;

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1)
                       - 5*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1))
                      * weight(iterQ);
    }

    if ( std::abs(partialSum + 23.0/120.0) > exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*y - x*y*y + y*y*y*y => exact value : 1/40
    partialSum=0.0;

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1)
                       - quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)
                       + quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1) )
                      * weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/40.0) > exactnessTol)
    {
        return 3;
    }

    // Degre 5: f= x^4*y - x*y^2 + y^5 => exact value : 1/84
    partialSum=0.0;

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += ( std::pow(quadPointCoor(iterQ,0),4.0)*quadPointCoor(iterQ,1)
                        -quadPointCoor(iterQ,0)*quadPointCoor(iterQ,1)*quadPointCoor(iterQ,1)
                        + std::pow(quadPointCoor(iterQ,1),5.0)
                      ) * weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/84.0) > exactnessTol)
    {
        return 4;
    }


    return 5;
}

UInt QuadRule::checkExactnessSegment() const
{
    // Degre 0: f=1 => exact value : 1
    Real partialSum(0.0);

    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0) > exactnessTol)
    {
        return 0;
    }

    // Degre 1: f=x => exact value: 1/2
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0))*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/2.0) > exactnessTol)
    {
        return 0;
    }

    // Degre 2: f=x*x => exact value: 1/3
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0) )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/3.0) > exactnessTol)
    {
        return 1;
    }

    // Degre 3: f=x*x*x => exact value: 1/4
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0) )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/4.0) > exactnessTol)
    {
        return 2;
    }

    // Degre 4: f=x*x*x*x => exact value: 1/5
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                      )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/5.0) > exactnessTol)
    {
        return 3;
    }

    // Degre 5: f=x*x*x*x*x => exact value: 1/6
    partialSum=0.0;
    for (UInt iterQ(0); iterQ<M_nbQuadPt; ++iterQ)
    {
        partialSum += (quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       *quadPointCoor(iterQ,0)*quadPointCoor(iterQ,0)
                       *quadPointCoor(iterQ,0)
                      )*weight(iterQ);
    }

    if ( std::abs(partialSum - 1.0/6.0) > exactnessTol)
    {
        return 4;
    }

    return 5;
}

void QuadRule::VTKexport( const std::string& filename) const
{
    std::ofstream output(filename.c_str());
    ASSERT(!output.fail(), " Unable to open the file for the export of the quadrature ");

    // Header
    output << "# vtk DataFile Version 3.0" << std::endl;
    output << "LifeV : Quadrature export" << std::endl;
    output << "ASCII" << std::endl;
    output << "DATASET POLYDATA" << std::endl;
    output << "POINTS " << M_nbQuadPt << " float" << std::endl;

    for (UInt i(0); i< M_nbQuadPt; ++i)
    {
        output << M_pt[i].coor(0) << " " << M_pt[i].coor(1) << " " << M_pt[i].coor(2) << std::endl;
    };

    output << "VERTICES " << M_nbQuadPt << " " << 2*M_nbQuadPt << std::endl;

    for (UInt i(0); i< M_nbQuadPt; ++i)
    {
        output << 1 << " " << i << std::endl;
    };

    output.close();
}

std::ostream& operator << ( std::ostream& c, const QuadRule& qr )
{
    c << " name: " << qr.M_name << std::endl;
    c << " shape:" << ( int ) qr.M_shape << std::endl;
    c << " nbQuadPt: " << qr.M_nbQuadPt << std::endl;
    c << " Points: \n";
    for ( UInt i (0); i < qr.M_nbQuadPt; ++i )
        c << qr.M_pt[ i ] << std::endl;
    return c;
}

void QuadRule::setPoints(const std::vector<QuadPoint>& pts)
{
    M_pt.clear();
    M_pt.resize(pts.size());
    for (UInt i(0); i<pts.size(); ++i)
    {
        M_pt[i]=QuadPoint(pts[i],M_dimension);
    }

    M_nbQuadPt=pts.size();
}

void QuadRule::setPoints(const std::vector<GeoVector>& coordinates, const std::vector<Real>& weights)
{
    ASSERT(coordinates.size()==weights.size(),"Non matching length of the arguments");

    M_pt.clear();
    M_pt.resize(coordinates.size());
    for (UInt i(0); i<M_pt.size(); ++i)
    {
        M_pt[i]=QuadPoint(coordinates[i],weights[i],M_dimension);
    }

    M_nbQuadPt=M_pt.size();
}

void QuadRule::setName(const std::string& newName)
{
    M_name = newName;
}

void QuadRule::setExactness(const UInt& exactness)
{
    M_degOfExact = exactness;
}

void QuadRule::setDimensionShape(const UInt& newDim, const ReferenceShapes& newShape)
{
    ASSERT(newDim >= getReferenceDimension(newShape)," Impossible shape-dimension combinaison ");
    M_dimension = newDim;
    M_shape = newShape;

    // Change also the dimension of the points if they are already set!
    for (UInt i(0); i<M_pt.size(); ++i)
    {
        M_pt[i] = QuadPoint(M_pt[i],newDim);
    }
}

}
