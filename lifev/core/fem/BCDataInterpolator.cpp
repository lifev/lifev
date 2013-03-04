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
 *  @file
 *  @brief File containing a class for interpolating boundary functions from scattered data
 *
 *  @date 09-06-2011
 *  @author Toni Lassila <toni.lassila@epfl.ch>
 *
 *  @maintainer Toni Lassila <toni.lassila@epfl.ch>
 */

#include <lifev/core/fem/BCDataInterpolator.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_array.hpp>
#include <boost/bind.hpp>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
BCDataInterpolator::BCDataInterpolator(  ) :
    M_interpolationMatrix(),
    M_rhs_x(),
    M_rhs_y(),
    M_rhs_z(),
    M_coeffs_x(),
    M_coeffs_y(),
    M_coeffs_z(),
    M_denseSolver(),
    M_interpolationMethod ( RBF_InverseMultiQuadric ),
    M_dataSites(), M_dataValues(),
    M_dataValues_timeSamples(),
    M_nofControlPoints(),
    M_filteringLevel ( 0 ),
    M_lastInterpolatedAtTime ( -1 ),
    M_timePeriod ( 0 ),
    M_timeInterval ( 0 ),
    M_flagInterpolated ( false ),
    M_verbose ( false )
{
    M_userDefinedFunction = boost::bind ( &BCDataInterpolator::interpolatedDataFunction, this, _1, _2, _3, _4, _5 );
}

BCDataInterpolator::~BCDataInterpolator()
{
}

// ===================================================
// Operators
// ===================================================
/*BCDataInterpolator&
 BCDataInterpolator::operator = (const BCDataInterpolator &BCdi)
 {
 if (this != &BCdi)
 {
 }

 return *this;
 }*/

// ===================================================
// Methods
// ===================================================
Real BCDataInterpolator::interpolatedDataFunction ( const Real& t,
                                                    const Real& x,
                                                    const Real& y,
                                                    const Real& z,
                                                    const ID& component )
{
    BCDataInterpolator_point vectEval;

    vectEval = interpolateVectorialFunction ( t, x, y, z );

    switch ( component )
    {
        case 0:
            return vectEval.x;
        case 1:
            return vectEval.y;
        case 2:
            return vectEval.z;
        default:
            std::ostringstream exception;
            exception << "Invalid component: " << component << std::endl;
            throw std::invalid_argument ( exception.str() );
    }
}

void BCDataInterpolator::showMe ( bool verbose, std::ostream& out ) const
{
    out << " Boundary Conditions Data Interpolator ====>" << std::endl;
    out << " Data sites:" << std::endl;

    for (UInt i = 0; i < M_nofControlPoints; i++)
    {
        out << "(" << M_dataSites[i].x << ", " << M_dataSites[i].y << ", " << M_dataSites[i].z << ")" << std::endl;
    }

    if (verbose)
    {
        for (Int t = 0; t < M_timePeriod / M_timeInterval; t++)
        {
            out << " Data values (at sites):" << std::endl;
            for (UInt i = 0; i < M_nofControlPoints; i++)
            {
                out << "(" << M_dataValues_timeSamples[indexInTime (i, t)].x << ", " << M_dataValues_timeSamples[indexInTime (i, t)].y << ", " << M_dataValues_timeSamples[indexInTime (i, t)].z << ")" << std::endl;
            }
        }
    }

    out << " Data values (interpolated):" << std::endl;
    for (UInt i = 0; i < M_nofControlPoints; i++)
    {
        out << "(" << M_dataValues[i].x << ", " << M_dataValues[i].y << ", " << M_dataValues[i].z << ")" << std::endl;
    }
}

void BCDataInterpolator::readData ( const std::string& fileName )
{
    std::ifstream fin;
    UInt RDIM;

    fin.open ( fileName.c_str() );

    if ( !fin.fail() )
    {
        fin.ignore ( 80, '\n' );

        fin >> M_nofControlPoints;
        fin >> RDIM;
        fin >> M_timeInterval;
        fin >> M_timePeriod;
        fin >> M_filteringLevel;

        if ( ( RDIM != 1 ) && ( RDIM != 3 ) )
        {
            std::ostringstream exception;
            exception << "Interpolated data must be either scalar or a 3-vector: " << std::endl;
            throw std::invalid_argument ( exception.str() );
        }

        if ( M_timeInterval > M_timePeriod)
        {
            std::ostringstream exception;
            exception << "Interpolation time interval and/or period are inconsistent." << std::endl;
            throw std::invalid_argument ( exception.str() );
        }

        Int n = static_cast< Int > (floor (M_timePeriod / M_timeInterval) );

        M_dataSites = boost::shared_array<BCDataInterpolator_point> (new BCDataInterpolator_point[M_nofControlPoints]);
        M_dataValues = boost::shared_array<BCDataInterpolator_point> (new BCDataInterpolator_point[M_nofControlPoints]);
        M_dataValues_timeSamples = boost::shared_array<BCDataInterpolator_point> (new BCDataInterpolator_point[M_nofControlPoints * (n + 1)]);

        fin.ignore ( 80, '\n' );
        fin.ignore ( 80, '\n' );

        // Read the data site locations
        for ( UInt i = 0; i < M_nofControlPoints; i++ )
        {
            fin >> M_dataSites[i].x;
            fin >> M_dataSites[i].y;
            fin >> M_dataSites[i].z;
        }

        fin.ignore ( 80, '\n' );
        fin.ignore ( 80, '\n' );

        // Read the data values at time instances
        for (Int t = 0; t <= n; t++)
        {

            for ( UInt i = 0; i < M_nofControlPoints; i++ )
            {
                fin >> M_dataValues_timeSamples[indexInTime (i, t)].x;
                fin >> M_dataValues_timeSamples[indexInTime (i, t)].y;
                fin >> M_dataValues_timeSamples[indexInTime (i, t)].z;
            }

            fin.ignore ( 80, '\n' );
            fin.ignore ( 80, '\n' );
        }

        fin.close();

        // Form the interpolation matrix now
        formRBFMatrix();
    }

}

void BCDataInterpolator::exportInterpolationMatrix() const
{
    // TODO: Write matrix exporting routine
}


// ===================================================
// Private Methods
// ===================================================
void BCDataInterpolator::formRBFMatrix()
{
    Int nofElementsPerRow = ( needSideConstraints() ? M_nofControlPoints + NDIM + 1 : M_nofControlPoints );

    M_interpolationMatrix.Shape ( static_cast< int > ( nofElementsPerRow ),
                                  static_cast< int > ( nofElementsPerRow ) );

    Real rbfEval;

    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {

        for ( UInt j = 0; j <= i; j++ )
        {
            // Symmetric interpolation matrix [A]_ij = phi(|x_i - x_j|)

            rbfEval = evaluateRBF ( M_dataSites[i], M_dataSites[j] );
            M_interpolationMatrix ( static_cast< int > ( i ), static_cast< int > ( j ) ) = static_cast< double > ( rbfEval );
            M_interpolationMatrix ( static_cast< int > ( j ), static_cast< int > ( i ) ) = static_cast< double > ( rbfEval );

        }

    }

    if ( needSideConstraints() )
    {

        // Add global linear polynomial to the interpolation basis

        for ( UInt i = 0; i < M_nofControlPoints; i++ )
        {

            // First order term of the polynomial
            M_interpolationMatrix ( i, M_nofControlPoints )     = M_dataSites[i].x;
            M_interpolationMatrix ( i, M_nofControlPoints + 1 ) = M_dataSites[i].y;
            M_interpolationMatrix ( i, M_nofControlPoints + 2 ) = M_dataSites[i].z;

            M_interpolationMatrix ( M_nofControlPoints,     i ) = M_dataSites[i].x;
            M_interpolationMatrix ( M_nofControlPoints + 1, i ) = M_dataSites[i].y;
            M_interpolationMatrix ( M_nofControlPoints + 2, i ) = M_dataSites[i].z;

            // Constant term of the polynomial
            M_interpolationMatrix ( i, M_nofControlPoints + 3 ) = 1;
            M_interpolationMatrix ( M_nofControlPoints + 3, i ) = 1;

        }
    }

    M_denseSolver.SetMatrix ( M_interpolationMatrix );
    M_flagInterpolated = false;

}

void BCDataInterpolator::solveInterpolationSystem()
{

#ifdef DEBUG
    // Check that condition number is "reasonable" (< 1e6), if not print a warning
    Real reciprocalCond = M_denseSolver.RCOND();

    if (reciprocalCond < 1e-6)
    {
        debugStream ( 5000 ) << "Radial basis interpolation matrix is ill-conditioned.\n";
        debugStream ( 5000 ) << "Estimated reciprocal of condition number: " << reciprocalCond << "\n";
    }
#endif

    // Solve the interpolation system for each component using LU-factorization

    M_denseSolver.SetVectors ( M_coeffs_x, M_rhs_x );
    M_denseSolver.Solve();

    M_denseSolver.SetVectors ( M_coeffs_y, M_rhs_y );
    M_denseSolver.Solve();

    M_denseSolver.SetVectors ( M_coeffs_z, M_rhs_z );
    M_denseSolver.Solve();
}

BCDataInterpolator::BCDataInterpolator_point BCDataInterpolator::interpolateVectorialFunction ( const Real& t,
        const Real& x,
        const Real& y,
        const Real& z )
{
    BCDataInterpolator_point rbfEval, evalPoint;
    Real rbfShape;

    rbfEval.x = 0;
    rbfEval.y = 0;
    rbfEval.z = 0;
    evalPoint.x = x;
    evalPoint.y = y;
    evalPoint.z = z;

    // Data values need to be interpolated in time, if not already done for this t
    if (M_lastInterpolatedAtTime != t)
    {
        interpolateDataValuesInTime ( t );
    }

    // If not yet solved the interpolation system, do so
    if ( !M_flagInterpolated )
    {
        solveInterpolationSystem();
        M_flagInterpolated = true;
    }

    // For each data site...
    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {

        rbfShape = evaluateRBF ( evalPoint, M_dataSites[i] );

        // Interpolated function value at point (x,y,z)
        rbfEval.x += M_coeffs_x ( i ) * rbfShape;
        rbfEval.y += M_coeffs_y ( i ) * rbfShape;
        rbfEval.z += M_coeffs_z ( i ) * rbfShape;

    } // ...for each data site


    if ( needSideConstraints() )
    {
        // Add the linear term of the polynomial...
        rbfEval.x += M_coeffs_x (     M_nofControlPoints ) * evalPoint.x;
        rbfEval.x += M_coeffs_x ( 1 + M_nofControlPoints ) * evalPoint.y;
        rbfEval.x += M_coeffs_x ( 2 + M_nofControlPoints ) * evalPoint.z;

        rbfEval.y += M_coeffs_y (     M_nofControlPoints ) * evalPoint.x;
        rbfEval.y += M_coeffs_y ( 1 + M_nofControlPoints ) * evalPoint.y;
        rbfEval.y += M_coeffs_y ( 2 + M_nofControlPoints ) * evalPoint.z;

        rbfEval.z += M_coeffs_z (     M_nofControlPoints ) * evalPoint.x;
        rbfEval.z += M_coeffs_z ( 1 + M_nofControlPoints ) * evalPoint.y;
        rbfEval.z += M_coeffs_z ( 2 + M_nofControlPoints ) * evalPoint.z;

        // ...and the constant term of the polynomial
        rbfEval.x += M_coeffs_x ( 3 + M_nofControlPoints );
        rbfEval.y += M_coeffs_y ( 3 + M_nofControlPoints );
        rbfEval.z += M_coeffs_z ( 3 + M_nofControlPoints );
    }

    return rbfEval;
}

Real BCDataInterpolator::evaluateRBF ( const BCDataInterpolator_point point1,
                                       const BCDataInterpolator_point point2 )
{
    Real r = std::sqrt ( std::pow ( point1.x - point2.x, 2 ) + std::pow ( point1.y - point2.y, 2 ) + std::pow ( point1.z - point2.z, 2 ) );

    switch ( M_interpolationMethod )
    {

        case RBF_ThinPlateSpline:
            return ( r < 1e-3 ? 0.0 : std::pow ( r, 2 ) * std::log ( r ) );

        case RBF_MultiQuadric:
            return std::sqrt ( std::pow ( r, 2 ) + 1 );

        case RBF_Cubic:
            return std::pow ( r, 3 );

        case RBF_Gaussian:
            return std::exp ( -std::pow ( r, 2 ) );

        case RBF_InverseMultiQuadric:
            return 1 / std::sqrt (std::pow ( r, 2 ) + 1);

        default:
            return 0.0;
    }
}

bool BCDataInterpolator::needSideConstraints() const
{
    switch ( M_interpolationMethod )
    {

        case RBF_MultiQuadric:
        case RBF_Gaussian:
        case RBF_InverseMultiQuadric:
            return false;

        default:
            return true;
    }
}

void BCDataInterpolator::formRBFvectors()
{
    Int nofElementsPerRow = ( needSideConstraints() ? 2 * M_nofControlPoints + NDIM : M_nofControlPoints );

    M_rhs_x.Size ( static_cast< int > ( nofElementsPerRow ) );
    M_rhs_y.Size ( static_cast< int > ( nofElementsPerRow ) );
    M_rhs_z.Size ( static_cast< int > ( nofElementsPerRow ) );

    M_coeffs_x.Size ( static_cast< int > ( nofElementsPerRow ) );
    M_coeffs_y.Size ( static_cast< int > ( nofElementsPerRow ) );
    M_coeffs_z.Size ( static_cast< int > ( nofElementsPerRow ) );

    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {
        // Form RHS for the values at the RBF sites

        M_rhs_x ( i ) = M_dataValues[i].x;
        M_rhs_y ( i ) = M_dataValues[i].y;
        M_rhs_z ( i ) = M_dataValues[i].z;
    }

    M_flagInterpolated = false;
}

void BCDataInterpolator::interpolateDataValuesInTime ( const Real t )
{
    Int n = static_cast< Int > (floor (0.5 + M_timePeriod / (2 * M_timeInterval) ) ); // 2n time instances per period

    // Fourier interpolate the values between two time instants

    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {
        BCDataInterpolator_point c1, c2;

        c1.x = 0;
        c1.y = 0;
        c1.z = 0;
        c2.x = 0;
        c2.y = 0;
        c2.z = 0;

        M_dataValues[i].x = 0;
        M_dataValues[i].y = 0;
        M_dataValues[i].z = 0;

        for (Int j = (M_filteringLevel - n); j <= (n - M_filteringLevel); j++)
        {
            for (Int k = 0; k <= 2 * n; k++)
            {
                Real tk = M_timePeriod / (2 * n + 1) * k;
                Int index_ik = indexInTime (i, k);

                c1.x += std::cos (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].x;
                c1.y += std::cos (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].y;
                c1.z += std::cos (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].z;

                c2.x += std::sin (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].x;
                c2.y += std::sin (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].y;
                c2.z += std::sin (-j * tk * 2 * M_PI / M_timePeriod) * M_dataValues_timeSamples[index_ik].z;
            }

            c1.x = c1.x / (2 * n + 1);
            c1.y = c1.y / (2 * n + 1);
            c1.z = c1.z / (2 * n + 1);
            c2.x = c2.x / (2 * n + 1);
            c2.y = c2.y / (2 * n + 1);
            c2.z = c2.z / (2 * n + 1);

            M_dataValues[i].x += c1.x * std::cos (j * t * 2 * M_PI / M_timePeriod) - c2.x * std::sin (j * t * 2 * M_PI / M_timePeriod);
            M_dataValues[i].y += c1.y * std::cos (j * t * 2 * M_PI / M_timePeriod) - c2.y * std::sin (j * t * 2 * M_PI / M_timePeriod);
            M_dataValues[i].z += c1.z * std::cos (j * t * 2 * M_PI / M_timePeriod) - c2.z * std::sin (j * t * 2 * M_PI / M_timePeriod);
        }
    }

    formRBFvectors();

    // Afterwards must solve again the interpolation weights because RHS changed
    M_flagInterpolated = false;
    M_lastInterpolatedAtTime = t;
}

} // namespace LifeV
