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
 @brief File containing BCDataInterpolator class for interpolating boundary data

 @author Toni Lassila <toni.lassila@epfl.ch>
 @maintainer Toni Lassila <toni.lassila@epfl.ch

 @date 23-05-2011

 Implements Radial Basis Function (RBF) interpolation of pointwise scalar or vectorial functions defined
 on a set of scattered interpolation points. Currently implements thin-plate splines and multiquadrics.
 For mathematical details see M.D. Buhmann. Radial basis functions: theory and implementations, Cambridge
 University Press, 2004.

 Inherits @c BCFunctionBase to facilitate use of interpolated data as boundary condition.

 The format of the data file passed to readData() should be as follows:\
\
     # header line before dimension definition\
     nof_data_sites nof_data_dimensions\
     # header line before control point definitions\
     data_site_1_x_coord data_site_1_y_coord data_site_1_z_coord\
     ...\
     data_site_n_x_coord data_site_n_y_coord data_site_n_z_coord\
     # header line before data definitions\
     data_value_1_x_coord data_value_1_y_coord data_value_1_z_coord\
     ...\
     data_value_n_x_coord data_value_n_y_coord data_value_n_z_coord\

 The variable nof_data_dimensions has to equal 1 or 3, depending on whether scalar or vectorial data
 is being interpolated. The variable nof_data_sites has to equal the number of rows passed in
 both the section involving the data_sites and the data values.

 Warning: in the current implementation the data sites are assumed fixed in time and they do not move
 with the mesh. Thus they should only be used in a Lagrangian frame of reference, i.e. with structural
 BC's.

 *///@HEADER

#include <sstream>
#include <stdexcept>
#include <boost/bind.hpp>
#include <life/lifecore/LifeV.hpp>
#include <life/lifefem/BCBase.hpp>
#include <life/lifefem/BCFunction.hpp>
#include <life/lifefem/BCDataInterpolator.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
BCDataInterpolator::BCDataInterpolator( BCInterpolationMethod interpolationMethod ) :
                M_interpolationMatrix(), M_rhs_x(), M_rhs_y(), M_rhs_z(), M_coeffs_x(), M_coeffs_y(),
                    M_coeffs_z(), M_denseSolver(), M_interpolationMethod ( interpolationMethod ), M_flagInterpolated( false )
{

    M_userDefinedFunction = boost::bind(&BCDataInterpolator::vectFct, this, _1, _2, _3, _4, _5);

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

Real BCDataInterpolator::vectFct( const Real& t,
                                  const Real& x,
                                  const Real& y,
                                  const Real& z,
                                  const ID& component )
{

    BCInterpolation_VectorialData vectEval;

    vectEval = interpolateVectorialFunction( x,
                                             y,
                                             z );

    switch ( component )
    {
        case 0:
            return vectEval.x;
        case 1:
            return vectEval.y;
        case 2:
            return vectEval.z;
        default:
            std::ostringstream __ex;
            __ex << "Invalid component: " << component << std::endl;
            throw std::invalid_argument( __ex.str() );
    }
}

void BCDataInterpolator::showMe( bool verbose,
                                 std::ostream& out ) const
{

    out << " Boundary Conditions Data Interpolator ====>" << std::endl;

}

void BCDataInterpolator::readData(const char *fileName)
{

    std::ifstream fin;
    UInt RDIM;

    fin.open(fileName);

    if ( !fin.fail() )
    {
        fin.ignore(80,'\n');

        fin >> M_nofControlPoints;
        fin >> RDIM;

        if ((RDIM != 1) && (RDIM != 3))
        {
            std::ostringstream __ex;
                        __ex << "Interpolated data must be either scalar or a 3-vector: " << std::endl;
                        throw std::invalid_argument( __ex.str() );
        }

        fin.ignore(80,'\n');
        fin.ignore(80,'\n');

        M_dataSites = new BCInterpolation_VectorialData[M_nofControlPoints];
        M_dataValues = new BCInterpolation_VectorialData[M_nofControlPoints];

        for (UInt i=0; i < M_nofControlPoints; i++)
        {
            fin >> M_dataSites[i].x;
            fin >> M_dataSites[i].y;
            fin >> M_dataSites[i].z;
        }

        fin.ignore(80,'\n');
        fin.ignore(80,'\n');

        for (UInt i=0; i < M_nofControlPoints; i++)
        {
            fin >> M_dataValues[i].x;
            fin >> M_dataValues[i].y;
            fin >> M_dataValues[i].z;
        }

        fin.close();

        // We can form the interpolation matrix now

        formRBFMatrix();
    }


}

void BCDataInterpolator::exportInterpolationMatrix()
{

    // TODO: Write matrix exporting routine

}

// ===================================================
// Set Methods
// ===================================================


// ===================================================
// Get Methods
// ===================================================
UInt&
BCDataInterpolator::nofControlPoints()
{

    return M_nofControlPoints;

}

// ===================================================
// Private Methods
// ===================================================

void BCDataInterpolator::formRBFMatrix()
{

    Int nofElementsPerRow = ( needSideConstraints() ? 2 * M_nofControlPoints + NDIM : M_nofControlPoints );

    M_interpolationMatrix.Shape( static_cast< int > ( nofElementsPerRow ),
                       static_cast< int > ( nofElementsPerRow ) );

    M_rhs_x.Size( static_cast< int > ( nofElementsPerRow ) );
    M_rhs_y.Size( static_cast< int > ( nofElementsPerRow ) );
    M_rhs_z.Size( static_cast< int > ( nofElementsPerRow ) );

    M_coeffs_x.Size( static_cast< int > ( nofElementsPerRow ) );
    M_coeffs_y.Size( static_cast< int > ( nofElementsPerRow ) );
    M_coeffs_z.Size( static_cast< int > ( nofElementsPerRow ) );

    Real rbfEval;

    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {

        for ( UInt j = 0; j <= i; j++ )
        {
            // Symmetric interpolation matrix [A]_ij = phi(|x_i - x_j|)

            rbfEval = evaluateRBF( M_dataSites[i],
                                   M_dataSites[j] );
            M_interpolationMatrix( static_cast< int >(i),
                         static_cast< int >(j) ) = static_cast< double >(rbfEval);
            M_interpolationMatrix( static_cast< int >(j),
                         static_cast< int >(i) ) = static_cast< double >(rbfEval);

        }

        // Form RHS for the values at the RBF sites

        M_rhs_x( i ) = M_dataValues[i].x;
        M_rhs_y( i ) = M_dataValues[i].y;
        M_rhs_z( i ) = M_dataValues[i].z;
    }

    if ( needSideConstraints() )
    {

        // Add global linear polynomial to the interpolation basis

        for ( UInt i = 0; i < M_nofControlPoints; i++ )
        {

            // Constant term of the polynomial
            M_interpolationMatrix( i,
                         M_nofControlPoints + i ) = 1;
            M_interpolationMatrix( M_nofControlPoints + i,
                         i ) = 1;

            // First order term of the polynomial
            M_interpolationMatrix( i,
                         2 * M_nofControlPoints ) = M_dataSites[i].x;
            M_interpolationMatrix( i,
                         2 * M_nofControlPoints + 1 ) = M_dataSites[i].y;
            M_interpolationMatrix( i,
                         2 * M_nofControlPoints + 2 ) = M_dataSites[i].z;

            M_interpolationMatrix( 2 * M_nofControlPoints,
                         i ) = M_dataSites[i].x;
            M_interpolationMatrix( 2 * M_nofControlPoints + 1,
                         i ) = M_dataSites[i].y;
            M_interpolationMatrix( 2 * M_nofControlPoints + 2,
                         i ) = M_dataSites[i].z;
        }
    }

    M_denseSolver.SetMatrix( M_interpolationMatrix );

}

void BCDataInterpolator::solveInterpolationSystem()
{

    #ifdef DEBUG
    // Check that condition number is "reasonable" (< 1e6), if not print a warning
    Real reciprocalCond = M_denseSolver.RCOND();

    if (reciprocalCond < 1e-6)
    {
        Debug( 5000 ) << "Radial basis interpolation matrix is ill-conditioned.\n";
        Debug( 5000 ) << "Estimated reciprocal of condition number: " << reciprocalCond << "\n";
    }
    #endif

    // Solve the interpolation system for each component using LU-factorization

    M_denseSolver.SetVectors( M_coeffs_x,
                               M_rhs_x );
    M_denseSolver.Solve();

    M_denseSolver.SetVectors( M_coeffs_x,
                               M_rhs_y );
    M_denseSolver.Solve();

    M_denseSolver.SetVectors( M_coeffs_z,
                               M_rhs_z );
    M_denseSolver.Solve();

}

BCDataInterpolator::BCInterpolation_VectorialData BCDataInterpolator::interpolateVectorialFunction( const Real& x,
                                                                                                    const Real& y,
                                                                                                    const Real& z )
{

    BCInterpolation_VectorialData rbfEval, evalPoint;
    Real rbfShape;

    rbfEval.x = 0;
    rbfEval.y = 0;
    rbfEval.z = 0;
    evalPoint.x = x;
    evalPoint.y = y;
    evalPoint.z = z;

    // If not yet solved the interpolation system, do so
    if ( !M_flagInterpolated )
    {
        solveInterpolationSystem();
        M_flagInterpolated = true;
    }

    // For each data site...
    for ( UInt i = 0; i < M_nofControlPoints; i++ )
    {

        rbfShape = evaluateRBF( evalPoint,
                                M_dataSites[i] );

        // Interpolated function value at point (x,y,z)
        rbfEval.x += M_coeffs_x( i ) * rbfShape;
        rbfEval.y += M_coeffs_y( i ) * rbfShape;
        rbfEval.z += M_coeffs_z( i ) * rbfShape;

    } // ...for each data site


    if (needSideConstraints())
    {
        // For each data site...
        for ( UInt i = 0; i < M_nofControlPoints; i++ )
        {
            // Add the constant term of the polynomial...
            rbfEval.x += M_coeffs_x( i + M_nofControlPoints );
            rbfEval.y += M_coeffs_y( i + M_nofControlPoints );
            rbfEval.z += M_coeffs_z( i + M_nofControlPoints );

            // ...and the linear term of the polynomial
            rbfEval.x += M_coeffs_x( 2 * M_nofControlPoints )    * M_dataSites[i].x;
            rbfEval.x += M_coeffs_x( 2 * M_nofControlPoints + 1) * M_dataSites[i].y;
            rbfEval.x += M_coeffs_x( 2 * M_nofControlPoints + 2) * M_dataSites[i].z;

            rbfEval.y += M_coeffs_y( 2 * M_nofControlPoints )    * M_dataSites[i].x;
            rbfEval.y += M_coeffs_y( 2 * M_nofControlPoints + 1) * M_dataSites[i].y;
            rbfEval.y += M_coeffs_y( 2 * M_nofControlPoints + 2) * M_dataSites[i].z;

            rbfEval.z += M_coeffs_z( 2 * M_nofControlPoints )    * M_dataSites[i].x;
            rbfEval.z += M_coeffs_z( 2 * M_nofControlPoints + 1) * M_dataSites[i].y;
            rbfEval.z += M_coeffs_z( 2 * M_nofControlPoints + 2) * M_dataSites[i].z;
        } // ...for each data site

    }

    return rbfEval;

}

Real BCDataInterpolator::evaluateRBF( const BCInterpolation_VectorialData point1,
                                      const BCInterpolation_VectorialData point2 )
{

    Real h = sqrt( pow( point1.x - point2.x,
                        2 ) + pow( point1.y - point2.y,
                                   2 ) + pow( point1.z - point2.z,
                                              2 ) );

    switch ( M_interpolationMethod )
    {

        case RBF_ThinPlateSpline:
            if ( h < 1e-3 )
            {
                return 0.0;
            }
            else
            {
                return pow( h,
                            2 ) * log( h );
            }

        case RBF_MultiQuadric:
            return sqrt( pow( h,
                              2 ) + 1 );

        default:
            return 0.0;
    }

}

bool BCDataInterpolator::needSideConstraints()
{
    switch ( M_interpolationMethod )
    {

        case RBF_MultiQuadric:
            return false;

        default:
            return true;
    }
}

} // namespace LifeV
