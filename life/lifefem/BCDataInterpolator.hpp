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

 @date 28-03-2011

 *///@HEADER

#ifndef BCDATAINTERPOLATOR_H
#define BCDATAINTERPOLATOR 1

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseVector.h"
#include <life/lifecore/LifeV.hpp>
#include <life/lifecore/Displayer.hpp>
#include <life/lifefem/BCBase.hpp>
#include <life/lifefem/BCFunction.hpp>
#include <life/lifefem/BCManage.hpp>

namespace LifeV {

//! BCDataInterpolator - class for interpolating boundary data
/*!
 @author Toni Lassila

 Implements Radial Basis Function (RBF) interpolation of pointwise scalar or vectorial functions defined
 on a set of scattered interpolation points. Currently implements thin-plate splines and multiquadrics.
 For mathematical details see M.D. Buhmann. Radial basis functions: theory and implementations, Cambridge
 University Press, 2004.

 Inherits @c BCFunctionBase to facilitate use of interpolated data as boundary condition.

 If the interpolated data depends on time, the user must pass the data values at 2n specific time instances
 uniformly sampled at interval M_timeInterval and with period M_timePeriod. The values of the data between
 these time instances is interpolated using Fourier interpolation, i.e. the interpolant is a trigonometric
 polynomial of order 2n and periodic with period M_timePeriod.

 The format of the data file passed to readData() should be as follows:\
\
     # header line before dimension definition\
     nof_data_sites nof_data_dimensions t_interval t_period\
     # header line before control point definitions\
     data_site_1_x_coord data_site_1_y_coord data_site_1_z_coord\
     ...\
     data_site_n_x_coord data_site_n_y_coord data_site_n_z_coord\
     # header line before data definitions\
     data_value_1_x_coord data_value_1_y_coord data_value_1_z_coord\
     ...\
     data_value_n_x_coord data_value_n_y_coord data_value_n_z_coord

 The variable nof_data_dimensions has to equal 1 or 3, depending on whether scalar or vectorial data
 is being interpolated. The variable nof_data_sites has to equal the number of rows passed in
 both the section involving the data_sites and the data values.

 Warning: in the current implementation the data sites are assumed fixed in time and they do not move
 with the mesh. Thus they should only be used in a Lagrangian frame of reference, i.e. with structural
 BC's.

 */

class BCDataInterpolator :
                           public BCFunctionBase
{
public:

    //! @name Public Types
    //@{

    typedef Epetra_SerialDenseSolver        solver_Type;
    typedef Epetra_SerialDenseMatrix        matrix_Type;
    typedef Epetra_SerialDenseVector        vector_Type;
    typedef boost::shared_ptr<matrix_Type>  matrix_ptrType;
    typedef boost::shared_ptr<vector_Type>  vector_ptrType;


    /*! @enum BCInterpolation_Type
     Interpolation methods: RBF_ThinPlateSpline
                            RBF_MultiQuadric
     */
    enum BCInterpolationMethod
    {
        RBF_ThinPlateSpline, /*!< Thin plate splines */
        RBF_MultiQuadric,    /*!< Multiquadrics */
    };

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Constructors for an data interpolator
    /*!
     \param comm  the Epetra_Comm to be used for communication
     \param localMap use localMap instead of M_FESpace.map()
     */
    BCDataInterpolator( BCInterpolationMethod interpolationMethod );

    //! Destructor
    ~BCDataInterpolator();

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
     @param bdFunctionDirectional The BCFunctionDirectional object
     @return Reference to a new BCFunctionDirectional object which is a copy of bcFunctionDirectional
     */
    BCFunctionDirectional&
    operator=( const BCDataInterpolator& bcDataInterpolator );

    //@}


    //! @name Methods
    //@{

    //! Evaluate the interpolated function
    /*!
     \param t Time
     \param x Coordinate
     \param y Coordinate
     \param z Coordinate
     \param component The component of the vector function
     \return The selected component of the vector function evaluated in (t,x,y,z)
     */
    Real vectFct( const Real& t,
                  const Real& x,
                  const Real& y,
                  const Real& z,
                  const ID& component );

    //! Display the content of the variables
    /*!
     @param verbose The verbosity (default: true)
     @param out The ostream output (default: std::cout)
     */
    void showMe( bool verbose = false,
                 std::ostream & out = std::cout ) const;

    //! Read control points and data from a file
    /*!

     */
    void readData(const char *fileName);

    //! Exports the interpolation matrix (for debugging purposes)
    /*!

     */
    void exportInterpolationMatrix();

    //@}

    //! @name Set Methods
    //@{



    //! @name Get Methods
    //@{

    //! Returns number of control points
    UInt& nofControlPoints();

    //@}

private:

    struct BCInterpolation_VectorialData
    {
        Real x;
        Real y;
        Real z;
    };

    struct BCInterpolation_ScalarData
    {
        Real s;
    };

    matrix_Type M_interpolationMatrix;

    vector_Type M_rhs_x, M_rhs_y, M_rhs_z;
    vector_Type M_coeffs_x, M_coeffs_y, M_coeffs_z;

    solver_Type M_denseSolver;

    BCInterpolationMethod M_interpolationMethod;

    BCInterpolation_VectorialData* M_dataSites;
    BCInterpolation_VectorialData* M_dataValues;
    BCInterpolation_VectorialData* M_dataValues_timeSamples;

    UInt M_nofControlPoints;

    Real M_lastInterpolatedAtTime;
    Real M_timePeriod;
    Real M_timeInterval;

    bool M_flagInterpolated;
    bool M_verbose;


    void formRBFMatrix();
    void solveInterpolationSystem();
    BCInterpolation_VectorialData interpolateVectorialFunction( const Real& t,
                                                                const Real& x,
                                                                const Real& y,
                                                                const Real& z );
    Real evaluateRBF( const BCInterpolation_VectorialData point1,
                      const BCInterpolation_VectorialData point2 );
    bool needSideConstraints();
    void formRBFvectors();
    void interpolateDataValuesInTime( const Real t );

    Int getIndexInTime(Int dataSite, Int timeInstant) const
    {
        return M_nofControlPoints * timeInstant + dataSite;
    }

};


} // namespace LifeV

#endif /* BCDATAINTERPOLATOR_H */
