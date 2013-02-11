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

#ifndef BCDATAINTERPOLATOR_H
#define BCDATAINTERPOLATOR 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/shared_array.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_SerialDenseVector.h"

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/fem/BCBase.hpp>
#include <lifev/core/fem/BCFunction.hpp>
#include <lifev/core/fem/BCManage.hpp>

namespace LifeV
{

//! BCDataInterpolator - Class for interpolating boundary functions from scattered data
/*!
  @author Toni Lassila
  @see Radial basis function interpolation \cite Buhmann2004

  Implements Radial Basis Function (RBF) interpolation of pointwise scalar or vectorial functions defined
  on a set of scattered interpolation points. Currently implements a variety of different basis functions
  for interpolation. Temporal interpolation done with trigonometric polynomials.

  Inherits @c BCFunctionBase to facilitate use of interpolated data as boundary condition.

  If the interpolated data depends on time, the user must pass the data values at 2n specific time instances,
  uniformly sampled at interval M_timeInterval and with period M_timePeriod. The values of the data between
  these time instances is interpolated using Fourier interpolation, i.e. the interpolant is a trigonometric
  polynomial of order 2n and periodic with period M_timePeriod.

  The format of the data file passed to readData() is the following: <BR>
 <BR>
     # HEADER LINE FOR PARAMETERS <BR>
     nof_data_sites nof_data_dimensions t_interval t_period filtering_level <BR>
     # HEADER LINE FOR DATA SITES <BR>
     data_site_1_x_coord data_site_1_y_coord data_site_1_z_coord <BR>
     ... <BR>
     data_site_n_x_coord data_site_n_y_coord data_site_n_z_coord <BR>
     # HEADER LINE FOR DATA VALUES <BR>
     data_value_1_x_coord data_value_1_y_coord data_value_1_z_coord <BR>
     ... <BR>
     data_value_n_x_coord data_value_n_y_coord data_value_n_z_coord <BR>
     # HEADER LINE FOR DATA VALUES <BR>
     data_value_1_x_coord data_value_1_y_coord data_value_1_z_coord <BR>
     ... <BR>
     data_value_n_x_coord data_value_n_y_coord data_value_n_z_coord <BR>
  <BR>
  The variable nof_data_dimensions has to equal 1 or 3, depending on whether scalar or vectorial data
  is being interpolated. The variable nof_data_sites has to equal the number of rows passed in
  both the section involving the data_sites and the data values. The data value section has to be
  repeated t_period / t_interval times. The value filtering_level >= 0 is used to drop the most
  oscillatory terms in the trigonometric polynomial, and should be an integer.

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
    typedef boost::shared_ptr<matrix_Type>  matrixPtr_Type;
    typedef boost::shared_ptr<vector_Type>  vectorPtr_Type;


    /*! @enum BCInterpolation_Type
     Interpolation methods: RBF_ThinPlateSpline
                            RBF_MultiQuadric
                            RBF_Cubic
                            RBF_Gaussian
                            RBF_InverseMultiQuadric
     */
    enum BCInterpolationMethod
    {
        RBF_InverseMultiQuadric, /*!< Inverse multiquadrics */
        RBF_Gaussian,            /*!< Gaussians */
        RBF_ThinPlateSpline,     /*!< Thin plate splines */
        RBF_MultiQuadric,        /*!< Multiquadrics */
        RBF_Cubic               /*!< Cubics */
    };

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructors for an data interpolator
    /*!

     */
    explicit BCDataInterpolator();

    //! Destructor
    virtual ~BCDataInterpolator();

    //@}


    //! @name Operators
    //@{

    //! Assignment operator
    /*!
     @param bdFunctionDirectional The BCFunctionDirectional object
     @return Reference to a new BCFunctionDirectional object which is a copy of bcFunctionDirectional
     */
    //BCFunctionDirectional& operator=( const BCDataInterpolator& bcDataInterpolator );

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
    Real interpolatedDataFunction ( const Real& t,
                                    const Real& x,
                                    const Real& y,
                                    const Real& z,
                                    const ID& component );

    //! Display the content of the variables
    /*!
     @param verbose The verbosity (default: false)
     @param out The ostream output (default: std::cout)
     */
    void showMe ( bool verbose = false,
                  std::ostream& out = std::cout ) const;

    //! Read control points and data from a file
    /*!
     @param filename The filename for the data sites and data values
     */
    void readData ( const std::string& fileName );

    //! Export the interpolation matrix for debugging purposes
    /*!

     */
    void exportInterpolationMatrix() const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the interpolation method
    /*!
     @param bcInterpolationMethod The interpolation method
     */
    void setInterpolationMethod ( const BCInterpolationMethod& bcInterpolationMethod )
    {
        M_interpolationMethod = bcInterpolationMethod;
    }

    //! Set the filtering level
    /*!
     @param filteringLevel The filtering level
     */
    void setFilteringLevel ( const Int& filteringLevel )
    {
        M_filteringLevel = filteringLevel;
    }


    //! @name Get Methods
    //@{

    //! Returns the number of control points
    /*!
     * @return The number of control points
     */
    const UInt& nofControlPoints() const
    {
        return M_nofControlPoints;
    }

    //@}

private:

    struct BCDataInterpolator_point
    {
        Real x;
        Real y;
        Real z;
    };

    //! @name Private Methods
    //@{

    void formRBFMatrix();
    void solveInterpolationSystem();
    BCDataInterpolator_point interpolateVectorialFunction ( const Real& t,
                                                            const Real& x,
                                                            const Real& y,
                                                            const Real& z );
    Real evaluateRBF ( const BCDataInterpolator_point point1,
                       const BCDataInterpolator_point point2 );
    bool needSideConstraints() const;
    void formRBFvectors();
    void interpolateDataValuesInTime ( const Real t );

    Int indexInTime (const Int dataSite, const Int timeInstant) const
    {
        return M_nofControlPoints * timeInstant + dataSite;
    }

    //@}

    matrix_Type M_interpolationMatrix;

    vector_Type M_rhs_x, M_rhs_y, M_rhs_z;
    vector_Type M_coeffs_x, M_coeffs_y, M_coeffs_z;

    solver_Type M_denseSolver;

    BCInterpolationMethod M_interpolationMethod;

    boost::shared_array<BCDataInterpolator_point> M_dataSites;
    boost::shared_array<BCDataInterpolator_point> M_dataValues;
    boost::shared_array<BCDataInterpolator_point> M_dataValues_timeSamples;

    UInt M_nofControlPoints;

    Int M_filteringLevel;

    Real M_lastInterpolatedAtTime;
    Real M_timePeriod;
    Real M_timeInterval;

    bool M_flagInterpolated;
    bool M_verbose;
};

} // namespace LifeV

#endif /* BCDATAINTERPOLATOR_H */
