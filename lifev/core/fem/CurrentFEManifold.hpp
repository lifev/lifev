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
    @brief  A class for a finite element living on a manifold

    @author Jean-Frederic Gerbeau
            Vincent Martin
    @date 00-09-2002

    @refactoring Luca Bertagna <lbertag@emory.edu>
    @date Sept 2012

    @contributor Luca Bertagna <lbertag@emory.edu>
    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Luca Bertagna <lbertag@emory.edu>
 */

#ifndef CURRENT_FE_MANIFOLD_HPP
#define CURRENT_FE_MANIFOLD_HPP 1

#include <boost/multi_array.hpp>
#include <lifev/core/fem/CurrentFE.hpp>

namespace LifeV
{

/* A few more update flags specific for the boundary elements.
   For more info about how update flags work, see CurrentFE.hpp
*/

const flag_Type UPDATE_ONLY_TANGENTS (16384);
const flag_Type UPDATE_ONLY_NORMALS (32768);
const flag_Type UPDATE_ONLY_METRIC (65536);
const flag_Type UPDATE_ONLY_DET_METRIC (131072);
const flag_Type UPDATE_ONLY_W_ROOT_DET_METRIC (262144);
const flag_Type UPDATE_ONLY_INV_METRIC (524288);

const flag_Type UPDATE_TANGENTS (UPDATE_ONLY_TANGENTS
                                 | UPDATE_ONLY_CELL_NODES);
const flag_Type UPDATE_NORMALS (UPDATE_ONLY_NORMALS
                                | UPDATE_ONLY_TANGENTS
                                | UPDATE_ONLY_CELL_NODES);
const flag_Type UPDATE_METRIC (UPDATE_ONLY_METRIC
                               | UPDATE_ONLY_TANGENTS
                               | UPDATE_ONLY_CELL_NODES);
const flag_Type UPDATE_INV_METRIC (UPDATE_ONLY_INV_METRIC
                                   | UPDATE_METRIC
                                   | UPDATE_ONLY_DET_METRIC);
const flag_Type UPDATE_W_ROOT_DET_METRIC (UPDATE_ONLY_W_ROOT_DET_METRIC
                                          | UPDATE_METRIC
                                          | UPDATE_ONLY_DET_METRIC);

/*!
  \class CurrentFEManifold
  \brief A class for a finite element on a manifold

  This class inherits from CurrentFE and adds some functionality related
  to the manifold, like normals, tangents, etc.
*/

class CurrentFEManifold : public CurrentFE
{

public:
    //! @name Constructor(s) & Destructor
    //@{

    //! Constructor without quadrature specification
    /*!
      If you use this constructor, you have to use the method CurrentFE::setQuadRule before
      using any update method!

      @param refFE Reference finite element used
      @param geoMap Geometric mapping used
     */
    CurrentFEManifold (const ReferenceFE& refFE, const GeometricMap& geoMap);

    //! Full constructor with default quadrature
    /*!
      @param refFE Reference finite element used
      @param geoMap Geometric mapping used
      @param qr Quadrature rule for the computations
     */
    CurrentFEManifold (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr);

    //! Constructor used to create a static reference FE for hybrid elements
    CurrentFEManifold (const ReferenceFE& refFE, const GeometricMap& geoMap, const QuadratureRule& qr,
                       const Real* refCoor, UInt currentId, Real invArea = 1. );

    //! Copy constructor
    CurrentFEManifold (const CurrentFEManifold& feManifold);

    // Virtual destructor (does not do anything)
    virtual ~CurrentFEManifold () {}
    //@}

    //! @name Methods
    //@{

    // The update method which takes the geometric element is inherited from the base class
    using CurrentFE::update;

    //! Update method using only point coordinates.
    /*!
        Overrides the method of base class CurrentFE. Actually, it calls the method of the base
        class and then performs some extra updates if upFlag contains also some boundary updates
        @param pts The coordinates of the points defining the current element
        @param upFlag   The update flag which determines what kind of update has to be done (see top of the file)
     */
    virtual void update (const std::vector<std::vector<Real> >& pts, flag_Type upFlag);

    //! Maps a point from the reference element to the current element
    /*!
      return the coordinate (x,y,z)= F(xi,eta), (F: geo mappping)
      where (xi,eta) are the coordinates in the reference element
      and (x,y,z) are the coor in the current element
      (if the code is compiled in 2D mode then z=0 and eta is disgarded)
    */
    void coorMap (Real& x, Real& y, Real& z, Real xi, Real eta) const;

    //! Compute the measure of the current element
    /*!
       Overrides the corresponding method in the mother class
     */
    Real measure () const;

    //! Compute the integral of a function f over the current element
    /*!
        @param f The function to be integrated. It must have the signature
                 Real f(Real x, Real y, Real z)
     */
    template <typename FunctorType>
    Real integral (const FunctorType& f) const;

    //! Compute the integral of the normal component of a vector function f over the current element
    /*!
        @param f The function to be integrated. It must have the signature
                 void f(Real x, Real y, Real z, Real* result)
                 where result is an array of size equal to the dimension of the ambient space.
                 The function should NOT create the array, but instead should assume that
                 the array is already created. The dimension of the ambient spase should be
                 equal to M_nbLocalCoor+1
     */
    template <typename FunctorType>
    Real normalIntegral (const FunctorType& f) const;
    //@}


    //! @name Get Methods
    //@{

    //! Values of the geometric basis functions on quadrature points (different from M_phi(iGeoNode,quadNode) only for Hybrid elements)
    Real phiGeo (UInt iGeoNode, UInt quadNode) const
    {
        return M_phiGeo[iGeoNode][quadNode];
    }

    //! Values of the tangents on the quadrature points
    Real tangent (UInt tangent, UInt coordinate, UInt quadNode) const
    {
        ASSERT (M_tangentsUpdated, "Tangents are not updated!\n");
        return M_tangents[tangent][coordinate][quadNode];
    }

    //! Values of the normal on the quadrature points
    Real normal (UInt coordinate, UInt quadNode) const
    {
        ASSERT (M_normalUpdated, "Normals are not updated!\n");
        return M_normal[coordinate][quadNode];
    }

    //! Metric tensor on the quadrature points
    Real metric (UInt iCoor, UInt jCoor, UInt quadNode) const
    {
        ASSERT (M_metricUpdated, "Metric is not updated!\n");
        return M_metric[iCoor][jCoor][quadNode];
    }

    //! Determinant of the metric tensor on the quadrature points
    Real detMetric (UInt quadNode)
    {
        ASSERT (M_detMetricUpdated, "Determinant of the metric is not updated!\n");
        return M_detMetric[quadNode];
    }

    //! Inverse of the metric tensor on the quadrature points
    Real inverseMetric (UInt iCoor, UInt jCoor, UInt quadNode) const
    {
        ASSERT (M_inverseMetricUpdated, "Inverse metric is not updated!\n");
        return M_inverseMetric[iCoor][jCoor][quadNode];
    }

    //! Square root of the determinant of the metric times the weight on the quadrature points
    Real wRootDetMetric (UInt quadNode) const
    {
        ASSERT (M_wRootDetMetricUpdated, "Weighted metric determinant is not updated!\n");
        return M_wRootDetMetric[quadNode];
    }

protected:

    //! Computes the tangent vectors in the quadrature nodes
    void computeTangents ();

    //! Computes the normal vectors in the quadrature nodes
    void computeNormal ();

    //! Computes the metric in the quadrature nodes
    void computeMetric ();

    //! Computes the determinant of the metric tensor in the quadrature nodes
    void computeDetMetric ();

    //! Computes the inverse of the metric tensor in the quadrature nodes
    void computeInverseMetric ();

    //! Computes the square root of the determinant of the metric times the weight in the quadrature nodes
    void computeWRootDetMetric ();

    //! Geometric map in the geometric nodes (different from M_phi only for Hybrid elements)
    boost::multi_array<Real, 2> M_phiGeo;

    //! Normal and tangent vectors, metric tensor and weighted measure in the quadrature nodes
    boost::multi_array<Real, 3> M_tangents;
    boost::multi_array<Real, 2> M_normal;
    boost::multi_array<Real, 3> M_metric;
    boost::multi_array<Real, 1> M_detMetric;
    boost::multi_array<Real, 3> M_inverseMetric;
    boost::multi_array<Real, 1> M_wRootDetMetric;

    //! Check variables
    bool M_tangentsUpdated;
    bool M_normalUpdated;
    bool M_metricUpdated;
    bool M_detMetricUpdated;
    bool M_inverseMetricUpdated;
    bool M_wRootDetMetricUpdated;
};

template <typename FunctorType>
Real CurrentFEManifold::integral (const FunctorType& f) const
{
    ASSERT (M_quadNodesUpdated && M_wRootDetMetricUpdated, "Error! Quadrature nodes and Jacobian Determinant have not been updated yet.\n");
    Real result = 0.0;
    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        result += f (M_quadNodes[iq][0], M_quadNodes[iq][1], M_quadNodes[iq][2]) * M_wRootDetMetric[iq];
    }

    return result;
}

template <typename FunctorType>
Real CurrentFEManifold::normalIntegral (const FunctorType& f) const
{
    ASSERT (M_quadNodesUpdated && M_normalUpdated && M_wRootDetMetricUpdated, "Error! Normal and Jacobian Determinant have not been updated yet.\n");
    Real result = 0.0;

    Real tmp;
    Real* returnValues = new Real[M_nbLocalCoor + 1];
    for (UInt iq (0); iq < M_nbQuadPt; ++iq)
    {
        tmp = 0;
        f (M_quadNodes[iq][0], M_quadNodes[iq][1], M_quadNodes[iq][2], returnValues);
        for (UInt iCoor (0); iCoor <= M_nbLocalCoor; ++iCoor)
        {
            tmp += returnValues[iCoor] * M_normal[iCoor][iq];
        }
        result += tmp * M_wRootDetMetric[iq];
    }
    delete[] returnValues;
    return result;
}

} // Namespace LifeV

#endif // CURRENT_FE_MANIFOLD_HPP
