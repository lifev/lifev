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
    @brief Base class for RefFE and GeoMap

    @author Jean-Frederic Gerbeau
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 00-04-2002

    @contributor
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef REFELE_H
#define REFELE_H 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <boost/numeric/ublas/vector.hpp>

#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <life/lifecore/life.hpp>

#include <life/lifemesh/basisElSh.hpp>


namespace LifeV
{
    typedef boost::numeric::ublas::vector<Real> GeoVector;

//! RefEle - The basis class for the geometric mapping and the reference finite elements.
/*!
  @author J.-F. Gerbeau
  @date 04/2002

  Implemented orginially by J.-F. Gerbeau (04/2002) but totally modified by S.Quinodoz (samuel.quinodoz@epfl.ch , 04/2010)

  This class contains all the basis functions, their derivatives and the reference coordinates. It is the basis class for the geometric map (LifeV::GeoMap) and the reference finite element (LifeV::RefFE).

  \todo Add Volume
  \todo Add the minimal dimension and checks
  \todo Incorporate vectorial FEs
  \todo Think dimensionless
  \todo change M_refCoor

*/
class RefEle
{

public:

    //! @name Public Types
    //@{

    // Some typedefs for functions
    typedef Real ( * function_Type ) ( const GeoVector& );

    //@}

    //! @name Constructor & Destructor
    //@{

    //! Full constructor
    /*!
      @param name Name of the reference element
      @param shape Shape related to this reference element
      @param nbDof Number of degrees of freedom
      @param nbCoor Number of local coordinates
      @param phi Array of the basis functions
      @param dPhi Array of the derivatives of the basis functions
      @param d2Phi Array of the second derivatives of the basis functions
      @param refCoor Array of the reference coordinates for this reference element
     */
    RefEle( std::string name, ReferenceShapes shape, UInt nbDof, UInt nbCoor, UInt feDim,
            const function_Type* phi, const function_Type* dPhi, const function_Type* d2Phi,
            const function_Type* divPhi, const Real* refCoor);

    //! Destructor
    virtual ~RefEle();

    //@}



    //! @name Methods
    //@{

    //! return the first local coordinate of the i-th node of the reference element
    Real xi( UInt i ) const
    {
        ASSERT_BD( i < M_nbDof )
        return M_refCoor[ 3 * i ];
    }
    //! return the second local coordinate of the i-th node of the reference element
    Real eta( UInt i ) const
    {
        ASSERT_BD( i < M_nbDof )
        return M_refCoor[ 3 * i + 1 ];
    }
    //! return the third local coordinate of the i-th node of the reference element
    Real zeta( UInt i ) const
    {
        ASSERT_BD( i < M_nbDof )
        return M_refCoor[ 3 * i + 2 ];
    }
    //! return the icoor-th local coordinate of the i-th node of the reference element
    Real refCoor( UInt i, UInt icoor ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_nbCoor )
        return M_refCoor[ 3 * i + icoor ];
    }
    //! return the coordinates of the reference element
    std::vector<GeoVector> refCoor() const;

    //! Return the value of the i-th basis function in the point v
    Real phi( UInt i, const GeoVector& v ) const
    {
        ASSERT_BD( i < M_nbDof )
        return M_phi[ i ] ( v );
    }

    //! return the value of the component icoor-th of the i-th basis function on point v.
    Real phi( UInt i, UInt icoor, const GeoVector& v ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_feDim )
        return M_phi[ i * M_feDim + icoor ] ( v );
    }

    //! return the value of the icoor-th derivative of the i-th basis function on point v
    Real dPhi( UInt i, UInt icoor, const GeoVector& v ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_nbCoor )
        return M_dPhi[ i * M_nbCoor + icoor ] ( v );
    }

    //!  return the value of the (icoor,jcoor)-th second derivative of the i-th basis function on point v
    Real d2Phi( UInt i, UInt icoor, UInt jcoor, const GeoVector& v ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_nbCoor && jcoor < M_nbCoor )
        return M_d2Phi[ ( i * M_nbCoor + icoor ) * M_nbCoor + jcoor ] ( v );
    }
    //! return the value of the divergence of the i-th basis function on point v.
    Real divPhi( UInt i, const GeoVector& v ) const
    {
        ASSERT_BD( i < M_nbDof )
        return M_divPhi[ i ] ( v );
    }


    //! Check if the refEle has phi functions
    bool hasPhi() const
    {
        return ( M_phi != static_cast<function_Type*>(NULL) );
    }
    //! Check if the refEle has dPhi functions
    bool hasDPhi() const
    {
        return ( M_dPhi != static_cast<function_Type*>(NULL) );
    }
    //! Check if the refEle has d2Phi functions
    bool hasD2Phi() const
    {
        return ( M_d2Phi != static_cast<function_Type*>(NULL) );
    }
    //! Check if the refEle has divPhi functions
    bool hasDivPhi() const
    {
        return ( M_divPhi != static_cast<function_Type*>(NULL) );
    }

    //! Method for transforming nodal values into FE values
    /*!
      This method can be used to retrieve the FE coefficients corresponding
      to the values given in the nodes (important for the interpolation
      procedures). For lagrangian elements (like P1 or P2),
      this method is just giving back the same values. However, for nodal but
      non-lagrangian finite elements (like P1Bubble), the values returned are
      different from the output.

      For example, using P1Bubble finite element, if one gives as input values
      that are all 1 (the finite element function is constant), this
      function will return 1 for the P1 nodes but 0 for the bubble. Indeed, by
      suming the P1 function, we already get the constant function over the
      whole element and there is no need for the bubble.

      Of course, this method is accessible only for nodal finite elements.
      If one tries to use this method with a non nodal finite element, an
      error will be displayed and the program will stop running.
     */
    virtual std::vector<Real> nodalToFEValues(const std::vector<Real>& /*nodalValues*/) const
    {
        //By default, it is not possible to use it.
        std::cerr << " Trying to access nodal values via nodalToFEValues function. " << std::endl;
        std::cerr << " This FE is not nodal, impossible operation! " << std::endl;
        abort();
    }

    //@}



    //! @name Get Methods
    //@{

    //! Return the name of the reference element.
    const std::string& name() const
    {
        return M_name;
    }

    //! Return the number of degrees of freedom for this reference element
    const UInt& nbDof() const
    {
        return M_nbDof;
    }

    //! Return the number of local coordinates
    const UInt& nbCoor() const
    {
        return M_nbCoor;
    }

    //! Return the dimension of the FE (scalar vs vectorial FE)
    const UInt& feDim() const
    {
        return M_feDim;
    }
/*    const UInt& __attribute__ (( __deprecated__)) FEDim() const
    {
        return feDim();
        }*/


    //! Return the shape of the element
    const ReferenceShapes& shape() const
    {
        return M_shape;
    }

    //@}

    //! @name Old Methods - Avoid using them!
    //@{

    //! return the value of the i-th basis function on point (x,y,z)
    inline Real phi( UInt i, const Real& x, const Real& y, const Real& z ) const
    {
        ASSERT_BD( i < M_nbDof )
        GeoVector v(3);
        v[0]=x;
        v[1]=y;
        v[2]=z;
        return M_phi[ i ] ( v );
    }

    //! return the value of the component icoor-th of the i-th basis function on point (x,y,z).
    inline Real phi( UInt i, UInt icoor, const Real& x, const Real& y, const Real& z ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_feDim )
        GeoVector v(3);
        v[0]=x;
        v[1]=y;
        v[2]=z;
        return M_phi[ i * M_feDim + icoor ] ( v );
    }

    //! return the value of the icoor-th derivative of the i-th basis function on point (x,y,z)
    inline Real dPhi( UInt i, UInt icoor, const Real& x, const Real& y, const Real& z ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_nbCoor )
        GeoVector v(3);
        v[0]=x;
        v[1]=y;
        v[2]=z;
        return M_dPhi[ i * M_nbCoor + icoor ] ( v );
    }
    //!  return the value of the (icoor,jcoor)-th second derivative of the i-th basis function on point (x,y,z)
    inline Real d2Phi( UInt i, UInt icoor, UInt jcoor, const Real& x, const Real& y, const Real& z ) const
    {
        ASSERT_BD( i < M_nbDof && icoor < M_nbCoor && jcoor < M_nbCoor )
        GeoVector v(3);
        v[0]=x;
        v[1]=y;
        v[2]=z;
        return M_d2Phi[ ( i * M_nbCoor + icoor ) * M_nbCoor + jcoor ] ( v );
    }
    //! return the value of the divergence of the i-th basis function on point (x,y,z).
    inline Real divPhi( UInt i, const Real& x, const Real& y, const Real& z ) const
    {
        ASSERT_BD( i < M_nbDof )
        GeoVector v(3);
        v[0]=x;
        v[1]=y;
        v[2]=z;
        return M_divPhi[ i ] ( v );
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No way to use the empty constructor
    RefEle();

    //! No way to use the copy constuctor
    RefEle(const RefEle&);

    //@}


    //! pointer on the basis functions
    const function_Type* M_phi;

    //! pointer on the derivatives of the basis functions
    const function_Type* M_dPhi;

    //! pointer on the second derivatives of the basis functions
    const function_Type* M_d2Phi;

    //! pointer on the divergence of the basis functions
    const function_Type* M_divPhi;

    //! reference coordinates. Order: xi_1,eta_1,zeta_1,xi_2,eta_2,zeta_2,...
    const Real* M_refCoor;



    //! name of the reference element
    const std::string M_name;

    //! geometrical shape of the element
    const ReferenceShapes M_shape;

    //! Total number of degrees of freedom
    const UInt M_nbDof;

    //! Number of local coordinates
    const UInt M_nbCoor;

    //! Number of dimension of the FE (1 for scalar FE, more for vectorial FE)
    const UInt M_feDim;

};






//======================================================================
//
//                            P0  (0D)
//
//======================================================================
/*
                           1
*/

Real fct1_P0_0D( const GeoVector& );
Real derfct1_P0_0D( const GeoVector & );
Real der2fct1_P0_0D( const GeoVector & );

static const Real refcoor_P0_0D[ 3 ] =
{
    1. , 0. , 0.
};

static const RefEle::function_Type fct_P0_0D[ 1 ] =
{
    fct1_P0_0D
};

static const RefEle::function_Type derfct_P0_0D[ 1 ] =
{
    derfct1_P0_0D
};

static const RefEle::function_Type der2fct_P0_0D[ 1 ] =
{
    der2fct1_P0_0D
};

//======================================================================
//
//                            P1  (1D)
//
//======================================================================
/*
                           1-----2
*/
Real fct1_P1_1D( const GeoVector& v );
Real fct2_P1_1D( const GeoVector& v );

Real derfct1_1_P1_1D( const GeoVector & );
Real derfct2_1_P1_1D( const GeoVector& );

Real der2fct1_P1_1D( const GeoVector& );

static const Real refcoor_P1_1D[ 6 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.
};

static const RefEle::function_Type fct_P1_1D[ 2 ] =
{
    fct1_P1_1D, fct2_P1_1D
};
static const RefEle::function_Type derfct_P1_1D[ 2 ] =
{
    derfct1_1_P1_1D, derfct2_1_P1_1D
};
static const RefEle::function_Type der2fct_P1_1D[ 2 ] =
{
    der2fct1_P1_1D, der2fct1_P1_1D
};

//======================================================================
//
//                            P2  (1D)
//
//======================================================================
/*
                           1--3--2
*/
Real fct1_P2_1D( const GeoVector& v );
Real fct2_P2_1D( const GeoVector& v );
Real fct3_P2_1D( const GeoVector& v );

Real derfct1_1_P2_1D( const GeoVector& v );
Real derfct2_1_P2_1D( const GeoVector& v );
Real derfct3_1_P2_1D( const GeoVector& v );

Real der2fct1_11_P2_1D( const GeoVector& v );
Real der2fct2_11_P2_1D( const GeoVector& v );
Real der2fct3_11_P2_1D( const GeoVector& v );

static const Real refcoor_P2_1D[ 9 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    0.5 , 0. , 0.
};
static const RefEle::function_Type fct_P2_1D[ 3 ] =
{
    fct1_P2_1D, fct2_P2_1D, fct3_P2_1D
};
static const RefEle::function_Type derfct_P2_1D[ 3 ] =
{
    derfct1_1_P2_1D, derfct2_1_P2_1D, derfct3_1_P2_1D
};
static const RefEle::function_Type der2fct_P2_1D[ 3 ] =
{
    der2fct1_11_P2_1D, der2fct2_11_P2_1D, der2fct3_11_P2_1D
};

//======================================================================
//
//                            P0  (2D)
//
//======================================================================
/*

                           |\
                           | \
                           | 1\
                            ---
*/
Real fct1_P0_2D( const GeoVector& v );
// First and Second derivatives are both equal (to 0).
Real derfct1_P0_2D( const GeoVector & );
Real der2fct1_P0_2D( const GeoVector & );

static const Real refcoor_P0_2D[ 3 ] =
{
    1. / 3. , 1. / 3. , 0.
}
;  // check this : gravity center??

static const RefEle::function_Type fct_P0_2D[ 1 ] =
{
    fct1_P0_2D
};

static const RefEle::function_Type derfct_P0_2D[ 2 ] =
{
    derfct1_P0_2D, derfct1_P0_2D
};

static const RefEle::function_Type der2fct_P0_2D[ 4 ] =
{
    der2fct1_P0_2D, der2fct1_P0_2D,
    der2fct1_P0_2D, der2fct1_P0_2D
};

//======================================================================
//
//                            P1  (2D)
//
//======================================================================
/*
                           3
                           |\
                           | \
                           |  \
                           1---2
*/
Real fct1_P1_2D( const GeoVector& v );
Real fct2_P1_2D( const GeoVector& v );
Real fct3_P1_2D( const GeoVector& v );

Real derfct1_1_P1_2D( const GeoVector& );
Real derfct1_2_P1_2D( const GeoVector& );
Real derfct2_1_P1_2D( const GeoVector& );
Real derfct2_2_P1_2D( const GeoVector& );
Real derfct3_1_P1_2D( const GeoVector& );
Real derfct3_2_P1_2D( const GeoVector& );

// Second derivatives
Real der2fctx_xx_P1_2D( const GeoVector& );

static const Real refcoor_P1_2D[ 9 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    0. , 1. , 0.
};

static const RefEle::function_Type fct_P1_2D[ 3 ] =
{
    fct1_P1_2D, fct2_P1_2D, fct3_P1_2D
};

static const RefEle::function_Type derfct_P1_2D[ 6 ] =
{
    derfct1_1_P1_2D, derfct1_2_P1_2D,
    derfct2_1_P1_2D, derfct2_2_P1_2D,
    derfct3_1_P1_2D, derfct3_2_P1_2D
};
static const RefEle::function_Type der2fct_P1_2D[ 12 ] =
{
    der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D,
    der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D,
    der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D, der2fctx_xx_P1_2D
};
//======================================================================
//
//                            P2  (2D)
//
//======================================================================
/*
                           3
                           |\
                           6 5
                           |  \
                           1-4-2
*/
Real fct1_P2_2D( const GeoVector& v );
Real fct2_P2_2D( const GeoVector& v );
Real fct3_P2_2D( const GeoVector& v );
Real fct4_P2_2D( const GeoVector& v );
Real fct5_P2_2D( const GeoVector& v );
Real fct6_P2_2D( const GeoVector& v );

Real derfct1_1_P2_2D( const GeoVector& v );
Real derfct1_2_P2_2D( const GeoVector& v );
Real derfct2_1_P2_2D( const GeoVector& v );
Real derfct2_2_P2_2D( const GeoVector& v );
Real derfct3_1_P2_2D( const GeoVector& v );
Real derfct3_2_P2_2D( const GeoVector& v );
Real derfct4_1_P2_2D( const GeoVector& v );
Real derfct4_2_P2_2D( const GeoVector& v );
Real derfct5_1_P2_2D( const GeoVector& v );
Real derfct5_2_P2_2D( const GeoVector& v );
Real derfct6_1_P2_2D( const GeoVector& v );
Real derfct6_2_P2_2D( const GeoVector& v );

Real der2fct1_11_P2_2D( const GeoVector& v );
Real der2fct1_12_P2_2D( const GeoVector& v );
Real der2fct1_21_P2_2D( const GeoVector& v );
Real der2fct1_22_P2_2D( const GeoVector& v );

Real der2fct2_11_P2_2D( const GeoVector& v );
Real der2fct2_12_P2_2D( const GeoVector& v );
Real der2fct2_21_P2_2D( const GeoVector& v );
Real der2fct2_22_P2_2D( const GeoVector& v );

Real der2fct3_11_P2_2D( const GeoVector& v );
Real der2fct3_12_P2_2D( const GeoVector& v );
Real der2fct3_21_P2_2D( const GeoVector& v );
Real der2fct3_22_P2_2D( const GeoVector& v );

Real der2fct4_11_P2_2D( const GeoVector& v );
Real der2fct4_12_P2_2D( const GeoVector& v );
Real der2fct4_21_P2_2D( const GeoVector& v );
Real der2fct4_22_P2_2D( const GeoVector& v );

Real der2fct5_11_P2_2D( const GeoVector& v );
Real der2fct5_12_P2_2D( const GeoVector& v );
Real der2fct5_21_P2_2D( const GeoVector& v );
Real der2fct5_22_P2_2D( const GeoVector& v );

Real der2fct6_11_P2_2D( const GeoVector& v );
Real der2fct6_12_P2_2D( const GeoVector& v );
Real der2fct6_21_P2_2D( const GeoVector& v );
Real der2fct6_22_P2_2D( const GeoVector& v );

static const Real refcoor_P2_2D[ 18 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    0. , 1. , 0.,
    0.5 , 0. , 0.,
    0.5 , 0.5 , 0.,
    0. , 0.5 , 0.
};

static const RefEle::function_Type fct_P2_2D[ 6 ] =
{
    fct1_P2_2D, fct2_P2_2D, fct3_P2_2D,
    fct4_P2_2D, fct5_P2_2D, fct6_P2_2D
};

static const RefEle::function_Type derfct_P2_2D[ 12 ] =
{
    derfct1_1_P2_2D, derfct1_2_P2_2D,
    derfct2_1_P2_2D, derfct2_2_P2_2D,
    derfct3_1_P2_2D, derfct3_2_P2_2D,
    derfct4_1_P2_2D, derfct4_2_P2_2D,
    derfct5_1_P2_2D, derfct5_2_P2_2D,
    derfct6_1_P2_2D, derfct6_2_P2_2D
};
static const RefEle::function_Type der2fct_P2_2D[ 24 ] =
{
    der2fct1_11_P2_2D, der2fct1_12_P2_2D, der2fct1_21_P2_2D, der2fct1_22_P2_2D,
    der2fct2_11_P2_2D, der2fct2_12_P2_2D, der2fct2_21_P2_2D, der2fct2_22_P2_2D,
    der2fct3_11_P2_2D, der2fct3_12_P2_2D, der2fct3_21_P2_2D, der2fct3_22_P2_2D,
    der2fct4_11_P2_2D, der2fct4_12_P2_2D, der2fct4_21_P2_2D, der2fct4_22_P2_2D,
    der2fct5_11_P2_2D, der2fct5_12_P2_2D, der2fct5_21_P2_2D, der2fct5_22_P2_2D,
    der2fct6_11_P2_2D, der2fct6_12_P2_2D, der2fct6_21_P2_2D, der2fct6_22_P2_2D
};


//======================================================================
//
//                            Q0  (2D)
//
//======================================================================
/*
                            -------
                           |       |
                           |   1   |
                           |       |
                            -------

*/
Real fct1_Q0_2D( const GeoVector& v );
Real derfct1_Q0_2D( const GeoVector& v );
// The second derivative is equal to the first : both = 0.
Real der2fct1_Q0_2D( const GeoVector& v );

static const Real refcoor_Q0_2D[ 3 ] =
{
    0.5, 0.5, 0.
};

static const RefEle::function_Type fct_Q0_2D[ 1 ] =
{
    fct1_Q0_2D
};

static const RefEle::function_Type derfct_Q0_2D[ 2 ] =
{
    derfct1_Q0_2D, derfct1_Q0_2D
};

static const RefEle::function_Type der2fct_Q0_2D[ 4 ] =
{
    der2fct1_Q0_2D, der2fct1_Q0_2D,
    der2fct1_Q0_2D, der2fct1_Q0_2D
};

//======================================================================
//
//                            Q1  (2D)
//
//======================================================================
/*
                           4-------3
                           |       |
                           |       |
                           |       |
                           1-------2
*/
Real fct1_Q1_2D( const GeoVector& v );
Real fct2_Q1_2D( const GeoVector& v );
Real fct3_Q1_2D( const GeoVector& v );
Real fct4_Q1_2D( const GeoVector& v );

Real derfct1_1_Q1_2D( const GeoVector& v );
Real derfct1_2_Q1_2D( const GeoVector& v );
Real derfct2_1_Q1_2D( const GeoVector& v );
Real derfct2_2_Q1_2D( const GeoVector& v );
Real derfct3_1_Q1_2D( const GeoVector& v );
Real derfct3_2_Q1_2D( const GeoVector& v );
Real derfct4_1_Q1_2D( const GeoVector& v );
Real derfct4_2_Q1_2D( const GeoVector& v );

// Second derivatives
Real der2fctx_xx_Q1_2D( const GeoVector& );

static const Real refcoor_Q1_2D[ 12 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    1. , 1. , 0.,
    0. , 1. , 0.
};

static const RefEle::function_Type fct_Q1_2D[ 4 ] =
{
    fct1_Q1_2D, fct2_Q1_2D, fct3_Q1_2D, fct4_Q1_2D
};

static const RefEle::function_Type derfct_Q1_2D[ 8 ] =
{
    derfct1_1_Q1_2D, derfct1_2_Q1_2D,
    derfct2_1_Q1_2D, derfct2_2_Q1_2D,
    derfct3_1_Q1_2D, derfct3_2_Q1_2D,
    derfct4_1_Q1_2D, derfct4_2_Q1_2D
};
static const RefEle::function_Type der2fct_Q1_2D[ 16 ] =
{
    der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D,
    der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D,
    der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D,
    der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D, der2fctx_xx_Q1_2D
};

//======================================================================
//
//                            Q2  (2D)
//
//======================================================================
/*
                           4---7---3
                           |       |
                           8   9   6
                           |       |
                           1---5---2
*/
Real fct1_Q2_2D( const GeoVector& v );
Real fct5_Q2_2D( const GeoVector& v );
Real fct2_Q2_2D( const GeoVector& v );
Real fct6_Q2_2D( const GeoVector& v );
Real fct3_Q2_2D( const GeoVector& v );
Real fct7_Q2_2D( const GeoVector& v );
Real fct4_Q2_2D( const GeoVector& v );
Real fct8_Q2_2D( const GeoVector& v );
Real fct9_Q2_2D( const GeoVector& v );

Real derfct1_1_Q2_2D( const GeoVector& v );
Real derfct1_2_Q2_2D( const GeoVector& v );
Real derfct5_1_Q2_2D( const GeoVector& v );
Real derfct5_2_Q2_2D( const GeoVector& v );
Real derfct2_1_Q2_2D( const GeoVector& v );
Real derfct2_2_Q2_2D( const GeoVector& v );
Real derfct6_1_Q2_2D( const GeoVector& v );
Real derfct6_2_Q2_2D( const GeoVector& v );
Real derfct3_1_Q2_2D( const GeoVector& v );
Real derfct3_2_Q2_2D( const GeoVector& v );
Real derfct7_1_Q2_2D( const GeoVector& v );
Real derfct7_2_Q2_2D( const GeoVector& v );
Real derfct4_1_Q2_2D( const GeoVector& v );
Real derfct4_2_Q2_2D( const GeoVector& v );
Real derfct8_1_Q2_2D( const GeoVector& v );
Real derfct8_2_Q2_2D( const GeoVector& v );
Real derfct9_1_Q2_2D( const GeoVector& v );
Real derfct9_2_Q2_2D( const GeoVector& v );

Real der2fct1_11_Q2_2D( const GeoVector& v );
Real der2fct1_12_Q2_2D( const GeoVector& v );
Real der2fct1_21_Q2_2D( const GeoVector& v );
Real der2fct1_22_Q2_2D( const GeoVector& v );

Real der2fct5_11_Q2_2D( const GeoVector& v );
Real der2fct5_12_Q2_2D( const GeoVector& v );
Real der2fct5_21_Q2_2D( const GeoVector& v );
Real der2fct5_22_Q2_2D( const GeoVector& v );

Real der2fct2_11_Q2_2D( const GeoVector& v );
Real der2fct2_12_Q2_2D( const GeoVector& v );
Real der2fct2_21_Q2_2D( const GeoVector& v );
Real der2fct2_22_Q2_2D( const GeoVector& v );

Real der2fct6_11_Q2_2D( const GeoVector& v );
Real der2fct6_12_Q2_2D( const GeoVector& v );
Real der2fct6_21_Q2_2D( const GeoVector& v );
Real der2fct6_22_Q2_2D( const GeoVector& v );

Real der2fct3_11_Q2_2D( const GeoVector& v );
Real der2fct3_12_Q2_2D( const GeoVector& v );
Real der2fct3_21_Q2_2D( const GeoVector& v );
Real der2fct3_22_Q2_2D( const GeoVector& v );

Real der2fct7_11_Q2_2D( const GeoVector& v );
Real der2fct7_12_Q2_2D( const GeoVector& v );
Real der2fct7_21_Q2_2D( const GeoVector& v );
Real der2fct7_22_Q2_2D( const GeoVector& v );

Real der2fct4_11_Q2_2D( const GeoVector& v );
Real der2fct4_12_Q2_2D( const GeoVector& v );
Real der2fct4_21_Q2_2D( const GeoVector& v );
Real der2fct4_22_Q2_2D( const GeoVector& v );

Real der2fct8_11_Q2_2D( const GeoVector& v );
Real der2fct8_12_Q2_2D( const GeoVector& v );
Real der2fct8_21_Q2_2D( const GeoVector& v );
Real der2fct8_22_Q2_2D( const GeoVector& v );

Real der2fct9_11_Q2_2D( const GeoVector& v );
Real der2fct9_12_Q2_2D( const GeoVector& v );
Real der2fct9_21_Q2_2D( const GeoVector& v );
Real der2fct9_22_Q2_2D( const GeoVector& v );


static const Real refcoor_Q2_2D[ 27 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    1. , 1. , 0.,
    0. , 1. , 0.,
    0.5 , 0. , 0.,
    1. , 0.5 , 0.,
    0.5 , 1. , 0.,
    0. , 0.5 , 0.,
    0.5 , 0.5 , 0.
};

static const RefEle::function_Type fct_Q2_2D[ 9 ] =
{
    fct1_Q2_2D, fct2_Q2_2D, fct3_Q2_2D, fct4_Q2_2D,
    fct5_Q2_2D, fct6_Q2_2D, fct7_Q2_2D, fct8_Q2_2D,
    fct9_Q2_2D
};


static const RefEle::function_Type derfct_Q2_2D[ 18 ] =
{
    derfct1_1_Q2_2D, derfct1_2_Q2_2D,
    derfct2_1_Q2_2D, derfct2_2_Q2_2D,
    derfct3_1_Q2_2D, derfct3_2_Q2_2D,
    derfct4_1_Q2_2D, derfct4_2_Q2_2D,
    derfct5_1_Q2_2D, derfct5_2_Q2_2D,
    derfct6_1_Q2_2D, derfct6_2_Q2_2D,
    derfct7_1_Q2_2D, derfct7_2_Q2_2D,
    derfct8_1_Q2_2D, derfct8_2_Q2_2D,
    derfct9_1_Q2_2D, derfct9_2_Q2_2D
};

static const RefEle::function_Type der2fct_Q2_2D[ 36 ] =
{
    der2fct1_11_Q2_2D, der2fct1_12_Q2_2D, der2fct1_21_Q2_2D, der2fct1_22_Q2_2D,
    der2fct2_11_Q2_2D, der2fct2_12_Q2_2D, der2fct2_21_Q2_2D, der2fct2_22_Q2_2D,
    der2fct3_11_Q2_2D, der2fct3_12_Q2_2D, der2fct3_21_Q2_2D, der2fct3_22_Q2_2D,
    der2fct4_11_Q2_2D, der2fct4_12_Q2_2D, der2fct4_21_Q2_2D, der2fct4_22_Q2_2D,
    der2fct5_11_Q2_2D, der2fct5_12_Q2_2D, der2fct5_21_Q2_2D, der2fct5_22_Q2_2D,
    der2fct6_11_Q2_2D, der2fct6_12_Q2_2D, der2fct6_21_Q2_2D, der2fct6_22_Q2_2D,
    der2fct7_11_Q2_2D, der2fct7_12_Q2_2D, der2fct7_21_Q2_2D, der2fct7_22_Q2_2D,
    der2fct8_11_Q2_2D, der2fct8_12_Q2_2D, der2fct8_21_Q2_2D, der2fct8_22_Q2_2D,
    der2fct9_11_Q2_2D, der2fct9_12_Q2_2D, der2fct9_21_Q2_2D, der2fct9_22_Q2_2D
};

//======================================================================
//
//                            P0  (3D)
//
//======================================================================
/*

               / .
              /  \.
             /  . \\
            / . 1  \\
           /.       \!
           ----------
*/
Real fct1_P0_3D( const GeoVector& v );

Real derfct1_P0_3D( const GeoVector& );

// Second derivatives
Real der2fct1_P0_3D( const GeoVector& );

static const Real refcoor_P0_3D[ 3 ] =
{
    0.25 , 0.25 , 0.25
};

static const RefEle::function_Type fct_P0_3D[ 1 ] =
{
    fct1_P0_3D
};

static const RefEle::function_Type derfct_P0_3D[ 3 ] =
{
    derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D
};
static const RefEle::function_Type der2fct_P0_3D[ 9 ] =
{
    derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D,
    derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D,
    derfct1_P0_3D, derfct1_P0_3D, derfct1_P0_3D
};

//======================================================================
//
//                            P1  (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / .    \\
           /.       \!
         1 ----------2
*/
Real fct1_P1_3D( const GeoVector& v );
Real fct2_P1_3D( const GeoVector& v );
Real fct3_P1_3D( const GeoVector& v );
Real fct4_P1_3D( const GeoVector& v );

Real derfct1_1_P1_3D( const GeoVector& );
Real derfct1_2_P1_3D( const GeoVector& );
Real derfct1_3_P1_3D( const GeoVector& );
Real derfct2_1_P1_3D( const GeoVector& );
Real derfct2_2_P1_3D( const GeoVector& );
Real derfct2_3_P1_3D( const GeoVector& );
Real derfct3_1_P1_3D( const GeoVector& );
Real derfct3_2_P1_3D( const GeoVector& );
Real derfct3_3_P1_3D( const GeoVector& );
Real derfct4_1_P1_3D( const GeoVector& );
Real derfct4_2_P1_3D( const GeoVector& );
Real derfct4_3_P1_3D( const GeoVector& );

// Second derivatives
Real der2fctx_xx_P1_3D( const GeoVector& );

static const Real refcoor_P1_3D[ 12 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    0. , 1. , 0.,
    0. , 0. , 1.
};

static const RefEle::function_Type fct_P1_3D[ 4 ] =
{
    fct1_P1_3D, fct2_P1_3D, fct3_P1_3D, fct4_P1_3D
};

static const RefEle::function_Type derfct_P1_3D[ 12 ] =
{
    derfct1_1_P1_3D, derfct1_2_P1_3D, derfct1_3_P1_3D,
    derfct2_1_P1_3D, derfct2_2_P1_3D, derfct2_3_P1_3D,
    derfct3_1_P1_3D, derfct3_2_P1_3D, derfct3_3_P1_3D,
    derfct4_1_P1_3D, derfct4_2_P1_3D, derfct4_3_P1_3D
};
static const RefEle::function_Type der2fct_P1_3D[ 36 ] =
{
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D,
    der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D, der2fctx_xx_P1_3D
};

//======================================================================
//
//                            P1 Bubble (3D)
//
//======================================================================
/*
                4
               / .
              /  \.3
             /  . \\
            / . .5 \\
           /.       \!
         1 ----------2
*/
Real fct1_P1bubble_3D( const GeoVector& v );
Real fct2_P1bubble_3D( const GeoVector& v );
Real fct3_P1bubble_3D( const GeoVector& v );
Real fct4_P1bubble_3D( const GeoVector& v );
Real fct5_P1bubble_3D( const GeoVector& v );

Real derfct1_1_P1bubble_3D( const GeoVector& );
Real derfct1_2_P1bubble_3D( const GeoVector& );
Real derfct1_3_P1bubble_3D( const GeoVector& );
Real derfct2_1_P1bubble_3D( const GeoVector& );
Real derfct2_2_P1bubble_3D( const GeoVector& );
Real derfct2_3_P1bubble_3D( const GeoVector& );
Real derfct3_1_P1bubble_3D( const GeoVector& );
Real derfct3_2_P1bubble_3D( const GeoVector& );
Real derfct3_3_P1bubble_3D( const GeoVector& );
Real derfct4_1_P1bubble_3D( const GeoVector& );
Real derfct4_2_P1bubble_3D( const GeoVector& );
Real derfct4_3_P1bubble_3D( const GeoVector& );
Real derfct5_1_P1bubble_3D( const GeoVector& );
Real derfct5_2_P1bubble_3D( const GeoVector& );
Real derfct5_3_P1bubble_3D( const GeoVector& );

// Second derivatives
Real der2fctx_xx_P1bubble_3D( const GeoVector& );
Real der2fct5_11_P1bubble_3D( const GeoVector& );
Real der2fct5_12_P1bubble_3D( const GeoVector& );
Real der2fct5_13_P1bubble_3D( const GeoVector& );
Real der2fct5_21_P1bubble_3D( const GeoVector& );
Real der2fct5_22_P1bubble_3D( const GeoVector& );
Real der2fct5_23_P1bubble_3D( const GeoVector& );
Real der2fct5_31_P1bubble_3D( const GeoVector& );
Real der2fct5_32_P1bubble_3D( const GeoVector& );
Real der2fct5_33_P1bubble_3D( const GeoVector& );

static const Real refcoor_P1bubble_3D[ 15 ] =
{
    0. , 0. , 0.,
    1. , 0. , 0.,
    0. , 1. , 0.,
    0. , 0. , 1.,
    0.25, 0.25, 0.25
};

static const RefEle::function_Type fct_P1bubble_3D[ 5 ] =
{
    fct1_P1bubble_3D, fct2_P1bubble_3D, fct3_P1bubble_3D, fct4_P1bubble_3D, fct5_P1bubble_3D
};

static const RefEle::function_Type derfct_P1bubble_3D[ 15 ] =
{
    derfct1_1_P1bubble_3D, derfct1_2_P1bubble_3D, derfct1_3_P1bubble_3D,
    derfct2_1_P1bubble_3D, derfct2_2_P1bubble_3D, derfct2_3_P1bubble_3D,
    derfct3_1_P1bubble_3D, derfct3_2_P1bubble_3D, derfct3_3_P1bubble_3D,
    derfct4_1_P1bubble_3D, derfct4_2_P1bubble_3D, derfct4_3_P1bubble_3D,
    derfct5_1_P1bubble_3D, derfct5_2_P1bubble_3D, derfct5_3_P1bubble_3D
};
static const RefEle::function_Type der2fct_P1bubble_3D[ 45 ] =
{
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D, der2fctx_xx_P1bubble_3D,
    der2fct5_11_P1bubble_3D, der2fct5_12_P1bubble_3D, der2fct5_13_P1bubble_3D,
    der2fct5_21_P1bubble_3D, der2fct5_22_P1bubble_3D, der2fct5_23_P1bubble_3D,
    der2fct5_31_P1bubble_3D, der2fct5_32_P1bubble_3D, der2fct5_33_P1bubble_3D
};

//======================================================================
//
//                            P2  (3D)
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7    \6
           /.       \!
         1 -----5----2
*/
Real fct1_P2_3D( const GeoVector& v );
Real fct2_P2_3D( const GeoVector& v );
Real fct3_P2_3D( const GeoVector& v );
Real fct4_P2_3D( const GeoVector& v );
Real fct5_P2_3D( const GeoVector& v );
Real fct6_P2_3D( const GeoVector& v );
Real fct7_P2_3D( const GeoVector& v );
Real fct8_P2_3D( const GeoVector& v );
Real fct9_P2_3D( const GeoVector& v );
Real fct10_P2_3D( const GeoVector& v );


Real derfct1_1_P2_3D( const GeoVector& v );
Real derfct1_2_P2_3D( const GeoVector& v );
Real derfct1_3_P2_3D( const GeoVector& v );

Real derfct2_1_P2_3D( const GeoVector& v );
Real derfct2_2_P2_3D( const GeoVector& v );
Real derfct2_3_P2_3D( const GeoVector& v );

Real derfct3_1_P2_3D( const GeoVector& v );
Real derfct3_2_P2_3D( const GeoVector& v );
Real derfct3_3_P2_3D( const GeoVector& v );

Real derfct4_1_P2_3D( const GeoVector& v );
Real derfct4_2_P2_3D( const GeoVector& v );
Real derfct4_3_P2_3D( const GeoVector& v );

Real derfct5_1_P2_3D( const GeoVector& v );
Real derfct5_2_P2_3D( const GeoVector& v );
Real derfct5_3_P2_3D( const GeoVector& v );

Real derfct6_1_P2_3D( const GeoVector& v );
Real derfct6_2_P2_3D( const GeoVector& v );
Real derfct6_3_P2_3D( const GeoVector& v );

Real derfct7_1_P2_3D( const GeoVector& v );
Real derfct7_2_P2_3D( const GeoVector& v );
Real derfct7_3_P2_3D( const GeoVector& v );

Real derfct8_1_P2_3D( const GeoVector& v );
Real derfct8_2_P2_3D( const GeoVector& v );
Real derfct8_3_P2_3D( const GeoVector& v );

Real derfct9_1_P2_3D( const GeoVector& v );
Real derfct9_2_P2_3D( const GeoVector& v );
Real derfct9_3_P2_3D( const GeoVector& v );

Real derfct10_1_P2_3D( const GeoVector& v );
Real derfct10_2_P2_3D( const GeoVector& v );
Real derfct10_3_P2_3D( const GeoVector& v );


Real der2fct1_11_P2_3D( const GeoVector& v );
Real der2fct1_12_P2_3D( const GeoVector& v );
Real der2fct1_13_P2_3D( const GeoVector& v );
Real der2fct1_21_P2_3D( const GeoVector& v );
Real der2fct1_22_P2_3D( const GeoVector& v );
Real der2fct1_23_P2_3D( const GeoVector& v );
Real der2fct1_31_P2_3D( const GeoVector& v );
Real der2fct1_32_P2_3D( const GeoVector& v );
Real der2fct1_33_P2_3D( const GeoVector& v );

Real der2fct2_11_P2_3D( const GeoVector& v );
Real der2fct2_12_P2_3D( const GeoVector& v );
Real der2fct2_13_P2_3D( const GeoVector& v );
Real der2fct2_21_P2_3D( const GeoVector& v );
Real der2fct2_22_P2_3D( const GeoVector& v );
Real der2fct2_23_P2_3D( const GeoVector& v );
Real der2fct2_31_P2_3D( const GeoVector& v );
Real der2fct2_32_P2_3D( const GeoVector& v );
Real der2fct2_33_P2_3D( const GeoVector& v );

Real der2fct3_11_P2_3D( const GeoVector& v );
Real der2fct3_12_P2_3D( const GeoVector& v );
Real der2fct3_13_P2_3D( const GeoVector& v );
Real der2fct3_21_P2_3D( const GeoVector& v );
Real der2fct3_22_P2_3D( const GeoVector& v );
Real der2fct3_23_P2_3D( const GeoVector& v );
Real der2fct3_31_P2_3D( const GeoVector& v );
Real der2fct3_32_P2_3D( const GeoVector& v );
Real der2fct3_33_P2_3D( const GeoVector& v );

Real der2fct4_11_P2_3D( const GeoVector& v );
Real der2fct4_12_P2_3D( const GeoVector& v );
Real der2fct4_13_P2_3D( const GeoVector& v );
Real der2fct4_21_P2_3D( const GeoVector& v );
Real der2fct4_22_P2_3D( const GeoVector& v );
Real der2fct4_23_P2_3D( const GeoVector& v );
Real der2fct4_31_P2_3D( const GeoVector& v );
Real der2fct4_32_P2_3D( const GeoVector& v );
Real der2fct4_33_P2_3D( const GeoVector& v );

Real der2fct5_11_P2_3D( const GeoVector& v );
Real der2fct5_12_P2_3D( const GeoVector& v );
Real der2fct5_13_P2_3D( const GeoVector& v );
Real der2fct5_21_P2_3D( const GeoVector& v );
Real der2fct5_22_P2_3D( const GeoVector& v );
Real der2fct5_23_P2_3D( const GeoVector& v );
Real der2fct5_31_P2_3D( const GeoVector& v );
Real der2fct5_32_P2_3D( const GeoVector& v );
Real der2fct5_33_P2_3D( const GeoVector& v );

Real der2fct6_11_P2_3D( const GeoVector& v );
Real der2fct6_12_P2_3D( const GeoVector& v );
Real der2fct6_13_P2_3D( const GeoVector& v );
Real der2fct6_21_P2_3D( const GeoVector& v );
Real der2fct6_22_P2_3D( const GeoVector& v );
Real der2fct6_23_P2_3D( const GeoVector& v );
Real der2fct6_31_P2_3D( const GeoVector& v );
Real der2fct6_32_P2_3D( const GeoVector& v );
Real der2fct6_33_P2_3D( const GeoVector& v );

Real der2fct7_11_P2_3D( const GeoVector& v );
Real der2fct7_12_P2_3D( const GeoVector& v );
Real der2fct7_13_P2_3D( const GeoVector& v );
Real der2fct7_21_P2_3D( const GeoVector& v );
Real der2fct7_22_P2_3D( const GeoVector& v );
Real der2fct7_23_P2_3D( const GeoVector& v );
Real der2fct7_31_P2_3D( const GeoVector& v );
Real der2fct7_32_P2_3D( const GeoVector& v );
Real der2fct7_33_P2_3D( const GeoVector& v );

Real der2fct8_11_P2_3D( const GeoVector& v );
Real der2fct8_12_P2_3D( const GeoVector& v );
Real der2fct8_13_P2_3D( const GeoVector& v );
Real der2fct8_21_P2_3D( const GeoVector& v );
Real der2fct8_22_P2_3D( const GeoVector& v );
Real der2fct8_23_P2_3D( const GeoVector& v );
Real der2fct8_31_P2_3D( const GeoVector& v );
Real der2fct8_32_P2_3D( const GeoVector& v );
Real der2fct8_33_P2_3D( const GeoVector& v );

Real der2fct9_11_P2_3D( const GeoVector& v );
Real der2fct9_12_P2_3D( const GeoVector& v );
Real der2fct9_13_P2_3D( const GeoVector& v );
Real der2fct9_21_P2_3D( const GeoVector& v );
Real der2fct9_22_P2_3D( const GeoVector& v );
Real der2fct9_23_P2_3D( const GeoVector& v );
Real der2fct9_31_P2_3D( const GeoVector& v );
Real der2fct9_32_P2_3D( const GeoVector& v );
Real der2fct9_33_P2_3D( const GeoVector& v );

Real der2fct10_11_P2_3D( const GeoVector& v );
Real der2fct10_12_P2_3D( const GeoVector& v );
Real der2fct10_13_P2_3D( const GeoVector& v );
Real der2fct10_21_P2_3D( const GeoVector& v );
Real der2fct10_22_P2_3D( const GeoVector& v );
Real der2fct10_23_P2_3D( const GeoVector& v );
Real der2fct10_31_P2_3D( const GeoVector& v );
Real der2fct10_32_P2_3D( const GeoVector& v );
Real der2fct10_33_P2_3D( const GeoVector& v );


static const Real refcoor_P2_3D[ 30 ] =
{
    0. , 0. , 0. ,
    1. , 0. , 0. ,
    0. , 1. , 0. ,
    0. , 0. , 1. ,
    0.5 , 0. , 0. ,
    0.5, 0.5 , 0. ,
    0. , 0.5 , 0. ,
    0. , 0. , 0.5,
    0.5, 0. , 0.5,
    0. , 0.5 , 0.5
};

static const RefEle::function_Type fct_P2_3D[ 10 ] =
{
    fct1_P2_3D, fct2_P2_3D, fct3_P2_3D, fct4_P2_3D,
    fct5_P2_3D, fct6_P2_3D, fct7_P2_3D, fct8_P2_3D,
    fct9_P2_3D, fct10_P2_3D
};

static const RefEle::function_Type derfct_P2_3D[ 30 ] =
{
    derfct1_1_P2_3D, derfct1_2_P2_3D, derfct1_3_P2_3D,
    derfct2_1_P2_3D, derfct2_2_P2_3D, derfct2_3_P2_3D,
    derfct3_1_P2_3D, derfct3_2_P2_3D, derfct3_3_P2_3D,
    derfct4_1_P2_3D, derfct4_2_P2_3D, derfct4_3_P2_3D,
    derfct5_1_P2_3D, derfct5_2_P2_3D, derfct5_3_P2_3D,
    derfct6_1_P2_3D, derfct6_2_P2_3D, derfct6_3_P2_3D,
    derfct7_1_P2_3D, derfct7_2_P2_3D, derfct7_3_P2_3D,
    derfct8_1_P2_3D, derfct8_2_P2_3D, derfct8_3_P2_3D,
    derfct9_1_P2_3D, derfct9_2_P2_3D, derfct9_3_P2_3D,
    derfct10_1_P2_3D, derfct10_2_P2_3D, derfct10_3_P2_3D
};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2\_3D, der2fct$i\_$j"."2\_P2\_3D, der2fct$i\_$j"."3_P2\_3D,\n";
  }
}
*/
static const RefEle::function_Type der2fct_P2_3D[ 90 ] =
{
    der2fct1_11_P2_3D, der2fct1_12_P2_3D, der2fct1_13_P2_3D,
    der2fct1_21_P2_3D, der2fct1_22_P2_3D, der2fct1_23_P2_3D,
    der2fct1_31_P2_3D, der2fct1_32_P2_3D, der2fct1_33_P2_3D,
    der2fct2_11_P2_3D, der2fct2_12_P2_3D, der2fct2_13_P2_3D,
    der2fct2_21_P2_3D, der2fct2_22_P2_3D, der2fct2_23_P2_3D,
    der2fct2_31_P2_3D, der2fct2_32_P2_3D, der2fct2_33_P2_3D,
    der2fct3_11_P2_3D, der2fct3_12_P2_3D, der2fct3_13_P2_3D,
    der2fct3_21_P2_3D, der2fct3_22_P2_3D, der2fct3_23_P2_3D,
    der2fct3_31_P2_3D, der2fct3_32_P2_3D, der2fct3_33_P2_3D,
    der2fct4_11_P2_3D, der2fct4_12_P2_3D, der2fct4_13_P2_3D,
    der2fct4_21_P2_3D, der2fct4_22_P2_3D, der2fct4_23_P2_3D,
    der2fct4_31_P2_3D, der2fct4_32_P2_3D, der2fct4_33_P2_3D,
    der2fct5_11_P2_3D, der2fct5_12_P2_3D, der2fct5_13_P2_3D,
    der2fct5_21_P2_3D, der2fct5_22_P2_3D, der2fct5_23_P2_3D,
    der2fct5_31_P2_3D, der2fct5_32_P2_3D, der2fct5_33_P2_3D,
    der2fct6_11_P2_3D, der2fct6_12_P2_3D, der2fct6_13_P2_3D,
    der2fct6_21_P2_3D, der2fct6_22_P2_3D, der2fct6_23_P2_3D,
    der2fct6_31_P2_3D, der2fct6_32_P2_3D, der2fct6_33_P2_3D,
    der2fct7_11_P2_3D, der2fct7_12_P2_3D, der2fct7_13_P2_3D,
    der2fct7_21_P2_3D, der2fct7_22_P2_3D, der2fct7_23_P2_3D,
    der2fct7_31_P2_3D, der2fct7_32_P2_3D, der2fct7_33_P2_3D,
    der2fct8_11_P2_3D, der2fct8_12_P2_3D, der2fct8_13_P2_3D,
    der2fct8_21_P2_3D, der2fct8_22_P2_3D, der2fct8_23_P2_3D,
    der2fct8_31_P2_3D, der2fct8_32_P2_3D, der2fct8_33_P2_3D,
    der2fct9_11_P2_3D, der2fct9_12_P2_3D, der2fct9_13_P2_3D,
    der2fct9_21_P2_3D, der2fct9_22_P2_3D, der2fct9_23_P2_3D,
    der2fct9_31_P2_3D, der2fct9_32_P2_3D, der2fct9_33_P2_3D,
    der2fct10_11_P2_3D, der2fct10_12_P2_3D, der2fct10_13_P2_3D,
    der2fct10_21_P2_3D, der2fct10_22_P2_3D, der2fct10_23_P2_3D,
    der2fct10_31_P2_3D, der2fct10_32_P2_3D, der2fct10_33_P2_3D
};
//======================================================================
//
//                            P2tilde  (3D)
// NAVIER-STOKES P2 Basis Oriented to the mass lumping
//
//======================================================================
/*
                4
               / .10
              /  \.3
             8  . 9\
            / 7 .11 6
           /.       \!
         1 -----5----2
*/
Real fct1_P2tilde_3D( const GeoVector& v );
Real fct2_P2tilde_3D( const GeoVector& v );
Real fct3_P2tilde_3D( const GeoVector& v );
Real fct4_P2tilde_3D( const GeoVector& v );
Real fct5_P2tilde_3D( const GeoVector& v );
Real fct6_P2tilde_3D( const GeoVector& v );
Real fct7_P2tilde_3D( const GeoVector& v );
Real fct8_P2tilde_3D( const GeoVector& v );
Real fct9_P2tilde_3D( const GeoVector& v );
Real fct10_P2tilde_3D( const GeoVector& v );
Real fct11_P2tilde_3D( const GeoVector& v );


Real derfct1_1_P2tilde_3D( const GeoVector& v );
Real derfct1_2_P2tilde_3D( const GeoVector& v );
Real derfct1_3_P2tilde_3D( const GeoVector& v );

Real derfct2_1_P2tilde_3D( const GeoVector& v );
Real derfct2_2_P2tilde_3D( const GeoVector& v );
Real derfct2_3_P2tilde_3D( const GeoVector& v );

Real derfct3_1_P2tilde_3D( const GeoVector& v );
Real derfct3_2_P2tilde_3D( const GeoVector& v );
Real derfct3_3_P2tilde_3D( const GeoVector& v );

Real derfct4_1_P2tilde_3D( const GeoVector& v );
Real derfct4_2_P2tilde_3D( const GeoVector& v );
Real derfct4_3_P2tilde_3D( const GeoVector& v );

Real derfct5_1_P2tilde_3D( const GeoVector& v );
Real derfct5_2_P2tilde_3D( const GeoVector& v );
Real derfct5_3_P2tilde_3D( const GeoVector& v );

Real derfct6_1_P2tilde_3D( const GeoVector& v );
Real derfct6_2_P2tilde_3D( const GeoVector& v );
Real derfct6_3_P2tilde_3D( const GeoVector& v );

Real derfct7_1_P2tilde_3D( const GeoVector& v );
Real derfct7_2_P2tilde_3D( const GeoVector& v );
Real derfct7_3_P2tilde_3D( const GeoVector& v );

Real derfct8_1_P2tilde_3D( const GeoVector& v );
Real derfct8_2_P2tilde_3D( const GeoVector& v );
Real derfct8_3_P2tilde_3D( const GeoVector& v );

Real derfct9_1_P2tilde_3D( const GeoVector& v );
Real derfct9_2_P2tilde_3D( const GeoVector& v );
Real derfct9_3_P2tilde_3D( const GeoVector& v );

Real derfct10_1_P2tilde_3D( const GeoVector& v );
Real derfct10_2_P2tilde_3D( const GeoVector& v );
Real derfct10_3_P2tilde_3D( const GeoVector& v );

Real derfct11_1_P2tilde_3D( const GeoVector& v );
Real derfct11_2_P2tilde_3D( const GeoVector& v );
Real derfct11_3_P2tilde_3D( const GeoVector& v );


Real der2fct1_11_P2tilde_3D( const GeoVector& v );
Real der2fct1_12_P2tilde_3D( const GeoVector& v );
Real der2fct1_13_P2tilde_3D( const GeoVector& v );
Real der2fct1_21_P2tilde_3D( const GeoVector& v );
Real der2fct1_22_P2tilde_3D( const GeoVector& v );
Real der2fct1_23_P2tilde_3D( const GeoVector& v );
Real der2fct1_31_P2tilde_3D( const GeoVector& v );
Real der2fct1_32_P2tilde_3D( const GeoVector& v );
Real der2fct1_33_P2tilde_3D( const GeoVector& v );

Real der2fct2_11_P2tilde_3D( const GeoVector& v );
Real der2fct2_12_P2tilde_3D( const GeoVector& v );
Real der2fct2_13_P2tilde_3D( const GeoVector& v );
Real der2fct2_21_P2tilde_3D( const GeoVector& v );
Real der2fct2_22_P2tilde_3D( const GeoVector& v );
Real der2fct2_23_P2tilde_3D( const GeoVector& v );
Real der2fct2_31_P2tilde_3D( const GeoVector& v );
Real der2fct2_32_P2tilde_3D( const GeoVector& v );
Real der2fct2_33_P2tilde_3D( const GeoVector& v );

Real der2fct3_11_P2tilde_3D( const GeoVector& v );
Real der2fct3_12_P2tilde_3D( const GeoVector& v );
Real der2fct3_13_P2tilde_3D( const GeoVector& v );
Real der2fct3_21_P2tilde_3D( const GeoVector& v );
Real der2fct3_22_P2tilde_3D( const GeoVector& v );
Real der2fct3_23_P2tilde_3D( const GeoVector& v );
Real der2fct3_31_P2tilde_3D( const GeoVector& v );
Real der2fct3_32_P2tilde_3D( const GeoVector& v );
Real der2fct3_33_P2tilde_3D( const GeoVector& v );

Real der2fct4_11_P2tilde_3D( const GeoVector& v );
Real der2fct4_12_P2tilde_3D( const GeoVector& v );
Real der2fct4_13_P2tilde_3D( const GeoVector& v );
Real der2fct4_21_P2tilde_3D( const GeoVector& v );
Real der2fct4_22_P2tilde_3D( const GeoVector& v );
Real der2fct4_23_P2tilde_3D( const GeoVector& v );
Real der2fct4_31_P2tilde_3D( const GeoVector& v );
Real der2fct4_32_P2tilde_3D( const GeoVector& v );
Real der2fct4_33_P2tilde_3D( const GeoVector& v );

Real der2fct5_11_P2tilde_3D( const GeoVector& v );
Real der2fct5_12_P2tilde_3D( const GeoVector& v );
Real der2fct5_13_P2tilde_3D( const GeoVector& v );
Real der2fct5_21_P2tilde_3D( const GeoVector& v );
Real der2fct5_22_P2tilde_3D( const GeoVector& v );
Real der2fct5_23_P2tilde_3D( const GeoVector& v );
Real der2fct5_31_P2tilde_3D( const GeoVector& v );
Real der2fct5_32_P2tilde_3D( const GeoVector& v );
Real der2fct5_33_P2tilde_3D( const GeoVector& v );

Real der2fct6_11_P2tilde_3D( const GeoVector& v );
Real der2fct6_12_P2tilde_3D( const GeoVector& v );
Real der2fct6_13_P2tilde_3D( const GeoVector& v );
Real der2fct6_21_P2tilde_3D( const GeoVector& v );
Real der2fct6_22_P2tilde_3D( const GeoVector& v );
Real der2fct6_23_P2tilde_3D( const GeoVector& v );
Real der2fct6_31_P2tilde_3D( const GeoVector& v );
Real der2fct6_32_P2tilde_3D( const GeoVector& v );
Real der2fct6_33_P2tilde_3D( const GeoVector& v );

Real der2fct7_11_P2tilde_3D( const GeoVector& v );
Real der2fct7_12_P2tilde_3D( const GeoVector& v );
Real der2fct7_13_P2tilde_3D( const GeoVector& v );
Real der2fct7_21_P2tilde_3D( const GeoVector& v );
Real der2fct7_22_P2tilde_3D( const GeoVector& v );
Real der2fct7_23_P2tilde_3D( const GeoVector& v );
Real der2fct7_31_P2tilde_3D( const GeoVector& v );
Real der2fct7_32_P2tilde_3D( const GeoVector& v );
Real der2fct7_33_P2tilde_3D( const GeoVector& v );

Real der2fct8_11_P2tilde_3D( const GeoVector& v );
Real der2fct8_12_P2tilde_3D( const GeoVector& v );
Real der2fct8_13_P2tilde_3D( const GeoVector& v );
Real der2fct8_21_P2tilde_3D( const GeoVector& v );
Real der2fct8_22_P2tilde_3D( const GeoVector& v );
Real der2fct8_23_P2tilde_3D( const GeoVector& v );
Real der2fct8_31_P2tilde_3D( const GeoVector& v );
Real der2fct8_32_P2tilde_3D( const GeoVector& v );
Real der2fct8_33_P2tilde_3D( const GeoVector& v );

Real der2fct9_11_P2tilde_3D( const GeoVector& v );
Real der2fct9_12_P2tilde_3D( const GeoVector& v );
Real der2fct9_13_P2tilde_3D( const GeoVector& v );
Real der2fct9_21_P2tilde_3D( const GeoVector& v );
Real der2fct9_22_P2tilde_3D( const GeoVector& v );
Real der2fct9_23_P2tilde_3D( const GeoVector& v );
Real der2fct9_31_P2tilde_3D( const GeoVector& v );
Real der2fct9_32_P2tilde_3D( const GeoVector& v );
Real der2fct9_33_P2tilde_3D( const GeoVector& v );

Real der2fct10_11_P2tilde_3D( const GeoVector& v );
Real der2fct10_12_P2tilde_3D( const GeoVector& v );
Real der2fct10_13_P2tilde_3D( const GeoVector& v );
Real der2fct10_21_P2tilde_3D( const GeoVector& v );
Real der2fct10_22_P2tilde_3D( const GeoVector& v );
Real der2fct10_23_P2tilde_3D( const GeoVector& v );
Real der2fct10_31_P2tilde_3D( const GeoVector& v );
Real der2fct10_32_P2tilde_3D( const GeoVector& v );
Real der2fct10_33_P2tilde_3D( const GeoVector& v );

Real der2fct11_11_P2tilde_3D( const GeoVector& v );
Real der2fct11_12_P2tilde_3D( const GeoVector& v );
Real der2fct11_13_P2tilde_3D( const GeoVector& v );
Real der2fct11_21_P2tilde_3D( const GeoVector& v );
Real der2fct11_22_P2tilde_3D( const GeoVector& v );
Real der2fct11_23_P2tilde_3D( const GeoVector& v );
Real der2fct11_31_P2tilde_3D( const GeoVector& v );
Real der2fct11_32_P2tilde_3D( const GeoVector& v );
Real der2fct11_33_P2tilde_3D( const GeoVector& v );


static const Real refcoor_P2tilde_3D[ 33 ] =
{
    0. , 0. , 0. ,
    1. , 0. , 0. ,
    0. , 1. , 0. ,
    0. , 0. , 1. ,
    0.5 , 0. , 0. ,
    0.5 , 0.5 , 0. ,
    0. , 0.5 , 0. ,
    0. , 0. , 0.5,
    0.5 , 0. , 0.5,
    0. , 0.5 , 0.5,
    0.25, 0.25, 0.25
};

static const RefEle::function_Type fct_P2tilde_3D[ 11 ] =
{
    fct1_P2tilde_3D, fct2_P2tilde_3D, fct3_P2tilde_3D, fct4_P2tilde_3D,
    fct5_P2tilde_3D, fct6_P2tilde_3D, fct7_P2tilde_3D, fct8_P2tilde_3D,
    fct9_P2tilde_3D, fct10_P2tilde_3D, fct11_P2tilde_3D
};

static const RefEle::function_Type derfct_P2tilde_3D[ 33 ] =
{
    derfct1_1_P2tilde_3D, derfct1_2_P2tilde_3D, derfct1_3_P2tilde_3D,
    derfct2_1_P2tilde_3D, derfct2_2_P2tilde_3D, derfct2_3_P2tilde_3D,
    derfct3_1_P2tilde_3D, derfct3_2_P2tilde_3D, derfct3_3_P2tilde_3D,
    derfct4_1_P2tilde_3D, derfct4_2_P2tilde_3D, derfct4_3_P2tilde_3D,
    derfct5_1_P2tilde_3D, derfct5_2_P2tilde_3D, derfct5_3_P2tilde_3D,
    derfct6_1_P2tilde_3D, derfct6_2_P2tilde_3D, derfct6_3_P2tilde_3D,
    derfct7_1_P2tilde_3D, derfct7_2_P2tilde_3D, derfct7_3_P2tilde_3D,
    derfct8_1_P2tilde_3D, derfct8_2_P2tilde_3D, derfct8_3_P2tilde_3D,
    derfct9_1_P2tilde_3D, derfct9_2_P2tilde_3D, derfct9_3_P2tilde_3D,
    derfct10_1_P2tilde_3D, derfct10_2_P2tilde_3D, derfct10_3_P2tilde_3D,
    derfct11_1_P2tilde_3D, derfct11_2_P2tilde_3D, derfct11_3_P2tilde_3D
};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2tilde\_3D, der2fct$i\_$j"."2\_P2tilde\_3D, der2fct$i\_$j"."3_P2tilde\_3D,\n";
  }
}
*/
static const RefEle::function_Type der2fct_P2tilde_3D[ 99 ] =
{
    der2fct1_11_P2tilde_3D, der2fct1_12_P2tilde_3D, der2fct1_13_P2tilde_3D,
    der2fct1_21_P2tilde_3D, der2fct1_22_P2tilde_3D, der2fct1_23_P2tilde_3D,
    der2fct1_31_P2tilde_3D, der2fct1_32_P2tilde_3D, der2fct1_33_P2tilde_3D,
    der2fct2_11_P2tilde_3D, der2fct2_12_P2tilde_3D, der2fct2_13_P2tilde_3D,
    der2fct2_21_P2tilde_3D, der2fct2_22_P2tilde_3D, der2fct2_23_P2tilde_3D,
    der2fct2_31_P2tilde_3D, der2fct2_32_P2tilde_3D, der2fct2_33_P2tilde_3D,
    der2fct3_11_P2tilde_3D, der2fct3_12_P2tilde_3D, der2fct3_13_P2tilde_3D,
    der2fct3_21_P2tilde_3D, der2fct3_22_P2tilde_3D, der2fct3_23_P2tilde_3D,
    der2fct3_31_P2tilde_3D, der2fct3_32_P2tilde_3D, der2fct3_33_P2tilde_3D,
    der2fct4_11_P2tilde_3D, der2fct4_12_P2tilde_3D, der2fct4_13_P2tilde_3D,
    der2fct4_21_P2tilde_3D, der2fct4_22_P2tilde_3D, der2fct4_23_P2tilde_3D,
    der2fct4_31_P2tilde_3D, der2fct4_32_P2tilde_3D, der2fct4_33_P2tilde_3D,
    der2fct5_11_P2tilde_3D, der2fct5_12_P2tilde_3D, der2fct5_13_P2tilde_3D,
    der2fct5_21_P2tilde_3D, der2fct5_22_P2tilde_3D, der2fct5_23_P2tilde_3D,
    der2fct5_31_P2tilde_3D, der2fct5_32_P2tilde_3D, der2fct5_33_P2tilde_3D,
    der2fct6_11_P2tilde_3D, der2fct6_12_P2tilde_3D, der2fct6_13_P2tilde_3D,
    der2fct6_21_P2tilde_3D, der2fct6_22_P2tilde_3D, der2fct6_23_P2tilde_3D,
    der2fct6_31_P2tilde_3D, der2fct6_32_P2tilde_3D, der2fct6_33_P2tilde_3D,
    der2fct7_11_P2tilde_3D, der2fct7_12_P2tilde_3D, der2fct7_13_P2tilde_3D,
    der2fct7_21_P2tilde_3D, der2fct7_22_P2tilde_3D, der2fct7_23_P2tilde_3D,
    der2fct7_31_P2tilde_3D, der2fct7_32_P2tilde_3D, der2fct7_33_P2tilde_3D,
    der2fct8_11_P2tilde_3D, der2fct8_12_P2tilde_3D, der2fct8_13_P2tilde_3D,
    der2fct8_21_P2tilde_3D, der2fct8_22_P2tilde_3D, der2fct8_23_P2tilde_3D,
    der2fct8_31_P2tilde_3D, der2fct8_32_P2tilde_3D, der2fct8_33_P2tilde_3D,
    der2fct9_11_P2tilde_3D, der2fct9_12_P2tilde_3D, der2fct9_13_P2tilde_3D,
    der2fct9_21_P2tilde_3D, der2fct9_22_P2tilde_3D, der2fct9_23_P2tilde_3D,
    der2fct9_31_P2tilde_3D, der2fct9_32_P2tilde_3D, der2fct9_33_P2tilde_3D,
    der2fct10_11_P2tilde_3D, der2fct10_12_P2tilde_3D, der2fct10_13_P2tilde_3D,
    der2fct10_21_P2tilde_3D, der2fct10_22_P2tilde_3D, der2fct10_23_P2tilde_3D,
    der2fct10_31_P2tilde_3D, der2fct10_32_P2tilde_3D, der2fct10_33_P2tilde_3D,
    der2fct11_11_P2tilde_3D, der2fct11_12_P2tilde_3D, der2fct11_13_P2tilde_3D,
    der2fct11_21_P2tilde_3D, der2fct11_22_P2tilde_3D, der2fct11_23_P2tilde_3D,
    der2fct11_31_P2tilde_3D, der2fct11_32_P2tilde_3D, der2fct11_33_P2tilde_3D
};

//======================================================================
//
//                            Q0  (3D)
//
//======================================================================
/*
                      ________
                     /.      /|
      / .     / |
     /_______/  |
     |  .  1 |  |
     |  .....|..|
     | .     | /
     |.      |/
     |_______|

*/
Real fct1_Q0_3D( const GeoVector& v );
Real derfct1_Q0_3D( const GeoVector& v );
// The second derivative is equal to the first : both = 0.
Real der2fct1_Q0_3D( const GeoVector& v );

static const Real refcoor_Q0_3D[ 3 ] =
{
    0.5 , 0.5 , 0.5
};


static const RefEle::function_Type fct_Q0_3D[ 1 ] =
{
    fct1_Q0_3D
};

static const RefEle::function_Type derfct_Q0_3D[ 3 ] =
{
    derfct1_Q0_3D, derfct1_Q0_3D, derfct1_Q0_3D
};

static const RefEle::function_Type der2fct_Q0_3D[ 9 ] =
{
    der2fct1_Q0_3D, der2fct1_Q0_3D, der2fct1_Q0_3D,
    der2fct1_Q0_3D, der2fct1_Q0_3D, der2fct1_Q0_3D,
    der2fct1_Q0_3D, der2fct1_Q0_3D, der2fct1_Q0_3D
};

//======================================================================
//
//                            Q1  (3D)
//
//======================================================================
/*
                      8-------7
                     /.      /|
      / .     / |
     5_______6  |
     |  .    |  |
     |  4....|..3
     | .     | /
     |.      |/
     1_______2
*/
Real fct1_Q1_3D( const GeoVector& v );
Real fct2_Q1_3D( const GeoVector& v );
Real fct3_Q1_3D( const GeoVector& v );
Real fct4_Q1_3D( const GeoVector& v );
Real fct5_Q1_3D( const GeoVector& v );
Real fct6_Q1_3D( const GeoVector& v );
Real fct7_Q1_3D( const GeoVector& v );
Real fct8_Q1_3D( const GeoVector& v );

Real derfct1_1_Q1_3D( const GeoVector& v );
Real derfct1_2_Q1_3D( const GeoVector& v );
Real derfct1_3_Q1_3D( const GeoVector& v );
Real derfct2_1_Q1_3D( const GeoVector& v );
Real derfct2_2_Q1_3D( const GeoVector& v );
Real derfct2_3_Q1_3D( const GeoVector& v );
Real derfct3_1_Q1_3D( const GeoVector& v );
Real derfct3_2_Q1_3D( const GeoVector& v );
Real derfct3_3_Q1_3D( const GeoVector& v );
Real derfct4_1_Q1_3D( const GeoVector& v );
Real derfct4_2_Q1_3D( const GeoVector& v );
Real derfct4_3_Q1_3D( const GeoVector& v );
Real derfct5_1_Q1_3D( const GeoVector& v );
Real derfct5_2_Q1_3D( const GeoVector& v );
Real derfct5_3_Q1_3D( const GeoVector& v );
Real derfct6_1_Q1_3D( const GeoVector& v );
Real derfct6_2_Q1_3D( const GeoVector& v );
Real derfct6_3_Q1_3D( const GeoVector& v );
Real derfct7_1_Q1_3D( const GeoVector& v );
Real derfct7_2_Q1_3D( const GeoVector& v );
Real derfct7_3_Q1_3D( const GeoVector& v );
Real derfct8_1_Q1_3D( const GeoVector& v );
Real derfct8_2_Q1_3D( const GeoVector& v );
Real derfct8_3_Q1_3D( const GeoVector& v );

Real der2fct1_11_Q1_3D( const GeoVector& v );
Real der2fct1_12_Q1_3D( const GeoVector& v );
Real der2fct1_13_Q1_3D( const GeoVector& v );
Real der2fct1_21_Q1_3D( const GeoVector& v );
Real der2fct1_22_Q1_3D( const GeoVector& v );
Real der2fct1_23_Q1_3D( const GeoVector& v );
Real der2fct1_31_Q1_3D( const GeoVector& v );
Real der2fct1_32_Q1_3D( const GeoVector& v );
Real der2fct1_33_Q1_3D( const GeoVector& v );

Real der2fct2_11_Q1_3D( const GeoVector& v );
Real der2fct2_12_Q1_3D( const GeoVector& v );
Real der2fct2_13_Q1_3D( const GeoVector& v );
Real der2fct2_21_Q1_3D( const GeoVector& v );
Real der2fct2_22_Q1_3D( const GeoVector& v );
Real der2fct2_23_Q1_3D( const GeoVector& v );
Real der2fct2_31_Q1_3D( const GeoVector& v );
Real der2fct2_32_Q1_3D( const GeoVector& v );
Real der2fct2_33_Q1_3D( const GeoVector& v );

Real der2fct3_11_Q1_3D( const GeoVector& v );
Real der2fct3_12_Q1_3D( const GeoVector& v );
Real der2fct3_13_Q1_3D( const GeoVector& v );
Real der2fct3_21_Q1_3D( const GeoVector& v );
Real der2fct3_22_Q1_3D( const GeoVector& v );
Real der2fct3_23_Q1_3D( const GeoVector& v );
Real der2fct3_31_Q1_3D( const GeoVector& v );
Real der2fct3_32_Q1_3D( const GeoVector& v );
Real der2fct3_33_Q1_3D( const GeoVector& v );

Real der2fct4_11_Q1_3D( const GeoVector& v );
Real der2fct4_12_Q1_3D( const GeoVector& v );
Real der2fct4_13_Q1_3D( const GeoVector& v );
Real der2fct4_21_Q1_3D( const GeoVector& v );
Real der2fct4_22_Q1_3D( const GeoVector& v );
Real der2fct4_23_Q1_3D( const GeoVector& v );
Real der2fct4_31_Q1_3D( const GeoVector& v );
Real der2fct4_32_Q1_3D( const GeoVector& v );
Real der2fct4_33_Q1_3D( const GeoVector& v );

Real der2fct5_11_Q1_3D( const GeoVector& v );
Real der2fct5_12_Q1_3D( const GeoVector& v );
Real der2fct5_13_Q1_3D( const GeoVector& v );
Real der2fct5_21_Q1_3D( const GeoVector& v );
Real der2fct5_22_Q1_3D( const GeoVector& v );
Real der2fct5_23_Q1_3D( const GeoVector& v );
Real der2fct5_31_Q1_3D( const GeoVector& v );
Real der2fct5_32_Q1_3D( const GeoVector& v );
Real der2fct5_33_Q1_3D( const GeoVector& v );

Real der2fct6_11_Q1_3D( const GeoVector& v );
Real der2fct6_12_Q1_3D( const GeoVector& v );
Real der2fct6_13_Q1_3D( const GeoVector& v );
Real der2fct6_21_Q1_3D( const GeoVector& v );
Real der2fct6_22_Q1_3D( const GeoVector& v );
Real der2fct6_23_Q1_3D( const GeoVector& v );
Real der2fct6_31_Q1_3D( const GeoVector& v );
Real der2fct6_32_Q1_3D( const GeoVector& v );
Real der2fct6_33_Q1_3D( const GeoVector& v );

Real der2fct7_11_Q1_3D( const GeoVector& v );
Real der2fct7_12_Q1_3D( const GeoVector& v );
Real der2fct7_13_Q1_3D( const GeoVector& v );
Real der2fct7_21_Q1_3D( const GeoVector& v );
Real der2fct7_22_Q1_3D( const GeoVector& v );
Real der2fct7_23_Q1_3D( const GeoVector& v );
Real der2fct7_31_Q1_3D( const GeoVector& v );
Real der2fct7_32_Q1_3D( const GeoVector& v );
Real der2fct7_33_Q1_3D( const GeoVector& v );

Real der2fct8_11_Q1_3D( const GeoVector& v );
Real der2fct8_12_Q1_3D( const GeoVector& v );
Real der2fct8_13_Q1_3D( const GeoVector& v );
Real der2fct8_21_Q1_3D( const GeoVector& v );
Real der2fct8_22_Q1_3D( const GeoVector& v );
Real der2fct8_23_Q1_3D( const GeoVector& v );
Real der2fct8_31_Q1_3D( const GeoVector& v );
Real der2fct8_32_Q1_3D( const GeoVector& v );
Real der2fct8_33_Q1_3D( const GeoVector& v );

static const Real refcoor_Q1_3D[ 24 ] =
{
    0. , 0. , 0. ,
    1. , 0. , 0. ,
    1. , 1. , 0. ,
    0. , 1. , 0. ,
    0. , 0. , 1. ,
    1. , 0. , 1. ,
    1. , 1. , 1. ,
    0. , 1. , 1.
};


static const RefEle::function_Type fct_Q1_3D[ 8 ] =
{
    fct1_Q1_3D, fct2_Q1_3D, fct3_Q1_3D, fct4_Q1_3D, fct5_Q1_3D,
    fct6_Q1_3D, fct7_Q1_3D, fct8_Q1_3D
};


static const RefEle::function_Type derfct_Q1_3D[ 24 ] =
{
    derfct1_1_Q1_3D, derfct1_2_Q1_3D, derfct1_3_Q1_3D,
    derfct2_1_Q1_3D, derfct2_2_Q1_3D, derfct2_3_Q1_3D,
    derfct3_1_Q1_3D, derfct3_2_Q1_3D, derfct3_3_Q1_3D,
    derfct4_1_Q1_3D, derfct4_2_Q1_3D, derfct4_3_Q1_3D,
    derfct5_1_Q1_3D, derfct5_2_Q1_3D, derfct5_3_Q1_3D,
    derfct6_1_Q1_3D, derfct6_2_Q1_3D, derfct6_3_Q1_3D,
    derfct7_1_Q1_3D, derfct7_2_Q1_3D, derfct7_3_Q1_3D,
    derfct8_1_Q1_3D, derfct8_2_Q1_3D, derfct8_3_Q1_3D
};
/*the perl-script:
  #!/usr/bin/perl
  for($i=1;$i<=10;$i++){
  for($j=1;$j<=3;$j++){
  printf "der2fct$i\_$j"."1\_P2\_3D, der2fct$i\_$j"."2\_P2\_3D, der2fct$i\_$j"."3_P2\_3D,\n";
  }
}
*/
static const RefEle::function_Type der2fct_Q1_3D[ 72 ] =
{
    der2fct1_11_Q1_3D, der2fct1_12_Q1_3D, der2fct1_13_Q1_3D,
    der2fct1_21_Q1_3D, der2fct1_22_Q1_3D, der2fct1_23_Q1_3D,
    der2fct1_31_Q1_3D, der2fct1_32_Q1_3D, der2fct1_33_Q1_3D,
    der2fct2_11_Q1_3D, der2fct2_12_Q1_3D, der2fct2_13_Q1_3D,
    der2fct2_21_Q1_3D, der2fct2_22_Q1_3D, der2fct2_23_Q1_3D,
    der2fct2_31_Q1_3D, der2fct2_32_Q1_3D, der2fct2_33_Q1_3D,
    der2fct3_11_Q1_3D, der2fct3_12_Q1_3D, der2fct3_13_Q1_3D,
    der2fct3_21_Q1_3D, der2fct3_22_Q1_3D, der2fct3_23_Q1_3D,
    der2fct3_31_Q1_3D, der2fct3_32_Q1_3D, der2fct3_33_Q1_3D,
    der2fct4_11_Q1_3D, der2fct4_12_Q1_3D, der2fct4_13_Q1_3D,
    der2fct4_21_Q1_3D, der2fct4_22_Q1_3D, der2fct4_23_Q1_3D,
    der2fct4_31_Q1_3D, der2fct4_32_Q1_3D, der2fct4_33_Q1_3D,
    der2fct5_11_Q1_3D, der2fct5_12_Q1_3D, der2fct5_13_Q1_3D,
    der2fct5_21_Q1_3D, der2fct5_22_Q1_3D, der2fct5_23_Q1_3D,
    der2fct5_31_Q1_3D, der2fct5_32_Q1_3D, der2fct5_33_Q1_3D,
    der2fct6_11_Q1_3D, der2fct6_12_Q1_3D, der2fct6_13_Q1_3D,
    der2fct6_21_Q1_3D, der2fct6_22_Q1_3D, der2fct6_23_Q1_3D,
    der2fct6_31_Q1_3D, der2fct6_32_Q1_3D, der2fct6_33_Q1_3D,
    der2fct7_11_Q1_3D, der2fct7_12_Q1_3D, der2fct7_13_Q1_3D,
    der2fct7_21_Q1_3D, der2fct7_22_Q1_3D, der2fct7_23_Q1_3D,
    der2fct7_31_Q1_3D, der2fct7_32_Q1_3D, der2fct7_33_Q1_3D,
    der2fct8_11_Q1_3D, der2fct8_12_Q1_3D, der2fct8_13_Q1_3D,
    der2fct8_21_Q1_3D, der2fct8_22_Q1_3D, der2fct8_23_Q1_3D,
    der2fct8_31_Q1_3D, der2fct8_32_Q1_3D, der2fct8_33_Q1_3D
};


//!======================================================================
//!
//!                           RT0  (3D)
//!
//!======================================================================
/*!

       8-------7
      /.      /|
     / .     / |
    5_______6  |
    |  .    |  |
    |  4....|..3
    | .     | /
    |.      |/
    1_______2

face 1: 1,4,3,2
face 2: 1,5,8,4
face 3: 1,2,6,5
face 4: 2,3,7,6
face 5: 3,4,8,7
face 6: 5,6,7,8

*/

Real fct1_RT0_1_HEXA_3D( const GeoVector& v );
Real fct1_RT0_2_HEXA_3D( const GeoVector& v );
Real fct1_RT0_3_HEXA_3D( const GeoVector& v );

Real fct2_RT0_1_HEXA_3D( const GeoVector& v );
Real fct2_RT0_2_HEXA_3D( const GeoVector& v );
Real fct2_RT0_3_HEXA_3D( const GeoVector& v );

Real fct3_RT0_1_HEXA_3D( const GeoVector& v );
Real fct3_RT0_2_HEXA_3D( const GeoVector& v );
Real fct3_RT0_3_HEXA_3D( const GeoVector& v );

Real fct4_RT0_1_HEXA_3D( const GeoVector& v );
Real fct4_RT0_2_HEXA_3D( const GeoVector& v );
Real fct4_RT0_3_HEXA_3D( const GeoVector& v );

Real fct5_RT0_1_HEXA_3D( const GeoVector& v );
Real fct5_RT0_2_HEXA_3D( const GeoVector& v );
Real fct5_RT0_3_HEXA_3D( const GeoVector& v );

Real fct6_RT0_1_HEXA_3D( const GeoVector& v );
Real fct6_RT0_2_HEXA_3D( const GeoVector& v );
Real fct6_RT0_3_HEXA_3D( const GeoVector& v );

Real fct1_DIV_RT0_HEXA_3D( const GeoVector& v );
Real fct2_DIV_RT0_HEXA_3D( const GeoVector& v );
Real fct3_DIV_RT0_HEXA_3D( const GeoVector& v );
Real fct4_DIV_RT0_HEXA_3D( const GeoVector& v );
Real fct5_DIV_RT0_HEXA_3D( const GeoVector& v );
Real fct6_DIV_RT0_HEXA_3D( const GeoVector& v );

static const Real refcoor_RT0_HEXA_3D[ 18 ] =
{
    0.5 , 0.5 , 0.  ,
    0.  , 0.5 , 0.5 ,
    0.5 , 0.  , 0.5 ,
    1.  , 0.5 , 0.5 ,
    0.5 , 1.  , 0.5 ,
    0.5 , 0.5 , 1.
};

static const RefEle::function_Type fct_RT0_HEXA_3D[ 18 ] =
{
    fct1_RT0_1_HEXA_3D, fct1_RT0_2_HEXA_3D, fct1_RT0_3_HEXA_3D,
    fct2_RT0_1_HEXA_3D, fct2_RT0_2_HEXA_3D, fct2_RT0_3_HEXA_3D,
    fct3_RT0_1_HEXA_3D, fct3_RT0_2_HEXA_3D, fct3_RT0_3_HEXA_3D,
    fct4_RT0_1_HEXA_3D, fct4_RT0_2_HEXA_3D, fct4_RT0_3_HEXA_3D,
    fct5_RT0_1_HEXA_3D, fct5_RT0_2_HEXA_3D, fct5_RT0_3_HEXA_3D,
    fct6_RT0_1_HEXA_3D, fct6_RT0_2_HEXA_3D, fct6_RT0_3_HEXA_3D
};

static const RefEle::function_Type fct_DIV_RT0_HEXA_3D[ 6 ] =
{
    fct1_DIV_RT0_HEXA_3D, fct2_DIV_RT0_HEXA_3D,
    fct3_DIV_RT0_HEXA_3D, fct4_DIV_RT0_HEXA_3D,
    fct5_DIV_RT0_HEXA_3D, fct6_DIV_RT0_HEXA_3D
};


//!======================================================================
//!
//!                           RT0  (3D)
//!
//!======================================================================
/*

    4
    | .
    |  \.3
    |  . \\
    | .    \\
    |.       \
    1 ---------2

SEE basisElSh.cc   for the ORIENTATION CONVENTIONS
point 1: 0, 0, 0
point 2: 1, 0, 0
point 3: 0, 1, 0
point 4: 0, 0, 1

face 1: 2, 3, 4
face 2: 1, 4, 3
face 3: 1, 2, 4
face 4: 1, 3, 2
*/

Real fct1_RT0_1_TETRA_3D( const GeoVector& v );
Real fct1_RT0_2_TETRA_3D( const GeoVector& v );
Real fct1_RT0_3_TETRA_3D( const GeoVector& v );

Real fct2_RT0_1_TETRA_3D( const GeoVector& v );
Real fct2_RT0_2_TETRA_3D( const GeoVector& v );
Real fct2_RT0_3_TETRA_3D( const GeoVector& v );

Real fct3_RT0_1_TETRA_3D( const GeoVector& v );
Real fct3_RT0_2_TETRA_3D( const GeoVector& v );
Real fct3_RT0_3_TETRA_3D( const GeoVector& v );

Real fct4_RT0_1_TETRA_3D( const GeoVector& v );
Real fct4_RT0_2_TETRA_3D( const GeoVector& v );
Real fct4_RT0_3_TETRA_3D( const GeoVector& v );


Real fct1_DIV_RT0_TETRA_3D( const GeoVector& v );
Real fct2_DIV_RT0_TETRA_3D( const GeoVector& v );
Real fct3_DIV_RT0_TETRA_3D( const GeoVector& v );
Real fct4_DIV_RT0_TETRA_3D( const GeoVector& v );


static const Real refcoor_RT0_TETRA_3D[ 12 ] =
{
    1. / 3  , 1. / 3. , 0.      ,
    1. / 3. , 0.      , 1. / 3. ,
    1. / 3. , 1. / 3. , 1. / 3. ,
    0.      , 1. / 3. , 1. / 3.
};

static const RefEle::function_Type fct_RT0_TETRA_3D[ 12 ] =
{
    fct1_RT0_1_TETRA_3D, fct1_RT0_2_TETRA_3D, fct1_RT0_3_TETRA_3D,
    fct2_RT0_1_TETRA_3D, fct2_RT0_2_TETRA_3D, fct2_RT0_3_TETRA_3D,
    fct3_RT0_1_TETRA_3D, fct3_RT0_2_TETRA_3D, fct3_RT0_3_TETRA_3D,
    fct4_RT0_1_TETRA_3D, fct4_RT0_2_TETRA_3D, fct4_RT0_3_TETRA_3D
};

static const RefEle::function_Type fct_DIV_RT0_TETRA_3D[ 4 ] =
{
    fct1_DIV_RT0_TETRA_3D, fct2_DIV_RT0_TETRA_3D,
    fct3_DIV_RT0_TETRA_3D, fct4_DIV_RT0_TETRA_3D
};

}
#endif
