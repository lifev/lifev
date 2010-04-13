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
#ifndef QUADRULE_H
#define QUADRULE_H
#include <life/lifecore/life.hpp>
#include <life/lifemesh/basisElSh.hpp>

namespace LifeV
{
/*----------------------------------------------------------------------
 *
 *                       Quadrature Rules
 *
 ----------------------------------------------------------------------*/

/*!
  \class QuadPoint
  \brief The class for the quadrature points i.e. (x,y,z,weight)
  \author J.-F. Gerbeau
  \date 04/2002
*/
class QuadPoint
{
public:

    //! Empty constructor (all zero data)
    QuadPoint() : M_weight( 0 )
    {
        M_coor[ 0 ] = 0;
        M_coor[ 1 ] = 0;
        M_coor[ 2 ] = 0;
    }
    
    //! Full constructor for 3D
    QuadPoint( Real x, Real y, Real z, Real weight ) : M_weight( weight )
    {
        M_coor[ 0 ] = x;
        M_coor[ 1 ] = y;
        M_coor[ 2 ] = z;
    }
    
    //! Full constructor for 2D
    QuadPoint( Real x, Real y, Real weight ) : M_weight( weight )
    {
        M_coor[ 0 ] = x;
        M_coor[ 1 ] = y;
        M_coor[ 2 ] = 0.;
    }
    
    //! Full constructor for 1D
    QuadPoint( Real x, Real weight ) : M_weight( weight )
    {
        M_coor[ 0 ] = x;
        M_coor[ 1 ] = 0.;
        M_coor[ 2 ] = 0.;
    }

    //! Copy constructor
    QuadPoint(const QuadPoint& qp) : M_weight(qp.M_weight)
    {
        M_coor[0]=qp.M_coor[0];
        M_coor[1]=qp.M_coor[1];
        M_coor[2]=qp.M_coor[2];
    };

    //! Destructor
    ~QuadPoint(){};
    
    //! Getter for the weight
    inline const Real& weight() const
    {
        return M_weight;
    }
    
    //! Getter for the first coordinate
    inline const Real& x() const
    {
        return M_coor[ 0 ];
    }

    //! Getter for the second coordinate
    inline const Real& y() const
    {
        return M_coor[ 1 ];
    }

    //! Getter for the third coordinate
    inline const Real& z() const
    {
        return M_coor[ 2 ];
    }

    //! Getter for the coordinate (0<=i<3)
    inline const Real& coor(const UInt& i ) const
    {
        ASSERT_BD(i < 3 );
        return M_coor[ i ];
    }

    //! Output operator
    friend std::ostream & operator << ( std::ostream& c, const QuadPoint& pt )
    {
        c << "coor = " << pt.M_coor[ 0 ] << "," << pt.M_coor[ 1 ] << "," << pt.M_coor[ 2 ]
        << "  weight = " << pt.M_weight ;
        return c;
    }

private:
    Real M_coor[ 3 ];
    Real M_weight;
};


/*!
  \class QuadRule
  \brief The class of quadrature rules
  \author J.-F. Gerbeau
  \date 04/2002

  It contains the quadrature points and the weights.

  \par How to add a new quadrature rule:

  Suppose you want to add the quadrature rule qr_Pipo on a tetrahedra:
  in the file quadRule.h, add
  \code
  extern const QuadRule quadRulePipo;
  \endcode
  in the file defQuadRuleFE.cc, add
  \code
  #define QUAD_RULE_TETRA_PIPO <label>
  \endcode
  where  <label> is an integer not yet attributed
  to another quadrature rule on a tetrahedra and copy/paste/modify an existing
  quadrature rule.
*/
class QuadRule
{
public:

    //! Full constructor
    QuadRule( const QuadPoint* pt, int _id, std::string _name, ReferenceShapes _shape,
              UInt _nbQuadPt, UInt _defOfExact );

    //! Copy constructor
    QuadRule( const QuadRule& qr);

    //! Destrutor
    ~QuadRule();
    
    //! quadPoint(ig) is the ig-th quadrature point
    const QuadPoint& quadPoint( const UInt& ig ) const
    {
        ASSERT_BD( ig < M_nbQuadPt );
        return M_pt[ ig ];
    }

    //! weight(ig) is the ig-th quadrature weight
    const Real& weight(const UInt& ig ) const
    {
        ASSERT_BD( ig < M_nbQuadPt );
        return M_pt[ ig ].weight();
    }

    //! quadPointCoor(ig,icoor) is the coordinate icoor of the quadrature point ig
    inline const Real& quadPointCoor( const UInt& ig, const UInt& icoor ) const
    {
        ASSERT_BD( ig < M_nbQuadPt );
        return M_pt[ ig ].coor( icoor );
    }

    //! Getter for the number of quadrature points
    inline const UInt& nbQuadPt() const
    {
        return M_nbQuadPt;
    };

    //! Output operator
    friend std::ostream& operator << ( std::ostream& c, const QuadRule& qr );

private:
    // Storage for the quadrature nodes
    const QuadPoint* M_pt;
    
    //! geometrical shape of the domain on which the quadrature rule can be used
    const ReferenceShapes M_shape;
    //! id of the quadrature rule. e.g. QR_GAUSS_LEGENDRE_1PT_1D
    const UInt M_id;
    //! name of the quadrature rule
    const std::string M_name;
    //! number of quadrature points
    const UInt M_nbQuadPt;
    //! degree of exactness
    const UInt M_degOfExact;

};


//======================================================================
extern const QuadRule quadRuleDummy;
extern const QuadRule quadRuleSeg1pt;
extern const QuadRule quadRuleSeg2pt;
extern const QuadRule quadRuleSeg3pt;

extern const QuadRule quadRuleTria1pt;
extern const QuadRule quadRuleTria3pt;
extern const QuadRule quadRuleTria4pt;
extern const QuadRule quadRuleTria6pt;
extern const QuadRule quadRuleTria7pt;

extern const QuadRule quadRuleQuad1pt;
extern const QuadRule quadRuleQuad4pt;
extern const QuadRule quadRuleQuad9pt;

extern const QuadRule quadRuleTetra1pt;
extern const QuadRule quadRuleTetra4pt;
extern const QuadRule quadRuleTetra5pt;
extern const QuadRule quadRuleTetra15pt;
extern const QuadRule quadRuleTetra64pt;

extern const QuadRule quadRuleHexa1pt;
extern const QuadRule quadRuleHexa8pt;
}
#endif
