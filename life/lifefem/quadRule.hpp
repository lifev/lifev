/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
#include "lifeV.hpp"
#include "basisElSh.hpp"

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
    Real _coor[3];
    Real _weight;
public:
    QuadPoint():_weight(0){
	_coor[0] = 0;_coor[1]=0;_coor[2]=0;
    }
    QuadPoint(Real x,Real y,Real z,Real weight):_weight(weight){
	_coor[0] = x;_coor[1] = y;_coor[2] = z;
    }
    QuadPoint(Real x,Real y,Real weight):_weight(weight){
	_coor[0] = x;_coor[1] = y;_coor[2] = 0.;
    }
    QuadPoint(Real x,Real weight):_weight(weight){
	_coor[0] = x;_coor[1] = 0.;_coor[2] = 0.;
    }
    inline Real weight() const {return _weight;}
    inline Real x() const {return _coor[0];}
    inline Real y() const {return _coor[1];}
    inline Real z() const {return _coor[2];}
    inline Real coor(int i) const{ASSERT_BD(i>=0 && i<3);return _coor[i];}
    friend std::ostream & operator << ( std::ostream& c, const QuadPoint& pt)
	{
	    c << "coor = " << pt._coor[0] << "," << pt._coor[1] << "," << pt._coor[2]
	      << "  weight = " << pt._weight ;
	    return c;
	}
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
  quadrature rule. Finally, increment NB_QUAD_RULE_TETRA,
  and add it to the array of all quadrature rule on a tetra, namely :

  \code
  static const QuadRule quad_rule_tetra[NB_QUAD_RULE_TETRA] = {
  quadRuleTetra1pt,
  quadRuleTetra4pt,
  quadRuleTetra5pt,
  quadRulePipo
};
  \endcode
*/
class QuadRule
{
  const  QuadPoint* _pt;
public:
  QuadRule(const QuadPoint* pt,int _id, std::string _name,ReferenceShapes _shape,
	   int _nbQuadPt,int _defOfExact);
  ~QuadRule();
  const ReferenceShapes shape;//!< geometrical shape of the domain on which the quadrature rule can be used
  const int id;//!< id of the quadrature rule. e.g. QR_GAUSS_LEGENDRE_1PT_1D
  const std::string name; //!< name of the quadrature rule
  const int nbQuadPt; //!< number of quadrature points
  const int degOfExact; //!< degree of exactness
  const QuadPoint& quadPoint(int ig) const//!< quadPoint(ig) is the ig-th quadrature point
  {ASSERT_BD(ig<nbQuadPt); return _pt[ig];}
  Real weight(int ig) const//!< weight(ig) is the ig-th quadrature weight
  {ASSERT_BD(ig<nbQuadPt); return _pt[ig].weight();}
  //! quadPointCoor(ig,icoor) is the coordinate icoor of the quadrature point ig
  inline const Real quadPointCoor(int ig,int icoor) const{
    ASSERT_BD(ig<nbQuadPt);
    return _pt[ig].coor(icoor);
  }
  friend std::ostream& operator << (std::ostream& c,const QuadRule& qr);
};

/*!
  \class SetOfQuadRule
  \author J.-F. Gerbeau
  \date 04/2002
  \brief Set of quadrature rule

  The purpose of the class SetOfQuadRule is to gather all the quadrature
  rules defined on a geometry (e.g. triangle, tetra, etc) in order to
  easily fill the array of a reference finite element with all avaible
  quadrature rule. This allows to switch very easily between various quadrature rules
  during the use of the finite element.
*/
class SetOfQuadRule
{
  const QuadRule* _qr; //!< array of quadrature rules
  int _totalNbQuadPoint;
  int _maxIdQuadRule;
  int* _posQuadRule;
public:
  const int nbQuadRule;
  SetOfQuadRule(const QuadRule* qr,int _nb);
  ~SetOfQuadRule();
  //! number of quadrature rules in the set
  int totalNbQuadPoint() const {return _totalNbQuadPoint;};
  //! the maximum of the id of quadrature rules in the set
  int maxIdQuadRule() const {return _maxIdQuadRule;};
  //! posQuadRule(t) = position in the set of the quadrature of id t
  const int posQuadRule(int t) const{
    ASSERT_BD(t>=0 && t<_maxIdQuadRule)
      return _posQuadRule[t];
  }
  //! return the k-th quadrature rule in the set
  const QuadRule& quadRule(int k) const
  {
    ASSERT_BD(k<nbQuadRule)
      return _qr[k];
  }
};

//======================================================================
extern const SetOfQuadRule allQuadRuleSeg;
extern const SetOfQuadRule allQuadRuleTria;
extern const SetOfQuadRule allQuadRuleQuad;
extern const SetOfQuadRule allQuadRuleHexa;
extern const SetOfQuadRule allQuadRuleTetra;
//======================================================================
extern const QuadRule quadRuleDummy;
extern const QuadRule quadRuleSeg1pt;
extern const QuadRule quadRuleSeg2pt;
extern const QuadRule quadRuleSeg3pt;

extern const QuadRule quadRuleTria1pt;
extern const QuadRule quadRuleTria3pt;
extern const QuadRule quadRuleTria4pt;

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
