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
    @brief File for the definition of the QuadratureRule class.

    @author Jean-Frederic Gerbeau
            Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01-06-2010

    @contributor Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @mantainer Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */


#ifndef QUADRULE_H
#define QUADRULE_H

#include <lifev/core/LifeV.hpp>

#include <lifev/core/mesh/ElementShapes.hpp>

#include <lifev/core/fem/QuadraturePoint.hpp>

#include <iostream>
#include <fstream>
#include <cstdarg>

namespace LifeV
{


//! QuadratureRule - The basis class for storing and accessing quadrature rules.
/*!

  <b> Definition of a quadrature rule </b>

  To define a quadrature rule, several constructors have been defined for this class. Besides the classical constructors, one is provided with an arbitrary number of arguements (indicated by "..."). Its usage is quite simple: first of all, one defines the standard arguements (name...) but also the number of quadrature points and their dimension. After that, one add the coordinates and weights of the quadrature points in the order: coordinates for the first point, weight for the first point, coordinates of the second point, weight of the second point,...

For example, if we want to define a Simpson rule in 1D, we can write:
 \code
  QuadratureRule simpson1D ("Simpson 1D",LINE,1,3,3,  // name,shape,dimension,exactness,number of points
                      0.0, 1.0/6.0,             // first point: x=0 with weight 1/6
                      0.5, 2.0/3.0,             // second point: x=0.5, weight 2/3
                      1.0, 1.0/6.0);            // third point: x=1, weight 1/6
  \endcode

If we want to see this quadrature as a 1D quadrature in a 2D space, the code would have been:

\code
  QuadratureRule simpson2D ("Simpson 2D",LINE,2,3,3,  // name,shape,dimension,exactness,number of points
                      0.0,0.0, 1.0/6.0,         // first point: (0,0) with weight 1/6
                      0.5,0.0, 2.0/3.0,         // second point: (0.5,0) weight 2/3
                      1.0,0.0, 1.0/6.0);        // third point: (1,0) weight 1/6
  \endcode

The following code would have produced the same quadrature:

\code
  QuadratureRule simpson2D (simpson1D,2);
  \endcode

  <b> Basic Use </b>

  A quadrature rule consists mainly in a container of quadrature points. The stored quadrature points can be accessed through accessors to the quadrature points or through specific accessors for the coordinates and the weight of the quadrature points.

  <b> Degree of exactness </b>

  The quadrature rules store the degree of exactness that they are supposed to achieve. It can be specified when building the quadrature rule. No test is performed automatically to check the order and when the quadrature rule is modified (with the addPoint method for example), the order of exactness is kept unchanged.

  However, the QuadratureRule class provides a test for the degree of exactness. When called, this test returns the degree of exactness that has been found by trying to perform several integrals, but it does not change the stored degree of exactness.

  <b> Dimension and Shape </b>

  For some problems, it makes sens to use quadrature rules outside their "original" space, e.g. one could need a quadrature rule for the triangles in a 3D space (when integrating on a surface). The QuadratureRule class provides the possibility of changing the dimension in which a quadrature rule in defined (i.e. how many coordinates are used for the coordinates of the quadrature points).

  When using a quadrature rule in a space with higher dimension, 0 coordinates are added to fill the missing coordinates. The shape of the quadrature (i.e. the geometric shape where it is defined) remains the same.

  It is currently not possible to move a quadrature to a lower space dimension (than the one defined by its shape): the problem resides in the area of the quadrature: a quadrature rule for tetrahedra has a total weight of 1/6 while a quadrature rule for triangles has a total weight of 1/2 (in general for simplexes, the area is 1/n!). If such change would be allowed, one would have to implement all the changes of weight depending on the shape (there is no such phenomena between hexahedra and quadrangles for example). Moreover, it makes very little sense to downgrade a quadrature (use a specific quadrature rule which will have very likely less quadrature nodes).

  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
  @date 1 June 2010
  @version 2.0
  @note The previous version was due to J.-F. Gerbeau (04/2002). Many ideas were taken from that version.

*/
class QuadratureRule
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Empty constructor
    QuadratureRule();

    //! Full constructor using pointers
    /*!
      With this constructor, the dimension is supposed to be 3 (old style).

      @param pt The set of quadrature points
      @param id Backward compatibility arguement
      @param name The name of the quadrature rule
      @param shape The shape of the element to be used with
      @param nbQuadPt The number of quadrature points defined
      @param degOfExact The degree of exactness of the quadrature rule
     */
    QuadratureRule ( const QuadraturePoint* pt, int id, std::string name, ReferenceShapes shape,
                     UInt nbQuadPt, UInt degOfExact );

    //! Full constructor
    /*!
      This constructor enables to use as many arguments as needed for the declaration of the
      quadrature rules (the number of quadrature points is not known a priori).

      The last arguements are the coordinates and the weights of the points, in the order:
      coordinates of the first point, weight of the first point, coordinates of the second
      point, weight of the second point,...

      @param name The name of the quadrature
      @param shape The shape were the quadrature is originally defined
      @param dimension The dimension of the space in which the quadrature is defined
      @param degreeOfExactness The degree of exactness of the quadrature
      @param nbQuadPt The number of quadrature points
     */
    QuadratureRule (std::string name, ReferenceShapes shape, UInt dimension, UInt degreeOfExactness, UInt nbQuadPt, ... );


    //! Copy constructor
    /*!
      @param qr The quadrature rule that we want to copy.
     */
    QuadratureRule ( const QuadratureRule& qr);

    //! Copy constructor using a different dimension
    /*!
      This can be used to export a quadrature rule from a space dimension
      to another one.

      @param qr The quadrature rule to export
      @param dim The new dimension to be used with the quadrature rule.
     */
    QuadratureRule ( const QuadratureRule& qr, const UInt dim);


    //! Destructor
    virtual ~QuadratureRule();

    //@}


    //! @name Operators
    //@{

    //! Output operator
    friend std::ostream& operator << ( std::ostream& c, const QuadratureRule& qr );

    //@}


    //! @name Methods
    //@{

    //! Method to add a point to an existing quadrature rule.
    /*!
      Beware to have set the dimension and the shape before calling
      this method.

      @note: the degree of exactness is not changed, it is up to the user
      to take care of it.

      @param qp The quadrature point to add.
     */
    void addPoint (const QuadraturePoint& qp)
    {
        M_pt.push_back (QuadraturePoint (qp, M_dimension) );
        M_nbQuadPt += 1;
    };

    //! ShowMe method
    void showMe ( std::ostream& output = std::cout) const;

    //! Check for the exactness of the quadrature
    /*!
      The quadrature rule is used to performed intergrals whose
      values are known. The degree of exactness is evaluated
      using these results.

      @note The degree of exactness stored internally is not
      changed when calling this method.
     */
    UInt checkExactness() const;

    //! VTK export for the quadrature
    /*!
      Creates a file with the positions of the
      quadrature points in a VTK format. No extension
      is added to the name of the file, please provide a name
      ending by ".vtk".

      @param filename The name of the file to be created.
     */
    void vtkExport ( const std::string& filename) const;

    //! Method for importing the quadrature rule from another class
    template <typename QRType>
    void import ( const QRType& qr);

    //@}


    //! @name Set methods
    //@{

    //! Change the quadrature points for the ones given here.
    /*!
      Beware to have set the dimension (default: 0!) before calling
      this method (use the QuadratureRule::setDimensionShape method for example).
     */
    void setPoints (const std::vector<QuadraturePoint>& pts);

    //! Change the quadrature points for the one given here
    /*!
      Use coordinates and weights separetly. The two vectors
      given in argument must have the same length. The quadrature
      points are then given by (coorindates[i],weights[i]).

      Beware to have set the dimension (default: 0!) before calling
      this method (use the QuadratureRule::setDimensionShape method for example).

     @param coordinates An array containing the coordinates of the points
     @param weights An array containing the weights of the points
     */
    void setPoints (const std::vector<GeoVector>& coordinates, const std::vector<Real>& weights);

    //! Change the name of the quadrature
    void setName (const std::string& newName);

    //! Change the degree of exactness
    void setExactness (const UInt& exactness);

    //! Change the dimension and the shape
    void setDimensionShape (const UInt& newDim, const ReferenceShapes& newShape);


    //@}


    //! @name Get methods
    //@{

    //! quadPoint(ig) is the ig-th quadrature point
    const QuadraturePoint& quadPoint ( const UInt& ig ) const
    {
        ASSERT_BD ( ig < M_nbQuadPt );
        return M_pt[ ig ];
    }

    //! weight(ig) is the ig-th quadrature weight
    const Real& weight (const UInt& ig ) const
    {
        ASSERT_BD ( ig < M_nbQuadPt );
        return M_pt[ ig ].weight();
    }

    //! quadPointCoor(ig,icoor) is the coordinate icoor of the quadrature point ig
    const Real& quadPointCoor ( const UInt& ig, const UInt& icoor ) const
    {
        ASSERT_BD ( ig < M_nbQuadPt );
        return M_pt[ ig ].coor ( icoor );
    }

    //! quadPointCoor(ig) is the full coordinates of the quadrature point ig
    const GeoVector& quadPointCoor ( const UInt& ig ) const
    {
        ASSERT_BD ( ig < M_nbQuadPt );
        return M_pt[ ig ].coor( );
    }

    //! Getter for the number of quadrature points
    const UInt& nbQuadPt() const
    {
        return M_nbQuadPt;
    };

    //! Getter for the name of the quadrature
    const std::string& name() const
    {
        return M_name;
    };

    //! Getter for the degree of exactness
    const UInt& degreeOfExactness() const
    {
        return M_degOfExact;
    };

    //@}


private:

    //! @name Private Methods
    //@{

    //! Check the exactness for quadrature rules on segments
    UInt checkExactnessSegment() const;

    //! Check the exactness for quadrature rules on triangles
    UInt checkExactnessTriangle() const;

    //! Check the exactness for quadrature rules on tetrahedra
    UInt checkExactnessTetra() const;

    //@}

    //! Tolerance for the test of exactness
    static Real S_exactnessTol;

    // Storage for the quadrature nodes
    std::vector<QuadraturePoint> M_pt;

    //! geometrical shape of the domain on which the quadrature rule can be used
    ReferenceShapes M_shape;

    //! name of the quadrature rule
    std::string M_name;

    //! number of quadrature points
    UInt M_nbQuadPt;

    //! degree of exactness
    UInt M_degOfExact;

    //! Dimension in which the quadrature is stored
    UInt M_dimension;
};


// Definition of the template method
template< typename QRType>
void
QuadratureRule::import (const QRType& qr)
{
    for (UInt i (0); i < qr.nbQuadPt(); ++i)
    {
        M_pt.push_back (qr.quadPoint (i) );
    }
    M_shape = qr.shape();
    M_name = "";
    M_nbQuadPt = qr.nbQuadPt();
    M_degOfExact = 0;
    M_dimension = qr.dimension();
}



//======================================================================
extern const QuadratureRule quadRuleDummy;

extern const QuadratureRule quadRuleNode1pt;

extern const QuadratureRule quadRuleSeg1pt;
extern const QuadratureRule quadRuleSeg2pt;
extern const QuadratureRule quadRuleSeg3pt;
extern const QuadratureRule quadRuleSeg4pt;

extern const QuadratureRule quadRuleTria1pt;
extern const QuadratureRule quadRuleTria3pt;
extern const QuadratureRule quadRuleTria4pt;
extern const QuadratureRule quadRuleTria6pt;
extern const QuadratureRule quadRuleTria7pt;

extern const QuadratureRule quadRuleQuad1pt;
extern const QuadratureRule quadRuleQuad4pt;
extern const QuadratureRule quadRuleQuad9pt;
extern const QuadratureRule quadRuleQuad16pt;

extern const QuadratureRule quadRuleTetra1pt;
extern const QuadratureRule quadRuleTetra4pt;
extern const QuadratureRule quadRuleTetra5pt;
extern const QuadratureRule quadRuleTetra15pt;
extern const QuadratureRule quadRuleTetra64pt;

extern const QuadratureRule quadRuleHexa1pt;
extern const QuadratureRule quadRuleHexa8pt;


}
#endif
