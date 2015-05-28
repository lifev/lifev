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
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 21 Feb 2012
 */

#ifndef QUADRATURERULEBOUNDARY_H
#define QUADRATURERULEBOUNDARY_H 1

#include <lifev/core/LifeV.hpp>
#include <lifev/core/fem/QuadratureRule.hpp>

#include <vector>

namespace LifeV
{

//! QuadratureRuleBoundary - Short description of the class
/*!
    @author Samuel Quinodoz
    @see Reference to papers (if available)
 */
class QuadratureRuleBoundary
{
public:

    QuadratureRuleBoundary() {};

    QuadratureRuleBoundary (const QuadratureRuleBoundary& QRBD)
    {
        for (UInt i (0); i < QRBD.M_quadrature.size(); ++i)
        {
            M_quadratures.push_back (QRBD.M_quadratures[i]);
        }
    }

    ~QuadratureRuleBoundary() {}

    QuadratureRule qr (const UInt i) const
    {
        return M_quadratures[i];
    }

    void setQR (const UInt i, QuadratureRule myQR)
    {
        while (M_quadratures.size() <= i)
        {
            M_quadratures.push_back (QuadratureRule() );
        }

        M_quadratures[i] = myQR;
    }

private:

    // Store the possible quadrature rules
    std::vector<QuadratureRule> M_quadratures

};


inline QuadratureRuleBoundary
createTetraBDQR ( const QuadratureRule& myQR)
{
    QuadratureRuleBoundary QRBD;

    QRBD.setQR (0, myQR);

    QuadratureRule F1;
    for (UInt i (0); i < myQR.nbQuadPt(); ++i)
    {
        Real x0 (myQR.quadPointCoor (i, 0) );
        Real y0 (myQR.quadPointCoor (i, 1) );

        Real x (x0);
        Real y (0);
        Real z (y0);

        F1.addPoint (QuadraturePoint (x, y, z, myQR.weight (i) ) );
    }
    QRBD.setQR (1, F1);

    QuadratureRule F2;
    for (UInt i (0); i < myQR.nbQuadPt(); ++i)
    {
        Real x0 (myQR.quadPointCoor (i, 0) );
        Real y0 (myQR.quadPointCoor (i, 1) );

        Real x (x0);
        Real y (y0);
        Real z (1 - x0 - y0);

        F2.addPoint (QuadraturePoint (x, y, z, myQR.weight (i) *std::sqrt (3) / (2.0 * sqrt (2) ) ) );
    }
    QRBD.setQR (2, F2);

    QuadratureRule F3;
    for (UInt i (0); i < myQR.nbQuadPt(); ++i)
    {
        Real x0 (myQR.quadPointCoor (i, 0) );
        Real y0 (myQR.quadPointCoor (i, 1) );

        Real x (0);
        Real y (y0);
        Real z (x0);

        F3.addPoint (QuadraturePoint (x, y, z, myQR.weight (i) ) );
    }
    QRBD.setQR (3, F3);

    return QRBD;
}


} // Namespace LifeV

#endif /* QUADRATURERULEBOUNDARY_H */
