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

 */

#ifndef QUADRATUREBOUNDARY_H
#define QUADRATUREBOUNDARY_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/fem/QuadratureRule.hpp>

#include <vector>

namespace LifeV
{

//! QuadratureBoundary - Short description of the class

class QuadratureBoundary
{
public:

    QuadratureBoundary() {};

    QuadratureBoundary ( const QuadratureBoundary& qrbd )
        : M_qrs (qrbd.M_qrs)
    {};


    void setQuadrature (const UInt i, const QuadratureRule& qr)
    {
        while (i >= M_qrs.size() )
        {
            M_qrs.push_back (QuadratureRule() );
        }

        M_qrs[i] = qr;
    }

    QuadratureRule qr (const UInt i) const
    {
        ASSERT (i < M_qrs.size(), "Invalid quadrature number");
        return M_qrs[i];
    }


    void showMe() const
    {
        for (UInt i (0); i < M_qrs.size(); ++i)
        {
            std::cout << " ## QR " << i << " ## " << std::endl;
            M_qrs[i].showMe();
        }
    }
private:

    std::vector< QuadratureRule > M_qrs;

};

inline
QuadratureBoundary
buildTetraBDQR (const QuadratureRule& my_qr)
{
    QuadratureBoundary qrbd;

    // Face 0
    qrbd.setQuadrature (0, my_qr);

    // Face 1
    QuadratureRule qf1 ("none", TRIANGLE, 3, 0, 0);
    for (UInt iq (0); iq < my_qr.nbQuadPt(); ++iq)
    {
        Real x (my_qr.quadPointCoor (iq, 0) );
        Real y (my_qr.quadPointCoor (iq, 1) );
        Real w (my_qr.weight (iq) );

        qf1.addPoint (QuadraturePoint (x, 0, y, w) );
    }
    qrbd.setQuadrature (1, qf1);

    // Face 2
    QuadratureRule qf2 ("none", TRIANGLE, 3, 0, 0);
    for (UInt iq (0); iq < my_qr.nbQuadPt(); ++iq)
    {
        Real x (my_qr.quadPointCoor (iq, 0) );
        Real y (my_qr.quadPointCoor (iq, 1) );
        Real w (my_qr.weight (iq) );

        qf2.addPoint (QuadraturePoint (x, y, 1 - x - y, w) );
    }
    qrbd.setQuadrature (2, qf2);

    // Face 3
    QuadratureRule qf3 ("none", TRIANGLE, 3, 0, 0);
    for (UInt iq (0); iq < my_qr.nbQuadPt(); ++iq)
    {
        Real x (my_qr.quadPointCoor (iq, 0) );
        Real y (my_qr.quadPointCoor (iq, 1) );
        Real w (my_qr.weight (iq) );

        qf3.addPoint (QuadraturePoint (0, y, x, w) );
    }
    qrbd.setQuadrature (3, qf3);

    return qrbd;

}

} // Namespace LifeV

#endif /* QUADRATUREBOUNDARY_H */
