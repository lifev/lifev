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
    @file Norm : distance between a given point
    @brief This file contains the computation of the distance between two points, used for the pseudo-ECG calculation

    @author Marie Dupraz <marie.dupraz@epfl.ch>
    @maintainer Marie Dupraz <marie.dupraz@epfl.ch>

    @date 2013-04-20

    Distance between two given points
 */


#include <lifev/core/LifeV.hpp>

namespace LifeV
{

class Norm
{
public:
    static Real f ( const Real& t, const Real& x, const Real& y,
                    const Real& z, const ID& i );
    //    {
    //        return  sqrt ( (x - M_xPosition) * (x - M_xPosition) + (y - M_yPosition) * (y - M_yPosition) + (z - M_zPosition) * (z - M_zPosition) ) ;
    //    }

    static inline void setPosition ( const Real& xPosition, const Real& yPosition, const Real& zPosition )
    {
        M_xPosition = xPosition;
        M_yPosition = yPosition;
        M_zPosition = zPosition;
    }

private:
    static Real  M_xPosition;
    static Real  M_yPosition;
    static Real  M_zPosition;

}; // class Laplacian

Real Norm::f ( const Real& /* t */, const Real& x, const Real& y, const Real& z, const ID& /* i */ )
{
    return  sqrt ( (x - M_xPosition) * (x - M_xPosition) + (y - M_yPosition) * (y - M_yPosition) + (z - M_zPosition) * (z - M_zPosition) ) ;
}

Real  Norm::M_xPosition = 1.;
Real  Norm::M_yPosition = 1.;
Real  Norm::M_zPosition = 0.5;

} // namespace LifeV

