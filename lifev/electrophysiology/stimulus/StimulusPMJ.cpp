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
 @brief Class for applying cardiac stimulus at Purkinje-muscle junctions

 @date 11-2013
 @author Toni Lassila <toni.lassila@epfl.ch>

 @last update 11-2013

 */

#include <lifev/electrophysiology/stimulus/StimulusPMJ.hpp>

namespace LifeV
{

// ===================================================
//! Constructors
// ===================================================
StimulusPMJ::StimulusPMJ() :
    M_radius ( 0 ),
    M_activationData ( * (new activationData_type() ) ),
    M_problemFolder ( "./" )
{

}

// ===================================================
//! Setters
// ===================================================
void StimulusPMJ::setPMJFromFile ( std::string fileName )
{
    std::ifstream fin;

    fin.open ( fileName.c_str() );
    while ( !fin.fail() )
    {
        StimulusPMJ_Activation junction;
        fin >> junction.x;
        fin >> junction.y;
        fin >> junction.z;
        fin >> junction.time;
        fin >> junction.duration;
        M_activationData.push_back ( junction );
    }
    fin.close();
}

void StimulusPMJ::setPMJAddJunction ( Real x, Real y, Real z, Real time, Real duration )
{
    StimulusPMJ_Activation junction;
    junction.x = x;
    junction.y = y;
    junction.z = z;
    junction.time = time;
    junction.duration = duration;
    M_activationData.push_back ( junction );
}

// ===================================================
//! Methods
// ===================================================
Real StimulusPMJ::appliedCurrent ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& /*i*/ )
{
    Real current = 0;
    const Real volumeOfBall = (4. / 3.) * M_PI * M_radius * M_radius * M_radius;

    for (activationData_type::iterator it = M_activationData.begin(); it != M_activationData.end(); ++it)
    {
        Real distance = std::sqrt ( (x - it->x) * (x - it->x) + (y - it->y) * (y - it->y) + (z - it->z) * (z - it->z) );


        if (distance <= M_radius && t >= (it->time) && t <= (it->time + it->duration) )
        {
            current += M_totalCurrent / volumeOfBall;
        }
    }

    return current;
}

void StimulusPMJ::showMe()
{
    std::cout << "\n\n\t\tPMJ activation Informations\n\n";

    std::cout << "\n\t\tList of parameters:\n\n";

    std::cout << "Radius current application: " << M_radius << std::endl;
    std::cout << "Total current: " << M_totalCurrent << std::endl;
    std::cout << "Problem folder: " << M_problemFolder << std::endl;
    std::cout << "\n\t\t End of PMJ activation Informations\n\n\n";
}

}
