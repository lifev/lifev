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
 *  @brief File containing the boundary conditions for the Monolithic Test
 *
 *  @date 2009-04-09
 *  @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  @contributor Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @maintainer Paolo Crosetto <crosetto@iacspc70.epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#ifndef FRES_HPP
#define FRES_HPP

// LifeV includes
#include <lifev/core/LifeV.hpp>

namespace LifeV
{
    
    Real exactVelocity (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i){
        
        Real a = 1.57;
        
        Real d = 0.78;
        
        Real nu = 0.01;
        
        Real e = std::exp (-d*d*t*nu);
        
        switch (i)
        {
            case 0:
                return -a * e * ( std::exp (a * x) * std::sin (a * y + d * z) + std::exp (a * z) * std::cos (a * x + d * y) );
            case 1:
                return -a * e * ( std::exp (a * y) * std::sin (a * z + d * x) + std::exp (a * x) * std::cos (a * y + d * z) );
            case 2:
                return -a * e * ( std::exp (a * z) * std::sin (a * x + d * y) + std::exp (a * y) * std::cos (a * z + d * x) );
            default:
                exit (1);
        }
        
    }
    
    Real gradientVelocity ( const UInt& icoor, const Real& t, const Real& x, const Real& y, const Real& z, const ID& i){
        
        Real a = 1.57;
        
        Real d = 0.78;
        
        Real nu = 0.01;
        
        Real e = std::exp (-d* d* t * nu);
        
        switch (icoor)
        {
            case 0:
                switch (i)
            {
                case 0:
                    return -a* e * ( a* std::exp (a* x) * std::sin (a* y + d* z) - a* std::exp (a* z) * std::sin (a* x + d* y) );
                case 1:
                    return -a* e * ( d* std::exp (a* y) * std::cos (a* z + d* x) + a* std::exp (a* x) * std::cos (a* y + d* z) );
                case 2:
                    return -a* e * ( a* std::exp (a* z) * std::cos (a* x + d* y) - d* std::exp (a* y) * std::sin (a* z + d* x) );
                default:
                    exit (1);
            }
            case 1:   // u_y
                switch (i)
            {
                case 0:
                    return -a* e * ( a* std::exp (a* x) * std::cos (a* y + d* z) - d* std::exp (a* z) * std::sin (a* x + d* y) );
                case 1:
                    return -a* e * ( a* std::exp (a* y) * std::sin (a* z + d* x) - a* std::exp (a* x) * std::sin (a* y + d* z) );
                case 2:
                    return -a* e * ( d* std::exp (a* z) * std::cos (a* x + d* y) + a* std::exp (a* y) * std::cos (a* z + d* x) );
                default:
                    exit (1);
            }
            case 2:
                switch (i)
            {
                case 0:
                    return -a* e * ( d* std::exp (a* x) * std::cos (a* y + d* z) + a* std::exp (a* z) * std::cos (a* x + d* y) );
                case 1:
                    return -a* e * ( a* std::exp (a* y) * std::cos (a* z + d* x) - d* std::exp (a* x) * std::sin (a* y + d* z) );
                case 2:
                    return -a* e * ( a* std::exp (a* z) * std::sin (a* x + d* y) - a* std::exp (a* y) * std::sin (a* z + d* x) );
                default:
                    exit (1);
            }
            default:
                exit (1);
        }
        return 1.;
        
        
    }
    
    Real exactPressure (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i){
        
        Real a = 1.58;
        
        Real d = 0.78;
        
        Real nu = 0.01;
        
        return - a*a/2.0*std::exp(-2.0*d*d*nu*t)*(std::exp (2.0*a*x)+std::exp(2.0*a*y)+std::exp(2.0*a*z)    +
                                                2.0*std::sin(a*x+d*y)*std::cos(a*z+d*x)*std::exp(a*(y+z)) +
                                                2.0*std::sin(a*y+d*z)*std::cos(a*x+d*y)*std::exp(a*(z+x)) +
                                                2.0*std::sin(a*z+d*x)*std::cos(a*y+d*z)*std::exp(a*(x+y)) );
        
    }
    
    Real normalStress (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i){
        
        Real dinamic_viscosity = 0.01;
        
        // Caution: accordingly to what is done in the driver navier_stokes.hpp, here I assume that the normal of the Neumann face
        //          is the one whose outward normal unit vector is n = [ 0; 0; 1 ].
        
        Real n[3] = {0.0, 0.0, -1.0};
        
        Real out = 0.0;
        
        for (UInt k = 0; k < 3; k++){
            
            out += dinamic_viscosity * gradientVelocity (k, t, x, y, z, i) * n[k];
        
        }
        
        for (UInt k = 0; k < 3; k++){
            
            out += dinamic_viscosity * gradientVelocity (i, t, x, y, z, k) * n[k];
        
        }
        
        out -= exactPressure (t, x, y, z, i) * n[i];
        
        return out;
        
    }
    
} // end Namespace LifeV

#endif // end of define FRES_HPP