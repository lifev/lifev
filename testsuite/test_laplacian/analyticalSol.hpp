/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

#ifndef ANALYTICALSOL_H_
#define ANALYTICALSOL_H_

#include <life/lifecore/life.hpp>

namespace LifeV
{

/*!
  \typedef The differential problem at hand
*/
enum ADRProblemSolution
{
    POISSON_POLYNOMIAL,       //!< A polynomial solution to the Poisson problem
    POISSON_TRIGONOMETRIC,    //!< A trigonometric solution to the Poisson problem
    ADR_STEADY_POLYNOMIAL,    //!< A polynomial solution to a general ADR problem
    ADR_UNSTEADY_POLYNOMIAL   //!< A polynomial solution to a general ADR problem
};


// Coefficient for Robin boundary condition
Real coefRobin(const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 2.0;
}

#ifdef TWODIM

class AnalyticalSol
{
public:
    // Give to the class a "functor" interface
    Real operator()(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return u_ex( t, x, y, z, ic );
    }
    // This method is required by the error calculation routines
    Real grad(UInt icoor, const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return grad_ex( icoor, t, x, y, z, ic );
    }
    // The exact solution to the problem
    static Real u_ex(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const UInt& /*ic=0*/)
    {
        //        return x*x + y*y; // Polynomial solution
        return sin(2*Pi*x)*sin(2*Pi*y); // Trigonometric solution
    }
    // The gradient of the exact solution
    static Real grad_ex(const UInt& icoor, const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const UInt& /*ic=0*/) {
        switch(icoor)
        {
        case 1: //der_x
            //            return 2*x; // Polynomial solution
            return 2*Pi*cos(2*Pi*x)*sin(2*Pi*y); // Trigonometric solution
        case 2: //der_y
            //            return 2*x; // Polynomial solution
            return 2*Pi*sin(2*Pi*x)*cos(2*Pi*y); // Trigonometric solution
        default:
            return 0;
        }
    }
    // The source term
    static Real source(const Real& /*t*/, const Real& x, const Real& y, const Real& /*z*/, const UInt& /*ic=0*/)
    {
        //        return -4; // Polynomial solution
        return 8*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y); // Trigonometric solution
    }
    // The normal derivative to the domain boundary
    // BEWARE: this is case dependent! (works for cubic domains)
    static Real fNeumann(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic)
    {
        Real n[2] = {0, 0}; Real out=0;
        if        ( x == xMin ) {
            n[0] = -1.;
        } else if ( x == xMax ) {
            n[0] =  1.;
        } else if ( y == yMin ) {
            n[1] = -1.;
        } else if ( y == yMax ) {
            n[1] =  1.;
        } else {
            std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                    << std::endl;
        }

        for (UInt k =0; k< 2; k++)  //mu gradu n
            out += grad_ex(k+1, t, x, y, z, ic)*n[k];

        return out;
    }
    // Robin boundary condition compatible with the exact solution
    static Real fRobin(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic)
    {
        return fNeumann(t,x,y,z,ic)+coefRobin(t,x,y,z,ic)*u_ex(t,x,y,z,ic);
    }
    // Initialization of the private members
    static void setup( const GetPot& dataFile, const std::string& section = "poisson" )
    {
        // We want a slash dividing the data file section from the variable name but only
        // when not provided by the user or when not looking in the root of the data file
        std::string corrected_section( section );
        if( ( ! section.empty() ) && ( section[section.length()-1] != '/' ) )
            corrected_section = section + '/';

        xMin = dataFile( (corrected_section+"problem/xMin").data(), -1. ) ;
        xMax = dataFile( (corrected_section+"problem/xMax").data(), 1. ) ;
        yMin = dataFile( (corrected_section+"problem/yMin").data(), -1. ) ;
        yMax = dataFile( (corrected_section+"problem/yMax").data(), 1. ) ;
    }

private:
    // The dimensions of the domain
    static Real xMin;
    static Real xMax;
    static Real yMin;
    static Real yMax;

};

// Declaration and default initialization of the static variables
Real AnalyticalSol::xMin(0.);
Real AnalyticalSol::xMax(1.);
Real AnalyticalSol::yMin(0.);
Real AnalyticalSol::yMax(1.);



Real beta( const Real& /* t */,
           const Real& ,
           const Real& ,
           const Real& ,
           const ID& icomp )
{
    switch(icomp)
    {
    case 1:
        return 0;
    case 2:
        return 0;
    default:
        return 0;
    }
}


#elif defined THREEDIM

class AnalyticalSol
{
public:
    // Give to the class a "functor" interface
    Real operator()(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return u_ex( t, x, y, z, ic );
    }
    // This method is required by the error calculation routines
    Real grad(UInt icoor, const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic) const
    {
        return grad_ex( icoor, t, x, y, z, ic );
    }
    // The exact solution to the problem
    static Real u_ex(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& /*ic*/)
    {
        switch( problemSolution )
        {
        case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
        case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            return x*x + y*y + z*z;
        case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
            return sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z);
        case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
            return t*t + x*x + y*y + z*z;
        default:
            return 0;
        }
    }
    // The gradient of the exact solution
    static Real grad_ex(const UInt& icoor,
                        const Real& /*t*/, const Real& x, const Real& y, const Real& z, const UInt& /*ic*/)
    {
        switch(icoor)
        {
        case 1: //der_x
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 2*x; // Polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 2*Pi*cos(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z); // Trigonometric solution
            }
        case 2: //der_y
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 2*y; // Polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 2*Pi*sin(2*Pi*x)*cos(2*Pi*y)*sin(2*Pi*z); // Trigonometric solution
            }
        case 3: //der_z
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 2*z; // Polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 2*Pi*sin(2*Pi*x)*sin(2*Pi*y)*cos(2*Pi*z); // Trigonometric solution
            }
        default:
            return 0;
        }
    }
    // The source term
    static Real fSource(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic)
    {
        Real value(0.);
        switch( problemSolution )
        {
        case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            return -diffusionCoeff*6.; // Polynomial solution
        case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
            return 12*Pi*Pi*sin(2*Pi*x)*sin(2*Pi*y)*sin(2*Pi*z); // Trigonometric solution
        case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
            value += 2*t;
        case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            value += -diffusionCoeff*6 + reactionCoeff*u_ex(t, x, y, z, ic);
            for(UInt icomp=1; icomp<=nDimensions; ++icomp)
                value += fAdvection(t, x, y, z, icomp)*grad_ex(icomp, t, x, y, z, ic);
            return value;
        default:
            return 0;
        }
    }
    // The advection term
    static Real fAdvection(const Real& /*t*/, const Real& x, const Real& y, const Real& z, const UInt& icomp)
    {
        switch(icomp)
        {
        case 1:
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 0.;
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 0.; //-0.1*x;
            }
        case 2:
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 0.;
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 0.; //-0.1*y;
            }
        case 3:
            switch( problemSolution )
            {
            case POISSON_POLYNOMIAL: // Poisson problem, polynomial solution
            case POISSON_TRIGONOMETRIC: // Poisson problem, trigonometric solution
                return 0.;
            case ADR_STEADY_POLYNOMIAL: // ADR problem, polynomial solution
            case ADR_UNSTEADY_POLYNOMIAL: // ADR problem, polynomial solution
                return 0.; //-0.1*z;
            }
        default:
            return 0.;
        }
    }
    // The normal derivative to the domain boundary
    // BEWARE: this is case dependent! (works for cubic domains)
    static Real fNeumann(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic)
    {
        Real n[3] = {0, 0, 0}; Real out=0;
        if        ( x == xMin ) {
            n[0] = -1.;
        } else if ( x == xMax ) {
            n[0] =  1.;
        } else if ( y == yMin ) {
            n[1] = -1.;
        } else if ( y == yMax ) {
            n[1] =  1.;
        } else if ( z == zMin ) {
            n[2] = -1.;
        } else if ( z == zMax ) {
            n[2] =  1.;
        } else {
            std::cout << "strange point: x=" << x << " y=" << y << " z=" << z
                    << std::endl;
        }

        for (UInt icomp =0; icomp< nDimensions; ++icomp)  //mu grad_u n
            out += diffusionCoeff*grad_ex(icomp+1, t, x, y, z, ic)*n[icomp];
        return out;
    }
    // Robin boundary condition compatible with the exact solution
    static Real fRobin(const Real& t, const Real& x, const Real& y, const Real& z, const UInt& ic)
    {
        return fNeumann(t,x,y,z,ic)+coefRobin(t,x,y,z,ic)*u_ex(t,x,y,z,ic);
    }
    // A generic diffusion coefficient
    //static Real diffusionCoeff(const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const UInt& /*ic*/)
    //{
    //    return t*t;
    //}
    // Initialization of the private members
    static void setup( const GetPot& dataFile, const std::string& section = "poisson" )
    {
        // We want a slash dividing the data file section from the variable name but only
        // when not provided by the user or when not looking in the root of the data file
        std::string corrected_section( section );
        if( ( ! section.empty() ) && ( section[section.length()-1] != '/' ) )
            corrected_section = section + '/';

        xMin = dataFile( (corrected_section+"problem/xMin").data(), -1. ) ;
        xMax = dataFile( (corrected_section+"problem/xMax").data(), 1. ) ;
        yMin = dataFile( (corrected_section+"problem/yMin").data(), -1. ) ;
        yMax = dataFile( (corrected_section+"problem/yMax").data(), 1. ) ;
        zMin = dataFile( (corrected_section+"problem/zMin").data(), -1. ) ;
        zMax = dataFile( (corrected_section+"problem/zMax").data(), 1. ) ;

        diffusionCoeff = dataFile( (corrected_section+"physics/diffusionCoefficient").data(), 1. ) ;
        reactionCoeff = dataFile( (corrected_section+"physics/reactionCoefficient").data(), 1. ) ;
    }
    // Output on screen
    static void showMe( std::ostream& c = std::cout )
    {
        c << "\n*** Analytical solution: values for user-defined data\n";

        c << "\n[/problem]" << std::endl;
        c << "xMin\t\t= " << xMin << std::endl;
        c << "xMax\t\t= " << xMax << std::endl;
        c << "yMin\t\t= " << yMin << std::endl;
        c << "yMax\t\t= " << yMax << std::endl;
        c << "zMin\t\t= " << zMin << std::endl;
        c << "zMax\t\t= " << zMax << std::endl;
        c << "solutionType\t\t= " << problemSolution << std::endl;

        c << "\n[/physics]" << std::endl;
        c << "diffusionCoeff\t\t= " << diffusionCoeff << std::endl;
        c << "reactionCoeff\t\t= " << reactionCoeff << std::endl;

        c << std::endl;
    }

    // The dimensions of the domain
    static Real xMin;
    static Real xMax;
    static Real yMin;
    static Real yMax;
    static Real zMin;
    static Real zMax;
    static Real diffusionCoeff;
    static Real reactionCoeff;

    static ADRProblemSolution problemSolution;

};

// Declaration and default initialization of the static variables
Real AnalyticalSol::xMin(0.);
Real AnalyticalSol::xMax(1.);
Real AnalyticalSol::yMin(0.);
Real AnalyticalSol::yMax(1.);
Real AnalyticalSol::zMin(0.);
Real AnalyticalSol::zMax(1.);
Real AnalyticalSol::diffusionCoeff(1.);
Real AnalyticalSol::reactionCoeff(1.);

ADRProblemSolution AnalyticalSol::problemSolution(POISSON_POLYNOMIAL);

#endif


} //namespace LifeV

#endif
