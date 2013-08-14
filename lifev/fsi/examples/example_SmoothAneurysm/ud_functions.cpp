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

#include "ud_functions.hpp"


namespace LifeV
{

Real f (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.;
}

Real u1 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real fZero (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real outerWallPressure (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& z, const ID& /*i*/)
{
    Real highestPressure( - ( 13330 - 113305 ) );
    Real pressure(0);
    Real totalTime = 0.8;
    Real halfTime = totalTime / 2.0;

    Real a = ( highestPressure / 2 ) * ( 1/ ( halfTime * halfTime ) );

    if ( t <= totalTime )
        pressure = ( highestPressure / totalTime ) * t;
    else
        pressure = highestPressure;

    // if ( t <= 0.8 )
    // {
    //     return ( value / ( 0.8 * 0.8 * 0.8 * 0.8 ) ) * ( t * t *t *t );
    // }
    // else
    // {
    //     return value;
    // }

    return -pressure;

}

Real epsilon (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 113487.36193;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real pressureInitial (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 113305;
}

// Initial velocity
Real u0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real p0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}


Real E (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    if( t <= 0.4 )
        return -10000000 * ( t - 0.4 );
    else
        return 0.0;

    return 0.0;

}


Real hydro (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1.33 * 1e5;
}

Real u2 (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 2:
            if ( t <= 0.003 )
            {
                return 1.5;
            }
            //      return 0.01;
            return 0.0;
            break;
        case 1:
            return 0.0;
            //      return 1.3332e4;
            //    else
            //      return 0.0;
            break;
    }
    return 0;
}


// Initial displacement and velocity
Real d0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{
    switch (i)
    {
        case 0:
            return 0.;
            break;
        case 1:
            return 0.;
            break;
        case 2:
            return 0.;
            break;
        default:
            ERROR_MSG ("This entry is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real w0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 10.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real fluxFunctionAneurysm (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{

    Real fluxFinal;
    Real rampAmpl (1.012);
    Real activeRamp ( rampAmpl / 1.0 );
    Real dt (0.001);

    if ( t <= activeRamp )
      {
        fluxFinal = ( 0.1747 / activeRamp ) * t; // 0.033 cm
      }
    else if ( t > activeRamp  && t <= rampAmpl )
      {
    	fluxFinal = ( 0.1747 );
      }
    else
      {


        // We change the flux for our geometry
        const Real pi   = 3.141592653589793;
        const Real area = 0.0033853;  // radius = 0.033 cm
        // const Real area = 0.013684; // radius = 0.066 cm

        const Real areaFactor = area / ( 0.195 * 0.195  * pi);
        //const Real Average = (48.21 * pow (area, 1.84) ) * 60; //Mean Cebral's Flux per minut

        // Unit conversion from ml/min to cm^3/s
        const Real unitFactor = 1. / 60.;

        // T is the period of the cardiac cycle
        const Real T          = 0.8;

        // a0 is the average VFR (the value is taken from Karniadakis p970)
        const Real a0         = 255;
        //const Real volumetric = Average / a0; //VolumetricFactor per minut

        // Fourrier
        const Int M (7);
        const Real a[M] = { -0.152001, -0.111619, 0.043304, 0.028871, 0.002098, -0.027237, -0.000557};
        const Real b[M] = { 0.129013, -0.031435, -0.086106, 0.028263, 0.010177, 0.012160, -0.026303};

        Real flux (0);
        const Real xi (2 * pi * (t - 0.8 + dt) / T);

        flux = a0;
        Int k (1);
        for (; k <= M ; ++k)
        {
            flux += a0 * (a[k - 1] * cos (k * xi) + b[k - 1] * sin (k * xi) );
        }

        fluxFinal =  (flux * areaFactor * unitFactor);
      }

    std::cout << "Flux that is imposed" << fluxFinal << std::endl;

    return fluxFinal;

}

Real aneurismFluxInVectorial (const Real&  t, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real n1 (0.0);
    Real n2 (0.0);
    Real n3 (1.0);

    Real x0 (0.0);
    Real y0 (0.0);


    Real flux (fluxFunctionAneurysm (t, x, y, z, i) );
    Real area (0.0033853); // 0.033 cm
    //Real area (0.013684);   // 0.066 cm

    //Parabolic profile
    Real radius(0.033);
    Real radiusSquared = radius * radius;
    Real peak(0);
    peak = ( 2 * flux ) / ( area );

    switch (i)
    {
    case 0:
        // Flat profile: flux / area;
        // return n1 * flux / area;
        return n1 * ( peak * ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ) ;
    case 1:
        // Flat profile: flux / area;
        //return n2 * flux / area;
        return n2 * ( peak * ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ) ;
    case 2:
        // Flat profile: flux / area;
        //return n3 * flux / area;
        return n3 * ( peak * ( (radiusSquared - ( (x-x0)*(x-x0) + (y-y0)*(y-y0)) )/radiusSquared) ) ;
    default:
        return 0.0;
    }
}


Real squareSinusoidalFluxFunction (const Real& t, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return - (t < (0.005 / 2) ) * std::sin (2 * M_PI * t / 0.005) * std::sin (2 * M_PI * t / 0.005);
}

//----------------------------------------------Fibers Directions--------------
Real Family1 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real theta =  0.9865; // value for anisotropic characterization taken from Robertson // ( PI / 6.0 );
    Real thetaChangeOfVariable = std::atan(  y / x );

    if( x < 0 )
    {
        // This is due to the periodicity of std::atan ( ref. official documentation )
        Real pi(3.141592653589793);
        thetaChangeOfVariable += pi;
    }

    switch (i)
    {
        case 0:
	    // Tube
            return - std::sin( thetaChangeOfVariable ) * std::cos( theta );
	    // Cube
            // return std::sin( theta );
            break;
        case 1:
	    // Tube
            return   std::cos( thetaChangeOfVariable ) * std::cos( theta );
	    // Cube
            // return std::cos( theta );
            break;
        case 2:
	    // Tube
            return std::sin( theta );
	    // Cube
            // return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family2 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
    Real theta = - 0.9865; //( - PI / 6.0 );
    Real thetaChangeOfVariable = std::atan( y / x );

    if( x < 0 )
    {
        // This is due to the periodicity of std::atan ( ref. official documentation )
        Real pi(3.141592653589793);
        thetaChangeOfVariable += pi;
    }

    switch (i)
    {
        case 0:
	    // Tube
            return - std::sin( thetaChangeOfVariable ) * std::cos( theta );
	    // Cube
            // return std::sin( theta );
            break;
        case 1:
	    // Tube
            return   std::cos( thetaChangeOfVariable ) * std::cos( theta );
	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
            return   std::sin( theta );
	    // Cube
	    // return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family3 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return -1.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real Family4 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}


Real Family5 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family6 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
        case 0:
            return 0.0;
            break;
        case 1:
            return 0.0;
            break;
        case 2:
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

// Method for the definition of the fibers
fibersDirectionList::fibersDirectionList() :
    M_mapNameDefinition( )
{}

fibersDirectionList::~fibersDirectionList()
{}

void fibersDirectionList::setupFiberDefinitions( const UInt nbFamilies )
{
    // At the moment the creation of the table of fiber functions is done
    // manually. There should be a way to make it automatically. Btw, only
    // the first nbFamilies that are set in the data file are taken into account

    ASSERT( nbFamilies < 6, "At the moment, a maximum number = 6 of families can be used! If you want more \n modifiy the file ud_functions.hpp in the application folder." );

    // Creation of the database of functions
    fiberFunctionPtr_Type pointerToFunction( new fiberFunction_Type( Family1 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family1", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family2 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family2", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family3 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family3", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family4 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family4", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family5 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family5", pointerToFunction ) );

    pointerToFunction.reset( new fiberFunction_Type( Family6 ) );
    M_mapNameDefinition.insert( std::pair<std::string, fiberFunctionPtr_Type>
                                  ( "Family6", pointerToFunction ) );


}

fibersDirectionList::fiberFunctionPtr_Type fibersDirectionList::fiberDefinition( const std::string nameFamily )
{

    mapNameDefinitionFiberFunction_Type::const_iterator IT;

    IT = M_mapNameDefinition.find ( nameFamily );

    if ( IT != M_mapNameDefinition.end() )
    {
        return IT->second;
    }
    else
    {
        std::cout << " Wrong identification of the fiber function! " << std::endl;
        fiberFunctionPtr_Type pointerToFunction( new fiberFunction_Type() );

        return pointerToFunction;
    }
}

}


