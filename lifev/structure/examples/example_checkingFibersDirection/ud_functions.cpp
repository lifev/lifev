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
 *  @author Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  @maintainer Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 *  Contains the functions to be assigned as boundary conditions, in the file boundaryConditions.hpp . The functions
 *  can depend on time and space, while they can take in input an ID specifying one of the three principal axis
 *  if the functions to assign is vectorial and the boundary condition is of type \c Full \c.
 */

#include "ud_functions.hpp"
#include<lifev/core/array/VectorSmall.hpp>
#include<lifev/core/array/MatrixSmall.hpp>

#define PI 3.14159265359

namespace LifeV
{

Real f (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
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


Real fzero_scalar (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return 0.0;
}

Real InternalPressure (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*i*/)
{
    return -1e+5;
    //return -260000*sin(80*3.141592*t);
}

// Initial displacement and velocity
Real d0 (const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    switch (i)
    {
    case 0:
        return  0.088002 * ( x + 0.5 );
        break;
    case 1:
        return - ( 0.02068 * 2.0 ) * ( y );
        break;
    case 2:
        return - ( 0.02068 * 2.0 ) * ( z );
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
            return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real a0 (const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& i)
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


//----------------------------------------------Boundary Conditions--------------
Real bcZero (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  0.;
}

Real bcNonZero (const Real& t, const Real&  X, const Real& Y, const Real& Z, const ID& i)
{

    return 19995;

	// Real highestPressure(6.666e+6);
	// Real totalTime = 20.0;
	// Real halfTime = totalTime / 2.0;

	// Real a = ( highestPressure / 2 ) * ( 1/ ((totalTime/2)*(totalTime/2)) );

	// if ( t <= halfTime )
	//     pressure = a * t*t;

	// if ( t > halfTime )
	//     pressure = - a * (t - totalTime)*(t - totalTime) + highestPressure;

    // switch (i)
    // {
    //     case 0:
    //         return 0.0;
    //         break;
    //     case 1:
    //         return pressure;
    //         break;
    //     case 2:
    //         return 0.0;
    //         break;
    //     default:
    //         ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
    //         return 0.;
    //         break;
    // }
}
 Real smoothPressure(const Real& t, const Real&  x, const Real& y, const Real& /*Z*/, const ID& i)
    {
        Real radius = std::sqrt( x*x + y*y);
        Real pressure(0);

        Real highestPressure(200000);
        Real totalTime = 4.5;
        Real halfTime = totalTime / 2.0;

        // Real a = ( highestPressure / 2 ) * ( 1/ ( halfTime*halfTime ) );

        // if ( t <= halfTime )
        //     pressure = a * t*t;

        // if ( t > halfTime )
        //     pressure = - a * (t - totalTime)*(t - totalTime) + highestPressure;

        pressure = 5000 * t;

        switch (i)
        {
        case 0:
            return  pressure *  ( x / radius ) ;
            break;
        case 1:
            return  pressure *  ( y / radius ) ;
            break;
        case 2:
            return 0.0;
            break;

        }
        return 0;

    }

Real bcNonZeroSecondOrderExponential (const Real& /*t*/, const Real&  /*X*/, const Real& /*Y*/, const Real& /*Z*/, const ID& /*i*/)
{
    return  19180.;
}

Real thetaFunction ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  return thetaPosition;
}

Real thetaRotationFunction ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  return thetaRotation;
}


//----------------------------------------------Fibers Directions--------------
Real Family1 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    // Tube
    // Real theta =  0.8426; // value for anisotropic characterization taken from Robertson // ( PI / 6.0 );
    // Real thetaChangeOfVariable = std::atan(  y / x );

    // if( x < 0 )
    // {
    //     // This is due to the periodicity of std::atan ( ref. official documentation )
    //     thetaChangeOfVariable += PI;
    // }

  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  Real xCenter;
  Real zCenter;
  xCenter = rT * std::cos( thetaPosition ) + xT;
  zCenter = rT * std::sin( thetaPosition );

  VectorSmall<3> positionCenter;
  positionCenter[0] = xCenter;
  positionCenter[1] = 0;
  positionCenter[2] = zCenter;

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  MatrixSmall<3,3> changeOfVariable;
  changeOfVariable(0,0) = std::cos( thetaRotation );  changeOfVariable(0,1) = 0.0; changeOfVariable(0,2) = std::sin( thetaRotation );
  changeOfVariable(1,0) = 0.0;                        changeOfVariable(1,1) = 1.0; changeOfVariable(1,2) = 0.0;
  changeOfVariable(2,0) = -std::sin( thetaRotation ); changeOfVariable(2,1) = 0.0; changeOfVariable(2,2) = std::cos( thetaRotation );

  VectorSmall<3> localPosition;
  localPosition = changeOfVariable * ( position - positionCenter );

  Real thetaPositionOnSection(0);
  thetaPositionOnSection = std::atan( localPosition[2] / localPosition[1] );

  if ( localPosition[1] < 0 )
    thetaPositionOnSection += PI;

  VectorSmall<3> localFibers;
  localFibers[0] = std::sin( thetaFiber );
  localFibers[1] = std::cos( thetaFiber ) *  std::sin( thetaPositionOnSection );
  localFibers[2] = -  std::cos( thetaFiber ) *  std::cos( thetaPositionOnSection );

  VectorSmall<3> originFibers;
  originFibers = changeOfVariable.transpose() * localFibers + positionCenter;

  originFibers.normalize();

    switch (i)
    {
        case 0:
	    // Tube
  	    // return - std::sin( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
	    return originFibers[0];

	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
   	    // return   std::cos( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
            return   originFibers[1];

	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
 	    //return   std::sin( theta );	    

	    // half thorus
            return  originFibers[2];

	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }

}

Real Family1Spherical ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  VectorSmall<3> centerSphere;
  centerSphere[0] = 0.138;
  centerSphere[1] = 0.0;
  centerSphere[2] = 0.185;

  VectorSmall<3> centerOnTorusAxis;
  centerOnTorusAxis[0] = rT * std::cos( thetaPosition ) + xT;
  centerOnTorusAxis[1] = 0.0;
  centerOnTorusAxis[2] = rT * std::sin( thetaPosition );

  Real xCenter;
  Real yCenter;
  Real zCenter;

  if( ( ( position - centerSphere ).norm() ) <= 0.05 ) // inside the dome
    {
      VectorSmall<3> tangentVectorToAxisTorus;
      tangentVectorToAxisTorus[0] = std::sin( thetaPosition ) ;
      tangentVectorToAxisTorus[1] = 0.0;
      tangentVectorToAxisTorus[2] = -std::cos( thetaPosition );
      
      VectorSmall<3> centerSphereOnSection;
      centerSphereOnSection = centerSphere - ( ( centerSphere - centerOnTorusAxis ).dot( tangentVectorToAxisTorus ) ) * tangentVectorToAxisTorus;

      xCenter = centerSphereOnSection[0];
      yCenter = centerSphereOnSection[1];
      zCenter = centerSphereOnSection[2];
    }
  else
    { 
      xCenter = centerOnTorusAxis[0];
      yCenter = centerOnTorusAxis[1];
      zCenter = centerOnTorusAxis[2];
    }

  // Changing reference system
  MatrixSmall<3,3> changeOfVariable;
  changeOfVariable(0,0) = std::cos( thetaRotation );  changeOfVariable(0,1) = 0.0; changeOfVariable(0,2) = std::sin( thetaRotation );
  changeOfVariable(1,0) = 0.0;                        changeOfVariable(1,1) = 1.0; changeOfVariable(1,2) = 0.0;
  changeOfVariable(2,0) = -std::sin( thetaRotation ); changeOfVariable(2,1) = 0.0; changeOfVariable(2,2) = std::cos( thetaRotation );

  VectorSmall<3> positionCenter;
  positionCenter[0] = xCenter;
  positionCenter[1] = yCenter;
  positionCenter[2] = zCenter;

  VectorSmall<3> localPosition;
  localPosition = changeOfVariable * ( position - positionCenter );

  Real thetaPositionOnSection(0);
  thetaPositionOnSection = std::atan( localPosition[2] / localPosition[1] );

  if ( localPosition[1] < 0 )
    thetaPositionOnSection += PI;

  VectorSmall<3> localFibers;
  localFibers[0] = std::sin( thetaFiber );
  localFibers[1] = std::cos( thetaFiber ) *  std::sin( thetaPositionOnSection );
  localFibers[2] = -  std::cos( thetaFiber ) *  std::cos( thetaPositionOnSection );

  VectorSmall<3> originFibers;
  originFibers = changeOfVariable.transpose() * localFibers + positionCenter;

  originFibers.normalize();

    switch (i)
    {
        case 0:
	    // Tube
  	    // return - std::sin( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
	    return originFibers[0];

	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
   	    // return   std::cos( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
            return   originFibers[1];

	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
 	    //return   std::sin( theta );	    

	    // half thorus
            return  originFibers[2];

	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }

}

Real positionCenter ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  Real xCenter;
  Real zCenter;
  xCenter = rT * std::cos( thetaPosition ) + xT;
  zCenter = rT * std::sin( thetaPosition );

    switch (i)
    {
        case 0:
	    // Tube
	  return xCenter;
	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
            return   0.0;
	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
            return  zCenter;
	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }


}

Real localPositionSpherical ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  VectorSmall<3> centerSphere;
  centerSphere[0] = 0.138;
  centerSphere[1] = 0.0;
  centerSphere[2] = 0.185;

  VectorSmall<3> centerOnTorusAxis;
  centerOnTorusAxis[0] = rT * std::cos( thetaPosition ) + xT;
  centerOnTorusAxis[1] = 0.0;
  centerOnTorusAxis[2] = rT * std::sin( thetaPosition );

  Real xCenter;
  Real yCenter;
  Real zCenter;

  if( ( ( position - centerSphere ).norm() ) <= 0.05 ) // inside the dome
    {
      VectorSmall<3> tangentVectorToAxisTorus;
      tangentVectorToAxisTorus[0] = std::sin( thetaPosition ) ;
      tangentVectorToAxisTorus[1] = 0.0;
      tangentVectorToAxisTorus[2] = -std::cos( thetaPosition );
      
      VectorSmall<3> centerSphereOnSection;
      centerSphereOnSection = centerSphere - ( ( centerSphere - centerOnTorusAxis ).dot( tangentVectorToAxisTorus ) ) * tangentVectorToAxisTorus;

      xCenter = centerSphereOnSection[0];
      yCenter = centerSphereOnSection[1];
      zCenter = centerSphereOnSection[2];
    }
  else
    { 
      xCenter = centerOnTorusAxis[0];
      yCenter = centerOnTorusAxis[1];
      zCenter = centerOnTorusAxis[2];
    }

  // Changing reference system
  MatrixSmall<3,3> changeOfVariable;
  changeOfVariable(0,0) = std::cos( thetaRotation );  changeOfVariable(0,1) = 0.0; changeOfVariable(0,2) = std::sin( thetaRotation );
  changeOfVariable(1,0) = 0.0;                        changeOfVariable(1,1) = 1.0; changeOfVariable(1,2) = 0.0;
  changeOfVariable(2,0) = -std::sin( thetaRotation ); changeOfVariable(2,1) = 0.0; changeOfVariable(2,2) = std::cos( thetaRotation );

  VectorSmall<3> positionCenter;
  positionCenter[0] = xCenter;
  positionCenter[1] = yCenter;
  positionCenter[2] = zCenter;

  VectorSmall<3> localPosition;
  localPosition = changeOfVariable * ( position - positionCenter );


    // returning local positions
    switch (i)
    {
        case 0:
	    // Tube
	  return localPosition[0];
	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
            return   localPosition[1];
	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
            return  localPosition[2];
	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }

}

Real localPosition ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    } 

  Real xCenter;
  Real zCenter;
  xCenter = rT * std::cos( thetaPosition ) + xT;
  zCenter = rT * std::sin( thetaPosition );

  VectorSmall<3> positionCenter;
  positionCenter[0] = xCenter;
  positionCenter[1] = 0;
  positionCenter[2] = zCenter;

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  MatrixSmall<3,3> changeOfVariable;
  changeOfVariable(0,0) = std::cos( thetaRotation );  changeOfVariable(0,1) = 0.0; changeOfVariable(0,2) = std::sin( thetaRotation );
  changeOfVariable(1,0) = 0.0;                        changeOfVariable(1,1) = 1.0; changeOfVariable(1,2) = 0.0;
  changeOfVariable(2,0) = -std::sin( thetaRotation ); changeOfVariable(2,1) = 0.0; changeOfVariable(2,2) = std::cos( thetaRotation );

  VectorSmall<3> localPosition;
  localPosition = changeOfVariable * ( position - positionCenter ); 

    switch (i)
    {
        case 0:
	    // Tube
	  return localPosition[0];
	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
            return  localPosition[1];
	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
            return  localPosition[2];
	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }


}

Real sphereIndicatorFunction ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  VectorSmall<3> centerSphere;
  centerSphere[0] = 0.138;
  centerSphere[1] = 0.0;
  centerSphere[2] = 0.185;

  if( ( ( position - centerSphere ).norm() ) <= ( 0.05 ) ) // inside the dome
    {
      return 1.0;
    }
  else
    { 
      return 0.0;
    }
  return 100000.0;
}


Real positionCenterSpherical ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
  // Half thorus
  Real thetaFiber ( 0.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  VectorSmall<3> centerSphere;
  centerSphere[0] = 0.138;
  centerSphere[1] = 0.0;
  centerSphere[2] = 0.185;

  VectorSmall<3> centerOnTorusAxis;
  centerOnTorusAxis[0] = rT * std::cos( thetaPosition ) + xT;
  centerOnTorusAxis[1] = 0.0;
  centerOnTorusAxis[2] = rT * std::sin( thetaPosition );

  Real xCenter;
  Real yCenter;
  Real zCenter;

  if( ( ( position - centerSphere ).norm() ) <= 0.05 ) // inside the dome
    {
      VectorSmall<3> tangentVectorToAxisTorus;
      tangentVectorToAxisTorus[0] = std::sin( thetaPosition ) ;
      tangentVectorToAxisTorus[1] = 0.0;
      tangentVectorToAxisTorus[2] = -std::cos( thetaPosition );
      
      VectorSmall<3> centerSphereOnSection;
      centerSphereOnSection = centerSphere - ( ( centerSphere - centerOnTorusAxis ).dot( tangentVectorToAxisTorus ) ) * tangentVectorToAxisTorus;

      xCenter = centerSphereOnSection[0];
      yCenter = centerSphereOnSection[1];
      zCenter = centerSphereOnSection[2];
    }
  else
    { 
      xCenter = centerOnTorusAxis[0];
      yCenter = centerOnTorusAxis[1];
      zCenter = centerOnTorusAxis[2];
    }

    switch (i)
    {
        case 0:
	    // Tube
	  return xCenter;
	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
            return   0.0;
	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
            return  zCenter;
	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }

}

Real Family2 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{
    // Tube
    // Real theta =  0.8426; // value for anisotropic characterization taken from Robertson // ( PI / 6.0 );
    // Real thetaChangeOfVariable = std::atan(  y / x );

    // if( x < 0 )
    // {
    //     // This is due to the periodicity of std::atan ( ref. official documentation )
    //     thetaChangeOfVariable += PI;
    // }

  // Half thorus
  Real thetaFiber (  - PI / 2.0 );
  Real xT (0.138);
  Real rT ( xT );
  Real thetaFree(0);
  Real thetaRotation(0);
  Real thetaPosition(0);
  thetaFree = std::atan( z / ( x - xT ) );
  thetaPosition = thetaFree;
  if ( x < xT )
    {
      thetaPosition = thetaFree + PI;
      thetaRotation = PI/2.0 - std::fabs( thetaFree );
    }
  else
    {
      thetaRotation = thetaFree - PI/2.0;
    }

  Real xCenter;
  Real zCenter;
  xCenter = rT * std::cos( thetaPosition ) + xT;
  zCenter = rT * std::sin( thetaPosition );

  VectorSmall<3> positionCenter;
  positionCenter[0] = xCenter;
  positionCenter[1] = 0;
  positionCenter[2] = zCenter;

  VectorSmall<3> position;
  position[0] = x;
  position[1] = y;
  position[2] = z;

  MatrixSmall<3,3> changeOfVariable;
  changeOfVariable(0,0) = std::cos( thetaRotation );  changeOfVariable(0,1) = 0.0; changeOfVariable(0,2) = std::sin( thetaRotation );
  changeOfVariable(1,0) = 0.0;                        changeOfVariable(1,1) = 1.0; changeOfVariable(1,2) = 0.0;
  changeOfVariable(2,0) = -std::sin( thetaRotation ); changeOfVariable(2,1) = 0.0; changeOfVariable(2,2) = std::cos( thetaRotation );

  VectorSmall<3> localPosition;
  localPosition = changeOfVariable * ( position - positionCenter );

  Real thetaPositionOnSection(0);
  thetaPositionOnSection = std::atan( localPosition[2] / localPosition[1] );

  if ( localPosition[1] < 0 )
    thetaPositionOnSection += PI;

  VectorSmall<3> localFibers;
  localFibers[0] = std::sin( thetaFiber );
  localFibers[1] = std::cos( thetaFiber ) *  std::sin( thetaPositionOnSection );
  localFibers[2] = -  std::cos( thetaFiber ) *  std::cos( thetaPositionOnSection );

  VectorSmall<3> originFibers;
  originFibers = changeOfVariable.transpose() * localFibers + positionCenter;

  originFibers.normalize();

    switch (i)
    {
        case 0:
	    // Tube
  	    // return - std::sin( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
	    return originFibers[0];

	    // Cube
            //return std::sin( theta );
            break;
        case 1:
	    // Tube
   	    // return   std::cos( thetaChangeOfVariable ) * std::cos( theta );

	    // half thorus
            return   originFibers[1];

	    // Cube
            //return std::cos( theta );
            break;
        case 2:
	    // Tube
 	    //return   std::sin( theta );	    

	    // half thorus
            return  originFibers[2];

	    // Cube
            //return 0.0;
            break;
        default:
            ERROR_MSG ("This entrie is not allowed: ud_functions.hpp");
            return 0.;
            break;
    }
}

Real Family3 ( const Real& /*t*/, const Real& x, const Real& y, const Real& z, const ID& i)
{

    Real theta = - 0.8426;
    Real thetaChangeOfVariable = std::atan( y / x );

    if( x < 0 )
    {
        // This is due to the periodicity of std::atan ( ref. official documentation )
        thetaChangeOfVariable += PI;
    }

    switch (i)
    {
        case 0:
	    // Tube
            return - std::sin( thetaChangeOfVariable ) * std::cos( theta );
	    // Cube
            //return std::sin( theta );
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
            //return 0.0;
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
    fiberFunctionPtr_Type pointerToFunction( new fiberFunction_Type( Family1Spherical ) );
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
