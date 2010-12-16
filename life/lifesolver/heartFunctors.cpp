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
    @brief  Heart Functors for the Luo-Rudy Kinetics

    @date 04−2010
    @author

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */




#include <life/lifesolver/heartFunctors.hpp>

using namespace LifeV;

// ===================================================
// Constructors & Destructor
// ===================================================
HeartFunctors::HeartFunctors():
    M_dataFile          ( ),
    M_stimulusSource    ( ),
    M_stimulusPeriod1   ( ),
    M_stimulusPeriod2   ( ),
    M_stimulusPeriod3   ( ),
    M_stimulusPeriod4   ( ),
    M_stimulusPeriod5   ( ),
    M_stimulusPeriod6   ( ),
    M_stimulusStart1    ( ),
    M_stimulusStop1     ( ),
    M_stimulusValue1    ( ),
    M_stimulusRadius1   ( ),
    M_stimulusCenter1   (3),
    M_stimulusStart2    ( ),
    M_stimulusStop2     ( ),
    M_stimulusValue2    ( ),
    M_stimulusRadius2   ( ),
    M_stimulusCenter2   (3),
    M_stimulusStart3    ( ),
    M_stimulusStop3     ( ),
    M_stimulusValue3    ( ),
    M_stimulusRadius3   ( ),
    M_stimulusCenter3   (3),
    M_stimulusStart4    ( ),
    M_stimulusStop4     ( ),
    M_stimulusValue4    ( ),
    M_stimulusRadius4   ( ),
    M_stimulusCenter4   (3),
    M_stimulusStart5    ( ),
    M_stimulusStop5     ( ),
    M_stimulusValue5    ( ),
    M_stimulusRadius5   ( ),
    M_stimulusCenter5   (3),
    M_stimulusStart6    ( ),
    M_stimulusStop6     ( ),
    M_stimulusValue6    ( ),
    M_stimulusRadius6   ( ),
    M_stimulusCenter6   (3),
    M_sphereX           ( ),
    M_sphereY           ( ),
    M_sphereZ           ( ),
    M_sphereR           ( ),
    M_sigmaReduction    (2),
    M_cylinderX         ( ),
    M_cylinderY         ( ),
    M_cylinderZ         ( ),
    M_cylinderA         ( ),
    M_cylinderB         ( ),
    M_cylinderC         ( ),
    M_cylinderR         ( ),
    M_minimumCylinderX  ( ),
    M_maximumCylinderX  ( ),
    M_minimumBoxX       ( ),
    M_minimumBoxY       ( ),
    M_minimumBoxZ       ( ),
    M_maximumBoxX       ( ),
    M_maximumBoxY       ( ),
    M_maximumBoxZ       ( ),
    //parametre Iapp IappZygote
    M_timePeriod                                ( ),
    M_appliedCurrentRightVentriculeAngle        ( ),
    M_appliedCurrentLeftVentriculeAngle         ( ),
    M_appliedCurrentStimulusTimeRightVentricule ( ),
    M_appliedCurrentStimulusTimeLeftVentrivcule ( ),
    M_ventricularFibrillation                   ( ),
    M_nb_fibrillationSources                    ( ),
    M_fibrillationSources                       ( ),
    M_restPotential                             ( )
{
}


HeartFunctors::HeartFunctors( GetPot& dataFile ):
    M_dataFile          (dataFile),
    M_stimulusSource    (dataFile("electric/physics/stim_source",1)),
    M_stimulusPeriod1   (dataFile("electric/physics/stim_period_1",200.)),
    M_stimulusPeriod2   (dataFile("electric/physics/stim_period_2",200.)),
    M_stimulusPeriod3   (dataFile("electric/physics/stim_period_3",200.)),
    M_stimulusPeriod4   (dataFile("electric/physics/stim_period_4",200.)),
    M_stimulusPeriod5   (dataFile("electric/physics/stim_period_5",200.)),
    M_stimulusPeriod6   (dataFile("electric/physics/stim_period_6",200.)),
    M_stimulusStart1    (dataFile("electric/physics/stim_start_1",0.)),
    M_stimulusStop1     (dataFile("electric/physics/stim_stop_1",0.)),
    M_stimulusValue1    (dataFile("electric/physics/stim_value_1",0.)),
    M_stimulusRadius1   (dataFile("electric/physics/stim_radius_1",0.)),
    M_stimulusCenter1   (3),
    M_stimulusStart2    (dataFile("electric/physics/stim_start_2",0.)),
    M_stimulusStop2     (dataFile("electric/physics/stim_stop_2",0.)),
    M_stimulusValue2    (dataFile("electric/physics/stim_value_2",0.)),
    M_stimulusRadius2   (dataFile("electric/physics/stim_radius_2",0.)),
    M_stimulusCenter2   (3),
    M_stimulusStart3    (dataFile("electric/physics/stim_start_3",0.)),
    M_stimulusStop3     (dataFile("electric/physics/stim_stop_3",0.)),
    M_stimulusValue3    (dataFile("electric/physics/stim_value_3",0.)),
    M_stimulusRadius3   (dataFile("electric/physics/stim_radius_3",0.)),
    M_stimulusCenter3   (3),
    M_stimulusStart4    (dataFile("electric/physics/stim_start_4",0.)),
    M_stimulusStop4     (dataFile("electric/physics/stim_stop_4",0.)),
    M_stimulusValue4    (dataFile("electric/physics/stim_value_4",0.)),
    M_stimulusRadius4   (dataFile("electric/physics/stim_radius_4",0.)),
    M_stimulusCenter4   (3),
    M_stimulusStart5    (dataFile("electric/physics/stim_start_5",0.)),
    M_stimulusStop5     (dataFile("electric/physics/stim_stop_5",0.)),
    M_stimulusValue5    (dataFile("electric/physics/stim_value_5",0.)),
    M_stimulusRadius5   (dataFile("electric/physics/stim_radius_5",0.)),
    M_stimulusCenter5   (3),
    M_stimulusStart6    (dataFile("electric/physics/stim_start_6",0.)),
    M_stimulusStop6     (dataFile("electric/physics/stim_stop_6",0.)),
    M_stimulusValue6    (dataFile("electric/physics/stim_value_6",0.)),
    M_stimulusRadius6   (dataFile("electric/physics/stim_radius_6",0.)),
    M_stimulusCenter6   (3),
    M_sphereX           (dataFile("electric/physics/sphere_center",0.,0)),
    M_sphereY           (dataFile("electric/physics/sphere_center",0.,1)),
    M_sphereZ           (dataFile("electric/physics/sphere_center",0.,2)),
    M_sphereR           (dataFile("electric/physics/sphere_radius",0.)),
    M_sigmaReduction    (2),
    M_cylinderX         (dataFile("electric/physics/x_cylinder",0.)),
    M_cylinderY         (dataFile("electric/physics/y_cylinder",0.)),
    M_cylinderZ         (dataFile("electric/physics/z_cylinder",0.)),
    M_cylinderA         (dataFile("electric/physics/a_cylinder",0.)),
    M_cylinderB         (dataFile("electric/physics/b_cylinder",0.)),
    M_cylinderC         (dataFile("electric/physics/c_cylinder",0.)),
    M_cylinderR         (dataFile("electric/physics/r_cylinder",0.)),
    M_minimumCylinderX  (dataFile("electric/physics/xmin_cylinder",0.)),
    M_maximumCylinderX  (dataFile("electric/physics/xmax_cylinder",0.)),
    M_minimumBoxX       (dataFile("electric/physics/box_vertex_min",0.,0)),
    M_minimumBoxY       (dataFile("electric/physics/box_vertex_min",0.,1)),
    M_minimumBoxZ       (dataFile("electric/physics/box_vertex_min",0.,2)),
    M_maximumBoxX       (dataFile("electric/physics/box_vertex_max",0.,0)),
    M_maximumBoxY       (dataFile("electric/physics/box_vertex_max",0.,1)),
    M_maximumBoxZ       (dataFile("electric/physics/box_vertex_max",0.,2)),
    //parametre Iapp IappZygote
    M_timePeriod                                (dataFile("electric/physics/Time_period",700.0)),
    M_appliedCurrentRightVentriculeAngle        (dataFile("electric/physics/Iapp_RV_angle",360.)),
    M_appliedCurrentLeftVentriculeAngle         (dataFile("electric/physics/Iapp_LV_angle",360.)),
    M_appliedCurrentStimulusTimeRightVentricule (dataFile("electric/physics/Iapp_stim_time_RV",6.)),
    M_appliedCurrentStimulusTimeLeftVentrivcule (dataFile("electric/physics/Iapp_stim_time_LV",10.)),
    M_ventricularFibrillation                   (dataFile("electric/physics/Ventricular_Fibrillation",0)),
    M_nb_fibrillationSources                    (dataFile("electric/physics/nb_fibrillation_sources",20)),
    M_fibrillationSources                       (dataFile("electric/physics/fibrillation_sources",0)),
    M_restPotential                             (dataFile("electric/physics/u0",0.))
{
    M_sigmaReduction(0)  = dataFile("electric/physics/sigma_reduction",1.,0);
    M_sigmaReduction(1)  = dataFile("electric/physics/sigma_reduction",1.,1);
    M_stimulusCenter1(0) = dataFile("electric/physics/stim_center_1",0.,0);
    M_stimulusCenter1(1) = dataFile("electric/physics/stim_center_1",0.,1);
    M_stimulusCenter1(2) = dataFile("electric/physics/stim_center_1",0.,2);
    M_stimulusCenter2(0) = dataFile("electric/physics/stim_center_2",0.,0);
    M_stimulusCenter2(1) = dataFile("electric/physics/stim_center_2",0.,1);
    M_stimulusCenter2(2) = dataFile("electric/physics/stim_center_2",0.,2);
    M_stimulusCenter3(0) = dataFile("electric/physics/stim_center_3",0.,0);
    M_stimulusCenter3(1) = dataFile("electric/physics/stim_center_3",0.,1);
    M_stimulusCenter3(2) = dataFile("electric/physics/stim_center_3",0.,2);
    M_stimulusCenter4(0) = dataFile("electric/physics/stim_center_4",0.,0);
    M_stimulusCenter4(1) = dataFile("electric/physics/stim_center_4",0.,1);
    M_stimulusCenter4(2) = dataFile("electric/physics/stim_center_4",0.,2);
    M_stimulusCenter5(0) = dataFile("electric/physics/stim_center_5",0.,0);
    M_stimulusCenter5(1) = dataFile("electric/physics/stim_center_5",0.,1);
    M_stimulusCenter5(2) = dataFile("electric/physics/stim_center_5",0.,2);
    M_stimulusCenter6(0) = dataFile("electric/physics/stim_center_6",0.,0);
    M_stimulusCenter6(1) = dataFile("electric/physics/stim_center_6",0.,1);
    M_stimulusCenter6(2) = dataFile("electric/physics/stim_center_6",0.,2);
}

// ===================================================
// Set Methods
// ===================================================
Real
HeartFunctors::Iapp ( const Real& x, const Real& y, const Real& z, const Real& t, const entityFlag_Type& id ) const
{
    Real appliedCurrent                 = 0.0;
    Real pi                             = std::acos( -1.0 );
    Real AngularVelocityLeftVentricule  = ( pi / 2 ) / M_appliedCurrentStimulusTimeLeftVentrivcule ; // rd/s
    Real AngularVelocityRightVentricule = ( pi / 2 ) / M_appliedCurrentStimulusTimeRightVentricule ;
    Real sumL1,sumL4,aL,bL,cL,e1,e4;
    //exitation istantane' pi=0;
    aL = 40;
    bL = 40;
    cL = 72;
    e4 = 0;
    e1  = 16*0.8;

    sumL1 = x * x / ( ( aL - e1 ) * ( aL - e1 ) ) * 100 +
            y * y / ( ( bL - e1 ) * ( bL - e1 ) ) * 100 +
            z * z / ( ( cL - e1 ) * ( cL - e1 ) ) * 100 -1;

    sumL4 = x * x / ( ( aL - e4 ) * ( aL - e4 ) ) * 100 +
            y * y / ( ( bL - e4 ) * ( bL - e4 ) ) * 100 +
            z * z / ( ( cL - e4 ) * ( cL - e4 ) ) * 100 -1;

    Real sumR2,aR,bR,cR,e2;
    aR = 78;
    bR = 40;
    cR = 72* 0.95;
    e2  = 10*0.8;

    sumR2 = x * x / ( ( aR - e2 ) * ( aR - e2 ) ) * 100 +
            y * y / ( ( bR - e2 ) * ( bR - e2 ) ) * 100 +
            z * z / ( ( cR - e2 ) * ( cR - e2 ) ) * 100 -1;

    //============ Coment if BBG
    if ( fmod(t, M_timePeriod) >= 0 &&
         fmod(t, M_timePeriod) <=  M_appliedCurrentStimulusTimeLeftVentrivcule )
    {
        if ( sumL1 < 0.15 &&
           ( std::atan( ( x + 1.5e+00) / ( 2.7586e+00 - z) ) < fmod(t, M_timePeriod) *  AngularVelocityLeftVentricule ) )      // -0.05<sumL1< 0.15
        {
            //      if (  ( atan(( x + 1.5e+00) /( 2.7586e+00 - z) ) < pi/4  ) ) // BBG Partiel angle pi/4
            if (  ( std::atan( ( x + 1.5e+00) / ( 2.7586e+00 - z) ) < M_appliedCurrentLeftVentriculeAngle * pi / 180  ) ) //
                appliedCurrent  = ( -5 / 2 * sumL1 + 0.375) * 8 ;
        }
    }

    //=====================septum droit
    if ( fmod(t, M_timePeriod) >= 0 &&
         fmod(t, M_timePeriod) <=  M_appliedCurrentStimulusTimeLeftVentrivcule &&
         sumL4 >= -0.2 && sumL4 <= 0.0 && x < -2.4 )
                appliedCurrent= (5 / 2 * sumL4 + 0.5) * 8 ;

    //============ Coment if BBD
    if ( fmod(t, M_timePeriod)>= 0  &&
         fmod(t, M_timePeriod) <=  M_appliedCurrentStimulusTimeLeftVentrivcule )
    {
        if (x < -2.586e+00)
            if (sumL4 > -0.2 )
            {
                if ( sumR2 >= 0.0 && sumR2 < 0.1 &&
                     std::atan( (- 3.7e+00 - x ) / ( 2.9586e+00 - z) ) < fmod(t,M_timePeriod) * AngularVelocityRightVentricule )
                    //              if ( ( atan((- 3.7e+00 - x) /( 2.9586e+00 - z) ) < pi/6 )) //BBD Partiel angle pi/6
                    if ( std::atan( (- 3.7e+00 - x) / ( 2.9586e+00 - z) ) < M_appliedCurrentRightVentriculeAngle * pi / 180  ) //
                            appliedCurrent = 0.5 * 8 ;//( -5 * sumR2 + 0.5) * 8 ;//0.3*4.;
            }
    }

    return appliedCurrent;

}


Real
HeartFunctors::IappZygote(const double& t,
                          const double& x,
                          const double& y,
                          const double& z,
                          const ID& i,
                          const entityFlag_Type& ref)
{
    // double pi = acos(-1.0);
    Real appliedCurrent = 0.0;
    Real x0 = 3.316424;
    Real y0 = 31.496351;
    Real z0 = 5.799850; //APEX, node number: 80185 : 3.316424 31.496351 5.799850

    if ( fmod(t,M_timePeriod) <= 25 )
    {
        if ( ref == 2 ) // ventricule gauche if (atan( d1/d2) < pi* t/20 )
        {
            if ( ( ( x - x0 ) * ( x - x0 ) +
                   ( y - y0 ) * ( y - y0 ) +
                   ( z - z0 ) * ( z - z0 ) >= fmod(t,M_timePeriod) * fmod(t,M_timePeriod) - 100 )
                            &&
                 ( ( x - x0 ) * ( x - x0 ) +
                   ( y - y0 ) * ( y - y0 ) +
                   ( z - z0 ) * ( z - z0 ) <= 3 * fmod(t,M_timePeriod) * fmod(t,M_timePeriod) ) )
                appliedCurrent = 2;
        }
        if ( ( ref == 1 ) || ( ref == 20 ) ) // if (ref==20) //BBD
        {
            //venticule droit    if (atan( d1/d2) < pi* t/20 )
            if ( ( ( x - x0 ) * ( x - x0 ) +
                   ( y - y0 ) * ( y - y0 ) +
                   ( z - z0 ) * ( z - z0 ) >= fmod(t,M_timePeriod) * fmod(t,M_timePeriod) - 100 )
                            &&
                 ( ( x - x0 ) * ( x - x0 ) +
                   ( y - y0 ) * ( y - y0 ) +
                   ( z - z0 ) * ( z - z0 ) <= 3 * fmod(t,M_timePeriod) * fmod(t,M_timePeriod) ) )
                appliedCurrent = 2;
        }
    }

    return appliedCurrent;
}


Real
HeartFunctors::stim( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const
{
    Real returnValue1;
    Real returnValue2;
    Real returnValue3;
    Real returnValue4;
    Real returnValue5;
    Real returnValue6;
    Real timeReset1(t - static_cast<int>(t/M_stimulusPeriod1) * M_stimulusPeriod1);
    Real timeReset2(t - static_cast<int>(t/M_stimulusPeriod2) * M_stimulusPeriod2);
    Real timeReset3(t - static_cast<int>(t/M_stimulusPeriod3) * M_stimulusPeriod3);
    Real timeReset4(t - static_cast<int>(t/M_stimulusPeriod4) * M_stimulusPeriod4);
    Real timeReset5(t - static_cast<int>(t/M_stimulusPeriod5) * M_stimulusPeriod5);
    Real timeReset6(t - static_cast<int>(t/M_stimulusPeriod6) * M_stimulusPeriod6);

    if ( (timeReset1 >= M_stimulusStart1 && timeReset1 <= M_stimulusStop1) &&
            ( ( ( x - M_stimulusCenter1(0) ) * ( x - M_stimulusCenter1(0) ) +
                ( y - M_stimulusCenter1(1) ) * ( y - M_stimulusCenter1(1) ) +
                ( z - M_stimulusCenter1(2) ) * ( z - M_stimulusCenter1(2) ) ) <= ( M_stimulusRadius1 * M_stimulusRadius1 ) ) )
    {
        returnValue1 = M_stimulusValue1;
    }
    else returnValue1 = 0.;

    if ( (timeReset2 >= M_stimulusStart2 && timeReset2 <= M_stimulusStop2) &&
            ( ( ( x - M_stimulusCenter2(0) ) * ( x - M_stimulusCenter2(0) ) +
                ( y - M_stimulusCenter2(1) ) * ( y - M_stimulusCenter2(1) ) +
                ( z - M_stimulusCenter2(2) ) * ( z - M_stimulusCenter2(2) ) ) <= ( M_stimulusRadius2 * M_stimulusRadius2 ) ) )
    {
        returnValue2 = M_stimulusValue2;
    }
    else returnValue2 = 0.;

    if ( (timeReset3 >= M_stimulusStart3 && timeReset3 <= M_stimulusStop3) &&
            ( ( ( x - M_stimulusCenter3(0) ) * ( x - M_stimulusCenter3(0) ) +
                ( y - M_stimulusCenter3(1) ) * ( y - M_stimulusCenter3(1) ) +
                ( z - M_stimulusCenter3(2) ) * ( z - M_stimulusCenter3(2) ) ) <= ( M_stimulusRadius3 * M_stimulusRadius3 ) ) )
    {
        returnValue3 = M_stimulusValue3;
    }
    else returnValue3 = 0.;

    if ( (timeReset4 >= M_stimulusStart4 && timeReset4 <= M_stimulusStop4) &&
            ( ( ( x - M_stimulusCenter4(0) ) * ( x - M_stimulusCenter4(0) ) +
                ( y - M_stimulusCenter4(1) ) * ( y - M_stimulusCenter4(1) ) +
                ( z - M_stimulusCenter4(2) ) * ( z - M_stimulusCenter4(2) ) ) <= ( M_stimulusRadius4 * M_stimulusRadius4 ) ) )
    {
        returnValue4 = M_stimulusValue4;
    }
    else returnValue4 = 0.;

    if ( (timeReset5 >= M_stimulusStart5 && timeReset5 <= M_stimulusStop5) &&
            ( ( ( x - M_stimulusCenter5(0) ) * ( x - M_stimulusCenter5(0) ) +
                ( y - M_stimulusCenter5(1) ) * ( y - M_stimulusCenter5(1) ) +
                ( z - M_stimulusCenter5(2) ) * ( z - M_stimulusCenter5(2) ) ) <= ( M_stimulusRadius5 * M_stimulusRadius5 ) ) )
    {
        returnValue5 = M_stimulusValue5;
    }
    else returnValue5 = 0.;

    if ( (timeReset6 >= M_stimulusStart6 && timeReset6 <= M_stimulusStop6) &&
            ( ( ( x - M_stimulusCenter6(0) ) * ( x - M_stimulusCenter6(0) ) +
                ( y - M_stimulusCenter6(1) ) * ( y - M_stimulusCenter6(1) ) +
                ( z - M_stimulusCenter6(2) ) * ( z - M_stimulusCenter6(2) ) ) <= ( M_stimulusRadius6 * M_stimulusRadius6 ) ) )
    {
        returnValue6 = M_stimulusValue6;
    }
    else returnValue6 = 0.;

    Real    returnValue = returnValue1 + returnValue2 + returnValue3 + returnValue4 + returnValue5 + returnValue6;

    if (id == 0)
        return returnValue;
    else return returnValue;

}








Real
HeartFunctors::reduced_sigma_sphere( const Real& x,
                                     const Real& y,
                                     const Real& z,
                                     const Real& /*t*/,
                                     const ID&   id,
                                     const Real& sigma) const
{
    if ( ( ( x - M_sphereX ) * ( x - M_sphereX ) +
           ( y - M_sphereY ) * ( y - M_sphereY ) +
           ( z - M_sphereZ ) * ( z - M_sphereZ ) ) < M_sphereR * M_sphereR )
        return sigma * M_sigmaReduction(id);
    else return sigma;
}



Real
HeartFunctors::reduced_sigma_cylinder( const Real& x,
                                       const Real& y,
                                       const Real& z,
                                       const Real& t,
                                       const ID&   id,
                                       const Real& sigma ) const
    {
        Real distance2, distanceX, distanceY, distanceZ;
        distanceX = ( ( M_cylinderB * M_cylinderB + M_cylinderC * M_cylinderC ) * ( M_cylinderX - x ) -
                      (                             M_cylinderA * M_cylinderC ) * ( M_cylinderZ - z ) -
                      (                             M_cylinderA * M_cylinderB ) * ( M_cylinderY - y ) ) /
                    (   M_cylinderA * M_cylinderA + M_cylinderB * M_cylinderB + M_cylinderC * M_cylinderC );

        distanceY = ( ( M_cylinderA * M_cylinderA + M_cylinderC * M_cylinderC ) * ( M_cylinderY - y ) -
                      (                             M_cylinderB * M_cylinderC ) * ( M_cylinderZ - z ) -
                      (                             M_cylinderA * M_cylinderB ) * ( M_cylinderX - x ) ) /
                    (   M_cylinderA * M_cylinderA + M_cylinderB * M_cylinderB + M_cylinderC * M_cylinderC );

        distanceZ = ( ( M_cylinderA * M_cylinderA + M_cylinderB * M_cylinderB ) * ( M_cylinderZ - z ) -
                      (                             M_cylinderA * M_cylinderC ) * ( M_cylinderX - x ) -
                      (                             M_cylinderC * M_cylinderB ) * ( M_cylinderY - y ) ) /
                    (   M_cylinderA * M_cylinderA + M_cylinderB * M_cylinderB + M_cylinderC * M_cylinderC );

        distance2 = pow(distanceX,2) + pow(distanceY,2) + pow(distanceZ,2);

        if ( ( distance2 < M_cylinderR * M_cylinderR ) && ( x < M_maximumCylinderX ) && ( x > M_minimumCylinderX ) )
        {
            return sigma * M_sigmaReduction(id);
        }
        else return sigma;
    }




Real
HeartFunctors::reduced_sigma_box( const Real& x,
                                  const Real& y,
                                  const Real& z,
                                  const Real& t,
                                  const ID&   id,
                                  const Real& sigma ) const
{
    if  ( ( x > M_minimumBoxX ) && ( x < M_maximumBoxX ) &&
          ( y > M_minimumBoxY ) && ( y < M_maximumBoxY ) &&
          ( z > M_minimumBoxZ ) && ( z < M_maximumBoxZ ) )
    {
        return sigma * M_sigmaReduction(id);
    }
    else return sigma;
}

// ===================================================
// Get Methods
// ===================================================
Real
HeartFunctors::initial_scalar( const Real& t ,
                               const Real& x,
                               const Real& y,
                               const Real& z,
                               const ID&   id )
{
    return M_restPotential;
}



Real
HeartFunctors::zero_scalar( const Real& t,
                            const Real& x,
                            const Real& y,
                            const Real& z,
                            const ID&   id )
{
    return 0.;
}


HeartFunctors::region1_Type
HeartFunctors::get_Iapp()
{
    region1_Type f;
    f = boost::bind(&HeartFunctors::Iapp, this, _1, _2, _3, _4, _5 );
    return f;
}

HeartFunctors::region1_Type
HeartFunctors::get_stim()
{
    region1_Type f;
    f = boost::bind(&HeartFunctors::stim, this, _1, _2, _3, _4, _5);
    return f;
}

const HeartFunctors::region_Type
HeartFunctors::get_reduced_sigma_sphere()
{
    region_Type f;
    f = boost::bind(&HeartFunctors::reduced_sigma_sphere, this, _1, _2, _3, _4, _5, _6);
    return f;
}

const HeartFunctors::region_Type
HeartFunctors::get_reduced_sigma_cylinder()
{
    region_Type f;
    f = boost::bind(&HeartFunctors::reduced_sigma_cylinder, this, _1, _2, _3, _4, _5, _6);
    return f;
}

const HeartFunctors::region_Type
HeartFunctors::get_reduced_sigma_box()
{
    region_Type f;
    f = boost::bind(&HeartFunctors::reduced_sigma_box, this, _1, _2, _3, _4, _5, _6);
    return f;
}

const HeartFunctors::region1_Type
HeartFunctors::get_initial_scalar()
{
    region1_Type f;
    f = boost::bind(&HeartFunctors::initial_scalar, this, _1, _2, _3, _4, _5);
    return f;
}

const HeartFunctors::region1_Type
HeartFunctors::get_zero_scalar()
{
    region1_Type f;
    f = boost::bind(&HeartFunctors::zero_scalar, this, _1, _2, _3, _4, _5);
    return f;
}
