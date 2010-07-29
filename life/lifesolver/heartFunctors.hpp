/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politecnico di Milano

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*!
  \file heartFunctors.hpp

  \brief Heart Functors for the Luo-Rudy Kinetics
   \date 04/10
*/

#ifndef _HEARTFUNCTORS_H_
#define _HEARTFUNCTORS_H_


#include <life/lifecore/application.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefilters/ensight.hpp>
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <fstream>

namespace LifeV
{
using namespace std;

class HeartFunctors
{
public:
    HeartFunctors( GetPot& dataFile );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&, Real const&)> fct_type;
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&)> fct_type1;
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&,const ID&, const EntityFlag& )> fct_typeREO;
    GetPot _dataFile;
    Epetra_Comm*   comm;
	int stim_source;
    Real stim_period_1;
    Real stim_period_2;
    Real stim_period_3;
    Real stim_period_4;
    Real stim_period_5;
    Real stim_period_6;

    Real stim_start_1;
    Real stim_stop_1;
    Real stim_value_1;
    Real stim_radius_1;
    KN<Real> stim_center_1;

    Real stim_start_2;
    Real stim_stop_2;
    Real stim_value_2;
    Real stim_radius_2;
    KN<Real> stim_center_2;

    Real stim_start_3;
    Real stim_stop_3;
    Real stim_value_3;
    Real stim_radius_3;
    KN<Real> stim_center_3;

    Real stim_start_4;
    Real stim_stop_4;
    Real stim_value_4;
    Real stim_radius_4;
    KN<Real> stim_center_4;

    Real stim_start_5;
    Real stim_stop_5;
    Real stim_value_5;
    Real stim_radius_5;
    KN<Real> stim_center_5;

    Real stim_start_6;
    Real stim_stop_6;
    Real stim_value_6;
    Real stim_radius_6;
    KN<Real> stim_center_6;

    Real x_sphere;
    Real y_sphere;
    Real z_sphere;
    Real r_sphere;
    KN<Real> sigma_reduction;

    Real x_cylinder;
    Real y_cylinder;
    Real z_cylinder;
    Real a_cylinder;
    Real b_cylinder;
    Real c_cylinder;
    Real r_cylinder;
    Real xmin_cylinder;
    Real xmax_cylinder;

    Real xmin_box;
    Real ymin_box;
    Real zmin_box;

    Real xmax_box;
    Real ymax_box;
    Real zmax_box;

    //parametres Iapp IappZygote REO
    Real G_Time_period;
    Real G_Iapp_RV_angle;
    Real G_Iapp_LV_angle;
    Real G_Iapp_stim_time_RV;
    Real G_Iapp_stim_time_LV;
    Real G_Ventricular_Fibrillation;
    Real G_nb_fibrillation_sources;
    Real G_fibrillation_sources;
    Real u0;

    /**
     * current volume source
     *
     * Define the stimulation current
     */

    Real Iapp ( const Real& x, const Real& y, const Real& z, const Real& t, const EntityFlag& ) const
    {
  		Real iapp=0.0;
        double pi = acos(-1.0);
        double Angular_velocity_LV = (pi/2)/G_Iapp_stim_time_LV ; // rd/s
        double Angular_velocity_RV = (pi/2)/G_Iapp_stim_time_RV ;
        Real sumL1,sumL4,a_L,b_L,c_L,e1,e4;
   		//exitation istantane' pi=0;
        a_L = 40;
        b_L = 40;
        c_L = 72;
        e4 = 0;
        e1  = 16*0.8;

        sumL1 = x*x/( (a_L-e1)*(a_L-e1) )*100
			+  y*y/( (b_L-e1)*(b_L-e1) )*100
            +  z*z/( (c_L-e1)*(c_L-e1) )*100 -1;
        sumL4 = x*x/( (a_L-e4)*(a_L-e4) ) *100
            +  y*y/( (b_L-e4)*(b_L-e4) ) *100
            +  z*z/( (c_L-e4)*(c_L-e4) )*100 -1;
        Real sumR2,a_R,b_R,c_R,e2;
        a_R = 78;
        b_R = 40;
        c_R = 72* 0.95;
        e2  = 10*0.8;
		sumR2 = x*x/( (a_R-e2)*(a_R-e2) )*100
            +  y*y/( (b_R-e2)*(b_R-e2) )*100
            +  z*z/( (c_R-e2)*(c_R-e2) )*100 -1;
        //============ Coment if BBG
        if ( (fmod(t,G_Time_period) >= 0) &&  (fmod(t,G_Time_period)<=  G_Iapp_stim_time_LV) )
        {
            if ( sumL1 < 0.15 && ( atan(( x + 1.5e+00) /( 2.7586e+00 - z) ) < fmod(t,G_Time_period) *  Angular_velocity_LV ) )      // -0.05<sumL1< 0.15
            {
                //      if (  ( atan(( x + 1.5e+00) /( 2.7586e+00 - z) ) < pi/4  ) ) // BBG Partiel angle pi/4
                if (  ( atan(( x + 1.5e+00) /( 2.7586e+00 - z) ) < G_Iapp_LV_angle*pi/180  ) ) //
                    iapp  = (-5/2 * sumL1 + 0.375) * 8 ;
            }
        }
        //=====================septum droit
        if (   (  fmod(t,G_Time_period) >= 0)  && ( fmod(t,G_Time_period) <=  G_Iapp_stim_time_LV)  &&  (sumL4 >= -0.2) && (sumL4 <= 0.0) && ( x< -2.4) )
            iapp= (5/2 * sumL4 + 0.5) * 8 ;
        //============ Coment if BBD
        if ( fmod(t,G_Time_period)>= 0  && fmod(t,G_Time_period)<=  G_Iapp_stim_time_LV)
        {
            if (x < -2.586e+00)
                if (sumL4 > -0.2 )
                {
                    if ( (sumR2 >= 0.0) && (sumR2 < 0.1) && ( atan((- 3.7e+00 - x) /( 2.9586e+00 - z) ) < (fmod(t,G_Time_period) )* Angular_velocity_RV  ) )
                        //              if ( ( atan((- 3.7e+00 - x) /( 2.9586e+00 - z) ) < pi/6 )) //BBD Partiel angle pi/6
                        if ( ( atan((- 3.7e+00 - x) /( 2.9586e+00 - z) ) < G_Iapp_RV_angle*pi/180 )) //
                            iapp = 0.5*8 ;//( -5 * sumR2 + 0.5) * 8 ;//0.3*4.;
                }
        }
        return iapp;
    }


//! To convert function in boost functor
	inline fct_type1 get_Iapp()
    {
        fct_type1 f;
        f = boost::bind(&HeartFunctors::Iapp, this, _1, _2, _3, _4, _5 );
        return f;
    }

    Real stim( const Real& t,
               const Real& x,
               const Real& y,
               const Real& z,
               const ID&   id) const
    {
		Real returnvalue1;
		Real returnvalue2;
        Real returnvalue3;
        Real returnvalue4;
        Real returnvalue5;
        Real returnvalue6;
        Real T_reset_1(t - static_cast<int>(t/stim_period_1) * stim_period_1);
        Real T_reset_2(t - static_cast<int>(t/stim_period_2) * stim_period_2);
        Real T_reset_3(t - static_cast<int>(t/stim_period_3) * stim_period_3);
        Real T_reset_4(t - static_cast<int>(t/stim_period_4) * stim_period_4);
        Real T_reset_5(t - static_cast<int>(t/stim_period_5) * stim_period_5);
		Real T_reset_6(t - static_cast<int>(t/stim_period_6) * stim_period_6);

		if( (T_reset_1>=stim_start_1 && T_reset_1<=stim_stop_1) &&
            ( ((x-stim_center_1(0))*(x-stim_center_1(0))+(y-stim_center_1(1))*(y-stim_center_1(1))+(z-stim_center_1(2))*(z-stim_center_1(2)))
              <= (stim_radius_1*stim_radius_1)) )
        {
		    returnvalue1=stim_value_1;
        }
		else returnvalue1=0.;
		if( (T_reset_2>=stim_start_2 && T_reset_2<=stim_stop_2) &&
            ( ((x-stim_center_2(0))*(x-stim_center_2(0))+(y-stim_center_2(1))*(y-stim_center_2(1))+(z-stim_center_2(2))*(z-stim_center_2(2)))
              <= (stim_radius_2*stim_radius_2)) )
        {
		    returnvalue2=stim_value_2;
        }
		else returnvalue2=0.;

		if( (T_reset_3>=stim_start_3 && T_reset_1<=stim_stop_3) &&
            ( ((x-stim_center_3(0))*(x-stim_center_3(0))+(y-stim_center_3(1))*(y-stim_center_3(1))+(z-stim_center_3(2))*(z-stim_center_3(2)))
              <= (stim_radius_3*stim_radius_3)) )
        {
            returnvalue3=stim_value_3;
        }
		else returnvalue3=0.;
		if( (T_reset_4>=stim_start_4 && T_reset_4<=stim_stop_4) &&
            ( ((x-stim_center_4(0))*(x-stim_center_4(0))+(y-stim_center_4(1))*(y-stim_center_4(1))+(z-stim_center_4(2))*(z-stim_center_4(2)))
              <= (stim_radius_4*stim_radius_4)) )
        {
            returnvalue4=stim_value_4;
        }
		else returnvalue4=0.;
		if( (T_reset_5>=stim_start_5 && T_reset_5<=stim_stop_5) &&
            ( ((x-stim_center_5(0))*(x-stim_center_5(0))+(y-stim_center_5(1))*(y-stim_center_5(1))+(z-stim_center_5(2))*(z-stim_center_5(2)))
              <= (stim_radius_5*stim_radius_5)) )
        {
            returnvalue5=stim_value_5;
        }
		else returnvalue5=0.;
		if( (T_reset_6>=stim_start_6 && T_reset_6<=stim_stop_6) &&
            ( ((x-stim_center_6(0))*(x-stim_center_6(0))+(y-stim_center_6(1))*(y-stim_center_6(1))+(z-stim_center_6(2))*(z-stim_center_6(2)))
              <= (stim_radius_6*stim_radius_6)) )
        {
            returnvalue6=stim_value_6;
        }
		else returnvalue6=0.;
        Real	returnvalue = returnvalue1+returnvalue2+returnvalue3+returnvalue4+returnvalue5+returnvalue6;
        if(id == 0)
	        return returnvalue;
        else return returnvalue;
    }

    //! To convert function in boost functor
    inline fct_type1 get_stim()
    {
        fct_type1 f;
        f = boost::bind(&HeartFunctors::stim, this, _1, _2, _3, _4, _5);
        return f;
    }

    Real IappZygote(const double& t,
                    const double& x,
                    const double& y,
                    const double& z,
                    const ID& i,
                    const EntityFlag& ref)
	{
        double pi = acos(-1.0);
    	Real iapp=0.0;
     	Real x_0= 3.316424,y_0 =  31.496351,z_0 = 5.799850;//APEX, node number: 80185 : 3.316424 31.496351 5.799850
     	if (fmod(t,G_Time_period)<=25)
       	{
            if (ref==2) // ventricule gauche if (atan( d1/d2) < pi* t/20 )
            {
                if ( ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) + (z-z_0)*(z-z_0)>=(fmod(t,G_Time_period))*(fmod(t,G_Time_period)) - 100 )   &&  ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) + (z-z_0)*(z-z_0)<=3*(fmod(t,G_Time_period))*(fmod(t,G_Time_period)) ))
                    iapp=2;
            }
            if ((ref==1)||(ref==20)) // if (ref==20) //BBD
            {
                //venticule droit    if (atan( d1/d2) < pi* t/20 )
	    		if ( ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) + (z-z_0)*(z-z_0)>=(fmod(t,G_Time_period))*(fmod(t,G_Time_period)) - 100 )   &&  ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) + (z-z_0)*(z-z_0)<=3*(fmod(t,G_Time_period))*(fmod(t,G_Time_period)) ))
                    iapp=2 ;
            }
        }
        return iapp;
    }
    //! To convert function in boost functor
    inline fct_typeREO get_IappZygote()
    {
        fct_typeREO f;
        f = boost::bind(&HeartFunctors::IappZygote, this, _1, _2, _3, _4, _5, _6);
        return f;
    }

    /**
     *
     * Reduces the conductivity in a sphere
     *
     */
    Real reduced_sigma_sphere( const Real& x,
                               const Real& y,
                               const Real& z,
                               const Real& /*t*/,
                               const ID&   id,
                               const Real& sigma) const
    {
        if (((x-x_sphere)*(x-x_sphere)+(y-y_sphere)*(y-y_sphere)+(z-z_sphere)*(z-z_sphere))<r_sphere*r_sphere)
            return sigma*sigma_reduction(id);
        else return sigma;
    }
    //! To convert function in boost functor
    inline const fct_type get_reduced_sigma_sphere()
    {
        fct_type f;
        f = boost::bind(&HeartFunctors::reduced_sigma_sphere, this, _1, _2, _3, _4, _5, _6);
        return f;
    }

    /**
     *
     * Reduces the conductivity in a cylinder
     *
     */
    Real reduced_sigma_cylinder( const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const Real& /*t*/,
                                 const ID&   id,
                                 const Real& sigma ) const
    {
    	Real distance2, distx, disty, distz;
    	distx=((b_cylinder*b_cylinder+c_cylinder*c_cylinder)*(x_cylinder-x)-a_cylinder*c_cylinder*(z_cylinder-z)-a_cylinder*b_cylinder*(y_cylinder-y))/(a_cylinder*a_cylinder+b_cylinder*b_cylinder+c_cylinder*c_cylinder);
    	disty=((a_cylinder*a_cylinder+c_cylinder*c_cylinder)*(y_cylinder-y)-b_cylinder*c_cylinder*(z_cylinder-z)-a_cylinder*b_cylinder*(x_cylinder-x))/(a_cylinder*a_cylinder+b_cylinder*b_cylinder+c_cylinder*c_cylinder);
    	distz=((a_cylinder*a_cylinder+b_cylinder*b_cylinder)*(z_cylinder-z)-a_cylinder*c_cylinder*(x_cylinder-x)-c_cylinder*b_cylinder*(y_cylinder-y))/(a_cylinder*a_cylinder+b_cylinder*b_cylinder+c_cylinder*c_cylinder);
    	distance2=pow(distx,2)+pow(disty,2)+pow(distz,2);
    	if ((distance2<r_cylinder*r_cylinder)&&(x<xmax_cylinder)&&(x>xmin_cylinder))
        {
            return sigma*sigma_reduction(id);
        }
    	else return sigma;
    }
    inline const fct_type get_reduced_sigma_cylinder()
    {
        fct_type f;
        f = boost::bind(&HeartFunctors::reduced_sigma_cylinder, this, _1, _2, _3, _4, _5, _6);
        return f;
    }

    /**
     *
     * Reduces the conductivity in a box
     *
     */
    Real reduced_sigma_box( const Real& x,
                            const Real& y,
                            const Real& z,
                            const Real& /*t*/,
                            const ID&   id,
                            const Real& sigma ) const
    {
    	if  ((x>xmin_box)&&(x<xmax_box)
             &&(y>ymin_box)&&(y<ymax_box)
             &&(z>zmin_box)&&(z<zmax_box))
        {
            return sigma*sigma_reduction(id);
        }
    	else return sigma;
    }

    inline const fct_type get_reduced_sigma_box()
    {
        fct_type f;
        f = boost::bind(&HeartFunctors::reduced_sigma_box, this, _1, _2, _3, _4, _5, _6);
        return f;
    }


    Real initial_scalar( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID& /* i */ )
    {
        return u0;
    }
    inline const fct_type1 get_initial_scalar()
    {
        fct_type1 f;
        f = boost::bind(&HeartFunctors::initial_scalar, this, _1, _2, _3, _4, _5);
        return f;
    }

    Real zero_scalar( const Real& /* t */,
            const Real& /* x */,
            const Real& /* y */,
            const Real& /* z */,
            const ID& /* i */ )
    {
        return 0.;
    }
    inline const fct_type1 get_zero_scalar()
    {
        fct_type1 f;
        f = boost::bind(&HeartFunctors::zero_scalar, this, _1, _2, _3, _4, _5);
        return f;
    }


};
}

#endif
