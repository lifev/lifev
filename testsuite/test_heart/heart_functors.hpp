#ifndef __HEART_FUNCTORS_H
#define __HEART_FUNCTORS_H


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


namespace LifeV
{
using namespace std;

class HeartFunctors
{
public:
    HeartFunctors( GetPot& dataFile );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&, Real const&)> fct_type;
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const&)> fct_type1;

    GetPot _dataFile;
    Epetra_Comm*   comm;
    
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
    
    
    
    /**
     * current volume source
     *
     * Define the stimulation current
     */
    Real stim( const Real& x,
              const Real& y,
              const Real& z,
              const Real& t,
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
//	        switch (stim_source)
//	        {
//		        case 0:
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
		          returnvalue1=stim_value_5;
		          }
		          else returnvalue5=0.;
		        
		        if( (T_reset_6>=stim_start_6 && T_reset_6<=stim_stop_6) && 
		        	( ((x-stim_center_6(0))*(x-stim_center_6(0))+(y-stim_center_6(1))*(y-stim_center_6(1))+(z-stim_center_6(2))*(z-stim_center_6(2)))
		        			<= (stim_radius_6*stim_radius_6)) )
		          {
		          returnvalue6=stim_value_6;
		          }
		          else returnvalue6=0.;

            Real returnvalue = returnvalue1+returnvalue2+returnvalue3+returnvalue4+returnvalue5+returnvalue6;
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

};


}

#endif
