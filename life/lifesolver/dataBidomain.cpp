/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file dataBidomain.cpp

*/


#include<life/lifesolver/dataBidomain.hpp>

namespace LifeV
{



//
// IMPLEMENTATION
//


// Constructors
#ifdef REO_CASE
//template <typename Mesh>
DataBidomain::
DataBidomain( boost::shared_ptr<HeartCaseBase> B_fct ) :
    DataMesh( B_fct->get_data_hdl(),  "electric/space_discretization" ),
    DataTime( B_fct->get_data_hdl() , "electric/time_discretization" )
{
    setup(B_fct->get_data_hdl());
}
#else
//template <typename Mesh>
DataBidomain::
DataBidomain( boost::shared_ptr<HeartFunctors> heart_fct ) :
    DataMesh( heart_fct->_dataFile, "electric/space_discretization" ),
    DataTime( heart_fct->_dataFile, "electric/time_discretization" ),
    red_sigma_sphere(heart_fct->get_reduced_sigma_sphere() ),
    red_sigma_cyl(heart_fct->get_reduced_sigma_cylinder() ),
    red_sigma_box(heart_fct->get_reduced_sigma_box() )
{
    setup(heart_fct->_dataFile);
}
#endif

//template <typename Mesh>
DataBidomain::
DataBidomain( const DataBidomain& dataBidomain ) :
    DataMesh                     ( dataBidomain ),
    DataTime                     ( dataBidomain ),
    M_uOrder(dataBidomain.M_uOrder),
    M_fibers_file(dataBidomain.M_fibers_file),
    M_Chi(dataBidomain.M_Chi),
    M_Cm(dataBidomain.M_Cm),
    M_D_i(dataBidomain.M_D_i),
    M_D_e(dataBidomain.M_D_e),
    M_sigmal_i(dataBidomain.M_sigmal_i),
    M_sigmat_i(dataBidomain.M_sigmat_i),
    M_sigmal_e(dataBidomain.M_sigmal_e),
    M_sigmat_e(dataBidomain.M_sigmat_e),
    M_heart_diff_fct(dataBidomain.M_heart_diff_fct),
    M_verbose(dataBidomain.M_verbose),
    M_post_dir(dataBidomain.M_post_dir),
    M_has_fibers(dataBidomain.M_has_fibers),
    M_fibers_format(dataBidomain.M_fibers_format),
    M_order_bdf(dataBidomain.M_order_bdf)
{
}




//template <typename Mesh>
void
DataBidomain::
setup(  const GetPot& dfile )
{
	M_Chi = dfile("electric/physics/Chi",1e3); 	// [1e-3 1/cm]   ColliPavarinoTaccardi2005
	M_Cm = dfile("electric/physics/Cm",1e-3);  	// [1e-3 mF/cm2]   ColliPavarinoTaccardi2005
    	if(dfile("electric/physics/ion_model",1) == 1)
		{
    		M_D_i = dfile("electric/physics/D_i" ,3.3e-2);       // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
    		M_D_e = dfile("electric/physics/D_e" ,4.29e-2);      // 4.29e-2 [1/Ohm/cm]	D_e_LR * D_RM/D_LR  see dataMonodomain
    		M_sigmal_i =  dfile("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
    		M_sigmat_i   = dfile("electric/physics/sigmat_i", 8.6e-3); 	// 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
    		M_sigmal_e =  dfile("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]	sigmal_e_LR * D_RM/D_LR
    		M_sigmat_e   = dfile("electric/physics/sigmat_e",3.69e-2); // 3.69e-2 [1/Ohm/cm]	sigmat_e_LR * D_RM/D_LR
		}
	else if(dfile("electric/physics/ion_model",1) == 2)
		{
		M_D_i = dfile("electric/physics/D_i" , 1.21e-3);  		// sigmal_i/3 + sigmat_i*2/3
		M_D_e = dfile("electric/physics/D_e" , 1.57e-3);  		// sigmal_e/3 + sigmat_e*2/3
		M_sigmal_i =  dfile("electric/physics/sigmal_i", 3e-3);  		// 3e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmat_i   = dfile("electric/physics/sigmat_i", 3.1525e-4); 	// 3.1525e-4 [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmal_e =  dfile("electric/physics/sigmal_e", 2e-3); 		// 2e-3      [1/Ohm/cm]   ColliPavarinoTaccardi2005
		M_sigmat_e   = dfile("electric/physics/sigmat_e",1.3514e-3); 	// 1.3514e-3 [1/Ohm/cm]   ColliPavarinoTaccardi2005
		}
	else if(dfile("electric/physics/ion_model",1) == 3)
                {
                M_D_i = dfile("electric/physics/D_i" ,3.3e-2);       // 3.3e-2  [1/Ohm/cm]   D_i_LR * D_RM/D_LR  see dataMonodomain
                M_D_e = dfile("electric/physics/D_e" ,4.29e-2);      // 4.29e-2 [1/Ohm/cm]      D_e_LR * D_RM/D_LR  see dataMonodomain
                M_sigmal_i =  dfile("electric/physics/sigmal_i", 8.19e-2); // 8.19e-2 [1/Ohm/cm]   sigmal_i_LR * D_RM/D_LR
                M_sigmat_i   = dfile("electric/physics/sigmat_i", 8.6e-3);      // 8.6e-3  [1/Ohm/cm]   sigmat_i_LR * D_RM/D_LR
                M_sigmal_e =  dfile("electric/physics/sigmal_e", 5.46e-2); // 5.46e-2 [1/Ohm/cm]        sigmal_e_LR * D_RM/D_LR
                M_sigmat_e   = dfile("electric/physics/sigmat_e",3.69e-2); // 3.69e-2 [1/Ohm/cm]        sigmat_e_LR * D_RM/D_LR
                }
    M_heart_diff_fct 	= dfile("electric/physics/heart_diff_fct",0);
    M_verbose 		= dfile( "electric/miscellaneous/verbose", 1 );
    M_post_dir  	= dfile("electric/miscellaneous/post_dir","./");
    M_uOrder 		= dfile( "electric/space_discretization/u_order", "P1");
    M_fibers_format 	= dfile("electric/space_discretization/fibers_format",0);
    M_has_fibers 	= dfile( "electric/space_discretization/has_fibers", 0);
    M_order_bdf       	= dfile("electric/time_discretization/BDF_order",1);
    if (M_has_fibers)
    {
        std::string fibers_dir = dfile( "electric/space_discretization/fibers_dir", this->meshDir().c_str() );
        std::string fibers_file = this->meshFile();
        fibers_file.replace(fibers_file.find(".mesh"), 5, "fibers");
        M_fibers_file = fibers_dir + dfile( "electric/space_discretization/fibers_file", fibers_file.c_str() );
        std::cout<<"Fibers File: "<<M_fibers_file<<std::endl;
    }
    else
    {
        M_fibers_file="";
        std::cout<<"Fibers not included!"<<std::endl;
    }
}

// Output
//template <typename Mesh>
void DataBidomain::
showMe( std::ostream& c )
{
    c << "\n*** Values for data [fluid/physics]\n\n";
    c << "endtime   = " << getEndTime() << std::endl;
    c << "\n*** Values for data [fluid/miscellaneous]\n\n";
    c << "verbose   = " << M_verbose << std::endl;
}





}
