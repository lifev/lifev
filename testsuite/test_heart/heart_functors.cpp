#include "heart_functors.hpp"

using namespace LifeV;


HeartFunctors::HeartFunctors( GetPot& dataFile ):	
	_dataFile(dataFile),
    stim_source (dataFile("electric/physics/stim_source",1)),
    stim_period_1 (dataFile("electric/physics/stim_period_1",200.)),
    stim_period_2 (dataFile("electric/physics/stim_period_2",200.)),	
    stim_period_3 (dataFile("electric/physics/stim_period_3",200.)),	
    stim_period_4 (dataFile("electric/physics/stim_period_4",200.)),	
    stim_period_5 (dataFile("electric/physics/stim_period_5",200.)),	
    stim_period_6 (dataFile("electric/physics/stim_period_6",200.)),	
    stim_start_1 (dataFile("electric/physics/stim_start_1",0.)),
    stim_stop_1 (dataFile("electric/physics/stim_stop_1",0.)),
    stim_value_1 (dataFile("electric/physics/stim_value_1",0.)),
    stim_radius_1 (dataFile("electric/physics/stim_radius_1",0.)),
    stim_center_1 (3),
    stim_start_2 (dataFile("electric/physics/stim_start_2",0.)),
    stim_stop_2 (dataFile("electric/physics/stim_stop_2",0.)),
    stim_value_2 (dataFile("electric/physics/stim_value_2",0.)),
    stim_radius_2 (dataFile("electric/physics/stim_radius_2",0.)),
    stim_center_2 (3),
    stim_start_3 (dataFile("electric/physics/stim_start_3",0.)),
    stim_stop_3 (dataFile("electric/physics/stim_stop_3",0.)),
    stim_value_3 (dataFile("electric/physics/stim_value_3",0.)),
    stim_radius_3 (dataFile("electric/physics/stim_radius_3",0.)),
    stim_center_3 (3),
    stim_start_4 (dataFile("electric/physics/stim_start_4",0.)),
    stim_stop_4 (dataFile("electric/physics/stim_stop_4",0.)),
    stim_value_4 (dataFile("electric/physics/stim_value_4",0.)),
    stim_radius_4 (dataFile("electric/physics/stim_radius_4",0.)),
    stim_center_4 (3),
    stim_start_5 (dataFile("electric/physics/stim_start_5",0.)),
    stim_stop_5 (dataFile("electric/physics/stim_stop_5",0.)),
    stim_value_5 (dataFile("electric/physics/stim_value_5",0.)),
    stim_radius_5 (dataFile("electric/physics/stim_radius_5",0.)),
    stim_center_5 (3),
    stim_start_6 (dataFile("electric/physics/stim_start_6",0.)),
    stim_stop_6 (dataFile("electric/physics/stim_stop_6",0.)),
    stim_value_6 (dataFile("electric/physics/stim_value_6",0.)),
    stim_radius_6 (dataFile("electric/physics/stim_radius_6",0.)),
    stim_center_6 (3),
    x_sphere(dataFile("electric/physics/sphere_center",0.,0)),
    y_sphere(dataFile("electric/physics/sphere_center",0.,1)),
    z_sphere(dataFile("electric/physics/sphere_center",0.,2)),
    r_sphere(dataFile("electric/physics/sphere_radius",0.)),
    sigma_reduction (2),
    x_cylinder(dataFile("electric/physics/x_cylinder",0.)),
    y_cylinder(dataFile("electric/physics/y_cylinder",0.)),
    z_cylinder(dataFile("electric/physics/z_cylinder",0.)),
    a_cylinder(dataFile("electric/physics/a_cylinder",0.)),
    b_cylinder(dataFile("electric/physics/b_cylinder",0.)),
    c_cylinder(dataFile("electric/physics/c_cylinder",0.)),
    r_cylinder(dataFile("electric/physics/r_cylinder",0.)),
    xmin_cylinder(dataFile("electric/physics/xmin_cylinder",0.)),
    xmax_cylinder(dataFile("electric/physics/xmax_cylinder",0.)),
    xmin_box(dataFile("electric/physics/box_vertex_min",0.,0)),
    ymin_box(dataFile("electric/physics/box_vertex_min",0.,1)),
    zmin_box(dataFile("electric/physics/box_vertex_min",0.,2)),
    xmax_box(dataFile("electric/physics/box_vertex_max",0.,0)),
    ymax_box(dataFile("electric/physics/box_vertex_max",0.,1)),
    zmax_box(dataFile("electric/physics/box_vertex_max",0.,2)),
    //parametre Iapp IappZygote REO
    G_Time_period (dataFile("electric/physics/Time_period",700.0)),
    G_Iapp_RV_angle (dataFile("electric/physics/Iapp_RV_angle",360.)),
    G_Iapp_LV_angle (dataFile("electric/physics/Iapp_LV_angle",360.)),
    G_Iapp_stim_time_RV (dataFile("electric/physics/Iapp_stim_time_RV",6.)),
    G_Iapp_stim_time_LV (dataFile("electric/physics/Iapp_stim_time_LV",10.)),
    G_Ventricular_Fibrillation (dataFile("electric/physics/Ventricular_Fibrillation",0)),
    G_nb_fibrillation_sources (dataFile("electric/physics/nb_fibrillation_sources",20)),
    G_fibrillation_sources (dataFile("electric/physics/fibrillation_sources",0))
    {
	sigma_reduction(0) = dataFile("electric/physics/sigma_reduction",1.,0);
	sigma_reduction(1) = dataFile("electric/physics/sigma_reduction",1.,1);
    stim_center_1(0)= dataFile("electric/physics/stim_center_1",0.,0);
    stim_center_1(1)= dataFile("electric/physics/stim_center_1",0.,1);
    stim_center_1(2)= dataFile("electric/physics/stim_center_1",0.,2);
    stim_center_2(0)= dataFile("electric/physics/stim_center_2",0.,0);
    stim_center_2(1)= dataFile("electric/physics/stim_center_2",0.,1);
    stim_center_2(2)= dataFile("electric/physics/stim_center_2",0.,2);
    stim_center_3(0)= dataFile("electric/physics/stim_center_3",0.,0);
    stim_center_3(1)= dataFile("electric/physics/stim_center_3",0.,1);
    stim_center_3(2)= dataFile("electric/physics/stim_center_3",0.,2);
    stim_center_4(0)= dataFile("electric/physics/stim_center_4",0.,0);
    stim_center_4(1)= dataFile("electric/physics/stim_center_4",0.,1);
    stim_center_4(2)= dataFile("electric/physics/stim_center_4",0.,2);
    stim_center_5(0)= dataFile("electric/physics/stim_center_5",0.,0);
    stim_center_5(1)= dataFile("electric/physics/stim_center_5",0.,1);
    stim_center_5(2)= dataFile("electric/physics/stim_center_5",0.,2);
    stim_center_6(0)= dataFile("electric/physics/stim_center_6",0.,0);
    stim_center_6(1)= dataFile("electric/physics/stim_center_6",0.,1);
    stim_center_6(2)= dataFile("electric/physics/stim_center_6",0.,2);
    }
