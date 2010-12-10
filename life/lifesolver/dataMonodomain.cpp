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
    @brief File containing a class for handling Monodomain data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#include <life/lifesolver/dataMonodomain.hpp>


namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
DataMonodomain::DataMonodomain(   boost::shared_ptr<HeartFunctors> heart ) :
    DataMesh                          ( heart-> M_dataFile, "electric/space_discretization" ),
    DataTime                          ( heart-> M_dataFile, "electric/time_discretization" ),
    M_reducedConductivityBox        ( heart-> get_reduced_sigma_box() ),
    M_reducedConductivityCylinder   ( heart-> get_reduced_sigma_cylinder() ),
    M_reducedConductivitySphere     ( heart-> get_reduced_sigma_sphere() )


{
    setup                           ( heart-> M_dataFile);
}

DataMonodomain::DataMonodomain() :
    DataMesh                        ( ),
    DataTime                        ( ),
    M_hasFibers                     ( ),
    M_verbose                       ( ),
    M_heartDiffusionFactor          ( ),
    M_diffusivity                   ( ),
    M_longitudinalConductivity      ( ),
    M_membraneCapacitance           ( ),
    M_transversalConductivity       ( ),
    M_volumeSurfaceRatio            ( ),
    M_fibersDirectory               ( ),
    M_fibersFile                    ( ),
    M_postProcessingDirectory       ( ),
    M_uOrder                        ( )

{
}

DataMonodomain::DataMonodomain( const DataMonodomain& dataMonodomain ) :
    DataMesh                        ( dataMonodomain ),
    DataTime                        ( dataMonodomain ),
    M_hasFibers                     ( dataMonodomain.M_hasFibers ),
    M_verbose                       ( dataMonodomain.M_verbose ),
    M_heartDiffusionFactor          ( dataMonodomain.M_heartDiffusionFactor ),
    M_diffusivity                   ( dataMonodomain.M_diffusivity ),
    M_longitudinalConductivity      ( dataMonodomain.M_longitudinalConductivity ),
    M_membraneCapacitance           ( dataMonodomain.M_membraneCapacitance ),
    M_transversalConductivity       ( dataMonodomain.M_transversalConductivity ),
    M_volumeSurfaceRatio            ( dataMonodomain.M_volumeSurfaceRatio ),
    M_fibersDirectory               ( dataMonodomain.M_fibersDirectory ),
    M_fibersFile                    ( dataMonodomain.M_fibersFile ),
    M_postProcessingDirectory       ( dataMonodomain.M_postProcessingDirectory ),
    M_uOrder                        ( dataMonodomain.M_uOrder )

{
}


// ===================================================
// Methods
// ===================================================
DataMonodomain&
DataMonodomain::operator=( const DataMonodomain& dataMonodomain )
{
    if ( this != &dataMonodomain )
    {
        M_hasFibers                     = dataMonodomain.M_hasFibers;
        M_heartDiffusionFactor          = dataMonodomain.M_heartDiffusionFactor;
        M_diffusivity                   = dataMonodomain.M_diffusivity;
        M_longitudinalConductivity      = dataMonodomain.M_longitudinalConductivity;
        M_membraneCapacitance           = dataMonodomain.M_membraneCapacitance;
        M_transversalConductivity       = dataMonodomain.M_transversalConductivity;
        M_volumeSurfaceRatio            = dataMonodomain.M_volumeSurfaceRatio;
        M_lambda                        = dataMonodomain.M_lambda;
        M_fibersDirectory               = dataMonodomain.M_fibersDirectory;
        M_fibersFile                    = dataMonodomain.M_fibersFile;
        M_postProcessingDirectory       = dataMonodomain.M_postProcessingDirectory;
        M_uOrder                        = dataMonodomain.M_uOrder;
        M_reducedConductivityBox        = dataMonodomain.M_reducedConductivityBox;
        M_reducedConductivityCylinder   = dataMonodomain.M_reducedConductivityCylinder;
        M_reducedConductivitySphere     = dataMonodomain.M_reducedConductivitySphere;
    }
    return *this;
}


void
DataMonodomain::setup(  const GetPot& dataFile )
{
    M_volumeSurfaceRatio           = dataFile("electric/physics/Chi", 1e3);   // [1e-3 1/cm]    ColliPavarinoTaccardi2005
    M_membraneCapacitance          = dataFile("electric/physics/Cm", 1e-3);	  // [1e-3 mF/cm2]  ColliPavarinoTaccardi2005

    if ( dataFile("electric/physics/ion_model",1) == 1)
    {
        M_diffusivity              = dataFile("electric/physics/D" , 0.0156);   // 0.0156   [1/Ohm/cm]    L^2/T*D,  L=0.099 cm, T=0.63 ms D=1,  //RogersMcCulloch1994
        M_longitudinalConductivity = dataFile("electric/physics/sigmal", 0.0328);    // 0.0328   [1/Ohm/cm]   sigmal_LR * D_RM/D_LR
        M_transversalConductivity  = dataFile("electric/physics/sigmat", 0.00699);    // 0.00699  [1/Ohm/cm]   sigmat_LR * D_RM/D_LR
    }
    else if ( dataFile("electric/physics/ion_model",1) == 2)
    {
        M_diffusivity = dataFile("electric/physics/D" , 5.7e-4);  // 5.7e-4 [1/Ohm/cm]              sigmal/3 + sigmat*2/3
        M_longitudinalConductivity = dataFile("electric/physics/sigmal", 1.2e-3);  // 1.2e-3  [1/Ohm/cm]   sigmal_i*sigmal_e/(sigmal_i+sigmal_e)    ColliPavarinoTaccardi2005
        M_transversalConductivity  = dataFile("electric/physics/sigmat", 2.56e-4); // 2.56e-4 [1/Ohm/cm]   sigmat_i*sigmat_e/(sigmat_i+sigmat_e)    ColliPavarinoTaccardi2005
    }
    M_lambda  =  dataFile("electric/physics/lambda", 0.66667); // 0.66667 [adim]       sigmal_e/sigmal_i
    M_heartDiffusionFactor         = dataFile("electric/physics/heartDiffusionFunctor",0);
    M_postProcessingDirectory      = dataFile("electric/miscellaneous/post_dir","./");
    M_uOrder                       = dataFile( "electric/space_discretization/u_order", "P1");
    M_hasFibers                    = dataFile( "electric/space_discretization/hasFibers", 0);

    if ( M_hasFibers )
    {
        std::string fibersDirectory = dataFile( "electric/space_discretization/fibers_dir", this->meshDir().c_str() );
        std::string fibersFile = this -> meshFile();
        fibersFile.replace( fibersFile.find(".mesh"), 5, "fibers" );
        M_fibersFile = fibersDirectory + dataFile( "electric/space_discretization/fibers_file", fibersFile.c_str() );
        std::cout << "Fibers File: " << M_fibersFile << std::endl;
    }
    else
    {
        M_fibersFile="";
        std::cout << "Fibers not included!" << std::endl;
    }
}


void
DataMonodomain::showMe( std::ostream& output )
{
    output << "\n*** Values for data [fluid/time_discretization]\n\n";
    output << "endtime   = " << getEndTime() << std::endl;
    output << "\n*** Values for data [fluid/miscellaneous]\n\n";
    output << "verbose   = " << M_verbose << std::endl;
}

}
