/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file dataElasticStructure.cpp
*/

#include <string>
#include <iostream>
//#include <>
#include <map>
#include <list>
#include <life/lifecore/util_string.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifesolver/dataElasticStructure.hpp>

namespace LifeV
{


//
// IMPLEMENTATION
//

// ===================================================
// Constructors
// ===================================================


// DataElasticStructure::
// DataElasticStructure( const GetPot& dfile ) :
//         //DataMesh<Mesh>( dfile, "solid/space_discretization" ),
//         DataTime( dfile, "solid/time_discretization" )
// {
//     // physics
//     _rho     = dfile( "solid/physics/density", 1. );
// //     _E       = dfile( "solid/physics/young" , 1. );
// //     _nu      = dfile( "solid/physics/poisson" , 0.25 );

//     // miscellaneous
//     _factor  = dfile( "solid/miscellaneous/factor", 1.0 );
// //    std::cout << "factor " << _factor << std::endl;
//     _verbose = dfile( "solid/miscellaneous/verbose", 1 );

//     M_order  = dfile( "solid/space_discretization/order", "P1");

//     M_thickness = dfile("solid/physics/thickness", 0.1);

//     // Lame coefficients
// //     _lambda  = _E * _nu / ( ( 1.0 + _nu ) * ( 1.0 - 2.0 * _nu ) );
// //     _mu      = _E / ( 2.0 * ( 1.0 + _nu ) );

//     std::string flagList;
//     flagList = dfile("solid/physics/material_flag", "0");
//     std::list<int> fList;
//     parseList(flagList, fList);

//     std::string youngList;
//     youngList =  dfile("solid/physics/young", "0.");
//     std::list<double> yList;
//     parseList(youngList, yList);

//     std::string poissonList;
//     poissonList = dfile("solid/physics/poisson", "0.");
//     std::list<double> pList;
//     parseList(poissonList, pList);

// //    if ((fList.size() != yList.size()) || (flist.size() != pList.size()))
//     ASSERT((fList.size() == yList.size()) && (fList.size() == pList.size()),"problem with the young modulus and poisson coef definiton : inconsistant sizes");

//     std::list<int>::iterator    fit;
//     std::list<double>::iterator yit = yList.begin();
//     std::list<double>::iterator pit = pList.begin();



//     for (fit = fList.begin(); fit != fList.end(); ++fit, ++yit, ++pit)
//     {
//         double young   = *yit;
//         double poisson = *pit;

//         M_young.insert(std::make_pair(*fit, young));
//         M_poisson.insert(std::make_pair(*fit, poisson));

//         double lambda = young*poisson/( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );
//         M_lambda.insert(std::make_pair(*fit, lambda));

//         double mu = young/( 2.0 * ( 1.0 + poisson ) );
//         M_mu.insert(std::make_pair(*fit, mu));

//         if (fList.size() == 1)
//             {
//                 _E      = young;
//                 _nu     = poisson;

//                 _lambda = lambda;
//                 _mu     = mu;
//             }
//     }

// }




DataElasticStructure::DataElasticStructure() :
        M_time                             ( ),
        M_density                          ( ),
        M_thickness                        ( ),
        M_poisson                          ( ),
        M_young                            ( ),
        M_order                            ( ),
        M_factor                           ( ),
        M_verbose                          ( )
{
}


DataElasticStructure::DataElasticStructure( const DataElasticStructure& dataElasticStructure ):
    DataTime                           ( dataElasticStructure ),
    M_time                             ( dataElasticStructure.M_time ),
    M_density                          ( dataElasticStructure.M_density ),
    M_thickness                        ( dataElasticStructure.M_thickness ),
    M_poisson                          ( dataElasticStructure.M_poisson ),
    M_young                            ( dataElasticStructure.M_young ),
    M_order                            ( dataElasticStructure.M_order ),
    M_factor                           ( dataElasticStructure.M_factor ),
    M_verbose                          ( dataElasticStructure.M_verbose )
{
}



// ===================================================
// Operators
// ===================================================



DataElasticStructure&
DataElasticStructure::operator=( const DataElasticStructure& dataElasticStructure )
{
    if ( this != &dataElasticStructure )
    {
        M_time                             = dataElasticStructure.M_time;
        M_density                          = dataElasticStructure.M_density;
        M_thickness                        = dataElasticStructure.M_thickness;
        M_poisson                          = dataElasticStructure.M_poisson;
        M_young                            = dataElasticStructure.M_young;
        M_order                            = dataElasticStructure.M_order;
        M_factor                           = dataElasticStructure.M_factor;
        M_verbose                          = dataElasticStructure.M_verbose;
    }

    return *this;
}




// ===================================================
// Methods
// ===================================================

void
DataElasticStructure::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // physics
    M_density   = dataFile( ( section + "/physics/density" ).data(), 1. );
    M_thickness = dataFile( ( section + "/physics/thickness" ).data(), 0.1 );

    UInt materialsNumber = dataFile.vector_variable_size( ( section + "/physics/material_flag" ).data() );
    if ( materialsNumber == 0 )
    {
        M_young[1]   = dataFile( ( section + "/physics/young" ).data(), 0. );
        M_poisson[1] = dataFile( ( section + "/physics/poisson" ).data(), 0. );
    }
    else
    {
        ASSERT( materialsNumber == dataFile.vector_variable_size( ( section + "/physics/young" ).data()),   "!!! ERROR: Inconsistent size for Young Modulus !!!");
        ASSERT( materialsNumber == dataFile.vector_variable_size( ( section + "/physics/poisson" ).data() ), "!!! ERROR: Inconsistent size for Poisson Coeff. !!!");

        UInt material(0);
        for ( UInt i(0) ; i < materialsNumber ; ++i )
        {
            material            = dataFile( ( section + "/physics/material_flag" ).data(), 0., i );
            M_young[material]   = dataFile( ( section + "/physics/young" ).data(), 0. );
            M_poisson[material] = dataFile( ( section + "/physics/poisson" ).data(), 0. );
        }
    }

    // space_discretization
    M_order     = dataFile( "solid/space_discretization/order", "P1" );

    // miscellaneous
    M_factor  = dataFile( "solid/miscellaneous/factor", 1.0 );
    M_verbose = dataFile( "solid/miscellaneous/verbose", 1 );
}


void
DataElasticStructure::showMe( std::ostream& output ) const
{
    // physics
    output << "\n*** Values for data [solid/physics]\n\n";
    output << "density                          = " << M_density << std::endl;
    output << "thickness                        = " << M_thickness << std::endl;
    for ( MaterialContainer_ConstIterator i = M_young.begin() ; i != M_young.end() ; ++i )
        output << "young[" << i->first << "]                         = " << i->second << std::endl;
    for ( MaterialContainer_ConstIterator i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
        output << "poisson[" << i->first << "]                       = " << i->second << std::endl;

    for ( MaterialContainer_ConstIterator i = M_poisson.begin() ; i != M_poisson.end() ; ++i )
    {
        output << "Lame - lambda[" << i->first << "]                 = " << lambda( i->first ) << std::endl;
        output << "Lame - mu[" << i->first << "]                     = " << mu( i->first ) << std::endl;
    }

    output << "\n*** Values for data [solid/miscellaneous]\n\n";
    output << "deformation factor               = " << M_factor << std::endl;
    output << "verbose                          = " << M_verbose << std::endl;

    output << "\n*** Values for data [solid/space_discretization]\n\n";
    output << "FE order                         = " << M_order << std::endl;

    output << "\n*** Values for data [solid/time_discretization]\n\n";
    M_time->showMe( output );
}

// ===================================================
// Set Method
// ===================================================

 void
DataElasticStructure::setDataTime( const Time_ptrType DataTime )
{
    M_time = DataTime;
}



 void
DataElasticStructure::setDensity( const Real& density )
{
    M_density = density;
}


 void
DataElasticStructure::setThickness( const Real& thickness )
{
    M_thickness = thickness;
}


 void
DataElasticStructure::setPoisson( const Real& poisson, const UInt& material )
{
    M_poisson[material] = poisson;
}


 void
DataElasticStructure::setYoung( const Real& young, const UInt& material )
{
    M_young[material] = young;
}

// ===================================================
// Get Method
// ===================================================

DataElasticStructure::Time_ptrType
DataElasticStructure::dataTime() const
{
    return M_time;
}


 const Real&
DataElasticStructure::rho() const
{
    return M_density;
}


 const Real&
DataElasticStructure::thickness() const
{
    return M_thickness;
}


 const Real&
DataElasticStructure::poisson( const UInt& material ) const
{
    return M_poisson.find( material )->second;
}


 const Real&
DataElasticStructure::young( const UInt& material ) const
{
    return M_young.find( material )->second;
}


 const Real
DataElasticStructure::lambda( const UInt& material ) const
{
    return M_young.find( material )->second * M_poisson.find( material )->second /
           ( ( 1.0 + M_poisson.find( material )->second ) * ( 1.0 - 2.0 * M_poisson.find( material )->second ) );
}


 const Real
DataElasticStructure::mu( const UInt& material ) const
{
    return M_young.find( material )->second/( 2.0 * ( 1.0 + M_poisson.find( material )->second ) );
}


 const std::string&
DataElasticStructure::order() const
{
    return M_order;
}


 const Real&
DataElasticStructure::factor() const
{
    return M_factor;
}


 const UInt&
DataElasticStructure::verbose() const
{
    return M_verbose;
}



}
