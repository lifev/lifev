/*
 This file is part of the LifeV library
 Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
  \file dataElasticStructure.h
  \author M.A. Fernandez
  \date 10/2003
  \version 1.0
  \brief
*/

#ifndef _DATAELASTICSTRUCTURE_H_
#define _DATAELASTICSTRUCTURE_H_
#include <string>
#include <iostream>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>
#include <map>
#include <list>

namespace LifeV
{
/*!
  \class DataElasticStructure
*/
template <typename Mesh>
class DataElasticStructure:
            public DataMesh<Mesh>,
            public DataTime
{
public:

    //! Constructor
    DataElasticStructure( const GetPot& dfile );
    DataElasticStructure( const DataElasticStructure &dataElasticStructure );

    //! Ouptut
    void showMe( std::ostream& c = std::cout ) const;

    //! End time
    Real endtime() const;

    //! getters

    const Real rho()     const {return _rho;}
//     const Real young()   const {return _E;}
//     const Real poisson() const {return _nu;}

//     const Real lambda()  const {return _lambda;}
//     const Real mu()      const {return _mu;}
    const Real factor()  const {return _factor;}

    const UInt verbose() const {return _verbose;}

    const double poisson(int mat) const
        {return M_poisson.find(mat)->second;}
    const double young(int mat) const
        {return M_young.find(mat)->second;}

    const double lambda(int mat) const
        {return M_lambda.find(mat)->second;}
    const double mu(int mat) const
        {return M_mu.find(mat)->second;}

    const double lambda() const {return _lambda;}
    const double mu() const {return _mu;}

    std::string order() const {return M_order;}
private:
    //! Physics
    Real _rho; // densisty
    Real _E;  // Young modulus
    Real _nu; // Poisson coeficient
    Real _lambda, _mu; // Lame coefficients
    Real _endtime; // end time

    std::map<int, double>  M_poisson;
    std::map<int, double>  M_young;

    std::map<int, double>  M_lambda;
    std::map<int, double>  M_mu;

    //! Miscellaneous
    Real _factor; // amplification factor for deformed mesh
    UInt _verbose; // temporal output verbose
    std::string M_order;


};



//
// IMPLEMENTATION
//


// Constructor
template <typename Mesh>
DataElasticStructure<Mesh>::
DataElasticStructure( const GetPot& dfile ) :
    DataMesh<Mesh>( dfile, "solid/discretization" ),
    DataTime( dfile, "solid/discretization" )
{
    // physics
    _rho     = dfile( "solid/physics/density", 1. );
    _E       = dfile( "solid/physics/young" , 1. );
    _nu      = dfile( "solid/physics/poisson" , 0.25 );
    _endtime = dfile( "solid/physics/endtime", 1. );

    // miscellaneous
    _factor  = dfile( "solid/miscellaneous/factor", 1.0 );
    std::cout << "factor " << _factor << std::endl;
    _verbose = dfile( "solid/miscellaneous/verbose", 1 );

    M_order  = dfile( "solid/discretization/order", "P1");

    // Lame coefficients
    _lambda  = _E * _nu / ( ( 1.0 + _nu ) * ( 1.0 - 2.0 * _nu ) );
    _mu      = _E / ( 2.0 * ( 1.0 + _nu ) );

    std::string flagList;
    flagList = dfile("solid/physics/material_flag", "0");
    std::list<int> fList;
    parseList(flagList, fList);

    std::string youngList;
    youngList =  dfile("solid/physics/young", "0.");
    std::list<double> yList;
    parseList(youngList, yList);

    std::string poissonList;
    poissonList = dfile("solid/physics/poisson", "0.");
    std::list<double> pList;
    parseList(poissonList, pList);

//    if ((fList.size() != yList.size()) || (flist.size() != pList.size()))
    ASSERT((fList.size() == yList.size()) && (fList.size() == pList.size()),"problem with the young modulus and poisson coef definiton : inconsistant sizes");

    std::list<int>::iterator    fit;
    std::list<double>::iterator yit = yList.begin();
    std::list<double>::iterator pit = pList.begin();

    std::cout << "flag       young       poisson" << std::endl;
    for (fit = fList.begin(); fit != fList.end(); ++fit, ++yit, ++pit)
    {
        double young   = *yit;
        double poisson = *pit;

        M_young.insert(make_pair(*fit, young));
        M_poisson.insert(make_pair(*fit, poisson));

        std::cout << *fit << " " << young << " " << poisson << std::endl;

        double lambda = young*poisson/( ( 1.0 + poisson ) * ( 1.0 - 2.0 * poisson ) );
        M_lambda.insert(make_pair(*fit, lambda));

        double mu = young/( 2.0 * ( 1.0 + poisson ) );
        M_mu.insert(make_pair(*fit, mu));

        std::cout << *fit << " " << lambda << " " << mu << std::endl;
    }

}


template <typename Mesh>
DataElasticStructure<Mesh>::
DataElasticStructure( const DataElasticStructure& dataElasticStructure ):
    DataMesh<Mesh>               ( dataElasticStructure ),
    DataTime                     ( dataElasticStructure ),
    _rho                         ( dataElasticStructure._rho ),
    _E                           ( dataElasticStructure._E ),
    _nu                          ( dataElasticStructure._nu ),
    _lambda                      ( dataElasticStructure._lambda ),
    _mu                          ( dataElasticStructure._mu ),
    _endtime                     ( dataElasticStructure._endtime ),
    M_young                      ( dataElasticStructure.M_young ),
    M_poisson                    ( dataElasticStructure.M_poisson ),
    M_lambda                     ( dataElasticStructure.M_lambda ),
    M_mu                         ( dataElasticStructure.M_mu ),
    _factor                      ( dataElasticStructure._factor ),
    _verbose                     ( dataElasticStructure._verbose ),
    M_order                      ( dataElasticStructure.M_order )
{
}




// Output
template <typename Mesh>
void DataElasticStructure<Mesh>::
showMe( std::ostream& c ) const
{
    // physics
    c << "\n*** Values for data [solid/physics]\n\n";
    c << "density                          = " << _rho << std::endl;
    c << "young                            = " << _E << std::endl;
    c << "poisson                          = " << _nu << std::endl;
    c << "lame constants (lambda, mu)      = " << _lambda << " " << _mu << std::endl;
    c << "endtime                          = " << _endtime << std::endl;

    c << "\n*** Values for data [solid/miscellaneous]\n\n";
    c << "deformation factor               = " << _factor << std::endl;
    c << "verbose                          = " << _verbose << std::endl;

    c << "\n*** Values for data [solid/discretization]\n\n";
    DataMesh<Mesh>::showMe();
    DataTime::showMe( c );
}


// The end time
template <typename Mesh>
Real DataElasticStructure<Mesh>::
endtime() const
{
    return _endtime;
}
}

#endif
