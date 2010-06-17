//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
 *  @file
 *  @brief DataElasticStructure - File containing a data container for solid problems with elastic structure
 *
 *  @version 1.0
 *  @author M.A. Fernandez
 *  @date 01-10-2003
 *
 *  @version 1.18
 *  @author Cristiano Malossi
 *  @date 10-06-2010
 */

#ifndef DATAELASTICSTRUCTURE_H
#define DATAELASTICSTRUCTURE_H

#include <string>
#include <iostream>
#include <map>
#include <life/lifecore/GetPot.hpp>
#include <life/lifecore/life.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>

namespace LifeV {

//! DataElasticStructure - Data container for solid problems with elastic structure
/*!
 *  @author M.A. Fernandez
 */
template <typename Mesh>
class DataElasticStructure
{
public:

    //! @name Type definitions
    //@{

    typedef DataTime                                                  Time_Type;
    typedef boost::shared_ptr< Time_Type >                            Time_ptrType;

    typedef DataMesh<Mesh>                                            Mesh_Type;
    typedef boost::shared_ptr< Mesh_Type >                            Mesh_ptrType;

    typedef std::map<UInt, Real>                                      MaterialContainer_Type;
    typedef MaterialContainer_Type::const_iterator                    MaterialContainer_ConstIterator;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Empty Constructor
    DataElasticStructure();

    //! Copy constructor
    /*!
     * @param dataElasticStructure - DataElasticStructure
     */
    DataElasticStructure( const DataElasticStructure& dataElasticStructure );

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param dataElasticStructure - DataElasticStructure
     */
    DataElasticStructure& operator=( const DataElasticStructure& dataElasticStructure );

    //@}


    //! @name Methods
    //@{

    //! Read the dataFile and set all the quantities
    /*!
     * @param dataFile data file
     * @param section section of the file
     */
    void setup( const GetPot& dataFile, const std::string& section = "solid" );

    //! Display the values
    void showMe( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set methods
    //@{

    //! Set data time container
    /*!
     * @param DataTime shared_ptr to dataTime container
     */
    inline void setDataTime( const Time_ptrType DataTime );

    //! Set mesh container
    /*!
     * @param DataMesh shared_ptr to dataMesh container
     */
    inline void setDataMesh( const Mesh_ptrType DataMesh );

    //! Set density
    /*!
     * @param density solid density value
     */
    inline void setDensity( const Real& density );

    //! Set thickness
    /*!
     * @param thickness solid thickness value
     */
    inline void setThickness( const Real& thickness );

    //! Set poisson
    /*!
     * @param poisson solid poisson value
     * @param material material ID (1 by default)
     */
    inline void setPoisson( const Real& poisson, const UInt& material = 1 );

    //! Set Young modulus
    /*!
     * @param Young solid young modulus value
     * @param material material ID (1 by default)
     */
    inline void setYoung( const Real& young, const UInt& material = 1 );

    //@}


    //! @name Get methods
    //@{

    //! Get data time container
    /*!
     * @return shared_ptr to dataTime container
     */
    Time_ptrType dataTime() const;

    //! Get mesh container
    /*!
     * @return shared_ptr to dataMesh container
     */
    Mesh_ptrType dataMesh() const;

    //! Get solid density
    /*!
     * @return Solid density
     */
    inline const Real& rho() const;

    //! Get solid thickness
    /*!
     * @return Solid thickness
     */
    inline const Real& thickness() const;

    //! Get solid poisson coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid poisson coefficient
     */
    inline const Real& poisson( const UInt& material = 1 ) const;

    //! Get solid young modulus
    /*!
     * @param material material ID (1 by default)
     * @return Solid young modulus
     */
    inline const Real& young( const UInt& material = 1 ) const;

    //! Get solid first lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid first Lame coefficient
     */
    inline Real lambda( const UInt& material = 1 ) const;

    //! Get solid second Lame coefficient
    /*!
     * @param material material ID (1 by default)
     * @return Solid second Lame coefficient
     */
    inline Real mu( const UInt& material = 1 ) const;

    //! Get FE order
    /*!
     * @return FE order
     */
    inline const std::string& order()     const;

    //! Get solid amplification factor
    /*!
     * @return Solid amplification factor
     */
    inline const Real& factor()    const;

    //! Get verbose level
    /*!
     * @return verbose level
     */
    inline const UInt& verbose()   const;

    //@}

private:

    //! Data containers for time and mesh
    Time_ptrType           M_time;
    Mesh_ptrType           M_mesh;

    //! Physics
    Real                   M_density;
    Real                   M_thickness;

    MaterialContainer_Type M_poisson;
    MaterialContainer_Type M_young;

    //! Space discretization
    std::string            M_order;

    //! Miscellaneous
    Real                   M_factor;  // amplification factor for deformed mesh
    UInt                   M_verbose; // temporal output verbose
};





// ===================================================
// Constructors
// ===================================================
template <typename Mesh>
DataElasticStructure<Mesh>::DataElasticStructure() :
        M_time                             ( ),
        M_mesh                             ( ),
        M_density                          ( ),
        M_thickness                        ( ),
        M_poisson                          ( ),
        M_young                            ( ),
        M_order                            ( ),
        M_factor                           ( ),
        M_verbose                          ( )
{
}

template <typename Mesh>
DataElasticStructure<Mesh>::DataElasticStructure( const DataElasticStructure& dataElasticStructure ):
    M_time                             ( dataElasticStructure.M_time ),
    M_mesh                             ( dataElasticStructure.M_mesh ),
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
template <typename Mesh>
DataElasticStructure<Mesh>&
DataElasticStructure<Mesh>::operator=( const DataElasticStructure& dataElasticStructure )
{
    if ( this != &dataElasticStructure )
    {
        M_time                             = dataElasticStructure.M_time;
        M_mesh                             = dataElasticStructure.M_mesh;
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
template <typename Mesh>
void
DataElasticStructure<Mesh>::setup( const GetPot& dataFile, const std::string& section )
{
    // If data time has not been set
    if ( !M_time.get() )
        M_time.reset( new Time_Type( dataFile, section + "/time_discretization" ) );

    // If data mesh has not been set
    if ( !M_mesh.get() )
        M_mesh.reset( new Mesh_Type( dataFile, section + "/space_discretization" ) );

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

template <typename Mesh>
void
DataElasticStructure<Mesh>::showMe( std::ostream& output ) const
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
    M_mesh->showMe( output );

    output << "\n*** Values for data [solid/time_discretization]\n\n";
    M_time->showMe( output );
}

// ===================================================
// Set Method
// ===================================================
template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setDataTime( const Time_ptrType DataTime )
{
    M_time = DataTime;
}

template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setDataMesh( const Mesh_ptrType DataMesh )
{
    M_mesh = DataMesh;
}

template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setDensity( const Real& density )
{
    M_density = density;
}

template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setThickness( const Real& thickness )
{
    M_thickness = thickness;
}

template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setPoisson( const Real& poisson, const UInt& material )
{
    M_poisson[material] = poisson;
}

template <typename Mesh>
inline void
DataElasticStructure<Mesh>::setYoung( const Real& young, const UInt& material )
{
    M_young[material] = young;
}

// ===================================================
// Get Method
// ===================================================
template <typename Mesh>
typename DataElasticStructure<Mesh>::Time_ptrType
DataElasticStructure<Mesh>::dataTime() const
{
    return M_time;
}

template <typename Mesh>
typename DataElasticStructure<Mesh>::Mesh_ptrType
DataElasticStructure<Mesh>::dataMesh() const
{
    return M_mesh;
}

template <typename Mesh>
inline const Real&
DataElasticStructure<Mesh>::rho() const
{
    return M_density;
}

template <typename Mesh>
inline const Real&
DataElasticStructure<Mesh>::thickness() const
{
    return M_thickness;
}

template <typename Mesh>
inline const Real&
DataElasticStructure<Mesh>::poisson( const UInt& material ) const
{
    return M_poisson.find( material )->second;
}

template <typename Mesh>
inline const Real&
DataElasticStructure<Mesh>::young( const UInt& material ) const
{
    return M_young.find( material )->second;
}

template <typename Mesh>
inline Real
DataElasticStructure<Mesh>::lambda( const UInt& material ) const
{
    return M_young.find( material )->second * M_poisson.find( material )->second /
           ( ( 1.0 + M_poisson.find( material )->second ) * ( 1.0 - 2.0 * M_poisson.find( material )->second ) );
}

template <typename Mesh>
inline Real
DataElasticStructure<Mesh>::mu( const UInt& material ) const
{
    return M_young.find( material )->second/( 2.0 * ( 1.0 + M_poisson.find( material )->second ) );
}

template <typename Mesh>
inline const std::string&
DataElasticStructure<Mesh>::order() const
{
    return M_order;
}

template <typename Mesh>
inline const Real&
DataElasticStructure<Mesh>::factor() const
{
    return M_factor;
}

template <typename Mesh>
inline const UInt&
DataElasticStructure<Mesh>::verbose() const
{
    return M_verbose;
}

} // end namespace LifeV

#endif // end DATAELASTICSTRUCTURE_H
