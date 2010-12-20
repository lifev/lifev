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
    @brief File containing a class for handling Ionic model data with GetPot

    @date 11−2007
    @author Lucia Mirabella <lucia.mirabella@gmail.com>
    @author Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, J.Castelneau (INRIA)
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#include <life/lifesolver/HeartIonicData.hpp>


namespace LifeV
{


// ===================================================
// Constructors & Destructor
// ===================================================
//! Constructors
HeartIonicData::HeartIonicData( const GetPot& dataFile ) :
        DataMesh( dataFile, "electric/space_discretization" ),
        DataTime( dataFile, "electric/time_discretization" )
{
    setup(dataFile);
}

HeartIonicData::HeartIonicData() :
        DataMesh                        ( ),
        DataTime                        ( ),
        M_verbose                       ( ),
        M_RMCParameterA                 ( ),
        M_RMCParameterB                 ( ),
        M_RMCParameterC1                ( ),
        M_RMCParameterC2                ( ),
        M_RMCParameterD                 ( ),
        M_RMCTimeUnit                   ( ),
        M_RMCPotentialAmplitude         ( ),
        M_RMCRestPotential              ( ),
        M_RMCInitialRepolarization      ( ),
        // Mitchell & Schaeffer
        M_MSTauIn                       ( ),
        M_MSTauOut                       ( ),
        M_MSTauOpen                     ( ),
        M_MSTauClose                    ( ),
        M_MSCriticalPotential           ( ),
        M_MSPotentialMinimum            ( ),
        M_MSPotentialMaximum            ( ),
        M_MSReactionAmplitude           ( ),
        M_MSInitialTime                 ( ),
        M_MSTend                        ( ),
        M_MSBDForder                    ( ),
        M_MSHasHeterogeneousTauClose    ( )
{
}

HeartIonicData::HeartIonicData( const HeartIonicData& dataIonic ) :
        DataMesh		                ( dataIonic ),
        DataTime                        ( dataIonic ),
        M_verbose                       ( dataIonic.M_verbose ),
        M_RMCParameterA                 ( dataIonic.M_RMCParameterA ),
        M_RMCParameterB                 ( dataIonic.M_RMCParameterB ),
        M_RMCParameterC1                ( dataIonic.M_RMCParameterC1 ),
        M_RMCParameterC2                ( dataIonic.M_RMCParameterC2 ),
        M_RMCParameterD                 ( dataIonic.M_RMCParameterD ),
        M_RMCTimeUnit                   ( dataIonic.M_RMCTimeUnit ),
        M_RMCPotentialAmplitude         ( dataIonic.M_RMCPotentialAmplitude ),
        M_RMCRestPotential              ( dataIonic.M_RMCRestPotential ),
        M_RMCInitialRepolarization      ( dataIonic.M_RMCInitialRepolarization ),
        // Mitchell & Schaeffer
        M_MSTauIn                       ( dataIonic.M_MSTauIn ),
        M_MSTauOut                       ( dataIonic.M_MSTauOut ),
        M_MSTauOpen                     ( dataIonic.M_MSTauOpen ),
        M_MSTauClose                    ( dataIonic.M_MSTauClose ),
        M_MSCriticalPotential           ( dataIonic.M_MSCriticalPotential ),
        M_MSPotentialMinimum            ( dataIonic.M_MSPotentialMinimum ),
        M_MSPotentialMaximum            ( dataIonic.M_MSPotentialMaximum ),
        M_MSReactionAmplitude           ( dataIonic.M_MSReactionAmplitude ),
        M_MSInitialTime                 ( dataIonic.M_MSInitialTime ),
        M_MSTend                        ( dataIonic.M_MSTend ),
        M_MSBDForder                    ( dataIonic.M_MSBDForder ),
        M_MSHasHeterogeneousTauClose    ( dataIonic.M_MSHasHeterogeneousTauClose )
{
}


// ===================================================
// Methods
// ===================================================
HeartIonicData&
HeartIonicData::operator=( const HeartIonicData& dataIonic )
{
    if( this != &dataIonic )
    {
        M_MSHasHeterogeneousTauClose    = dataIonic.M_MSHasHeterogeneousTauClose;
        M_RMCParameterA                           = dataIonic.M_RMCParameterA;
        M_RMCParameterB                           = dataIonic.M_RMCParameterB;
        M_RMCParameterC1                          = dataIonic.M_RMCParameterC1;
        M_RMCParameterC2                          = dataIonic.M_RMCParameterC2;
        M_RMCParameterD                           = dataIonic.M_RMCParameterD;
        M_RMCTimeUnit                    = dataIonic.M_RMCTimeUnit;
        M_RMCPotentialAmplitude          = dataIonic.M_RMCPotentialAmplitude;
        M_RMCRestPotential               = dataIonic.M_RMCRestPotential;
        M_RMCInitialRepolarization       = dataIonic.M_RMCInitialRepolarization;
        // Mitchell & Schaeffer
        M_MSTauIn                      = dataIonic.M_MSTauIn;
        M_MSTauOut                     = dataIonic.M_MSTauOut;
        M_MSTauOpen                    = dataIonic.M_MSTauOpen;
        M_MSTauClose                   = dataIonic.M_MSTauClose;
        M_MSCriticalPotential           = dataIonic.M_MSCriticalPotential;
        M_MSPotentialMinimum            = dataIonic.M_MSPotentialMinimum;
        M_MSPotentialMaximum            = dataIonic.M_MSPotentialMaximum;
        M_MSReactionAmplitude           = dataIonic.M_MSReactionAmplitude;
        M_MSInitialTime                 = dataIonic.M_MSInitialTime;
        M_MSTend                        = dataIonic.M_MSTend;
        M_MSBDForder                    = dataIonic.M_MSBDForder;

    }
    return *this;
}

void
HeartIonicData::setup(  const GetPot& dataFile )
{
    M_MSHasHeterogeneousTauClose = dataFile( "electric/physics/hasHeteroTauClose",1 );
    M_RMCParameterA   		               = dataFile( "electric/physics/a",0.13 );   // 0.13  adim  //RogersMcCulloch1994
    M_RMCParameterB   		               = dataFile( "electric/physics/b",0.013 );  // 0.013 adim //RogersMcCulloch1994
    M_RMCParameterC1   		               = dataFile( "electric/physics/c1",0.26 );  // 0.26  adim //RogersMcCulloch1994
    M_RMCParameterC2   		               = dataFile( "electric/physics/c2",0.1 );   //0.1    adim //RogersMcCulloch1994
    M_RMCParameterD   		               = dataFile( "electric/physics/d",1 );      //1      adim //RogersMcCulloch1994
    M_RMCTimeUnit   		       = dataFile( "electric/physics/T",0.63 );    //0.63ms    //RogersMcCulloch1994
    M_RMCPotentialAmplitude   	   = dataFile( "electric/physics/A",110 );    //130mV    //RogersMcCulloch1994
    M_RMCRestPotential            = dataFile( "electric/physics/u0",-84.0 );	  //-84mV    //RogersMcCulloch1994
    M_RMCInitialRepolarization    = dataFile( "electric/physics/winit", 0 );
    // Mitchell & Schaeffer
    M_MSTauIn                   = dataFile( "electric/physics/tau_in",0.8 );
    M_MSPotentialMinimum         = dataFile( "electric/physics/v_min",-80.0 );
    M_MSPotentialMaximum         = dataFile( "electric/physics/v_max", 20.0 );
    M_MSReactionAmplitude        = dataFile( "electric/physics/reac_amp", 0.2 );
    M_MSTauOut                  = dataFile( "electric/physics/tau_out",18.0 );
    M_MSTauOpen                 = dataFile( "electric/physics/tau_open",100.0 );
    M_MSTauClose                = dataFile( "electric/physics/tau_close",100.0 );
    M_MSCriticalPotential        = dataFile( "electric/physics/vcrit",-67.0 );
    M_MSInitialTime              = dataFile( "electric/physics/init_time",0.0 );
    M_MSTend                     = dataFile( "electric/physics/end_time",1000.0 );
    M_MSBDForder                 = dataFile( "electric/time_discretization/BDF_order",1 );

}

// Output
void HeartIonicData::showMe( std::ostream& output )
{
 output << " The output is still to be coded! \n" << std::endl;
}

}

