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

#include <life/lifesolver/dataIonic.hpp>


namespace LifeV
{


// ===================================================
// Constructors & Destructor
// ===================================================
//! Constructors
DataIonic::DataIonic( const GetPot& dataFile ) :
        DataMesh( dataFile, "electric/space_discretization" ),
        DataTime( dataFile, "electric/time_discretization" )
{
    setup(dataFile);
}

DataIonic::DataIonic() :
        DataMesh                      ( ),
        DataTime                      ( ),
        M_verbose                     ( ),
        M_hasHeterogeneousTauClose    ( ),
        M_a                           ( ),
        M_b                           ( ),
        M_c1                          ( ),
        M_c2                          ( ),
        M_d                           ( ),
        M_timeUnit                    ( ),
        M_potentialAmplitude          ( ),
        M_restPotential               ( ),
        M_initialRepolarization       ( ),
        // Mitchell & Schaeffer
        M_tau_in                      ( ),
        M_tau_out                     ( ),
        M_tau_open                    ( ),
        M_tau_close                   ( ),
        M_criticalPotential           ( ),
        M_potentialMinimum            ( ),
        M_potentialMaximum            ( ),
        M_reactionAmplitude           ( ),
        M_initialTime                 ( ),
        M_tend                        ( ),
        M_BDForder                    ( )

{
}

DataIonic::DataIonic( const DataIonic& dataIonic ) :
        DataMesh		              ( dataIonic ),
        DataTime                      ( dataIonic ),
        M_hasHeterogeneousTauClose    ( dataIonic.M_hasHeterogeneousTauClose ),
        M_verbose                     ( dataIonic.M_verbose ),
        M_a                           ( dataIonic.M_a ),
        M_b                           ( dataIonic.M_b ),
        M_c1                          ( dataIonic.M_c1 ),
        M_c2                          ( dataIonic.M_c2 ),
        M_d                           ( dataIonic.M_d ),
        M_timeUnit                    ( dataIonic.M_timeUnit ),
        M_potentialAmplitude          ( dataIonic.M_potentialAmplitude ),
        M_restPotential               ( dataIonic.M_restPotential ),
        M_initialRepolarization       ( dataIonic.M_initialRepolarization ),
        // Mitchell & Schaeffer
        M_tau_in                      ( dataIonic.M_tau_in ),
        M_tau_out                     ( dataIonic.M_tau_out ),
        M_tau_open                    ( dataIonic.M_tau_open ),
        M_tau_close                   ( dataIonic.M_tau_close ),
        M_criticalPotential           ( dataIonic.M_criticalPotential ),
        M_potentialMinimum            ( dataIonic.M_potentialMinimum ),
        M_potentialMaximum            ( dataIonic.M_potentialMaximum ),
        M_reactionAmplitude           ( dataIonic.M_reactionAmplitude ),
        M_initialTime                 ( dataIonic.M_initialTime ),
        M_tend                        ( dataIonic.M_tend ),
        M_BDForder                    ( dataIonic.M_BDForder )

{
}


// ===================================================
// Methods
// ===================================================
DataIonic&
DataIonic::operator=( const DataIonic& dataIonic )
{
    if( this != &dataIonic )
    {
        M_hasHeterogeneousTauClose    = dataIonic.M_hasHeterogeneousTauClose;
        M_a                           = dataIonic.M_a;
        M_b                           = dataIonic.M_b;
        M_c1                          = dataIonic.M_c1;
        M_c2                          = dataIonic.M_c2;
        M_d                           = dataIonic.M_d;
        M_timeUnit                    = dataIonic.M_timeUnit;
        M_potentialAmplitude          = dataIonic.M_potentialAmplitude;
        M_restPotential               = dataIonic.M_restPotential;
        M_initialRepolarization       = dataIonic.M_initialRepolarization;
        // Mitchell & Schaeffer
        M_tau_in                      = dataIonic.M_tau_in;
        M_tau_out                     = dataIonic.M_tau_out;
        M_tau_open                    = dataIonic.M_tau_open;
        M_tau_close                   = dataIonic.M_tau_close;
        M_criticalPotential           = dataIonic.M_criticalPotential;
        M_potentialMinimum            = dataIonic.M_potentialMinimum;
        M_potentialMaximum            = dataIonic.M_potentialMaximum;
        M_reactionAmplitude           = dataIonic.M_reactionAmplitude;
        M_initialTime                 = dataIonic.M_initialTime;
        M_tend                        = dataIonic.M_tend;
        M_BDForder                    = dataIonic.M_BDForder;

    }
    return *this;
}

void
DataIonic::setup(  const GetPot& dataFile )
{
    M_hasHeterogeneousTauClose = dataFile( "electric/physics/hasHeteroTauClose",1 );
    M_a   		               = dataFile( "electric/physics/a",0.13 );   // 0.13  adim  //RogersMcCulloch1994
    M_b   		               = dataFile( "electric/physics/b",0.013 );  // 0.013 adim //RogersMcCulloch1994
    M_c1   		               = dataFile( "electric/physics/c1",0.26 );  // 0.26  adim //RogersMcCulloch1994
    M_c2   		               = dataFile( "electric/physics/c2",0.1 );   //0.1    adim //RogersMcCulloch1994
    M_d   		               = dataFile( "electric/physics/d",1 );      //1      adim //RogersMcCulloch1994
    M_timeUnit   		       = dataFile( "electric/physics/T",0.63 );    //0.63ms    //RogersMcCulloch1994
    M_potentialAmplitude   	   = dataFile( "electric/physics/A",110 );    //130mV    //RogersMcCulloch1994
    M_restPotential            = dataFile( "electric/physics/u0",-84.0 );	  //-84mV    //RogersMcCulloch1994
    M_initialRepolarization    = dataFile( "electric/physics/winit", 0 );
    // Mitchell & Schaeffer
    M_tau_in                   = dataFile( "electric/physics/tau_in",0.8 );
    M_potentialMinimum         = dataFile( "electric/physics/v_min",-80.0 );
    M_potentialMaximum         = dataFile( "electric/physics/v_max", 20.0 );
    M_reactionAmplitude        = dataFile( "electric/physics/reac_amp", 0.2 );
    M_tau_out                  = dataFile( "electric/physics/tau_out",18.0 );
    M_tau_open                 = dataFile( "electric/physics/tau_open",100.0 );
    M_tau_close                = dataFile( "electric/physics/tau_close",100.0 );
    M_criticalPotential        = dataFile( "electric/physics/vcrit",-67.0 );
    M_initialTime              = dataFile( "electric/physics/init_time",0.0 );
    M_tend                     = dataFile( "electric/physics/end_time",1000.0 );
    M_BDForder                 = dataFile( "electric/time_discretization/BDF_order",1 );

}

// Output
void DataIonic::showMe( std::ostream& output )
{

}

}

