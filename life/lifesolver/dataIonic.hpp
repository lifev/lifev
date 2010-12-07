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
    @author Lucia Mirabella <lucia.mirabella@gmail.com>, Mauro Perego <perego.mauro@gmail.com>

    @contributor Simone Rossi <simone.rossi@epfl.ch>, Ricardo Ruiz-Baier <ricardo.ruiz@epfl.ch>, J.Castelneau (INRIA)
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef _DATAIONIC_H_
#define _DATAIONIC_H_

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/dataTime.hpp>

namespace LifeV
{
/*!
  \class DataIonic

  Base class which holds usual data for the ionic model solvers

*/
class DataIonic:
        public DataMesh,
        public DataTime
{
public:

    //! @name Constructor & Destructor
    //@{

    //! Constructors
    DataIonic();

    DataIonic( const GetPot& dataFile );

    DataIonic( const DataIonic& dataIonic );

    //@}


    //! @name Operators
    //@{



    DataIonic& operator=( const DataIonic& dataIonic );

    //@}

    //! @name Methods
    //@{

    //! output: show the data used for the simulation
    void showMe( std::ostream& output = std::cout );

    //@}


    //! @name Set Methods
    //@{

    //!external setup: set all the data for the simulation
    void setup( const GetPot& dataFile );

    //@}

    /*//! End time
    Real endtime() const;

    //! FE space order
    std::string wOrder() const;
*/

    bool        M_hasHeterogeneousTauClose;

    UInt        M_verbose;

    //! RogersMcCulloch 1994 Ionic Model parameters
    Real        M_a;
    Real        M_b;
    Real        M_c1;
    Real        M_c2;
    Real        M_d;
    Real        M_timeUnit;
    Real        M_potentialAmplitude;
    Real        M_restPotential;
    Real        M_initialRepolarization;
    //Mitchell & Schaeffer
    Real        M_tau_in;   // = 0.8
    Real        M_tau_out;  // = 18.0
    Real        M_tau_open; // = 300.0
    Real        M_tau_close;// = 100.0
    Real        M_criticalPotential;    // =  -67.0
    Real        M_potentialMinimum;
    Real        M_potentialMaximum;
    Real        M_reactionAmplitude;
    Real        M_initialTime;
    Real        M_tend;
    Real        M_BDForder;       //= 1

    std::string M_meshFile;
    std::string M_meshDirectory;

private:


};

}
#endif
