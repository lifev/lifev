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
    @brief Set of utility for the electrophysiology benchmark

    @date 03 - 2014
    @author Simone Rossi <simone.rossi@epfl.ch>

    @contributor
    @mantainer Simone Rossi <simone.rossi@epfl.ch>
 */

#ifndef BENCHMARKUTILITY_HPP_
#define BENCHMARKUTILITY_HPP_


#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


// ---------------------------------------------------------------
// Include the ionic models
// ---------------------------------------------------------------

#include <lifev/electrophysiology/solver/IonicModels/IonicMinimalModel.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicLuoRudyI.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicTenTusscher06.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicHodgkinHuxley.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicNoblePurkinje.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicFox.hpp>
#include <lifev/electrophysiology/solver/IonicModels/IonicAlievPanfilov.hpp>

// ---------------------------------------------------------------
// Include LifeV
// ---------------------------------------------------------------

#include <lifev/core/LifeV.hpp>

// ---------------------------------------------------------------
// In LifeV namespace we create the benchmarkUtility namespace
// ---------------------------------------------------------------

namespace LifeV{

namespace BenchmarkUtility{


// ---------------------------------------------------------------
// We define some useful typedefs
// ---------------------------------------------------------------

typedef ElectroIonicModel                                        ionicModel_Type;
typedef boost::shared_ptr<ionicModel_Type>                       ionicModelPtr_Type;
typedef boost::function < Real (const Real& t,
                                 const Real& x,
                                 const Real& y,
                                 const Real& z,
                                 const ID&   i ) >   function_Type;

// ---------------------------------------------------------------
// This function is used to choose among the ionic models with an
// if among the ionic model names.
// ---------------------------------------------------------------

Real chooseIonicModel(ionicModelPtr_Type& model, std::string ionic_model, Epetra_Comm& Comm )
{
	Real activationThreshold(0.95);
	bool ionicModelSet = false;
    if ( Comm.MyPID() == 0 )
    {
        std::cout << "Constructing the ionic model ... !"; // << std::endl;
    }

    if ( ionic_model == "AlievPanfilov" )
    {
        model.reset ( new IonicAlievPanfilov() );
        ionicModelSet = true;
    }
    if ( ionic_model == "LuoRudyI" )
    {
        model.reset ( new IonicLuoRudyI() );
        ionicModelSet = true;
    }
    if ( ionic_model == "TenTusscher06")
    {
        model.reset (new IonicTenTusscher06() );
        ionicModelSet = true;
    }
    if ( ionic_model == "HodgkinHuxley")
    {
        model.reset (new IonicHodgkinHuxley() );
        activationThreshold = 10.0;
        ionicModelSet = true;
    }
    if ( ionic_model == "NoblePurkinje")
    {
        model.reset (new IonicNoblePurkinje() );
        ionicModelSet = true;
    }
    if ( ionic_model == "MinimalModel")
    {
        model.reset ( new IonicMinimalModel() );
        ionicModelSet = true;
    }
    if ( ionic_model == "Fox")
    {
        model.reset ( new IonicFox() );
        if ( Comm.MyPID() == 0 )
        {
//            assert(0 && "Fox model is not properly working in 3D."); //TO DO: Fix It!
        }
    }

    if ( Comm.MyPID() == 0 )
    {
        std::cout << " Done!" << std::endl;
    }

    assert(ionicModelSet  && "Ionic model not specified.");

    model -> showMe();

    return activationThreshold;
}

// ---------------------------------------------------------------
// This function is used to pace the minimal model (and the
// Aliev Panfilov model )
// ---------------------------------------------------------------

Real PacingProtocolMM ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 10;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

// ---------------------------------------------------------------
// This function is used to pace the Hodgkin-Huxley model
// ---------------------------------------------------------------

Real PacingProtocolHH ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 500.;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            &&
            t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

// ---------------------------------------------------------------
// This function is used to pace the cardiac models
// ---------------------------------------------------------------

Real PacingProtocol ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   /*id*/)
{

    Real pacingSite_X = 0.0;
    Real pacingSite_Y = 0.0;
    Real pacingSite_Z = 0.0;
    Real stimulusRadius = 0.15;
    Real stimulusValue = 40.0;

    Real returnValue;

    if ( std::abs ( x - pacingSite_X ) <= stimulusRadius
            &&
            std::abs ( z - pacingSite_Z ) <= stimulusRadius
            &&
            std::abs ( y - pacingSite_Y ) <= stimulusRadius
            && t >= 0.1 && t <= 2)
    {
        returnValue = stimulusValue;
    }
    else
    {
        returnValue = 0.;
    }

    return returnValue;
}

// ---------------------------------------------------------------
// This function is used to set the external stimulus
// ---------------------------------------------------------------

void setStimulus(function_Type& f, std::string ionic_model)
{
    if (ionic_model == "MinimalModel" || ionic_model == "AlievPanfilov")
    {
        f = &BenchmarkUtility::PacingProtocolMM;
    }
    else if (ionic_model == "HodgkinHuxley" )
    {
        f = &BenchmarkUtility::PacingProtocolHH;
    }
    else
    {
        f = &BenchmarkUtility::PacingProtocol;
    }
}

}//BenchmarkUtility

} //LifeV


#endif /* BENCHMARKUTILITY_HPP_ */
