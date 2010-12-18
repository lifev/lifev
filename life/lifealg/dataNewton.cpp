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
 *  @file
 *  @brief File containing a class for handling temporal discretization with GetPot.
 *
 *  @date 01-10-2003
 *  @author Miguel Angel Fernandez
 *
 *  @contributor Paolo Tricerri <paolo.tricerri@epfl.ch>
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#include <life/lifecore/life.hpp>
#include <life/lifealg/dataNewton.hpp>

namespace LifeV
{
  //====================================
  // Constructor & Destructor
  //===================================

DataNewton::DataNewton( const GetPot& dfile, const std::string& section )
{
    M_maxiter = dfile( ( section + "/maxiter" ).data(), 100 );
    M_abstol = dfile( ( section + "/abstol" ).data(), 0.0 );
    M_reltol = dfile( ( section + "/reltol" ).data(), 0.0 );
    M_etamax = dfile( ( section + "/etamax" ).data(), 1.e-3 );
    M_linesearch = dfile( ( section + "/linesearch" ).data(), 2 );
}

DataNewton::~DataNewton()
{}

//=====================================
// Get Methods
//=====================================

// The max number of interations
const UInt DataNewton::getMaxiter() const
{
    return M_maxiter;
}

// The absolute tolerance
const Real DataNewton::getAbstol() const
{
    return M_abstol;
}

// The relative tolerance
const Real DataNewton::getReltol() const
{
    return M_reltol;
}

// The relative tolerance
const Real DataNewton::getEtamax() const
{
    return M_etamax;
}

// The linesearch option
const UInt DataNewton::getLinesearch() const
{
    return M_linesearch;
}

// Output
void DataNewton::showMe( std::ostream& c ) const
{
    //
    c << "maxiter        = " << M_maxiter << std::endl;
    c << "abstol         = " << M_abstol << std::endl;
    c << "reltol         = " << M_reltol << std::endl;
    c << "etamax         = " << M_reltol << std::endl;
    c << "linesearch     = " << M_linesearch << std::endl;
}
}
