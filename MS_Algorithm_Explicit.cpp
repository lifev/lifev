//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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
 *  @brief MultiScale Explicit Algorithm
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 26-10-2009
 */

#include <lifemc/lifesolver/MS_Algorithm_Explicit.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
MS_Algorithm_Explicit::MS_Algorithm_Explicit() :
        super               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Explicit::MS_Algorithm_Explicit() \n";
#endif

    M_type = Explicit;
}

MS_Algorithm_Explicit::MS_Algorithm_Explicit( const MS_Algorithm_Explicit& algorithm ) :
        super               ( algorithm )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8011 ) << "MS_Algorithm_Explicit::MS_Algorithm_Explicit( algorithm ) \n";
#endif

}

// ===================================================
// Operators
// ===================================================
MS_Algorithm_Explicit&
MS_Algorithm_Explicit::operator=( const MS_Algorithm_Explicit& algorithm )
{
    if ( this != &algorithm )
    {
        super::operator=( algorithm );
    }
    return *this;
}

// ===================================================
// MultiScale Algorithm Virtual Methods
// ===================================================
void
MS_Algorithm_Explicit::SetupData( const std::string& FileName )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MS_Algorithm_Explicit::SetupData( DataFile ) \n";
#endif

    super::SetupData( FileName );
}

void
MS_Algorithm_Explicit::SubIterate()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 8012 ) << "MS_Algorithm_Explicit::SubIterate( tolerance, subITMax ) \n";
#endif

    ToleranceSatisfied();
}

void
MS_Algorithm_Explicit::UpdateCouplingVariables()
{
    // We use the initialize method for updating the coupling
    M_multiscale->InitializeCouplingVariables();
}

void
MS_Algorithm_Explicit::ShowMe()
{
    if ( M_displayer->isLeader() )
    {
        super::ShowMe();
    }
}

} // Namespace LifeV
