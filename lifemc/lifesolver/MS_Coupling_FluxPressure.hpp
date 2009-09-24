/* -*- mode: c++ -*-

  This file is part of the LifeV Applications.

  Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
       Date: 2009-08-24

  Copyright (C) 2009 EPFL

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2.1 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA
*/
/**
   \file MS_Coupling_FluxPressure.hpp
   \author Cristiano Malossi <cristiano.malossi@epfl.ch>
   \date 2009-08-24
 */

#ifndef __MS_Coupling_FluxPressure_H
#define __MS_Coupling_FluxPressure_H 1

#include <lifemc/lifefem/BCInterface.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

//! MS_Coupling_FluxPressure - Flux Pressure coupling condition
/*!
 *  The MS_Coupling_FluxPressure class is an implementation of the MS_PhysicalCoupling
 *  for applying Flux Pressure coupling conditions on the models.
 *
 *  The coupling equations are:
 *  Q_j = \sum Q_i
 *  P_i = P_j
 *  where Q is the flux and P is the total pressure.
 *
 *  @author Cristiano Malossi
 */
class MS_Coupling_FluxPressure : public MS_PhysicalCoupling
{
public:

	typedef MS_PhysicalCoupling								super;

	//! @name Constructors & Destructor
    //@{

    //! Constructor
	MS_Coupling_FluxPressure();

	//! Copy constructor
	/*!
	 * \param FluxPressure		- MS_Coupling_FluxPressure
	 */
	MS_Coupling_FluxPressure( const MS_Coupling_FluxPressure& FluxPressure );

    //! Destructor
    ~MS_Coupling_FluxPressure() {}

    //@}



    //! @name Methods
    //@{

	//! Operator=
	/*!
	 * \param fluxPressure - MS_Coupling_FluxPressure
	 */
    MS_Coupling_FluxPressure& operator=( const MS_Coupling_FluxPressure& fluxPressure );

    //@}



    //! @name MultiScale Physical Coupling
    //@{

	//! Setup the data of the coupling
    void SetupData( void );

	//! Setup the coupling
    void SetupCoupling( void );

	//! Update the values of the couplings variables (Q and P)
    void UpdateCoupling( void );

	//! Display some information about the coupling
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
	//@{

    template< class model >
    inline void imposeFlux		( const PhysicalModel_ptr& physicalModel );

    template< class model >
    inline void imposePressure	( const PhysicalModel_ptr& physicalModel, const UInt& i );

    Real functionFlux		( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ );
    Real functionPressure	( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/ );

	//@}

	BCFunctionBase 										M_baseFlux;
	BCFunctionBase 										M_basePressure;

	Real												M_flux;
	Real												M_pressure;

};

//! Factory create function
inline MS_PhysicalCoupling* createFluxPressure()
{
	return new MS_Coupling_FluxPressure();
}



// ===================================================
//! Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_FluxPressure::imposeFlux( const PhysicalModel_ptr& physicalModel )
{
	model *Model = dynamic_cast<model *>( &( *physicalModel ) );

	Model->GetBCInterface().addBC( "imposeFlux_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[0] ),
									M_flags[0], Flux, Full, M_baseFlux, 3 );
}

template< class model >
inline void
MS_Coupling_FluxPressure::imposePressure( const PhysicalModel_ptr& physicalModel, const UInt& i )
{
	model *Model = dynamic_cast<model *>( &( *physicalModel ) );

	Model->GetBCInterface().addBC( "imposePressure_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ),
									M_flags[i], Natural, Normal, M_basePressure );
}

} // Namespace LifeV

#endif /* __MS_Coupling_FluxPressure_H */
