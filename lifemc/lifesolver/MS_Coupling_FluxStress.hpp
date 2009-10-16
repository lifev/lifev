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
 \file MS_Coupling_FluxStress.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-08-24
 */

#ifndef __MS_Coupling_FluxStress_H
#define __MS_Coupling_FluxStress_H 1
#define DEBUG 1;
#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

//! MS_Coupling_FluxStress - Flux-Stress coupling condition
/*!
 *  The MS_Coupling_FluxStress class is an implementation of the MS_PhysicalCoupling
 *  for applying Flux-Stress coupling conditions on the models.
 *
 *  The coupling equations are:
 *  Q_j = - \sum Q_i
 *  \sigma_i = -P_j
 *  where Q is the flux and P is the pressure (or the total pressure).
 *
 *  @author Cristiano Malossi
 */
class MS_Coupling_FluxStress: public MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_FluxStress();

    //! Copy constructor
    /*!
     * \param FluxStress - MS_Coupling_FluxStress
     */
    MS_Coupling_FluxStress( const MS_Coupling_FluxStress& FluxStress );

    //! Destructor
    ~MS_Coupling_FluxStress() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param FluxStress - MS_Coupling_FluxStress
     */
    MS_Coupling_FluxStress& operator=( const MS_Coupling_FluxStress& FluxStress );

    //@}


    //! @name MultiScale Physical Coupling
    //@{

    //! Setup the data of the coupling
    void SetupData( void );

    //! Setup the coupling
    void SetupCoupling( void );

    //! Setup parameters for the implicit coupling
    void SetupImplicitCoupling( ContainerOfVectors< EpetraVector >& couplingVariables,
                                ContainerOfVectors< EpetraVector >& couplingResiduals );

    //! Update the values of the coupling variables (Flux and Normal Stress i.e. -P)
    void UpdateCouplingVariables( void );

    //! Update the values of the coupling residual: R = C - F(C)
    void UpdateCouplingResiduals( void );

    //! Display some information about the coupling
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
    //@{

    EpetraVector CouplingFunction( void );

    template< class model >
    inline void imposeFlux( const UInt& i = 0 );

    template< class model >
    inline void imposeStress( const UInt& i );

    template< class model >
    inline Real computeFlux( const UInt& i );

    template< class model >
    inline Real computeStress( const UInt& i = 0 );

    template< class model >
    inline Real computeStaticPressure( const UInt& i );

    template< class model >
    inline Real computeDynamicPressure( const UInt& i );

    Real functionFlux  ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);
    Real functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);

    //@}

    enum pressureTypes
    {
        Static,
        Total
    };

    BCFunctionBase M_baseFlux;
    BCFunctionBase M_baseStress;

    static std::map< std::string, pressureTypes > M_mapPressure;
    pressureTypes                                 M_pressureType;
};

//! Factory create function
inline MS_PhysicalCoupling* createFluxStress()
{
    return new MS_Coupling_FluxStress();
}

// ===================================================
//! Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_FluxStress::imposeFlux( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    Model->GetBCInterface().addBC( "imposeFlux_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Flux, Full, M_baseFlux, 3 );
}

template< class model >
inline void
MS_Coupling_FluxStress::imposeStress( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    Model->GetBCInterface().addBC( "imposeStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Natural, Normal, M_baseStress );
}

template< class model >
inline Real
MS_Coupling_FluxStress::computeFlux( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return Model->GetFlux( M_flags[i] );
}

template< class model >
inline Real
MS_Coupling_FluxStress::computeStress( const UInt& i )
{
    switch ( M_pressureType )
    {
        case Static:
        {
            return -computeStaticPressure< model > ( i );
        }

        case Total:
        {
            return -computeStaticPressure< model > ( i ) + ( computeFlux< model >( i ) > 0.0 ) * computeDynamicPressure< model > ( i );
        }

        default:

            std::cout << "ERROR: Invalid pressure type [" << Enum2String( M_pressureType, M_mapPressure ) << "]" << std::endl;

            return 0.0;
    }
}

template< class model >
inline Real
MS_Coupling_FluxStress::computeStaticPressure( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return Model->GetPressure( M_flags[i] );
}

template< class model >
inline Real
MS_Coupling_FluxStress::computeDynamicPressure( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return 0.5 * Model->GetDensity( M_flags[i] ) * Model->GetFlux( M_flags[i] ) * Model->GetFlux( M_flags[i] ) / ( Model->GetArea( M_flags[i] ) * Model->GetArea( M_flags[i] ) );
}

} // Namespace LifeV

#endif /* __MS_Coupling_FluxStress_H */
