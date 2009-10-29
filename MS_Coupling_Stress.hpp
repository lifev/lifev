/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-10-20

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
 \file MS_Coupling_Stress.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-10-20
 */

#ifndef __MS_Coupling_Stress_H
#define __MS_Coupling_Stress_H 1

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

//! MS_Coupling_Stress - Stress coupling condition
/*!
 *  The MS_Coupling_Stress class is an implementation of the MS_PhysicalCoupling
 *  for applying Stress coupling conditions on the models.
 *
 *  @author Cristiano Malossi
 */
class MS_Coupling_Stress: public virtual MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_Stress();

    //! Copy constructor
    /*!
     * \param Stress - MS_Coupling_Stress
     */
    MS_Coupling_Stress( const MS_Coupling_Stress& Stress );

    //! Destructor
    ~MS_Coupling_Stress() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param Stress - MS_Coupling_Stress
     */
    MS_Coupling_Stress& operator=( const MS_Coupling_Stress& Stress );

    //@}


    //! @name MultiScale Physical Coupling
    //@{

    //! Setup the data of the coupling
    void SetupData( void );

    //! Setup the coupling
    void SetupCoupling( void );

    //! Initialize the values of the coupling variables
    void InitializeCouplingVariables( void );

    //! Update the values of the coupling residuals
    void ExportCouplingResiduals( VectorType& CouplingResiduals );

    //! Compute Jacobian product: J(X) * X
    /*!
     * \param deltaCouplingVariables - variation of the coupling variables
     */
    void ComputeJacobianProduct( const VectorType& deltaCouplingVariables );

    //! Display some information about the coupling
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
    //@{

    EpetraVector CouplingFunction( void );

    //! Print the couplings quantities
    void PrintQuantities( void );



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

    Real functionStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);



    template< class model >
    inline void imposeDeltaStress( const UInt& i );

    template< class model >
    inline Real computeDeltaFlux( const UInt& i );

    Real functionDeltaStress( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);

    //@}

    enum stressTypes
    {
        StaticPressure,
        TotalPressure
    };

    BCFunctionBase M_baseStress;

    BCFunctionBase M_baseDeltaStress;

    static std::map< std::string, stressTypes > M_mapStress;
    stressTypes                                 M_stressType;
};

//! Factory create function
inline MS_PhysicalCoupling* createStress()
{
    return new MS_Coupling_Stress();
}

// ===================================================
//! Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_Stress::imposeStress( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    Model->GetBC().addBC( "imposeStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Natural, Normal, M_baseStress );
}

template< class model >
inline Real
MS_Coupling_Stress::computeFlux( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return Model->GetFlux( M_flags[i] );
}

template< class model >
inline Real
MS_Coupling_Stress::computeStress( const UInt& i )
{
    switch ( M_stressType )
    {
        case StaticPressure:
        {
            return -computeStaticPressure< model > ( i );
        }

        case TotalPressure:
        {
            return -computeStaticPressure< model > ( i ) + ( computeFlux< model >( i ) > 0.0 ) * computeDynamicPressure< model > ( i );
        }

        default:

            std::cout << "ERROR: Invalid pressure type [" << Enum2String( M_stressType, M_mapStress ) << "]" << std::endl;

            return 0.0;
    }
}

template< class model >
inline Real
MS_Coupling_Stress::computeStaticPressure( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return Model->GetPressure( M_flags[i] );
}

template< class model >
inline Real
MS_Coupling_Stress::computeDynamicPressure( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    return 0.5 * Model->GetDensity( M_flags[i] ) * Model->GetFlux( M_flags[i] ) * Model->GetFlux( M_flags[i] ) / ( Model->GetArea( M_flags[i] ) * Model->GetArea( M_flags[i] ) );
}

template< class model >
inline void
MS_Coupling_Stress::imposeDeltaStress( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    Model->GetLinearBC().addBC( "imposeDeltaStress_Model_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ), M_flags[i], Natural, Normal, M_baseDeltaStress );
}

template< class model >
inline Real
MS_Coupling_Stress::computeDeltaFlux( const UInt& i )
{
    model *Model = dynamic_cast< model * > ( &( *M_models[i] ) );

    Model->UpdateLinearSystem();
    Model->SolveLinearSystem();

    return Model->GetLinearFlux( M_flags[i] );
}

} // Namespace LifeV

#endif /* __MS_Coupling_Stress_H */
