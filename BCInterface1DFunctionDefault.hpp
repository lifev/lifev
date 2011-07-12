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
 *  @brief File containing the BCInterface1DFunctionDefault class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1DFunctionDefault_H
#define BCInterface1DFunctionDefault_H 1

#include <lifemc/lifesolver/BCInterfaceDefinitions.hpp>
#include <lifemc/lifesolver/BCInterfaceData.hpp>

#include <lifemc/lifesolver/OneDimensionalSolver.hpp>

namespace LifeV
{

//! BCInterface1DFunctionDefault - LifeV boundary condition function wrapper for \c BCInterface1D
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1DFunctionDefault class provides a general interface between the
 *  \c BCInterface1D and the default boundary condition for the \c OneDimensionalSolver.
 *
 *  <b>DETAILS:</b> <BR>
 *  The list of available conditions is described by the \c defaultFunctions enum type.
 *
 *  They are:
 *  <ol>
 *      <li> Riemann;
 *      <li> Compatibility;
 *	    <li> Absorbing;
 *	    <li> Resistance;
 *	</ol>
 */
template< class PhysicalSolverType >
class BCInterface1DFunctionDefault
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterfaceData                                                       data_Type;

    typedef OneDimensionalBC                                                      bc_Type;
    typedef bc_Type::bcFunctionDefaultPtr_Type                                    bcFunction_Default_PtrType;

    typedef bc_Type::vectorPtrContainer_Type                                      vectorPtrContainer_Type;

    typedef bc_Type::fluxPtr_Type                                                 fluxPtr_Type;
    typedef bc_Type::sourcePtr_Type                                               sourcePtr_Type;
    typedef bc_Type::solutionPtr_Type                                             solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface1DFunctionDefault();

    //! Constructor
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    explicit BCInterface1DFunctionDefault( const data_Type& data );

    //! Destructor
    virtual ~BCInterface1DFunctionDefault() {}

    //@}


    //! @name Methods
    //@{

    //! Assign the function to the base
    /*!
     * @param base base of the bc
     */
    void assignFunction( OneDimensionalBCFunction& base );

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void setData( const data_Type& data );

    //! Set flux and source
    /*!
     * @param flux flux object of the 1D model
     * @param source source object of the 1D model
     */
    void setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source ) { M_function->setFluxSource( flux, source ); }

    //! Set solution
    /*!
     * @param solution solution container of the 1D model
     */
    void setSolution( const solutionPtr_Type& solution ) { M_function->setSolution( solution ); }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface1DFunctionDefault( const BCInterface1DFunctionDefault& defaultFunctions );

    BCInterface1DFunctionDefault& operator=( const BCInterface1DFunctionDefault& defaultFunctions );

    //@}

    enum defaultFunctions
    {
        Riemann,
        Compatibility,
        Absorbing,
        Resistance
    };

    defaultFunctions            M_defaultFunction;
    bcFunction_Default_PtrType  M_function;
};

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface1DFunctionDefault< PhysicalSolverType >::BCInterface1DFunctionDefault() :
        M_defaultFunction (),
        M_function        ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::BCInterface1DFunctionDefault()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface1DFunctionDefault< PhysicalSolverType >::BCInterface1DFunctionDefault( const data_Type& data ) :
        M_defaultFunction (),
        M_function        ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::BCInterface1DFunctionDefault( data )" << "\n";
#endif

    this->setData( data );
}


// ===================================================
// Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface1DFunctionDefault< PhysicalSolverType >::assignFunction( OneDimensionalBCFunction& base )
{
    switch ( M_defaultFunction )
    {
    case Riemann:

        base.setFunction( boost::bind( &OneDimensionalBCFunctionRiemann::operator(),
                                       dynamic_cast<OneDimensionalBCFunctionRiemann *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Compatibility:

        base.setFunction( boost::bind( &OneDimensionalBCFunctionCompatibility::operator(),
                                       dynamic_cast<OneDimensionalBCFunctionCompatibility *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Absorbing:

        base.setFunction( boost::bind( &OneDimensionalBCFunctionAbsorbing::operator(),
                                       dynamic_cast<OneDimensionalBCFunctionAbsorbing *> ( &( *M_function ) ), _1, _2 ) );

        break;

    case Resistance:

        base.setFunction( boost::bind( &OneDimensionalBCFunctionResistance::operator(),
                                       dynamic_cast<OneDimensionalBCFunctionResistance *> ( &( *M_function ) ), _1, _2 ) );

        break;
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface1DFunctionDefault< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::setData( data )" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, defaultFunctions > mapFunction;
    mapFunction["Riemann"]             = Riemann;
    mapFunction["Compatibility"]       = Compatibility;
    mapFunction["Absorbing"]           = Absorbing;
    mapFunction["Resistance"]          = Resistance;

    M_defaultFunction = mapFunction[data.baseString()];

    switch ( M_defaultFunction )
    {
    case Riemann:

        M_function.reset( new OneDimensionalBCFunctionRiemann( data.side(), data.quantity() ) );

        break;

    case Compatibility:

        M_function.reset( new OneDimensionalBCFunctionCompatibility( data.side(), data.quantity() ) );

        break;

    case Absorbing:

        M_function.reset( new OneDimensionalBCFunctionAbsorbing( data.side(), data.quantity() ) );

        break;

    case Resistance:

        M_function.reset( new OneDimensionalBCFunctionResistance( data.side(), data.quantity(), data.resistance()[0] ) );

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface1DFunctionDefault_H */
