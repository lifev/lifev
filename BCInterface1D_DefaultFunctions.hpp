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
 *  @brief File containing the BCInterface1D_DefaultFunctions class
 *
 *  @date 10-05-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface1D_DefaultFunctions_H
#define BCInterface1D_DefaultFunctions_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface1D_Data.hpp>

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV
{

//! BCInterface1D_DefaultFunctions - Interface with the default boundary conditions of the 1D model
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1D_DefaultFunctions class provides a general interface between the
 *  BCInterface1D and the default BC for the 1D Model.
 *
 *  <b>DETAILS:</b>
 *
 *  The list of available conditions is described by the \c defaultFunctions enum type.
 *  They are:
 *
 *	- Riemann,
 *	- Compatibility,
 *	- Absorbing,
 *	- Resistance
 *
 *	To get the base for the boundary condition, call the \c base() method.
 */
template< class PhysicalSolver >
class BCInterface1D_DefaultFunctions
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_Data                                                    data_Type;

    typedef OneDimensionalModel_BC::BCFunction_Type                               bcFunction_Type;
    typedef OneDimensionalModel_BC::BCFunction_PtrType                            bcFunction_PtrType;
    typedef OneDimensionalModel_BC::BCFunction_Default_PtrType                    bcFunction_Default_PtrType;

    typedef OneDimensionalModel_BC::Flux_PtrType                                  fluxPtr_Type;
    typedef OneDimensionalModel_BC::Source_PtrType                                sourcePtr_Type;
    typedef OneDimensionalModel_BC::Solution_PtrType                              solutionPtr_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_DefaultFunctions();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_DefaultFunctions( const data_Type& data );

    //! Copy constructor
    /*!
     * @param Default BCInterface1D_DefaultFunctions
     */
    BCInterface1D_DefaultFunctions( const BCInterface1D_DefaultFunctions& defaultFunctions );

    //! Destructor
    ~BCInterface1D_DefaultFunctions() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param Default BCInterface1D_DefaultFunctions
     * @return reference to a copy of the class
     */
    BCInterface1D_DefaultFunctions& operator=( const BCInterface1D_DefaultFunctions& defaultFunctions );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void setData( const data_Type& data );

    //! Set solution
    /*!
     * @param solution solution container of the 1D model
     */
    void setSolution( const solutionPtr_Type solution ) { M_defaultFunction->setSolution( solution ); }

    //! Set flux and source
    /*!
     * @param flux flux object of the 1D model
     * @param source source object of the 1D model
     */
    void setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source ) { M_defaultFunction->setFluxSource( flux, source ); }

    //@}


    //! @name Get functions
    //@{

    //! Get the base of the boundary condition
    /*!
     * @return boundary condition base
     */
    bcFunction_Type& base() { return *M_base; }

    //@}

private:

    enum defaultFunctions
    {
        Riemann,
        Compatibility,
        Absorbing,
        Resistance
    };

    bcFunction_PtrType          M_base;
    bcFunction_Default_PtrType  M_defaultFunction;
};

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolver >
BCInterface1D_DefaultFunctions< PhysicalSolver >::BCInterface1D_DefaultFunctions() :
        M_base              ( new bcFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::BCInterface1D_DefaultFunctions()" << "\n";
#endif

}

template< class PhysicalSolver >
BCInterface1D_DefaultFunctions< PhysicalSolver >::BCInterface1D_DefaultFunctions( const data_Type& data ) :
        M_base              ( new bcFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::BCInterface1D_DefaultFunctions( data )" << "\n";
#endif

    this->setData( data );
}

template< class PhysicalSolver >
BCInterface1D_DefaultFunctions< PhysicalSolver >::BCInterface1D_DefaultFunctions( const BCInterface1D_DefaultFunctions& defaultFunctions ) :
        M_base              ( defaultFunctions.M_base ),
        M_defaultFunction   ()
{
}

// ===================================================
// Methods
// ===================================================
template< class PhysicalSolver >
BCInterface1D_DefaultFunctions< PhysicalSolver >&
BCInterface1D_DefaultFunctions< PhysicalSolver >::operator=( const BCInterface1D_DefaultFunctions& defaultFunctions )
{
    if ( this != &defaultFunctions )
    {
        M_base            = defaultFunctions.M_base;
        M_defaultFunction = defaultFunctions.M_defaultFunctions;
    }

    return *this;
}

template< class PhysicalSolver >
inline void
BCInterface1D_DefaultFunctions< PhysicalSolver >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, defaultFunctions > mapFunction;
    mapFunction["Riemann"]             = Riemann;
    mapFunction["Compatibility"]       = Compatibility;
    mapFunction["Absorbing"]           = Absorbing;
    mapFunction["Resistance"]          = Resistance;

    switch ( mapFunction[data.baseString()] )
    {
    case Riemann:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Riemann( data.side(), data.quantity() ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Compatibility:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Compatibility( data.side(), data.quantity() ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Absorbing:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface1D_DefaultFunctions::checkFunction                          Absorbing" << "\n";
#endif

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Absorbing( data.side(), data.quantity() ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Absorbing::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Absorbing *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Resistance:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Resistance( data.side(), data.quantity(), data.resistance()[0] ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Resistance::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Resistance *> ( &( *M_defaultFunction ) ), _1, _2 ) );


        break;

    }
}

} // Namespace LifeV

#endif /* BCInterface1D_DefaultFunctions_H */
