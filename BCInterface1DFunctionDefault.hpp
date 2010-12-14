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

#include <lifemc/lifesolver/BCInterface1DDefinitions.hpp>
#include <lifemc/lifesolver/BCInterface1DData.hpp>

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV
{

//! BCInterface1DFunctionDefault - Interface with the default boundary conditions of the 1D model
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1DFunctionDefault class provides a general interface between the
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
template< class PhysicalSolverType >
class BCInterface1DFunctionDefault
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface1DData                                                     data_Type;

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
    explicit BCInterface1DFunctionDefault();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface1DFunctionDefault( const data_Type& data );

    //! Destructor
    virtual ~BCInterface1DFunctionDefault() {}

    //@}


    //! @name Set Methods
    //@{

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


    //! @name Get Methods
    //@{

    //! Get the base of the boundary condition
    /*!
     * @return boundary condition base
     */
    bcFunction_Type& base() { return *M_base; }

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

    bcFunction_PtrType          M_base;
    bcFunction_Default_PtrType  M_defaultFunction;
};

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface1DFunctionDefault< PhysicalSolverType >::BCInterface1DFunctionDefault() :
        M_base              ( new bcFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::BCInterface1DFunctionDefault()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface1DFunctionDefault< PhysicalSolverType >::BCInterface1DFunctionDefault( const data_Type& data ) :
        M_base              ( new bcFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::BCInterface1DFunctionDefault( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface1DFunctionDefault< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1DFunctionDefault::setData" << "\n";
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

#endif /* BCInterface1DFunctionDefault_H */
