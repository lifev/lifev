//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @brief BCInterface_DefaultFunctions
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-05-2010
 */

#ifndef BCInterface1D_DefaultFunctions_H
#define BCInterface1D_DefaultFunctions_H 1

#include <lifemc/lifesolver/BCInterface1D_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface1D_Data.hpp>

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

namespace LifeV
{

//! BCInterface1D_DefaultFunctions - Interface with the default boundary conditions of the OneDimensionalModel
/*!
 *  @author Cristiano Malossi
 *
 *  The BCInterface1D_DefaultFunctions class provides a general interface between the
 *  BCInterface1D and the default BC for the One Dimensional Model.
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the ID of the BC to impose,
 *  and the Operator. The list of available conditions is described by DefaultFunctions. They are:
 *
 *	- Riemann,
 *	- Compatibility,
 *	- Absorbing,
 *	- Resistance
 *
 *	To get the base for the boundary condition, call the GetBase function.
 */
template< class Operator >
class BCInterface1D_DefaultFunctions
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_Data                                                    Data_Type;

    typedef OneDimensionalModel_BC::BCFunction_Type                               BCFunction_Type;
    typedef OneDimensionalModel_BC::BCFunction_PtrType                            BCFunction_PtrType;
    typedef OneDimensionalModel_BC::BCFunction_Default_PtrType                    BCFunction_Default_PtrType;

    typedef OneDimensionalModel_BC::Flux_PtrType                                  Flux_PtrType;
    typedef OneDimensionalModel_BC::Source_PtrType                                Source_PtrType;
    typedef OneDimensionalModel_BC::Solution_PtrType                              Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_DefaultFunctions();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_DefaultFunctions( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param Default BCInterface1D_DefaultFunctions
     */
    BCInterface1D_DefaultFunctions( const BCInterface1D_DefaultFunctions& Default );

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
    BCInterface1D_DefaultFunctions& operator=( const BCInterface1D_DefaultFunctions& Default );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void SetData( const Data_Type& data );

    //! Set solution
    /*!
     * @param solution solution container of the 1D model
     */
    inline void SetSolution( const Solution_PtrType solution );

    //! Set flux and source
    /*!
     * @param flux flux object of the 1D model
     * @param source source object of the 1D model
     */
    inline void SetFluxSource( const Flux_PtrType& flux, const Source_PtrType& source );

    //@}


    //! @name Get functions
    //@{

    //! Get the base of the boundary condition
    BCFunction_Type& GetBase()
    {
        return *M_base;
    }

    //@}

private:

    enum DefaultFunctions
    {
        Riemann,
        Compatibility,
        Absorbing,
        Resistance
    };

    BCFunction_PtrType          M_base;
    BCFunction_Default_PtrType  M_defaultFunction;
};

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface1D_DefaultFunctions< Operator >::BCInterface1D_DefaultFunctions() :
        M_base              ( new BCFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::BCInterface1D_DefaultFunctions()" << "\n";
#endif

}

template< class Operator >
BCInterface1D_DefaultFunctions< Operator >::BCInterface1D_DefaultFunctions( const Data_Type& data ) :
        M_base              ( new BCFunction_Type() ),
        M_defaultFunction   ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::BCInterface1D_DefaultFunctions( data )" << "\n";
#endif

    this->SetData( data );
}

template< class Operator >
BCInterface1D_DefaultFunctions< Operator >::BCInterface1D_DefaultFunctions( const BCInterface1D_DefaultFunctions& Default ) :
        M_base              ( Default.M_base ),
        M_defaultFunction   ()
{
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
BCInterface1D_DefaultFunctions< Operator >&
BCInterface1D_DefaultFunctions< Operator >::operator=( const BCInterface1D_DefaultFunctions& Default )
{
    if ( this != &Default )
    {
        M_base        = Default.M_base;
    }

    return *this;
}

template< class Operator >
void BCInterface1D_DefaultFunctions< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::setData" << "\n";
#endif

    //Set mapFunction
    std::map< std::string, DefaultFunctions > mapFunction;
    mapFunction["Riemann"]             = Riemann;
    mapFunction["Compatibility"]       = Compatibility;
    mapFunction["Absorbing"]           = Absorbing;
    mapFunction["Resistance"]          = Resistance;

    switch ( mapFunction[data.GetBaseString()] )
    {
    case Riemann:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Riemann( data.GetSide(),
                                                                             data.GetQuantity()
                                                                           ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Compatibility:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Compatibility( data.GetSide(),
                                                                                   data.GetQuantity()
                                                                                 ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Absorbing:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface1D_DefaultFunctions::checkFunction                          Absorbing" << "\n";
#endif

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Absorbing( data.GetSide(),
                                                                               data.GetQuantity()
                                                                             ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Absorbing::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Absorbing *> ( &( *M_defaultFunction ) ), _1, _2 ) );

        break;

    case Resistance:

        M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Resistance( data.GetSide(),
                                                                                data.GetQuantity(),
                                                                                data.GetResistance()[0]
                                                                              ) );

        M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Resistance::operator(),
                                          dynamic_cast<OneDimensionalModel_BCFunction_Resistance *> ( &( *M_defaultFunction ) ), _1, _2 ) );


        break;

    }
}

template< class Operator >
inline void
BCInterface1D_DefaultFunctions< Operator >::SetSolution( const Solution_PtrType solution )
{
    M_defaultFunction->setSolution( solution );
}

template< class Operator >
inline void
BCInterface1D_DefaultFunctions< Operator >::SetFluxSource( const Flux_PtrType& flux, const Source_PtrType& source )
{
    M_defaultFunction->setFluxSource( flux, source );
}

} // Namespace LifeV

#endif /* BCInterface1D_DefaultFunctions_H */
