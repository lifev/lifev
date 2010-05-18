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

namespace LifeV {

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

    typedef BCInterface1D_Data< Operator >                                        Data_Type;

    typedef OneDimensionalModel_BC::BCFunction_Type                               BCFunction_Type;
    typedef OneDimensionalModel_BC::BCFunction_PtrType                            BCFunction_PtrType;
    typedef OneDimensionalModel_BC::BCFunction_Default_PtrType                    BCFunction_Default_PtrType;

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

#ifdef DEBUG
    Debug( 5025 ) << "BCInterface1D_DefaultFunctions::BCInterface1D_DefaultFunctions()" << "\n";
#endif

}

template< class Operator >
BCInterface1D_DefaultFunctions< Operator >::BCInterface1D_DefaultFunctions( const Data_Type& data ) :
    M_base              ( new BCFunction_Type() ),
    M_defaultFunction   ()
{

#ifdef DEBUG
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

#ifdef DEBUG
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

            M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Riemann( data.GetOperator()->Flux(),
                                                                                 data.GetOperator()->Source(),
                                                                                 data.GetOperator()->U_thistime(),
                                                                                 data.GetSide(),
                                                                                 data.GetType()
                             ) );

            M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                              dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunction ) ), _1 ) );

            break;

        case Compatibility:

            M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Compatibility( data.GetOperator()->Flux(),
                                                                                       data.GetOperator()->Source(),
                                                                                       data.GetOperator()->U_thistime(),
                                                                                       data.GetSide(),
                                                                                       data.GetType()
                             ) );

            M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                              dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunction ) ), _1 ) );

            break;

        case Absorbing:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface1D_DefaultFunctions::checkFunction                          Absorbing" << "\n";
#endif

            M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Absorbing( data.GetOperator()->Flux(),
                                                                                   data.GetOperator()->Source(),
                                                                                   data.GetOperator()->U_thistime(),
                                                                                   data.GetSide(),
                                                                                   data.GetType()
                             ) );

            M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Absorbing::operator(),
                                              dynamic_cast<OneDimensionalModel_BCFunction_Absorbing *> ( &( *M_defaultFunction ) ), _1 ) );

            break;

        case Resistance:

            M_defaultFunction.reset( new OneDimensionalModel_BCFunction_Resistance( data.GetOperator()->Flux(),
                                                                                    data.GetOperator()->Source(),
                                                                                    data.GetOperator()->U_thistime(),
                                                                                    data.GetSide(),
                                                                                    data.GetType(),
                                                                                    0. //Resistance value (add)
                             ) );

            M_base->setFunction( boost::bind( &OneDimensionalModel_BCFunction_Resistance::operator(),
                                              dynamic_cast<OneDimensionalModel_BCFunction_Resistance *> ( &( *M_defaultFunction ) ), _1 ) );


            break;

    }
}

} // Namespace LifeV

#endif /* BCInterface1D_DefaultFunctions_H */
