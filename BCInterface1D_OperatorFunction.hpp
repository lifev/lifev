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
 *  @brief BCInterface1D_OperatorFunction
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 10-05-2010
 */

#ifndef BCInterface1D_OperatorFunction_H
#define BCInterface1D_OperatorFunction_H 1

#include <lifemc/lifesolver/OneDimensionalModel_Solver.hpp>

#include <lifemc/lifesolver/BCInterface1D_Function.hpp>

namespace LifeV
{

//! BCInterface1D_OperatorFunction - LifeV bcFunction wrapper for BCInterface1D (with operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface1D, SpiritParser and a general
 *  LifeV operator (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from
 *  a GetPot file in which are present some operator parameters.
 *
 *  The class can be used in two ways:
 *
 *  1) hereditating it and implementing the template specialization of createAccessList() and UpdateOperatorVariables();
 *  2) manually setting the variables by using the setVariable() function.
 *
 *	<b>AVAILABLE OPERATORS</b>
 *
 *	Available operators are:
 *
 *	f_area
 *	f_density
 *	f_flux
 *	f_pressure
 *	f_viscosity
 *	s_density
 *	s_poisson
 *	s_thickness
 *	s_young
 */
template< class Operator >
class BCInterface1D_OperatorFunction: public virtual BCInterface1D_Function< Operator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface1D_Function< Operator >                      super;
    typedef BCInterface1D_Data                                      Data_Type;
    typedef OneDimensionalModel_Solver::Solution_PtrType            Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface1D_OperatorFunction();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface1D_OperatorFunction( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param function BCInterface1D_OperatorFunction
     */
    BCInterface1D_OperatorFunction( const BCInterface1D_OperatorFunction& function );

    //! Destructor
    virtual ~BCInterface1D_OperatorFunction() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param function BCInterface1D_OperatorFunction
     * @return reference to a copy of the class
     */
    virtual BCInterface1D_OperatorFunction& operator=( const BCInterface1D_OperatorFunction& function );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void SetData( const Data_Type& data );

    //! Set an operator
    /*!
     * @param Oper operator
     */
    inline void SetOperator( const boost::shared_ptr< Operator >& Oper );

    //! Set solution
    /*!
     * @param solution The solution container of the 1D problem
     */
    inline void SetSolution( const Solution_PtrType solution );

    //! Set variable function
    /*!
     * @param name name of the variable
     * @param value value of the variable
     */
    inline void SetVariable( const std::string& name, const Real& value );

    //! Update operator variables
    inline void UpdateOperatorVariables();

    //@}

protected:

    //! @name Protected Methods
    //@{

    inline void CreateAccessList( const Data_Type& data );

    //@}

    //List of all available operators
    enum operatorList
    {
        f_area,
        f_flux,
        f_density,
        f_pressure,
        f_viscosity,
        s_density,
        s_poisson,
        s_thickness,
        s_young,
    };

    boost::shared_ptr< Operator >         M_operator;
    Solution_PtrType                      M_solution;

    OneD_BCSide                           M_side;
    std::set< operatorList >              M_list;

private:

    //! @name Private Methods
    //@{

    inline void CreateFluidMap( std::map< std::string, operatorList >& mapList );
    inline void CreateSolidMap( std::map< std::string, operatorList >& mapList );
    inline void CreateList( const std::map< std::string, operatorList >& mapList, const Data_Type& data );

    inline void SwitchErrorMessage( const std::string& operatorType );

    //@}

};

//! Factory create function
template< typename Operator >
inline BCInterface1D_Function< Operator >* BCInterface1D_CreateOperatorFunction()
{
    return new BCInterface1D_OperatorFunction< Operator > ();
}

// ===================================================
// Constructors
// ===================================================
template< class Operator >
BCInterface1D_OperatorFunction< Operator >::BCInterface1D_OperatorFunction() :
        super                            (),
        M_operator                       (),
        M_solution                       (),
        M_side                           (),
        M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface1D_OperatorFunction::BCInterface1D_OperatorFunction()" << "\n";
#endif

}

template< class Operator >
BCInterface1D_OperatorFunction< Operator >::BCInterface1D_OperatorFunction( const Data_Type& data ) :
        super                            (),
        M_operator                       (),
        M_solution                       (),
        M_side                           (),
        M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface1D_OperatorFunction::BCInterface1D_OperatorFunction( data )" << "\n";
#endif

    this->SetData( data );
}

template< class Operator >
BCInterface1D_OperatorFunction< Operator >::BCInterface1D_OperatorFunction( const BCInterface1D_OperatorFunction& function ) :
        super                            ( function ),
        M_operator                       ( function.M_operator ),
        M_solution                       ( function.M_solution ),
        M_side                           ( function.M_side ),
        M_list                           ( function.M_list )
{
}

// ===================================================
// Methods
// ===================================================
template< class Operator >
BCInterface1D_OperatorFunction< Operator >&
BCInterface1D_OperatorFunction< Operator >::operator=( const BCInterface1D_OperatorFunction& function )
{
    if ( this != &function )
    {
        super::operator=( function );

        M_operator = function.M_operator;
        M_solution = function.M_solution;
        M_side     = function.M_side;
        M_list     = function.M_list;
    }

    return *this;
}

template< class Operator >
void
BCInterface1D_OperatorFunction< Operator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface1D_OperatorFunction::setData" << "\n";
#endif

    M_side     = data.GetSide();

    super::SetData( data );

    CreateAccessList( data );
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::SetOperator( const boost::shared_ptr< Operator >& Oper )
{
    M_operator = Oper;
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::SetSolution(  const Solution_PtrType solution )
{
    M_solution = solution;
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::SetVariable( const std::string& name, const Real& value )
{
    super::M_parser->SetVariable( name, value );
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::UpdateOperatorVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface1D_OperatorFunction<FSIOperator>::UpdateOperatorVariables  " << "\n";
#endif

    // Create/Update variables for 1D problem
    for ( typename std::set< operatorList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
            // f_ -> FLUID
        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_side) << "): " << M_operator->BoundaryValue( OneD_A, M_side ) << "\n";
#endif
            SetVariable( "f_area", M_operator->BoundaryValue( *M_solution, OneD_A, M_side ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                  f_density: " << M_operator->Physics()->Data()->DensityRho() << "\n";
#endif
            SetVariable( "f_density", M_operator->Physics()->Data()->DensityRho() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_side) << "): " << M_operator->BoundaryValue( OneD_Q, M_side ) << "\n";
#endif

            SetVariable( "f_flux", M_operator->BoundaryValue( *M_solution, OneD_Q, M_side ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_side) << "): " << M_operator->BoundaryValue( OneD_P, M_side ) << "\n";
#endif

            SetVariable( "f_pressure", M_operator->BoundaryValue( *M_solution, OneD_P, M_side ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                f_viscosity: " << M_operator->fluid().viscosity() << "\n";
#endif
            SetVariable( "f_viscosity", M_operator->Physics()->Data()->Viscosity() );

            break;

            // s_ -> SOLID
        case s_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   s_density: " << M_operator->solid().rho() << "\n";
#endif

            SetVariable( "s_density", M_operator->Physics()->Data()->DensityWall() );

            break;

        case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   s_poisson: " << M_operator->solid().poisson() << "\n";
#endif

            SetVariable( "s_poisson", M_operator->Physics()->Data()->Poisson() );

            break;

        case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                 s_thickness: " << M_operator->solid().thickness() << "\n";
#endif

            SetVariable( "s_thickness", M_operator->Physics()->Data()->Thickness() );

            break;

        case s_young:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                     s_young: " << M_operator->solid().young() << "\n";
#endif

            SetVariable( "s_young", M_operator->Physics()->Data()->Young() );

            break;

        default:
            SwitchErrorMessage( "OneDimensionalModel_Solver" );
        }
}

// ===================================================
// Protected Methods
// ===================================================
template< class Operator>
inline void
BCInterface1D_OperatorFunction< Operator >::CreateAccessList( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface1D_OperatorFunction<FSIOperator>::createAccessList" << "\n";
#endif

    std::map< std::string, operatorList > mapList;

    CreateFluidMap( mapList );
    CreateSolidMap( mapList );
    CreateList( mapList, data );

    //if ( M_operator.get() )
    //    UpdateOperatorVariables();
}

// ===================================================
// Private Methods
// ===================================================
template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::CreateFluidMap( std::map< std::string, operatorList >& mapList )
{
    mapList["f_area"]      = f_area;
    mapList["f_density"]   = f_density;
    mapList["f_flux"]      = f_flux;
    mapList["f_pressure"]  = f_pressure;
    mapList["f_viscosity"] = f_viscosity;
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::CreateSolidMap( std::map< std::string, operatorList >& mapList )
{
    mapList["s_density"]   = s_density;
    mapList["s_poisson"]   = s_poisson;
    mapList["s_thickness"] = s_thickness;
    mapList["s_young"]     = s_young;
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::CreateList( const std::map< std::string, operatorList >& mapList, const Data_Type& data )
{
    M_list.clear();
    for ( typename std::map< std::string, operatorList >::const_iterator j = mapList.begin(); j != mapList.end(); ++j )
        if ( boost::find_first( data.GetBaseString(), j->first ) )
            M_list.insert( j->second );
}

template< class Operator >
inline void
BCInterface1D_OperatorFunction< Operator >::SwitchErrorMessage( const std::string& operatorType )
{
    std::cout << "ERROR: Invalid variable type for " << operatorType << "OperatorFunction" << std::endl;
}

} // Namespace LifeV

#endif /* BCInterface1D_OperatorFunction_H */
