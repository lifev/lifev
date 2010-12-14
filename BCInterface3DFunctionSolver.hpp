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
 *  @brief File containing the BCInterface_OperatorFunction class
 *
 *  @date 24-08-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface_OperatorFunction_H
#define BCInterface_OperatorFunction_H 1

#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifesolver/OseenShapeDerivative.hpp>

#include <lifemc/lifesolver/BCInterface_Function.hpp>

namespace LifeV
{

//! BCInterface_OperatorFunction - LifeV bcFunction wrapper for BCInterface (with operators)
/*!
 *  @author Cristiano Malossi
 *
 *  This class is an interface between BCInterface, the grammar parser and a general
 *  LifeV solver (such as Oseen or FSIOperator). It allows to construct LifeV
 *  functions type for boundary conditions, using a functions string loaded from
 *  a GetPot file in which are present some "solver" parameters.
 *
 *  The class can be used in two ways:
 *
 *  1) hereditating it and implementing the template specialization of createAccessList() and updatePhysicalSolverVariables();
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
template< class PhysicalSolverType >
class BCInterface_OperatorFunction: public virtual BCInterface_Function< PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface_Data                                                      data_Type;
    typedef BCInterface_Function< physicalSolver_Type >                           function_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface_OperatorFunction();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface_OperatorFunction( const data_Type& data );

    //! Destructor
    virtual ~BCInterface_OperatorFunction() {}

    //@}


    //! @name Methods
    //@{

    //! Update operator variables
    void updatePhysicalSolverVariables() {}

    //@}


    //! @name Set Methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    virtual void setData( const data_Type& data );

    //! Set an operator
    /*!
     * @param physicalSolver operator
     */
    void setPhysicalSolver( const boost::shared_ptr< PhysicalSolverType >& physicalSolver ) { M_physicalSolver = physicalSolver; }

    //! Set variable function
    /*!
     * @param name name of the variable
     * @param value value of the variable
     */
    void setVariable( const std::string& name, const Real& value ) { function_Type::M_parser->setVariable( name, value ); }

    //@}

protected:

    //! @name Protected Methods
    //@{

    void createAccessList( const data_Type& data );

    //@}

    //List of all available operators
    enum physicalSolverList
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

    boost::shared_ptr< PhysicalSolverType >    M_physicalSolver;
    BCFlag                                     M_flag;
    std::set< physicalSolverList >             M_list;

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface_OperatorFunction( const BCInterface_OperatorFunction& function );

    BCInterface_OperatorFunction& operator=( const BCInterface_OperatorFunction& function );

    //@}


    //! @name Private Methods
    //@{

    void createFluidMap( std::map< std::string, physicalSolverList >& mapList );
    void createSolidMap( std::map< std::string, physicalSolverList >& mapList );
    void createList( const std::map< std::string, physicalSolverList >& mapList, const data_Type& data );

    void switchErrorMessage( const std::string& operatorType ) { std::cout << "ERROR: Invalid variable type for " << operatorType << "OperatorFunction" << std::endl; }

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterface_Function< PhysicalSolverType >* createBCInterface_OperatorFunction()
{
    return new BCInterface_OperatorFunction< PhysicalSolverType > ();
}

// ===================================================
// Constructors
// ===================================================
template< class PhysicalSolverType >
BCInterface_OperatorFunction< PhysicalSolverType >::BCInterface_OperatorFunction() :
        function_Type                    (),
        M_physicalSolver                 (),
        M_flag                           (),
        M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction::BCInterface_OperatorFunction()" << "\n";
#endif

}

template< class PhysicalSolverType >
BCInterface_OperatorFunction< PhysicalSolverType >::BCInterface_OperatorFunction( const data_Type& data ) :
        function_Type                    (),
        M_physicalSolver                 (),
        M_flag                           (),
        M_list                           ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction::BCInterface_OperatorFunction( data )" << "\n";
#endif

    this->setData( data );
}

// ===================================================
// Methods
// ===================================================
template< >
inline void
BCInterface_OperatorFunction< FSIOperator >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<FSIOperator>::UpdateOperatorVariables  " << "\n";
#endif

    // Create/Update variables for FSI problem
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
            // f_ -> FLUID
        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->fluid().area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                  f_density: " << M_physicalSolver->fluid().density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->fluid().density() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().flux( M_flag ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->fluid().flux( M_flag, *M_physicalSolver->un() ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->fluid().pressure( M_flag ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->fluid().pressure( M_flag, *M_physicalSolver->un() ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                f_viscosity: " << M_physicalSolver->fluid().viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->fluid().viscosity() );

            break;

            // s_ -> SOLID
        case s_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   s_density: " << M_physicalSolver->solid().rho() << "\n";
#endif

            setVariable( "s_density", M_physicalSolver->solid().rho() );

            break;

        case s_poisson:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   s_poisson: " << M_physicalSolver->solid().poisson() << "\n";
#endif

            setVariable( "s_poisson", M_physicalSolver->solid().poisson() );

            break;

        case s_thickness:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                 s_thickness: " << M_physicalSolver->solid().thickness() << "\n";
#endif

            setVariable( "s_thickness", M_physicalSolver->solid().thickness() );

            break;

        case s_young:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                     s_young: " << M_physicalSolver->solid().young() << "\n";
#endif

            setVariable( "s_young", M_physicalSolver->solid().young() );

            break;

        default:
            switchErrorMessage( "FSI" );
        }
}

template< >
inline void
BCInterface_OperatorFunction< Oseen< RegionMesh3D< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<Oseen>::UpdateOperatorVariables  " << "\n";
#endif

    // Create/Update variables for Oseen problem
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
            // f_ -> FLUID
        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                  f_density: " << M_physicalSolver->density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->density() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux( M_flag ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->flux( M_flag ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure( M_flag ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->pressure( M_flag ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                f_viscosity: " << M_physicalSolver->viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->viscosity() );

            break;

        default:

            switchErrorMessage( "OSEEN" );
        }
}

template< >
inline void
BCInterface_OperatorFunction< OseenShapeDerivative< RegionMesh3D< LinearTetra > > >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<OseenShapeDerivative>::UpdateOperatorVariables  " << "\n";
#endif

    // Create/Update variables for OseenShapeDerivative problem
    for ( std::set< physicalSolverList >::iterator j = M_list.begin(); j != M_list.end(); ++j )
        switch ( *j )
        {
            // f_ -> FLUID
        case f_area:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_area(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->area( M_flag ) << "\n";
#endif
            setVariable( "f_area", M_physicalSolver->area( M_flag ) );

            break;

        case f_density:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                f_density(): " << M_physicalSolver->density() << "\n";
#endif
            setVariable( "f_density", M_physicalSolver->density() );

            break;

        case f_flux:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                                   f_flux(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->flux( M_flag ) << "\n";
#endif

            setVariable( "f_flux", M_physicalSolver->flux( M_flag ) );

            break;

        case f_pressure:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                               f_pressure(" << static_cast<Real> (M_flag) << "): " << M_physicalSolver->pressure( M_flag ) << "\n";
#endif

            setVariable( "f_pressure", M_physicalSolver->pressure( M_flag ) );

            break;

        case f_viscosity:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5023 ) << "                                              f_viscosity(): " << M_physicalSolver->viscosity() << "\n";
#endif
            setVariable( "f_viscosity", M_physicalSolver->viscosity() );

            break;

        default:

            switchErrorMessage( "OSEENSHAPEDERIVATIVE" );
        }
}

// ===================================================
// Set Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface_OperatorFunction< PhysicalSolverType >::setData( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction::setData" << "\n";
#endif

    M_flag     = data.flag();

    function_Type::setData( data );

    createAccessList( data );
}

// ===================================================
// Protected Methods
// ===================================================
template< >
inline void
BCInterface_OperatorFunction< FSIOperator >::createAccessList( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<FSIOperator>::createAccessList" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createSolidMap( mapList );
    createList( mapList, data );

    //if ( M_physicalSolver.get() )
    //    updatePhysicalSolverVariables();
}

template< >
inline void
BCInterface_OperatorFunction< Oseen< RegionMesh3D< LinearTetra > > >::createAccessList( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<Oseen>::createAccessList" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createList( mapList, data );

    //if ( M_physicalSolver.get() )
    //    updatePhysicalSolverVariables();
}

template< >
inline void
BCInterface_OperatorFunction< OseenShapeDerivative< RegionMesh3D< LinearTetra > > >::createAccessList( const data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5023 ) << "BCInterface_OperatorFunction<OseenShapeDerivative>::createAccessList" << "\n";
#endif

    std::map< std::string, physicalSolverList > mapList;

    createFluidMap( mapList );
    createList( mapList, data );

    //if ( M_physicalSolver.get() )
    //    updatePhysicalSolverVariables();
}

// ===================================================
// Private Methods
// ===================================================
template< class PhysicalSolverType >
inline void
BCInterface_OperatorFunction< PhysicalSolverType >::createFluidMap( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["f_area"]      = f_area;
    mapList["f_density"]   = f_density;
    mapList["f_flux"]      = f_flux;
    mapList["f_pressure"]  = f_pressure;
    mapList["f_viscosity"] = f_viscosity;
}

template< class PhysicalSolverType >
inline void
BCInterface_OperatorFunction< PhysicalSolverType >::createSolidMap( std::map< std::string, physicalSolverList >& mapList )
{
    mapList["s_density"]   = s_density;
    mapList["s_poisson"]   = s_poisson;
    mapList["s_thickness"] = s_thickness;
    mapList["s_young"]     = s_young;
}

template< class PhysicalSolverType >
inline void
BCInterface_OperatorFunction< PhysicalSolverType >::createList( const std::map< std::string, physicalSolverList >& mapList, const data_Type& data )
{
    M_list.clear();
    for ( typename std::map< std::string, physicalSolverList >::const_iterator j = mapList.begin(); j != mapList.end(); ++j )
        if ( boost::find_first( data.baseString(), j->first ) )
            M_list.insert( j->second );
}

} // Namespace LifeV

#endif /* BCInterface_OperatorFunction_H */
