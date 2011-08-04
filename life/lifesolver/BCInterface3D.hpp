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
 *  @brief File containing the BCInterface3D class
 *
 *  @date 01-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3D_H
#define BCInterface3D_H 1

#include <life/lifefem/BCDataInterpolator.hpp>

#include <life/lifesolver/BCInterface.hpp>

#include <life/lifesolver/BCInterface3DFSI.hpp>

namespace LifeV
{

//! BCInterface3D - LifeV interface to load boundary conditions for 3D problems completely from a \c GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions for a 3D problem completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b> <BR>
 *  In the GetPot data file, \c BCInterface reads a new section: <CODE> [boundary_conditions] </CODE>.
 *
 *  Inside the new section there is a list of boundary conditions which correspond to other sub-section
 *  with the same name, for example: <CODE> list = 'InFlow OutFlow' </CODE>
 *
 *  Each boundary condition has a similar structure. The list of properties depends from the type of the
 *  boundary condition. For example:
 *
 *  <CODE>
 *  [InFlow]                             <BR>
 *  type       = Essential               <BR>
 *  flag       = 2                       <BR>
 *  mode       = Full                    <BR>
 *  component  = 3                       <BR>
 *  function   = '[0, 0, 3*0.03*(1/4-(x^2+y^2)]' <BR>
 *
 *  [OutFlow]        <BR>
 *  type       = Essential               <BR>
 *  flag       = 3                       <BR>
 *  mode       = Full                    <BR>
 *  component  = 3                       <BR>
 *  function   = '0'                     <BR>
 *  </CODE>
 *
 *  where \c type, \c flag, \c mode, \c component are the classical parameters for a 3D boundary condition.
 *  The string \c function represents the base module and can be replaced by other derived/alternative modules.
 *  The following functions are available (see the related classes for more information):
 *
 *  <ol>
 *      <li> \c function, which is implemented in \c BCInterfaceFunction;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionFileSolver;
 *      <li> \c FSI, which is implemented in \c BCInterface3DFSI;
 *      <li> \c dataInterpolator, which is implemented in\c BCDataInterpolator;
 *  </ol>
 *
 *  All the parameters are case sensitive.
 *
 *  See \c BCInterface base class for more details.
 */
template< class BcHandler, class PhysicalSolverType >
class BCInterface3D : public virtual BCInterface< BcHandler, PhysicalSolverType >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface< BcHandler, PhysicalSolverType >          bcInterface_Type;

    typedef typename bcInterface_Type::bcHandler_Type             bcHandler_Type;
    typedef typename bcInterface_Type::bcHandlerPtr_Type          bcHandlerPtr_Type;

    typedef typename bcInterface_Type::physicalSolver_Type        physicalSolver_Type;
    typedef typename bcInterface_Type::physicalSolverPtr_Type     physicalSolverPtr_Type;

    typedef typename bcInterface_Type::factory_Type               factory_Type;
    typedef typename bcInterface_Type::bcFunctionPtr_Type         bcFunctionPtr_Type;

    typedef typename bcInterface_Type::data_Type                  data_Type;

    typedef typename bcInterface_Type::vectorFunction_Type        vectorFunction_Type;

    typedef BCFunctionRobin                                       bcFunctionRobin_Type;
    typedef boost::shared_ptr< bcFunctionRobin_Type >             bcFunctionRobinPtr_Type;
    typedef std::vector< bcFunctionRobinPtr_Type >                vectorFunctionRobin_Type;

    typedef BCFunctionDirectional                                 bcFunctionDirectional_Type;
    typedef boost::shared_ptr< bcFunctionDirectional_Type >       bcFunctionDirectionalPtr_Type;
    typedef std::vector< bcFunctionDirectionalPtr_Type >          vectorFunctionDirectional_Type;

    typedef BCDataInterpolator                                    bcFunctionDataInterpolator_Type;
    typedef boost::shared_ptr< bcFunctionDataInterpolator_Type >  bcFunctionDataInterpolatorPtr_Type;
    typedef std::vector< bcFunctionDataInterpolatorPtr_Type >     vectorDataInterpolator_Type;

    typedef BCInterface3DFSI< physicalSolver_Type >               bcFunctionFSI_Type;
    typedef boost::shared_ptr< bcFunctionFSI_Type >               bcFunctionFSIPtr_Type;
    typedef std::vector< bcFunctionFSIPtr_Type >                  vectorFSI_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3D();

    //! Destructor
    virtual ~BCInterface3D() {}

    //@}


    //! @name Methods
    //@{

    //! Insert the current boundary condition in the BChandler
    void insertBC();

    //! Update the variables inside the physical solver
    void updatePhysicalSolverVariables();

    //@}


    //! @name External interface for BCHandler functions
    //@{

    //! Add a Boundary Condition using the standard interface of the BCHandler
    /*!
     * @param name name of the condition
     * @param flag list of flags
     * @param type type of the condition
     * @param mode mode of the condition
     * @param base base of the condition
     */
    template< class BCBaseType >
    void addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base )
    {
        this->M_handler->addBC( name, flag, type, mode, base );
    }

    //! Add a Boundary Condition with component using the standard interface of the BCHandler
    /*!
     * @param name name of the condition
     * @param flag list of flags
     * @param type type of the condition
     * @param mode mode of the condition
     * @param base base of the condition
     * @param comp component of the condition
     */
    template< class BCBaseType, class BCCompType >
    void addBC( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base, const BCCompType& comp )
    {
        this->M_handler->addBC( name, flag, type, mode, base, comp );
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set a physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3D( const BCInterface3D& bcInterface3D );

    BCInterface3D& operator=( const BCInterface3D& bcInterface3D );

    //@}


    //! @name Private Methods
    //@{

    template< class BCBaseType >
    void createFunctionRobin( BCBaseType& base );

    template< class BCBaseType >
    void createFunctionDirectional( BCBaseType& base );

    void createFunctionDataInterpolator();

    void createFunctionFSI();

    template< class BCBaseType >
    void addBcToHandler( BCBaseType& base );

    //@}

    // Functions Robin
    vectorFunctionRobin_Type        M_vectorFunctionRobin;

    // Functions Directions
    vectorFunctionDirectional_Type  M_vectorFunctionDirection;

    // Data Interpolator Functions
    vectorDataInterpolator_Type     M_vectorDataInterpolator;

    // FSI Functions
    vectorFSI_Type                  M_vectorFSI;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class BcHandler, class PhysicalSolverType >
BCInterface3D< BcHandler, PhysicalSolverType >::BCInterface3D() :
        bcInterface_Type          (),
        M_vectorFunctionRobin     (),
        M_vectorFunctionDirection (),
        M_vectorDataInterpolator  (),
        M_vectorFSI               ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::BCInterface3D" << "\n";
#endif

}

// ===================================================
// Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface3D< BcHandler, PhysicalSolverType >::insertBC()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::insertBC\n";
#endif

    switch ( this->M_data.base().second )
    {
    case BCIFunction:
    case BCIFunctionFile:
    case BCIFunctionSolver:
    case BCIFunctionFileSolver:
    {
        factory_Type factory;
        this->M_vectorFunction.push_back( factory.createFunction( this->M_data ) );

        BCFunctionBase base;
        this->M_vectorFunction.back()->assignFunction( base );

        // Directional BC
        if ( this->M_data.mode() == Directional )
        {
            createFunctionDirectional( base );
            addBcToHandler( *M_vectorFunctionDirection.back() );

            return;
        }

        // Robin BC
        if ( this->M_data.type() == Robin )
        {
            createFunctionRobin( base );
            addBcToHandler( *M_vectorFunctionRobin.back() );

            return;
        }

        // All the other type of BC
        addBcToHandler( base );

        return;
    }
    case BCI3DDataInterpolator:

        createFunctionDataInterpolator();
        addBcToHandler( *M_vectorDataInterpolator.back() );

        return;

    case BCI3DFSI:

        createFunctionFSI();

        return;

    default:

        std::cout << " !!! Error: " << this->M_data.base().first << " is not valid in BCInterface3D !!!" << std::endl;

        break;
    }
}

template< class BcHandler, class PhysicalSolverType >
void
BCInterface3D< BcHandler, PhysicalSolverType >::updatePhysicalSolverVariables()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5020 ) << "BCInterface3D::updatePhysicalSolverVariables\n";
#endif

    bcInterface_Type::updatePhysicalSolverVariables();

    for ( typename vectorFSI_Type::const_iterator i = M_vectorFSI.begin() ; i < M_vectorFSI.end() ; ++i )
        ( *i )->updatePhysicalSolverVariables();
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface3D< BcHandler, PhysicalSolverType >::setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver )
{
    bcInterface_Type::setPhysicalSolver( physicalSolver );

    for ( typename vectorFSI_Type::const_iterator i = M_vectorFSI.begin() ; i < M_vectorFSI.end() ; ++i )
    {
        ( *i )->exportData( this->M_data );

        // Robin BC
        if ( this->M_data.type() == Robin )
        {
            BCVector base;

            ( *i )->assignFunction( physicalSolver, base );
            addBcToHandler( base );
        }
        else
        {
            BCVectorInterface base;

            ( *i )->assignFunction( physicalSolver, base );
            addBcToHandler( base );
        }
    }
}

// ===================================================
// Private Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >  template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionRobin( BCBaseType& base )
{
    // Parameters for direction BC
    this->M_data.setName( this->M_data.name() + "_robinMassTerm" );
    this->M_data.setRobinBaseAlpha();

    // Create the mass term function
    factory_Type factory;
    this->M_vectorFunction.push_back( factory.createFunction( this->M_data ) );

    BCFunctionBase baseRobin;
    this->M_vectorFunction.back()->assignFunction( baseRobin );

    // Robin base
    bcFunctionRobinPtr_Type robinBase( new bcFunctionRobin_Type( base.Function(), baseRobin.Function() ) );
    M_vectorFunctionRobin.push_back( robinBase );
}

template< class BcHandler, class PhysicalSolverType >  template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionDirectional( BCBaseType& base )
{
    // Parameters for direction BC
    this->M_data.setName( this->M_data.name() + "_directionalField" );
    this->M_data.setDirectionalBase();

    // Create the directional field
    factory_Type factory;
    this->M_vectorFunction.push_back( factory.createFunction( this->M_data ) );

    BCFunctionBase baseDirectional;
    this->M_vectorFunction.back()->assignFunction( baseDirectional );

    // Directional base
    bcFunctionDirectionalPtr_Type directionalBase( new bcFunctionDirectional_Type( base.Function(), baseDirectional.Function() ) );
    M_vectorFunctionDirection.push_back( directionalBase );
}

template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionDataInterpolator()
{
    // Directional base
    bcFunctionDataInterpolatorPtr_Type dataInterpolatorBase( new bcFunctionDataInterpolator_Type() );
    dataInterpolatorBase->readData( this->M_data.baseString() );
    M_vectorDataInterpolator.push_back( dataInterpolatorBase );
}

template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionFSI()
{
    bcFunctionFSIPtr_Type function( new bcFunctionFSI_Type( this->M_data ) );
    M_vectorFSI.push_back( function );
}

template< class BcHandler, class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::addBcToHandler( BCBaseType& base )
{
    if ( !this->M_handler.get() ) // If BCHandler has not been created yet, we do it now
        this->createHandler();

    switch ( this->M_data.mode() )
    {
    case Scalar:
    case Normal:
    case Tangential:
    case Directional:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBcToHandler                            Scalar, Normal, Tangential, Directional" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base );

        break;

    case Full:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBcToHandler                            Full" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base, this->M_data.comN() );

        break;

    case Component:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBcToHandler                            Component" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base, this->M_data.comV() );

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3D_H */
