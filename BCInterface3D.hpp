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

#include <lifemc/lifesolver/BCInterface.hpp>

#include <lifemc/lifesolver/BCInterface3DFSI.hpp>

namespace LifeV
{

//! BCInterface3D - LifeV Interface to load Boundary Conditions completely from a GetPot file
/*!
 *  @author Cristiano Malossi
 *
 *  This class allows to impose boundary conditions completely from a file.
 *
 *  <b>EXAMPLE - DATA FILE</b>
 *
 *  In the GetPot data file, inside each subsection ([fluid], [solid], etc...) BCInterface3D reads a new section: [boundary_conditions].
 *
 *  Inside the new section there is a list of condition which correspond to other sub-section
 *  with the same name. The list must be inside the apex ' '.
 *
 *  Each condition has a similar structure; here there is an example:
 *
 *  [InFlow]               <br>
 *  type       = Essential <br>
 *  flag       = 2         <br>
 *  mode       = Full      <br>
 *  component  = 3         <br>
 *  function   = '[0, 0, 3*0.03*(1/4-(x^2+y^2)]' <br>
 *
 *  NOTE: All the parameters are case sensitive.
 *
 *  type - can be: Essential Natural Robin Flux
 *  flag - contains the flag
 *  mode - can be: Full Component Scalar Tangential Normal.
 *  component - if mode is Scalar, Tangential or Normal it is missing.
 *              if mode is Component it contains the ID of the component (or of the components list inside apex)
 *              if mode is Full it contains the total number of components
 *  function - contains the function. See BCInterfaceFunction for more details about the syntax.
 *
 *  <b>NOTE:</b>
 *
 *  The string "function" represent the base module and can be replaced by other expanded modules.
 *  Up to now we have the following modules for function:
 *
 *  - function
 *  - functionFile
 *  - OperatorFunction
 *  - OperatorFunctionFile
 *  - FSI
 *
 *  To see some example look at test_fsi.
 *
 *  <b>EXAMPLE - HOW TO USE</b>
 *
 *  Here there is a short example on how to use it.
 *
 *  1) You can define your BCInterface class in a shared pointer:
 *     boost::shared_ptr<BCInterface> 	M_fluidBC;
 *
 *  2) Build the BCInterface using empty constructor;
 *
 *  3) If you have operator conditions you have to give the operator to access variables
 *     M_fluidBC->setPhysicalSolver( M_fsi->FSIOper() );
 *
 *  4) Then you can fill the handler from a file and a section (this can be done for multiple files & sections)
 *     M_fluidBC->fillHandler( "fileName.dat", "fluid" );
 *
 *  5) Finally, to get the handler you can use:
 *     M_fluidBC->handler();
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

    typedef BCFunctionDirectional                                 bcFunctionDirectional_Type;
    typedef boost::shared_ptr< bcFunctionDirectional_Type >       bcFunctionDirectionalPtr_Type;
    typedef std::vector< bcFunctionDirectionalPtr_Type >          vectorFunctionDirectional_Type;

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

    template< class BCInterfaceBaseType >
    void createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector );

    // This method should be removed: it is a workaround due to legacy of LifeV BC.
    void addBcToHandler( BCVectorInterface& base );

    template< class BCBaseType >
    void addBcToHandler( BCBaseType& base );

    //@}

    // Functions Directions
    vectorFunctionDirectional_Type  M_vectorFunctionDirection;

    // FSI Functions
    vectorFSI_Type                  M_vectorFSI;
};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class BcHandler, class PhysicalSolverType >
BCInterface3D< BcHandler, PhysicalSolverType >::BCInterface3D() :
        M_vectorFunctionDirection (),
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

        addBcToHandler( base );

        break;
    }
    case BCI3DFSI:

        createFunction( M_vectorFSI );

        break;

    default:

        std::cout << " !!! Error: " << this->M_data.base().first << " is not valid in BCInterface3D !!!" << std::endl;
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface3D< BcHandler, PhysicalSolverType >::setPhysicalSolver( const boost::shared_ptr< physicalSolver_Type >& physicalSolver )
{
    //for ( typename vectorFunction_Type::const_iterator i = this->M_vectorFunction.begin() ; i < this->M_vectorFunction.end() ; ++i )
    for ( UInt i( 0 ); i < this->M_vectorFunction.size(); ++i )
    {
        BCInterfaceFunctionSolver< physicalSolver_Type > *castedOperator =
            dynamic_cast < BCInterfaceFunctionSolver< physicalSolver_Type >* > ( &( *this->M_vectorFunction[i] ) );

        if ( castedOperator != 0 )
            castedOperator->setPhysicalSolver( physicalSolver );
    }

    for ( typename vectorFSI_Type::const_iterator i = M_vectorFSI.begin() ; i < M_vectorFSI.end() ; ++i )
    {
        BCVectorInterface base;

        ( *i )->assignFunction( physicalSolver, base );
        ( *i )->exportData( this->M_data );

        addBcToHandler( base );
    }
}


// ===================================================
// Private Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType > template< class BCInterfaceBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunction( std::vector< boost::shared_ptr< BCInterfaceBaseType > >& baseVector )
{
    boost::shared_ptr< BCInterfaceBaseType > function( new BCInterfaceBaseType( this->M_data ) );
    baseVector.push_back( function );
}

template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::addBcToHandler( BCVectorInterface& base )
{
    if ( !this->M_handler.get() ) // If BCHandler has not been created yet, we do it now
        this->createHandler();

    switch ( this->M_data.mode() )
    {
    case Scalar:
    case Normal:
    case Tangential:
    case Directional:
    case Component:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              others" << "\n\n";
#endif

        std::cout << "ERROR: Scalar, Normal, Tangential, Directional, Component NOT AVAILABLE FOR FSI BC" << std::endl;
        break;

    case Full:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              Full" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base, this->M_data.comN() );

        break;
    }
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

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              Scalar, Normal, Tangential" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base );

        break;

    case Directional:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              Directional" << "\n\n";
#endif
        {
            // Parameters for direction BC
            this->M_data.setName( this->M_data.name() + "_direction" );
            this->M_data.setBase( make_pair( "function", BCIFunction ) );
            this->M_data.setBaseString( this->M_data.direction() );

            // Directional field
            factory_Type factory;
            this->M_vectorFunction.push_back( factory.createFunction( this->M_data ) );

            BCFunctionBase baseDirectional;
            this->M_vectorFunction.back()->assignFunction( baseDirectional );

            // Directional base
            boost::shared_ptr< BCFunctionDirectional > directionalBase( new BCFunctionDirectional( base.Function(), baseDirectional.Function() ) );
            M_vectorFunctionDirection.push_back( directionalBase );

            this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), *M_vectorFunctionDirection.back() );
        }

        break;

    case Full:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              Full" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base, this->M_data.comN() );

        break;

    case Component:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5020 ) << "BCInterface3D::addBCManager                              Component" << "\n\n";
#endif

        this->M_handler->addBC( this->M_data.name(), this->M_data.flag(), this->M_data.type(), this->M_data.mode(), base, this->M_data.comV() );

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3D_H */
