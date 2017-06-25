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

// BCInterface includes
#include <lifev/bc_interface/core/bc/BCInterface.hpp>
#include <lifev/bc_interface/3D/bc/BCInterfaceData3D.hpp>

/** When using BCInterface3D, you have to include beforehand the relevant functions.
 * In particular:
 *
 * for Fluid:
#include <lifev/bc_interface/3D/function/fluid/BCInterfaceFunctionParserFluid3D.hpp>
#include <lifev/bc_interface/3D/function/fluid/BCInterfaceFunctionParserSolverFluid3D.hpp>
#include <lifev/bc_interface/3D/function/fluid/BCInterfaceFunctionUserDefinedFluid3D.hpp>

 * for Solid:
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionParserSolid3D.hpp>
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionParserSolverSolid3D.hpp>
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionSolverDefinedSolid3D.hpp>
#include <lifev/bc_interface/3D/function/solid/BCInterfaceFunctionUserDefinedSolid3D.hpp>

 * For FSI:
#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionParserFSI3D.hpp>
#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionParserSolverFSI3D.hpp>
#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionSolverDefinedFSI3D.hpp>
#include <lifev/bc_interface/3D/function/fsi/BCInterfaceFunctionUserDefinedFSI3D.hpp>
*/

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
 *      <li> \c function, which is implemented in \c BCInterfaceFunctionParser;
 *      <li> \c functionFile, which is implemented in \c BCInterfaceFunctionParserFile;
 *      <li> \c functionSolver, which is implemented in \c BCInterfaceFunctionParserSolver;
 *      <li> \c functionFileSolver, which is implemented in \c BCInterfaceFunctionParserFileSolver;
 *      <li> \c functionUD, which is implemented in \c BCInterfaceFunctionUserDefined;
 *      <li> \c functionSD, which is implemented in \c BCInterfaceFunctionSolverDefined;
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

    typedef BCInterface< BcHandler, PhysicalSolverType >                bcInterface_Type;

    typedef typename bcInterface_Type::bcHandler_Type                   bcHandler_Type;
    typedef typename bcInterface_Type::bcHandlerPtr_Type                bcHandlerPtr_Type;

    typedef typename bcInterface_Type::physicalSolver_Type              physicalSolver_Type;
    typedef typename bcInterface_Type::physicalSolverPtr_Type           physicalSolverPtr_Type;

    typedef typename bcInterface_Type::factory_Type                     factory_Type;

    typedef typename bcInterface_Type::bcFunctionPtr_Type               bcFunctionPtr_Type;
    typedef typename bcInterface_Type::vectorFunction_Type              vectorFunction_Type;

    typedef typename bcInterface_Type::bcFunctionSolverDefinedPtr_Type  bcFunctionSolverDefinedPtr_Type;
    typedef typename bcInterface_Type::vectorFunctionSolverDefined_Type vectorFunctionSolverDefined_Type;

    typedef BCInterfaceData3D                                           data_Type;
    typedef std::shared_ptr< data_Type >                              dataPtr_Type;

    typedef BCFunctionRobin                                             bcFunctionRobin_Type;
    typedef std::shared_ptr< bcFunctionRobin_Type >                   bcFunctionRobinPtr_Type;
    typedef std::vector< bcFunctionRobinPtr_Type >                      vectorFunctionRobin_Type;

    typedef BCFunctionDirectional                                       bcFunctionDirectional_Type;
    typedef std::shared_ptr< bcFunctionDirectional_Type >             bcFunctionDirectionalPtr_Type;
    typedef std::vector< bcFunctionDirectionalPtr_Type >                vectorFunctionDirectional_Type;

    typedef BCDataInterpolator                                          bcFunctionDataInterpolator_Type;
    typedef std::shared_ptr< bcFunctionDataInterpolator_Type >        bcFunctionDataInterpolatorPtr_Type;
    typedef std::vector< bcFunctionDataInterpolatorPtr_Type >           vectorDataInterpolator_Type;

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

    //! Read a specific boundary condition from a file and add it to the data container
    /*!
     * @param fileName Name of the data file
     * @param dataSection section in the data file
     * @param name name of the boundary condition
     */
    void readBC ( const std::string& fileName, const std::string& dataSection, const std::string& name )
    {
        M_data->readBC ( fileName, dataSection, name );
    }

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
    void addBC ( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base )
    {
        this->M_handler->addBC ( name, flag, type, mode, base );
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
    void addBC ( const bcName_Type& name, const bcFlag_Type& flag, const bcType_Type& type, const bcMode_Type& mode, BCBaseType& base, const BCCompType& comp )
    {
        this->M_handler->addBC ( name, flag, type, mode, base, comp );
    }

    //@}


    //! @name Set Methods
    //@{

    //! Set a physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver );

    //@}


    //! @name Get Methods
    //@{

    //! Get the data container
    /*!
     * @return the data container
     */
    data_Type& dataContainer()
    {
        return *M_data;
    }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3D ( const BCInterface3D& bcInterface3D );

    BCInterface3D& operator= ( const BCInterface3D& bcInterface3D );

    //@}


    //! @name Private Methods
    //@{

    template< class BCBaseType >
    void createFunctionRobin ( BCBaseType& base );

    template< class BCBaseType >
    void createFunctionDirectional ( BCBaseType& base );

    void createFunctionDataInterpolator();

    template< class BCBaseType >
    void addBcToHandler ( BCBaseType& base );

    //@}

    // Data
    dataPtr_Type                    M_data;

    // Functions Robin
    vectorFunctionRobin_Type        M_vectorFunctionRobin;

    // Functions Directions
    vectorFunctionDirectional_Type  M_vectorFunctionDirection;

    // Data Interpolator Functions
    vectorDataInterpolator_Type     M_vectorDataInterpolator;

};

// ===================================================
// Constructors & Destructor
// ===================================================
template< class BcHandler, class PhysicalSolverType >
BCInterface3D< BcHandler, PhysicalSolverType >::BCInterface3D() :
    bcInterface_Type          (),
    M_data                    ( new data_Type() ),
    M_vectorFunctionRobin     (),
    M_vectorFunctionDirection (),
    M_vectorDataInterpolator  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 5020 ) << "BCInterface3D::BCInterface3D" << "\n";
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
    debugStream ( 5020 ) << "BCInterface3D::insertBC\n";
#endif

    switch ( M_data->base().second )
    {
        case BCIFunctionParser:
        case BCIFunctionParserFile:
        case BCIFunctionParserSolver:
        case BCIFunctionParserFileSolver:
        case BCIFunctionUserDefined:
        {
            factory_Type factory;
            this->M_vectorFunction.push_back ( factory.createFunctionParser ( M_data ) );

            BCFunctionBase base;
            this->M_vectorFunction.back()->assignFunction ( base );

            // Directional BC
            if ( M_data->mode() == Directional )
            {
                createFunctionDirectional ( base );
                addBcToHandler ( *M_vectorFunctionDirection.back() );

                return;
            }

            // Robin BC
            if ( M_data->type() == Robin )
            {
                createFunctionRobin ( base );
                addBcToHandler ( *M_vectorFunctionRobin.back() );

                return;
            }

            // All the other type of BC
            addBcToHandler ( base );

            return;
        }
        case BCIFunctionSolverDefined:
        {
            factory_Type factory;
            this->M_vectorFunctionSolverDefined.push_back ( factory.createFunctionSolverDefined ( M_data ) );

            return;
        }
        case BCI3DDataInterpolator:

            createFunctionDataInterpolator();
            addBcToHandler ( *M_vectorDataInterpolator.back() );

            return;

        default:

            std::cout << " !!! Error: " << M_data->base().first << " is not valid in BCInterface3D !!!" << std::endl;

            break;
    }
}

// ===================================================
// Set Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >
void
BCInterface3D< BcHandler, PhysicalSolverType >::setPhysicalSolver ( const physicalSolverPtr_Type& physicalSolver )
{
    bcInterface_Type::setPhysicalSolver ( physicalSolver );

    for ( typename vectorFunctionSolverDefined_Type::const_iterator i = this->M_vectorFunctionSolverDefined.begin() ; i < this->M_vectorFunctionSolverDefined.end() ; ++i )
    {
        ( *i )->exportData ( M_data );

        // Robin BC
        if ( M_data->type() == Robin )
        {
            BCVector base;

            ( *i )->assignFunction ( base );
            addBcToHandler ( base );
        }
        else
        {
            BCVectorInterface base;

            ( *i )->assignFunction ( base );
            addBcToHandler ( base );
        }
    }
}

// ===================================================
// Private Methods
// ===================================================
template< class BcHandler, class PhysicalSolverType >  template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionRobin ( BCBaseType& base )
{
    // Parameters for direction BC
    M_data->setName ( M_data->name() + "_robinMassTerm" );
    M_data->setRobinBaseAlpha();

    // Create the mass term function
    factory_Type factory;
    this->M_vectorFunction.push_back ( factory.createFunctionParser ( M_data ) );

    BCFunctionBase baseRobin;
    this->M_vectorFunction.back()->assignFunction ( baseRobin );

    // Robin base
    bcFunctionRobinPtr_Type robinBase ( new bcFunctionRobin_Type ( base.Function(), baseRobin.Function() ) );
    M_vectorFunctionRobin.push_back ( robinBase );
}

template< class BcHandler, class PhysicalSolverType >  template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionDirectional ( BCBaseType& base )
{
    // Parameters for direction BC
    M_data->setName ( M_data->name() + "_directionalField" );
    M_data->setDirectionalBase();

    // Create the directional field
    factory_Type factory;
    this->M_vectorFunction.push_back ( factory.createFunctionParser ( M_data ) );

    BCFunctionBase baseDirectional;
    this->M_vectorFunction.back()->assignFunction ( baseDirectional );

    // Directional base
    bcFunctionDirectionalPtr_Type directionalBase ( new bcFunctionDirectional_Type ( base.Function(), baseDirectional.Function() ) );
    M_vectorFunctionDirection.push_back ( directionalBase );
}

template< class BcHandler, class PhysicalSolverType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::createFunctionDataInterpolator()
{
    // Directional base
    bcFunctionDataInterpolatorPtr_Type dataInterpolatorBase ( new bcFunctionDataInterpolator_Type() );
    dataInterpolatorBase->readData ( M_data->baseString() );
    dataInterpolatorBase->setInterpolationMethod ( LifeV::BCDataInterpolator::RBF_InverseMultiQuadric);
    M_vectorDataInterpolator.push_back ( dataInterpolatorBase );
}

template< class BcHandler, class PhysicalSolverType > template< class BCBaseType >
inline void
BCInterface3D< BcHandler, PhysicalSolverType >::addBcToHandler ( BCBaseType& base )
{
    if ( !this->M_handler.get() ) // If BCHandler has not been created yet, we do it now
    {
        this->createHandler();
    }

    switch ( M_data->mode() )
    {
        case Scalar:
        case Normal:
        case Tangential:
        case Directional:

#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 5020 ) << "BCInterface3D::addBcToHandler                            Scalar, Normal, Tangential, Directional" << "\n\n";
#endif

            this->M_handler->addBC ( M_data->name(), M_data->flag(), M_data->type(), M_data->mode(), base );

            break;

        case Full:

#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 5020 ) << "BCInterface3D::addBcToHandler                            Full" << "\n\n";
#endif

            this->M_handler->addBC ( M_data->name(), M_data->flag(), M_data->type(), M_data->mode(), base, M_data->componentsNumber() );

            break;

        case Component:

#ifdef HAVE_LIFEV_DEBUG
            debugStream ( 5020 ) << "BCInterface3D::addBcToHandler                            Component" << "\n\n";
#endif

            this->M_handler->addBC ( M_data->name(), M_data->flag(), M_data->type(), M_data->mode(), base, M_data->componentsVector() );

            break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3D_H */
