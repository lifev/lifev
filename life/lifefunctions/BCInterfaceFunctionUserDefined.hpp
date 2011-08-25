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
 *  @brief File containing the BCInterfaceFunctionSolverDefined class
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterfaceFunctionSolverDefined_H
#define BCInterfaceFunctionSolverDefined_H 1

#include <life/lifesolver/FSIExactJacobian.hpp>
#include <life/lifesolver/FSIFixedPoint.hpp>

#include <life/lifesolver/FSIMonolithicGE.hpp>
#include <life/lifesolver/FSIMonolithicGI.hpp>

#include <life/lifesolver/BCInterfaceDefinitions.hpp>
#include <life/lifesolver/BCInterfaceData3D.hpp>
#include <life/lifefunctions/BCInterfaceFactory.hpp>

namespace LifeV
{

//! BCInterfaceFunctionSolverDefined - Empty class to be used for solver user defined specializations.
/*!
 *  @author Cristiano Malossi
 *
 *  This class provide the interfaces for the user defined functions.
 */
template< class PhysicalSolverType >
class BCInterfaceFunctionSolverDefined
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                    physicalSolver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit BCInterfaceFunctionSolverDefined() {}

    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

#ifdef MULTISCALE_IS_IN_LIFEV
    //! Copy the stored parameters in the 0D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData( BCInterfaceData0D& /*data*/ ) {}
#endif

    //! Copy the stored parameters in the 1D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData( BCInterfaceData1D& /*data*/ ) {}

    //! Copy the stored parameters in the 3D data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData( BCInterfaceData3D& /*data*/ ) {}

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition base
     */
    template< class BCBaseType >
    void assignFunction( BCBaseType& /*base*/ ) {}

    //! Update the solver variables
    void updatePhysicalSolverVariables() {}


    //@}


    //! @name Set Methods
    //@{

#ifdef MULTISCALE_IS_IN_LIFEV
    //! Set data for 0D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData( const BCInterfaceData0D& /*data*/ ) {}
#endif

    //! Set data for 1D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData( const BCInterfaceData1D& /*data*/ ) {}

    //! Set data for 3D boundary conditions
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData( const BCInterfaceData3D& /*data*/ ) {}

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver( const boost::shared_ptr< PhysicalSolverType >& /*physicalSolver*/ ) {}

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined( const BCInterfaceFunctionSolverDefined& function);

    BCInterfaceFunctionSolverDefined& operator=( const BCInterfaceFunctionSolverDefined& function );

    //@}

};

// ===================================================
// Factory
// ===================================================
//! Factory create function
template< typename PhysicalSolverType >
inline BCInterfaceFunctionSolverDefined< PhysicalSolverType >* createBCInterfaceFunctionSolverDefined()
{
    return new BCInterfaceFunctionSolverDefined< PhysicalSolverType > ();
}



//! BCInterfaceFunctionSolverDefined - Template specialization of BCInterfaceFunctionSolverDefined for 3D FSI problems
/*!
 *  @author Cristiano Malossi
 *
 *
 *  The BCInterfaceFunctionSolverDefined class provides a general interface between the
 *  \c BCInterface3D and the default boundary condition for the \c FSIOperator.
 *
 *  <b>DETAILS:</b> <BR>
 *  The constructor of the class takes a string contains the ID of the interface condition to impose,
 *  and the FSI. The list of available conditions is the FSIFunction variable. These are:
 *
 *  <ol>
 *      <li> DerFluidLoadToFluid,                        (not implemented)
 *        <li> DerFluidLoadToStructure,
 *        <li> DerHarmonicExtensionVelToFluid,
 *        <li> DerStructureDispToSolid,                  (not implemented)
 *        <li> FluidInterfaceDisp,                       (not working)
 *        <li> FluidLoadToStructure,
 *        <li> HarmonicExtensionVelToFluid,
 *        <li> SolidLoadToStructure,
 *        <li> StructureDispToHarmonicExtension,
 *        <li> StructureDispToSolid,                     (not implemented)
 *        <li> StructureToFluid
 *        <li> RobinWall
 *  </ol>
 *
 *  The class automatically recognize which FSI algorithm is used among:
 *  <ol>
 *      <li> EXACTJACOBIAN;
 *      <li> FIXEDPOINT;
 *      <li> MONOLITHIC_GE;
 *      <li> MONOLITHIC_GI;
 *  </ol>
 */
template< >
class BCInterfaceFunctionSolverDefined< FSIOperator >
{
public:

    //! @name Type definitions
    //@{

    typedef FSIOperator                                           physicalSolver_Type;
    typedef boost::shared_ptr< physicalSolver_Type >              physicalSolverPtr_Type;

    typedef BCInterfaceFactory< FSIOperator >                     factory_Type;

    typedef BCInterfaceFunctionParser< physicalSolver_Type >      bcFunctionParser_Type;
    typedef boost::shared_ptr< bcFunctionParser_Type >            bcFunctionParserPtr_Type;
    typedef std::vector< bcFunctionParserPtr_Type >               vectorFunctionParser_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterfaceFunctionSolverDefined();

    //! Destructor
    virtual ~BCInterfaceFunctionSolverDefined() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData( BCInterfaceData3D& data );

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition base
     */
    template< class BCBaseType >
    void assignFunction( BCBaseType& base );

    //! Update the solver variables
    void updatePhysicalSolverVariables();

    //@}


    //! @name Set methods
    //@{

    //! Set data
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData( const BCInterfaceData3D& data );

    //! Set the physical solver
    /*!
     * @param physicalSolver physical solver
     */
    void setPhysicalSolver( const physicalSolverPtr_Type& physicalSolver ) { M_physicalSolver = physicalSolver; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterfaceFunctionSolverDefined( const BCInterfaceFunctionSolverDefined& fsi);

    BCInterfaceFunctionSolverDefined& operator=( const BCInterfaceFunctionSolverDefined& fsi );

    //@}


    //! @name Private Methods
    //@{

    template< class MethodType >
    void checkFunction( BCVectorInterface& base );

    template< class MethodType >
    void checkFunction( BCVector& base );

    //@}

    enum FSIMethod
    {
        EXACTJACOBIAN, FIXEDPOINT, MONOLITHIC_GE, MONOLITHIC_GI
    };

    enum FSIFunction
    {
        DerFluidLoadToFluid,
        DerFluidLoadToStructure,
        DerHarmonicExtensionVelToFluid,
        DerStructureDispToSolid,
        FluidInterfaceDisp,
        FluidLoadToStructure,
        HarmonicExtensionVelToFluid,
        SolidLoadToStructure,
        StructureDispToHarmonicExtension,
        StructureDispToSolid,
        StructureToFluid,
        RobinWall
    };

    FSIFunction                                    M_FSIFunction;

    physicalSolverPtr_Type                         M_physicalSolver;

    // The following members are required since the FSI BC are applied
    // a posteriori, when setPhysicalSolver() is called.

    // Classical parameters
    bcName_Type                                    M_name;
    bcFlag_Type                                    M_flag;
    bcType_Type                                    M_type;
    bcMode_Type                                    M_mode;
    bcComponentsVec_Type                           M_componentsVector;

    // RobinViscoelastic
    vectorFunctionParser_Type                      M_vectorFunctionRobin;
    FSIOperator::vectorPtr_Type                    M_robinRHS;
    FSIOperator::vectorPtr_Type                    M_robinAlphaCoefficient;
    FSIOperator::vectorPtr_Type                    M_robinBetaCoefficient;
};

// ===================================================
// Methods
// ===================================================
template< class BCBaseType >
inline void
BCInterfaceFunctionSolverDefined< FSIOperator >::assignFunction( BCBaseType& base )
{
    //Set mapMethod
    std::map< std::string, FSIMethod > mapMethod;

    mapMethod["exactJacobian"] = EXACTJACOBIAN;
    mapMethod["fixedPoint"]    = FIXEDPOINT;
    mapMethod["monolithicGE"]  = MONOLITHIC_GE;
    mapMethod["monolithicGI"]  = MONOLITHIC_GI;

    switch ( mapMethod[M_physicalSolver->data().method()] )
    {
    case EXACTJACOBIAN:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkMethod                            exactJacobian" << "\n";
#endif

        checkFunction< FSIExactJacobian > ( base );

        break;

    case FIXEDPOINT:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkMethod                            fixedPoint" << "\n";
#endif

        checkFunction< FSIFixedPoint > ( base );

        break;

    case MONOLITHIC_GE:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkMethod                            monolithicGE" << "\n";
#endif

        checkFunction< FSIMonolithicGE >( base );

        break;

    case MONOLITHIC_GI:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkMethod                            monolithicGI" << "\n";
#endif

        checkFunction< FSIMonolithicGI >( base );

        break;

    default:

        std::cout << " !!! Warning:" << mapMethod[M_physicalSolver->data().method()] << " not assigned !!!" << std::endl;

        break;

    }
}

// ===================================================
// Private functions
// ===================================================
template< class MethodType >
inline void BCInterfaceFunctionSolverDefined< FSIOperator >::checkFunction( BCVectorInterface& base )
{
    boost::shared_ptr< MethodType > operMethod = boost::dynamic_pointer_cast< MethodType > ( M_physicalSolver );

    switch ( M_FSIFunction )
    {
    case DerFluidLoadToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          DerFluidLoadToFluid" << "\n";
#endif

        break;

    case DerFluidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          DerFluidLoadToStructure" << "\n";
#endif
        if ( !operMethod->isSolid() )
            return;

        operMethod->setDerFluidLoadToStructure( operMethod->sigmaSolidRepeated() );

        base = *operMethod->bcvDerFluidLoadToStructure();

        break;

    case DerHarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          DerHarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setDerHarmonicExtensionVelToFluid( operMethod->derVeloFluidMesh() );

        base = *operMethod->bcvDerHarmonicExtensionVelToFluid();

        break;

    case DerStructureDispToSolid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          DerStructureDispToSolid" << "\n";
#endif

        break;

    case FluidInterfaceDisp:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          FluidInterfaceDisp" << "\n";
#endif

        //operMethod->FluidInterfaceDisp( (LifeV::Vector&) operMethod->lambdaFluidRepeated() );

        //base = *operMethod->bcvFluidInterfaceDisp();

        break;

    case FluidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          FluidLoadToStructure" << "\n";
#endif

        if ( !operMethod->isSolid() )
            return;

        operMethod->setFluidLoadToStructure( operMethod->sigmaSolidRepeated() );

        base = *operMethod->bcvFluidLoadToStructure();

        break;

    case HarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          HarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setHarmonicExtensionVelToFluid( operMethod->veloFluidMesh() );

        base = *operMethod->bcvHarmonicExtensionVelToFluid();

        break;

    case SolidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          SolidLoadToStructure" << "\n";
#endif
        if ( !operMethod->isFluid() )
            return;

        operMethod->setSolidLoadToStructure( operMethod->minusSigmaFluidRepeated() );

        base = *operMethod->bcvSolidLoadToStructure();

        break;

    case StructureDispToHarmonicExtension:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          StructureDispToHarmonicExtension" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setStructureDispToHarmonicExtension( operMethod->lambdaFluidRepeated() );

        base = *operMethod->bcvStructureDispToHarmonicExtension();

        break;

    case StructureDispToSolid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          StructureDispToSolid" << "\n";
#endif

        break;

    case StructureToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          StructureToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setStructureToFluid( operMethod->veloFluidMesh() );
        operMethod->setStructureToFluidParametres();

        base = *operMethod->bcvStructureToFluid();

        break;

    default:

        std::cout << " !!! Error: " << M_FSIFunction << " is not available as a BCVectorInterface !!!" << std::endl;

        break;
    }
}

template< class MethodType >
inline void BCInterfaceFunctionSolverDefined< FSIOperator >::checkFunction( BCVector& base )
{
    boost::shared_ptr< MethodType > operMethod = boost::dynamic_pointer_cast< MethodType > ( M_physicalSolver );

    switch ( M_FSIFunction )
    {
    case RobinWall:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterfaceFunctionSolverDefined::checkFunction                          RobinWall" << "\n";
#endif

        if ( !operMethod->isSolid() )
            return;

        // Define the vectors
        M_robinRHS.reset( new physicalSolver_Type::vector_Type( operMethod->dFESpace().map(), Repeated, Zero ) );
        M_robinAlphaCoefficient.reset( new physicalSolver_Type::vector_Type( operMethod->dFESpace().map(), Repeated, Zero ) );
        M_robinBetaCoefficient.reset( new physicalSolver_Type::vector_Type( operMethod->dFESpace().map(), Repeated, Zero ) );

        // Set the vectors (still empty)
        base.setRhsVector( *M_robinRHS, operMethod->dFESpace().dof().numTotalDof(), 0 );
        base.setRobinCoeffVector( *M_robinAlphaCoefficient );
        base.setBetaCoeffVector( *M_robinBetaCoefficient );

        // Set the physical solver in the Robin functions for alpha and beta
        for ( UInt i( 0 ); i < M_vectorFunctionRobin.size(); ++i )
        {
            boost::shared_ptr< BCInterfaceFunctionParserSolver< physicalSolver_Type > > castedFunctionSolver =
                boost::dynamic_pointer_cast< BCInterfaceFunctionParserSolver< physicalSolver_Type > > ( M_vectorFunctionRobin[i] );

            if ( castedFunctionSolver != 0 )
                castedFunctionSolver->setPhysicalSolver( M_physicalSolver );
        }

        break;

    default:

        std::cout << " !!! Error: " << M_FSIFunction << " is not available as a BCVector !!!" << std::endl;

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterfaceFunctionSolverDefined_H */
