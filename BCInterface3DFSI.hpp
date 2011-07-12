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
 *  @brief File containing the BCInterface3DFSI class
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3DFSI_H
#define BCInterface3DFSI_H 1

#include <life/lifesolver/FSIExactJacobian.hpp>
#include <life/lifesolver/FSIFixedPoint.hpp>

#include <lifemc/lifesolver/FSIMonolithicGE.hpp>
#include <lifemc/lifesolver/FSIMonolithicGI.hpp>

#include <lifemc/lifesolver/BCInterfaceDefinitions.hpp>
#include <lifemc/lifesolver/BCInterfaceData.hpp>

namespace LifeV
{

//! BCInterface3DFSI Fake class for non-FSI problems.
/*!
 *  @author Cristiano Malossi
 */
template< class PhysicalSolverType >
class BCInterface3DFSI
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterfaceData                                                       data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit BCInterface3DFSI() {}
    explicit BCInterface3DFSI( const data_Type& /*data*/) {}

    virtual ~BCInterface3DFSI() {}

    //@}


    //! @name Methods
    //@{

    void exportData( data_Type& /*data*/ ) {}
    void assignFunction( const boost::shared_ptr< physicalSolver_Type >& /*physicalSolver*/, BCVectorInterface& /*base*/ ) {}

    //@}


    //! @name Set Methods
    //@{

    void setData( const data_Type& /*data*/) {}

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3DFSI( const BCInterface3DFSI& fsi);

    BCInterface3DFSI& operator=( const BCInterface3DFSI& fsi );

    //@}

};

//! BCInterface3DFSI - LifeV boundary condition function wrapper for \c BCInterface3D and FSI problems
/*!
 *  @author Cristiano Malossi
 *
 *
 *  The BCInterface3DFSI class provides a general interface between the
 *  \c BCInterface3D and the default boundary condition for the \c FSIOperator.
 *
 *  <b>DETAILS:</b> <BR>
 *  The constructor of the class takes a string contains the ID of the interface condition to impose,
 *  and the FSI. The list of available conditions is the FSIFunction variable. These are:
 *
 *	<ol>
 *      <li> DerFluidLoadToFluid,					(not implemented)
 *	    <li> DerFluidLoadToStructure,
 *	    <li> DerHarmonicExtensionVelToFluid,
 *	    <li> DerStructureDispToSolid,				(not implemented)
 *	    <li> FluidInterfaceDisp,					(not working)
 *	    <li> FluidLoadToStructure,
 *	    <li> HarmonicExtensionVelToFluid,
 *	    <li> SolidLoadToStructure,
 *	    <li> StructureDispToHarmonicExtension,
 *	    <li> StructureDispToSolid, 				    (not implemented)
 *	    <li> StructureToFluid
 *  <ol>
 *
 *	The class automatically recognize which FSI algorithm is used among:
 *  <ol>
 *      <li> EXACTJACOBIAN;
 *      <li> FIXEDPOINT;
 *      <li> MONOLITHIC (both GE and GI);
 *  </ol>
 */
template< >
class BCInterface3DFSI< FSIOperator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterfaceData                                                      data_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3DFSI();

    //! Constructor
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    explicit BCInterface3DFSI( const data_Type& data );

    //! Destructor
    virtual ~BCInterface3DFSI() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the data container
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void exportData( data_Type& data );

    //! Assign a boundary function to the boundary condition vector base
    /*!
     * @param physicalSolver FSI physical solver,
     * @param base boundary condition vector base
     */
    void assignFunction( const boost::shared_ptr< FSIOperator >& physicalSolver, BCVectorInterface& base );

    //@}


    //! @name Set methods
    //@{

    //! Set data
    /*!
     * @param data boundary condition data loaded from \c GetPot file
     */
    void setData( const data_Type& data );

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3DFSI( const BCInterface3DFSI& fsi);

    BCInterface3DFSI& operator=( const BCInterface3DFSI& fsi );

    //@}


    //! @name Private Methods
    //@{

    template< class method >
    void checkFunction( const boost::shared_ptr< FSIOperator >& physicalSolver, BCVectorInterface& base );

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
        StructureToFluid
    };

    FSIFunction                                    M_FSIFunction;

    // These are required since FSI BC are applied a posteriori, when setPhysicalSolver is called
    bcName_Type                                    M_name;
    bcFlag_Type                                    M_flag;
    bcType_Type                                    M_type;
    bcMode_Type                                    M_mode;
    bcComponentsVec_Type                           M_comV;
};

// ===================================================
// Private functions
// ===================================================
template< class method >
inline void BCInterface3DFSI< FSIOperator >::checkFunction( const boost::shared_ptr< FSIOperator >& physicalSolver, BCVectorInterface& base )
{
    method *operMethod = dynamic_cast< method * > ( &*physicalSolver );

    switch ( M_FSIFunction )
    {
    case DerFluidLoadToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          DerFluidLoadToFluid" << "\n";
#endif

        break;

    case DerFluidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          DerFluidLoadToStructure" << "\n";
#endif
        if ( !operMethod->isSolid() )
            return;

        operMethod->setDerFluidLoadToStructure( operMethod->sigmaSolidRepeated() );

        base = *operMethod->bcvDerFluidLoadToStructure();

        break;

    case DerHarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          DerHarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setDerHarmonicExtensionVelToFluid( operMethod->derVeloFluidMesh() );

        base = *operMethod->bcvDerHarmonicExtensionVelToFluid();

        break;

    case DerStructureDispToSolid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          DerStructureDispToSolid" << "\n";
#endif

        break;

    case FluidInterfaceDisp:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          FluidInterfaceDisp" << "\n";
#endif

        //operMethod->FluidInterfaceDisp( (LifeV::Vector&) operMethod->lambdaFluidRepeated() );

        //base = *operMethod->bcvFluidInterfaceDisp();

        break;

    case FluidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          FluidLoadToStructure" << "\n";
#endif

        if ( !operMethod->isSolid() )
            return;

        operMethod->setFluidLoadToStructure( operMethod->sigmaSolidRepeated() );

        base = *operMethod->bcvFluidLoadToStructure();

        break;

    case HarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          HarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setHarmonicExtensionVelToFluid( operMethod->veloFluidMesh() );

        base = *operMethod->bcvHarmonicExtensionVelToFluid();

        break;

    case SolidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          SolidLoadToStructure" << "\n";
#endif
        if ( !operMethod->isFluid() )
            return;

        operMethod->setSolidLoadToStructure( operMethod->minusSigmaFluidRepeated() );

        base = *operMethod->bcvSolidLoadToStructure();

        break;

    case StructureDispToHarmonicExtension:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          StructureDispToHarmonicExtension" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setStructureDispToHarmonicExtension( operMethod->lambdaFluidRepeated() );

        base = *operMethod->bcvStructureDispToHarmonicExtension();

        break;

    case StructureDispToSolid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          StructureDispToSolid" << "\n";
#endif

        break;

    case StructureToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          StructureToFluid" << "\n";
#endif

        if ( !operMethod->isFluid() )
            return;

        operMethod->setStructureToFluid( operMethod->veloFluidMesh() );
        operMethod->setStructureToFluidParametres();

        base = *operMethod->bcvStructureToFluid();

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3DFSI_H */
