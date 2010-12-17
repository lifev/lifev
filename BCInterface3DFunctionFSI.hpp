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
 *  @brief File containing the BCInterface3DFunctionFSI class
 *
 *  @date 23-04-2009
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#ifndef BCInterface3DFunctionFSI_H
#define BCInterface3DFunctionFSI_H 1

#include <lifemc/lifesolver/BCInterface3DDefinitions.hpp>
#include <lifemc/lifesolver/BCInterface3DData.hpp>

#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
//#include <life/lifesolver/steklovPoincareBase.hpp>
#include <lifemc/lifesolver/MonolithicGE.hpp>
#include <lifemc/lifesolver/MonolithicGI.hpp>

namespace LifeV
{

//! BCInterface3DFunctionFSI Fake class for non-FSI problems.
/*!
 *  @author Cristiano Malossi
 */
template< class PhysicalSolverType >
class BCInterface3DFunctionFSI
{
public:

    //! @name Type definitions
    //@{

    typedef PhysicalSolverType                                                    physicalSolver_Type;
    typedef BCInterface3DData                                                     data_Type;
    typedef BCVectorInterface                                                     bcFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    explicit BCInterface3DFunctionFSI() {}
    explicit BCInterface3DFunctionFSI( const data_Type& /*data*/) {}

    virtual ~BCInterface3DFunctionFSI() {}

    //@}


    //! @name Methods
    //@{

    void exportData( data_Type& /*data*/ ) {}
    void checkMethod( const boost::shared_ptr< physicalSolver_Type >& /*physicalSolver*/ ) {}

    //@}


    //! @name Set Methods
    //@{

    void setData( const data_Type& /*data*/) {}

    //@}


    //! @name Get Methods
    //@{

    bcFunction_Type& base() { return *M_base; }

    //@}

private:

    //! @name Unimplemented Methods
    //@{

    BCInterface3DFunctionFSI( const BCInterface3DFunctionFSI& fsi);

    BCInterface3DFunctionFSI& operator=( const BCInterface3DFunctionFSI& fsi );

    //@}

    boost::shared_ptr< bcFunction_Type > M_base;
};

//! BCInterface3DFunctionFSI - specialized template implementation for FSI problems.
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_PhysicalCoupling class provides a general interface between the
 *  MS_Algorithm and all the coupling conditions.
 *
 *  This class allows to use impose interface conditions for FSI problems.
 *
 *  <b>DETAILS:</b>
 *
 *  The constructor of the class takes a string contains the ID of the interface condition to impose,
 *  and the FSIOperator. The list of available conditions is the FSIFunction variable. These are:
 *
 *	- DerFluidLoadToFluid,					(not implemented)
 *	- DerFluidLoadToStructure,
 *	- DerHarmonicExtensionVelToFluid,
 *	- DerStructureDispToSolid,				(not implemented)
 *	- FluidInterfaceDisp,					(not working)
 *	- FluidLoadToStructure,
 *	- HarmonicExtensionVelToFluid,
 *	- SolidLoadToStructure,
 *	- StructureDispToHarmonicExtension,
 *	- StructureDispToSolid, 				(not implemented)
 *	- StructureToFluid
 *
 *	The class automatically recognize which method is used among:
 *
 *	- EXACTJACOBIAN
 *	- FIXEDPOINT
 *	- MONOLITHIC
 *	- STEKLOVPOINCARE 	(not working)
 *
 *	To get the base for the boundary condition call the getBase function.
 */
template< >
class BCInterface3DFunctionFSI< FSIOperator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface3DData                                                      data_Type;
    typedef BCVectorInterface                                                     bcFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    explicit BCInterface3DFunctionFSI();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    explicit BCInterface3DFunctionFSI( const data_Type& data );

    //! Destructor
    virtual ~BCInterface3DFunctionFSI() {}

    //@}


    //! @name Methods
    //@{

    //! Copy the stored parameters in the data container
    /*!
     * @param data BC data loaded from GetPot file
     */
    void exportData( data_Type& data );

    //! Check method passing the operator
    /*!
     * @param physicalSolver FSIOperator
     */
    void checkMethod( const boost::shared_ptr< FSIOperator >& physicalSolver );

    //@}


    //! @name Set methods
    //@{

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void setData( const data_Type& data );

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

    BCInterface3DFunctionFSI( const BCInterface3DFunctionFSI& fsi);

    BCInterface3DFunctionFSI& operator=( const BCInterface3DFunctionFSI& fsi );

    //@}


    //! @name Private Methods
    //@{

    template< class method >
    void checkFunction( const boost::shared_ptr< FSIOperator >& physicalSolver );

    //@}

    enum FSIMethod
    {
        EXACTJACOBIAN, FIXEDPOINT, MONOLITHIC_GE, MONOLITHIC_GI, STEKLOVPOINCARE
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

    FSIFunction                               M_FSIFunction;

    // These are required since FSI BC are applied a posteriori
    bcName_Type                                    M_name;
    BCFlag                                    M_flag;
    bcType_Type                                    M_type;
    bcMode_Type                                    M_mode;
    BCComV                                    M_comV;

    boost::shared_ptr< bcFunction_Type >      M_base;
};

// ===================================================
// Private functions
// ===================================================
template< class method >
inline void BCInterface3DFunctionFSI< FSIOperator >::checkFunction( const boost::shared_ptr< FSIOperator >& physicalSolver )
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
        if ( !physicalSolver->isSolid() )
            return;

        operMethod->setDerFluidLoadToStructure( physicalSolver->sigmaSolidRepeated() );

        M_base = operMethod->bcvDerFluidLoadToStructure();

        break;

    case DerHarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          DerHarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !physicalSolver->isFluid() )
            return;

        operMethod->setDerHarmonicExtensionVelToFluid( physicalSolver->derVeloFluidMesh() );

        M_base = operMethod->bcvDerHarmonicExtensionVelToFluid();

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

        //operMethod->FluidInterfaceDisp( (LifeV::Vector&) physicalSolver->lambdaFluidRepeated() );

        //M_base = operMethod->bcvFluidInterfaceDisp();

        break;

    case FluidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          FluidLoadToStructure" << "\n";
#endif

        if ( !physicalSolver->isSolid() )
            return;

        operMethod->setFluidLoadToStructure( physicalSolver->sigmaSolidRepeated() );

        M_base = operMethod->bcvFluidLoadToStructure();

        break;

    case HarmonicExtensionVelToFluid:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          HarmonicExtensionVelToFluid" << "\n";
#endif

        if ( !physicalSolver->isFluid() )
            return;

        physicalSolver->setHarmonicExtensionVelToFluid( physicalSolver->veloFluidMesh() );

        M_base = physicalSolver->bcvHarmonicExtensionVelToFluid();

        break;

    case SolidLoadToStructure:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          SolidLoadToStructure" << "\n";
#endif
        if ( !physicalSolver->isFluid() )
            return;

        physicalSolver->setSolidLoadToStructure( physicalSolver->minusSigmaFluidRepeated() );

        M_base = physicalSolver->bcvSolidLoadToStructure();

        break;

    case StructureDispToHarmonicExtension:

#ifdef HAVE_LIFEV_DEBUG
        Debug( 5025 ) << "BCInterface3DFunctionFSI::checkFunction                          StructureDispToHarmonicExtension" << "\n";
#endif

        if ( !physicalSolver->isFluid() )
            return;

        operMethod->setStructureDispToHarmonicExtension( physicalSolver->lambdaFluidRepeated() );

        M_base = operMethod->bcvStructureDispToHarmonicExtension();

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

        if ( !physicalSolver->isFluid() )
            return;

        physicalSolver->setStructureToFluid( physicalSolver->veloFluidMesh() );
        physicalSolver->setStructureToFluidParametres();

        M_base = physicalSolver->bcvStructureToFluid();

        break;
    }
}

} // Namespace LifeV

#endif /* BCInterface3DFunctionFSI_H */
