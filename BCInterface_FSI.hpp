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
 *  @brief BCInterface_FSI
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 23-04-2009
 */

#ifndef BCInterface_FSI_H
#define BCInterface_FSI_H 1

#include <lifemc/lifesolver/BCInterface_Definitions.hpp>
#include <lifemc/lifesolver/BCInterface_Data.hpp>

#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
//#include <life/lifesolver/steklovPoincareBase.hpp>
#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV {

//! BCInterface_FSI Fake class for non-FSI problems.
/*!
 *  @author Cristiano Malossi
 */
template< class Operator >
class BCInterface_FSI
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_Data< Operator >                                          Data_Type;
    typedef BCVectorInterface                                                     BCFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    BCInterface_FSI() {}
    BCInterface_FSI( const Data_Type& /*data*/) {}
    BCInterface_FSI( const BCInterface_FSI& /*fsi*/) {}

    ~BCInterface_FSI() {}

    //@}


    //! @name Methods
    //@{

    BCInterface_FSI& operator=( const BCInterface_FSI& /*fsi*/) {}
    void SetData( const Data_Type& /*data*/) {}

    //@}


    //! @name Get functions
    //@{

    BCFunction_Type& GetBase()
    {
        return *M_base;
    }

    //@}

private:

    boost::shared_ptr< BCFunction_Type > M_base;
};

//! BCInterface_FSI - specialized template implementation for FSI problems.
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
class BCInterface_FSI< FSIOperator >
{
public:

    //! @name Type definitions
    //@{

    typedef BCInterface_Data< FSIOperator >                                       Data_Type;
    typedef BCVectorInterface                                                     BCFunction_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_FSI();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_FSI( const Data_Type& data );

    //! Copy constructor
    /*!
     * @param fsiOperator BCInterface_FSI
     */
    BCInterface_FSI( const BCInterface_FSI& fsi );

    //! Destructor
    ~BCInterface_FSI() {}

    //@}


    //! @name Methods
    //@{

    //! Operator =
    /*!
     * @param fsiOperator BCInterface_FSI
     * @return reference to a copy of the class
     */
    BCInterface_FSI& operator=( const BCInterface_FSI& fsi );

    //! Set data
    /*!
     * @param data BC data loaded from GetPot file
     */
    void SetData( const Data_Type& data );

    //@}


    //! @name Get functions
    //@{

    //! Get the base of the boundary condition
    BCFunction_Type& GetBase();

    //@}

private:

    //! @name Private functions
    //@{

    template< class method >
    inline void CheckFunction( const Data_Type& data );

    //@}

    enum FSIMethod
    {
        EXACTJACOBIAN, FIXEDPOINT, MONOLITHIC, STEKLOVPOINCARE
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

    boost::shared_ptr< FSIOperator >       M_operator;
    boost::shared_ptr< BCFunction_Type >   M_base;
};

// ===================================================
// Private functions
// ===================================================
template< class method >
inline void BCInterface_FSI< FSIOperator >::CheckFunction( const Data_Type& data )
{
    method *operMethod = dynamic_cast< method * > ( &*M_operator );

    //Set mapFunction
    std::map< std::string, FSIFunction > mapFunction;
    mapFunction["DerFluidLoadToFluid"]              = DerFluidLoadToFluid;
    mapFunction["DerFluidLoadToStructure"]          = DerFluidLoadToStructure;
    mapFunction["DerHarmonicExtensionVelToFluid"]   = DerHarmonicExtensionVelToFluid;
    mapFunction["DerStructureDispToSolid"]          = DerStructureDispToSolid;
    mapFunction["FluidInterfaceDisp"]               = FluidInterfaceDisp;
    mapFunction["FluidLoadToStructure"]             = FluidLoadToStructure;
    mapFunction["HarmonicExtensionVelToFluid"]      = HarmonicExtensionVelToFluid;
    mapFunction["SolidLoadToStructure"]             = SolidLoadToStructure;
    mapFunction["StructureDispToHarmonicExtension"] = StructureDispToHarmonicExtension;
    mapFunction["StructureDispToSolid"]             = StructureDispToSolid;
    mapFunction["StructureToFluid"]                 = StructureToFluid;

    switch ( mapFunction[ data.GetBaseString() ] )
    {
        case DerFluidLoadToFluid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          DerFluidLoadToFluid" << "\n";
#endif

            break;

        case DerFluidLoadToStructure:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          DerFluidLoadToStructure" << "\n";
#endif
            if ( !M_operator->isSolid() )
                return;

            operMethod->setDerFluidLoadToStructure( M_operator->sigmaSolidRepeated() );

            M_base = operMethod->bcvDerFluidLoadToStructure();

            break;

        case DerHarmonicExtensionVelToFluid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          DerHarmonicExtensionVelToFluid" << "\n";
#endif

            if ( !M_operator->isFluid() )
                return;

            operMethod->setDerHarmonicExtensionVelToFluid( M_operator->derVeloFluidMesh() );

            M_base = operMethod->bcvDerHarmonicExtensionVelToFluid();

            break;

        case DerStructureDispToSolid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          DerStructureDispToSolid" << "\n";
#endif

            break;

        case FluidInterfaceDisp:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          FluidInterfaceDisp" << "\n";
#endif

            //operMethod->FluidInterfaceDisp( (LifeV::Vector&) M_operator->lambdaFluidRepeated() );

            //M_base = operMethod->bcvFluidInterfaceDisp();

            break;

        case FluidLoadToStructure:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          FluidLoadToStructure" << "\n";
#endif

            if ( !M_operator->isSolid() )
                return;

            operMethod->setFluidLoadToStructure( M_operator->sigmaSolidRepeated() );

            M_base = operMethod->bcvFluidLoadToStructure();

            break;

        case HarmonicExtensionVelToFluid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          HarmonicExtensionVelToFluid" << "\n";
#endif

            if ( !M_operator->isFluid() )
                return;

            M_operator->setHarmonicExtensionVelToFluid( M_operator->veloFluidMesh() );

            M_base = M_operator->bcvHarmonicExtensionVelToFluid();

            break;

        case SolidLoadToStructure:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          SolidLoadToStructure" << "\n";
#endif
            if ( !M_operator->isFluid() )
                return;

            M_operator->setSolidLoadToStructure( M_operator->minusSigmaFluidRepeated() );

            M_base = M_operator->bcvSolidLoadToStructure();

            break;

        case StructureDispToHarmonicExtension:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          StructureDispToHarmonicExtension" << "\n";
#endif

            if ( !M_operator->isFluid() )
                return;

            operMethod->setStructureDispToHarmonicExtension( M_operator->lambdaFluidRepeated() );

            M_base = operMethod->bcvStructureDispToHarmonicExtension();

            break;

        case StructureDispToSolid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          StructureDispToSolid" << "\n";
#endif

            break;

        case StructureToFluid:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkFunction                          StructureToFluid" << "\n";
#endif

            if ( !M_operator->isFluid() )
                return;

            M_operator->setStructureToFluid( M_operator->veloFluidMesh() );
            M_operator->setStructureToFluidParametres();

            M_base = M_operator->bcvStructureToFluid();

            break;
    }
}

} // Namespace LifeV

#endif /* BCInterface_FSI_H */
