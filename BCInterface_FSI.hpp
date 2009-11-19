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

#include <lifemc/lifefem/BCInterface_Definitions.hpp>
#include <lifemc/lifefem/BCInterface_Data.hpp>

#include <life/lifesolver/exactJacobianBase.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
//#include <life/lifesolver/steklovPoincareBase.hpp>

namespace LifeV {

//! BCInterface_FSI Fake class for non-FSI problems.
/*!
 *  @author Cristiano Malossi
 */
template< class Operator >
class BCInterface_FSI
{
public:

    //! @name Constructors & Destructor
    //@{

    BCInterface_FSI() {}
    BCInterface_FSI( const BCInterface_Data< Operator >& /*data*/) {}
    BCInterface_FSI( const BCInterface_FSI& /*fsi*/) {}

    ~BCInterface_FSI() {}

    //@}


    //! @name Methods
    //@{

    BCInterface_FSI& operator=( const BCInterface_FSI& /*fsi*/) {}
    void SetData( const BCInterface_Data< Operator >& /*data*/) {}

    bool Compare( const BCInterface_Data< Operator >& /*data*/)
    {
        return true;
    }

    //@}


    //! @name Get functions
    //@{

    BCVectorInterface& GetBase()
    {
        return *M_base;
    }

    //@}

private:

    boost::shared_ptr< BCVectorInterface > M_base;
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
 *	- MONOLITHIC		(not working)
 *	- STEKLOVPOINCARE 	(not working)
 *
 *	To get the base for the boundary condition call the getBase function.
 */
template< >
class BCInterface_FSI< FSIOperator >
//     :
//     public LifeV::Application
{
public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    BCInterface_FSI();

    //! Constructor
    /*!
     * @param data BC data loaded from GetPot file
     */
    BCInterface_FSI( const BCInterface_Data< FSIOperator >& data );

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
    void SetData( const BCInterface_Data< FSIOperator >& data );

    //! Compare function
    /*!
     * @param data BC data loaded from GetPot file
     * @return true if the functions are equal, false if they aren't
     */
    bool Compare( const BCInterface_Data< FSIOperator >& data );

    //@}


    //! @name Get functions
    //@{

    //! Get the base of the boundary condition
    BCVectorInterface& GetBase()
    {
        return *M_base;
    }

    //@}

private:

    //! @name Private functions
    //@{

    inline void CheckMethod();

    template< class method >
    inline void CheckFunction();

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
    std::string                            M_baseString;
    boost::shared_ptr< BCVectorInterface > M_base;

    std::map< std::string, FSIMethod >     M_mapMethod;
    std::map< std::string, FSIFunction >   M_mapFunction;
};

// ===================================================
// Private functions
// ===================================================
template< class method >
inline void BCInterface_FSI< FSIOperator >::CheckFunction()
{
    method *operMethod = dynamic_cast< method * > ( &*M_operator );

    //Set mapFunction
    M_mapFunction["DerFluidLoadToFluid"]              = DerFluidLoadToFluid;
    M_mapFunction["DerFluidLoadToStructure"]          = DerFluidLoadToStructure;
    M_mapFunction["DerHarmonicExtensionVelToFluid"]   = DerHarmonicExtensionVelToFluid;
    M_mapFunction["DerStructureDispToSolid"]          = DerStructureDispToSolid;
    M_mapFunction["FluidInterfaceDisp"]               = FluidInterfaceDisp;
    M_mapFunction["FluidLoadToStructure"]             = FluidLoadToStructure;
    M_mapFunction["HarmonicExtensionVelToFluid"]      = HarmonicExtensionVelToFluid;
    M_mapFunction["SolidLoadToStructure"]             = SolidLoadToStructure;
    M_mapFunction["StructureDispToHarmonicExtension"] = StructureDispToHarmonicExtension;
    M_mapFunction["StructureDispToSolid"]             = StructureDispToSolid;
    M_mapFunction["StructureToFluid"]                 = StructureToFluid;

    switch ( M_mapFunction[M_baseString] )
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
