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

#include <lifemc/lifesolver/BCInterface_FSI.hpp>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
BCInterface_FSI< FSIOperator >::BCInterface_FSI() :
    M_operator      (),
    M_base          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI()" << "\n";
#endif

}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const Data_Type& data ) :
    M_operator      (),
    M_base          ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI( data )" << "\n";
#endif

    this->SetData( data );
}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const BCInterface_FSI& fsi ) :
    M_operator    ( fsi.M_operator ),
    M_base        ( fsi.M_base )
{
}

// ===================================================
// Methods
// ===================================================
BCInterface_FSI< FSIOperator >&
BCInterface_FSI< FSIOperator >::operator=( const BCInterface_FSI& fsi )
{
    if ( this != &fsi )
    {
        M_operator    = fsi.M_operator;
        M_base        = fsi.M_base;
    }

    return *this;
}

void
BCInterface_FSI< FSIOperator >::SetData( const Data_Type& data )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 5025 ) << "BCInterface_FSIFunctionFile::setData" << "\n";
#endif

    M_operator   = data.GetOperator();

    //Set mapMethod
    std::map< std::string, FSIMethod > mapMethod;

    mapMethod["exactJacobian"]   = EXACTJACOBIAN;
    mapMethod["fixedPoint"]      = FIXEDPOINT;
    mapMethod["monolithic"]      = MONOLITHIC;
    mapMethod["steklovPoincare"] = STEKLOVPOINCARE;

    switch ( mapMethod[M_operator->data().method()] )
    {
        case EXACTJACOBIAN:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            exactJacobian" << "\n";
#endif

            CheckFunction< exactJacobian > ( data );

            break;

        case FIXEDPOINT:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            fixedPoint" << "\n";
#endif

            CheckFunction< fixedPoint > ( data );

            break;

        case MONOLITHIC:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            monolithic" << "\n";
#endif

            CheckFunction< Monolithic >( data );

            break;

        case STEKLOVPOINCARE:

#ifdef HAVE_LIFEV_DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            steklovPoincare" << "\n";
#endif

            //CheckFunction< steklovPoincare >( data );

            break;
    }
}

// ===================================================
// Get Methods
// ===================================================
BCInterface_FSI< FSIOperator >::BCFunction_Type&
BCInterface_FSI< FSIOperator >::GetBase()
{
    return *M_base;
}

} // Namespace LifeV
