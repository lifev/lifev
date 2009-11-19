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

#include <lifemc/lifefem/BCInterface_FSI.hpp>

namespace LifeV {

// ===================================================
// Constructors
// ===================================================
BCInterface_FSI< FSIOperator >::BCInterface_FSI() :
    M_operator      (),
    M_baseString    (),
    M_base          (),
    M_mapMethod     (),
    M_mapFunction   ()
{

#ifdef DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI()" << "\n";
#endif

}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const BCInterface_Data< FSIOperator >& data ) :
    M_operator      (),
    M_baseString    (),
    M_base          (),
    M_mapMethod     (),
    M_mapFunction   ()
{

#ifdef DEBUG
    Debug( 5025 ) << "BCInterface_FSI::BCInterface_FSI( data )" << "\n";
#endif

    this->SetData( data );
}

BCInterface_FSI< FSIOperator >::BCInterface_FSI( const BCInterface_FSI& fsi ) :
    M_operator    ( fsi.M_operator ),
    M_baseString  ( fsi.M_baseString ),
    M_base        ( fsi.M_base ),
    M_mapMethod   ( fsi.M_mapMethod ),
    M_mapFunction ( fsi.M_mapFunction )
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
        M_baseString  = fsi.M_baseString;
        M_base        = fsi.M_base;
        M_mapMethod   = fsi.M_mapMethod;
        M_mapFunction = fsi.M_mapFunction;
    }

    return *this;
}

void BCInterface_FSI< FSIOperator >::SetData( const BCInterface_Data< FSIOperator >& data )
{

#ifdef DEBUG
    Debug( 5025 ) << "BCInterface_FSIFunctionFile::setData" << "\n";
#endif

    M_operator   = data.GetOperator();
    M_baseString = data.GetBaseString();

    this->CheckMethod();
}

bool BCInterface_FSI< FSIOperator >::Compare( const BCInterface_Data< FSIOperator >& data )
{
    return M_baseString.compare( data.GetBaseString() ) == 0; //&& add compare for Operator!
}

// ===================================================
// Private functions
// ===================================================
inline void BCInterface_FSI< FSIOperator >::CheckMethod()
{
    //Set mapMethod
    M_mapMethod["exactJacobian"]   = EXACTJACOBIAN;
    M_mapMethod["fixedPoint"]      = FIXEDPOINT;
    M_mapMethod["monolithic"]      = MONOLITHIC;
    M_mapMethod["steklovPoincare"] = STEKLOVPOINCARE;

    switch ( M_mapMethod[M_operator->method()] )
    {
        case EXACTJACOBIAN:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            exactJacobian" << "\n";
#endif

            CheckFunction< exactJacobian > ();

            break;

        case FIXEDPOINT:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            fixedPoint" << "\n";
#endif

            CheckFunction< fixedPoint > ();

            break;

        case MONOLITHIC:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            monolithic" << "\n";
#endif

            //CheckFunction<monolithic>();

            break;

        case STEKLOVPOINCARE:

#ifdef DEBUG
            Debug( 5025 ) << "BCInterface_FSI::checkMethod                            steklovPoincare" << "\n";
#endif

            //CheckFunction<steklovPoincare>();

            break;
    }
}

} // Namespace LifeV
