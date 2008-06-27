/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-11-09

  Copyright (C) 2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file EpetraPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2006-11-09
 */

#include "EpetraPreconditioner.hpp"


namespace LifeV
{
namespace Epetra
{

Preconditioner::Preconditioner():
  M_Prec(),
  M_Oper()
{
}

Preconditioner::Preconditioner(operator_type& oper):
  M_Prec(),
  M_Oper()
{
    buildPreconditioner(oper);
}


void Preconditioner::setDataFromGetPot( const GetPot& dataFile, const std::string& section )
{
    M_precType = dataFile((section + "/prectype").data(),"Amesos");

    M_overlapLevel       = dataFile((section + "/overlap").data(),     4);
    double dropTolerance = dataFile((section + "/droptol").data(),     1e-5);
    double levelOfFill   = dataFile((section + "/fill").data(),        4.);
    double athr   = dataFile((section + "/athr").data(),        0.);
    double rthr   = dataFile((section + "/rthr").data(),        1.);
    //double relax_value   = dataFile((section + "/relax_value").data(), 0.);
    //int    localParts    = dataFile((section + "/localparts").data(),  4);

    M_List.set("fact: drop tolerance",     dropTolerance);
    M_List.set("fact: ilut level-of-fill", levelOfFill);
    M_List.set("fact: absolute threshold", athr);
    M_List.set("fact: relative threshold", rthr);
    //M_List.set("fact: level-of-fill",      levelOfFill);
    //M_List.set("fact: relax value",      relax_value);

}

int Preconditioner::buildPreconditioner(operator_type& oper)
{
    M_Oper = oper;

    //List.set("schwarz: combine mode", "Zero"); //
    M_List.set("schwarz: filter singletons", true);

//    List.set("amesos: solver type", "Amesos_Lapack");

//     M_List.set("PrintTiming", true);
//     M_List.set("PrintStatus", true);

    Ifpack factory;

    M_Prec.reset(factory.Create(M_precType, &M_Oper->getEpetraMatrix(), M_overlapLevel));
//    M_Prec.reset(new prec_type(&A.getEpetraMatrix(), OverlapLevel));
    if ( !M_Prec.get() )
    { //! if not filled, I do not know how to diagonalize.
      ERROR_MSG( "Preconditioner not setted, something went wrong in its computation\n" );
    }

    IFPACK_CHK_ERR(M_Prec->SetParameters(M_List));
    IFPACK_CHK_ERR(M_Prec->Initialize());
    IFPACK_CHK_ERR(M_Prec->Compute());
    return EXIT_SUCCESS;
}

double Preconditioner::Condest()
{
    return M_Prec->Condest();
}

Preconditioner::prec_raw_type* Preconditioner::getPrec()
{
    return M_Prec.get();
}

void
Preconditioner::precReset()
{
    M_Oper.reset();
    M_Prec.reset();
}

void
Preconditioner::createList( const GetPot& /*dataFile*/ )
{
}




} // namespace Epetra
} // namespace LifeV
