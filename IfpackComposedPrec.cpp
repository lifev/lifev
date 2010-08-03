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
   \file IfpackComposedPrec.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2009-05-07
 */

#include <lifemc/lifealg/IfpackComposedPrec.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>



namespace LifeV
{


IfpackComposedPrec::IfpackComposedPrec(const Epetra_Comm* comm):
  super (comm),
  M_Prec(new prec_raw_type()),
  M_OperVector(0),
  M_precType()
{
}

IfpackComposedPrec::IfpackComposedPrec(IfpackComposedPrec& P):
    super(P, &P.getPrec()->Comm()),
    M_Prec(new prec_raw_type(*boost::dynamic_pointer_cast<prec_raw_type>(P.getPrecPtr()))),
    M_OperVector(P.getOperVector()),
    M_precType(P.precType())
{
    //    *M_Prec=*P.getPrec();
}
IfpackComposedPrec::~IfpackComposedPrec()
{}


void
IfpackComposedPrec::setDataFromGetPot( const GetPot& dataFile,
                                       const std::string& section )
{

    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters
    createIfpackList(dataFile, section, this->M_List);

    M_overlapLevel = this->M_List.get("overlap level", -1);
    M_precType     = this->M_List.get("prectype", "Amesos");

}

int
IfpackComposedPrec::buildPreconditioner(operator_type& oper)
{
    M_Prec.reset(new prec_raw_type());
    return push_back(oper);
}

int
IfpackComposedPrec::buildPreconditioner(operator_type& oper,
                                        const bool useInverse,
                                        const bool useTranspose)
{
    M_Prec.reset(new prec_raw_type());
    return push_back(oper, useInverse, useTranspose);
}

int
IfpackComposedPrec::createIfpackPrec(operator_type& oper,
                                     prec_raw_type::prec_type& prec)
{
    M_overlapLevel = this->M_List.get("overlap level", -1);
    M_precType     = this->M_List.get("prectype", "Amesos");

    Ifpack factory;


    prec.reset(factory.Create(M_precType, oper->getMatrixPtr().get(), M_overlapLevel));

    if ( !prec.get() )
    {
        ERROR_MSG( "Preconditioner not set, something went wrong in its computation\n" );
    }

    IFPACK_CHK_ERR(prec->SetParameters(this->M_List));
    IFPACK_CHK_ERR(prec->Initialize());
    IFPACK_CHK_ERR(prec->Compute());

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

int
IfpackComposedPrec::push_back(operator_type& oper,
                              const bool useInverse,
                              const bool useTranspose)
{
    if(!M_Prec.get())
        M_Prec.reset(new prec_raw_type());
    M_OperVector.push_back(oper);
    Chrono chrono;
    prec_raw_type::prec_type prec;

    this->M_displayer.leaderPrint("ICP-  Computing prec. factorization ...        ");
    chrono.start();
    createIfpackPrec(oper, prec);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());
    M_Prec->push_back(prec, useInverse, useTranspose);

    return EXIT_SUCCESS;
}

int
IfpackComposedPrec::replace(operator_type& oper,
                            const UInt index,
                            const bool useInverse,
                            const bool useTranspose)
{
    ASSERT(index <= M_OperVector.size(), "IfpackComposedPrec::replace: index too large");

    M_OperVector[index] = oper;
    Chrono chrono;
    prec_raw_type::prec_type prec;
    this->M_displayer.leaderPrint("ICP-  Computing prec. factorization ...        ");
    chrono.start();
    createIfpackPrec(oper, prec);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());

    M_Prec->replace(prec, index, useInverse, useTranspose);

    return EXIT_SUCCESS;
}


double
IfpackComposedPrec::Condest()
{
    return M_Prec->Condest();
}

EpetraPreconditioner::prec_raw_type*
IfpackComposedPrec::getPrec()
{
    return M_Prec.get();
}

void
IfpackComposedPrec::precReset()
{
    M_Oper.reset();
    M_Prec.reset();

    this->M_preconditionerCreated = false;
}

} // namespace LifeV
