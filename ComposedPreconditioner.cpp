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
   \file ComposedPreconditioner.cpp
   \author Simone Deparis <simone.deparis@epfl.ch>
   \date 2009-05-07
 */

#include <lifemc/lifealg/ComposedPreconditioner.hpp>
#include <life/lifealg/IfpackPreconditioner.hpp>
#include <life/lifealg/MLPreconditioner.hpp>



namespace LifeV
{

ComposedPreconditioner::ComposedPreconditioner( boost::shared_ptr<Epetra_Comm> comm):
    super (comm ),
    M_Prec(new prec_raw_type(comm)),
    M_OperVector(0)
    //M_precType()
{
}

ComposedPreconditioner::ComposedPreconditioner(ComposedPreconditioner& P):
    super(P, boost::dynamic_pointer_cast<ComposedOperator<Ifpack_Preconditioner> >(P.getPrecPtr())->getCommPtr()),
    M_Prec(new prec_raw_type(*boost::dynamic_pointer_cast<prec_raw_type>(P.getPrecPtr()))),
    M_OperVector(P.getOperVector())
    //M_precType(P.precType())
{
    //    *M_Prec=*P.getPrec();
}
ComposedPreconditioner::~ComposedPreconditioner()
{}


void
ComposedPreconditioner::setDataFromGetPot( const GetPot&      dataFile,
                                           const std::string& section )
{
    createList( M_List, dataFile, section, "Composed" );
}

void
ComposedPreconditioner::createList(       list_Type& /*list*/,
                                    const GetPot&      dataFile,
                                    const std::string& section,
                                    const std::string& subSection )
{
    ASSERT( !M_Prec->getNumber(), "Error, when initializing the preconditioner, it must be empty" );

    for( UInt i(0); i < dataFile.vector_variable_size( ( section + "/" + subSection + "/list" ).data() ); ++i )
    {
        epetra_prec_type tmp( PRECFactory::instance().createObject( dataFile( ( section + "/" + subSection + "/list" ).data(), "ML", i ) ) );
        M_Prec->push_back(tmp);
        M_Prec->P()[i]->createList(M_Prec->P()[i]->list(), dataFile, section, dataFile( ( section + "/" + subSection + "/sections" ).data(), "ML", i ));
    }
}

int
ComposedPreconditioner::buildPreconditioner(operator_type& oper)
{
    //M_Prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back(oper);
}

int
ComposedPreconditioner::buildPreconditioner(operator_type& oper,
                                        const bool useInverse,
                                        const bool useTranspose)
{
    //M_Prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back(oper, useInverse, useTranspose);
}

int
ComposedPreconditioner::createPrec(operator_type& oper,
                               boost::shared_ptr<EpetraPreconditioner> & prec )
{
    prec->buildPreconditioner( oper );
}


int
ComposedPreconditioner::push_back(operator_type& oper,
                              const bool useInverse,
                              const bool useTranspose
                              )
{
    if(!M_Prec.get())
        M_Prec.reset(new prec_raw_type(M_displayer.comm()));
    M_OperVector.push_back(oper);
    Chrono chrono;
    epetra_prec_type prec;

    this->M_displayer.leaderPrint(std::string("ICP-  Computing prec. factorization, type:") + M_Prec->P()[M_OperVector.size()-1]->precType() + " ...        ");
    chrono.start();
    createPrec(oper, M_Prec->P()[M_OperVector.size()-1]);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());
    M_Prec->replace(prec, useInverse, useTranspose);// \TODO to reset as push_back
    if( M_Prec->P().size() == M_OperVector.size() )
        this->M_preconditionerCreated=true;
    return EXIT_SUCCESS;
}

int
ComposedPreconditioner::replace(operator_type& oper,
                            const UInt index,
                            const bool useInverse,
                            const bool useTranspose)
{
    ASSERT(index <= M_OperVector.size(), "ComposedPreconditioner::replace: index too large");

    M_OperVector[index] = oper;
    Chrono chrono;
    //ifpack_prec_type prec;
    this->M_displayer.leaderPrint(std::string("ICP-  Computing prec. factorization, type:") + (M_Prec->P()[index]->precType()) + " ...        ");
    chrono.start();
    createPrec(oper, M_Prec->P()[index]);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());

    M_Prec->replace(M_Prec->P()[index], index, useInverse, useTranspose);

    return EXIT_SUCCESS;
}


double
ComposedPreconditioner::Condest()
{
    return M_Prec->Condest();
}

EpetraPreconditioner::prec_raw_type*
ComposedPreconditioner::getPrec()
{
    return M_Prec.get();
}

void
ComposedPreconditioner::precReset()
{
    //M_OperVector.reset();
    M_Prec.reset();

    this->M_preconditionerCreated = false;
}

bool ComposedPreconditioner::registerComposed = PRECFactory::instance().registerProduct( "Composed", &ComposedPreconditioner::createComposedPreconditioner );

} // namespace LifeV
