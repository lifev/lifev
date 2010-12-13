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


// ===================================================
//! Constructurs and Destructor
// ===================================================

ComposedPreconditioner::ComposedPreconditioner( boost::shared_ptr<Epetra_Comm> comm):
    super_Type (comm ),
    M_prec(new prec_raw_type(comm)),
    M_operVector(0)
    //M_precType()
{
}

ComposedPreconditioner::ComposedPreconditioner(ComposedPreconditioner& P):
    super(P, boost::dynamic_pointer_cast<ComposedOperator<Ifpack_Preconditioner> >(P.getPrecPtr())->getCommPtr()),
    M_prec(new prec_raw_type(*boost::dynamic_pointer_cast<prec_raw_type>(P.getPrecPtr()))),
    M_operVector(P.getOperVector())
    //M_precType(P.precType())
{
    //    *M_prec=*P.getPrec();
}
ComposedPreconditioner::~ComposedPreconditioner()
{}


// ===================================================
//! Public Methods
// ===================================================


void
ComposedPreconditioner::setDataFromGetPot( const GetPot&      dataFile,
                                           const std::string& section )
{
	list_Type uselessList __attribute__ ((deprecated));
    createList( uselessList, dataFile, section, "Composed" );
}

void
ComposedPreconditioner::createList(       list_Type& /*list*/,
                                          const GetPot&      dataFile,
                                          const std::string& section,
                                          const std::string& subSection )
{
    //M_list.resize(listNumber);
    //! See http://trilinos.sandia.gov/packages/docs/r9.0/packages/ifpack/doc/html/index.html
    //! for more informations on the parameters
    ASSERT( !M_prec->getNumber(), "Error, when initializing the preconditioner, it must be empty" );
    for ( UInt i(0); i < dataFile.vector_variable_size( ( section + "/" + subSection + "/list" ).data() ); ++i )
    {
        epetra_prec_type tmp( PRECFactory::instance().createObject( dataFile( ( section + "/" + subSection + "/list" ).data(), "ML", i ) ) );
        M_prec->push_back(tmp);
        M_prec->P()[i]->createList(M_prec->P()[i]->list(), dataFile, section, dataFile( ( section + "/" + subSection + "/sections" ).data(), "ML", i ));
    }
}

double
ComposedPreconditioner::Condest()
{
    return M_prec->Condest();
}

EpetraPreconditioner::prec_raw_type*
ComposedPreconditioner::getPrec()
{
    return M_prec.get();
}

int
ComposedPreconditioner::buildPreconditioner(operatorPtr_Type& oper)
{
    //M_prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back(oper);
}

int
ComposedPreconditioner::buildPreconditioner(operatorPtr_Type& oper,
                                        const bool useInverse,
                                        const bool useTranspose)
{
    //M_prec.reset(new prec_raw_type(M_displayer.comm()));
    return push_back(oper, useInverse, useTranspose);
}

int
ComposedPreconditioner::push_back(operatorPtr_Type& oper,
                              const bool useInverse,
                              const bool useTranspose
                              )
{
    if(!M_prec.get())
        M_prec.reset(new prec_raw_type(M_displayer.comm()));
    M_operVector.push_back(oper);
    Chrono chrono;
    epetraPrecPtr_Type prec;

    this->M_displayer.leaderPrint(std::string("ICP-  Computing prec. factorization, type:") + M_prec->P()[M_operVector.size()-1]->precType() + " ...        ");
    chrono.start();
    createPrec(oper, M_prec->P()[M_operVector.size()-1]);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());
    M_prec->replace(prec, useInverse, useTranspose);// \TODO to reset as push_back
    if( M_prec->P().size() == M_operVector.size() )
        this->M_preconditionerCreated=true;
    return EXIT_SUCCESS;
}

int
ComposedPreconditioner::replace(operatorPtr_Type& oper,
                            const UInt index,
                            const bool useInverse,
                            const bool useTranspose)
{
    ASSERT(index <= M_operVector.size(), "ComposedPreconditioner::replace: index too large");

    M_operVector[index] = oper;
    Chrono chrono;
    //ifpack_prec_type prec;
    this->M_displayer.leaderPrint(std::string("ICP-  Computing prec. factorization, type:") + (M_prec->P()[index]->precType()) + " ...        ");
    chrono.start();
    createPrec(oper, M_prec->P()[index]);
    chrono.stop();
    this->M_displayer.leaderPrintMax("done in ", chrono.diff());

    M_prec->replace(M_prec->P()[index], index, useInverse, useTranspose);

    return EXIT_SUCCESS;
}

void
ComposedPreconditioner::precReset()
{
    //M_operVector.reset();
    M_prec.reset();

    this->M_preconditionerCreated = false;
}


// ===================================================
//! Private Methods
// ===================================================

int
ComposedPreconditioner::createPrec(operator_type& oper,
                                   boost::shared_ptr<EpetraPreconditioner> & prec )
{
    prec->buildPreconditioner( oper );
}



bool ComposedPreconditioner::registerComposed = PRECFactory::instance().registerProduct( "Composed", &ComposedPreconditioner::createComposedPreconditioner );

} // namespace LifeV
