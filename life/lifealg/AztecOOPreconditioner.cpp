//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
 *  @brief AztecOO preconditioner
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 17-11-2009
 */

#include "AztecOOPreconditioner.hpp"

namespace LifeV {

AztecOOPreconditioner::AztecOOPreconditioner():
        super                   ( ),
        M_solver                ( )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::AztecOOPreconditioner() \n";
#endif

}

void
AztecOOPreconditioner::setDataFromGetPot( const GetPot&      dataFile,
                                          const std::string& section,
                                          const UInt&        /*listNumber*/ )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::setDataFromGetPot(dataFile, section) \n";
#endif

    // Preconditioner
    M_solver->getParameterList().set("precond", dataFile( ( section + "/AztecOO/precond" ).data(), "dom_decomp" ));

    // Compute the preconditioner
    M_solver->getParameterList().set("pre_calc", dataFile( ( section + "/AztecOO/pre_calc" ).data(), "calc" ));

    // Reordering
    M_solver->getParameterList().set("reorder", dataFile( ( section + "/AztecOO/reorder" ).data(), 1 ));

    // Keep the information
    M_solver->getParameterList().set("keep_info", dataFile( ( section + "/AztecOO/keep_info" ).data(), 1 ));

    // Polynomial order when using Jacobi and GS preconditioners
    //M_solver->getParameterList().set("poly_ord", dataFile( ( section + "/AztecOO/poly_ord" ).data(), 3 ));

    // Overlap level
    M_solver->getParameterList().set("overlap", dataFile( ( section + "/AztecOO/overlap" ).data(), AZ_none ));

    // Overlap type
    M_solver->getParameterList().set("type_overlap", dataFile( ( section + "/AztecOO/type_overlap" ).data(), AZ_standard ));



    // SUBDOMAIN SOLVER

    M_solver->getParameterList().set("subdomain_solve", dataFile( ( section + "/AztecOO/subdomain_solve" ).data(), "ILUT" ));

    M_solver->getParameterList().set("drop", dataFile( ( section + "/AztecOO/drop" ).data(), 1.e-5 ));

    //M_solver->getParameterList().set("graph_fill", dataFile( ( section + "/AztecOO/graph_fill" ).data(), 6. ));

    M_solver->getParameterList().set("ilut_fill", dataFile( ( section + "/AztecOO/ilut_fill" ).data(), 4. ));

    //M_solver->getParameterList().set("init_guess", dataFile( ( section + "/AztecOO/init_guess" ).data(), AZ_NOT_ZERO ));

    //M_solver->getParameterList().set("keep_kvecs", dataFile( ( section + "/AztecOO/keep_kvecs" ).data(), 0 ));

    //M_solver->getParameterList().set("apply_kvecs", dataFile( ( section + "/AztecOO/apply_kvecs" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("orth_kvecs", dataFile( ( section + "/AztecOO/orth_kvecs" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("ignore_scaling", dataFile( ( section + "/AztecOO/ignore_scaling" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("check_update_size", dataFile( ( section + "/AztecOO/check_update_size" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("omega", dataFile( ( section + "/AztecOO/omega" ).data(), 1. ));

    //M_solver->getParameterList().set("update_reduction", dataFile( ( section + "/AztecOO/update_reduction" ).data(), 10e10 ));


    // SET PARAMETERS
    M_solver->setParameters( true );

    // DISPLAY LIST
    if ( dataFile( (section + "/displayList").data(), false ) )
        M_solver->getParameterList().print(std::cout);
}

int
AztecOOPreconditioner::buildPreconditioner( operator_type& Operator )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::buildPreconditioner( Operator ) \n";
#endif

    if ( this->M_preconditionerCreated )
        precReset();

    M_solver->getSolver().SetPrecMatrix( Operator->getMatrixPtr().get() );

    M_solver->getSolver().SetAztecOption( AZ_pre_calc, AZ_calc );
    M_solver->getSolver().SetAztecOption( AZ_keep_info, 1 );

    Real estimateConditionNumber;
    M_solver->getSolver().ConstructPreconditioner( estimateConditionNumber );

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

Real
AztecOOPreconditioner::Condest()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::Condest() \n";
#endif

    return M_solver->getSolver().Condest();
}

EpetraPreconditioner::prec_raw_type*
AztecOOPreconditioner::getPrec()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::getPrec() \n";
#endif

    if ( this->M_preconditionerCreated )
        return M_solver->getSolver().GetPrecMatrix();

    return 0;
}

EpetraPreconditioner::prec_type
AztecOOPreconditioner::getPrecPtr()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::getPrec() \n";
#endif

    boost::shared_ptr<Epetra_RowMatrix> prec;
    if ( this->M_preconditionerCreated )
    {
        prec.reset(M_solver->getSolver().GetPrecMatrix());
        return  prec;
    }
    return EpetraPreconditioner::prec_type();
}

void
AztecOOPreconditioner::precReset()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::precReset() \n";
#endif

    M_solver->getSolver().SetAztecOption( AZ_keep_info, 0 );
    M_solver->getSolver().SetAztecOption( AZ_pre_calc, AZ_calc );

    // Perform one "fake" iteration to delete the preconditioner
    int AZoutputOption = M_solver->getSolver().GetAztecOption( AZ_output );
    M_solver->getSolver().SetAztecOption( AZ_output, AZ_none );
    //M_solver->getSolver().GetRHS()->PutScalar( 1.0 );
    //M_solver->getSolver().GetLHS()->PutScalar( 0.0 );
    M_solver->getSolver().Iterate( 0, 1.e14 );
    M_solver->getSolver().SetAztecOption( AZ_output, AZoutputOption );

    //M_solver->getSolver().DestroyPreconditioner();

    this->M_preconditionerCreated = false;
}

} // namespace LifeV
