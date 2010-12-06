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
    @file
    @brief AztecOO preconditioner

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-11-2009
 */

#include "AztecOOPreconditioner.hpp"

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
AztecOOPreconditioner::AztecOOPreconditioner():
        super                   ( ),
        M_solver                ( )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::AztecOOPreconditioner() \n";
#endif

}


// ===================================================
// Methods
// ===================================================
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

void
AztecOOPreconditioner::createList( list_Type& /*list*/, const GetPot& dataFile, const std::string& section, const std::string& subSection )
{
    // Preconditioner
    M_solver->getParameterList().set("precond", dataFile( ( section + "/" + subSection + "/precond" ).data(), "dom_decomp" ));

    // Compute the preconditioner
    M_solver->getParameterList().set("pre_calc", dataFile( ( section + "/" + subSection + "/pre_calc" ).data(), "calc" ));

    // Reordering
    M_solver->getParameterList().set("reorder", dataFile( ( section + "/" + subSection + "/reorder" ).data(), 1 ));

    // Keep the information
    M_solver->getParameterList().set("keep_info", dataFile( ( section + "/" + subSection + "/keep_info" ).data(), 1 ));

    // Polynomial order when using Jacobi and GS preconditioners
    //M_solver->getParameterList().set("poly_ord", dataFile( ( section + "/" + subSection + "/poly_ord" ).data(), 3 ));

    // Overlap level
    M_solver->getParameterList().set("overlap", dataFile( ( section + "/" + subSection + "/overlap" ).data(), AZ_none ));

    // Overlap type
    M_solver->getParameterList().set("type_overlap", dataFile( ( section + "/" + subSection + "/type_overlap" ).data(), AZ_standard ));


    // SUBDOMAIN SOLVER

    M_solver->getParameterList().set("subdomain_solve", dataFile( ( section + "/" + subSection + "/subdomain_solve" ).data(), "ILUT" ));

    M_solver->getParameterList().set("drop", dataFile( ( section + "/" + subSection + "/drop" ).data(), 1.e-5 ));

    //M_solver->getParameterList().set("graph_fill", dataFile( ( section + "/" + subSection + "/graph_fill" ).data(), 6. ));

    M_solver->getParameterList().set("ilut_fill", dataFile( ( section + "/" + subSection + "/ilut_fill" ).data(), 4. ));

    //M_solver->getParameterList().set("init_guess", dataFile( ( section + "/" + subSection + "/init_guess" ).data(), AZ_NOT_ZERO ));

    //M_solver->getParameterList().set("keep_kvecs", dataFile( ( section + "/" + subSection + "/keep_kvecs" ).data(), 0 ));

    //M_solver->getParameterList().set("apply_kvecs", dataFile( ( section + "/" + subSection + "/apply_kvecs" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("orth_kvecs", dataFile( ( section + "/" + subSection + "/orth_kvecs" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("ignore_scaling", dataFile( ( section + "/" + subSection + "/ignore_scaling" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("check_update_size", dataFile( ( section + "/" + subSection + "/check_update_size" ).data(), AZ_FALSE ));

    //M_solver->getParameterList().set("omega", dataFile( ( section + "/" + subSection + "/omega" ).data(), 1. ));

    //M_solver->getParameterList().set("update_reduction", dataFile( ( section + "/" + subSection + "/update_reduction" ).data(), 10e10 ));

    // SET PARAMETERS
    M_solver->setParameters( true );

    // DISPLAY LIST
    if ( dataFile( (section + "/displayList").data(), false ) )
        M_solver->getParameterList().print(std::cout);
}

// ===================================================
// Set Methods
// ===================================================
void
AztecOOPreconditioner::setDataFromGetPot( const GetPot&      dataFile,
                                          const std::string& section )
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 7100 ) << "AztecOOPreconditioner::setDataFromGetPot(dataFile, section) \n";
#endif

    createList( M_List, dataFile, section, "AztecOO" );
}


// ===================================================
// Get Methods
// ===================================================
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

} // namespace LifeV
