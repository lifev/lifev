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

#include "PreconditionerAztecOO.hpp"
#include <lifev/core/LifeV.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerAztecOO::PreconditionerAztecOO() :
    super   (),
    M_solver()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::PreconditionerAztecOO() \n";
#endif

}


// ===================================================
// Methods
// ===================================================
Int
PreconditionerAztecOO::buildPreconditioner ( operator_type& matrix )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::buildPreconditioner( Operator ) \n";
#endif

    if ( this->M_preconditionerCreated )
    {
        resetPreconditioner();
    }

    M_solver->solver().SetPrecMatrix ( matrix->matrixPtr().get() );

    M_solver->solver().SetAztecOption ( AZ_pre_calc, AZ_calc );
    M_solver->solver().SetAztecOption ( AZ_keep_info, 1 );

    Real estimateConditionNumber;
    M_solver->solver().ConstructPreconditioner ( estimateConditionNumber );

    this->M_preconditionerCreated = true;

    return ( EXIT_SUCCESS );
}

void
PreconditionerAztecOO::resetPreconditioner()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::precReset() \n";
#endif

    M_solver->solver().SetAztecOption ( AZ_keep_info, 0 );
    M_solver->solver().SetAztecOption ( AZ_pre_calc, AZ_calc );

    // Perform one "fake" iteration to delete the preconditioner
    Int AZoutputOption = M_solver->solver().GetAztecOption ( AZ_output );
    M_solver->solver().SetAztecOption ( AZ_output, AZ_none );
    M_solver->solver().Iterate ( 0, 1.e14 );
    M_solver->solver().SetAztecOption ( AZ_output, AZoutputOption );

    this->M_preconditionerCreated = false;
}

void
PreconditionerAztecOO::createParametersList ( list_Type& /*list*/, const GetPot& dataFile, const std::string& section, const std::string& subSection )
{
    // Preconditioner
    M_solver->getParametersList().set ( "precond", dataFile ( ( section + "/" + subSection + "/precond" ).data(), "dom_decomp" ) );

    // Compute the preconditioner
    M_solver->getParametersList().set ( "pre_calc", dataFile ( ( section + "/" + subSection + "/pre_calc" ).data(), "calc" ) );

    // Reordering
    M_solver->getParametersList().set ( "reorder", dataFile ( ( section + "/" + subSection + "/reorder" ).data(), 1 ) );

    // Keep the information
    M_solver->getParametersList().set ( "keep_info", dataFile ( ( section + "/" + subSection + "/keep_info" ).data(), 1 ) );

    // Polynomial order when using Jacobi and GS preconditioners
    //M_solver->getParametersList().set( "poly_ord", dataFile( ( section + "/" + subSection + "/poly_ord" ).data(), 3 ) );

    // Overlap level
    M_solver->getParametersList().set ( "overlap", dataFile ( ( section + "/" + subSection + "/overlap" ).data(), AZ_none ) );

    // Overlap type
    M_solver->getParametersList().set ( "type_overlap", dataFile ( ( section + "/" + subSection + "/type_overlap" ).data(), AZ_standard ) );


    // SUBDOMAIN SOLVER

    M_solver->getParametersList().set ( "subdomain_solve", dataFile ( ( section + "/" + subSection + "/subdomain_solve" ).data(), "ILUT" ) );

    M_solver->getParametersList().set ( "drop", dataFile ( ( section + "/" + subSection + "/drop" ).data(), 1.e-5 ) );

    //M_solver->getParametersList().set( "graph_fill", dataFile( ( section + "/" + subSection + "/graph_fill" ).data(), 6. ) );

    M_solver->getParametersList().set ( "ilut_fill", dataFile ( ( section + "/" + subSection + "/ilut_fill" ).data(), 4. ) );

    //M_solver->getParametersList().set( "init_guess", dataFile( ( section + "/" + subSection + "/init_guess" ).data(), AZ_NOT_ZERO ) );

    //M_solver->getParametersList().set( "keep_kvecs", dataFile( ( section + "/" + subSection + "/keep_kvecs" ).data(), 0 ) );

    //M_solver->getParametersList().set( "apply_kvecs", dataFile( ( section + "/" + subSection + "/apply_kvecs" ).data(), AZ_FALSE ) );

    //M_solver->getParametersList().set( "orth_kvecs", dataFile( ( section + "/" + subSection + "/orth_kvecs" ).data(), AZ_FALSE ) );

    //M_solver->getParametersList().set( "ignore_scaling", dataFile( ( section + "/" + subSection + "/ignore_scaling" ).data(), AZ_FALSE ) );

    //M_solver->getParametersList().set( "check_update_size", dataFile( ( section + "/" + subSection + "/check_update_size" ).data(), AZ_FALSE ) );

    //M_solver->getParametersList().set( "omega", dataFile( ( section + "/" + subSection + "/omega" ).data(), 1. ) );

    //M_solver->getParametersList().set( "update_reduction", dataFile( ( section + "/" + subSection + "/update_reduction" ).data(), 10e10 ) );

    // SET PARAMETERS
    M_solver->setParameters ( true );

    // DISPLAY LIST
    if ( dataFile ( (section + "/displayList").data(), false ) )
    {
        M_solver->getParametersList().print ( std::cout );
    }
}

void
PreconditionerAztecOO::showMe ( std::ostream& output ) const
{
    output << "showMe must be implemented for the PreconditionerAztecOO class" << std::endl;
}

Real
PreconditionerAztecOO::condest()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::Condest() \n";
#endif

    return M_solver->solver().Condest();
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerAztecOO::setDataFromGetPot ( const GetPot&      dataFile,
                                           const std::string& section )
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::setDataFromGetPot(dataFile, section) \n";
#endif

    createParametersList ( M_list, dataFile, section, "AztecOO" );
}


// ===================================================
// Get Methods
// ===================================================
Preconditioner::prec_raw_type*
PreconditionerAztecOO::preconditioner()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::getPrec() \n";
#endif

    if ( this->M_preconditionerCreated )
    {
        return M_solver->solver().GetPrecMatrix();
    }

    return 0;
}

Preconditioner::prec_type
PreconditionerAztecOO::preconditionerPtr()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 7100 ) << "PreconditionerAztecOO::getPrec() \n";
#endif

    boost::shared_ptr<Epetra_RowMatrix> prec;
    if ( this->M_preconditionerCreated )
    {
        prec.reset (M_solver->solver().GetPrecMatrix() );
        return  prec;
    }
    return Preconditioner::prec_type();
}

} // namespace LifeV
