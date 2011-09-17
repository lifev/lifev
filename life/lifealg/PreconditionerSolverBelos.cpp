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
    @brief Solver Belos Operator

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-09-2011
 */

#include <life/lifealg/PreconditionerSolverBelos.hpp>
#include <life/lifecore/LifeDebug.hpp>
#include <life/lifecore/LifeChrono.hpp>
#include <life/lifefilters/GetPot.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
PreconditionerSolverBelos::PreconditionerSolverBelos( boost::shared_ptr<Epetra_Comm> comm ) :
        Preconditioner( comm )
{

}

PreconditionerSolverBelos::~PreconditionerSolverBelos()
{
    M_prec.reset();
}

// ===================================================
// Methods
// ===================================================

void
PreconditionerSolverBelos::createParametersList( list_Type&         list,
                                                  const GetPot&      dataFile,
                                                  const std::string& section,
                                                  const std::string& subSection )
{
    createSolverAmesosList( list, dataFile, section, subSection );
}

void
PreconditionerSolverBelos::createSolverBelosList( list_Type&         list,
                                                  const GetPot&      dataFile,
                                                  const std::string& section,
                                                  const std::string& subsection )
{
    bool displayList = dataFile( ( section + "/displayList" ).data(), false);

    // Status parameters
    list.set( "OutputLevel",  dataFile( ( section + "/" + subsection + "/outputlevel").data(), 0 ) );
    list.set( "PrintStatus",  dataFile( ( section + "/" + subsection + "/print_status").data(), false ) );
    list.set( "PrintTiming",  dataFile( ( section + "/" + subsection + "/print_timing").data(), false ) );
    list.set( "ComputeVectorNorms", dataFile( ( section + "/" + subsection + "/computevectornorms").data(), false ) );
    list.set( "ComputeTrueResidual", dataFile( ( section + "/" + subsection + "/computeresidual").data(), false ) );

    // Control parameters
    list.set( "AddZeroToDiag",  dataFile( ( section + "/" + subsection + "/addzerotodiag").data(), false ) );
    list.set( "Refactorize", dataFile( ( section + "/" + subsection + "/refactorize").data(), false ) );
    list.set( "RcondThreshold", dataFile( ( section + "/" + subsection + "/rcondthreshold").data(), 1.e-2) );
    list.set( "Redistribute", dataFile( ( section + "/" + subsection + "/redistribute").data(), true ) ); // SuperLU
    list.set( "MaxProcs", dataFile( ( section + "/" + subsection + "/maxprocs").data(), -1) ); // ScalaPack
    list.set( "ScaleMethod", dataFile( ( section + "/" + subsection + "/amesos/scalemethod").data(), 1) );

    // Type of the matrix: symmetric, SDP, general
    list.set( "MatrixProperty", dataFile( ( section + "/" + subsection + "/matrixproperty").data(), "general" ) );

    // Type of the solver
    list.set( "SolverType", dataFile( ( section + "/" + subsection + "/solvertype"  ).data(), "Klu" ) );

    if ( displayList ) list.print( std::cout );
}

Int
PreconditionerSolverBelos::buildPreconditioner( operator_type& matrix )
{
    M_prec.reset( new precOperator_Type( this->M_displayer.comm() ) );
    M_prec->buildPreconditioner( matrix, M_list );

    this->M_preconditionerCreated = true;

    return 0;
}

void
PreconditionerSolverBelos::resetPreconditioner()
{
    M_prec.reset();
    this->M_preconditionerCreated = false;
}

Real
PreconditionerSolverBelos::condest()
{
    return 0.0;
}

void
PreconditionerSolverBelos::showMe( std::ostream& output ) const
{
    M_prec->showMe(output);
}

// ===================================================
// Epetra Operator Interface Methods
// ===================================================
Int
PreconditionerSolverBelos::SetUseTranspose( const bool useTranspose )
{
    return M_prec->SetUseTranspose(useTranspose);
}

bool
PreconditionerSolverBelos::UseTranspose()
{
    return M_prec->UseTranspose();
}

Int
PreconditionerSolverBelos::Apply( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    return M_prec->Apply( X, Y );
}

Int
PreconditionerSolverBelos::ApplyInverse( const Epetra_MultiVector& X, Epetra_MultiVector& Y ) const
{
    if( M_prec )
    {
        M_prec->ApplyInverse( X, Y );
    }
    return 0;
}

const Epetra_Map&
PreconditionerSolverBelos::OperatorRangeMap() const
{
    return M_prec->OperatorRangeMap();
}

const Epetra_Map&
PreconditionerSolverBelos::OperatorDomainMap() const
{
    return M_prec->OperatorDomainMap();
}

// ===================================================
// Set Methods
// ===================================================
void
PreconditionerSolverBelos::setDataFromGetPot ( const GetPot& dataFile, const std::string& section )
{
    createSolverAmesosList( M_list, dataFile, section, "SolverAmesos" );
}

void
PreconditionerSolverBelos::setSolver( SolverAztecOO& /*solver*/ )
{

}

// ===================================================
// Get Methods
// ===================================================
//! Return true if the preconditioner is set
bool
PreconditionerSolverBelos::isPreconditionerSet() const
{
    return M_prec;
}

PreconditionerSolverBelos::prec_raw_type*
PreconditionerSolverBelos::preconditioner()
{
    return M_prec.get();
}

PreconditionerSolverBelos::prec_type
PreconditionerSolverBelos::preconditionerPtr()
{
    return M_prec;
}

std::string
PreconditionerSolverBelos::preconditionerType()
{
    return "SolverBelos preconditioner";
}

// ===================================================
// Private Methods
// ===================================================


} // namespace LifeV

