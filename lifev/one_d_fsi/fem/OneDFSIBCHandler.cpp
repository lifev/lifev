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
 *  @file
 *  @brief File containing a class for the boundary conditions handling of the 1D model.
 *
 *  @version 1.0
 *  @date 01-28-2006
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *
 *  @version 2.0
 *  @date 20-04-2010
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *
 *  @maintainer Cristiano Malossi <cristiano.malossi@epfl.ch>
 */

#include <lifev/one_d_fsi/fem/OneDFSIBCHandler.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDFSIBCHandler::OneDFSIBCHandler() :
    M_boundary          (),
    M_boundarySet       (),
    M_defaultFunctions  ()
{

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::OneDFSIModel_BCHandler] Creating OneDFSIModel_BC classes.\n";
#endif

    M_boundary[ OneDFSI::left ].reset (  new bc_Type ( OneDFSI::left ) );
    M_boundary[ OneDFSI::right ].reset ( new bc_Type ( OneDFSI::right ) );

    M_boundarySet[ OneDFSI::left ].insert (  std::make_pair ( OneDFSI::first,  false ) );
    M_boundarySet[ OneDFSI::left ].insert (  std::make_pair ( OneDFSI::second, false ) );
    M_boundarySet[ OneDFSI::right ].insert ( std::make_pair ( OneDFSI::first,  false ) );
    M_boundarySet[ OneDFSI::right ].insert ( std::make_pair ( OneDFSI::second, false ) );
}

OneDFSIBCHandler::OneDFSIBCHandler ( const OneDFSIBCHandler& bcHandler ) :
    M_boundary          (),
    M_boundarySet       ( bcHandler.M_boundarySet ),
    M_defaultFunctions  ()
{
    M_boundary[ OneDFSI::left ].reset (  new bc_Type ( *bcHandler.M_boundary.find ( OneDFSI::left )->second ) );
    M_boundary[ OneDFSI::right ].reset ( new bc_Type ( *bcHandler.M_boundary.find ( OneDFSI::right )->second ) );

    for ( std::vector < bcFunctionSolverDefinedPtr_Type >::const_iterator i = bcHandler.M_defaultFunctions.begin();
            i != bcHandler.M_defaultFunctions.end() ; ++i )
    {
        bcFunctionSolverDefinedPtr_Type BCDefaultFunction ( new bcFunctionSolverDefined_Type ( *i->get() ) );
        M_defaultFunctions.push_back ( BCDefaultFunction );
    }

    // NOTE: The copy constructor is not working correctly. All the members of the class are true copy, but
    // the BCFunctions inside M_boundary are still pointing to the original M_defaultFunction, instead
    // of to the copy (which remain unused!). This is because the link between M_boundary and M_defaultFunction
    // is provided by std::bind and, for now, we have no solution for this.
    std::cerr << "!!! WARNING: COPY CONSTRUCTOR DOES NOT CREATE A TRUE COPY !!!" << std::endl;
    std::exit ( EXIT_FAILURE );
}

// ===================================================
// Methods
// ===================================================
void
OneDFSIBCHandler::applyBC ( const Real& time, const Real& timeStep, const solution_Type& solution,
                            const fluxPtr_Type& fluxPtr, vectorPtrContainer_Type& rhs )
{
    M_boundary[ OneDFSI::left  ]->applyBC ( time, timeStep, solution, fluxPtr, rhs );
    M_boundary[ OneDFSI::right ]->applyBC ( time, timeStep, solution, fluxPtr, rhs );
}

void
OneDFSIBCHandler::applyViscoelasticBC ( const fluxPtr_Type& fluxPtr, matrix_Type& matrix, vector_Type& rhs )
{
    M_boundary[ OneDFSI::left  ]->applyViscoelasticBC ( fluxPtr, matrix, rhs );
    M_boundary[ OneDFSI::right ]->applyViscoelasticBC ( fluxPtr, matrix, rhs );
}

// ===================================================
// Set Methods
// ===================================================
void
OneDFSIBCHandler::setBC ( const bcSide_Type& bcSide, const bcLine_Type& bcLine,
                          const bcType_Type& bcType, const bcFunction_Type& bcFunction )
{
    M_boundarySet[bcSide][bcLine] = true;
    M_boundary[bcSide]->setType ( bcLine, bcType );
    M_boundary[bcSide]->setBCFunction ( bcLine, bcFunction );

#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::setBC] imposing function at "
                         << bcSide << " boundary (" << bcLine << " bcLine), variable " << bcType << ".\n";
#endif
}

void
OneDFSIBCHandler::setDefaultBC()
{
#ifdef HAVE_LIFEV_DEBUG
    debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::OneDFSIModel_BCHandler] Set Default BC ... \n";
#endif

    if ( !M_boundarySet[ OneDFSI::left ][OneDFSI::first] )
    {
        //bcFunctionSolverDefinedPtr_Type bcFunction( new OneDFSIFunctionSolverDefined( OneDFSI::OneD_left, OneDFSI::OneD_W1 ) );
        bcFunctionSolverDefinedPtr_Type bcDefaultFunction ( new OneDFSIFunctionSolverDefinedRiemann ( OneDFSI::left, OneDFSI::W1 ) );
        M_defaultFunctions.push_back ( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction ( std::bind ( &OneDFSIFunctionSolverDefinedRiemann::operator(),
                                               dynamic_cast<OneDFSIFunctionSolverDefinedRiemann*> ( & ( *M_defaultFunctions.back() ) ), std::placeholders::_1, std::placeholders::_2 ) );

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::setDefaultBC] left-first-W1 Invoking setBC.\n";
#endif
        setBC ( OneDFSI::left, OneDFSI::first, OneDFSI::W1, bcFunction );
    }

    if ( !M_boundarySet[ OneDFSI::left ][OneDFSI::second] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDFSIFunctionSolverDefinedCompatibility( OneDFSI::OneD_left, OneDFSI::OneD_W2 ) );
        bcFunctionSolverDefinedPtr_Type bcDefaultFunction ( new OneDFSIFunctionSolverDefinedCompatibility ( OneDFSI::left, OneDFSI::W2 ) );
        M_defaultFunctions.push_back ( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction ( std::bind ( &OneDFSIFunctionSolverDefinedCompatibility::operator(),
                                               dynamic_cast<OneDFSIFunctionSolverDefinedCompatibility*> ( & ( *M_defaultFunctions.back() ) ), std::placeholders::_1, std::placeholders::_2 ) );

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::setDefaultBC] left-second-W2 Invoking setBC.\n";
#endif
        setBC ( OneDFSI::left, OneDFSI::second, OneDFSI::W2, bcFunction );
    }

    if ( !M_boundarySet[ OneDFSI::right ][ OneDFSI::first ] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDFSIFunctionSolverDefinedRiemann( OneDFSI::OneD_right, OneDFSI::OneD_W2 ) );
        bcFunctionSolverDefinedPtr_Type bcDefaultFunction ( new OneDFSIFunctionSolverDefinedRiemann ( OneDFSI::right, OneDFSI::W2 ) );
        M_defaultFunctions.push_back ( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction ( std::bind ( &OneDFSIFunctionSolverDefinedRiemann::operator(),
                                               dynamic_cast<OneDFSIFunctionSolverDefinedRiemann*> ( & ( *M_defaultFunctions.back() ) ), std::placeholders::_1, std::placeholders::_2 ) );

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::setDefaultBC] right-first-W2 Invoking setBC.\n";
#endif
        setBC ( OneDFSI::right, OneDFSI::first, OneDFSI::W2, bcFunction );
    }

    if ( !M_boundarySet[ OneDFSI::right ][ OneDFSI::second ] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDFSIFunctionSolverDefinedCompatibility( OneDFSI::OneD_right, OneDFSI::OneD_W1 ) );
        bcFunctionSolverDefinedPtr_Type bcDefaultFunction ( new OneDFSIFunctionSolverDefinedCompatibility ( OneDFSI::right, OneDFSI::W1 ) );
        M_defaultFunctions.push_back ( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction ( std::bind ( &OneDFSIFunctionSolverDefinedCompatibility::operator(),
                                               dynamic_cast<OneDFSIFunctionSolverDefinedCompatibility*> ( & ( *M_defaultFunctions.back() ) ), std::placeholders::_1, std::placeholders::_2 ) );

#ifdef HAVE_LIFEV_DEBUG
        debugStream ( 6311 ) << "[OneDFSIModel_BCHandler::setDefaultBC] right-second-W1 Invoking setBC.\n";
#endif
        setBC ( OneDFSI::right, OneDFSI::second, OneDFSI::W1, bcFunction );
    }
}

void
OneDFSIBCHandler::setFluxSource ( const fluxPtr_Type& fluxPtr, const sourcePtr_Type& sourcePtr )
{
    for ( std::vector < bcFunctionSolverDefinedPtr_Type >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
    {
        ( *i )->setFluxSource ( fluxPtr, sourcePtr );
    }
}

void
OneDFSIBCHandler::setSolution ( const solutionPtr_Type& solution )
{
    for ( std::vector < bcFunctionSolverDefinedPtr_Type >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
    {
        ( *i )->setSolution ( solution );
    }
}

}
