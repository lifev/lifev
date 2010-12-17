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

#include <lifemc/lifefem/OneDimensionalModel_BCHandler.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler() :
        M_boundary          (),
        M_boundarySet       (),
        M_defaultFunctions  ()
{

#ifdef HAVE_LIFEV_DEBUG
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Creating OneDimensionalModel_BC classes.\n";
#endif

    M_boundary[ OneDimensional::left ].reset(  new bc_Type( OneDimensional::left ) );
    M_boundary[ OneDimensional::right ].reset( new bc_Type( OneDimensional::right ) );

    M_boundarySet[ OneDimensional::left ].insert(  std::make_pair( OneDimensional::first,  false ) );
    M_boundarySet[ OneDimensional::left ].insert(  std::make_pair( OneDimensional::second, false ) );
    M_boundarySet[ OneDimensional::right ].insert( std::make_pair( OneDimensional::first,  false ) );
    M_boundarySet[ OneDimensional::right ].insert( std::make_pair( OneDimensional::second, false ) );
}

OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler( const OneDimensionalModel_BCHandler& bcHandler ) :
        M_boundary          (),
        M_boundarySet       ( bcHandler.M_boundarySet ),
        M_defaultFunctions  ()
{
    M_boundary[ OneDimensional::left ].reset(  new bc_Type( *bcHandler.M_boundary.find( OneDimensional::left )->second ) );
    M_boundary[ OneDimensional::right ].reset( new bc_Type( *bcHandler.M_boundary.find( OneDimensional::right )->second ) );

    for ( std::vector < bcFunctionDefaultPtr_Type >::const_iterator i = bcHandler.M_defaultFunctions.begin();
            i != bcHandler.M_defaultFunctions.end() ; ++i )
    {
        bcFunctionDefaultPtr_Type BCDefaultFunction( new bcFunctionDefault_Type( *i->get() ) );
        M_defaultFunctions.push_back( BCDefaultFunction );
    }

    // NOTE: The copy constructor is not working correctly. All the members of the class are true copy, but
    // the BCFunctions inside M_boundary are still pointing to the original M_defaultFunction, instead
    // of to the copy (which remain unused!). This is because the link between M_boundary and M_defaultFunction
    // is provided by boost::bind and, for now, we have no solution for this.
    std::cout << "!!! WARNING: COPY CONSTRUCTOR DOES NOT CREATE A TRUE COPY !!!" << std::endl;
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_BCHandler::applyBC( const Real& time, const Real& timeStep, const solution_Type& solution, const fluxPtr_Type& flux,
                                        container2D_Type& leftBC, container2D_Type& rightBC )
{
#ifdef HAVE_LIFEV_DEBUG
    ASSERT_PRE( leftBC.size() == 2 && rightBC.size() == 2, "applyBC works only for 2D vectors" );
#endif

    M_boundary[ OneDimensional::left  ]->applyBC( time, timeStep, solution, flux, leftBC  );
    M_boundary[ OneDimensional::right ]->applyBC( time, timeStep, solution, flux, rightBC );

#ifdef HAVE_LIFEV_DEBUG
    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at left "
    << " imposing [ A, Q ] = [ " << leftBC[0]
    << ", " << leftBC[1] << " ]\n";
    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at right "
    << " imposing [ A, Q ] = [ " << rightBC[0]
    << ", " << rightBC[1] << " ]\n";
#endif
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCHandler::setBC( const bcSide_Type& bcSide, const bcLine_Type& bcLine,
                                      const bcTypeOneD_Type& bcType, const bcFunction_Type& BCfunction )
{
    M_boundarySet[bcSide][bcLine] = true;
    M_boundary[bcSide]->setType( bcLine, bcType );
    M_boundary[bcSide]->setBCFunction( bcLine, BCfunction );

#ifdef HAVE_LIFEV_DEBUG
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setBC] imposing function at "
                  << bcSide << " boundary (" << bcLine << " bcLine), variable " << bcType << ".\n";
#endif
}

void
OneDimensionalModel_BCHandler::setDefaultBC()
{
#ifdef HAVE_LIFEV_DEBUG
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Set Default BC ... \n";
#endif

    if ( !M_boundarySet[ OneDimensional::left ][OneDimensional::first] )
    {
        //bcFunctionDefaultPtr_Type bcFunction( new OneDimensionalModel_BCFunction_Default( OneDimensional::OneD_left, OneDimensional::OneD_W1 ) );
        bcFunctionDefaultPtr_Type bcDefaultFunction( new OneDimensionalModel_BCFunction_Riemann( OneDimensional::left, OneDimensional::W1 ) );
        M_defaultFunctions.push_back( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-first-W1 Invoking setBC.\n";
#endif
        setBC( OneDimensional::left, OneDimensional::first, OneDimensional::W1, bcFunction );
    }

    if ( !M_boundarySet[ OneDimensional::left ][OneDimensional::second] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDimensionalModel_BCFunction_Compatibility( OneDimensional::OneD_left, OneDimensional::OneD_W2 ) );
        bcFunctionDefaultPtr_Type bcDefaultFunction( new OneDimensionalModel_BCFunction_Compatibility( OneDimensional::left, OneDimensional::W2 ) );
        M_defaultFunctions.push_back( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-second-W2 Invoking setBC.\n";
#endif
        setBC( OneDimensional::left, OneDimensional::second, OneDimensional::W2, bcFunction );
    }

    if ( !M_boundarySet[ OneDimensional::right ][ OneDimensional::first ] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDimensionalModel_BCFunction_Riemann( OneDimensional::OneD_right, OneDimensional::OneD_W2 ) );
        bcFunctionDefaultPtr_Type bcDefaultFunction( new OneDimensionalModel_BCFunction_Riemann( OneDimensional::right, OneDimensional::W2 ) );
        M_defaultFunctions.push_back( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-first-W2 Invoking setBC.\n";
#endif
        setBC( OneDimensional::right, OneDimensional::first, OneDimensional::W2, bcFunction );
    }

    if ( !M_boundarySet[ OneDimensional::right ][ OneDimensional::second ] )
    {
        //bcFunctionPtr_Type bcFunction( new OneDimensionalModel_BCFunction_Compatibility( OneDimensional::OneD_right, OneDimensional::OneD_W1 ) );
        bcFunctionDefaultPtr_Type bcDefaultFunction( new OneDimensionalModel_BCFunction_Compatibility( OneDimensional::right, OneDimensional::W1 ) );
        M_defaultFunctions.push_back( bcDefaultFunction );

        //bcFunctionPtr_Type bcFunction ( new bcFunction_Type() );
        bcFunction_Type bcFunction;
        bcFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

#ifdef HAVE_LIFEV_DEBUG
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-second-W1 Invoking setBC.\n";
#endif
        setBC( OneDimensional::right, OneDimensional::second, OneDimensional::W1, bcFunction );
    }
}

void
OneDimensionalModel_BCHandler::setSolution( const solutionPtr_Type& solution )
{
    for ( std::vector < bcFunctionDefaultPtr_Type >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
        ( *i )->setSolution( solution );
}

void
OneDimensionalModel_BCHandler::setFluxSource( const fluxPtr_Type& flux, const sourcePtr_Type& source )
{
    for ( std::vector < bcFunctionDefaultPtr_Type >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
        ( *i )->setFluxSource( flux, source );
}

}
