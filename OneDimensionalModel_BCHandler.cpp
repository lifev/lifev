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
 *  @author Lucia Mirabella <lucia@mathcs.emory.edu>
 *  @author Tiziano Passerini <tiziano@mathcs.emory.edu>
 *  @date 01-28-2006
 *
 *  @version 2.0
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 20-04-2010
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
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Creating OneDimensionalModel_BC classes.\n";

    M_boundary[ OneD_left ].reset(  new BC_Type( OneD_left ) );
    M_boundary[ OneD_right ].reset( new BC_Type( OneD_right ) );

    M_boundarySet[ OneD_left ].insert(  std::make_pair( OneD_first,  false ) );
    M_boundarySet[ OneD_left ].insert(  std::make_pair( OneD_second, false ) );
    M_boundarySet[ OneD_right ].insert( std::make_pair( OneD_first,  false ) );
    M_boundarySet[ OneD_right ].insert( std::make_pair( OneD_second, false ) );
}

OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler( const OneDimensionalModel_BCHandler& BCH ) :
        M_boundary          (),
        M_boundarySet       ( BCH.M_boundarySet ),
        M_defaultFunctions  ()
{
    M_boundary[ OneD_left ].reset(  new BC_Type( *BCH.M_boundary.find( OneD_left )->second ) );
    M_boundary[ OneD_right ].reset( new BC_Type( *BCH.M_boundary.find( OneD_right )->second ) );

    for ( std::vector < BCFunction_Default_PtrType >::const_iterator i = BCH.M_defaultFunctions.begin();
            i != BCH.M_defaultFunctions.end() ; ++i )
    {
        BCFunction_Default_PtrType BCDefaultFunction( new BCFunction_Default_Type( *i->get() ) );
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
OneDimensionalModel_BCHandler::applyBC( const Real&              time,
                                        const Real&              timeStep,
                                        const Solution_Type&     solution,
                                        const Flux_PtrType&      flux,
                                        container2D_Type&  leftBC,
                                        container2D_Type&  rightBC )
{
    ASSERT_PRE( leftBC.size() == 2 && rightBC.size() == 2, "applyBC works only for 2D vectors" );

    M_boundary[ OneD_left  ]->applyBC( time, timeStep, solution, flux, leftBC  );
    M_boundary[ OneD_right ]->applyBC( time, timeStep, solution, flux, rightBC );

    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at left "
    << " imposing [ A, Q ] = [ " << leftBC[0]
    << ", " << leftBC[1] << " ]\n";
    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at right "
    << " imposing [ A, Q ] = [ " << rightBC[0]
    << ", " << rightBC[1] << " ]\n";
}

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCHandler::setBC( const OneD_BCSide&     side,
                                      const OneD_BCLine&     line,
                                      const OneD_BC&         bcType,
                                      const BCFunction_Type& BCfunction )
{
    M_boundarySet[side][line] = true;
    M_boundary[side]->setType( line, bcType );
    M_boundary[side]->setBCFunction( line, BCfunction );
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setBC] imposing function at "
                  << side << " boundary ("
                  << line << " line), variable "
                  << bcType << ".\n";
}

void
OneDimensionalModel_BCHandler::setDefaultBC()
{
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Set Default BC ... \n";

    if ( !M_boundarySet[ OneD_left ][OneD_first] )
    {
        //BCFunction_Default_PtrType BCFunction ( new OneDimensionalModel_BCFunction_Default( OneD_left, OneD_W1 ) );
        BCFunction_Default_PtrType BCDefaultFunction ( new OneDimensionalModel_BCFunction_Riemann( OneD_left, OneD_W1 ) );
        M_defaultFunctions.push_back( BCDefaultFunction );

        //BCFunction_PtrType BCFunction ( new BCFunction_Type() );
        BCFunction_Type BCFunction;
        BCFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-first-W1 Invoking setBC.\n";
        setBC( OneD_left, OneD_first, OneD_W1, BCFunction );
    }

    if ( !M_boundarySet[ OneD_left ][OneD_second] )
    {
        //BCFunction_PtrType BCFunction ( new OneDimensionalModel_BCFunction_Compatibility( OneD_left, OneD_W2 ) );
        BCFunction_Default_PtrType BCDefaultFunction ( new OneDimensionalModel_BCFunction_Compatibility( OneD_left, OneD_W2 ) );
        M_defaultFunctions.push_back( BCDefaultFunction );

        //BCFunction_PtrType BCFunction ( new BCFunction_Type() );
        BCFunction_Type BCFunction;
        BCFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );


        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-second-W2 Invoking setBC.\n";
        setBC( OneD_left, OneD_second, OneD_W2, BCFunction );
    }

    if ( !M_boundarySet[ OneD_right ][ OneD_first ] )
    {
        //BCFunction_PtrType BCFunction ( new OneDimensionalModel_BCFunction_Riemann( OneD_right, OneD_W2 ) );
        BCFunction_Default_PtrType BCDefaultFunction ( new OneDimensionalModel_BCFunction_Riemann( OneD_right, OneD_W2 ) );
        M_defaultFunctions.push_back( BCDefaultFunction );

        //BCFunction_PtrType BCFunction ( new BCFunction_Type() );
        BCFunction_Type BCFunction;
        BCFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Riemann::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Riemann *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-first-W2 Invoking setBC.\n";
        setBC( OneD_right, OneD_first, OneD_W2, BCFunction );
    }

    if ( !M_boundarySet[ OneD_right ][ OneD_second ] )
    {
        //BCFunction_PtrType BCFunction ( new OneDimensionalModel_BCFunction_Compatibility( OneD_right, OneD_W1 ) );
        BCFunction_Default_PtrType BCDefaultFunction ( new OneDimensionalModel_BCFunction_Compatibility( OneD_right, OneD_W1 ) );
        M_defaultFunctions.push_back( BCDefaultFunction );

        //BCFunction_PtrType BCFunction ( new BCFunction_Type() );
        BCFunction_Type BCFunction;
        BCFunction.setFunction( boost::bind( &OneDimensionalModel_BCFunction_Compatibility::operator(),
                                             dynamic_cast<OneDimensionalModel_BCFunction_Compatibility *> ( &( *M_defaultFunctions.back() ) ), _1, _2 ) );

        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-second-W1 Invoking setBC.\n";
        setBC( OneD_right, OneD_second, OneD_W1, BCFunction );
    }
}


void
OneDimensionalModel_BCHandler::setSolution( const Solution_PtrType& solution )
{
    for ( std::vector < BCFunction_Default_PtrType >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
        ( *i )->setSolution( solution );
}

void
OneDimensionalModel_BCHandler::setFluxSource( const Flux_PtrType& flux, const Source_PtrType& source )
{
    for ( std::vector < BCFunction_Default_PtrType >::const_iterator i = M_defaultFunctions.begin() ; i < M_defaultFunctions.end() ; ++i )
        ( *i )->setFluxSource( flux, source );
}

void
OneDimensionalModel_BCHandler::setInternalNode( const OneD_BCSide& side )
{
    M_boundary[ side ]->setInternalFlag( true );
}


}
