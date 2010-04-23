//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

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

namespace LifeV {

// ===================================================
// Constructors & Destructor
// ===================================================
OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler( const std::vector<Vector_Type>& U_thistime,
                                                              const Flux_PtrType              fluxFun,
                                                              const Real&                     dimDof ) :
    M_U_thistime ( U_thistime ),
    M_fluxFun    ( fluxFun )
{
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Creating OneDimensionalModel_BC classes.\n";

    M_boundary["left"].reset( new OneDimensionalModel_BC( U_thistime, /*W_thistime,*/ fluxFun, dimDof, "left" ) );
    M_boundary["right"].reset( new OneDimensionalModel_BC( U_thistime, /*W_thistime,*/ fluxFun, dimDof, "right" ) );

    M_boundarybool["left"].insert( make_pair("first", false));
    M_boundarybool["left"].insert( make_pair("second", false));
    M_boundarybool["right"].insert( make_pair("first", false));
    M_boundarybool["right"].insert( make_pair("second", false));

    M_OneDimensionalModel_BCHandlerMapStringValues["left"]   = OneDBCLeftBoundary;
    M_OneDimensionalModel_BCHandlerMapStringValues["right"]  = OneDBCRightBoundary;
    M_OneDimensionalModel_BCHandlerMapStringValues["first"]  = OneDBCFirstRHS;
    M_OneDimensionalModel_BCHandlerMapStringValues["second"] = OneDBCSecondRHS;
}

// ===================================================
// Methods
// ===================================================
void
OneDimensionalModel_BCHandler::applyBC(const Real&   time_val,
                                             Vec2D&  left_BC_dir,
                                             Vec2D&  right_BC_dir)
{

    ASSERT_PRE( left_BC_dir.size() == 2 && right_BC_dir.size() == 2,
                "applyBC works only for 2D vectors");

    M_boundary["left" ]->applyBC( time_val, left_BC_dir  );
    M_boundary["right"]->applyBC( time_val, right_BC_dir );

    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at left "
                << " imposing [ A, Q ] = [ " << left_BC_dir[0]
                << ", " << left_BC_dir[1] << " ]\n";
    Debug(6311) << "[OneDimensionalModel_BCHandler::applyBC] at right "
                << " imposing [ A, Q ] = [ " << right_BC_dir[0]
                << ", " << right_BC_dir[1] << " ]\n";
};

// ===================================================
// Set Methods
// ===================================================
void
OneDimensionalModel_BCHandler::setDefaultBC( const FESpace_Type&  fespace,
                                             const Source_PtrType sourceFun,
                                             const Real&          dt )
{
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::OneDimensionalModel_BCHandler] Set Default BC ... \n";
    std::string border;
    std::string var;

    if(!M_boundarybool["left"].operator[]("first"))
    {
        border = "left";
        var    = "W1";
        OneDBCFunction_PtrType point ( new Riemann( fespace,
                                                    M_fluxFun,
                                                    M_U_thistime,
                                                    border,
                                                    var ) );
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-first-W1 Invoking setBC.\n";
        setBC(point, "left", "first", "W1");
    }
    if(!M_boundarybool["left"].operator[]("second"))
    {
        border = "left";
        var    = "W2";
        OneDBCFunction_PtrType point ( new Compatibility( fespace,
                                                          M_fluxFun,
                                                          sourceFun,
                                                          M_U_thistime,
                                                          dt,
                                                          border,
                                                          var ) );
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] left-second-W2 Invoking setBC.\n";
        setBC(point, "left", "second", "W2");
    }
    if(!M_boundarybool["right"]["first"])
    {

        border = "right";
        var    = "W2";

        OneDBCFunction_PtrType point ( new Riemann( fespace,
                                                    M_fluxFun,
                                                    M_U_thistime,
                                                    border, var
                                                    ) );
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-first-W2 Invoking setBC.\n";
        setBC(point, "right", "first", "W2");
    }
    if(!M_boundarybool["right"]["second"])
    {
        border = "right";
        var    = "W1";

        OneDBCFunction_PtrType point ( new Compatibility( fespace,
                                                          M_fluxFun,
                                                          sourceFun,
                                                          M_U_thistime,
                                                          dt,
                                                          border,
                                                          var ) );
        Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setDefaultBC] right-second-W1 Invoking setBC.\n";
        setBC(point, "right", "second", "W1");
    }
}

void
OneDimensionalModel_BCHandler::setBC( const OneDBCFunction_PtrType& funptr,
                                      const std::string& border,
                                      const std::string& line,
                                      const std::string& var )
{
    M_boundarybool[border][line] = true;
    M_boundary[border]->variable(line)=var;
    M_boundary[border]->rhs(line) = funptr;
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setBC] imposing function at "
                  << border << " boundary ("
                  << line << " line), variable "
                  << var << ".\n";
}

void
OneDimensionalModel_BCHandler::setBC( const OneDBCFunction_PtrType& funptr,
                                      const std::string& border,
                                      const std::string& line,
                                      const std::string& var,
                                      const Vec2D& matrixrow )
{
    ASSERT_PRE( matrixrow.size() == 2,
                "setBC works only for 2D vectors");

    setBC(funptr, border, line, var);
    Debug( 6311 ) << "[OneDimensionalModel_BCHandler::setBC] imposing matrix row as well.\n";
    M_boundary[border]->matrixrow(line)=matrixrow;
}

inline void
OneDimensionalModel_BCHandler::setBCLeft_internalnode()
{
    M_boundary["left"]->isInternal()=true;
}

inline void
OneDimensionalModel_BCHandler::setBCRight_internalnode()
{
    M_boundary["right"]->isInternal()=true;
}

// ===================================================
// Get Methods
// ===================================================
inline OneDimensionalModel_BC&
OneDimensionalModel_BCHandler::BC( const std::string& bound )
{
    return *(M_boundary[bound]);
}

OneDimensionalModel_BCHandler::OneDBCFunction_PtrType&
OneDimensionalModel_BCHandler::leftBCFunction( const std::string& line )
{
    return M_boundary["left"]->rhs(line);
}

OneDimensionalModel_BCHandler::OneDBCFunction_PtrType&
OneDimensionalModel_BCHandler::rightBCFunction( const std::string& line )
{
    return M_boundary["right"]->rhs(line);
}

inline bool&
OneDimensionalModel_BCHandler::leftBCReady( const std::string& line )
{
    return M_boundarybool["left"][line];
}

inline bool&
OneDimensionalModel_BCHandler::rightBCReady( const std::string& line )
{
    return M_boundarybool["right"][line];
}

}

