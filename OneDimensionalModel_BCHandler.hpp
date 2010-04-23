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

#ifndef ONEDIMENSIONALMODEL_BCHANDLER_H
#define ONEDIMENSIONALMODEL_BCHANDLER_H

// LIFEV - MATHCARD
#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BCFunction.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifesolver/OneDimensionalModel_Data.hpp>

namespace LifeV {

//! OneDimensionalModel_BCHandler - Class featuring methods to handle boundary conditions.
/*!
 *  @author Lucia Mirabella, Tiziano Passerini
 */
class OneDimensionalModel_BCHandler
{
public:

    //! @name Type definitions
    //@{

    typedef boost::shared_ptr<OneDimensionalModel_BCFunction>    OneDBCFunction_PtrType;

    typedef OneDimensionalModel_BCFunction::Flux_Type            Flux_Type;
    typedef OneDimensionalModel_BCFunction::Flux_PtrType         Flux_PtrType;

    typedef OneDimensionalModel_BCFunction::Source_Type          Source_Type;
    typedef OneDimensionalModel_BCFunction::Source_PtrType       Source_PtrType;

    typedef OneDimensionalModel_BCFunction::Data_Type            Data_Type;
    typedef OneDimensionalModel_BCFunction::Mesh_Type            Mesh_Type;

    typedef OneDimensionalModel_BCFunction::FESpace_Type         FESpace_Type;

    typedef OneDimensionalModel_BCFunction::LinearSolver_Type    LinearSolver_Type;
    typedef OneDimensionalModel_BCFunction::Vector_Type          Vector_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCHandler( const std::vector<Vector_Type>& U_thistime,
                                   const Flux_PtrType              fluxFun,
                                   const Real&                     dimDof );

    //! Destructor
    ~OneDimensionalModel_BCHandler() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions
    void applyBC (const Real&  time_val,
                        Vec2D& left_BC_dir,
                        Vec2D& right_BC_dir );

    //@}


    //! @name Set Methods
    //@{

    void setBC( const OneDBCFunction_PtrType& funptr, const std::string& border,
                const std::string& line,              const std::string& var );

    void setBC( const OneDBCFunction_PtrType& funptr, const std::string& border,
                const std::string& line,              const std::string& var,
                const Vec2D & matrixrow );

    inline void setBCLeft_internalnode();
    inline void setBCRight_internalnode();

    void setDefaultBC( const FESpace_Type&  fespace,
                       const Source_PtrType sourceFun,
                       const Real&          dt);

    //@}


    //! @name Get Methods
    //@{

    inline OneDimensionalModel_BC& BC( const std::string& bound );

    OneDBCFunction_PtrType& leftBCFunction( const std::string& line );
    OneDBCFunction_PtrType& rightBCFunction( const std::string& line );

    inline bool& leftBCReady( const std::string& line );
    inline bool& rightBCReady( const std::string& line );

    //@}

private:

    //! Trick to use strings in C++ switch construct
    std::map<std::string, OneDBCStringValue>                 M_OneDimensionalModel_BCHandlerMapStringValues;

    //! Reference to the solver current unknowns (U)
    const std::vector<Vector_Type>&                          M_U_thistime;

    //! Reference to the solver non linear flux functions
    Flux_PtrType                                             M_fluxFun;

    std::map<std::string, boost::shared_ptr<OneDimensionalModel_BC > > M_boundary;

    std::map<std::string, std::map< std::string, bool > >    M_boundarybool;
};

}

#endif //ONEDIMENSIONALMODEL_BCHANDLER_H
