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
#include <lifemc/lifesolver/OneDimensionalModel_Definitions.hpp>
#include <lifemc/lifefem/OneDimensionalModel_BC.hpp>

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

    typedef OneDimensionalModel_BC::BCFunction_Type         BCFunction_Type;
    typedef OneDimensionalModel_BC::BCFunction_PtrType         BCFunction_PtrType;
    typedef OneDimensionalModel_BC::BCFunction_Default_PtrType BCFunction_Default_PtrType;

    typedef OneDimensionalModel_BC::Flux_PtrType               Flux_PtrType;
    typedef OneDimensionalModel_BC::Source_PtrType             Source_PtrType;
    typedef OneDimensionalModel_BC::Solution_PtrType           Solution_PtrType;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    OneDimensionalModel_BCHandler();

    //! Destructor
    ~OneDimensionalModel_BCHandler() {}

    //@}


    //! @name Methods
    //@{

    //! Apply boundary conditions
    void applyBC ( const Real&             time,
                   const Solution_PtrType& solution,
                   const Flux_PtrType&     flux,
                         Container2D_Type& left_BC_dir,
                         Container2D_Type& right_BC_dir );

    //@}


    //! @name Set Methods
    //@{

    void setBC( const BCFunction_Type& BCfunction, const OneD_BCSide& side,
                const OneD_BCLine& line,           const OneD_BC& bcType );

    void setBC( const BCFunction_Type& BCfunction, const OneD_BCSide& side,
                const OneD_BCLine& line,           const OneD_BC& bcType,
                const Container2D_Type& matrixrow );

    void setDefaultBC( const Flux_PtrType     fluxFun,
                       const Source_PtrType   sourceFun,
                       const Solution_PtrType solution );

    //@}


    //! @name Get Methods
    //@{

    const OneDimensionalModel_BC& BC( const OneD_BCSide& side );

    const bool& BCReady( const OneD_BCSide& side, const OneD_BCLine& line );

    //@}

private:

    std::map< OneD_BCSide, boost::shared_ptr< OneDimensionalModel_BC > > M_boundary;
    std::map< OneD_BCSide, std::map< OneD_BCLine, bool > >               M_boundarySet;

    std::vector < BCFunction_Default_PtrType >                           M_defaultBC;
};

}

#endif //ONEDIMENSIONALMODEL_BCHANDLER_H
