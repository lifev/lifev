/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-09-02

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/**
 \file MS_Coupling_BoundaryCondition.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-09-02
 */

#ifndef __MS_Coupling_BoundaryCondition_H
#define __MS_Coupling_BoundaryCondition_H 1

#include <lifemc/lifefem/BCInterface.hpp>

#include <lifemc/lifesolver/MS_PhysicalCoupling.hpp>
#include <lifemc/lifesolver/MS_Model_Fluid3D.hpp>

namespace LifeV {

//! MS_Coupling_BoundaryCondition - Coupling condition for standard boundary conditions
/*!
 *  The MS_Coupling_BoundaryCondition class is an implementation of the MS_PhysicalCoupling
 *  for applying standard boundary conditions on the models.
 *
 *  @author Cristiano Malossi
 */
class MS_Coupling_BoundaryCondition: public MS_PhysicalCoupling
{
public:

    typedef MS_PhysicalCoupling super;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Coupling_BoundaryCondition();

    //! Copy constructor
    /*!
     * \param BoundaryCondition	- MS_Coupling_BoundaryCondition
     */
    MS_Coupling_BoundaryCondition( const MS_Coupling_BoundaryCondition& BoundaryCondition );

    //! Destructor
    ~MS_Coupling_BoundaryCondition() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param boundaryCondition - MS_Coupling_BoundaryCondition
     */
    MS_Coupling_BoundaryCondition& operator=( const MS_Coupling_BoundaryCondition& boundaryCondition );

    //@}


    //! @name MultiScale Physical Coupling
    //@{

    //! Setup the data of the coupling
    void SetupData( void );

    //! Setup the coupling
    void SetupCoupling( void );

    //! Setup parameters for the implicit coupling (DO NOTHING)
    void SetupImplicitCoupling( ContainerOfVectors< EpetraVector >& /*couplingVariables*/,
                                ContainerOfVectors< EpetraVector >& /*couplingResiduals*/) {}

    //! Update the values of the coupling variables (DO NOTHING)
    void UpdateCouplingVariables( void ) {}

    //! Update the values of the coupling residual (DO NOTHING)
    void UpdateCouplingResiduals( void ) {}

    //! Display some information about the coupling
    void ShowMe( void );

    //@}

private:

    //! @name Private Methods
    //@{


    //! Apply the boundary condition to the specific model
    /*!
     * \param physicalModel - shared_ptr to the specific model
     */
    template< class model >
    inline void ApplyBoundaryConditions( const PhysicalModel_ptr& physicalModel );

    //@}

    std::vector< BCName > M_list;
    UInt                  M_listSize;

};

//! Factory create function
inline MS_PhysicalCoupling* createBoundaryCondition()
{
    return new MS_Coupling_BoundaryCondition();
}

// ===================================================
//! Template implementation
// ===================================================
template< class model >
inline void
MS_Coupling_BoundaryCondition::ApplyBoundaryConditions( const PhysicalModel_ptr& physicalModel )
{
    model *Model = dynamic_cast< model * > ( &( *physicalModel ) );

    for ( UInt i( 0 ); i < M_listSize; ++i )
    {
        Model->GetBCInterface().ReadExternalBC( M_list[i], "boundary_conditions/", M_dataFile );

        Model->GetBCInterface().GetDataContainer().SetName( "BC_" + number2string( Model->GetID() ) + "_Flag_" + number2string( M_flags[i] ) + "_" + M_list[i] );
        Model->GetBCInterface().GetDataContainer().SetFlag( M_flags[i] );

        Model->GetBCInterface().InsertExternalBC();
    }

    //MPI Barrier
    M_comm->Barrier();
}

} // Namespace LifeV

#endif /* __MS_Coupling_BoundaryCondition_H */
