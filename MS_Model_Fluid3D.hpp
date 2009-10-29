/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Cristiano Malossi <cristiano.malossi@epfl.ch>
 Date: 2009-03-12

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
 \file MS_Model_Fluid3D.hpp
 \author Cristiano Malossi <cristiano.malossi@epfl.ch>
 \date 2009-03-12
 */

#ifndef __MS_Model_Fluid3D_H
#define __MS_Model_Fluid3D_H 1

#include <life/lifearray/EpetraMatrix.hpp>
#include <life/lifealg/EpetraMap.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#include <life/lifesolver/OseenShapeDerivative.hpp>

#include <lifemc/lifefem/BCInterface.hpp>

#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

namespace LifeV {

//! MS_Model_Fluid3D - MultiScale model for 3D Fluid simulations
/*!
 *  The MS_Model_Fluid3D class is an implementation of the MS_PhysicalModel
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 *
 *  @author Cristiano Malossi
 */
class MS_Model_Fluid3D: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                  super;

    typedef RegionMesh3D< LinearTetra >       MeshType;
    typedef Ensight< MeshType >               OutputType;
    typedef OseenShapeDerivative< MeshType >  FluidType;

    typedef FluidType::vector_type            FluidVectorType;

    typedef BCInterface< FluidType >          FluidBCType;
    typedef BdfTNS< FluidVectorType >         FluidBDFType;
    typedef DataNavierStokes< MeshType >      FluidDataType;
    typedef partitionMesh< MeshType >         FluidMeshType;

    typedef FESpace< MeshType, EpetraMap >    FESpaceType;

    //! @name Constructors & Destructor
    //@{

    //! Constructor
    MS_Model_Fluid3D();

    //! Copy constructor
    /*!
     * \param Fluid3D - MS_Model_Fluid3D
     */
    MS_Model_Fluid3D( const MS_Model_Fluid3D& Fluid3D );

    //! Destructor
    ~MS_Model_Fluid3D() {}

    //@}


    //! @name Methods
    //@{

    //! Operator=
    /*!
     * \param fluid3D - MS_Model_Fluid3D
     */
    MS_Model_Fluid3D& operator=( const MS_Model_Fluid3D& fluid3D );

    //! Setup the data of the linear model
    void SetupLinearData( void );

    //! Setup the linear model
    void SetupLinearModel( void );

    //! Build the linear system matrix and vectors
    void BuildLinearSystem( void );

    //! Update the linear system matrix and vectors
    void UpdateLinearSystem( void );

    //! Solve the linear problem
    void SolveLinearSystem( void );

    //@}


    //! @name MultiScale Physical Model
    //@{

    //! Setup the data of the model
    void SetupData( void );

    //! Setup the model
    void SetupModel( void );

    //! Build the system matrix and vectors
    void BuildSystem( void );

    //! Update the system matrix and vectors
    void UpdateSystem( void );

    //! Solve the problem
    void SolveSystem( void );

    //! Save the solution
    void SaveSolution( void );

    //! Display some information about the model
    void ShowMe( void );

    //@}


    //! @name Get functions
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    FluidBCType& GetBC( void )
    {
        return *M_FluidBC;
    }

    //! Get the BCInterface container of the boundary conditions of the linear model
    FluidBCType& GetLinearBC( void )
    {
        return *M_LinearFluidBC;
    }

    //! Get the area on a specific boundary face of the model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetArea( const BCFlag& flag ) const
    {
        return M_Fluid->area( flag );
    }

    //! Get the flux on a specific boundary face of the model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetFlux( const BCFlag& flag ) const
    {
        return M_Fluid->flux( flag );
    }

    //! Get the flux on a specific boundary face of the linear model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetLinearFlux( const BCFlag& flag ) const
    {
        return M_Fluid->GetLinearFlux( flag );
    }

    //! Get the pressure on a specific boundary face of the model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetPressure( const BCFlag& flag ) const
    {
        return M_Fluid->pressure( flag );
    }

    //! Get the pressure on a specific boundary face of the linear model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetLinearPressure( const BCFlag& flag ) const
    {
        return M_Fluid->GetLinearPressure( flag );
    }

    //! Get the density on a specific boundary face of the model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetDensity( const BCFlag& /*flag*/) const
    {
        return M_FluidData->density();
    }

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * \param flag - flag of the boundary face
     */
    Real GetViscosity( const BCFlag& /*flag*/) const
    {
        return M_FluidData->viscosity();
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the FE space for pressure and velocity
    void setupFEspace( void );

    //! Setup the DOF of the model
    void setupDOF( void );

    //! Compute the number fluxes applied to the model
    UInt imposeFluxes( const boost::shared_ptr< FluidBCType >& BC );

    //@}

    boost::shared_ptr< OutputType >       M_output;

    // Fluid problem
    boost::shared_ptr< FluidType >        M_Fluid;
    boost::shared_ptr< FluidBCType >      M_FluidBC;
    boost::shared_ptr< FluidBDFType >     M_FluidBDF;
    boost::shared_ptr< FluidDataType >    M_FluidData;
    boost::shared_ptr< FluidMeshType >    M_FluidMesh;
    boost::shared_ptr< EpetraMap >        M_FluidFullMap;
    boost::shared_ptr< FluidVectorType >  M_FluidSolution;

    // Linear Fluid problem
    boost::shared_ptr< FluidBCType >      M_LinearFluidBC;

    // FE spaces
    boost::shared_ptr< FESpaceType >      M_uFESpace;
    boost::shared_ptr< FESpaceType >      M_pFESpace;

    // Degrees of freedom of the problem
    UInt                                  M_uDOF;
    UInt                                  M_pDOF;

    // Problem coefficients
    Real                                  M_alpha;
    boost::shared_ptr< FluidVectorType >  M_beta;
    boost::shared_ptr< FluidVectorType >  M_RHS;
};

//! Factory create function
inline MS_PhysicalModel* createFluid3D()
{
    return new MS_Model_Fluid3D();
}

} // Namespace LifeV

#endif /* __MS_Model_Fluid3D_H */
