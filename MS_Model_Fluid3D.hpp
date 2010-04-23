//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2009 EPFL, Politecnico di Milano, INRIA

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
 *  @brief MultiScale Model Fluid3D
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 12-03-2009
 */


#ifndef MS_Model_Fluid3D_H
#define MS_Model_Fluid3D_H 1

// Mathcard includes
#include <lifemc/lifealg/AztecOOPreconditioner.hpp>
#include <lifemc/lifefem/BCInterface.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

// LifeV includes
#include <life/lifefilters/ensight.hpp>
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/exporter.hpp>
#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
#else
    #include <life/lifefilters/ensight.hpp>
#endif
#include <life/lifesolver/OseenShapeDerivative.hpp>

namespace LifeV {

//! MS_Model_Fluid3D - MultiScale model for 3D Fluid simulations
/*!
 *  @author Cristiano Malossi
 *
 *  The MS_Model_Fluid3D class is an implementation of the MS_PhysicalModel
 *  for 3D Fluid problem (in particular Oseen with Shape Derivatives).
 */
class MS_Model_Fluid3D: public virtual MS_PhysicalModel
{
public:

    typedef MS_PhysicalModel                  super;

    typedef RegionMesh3D< LinearTetra >       MeshType;
    typedef OseenShapeDerivative< MeshType >  FluidType;

#ifdef HAVE_HDF5
    typedef Hdf5exporter< MeshType >          IOFileType;
#else
    typedef Ensight< MeshType >               IOFileType;
#endif

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
     * @param Fluid3D MS_Model_Fluid3D
     */
    MS_Model_Fluid3D( const MS_Model_Fluid3D& Fluid3D );

    //! Destructor
    ~MS_Model_Fluid3D() {}

    //@}


    //! @name Operators
    //@{

    //! Operator=
    /*!
     * @param fluid3D MS_Model_Fluid3D
     * @return reference to a copy of the class
     */
    MS_Model_Fluid3D& operator=( const MS_Model_Fluid3D& fluid3D );

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the model.
    void SetupModel();

    //! Build the initial system (matrix and vectors).
    void BuildSystem();

    //! Update the system for (matrix and vectors).
    void UpdateSystem();

    //! Solve the problem.
    void SolveSystem();

    //! Save the solution
    void SaveSolution();

    //! Display some information about the model.
    void ShowMe();

    //@}


    //! @name Methods
    //@{

    //! Setup the data of the linear model
    /*!
     * @param FileName Name of data file.
     */
    void SetupLinearData( const std::string& FileName );

    //! Setup the linear model
    void SetupLinearModel();

    //! Update the linear system matrix and vectors
    void UpdateLinearModel();

    //! Solve the linear problem
    void SolveLinearModel( bool& SolveLinearSystem );

    //@}


    //! @name Get Methods
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    FluidBCType& GetBC();

    //! Get the BCInterface container of the boundary conditions of the linear model
    /*!
     * @return BCInterface container
     */
    FluidBCType& GetLinearBC();

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real GetDensity( const BCFlag& /*flag*/) const;

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real GetViscosity( const BCFlag& /*flag*/) const;

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real GetArea( const BCFlag& Flag ) const;

    //! Get the flux on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flux value
     */
    Real GetFlux( const BCFlag& Flag ) const;

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real GetPressure( const BCFlag& Flag ) const;

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real GetDynamicPressure( const BCFlag& Flag ) const;

    //! Get the value of the Lagrange multiplier associated to a specific boundary face
    /*!
     * @param Flag flag of the boundary face
     * @return Lagrange multiplier value
     */
    Real GetLagrangeMultiplier( const BCFlag& Flag ) const;

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param StressType Type of approximation for the stress
     * @return stress value
     */
    Real GetStress( const BCFlag& Flag, const stressTypes& StressType = StaticPressure ) const;

    //! Get the variation of the flux (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flux
     */
    Real GetDeltaFlux( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real GetDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the total pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the dynamic pressure
     */
    Real GetDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the Lagrange multiplier associated to a specific boundary face, using the linear model
    /*!
     * @param Flag flag of the boundary face
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return Lagrange multiplier value
     */
    Real GetDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @param StressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real GetDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType = StaticPressure );

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the FE space for pressure and velocity
    void SetupFEspace();

    //! Setup the DOF of the model
    void SetupDOF();

    //! Setup the offset for fluxes boundary conditions
    void SetupBCOffset( const boost::shared_ptr< FluidBCType >& BC );

    //! Setup the solution.
    void SetupSolution();

    //@}

    boost::shared_ptr< IOFileType >       M_exporter;
    boost::shared_ptr< IOFileType >       M_importer;

    std::string                           M_FileName; //TODO Temporary - To be removed

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
    bool                                  M_UpdateLinearModel;

    // FE spaces
    boost::shared_ptr< FESpaceType >      M_uFESpace;
    boost::shared_ptr< FESpaceType >      M_pFESpace;

    // Degrees of freedom of the problem
    UInt                                  M_uDOF;  // Velocity degrees of freedom (one component)
    UInt                                  M_pDOF;  // Pressure degrees of freedom
    UInt                                  M_lmDOF; // Lagrange Multipliers degrees of freedom

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

#endif /* MS_Model_Fluid3D_H */
