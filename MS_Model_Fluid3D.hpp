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
#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

// LifeV includes
#include <life/lifemesh/partitionMesh.hpp>
#include <life/lifesolver/dataNavierStokes.hpp>
#include <life/lifemesh/dataMesh.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifefem/bdfNS_template.hpp>
#include <life/lifefilters/ensight.hpp>
#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
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

    //! @name Public Types
    //@{

    typedef MS_PhysicalModel                      super;

    typedef RegionMesh3D< LinearTetra >           Mesh_Type;
    typedef partitionMesh< Mesh_Type >            PartitionMesh_Type;

    typedef OseenShapeDerivative< Mesh_Type >     Fluid_Type;
    typedef Fluid_Type::vector_type               FluidVector_Type;
    typedef boost::shared_ptr< FluidVector_Type > FluidVector_PtrType;

    typedef Exporter< Mesh_Type >                 IOFile_Type;
    typedef boost::shared_ptr< IOFile_Type >      IOFile_PtrType;

    typedef Ensight< Mesh_Type >                  EnsightIOFile_Type;
#ifdef HAVE_HDF5
    typedef Hdf5exporter< Mesh_Type >             HDF5IOFile_Type;
#endif

    typedef BCHandler                             BC_Type;
    typedef boost::shared_ptr< BC_Type >          BC_PtrType;
    typedef BCInterface< Fluid_Type >             BCInterface_Type;
    typedef boost::shared_ptr< BCInterface_Type > BCInterface_PtrType;

    typedef BdfTNS< FluidVector_Type >            BDF_Type;
    typedef DataNavierStokes		              Data_Type;

    typedef FESpace< Mesh_Type, EpetraMap >       FESpace_Type;
    typedef boost::shared_ptr< FESpace_Type >     FESpace_PtrType;

    //@}


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

    //! Setup the linear model
    void SetupLinearModel();

    //! Update the linear system matrix and vectors
    void UpdateLinearModel();

    //! Solve the linear problem
    void SolveLinearModel( bool& SolveLinearSystem );

    //@}


    //! @name Set Methods
    //@{

    //! Set the solution vector
    /*!
     * @param Solution Solution vector
     */
    void SetSolution( const FluidVector_PtrType& Solution );

    //@}


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BCInterface_Type& GetBCInterface();

    //! Get the density on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return density value
     */
    Real GetBoundaryDensity( const BCFlag& /*flag*/) const;

    //! Get the viscosity on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return viscosity value
     */
    Real GetBoundaryViscosity( const BCFlag& /*flag*/) const;

    //! Get the area on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return area value
     */
    Real GetBoundaryArea( const BCFlag& Flag ) const;

    //! Get the flow rate on a specific boundary face of the model
    /*!
     * @param flag flag of the boundary face
     * @return flow rate value
     */
    Real GetBoundaryFlowRate( const BCFlag& Flag ) const;

    //! Get the integral of the pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return pressure value
     */
    Real GetBoundaryPressure( const BCFlag& Flag ) const;

    //! Get the integral of the dynamic pressure (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @return dynamic pressure value
     */
    Real GetBoundaryDynamicPressure( const BCFlag& Flag ) const;

    //! Get the value of the Lagrange multiplier associated to a specific boundary face
    /*!
     * @param Flag flag of the boundary face
     * @return Lagrange multiplier value
     */
    Real GetBoundaryLagrangeMultiplier( const BCFlag& Flag ) const;

    //! Get the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param StressType Type of approximation for the stress
     * @return stress value
     */
    Real GetBoundaryStress( const BCFlag& Flag, const stressTypes& StressType = StaticPressure ) const;

    //! Get the variation of the flow rate (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flow rate
     */
    Real GetBoundaryDeltaFlowRate( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the pressure
     */
    Real GetBoundaryDeltaPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the total pressure (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the dynamic pressure
     */
    Real GetBoundaryDeltaDynamicPressure( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the Lagrange multiplier associated to a specific boundary face, using the linear model
    /*!
     * @param Flag flag of the boundary face
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return Lagrange multiplier value
     */
    Real GetBoundaryDeltaLagrangeMultiplier( const BCFlag& Flag, bool& SolveLinearSystem );

    //! Get the variation of the integral of the normal stress (on a specific boundary face)
    /*!
     * @param flag flag of the boundary face
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @param StressType Type of approximation for the stress
     * @return variation of the stress
     */
    Real GetBoundaryDeltaStress( const BCFlag& Flag, bool& SolveLinearSystem, const stressTypes& StressType = StaticPressure );

    //@}


    //! @name Get Methods
    //@{

    //! Get the data container
    /*!
     * @return Data container
     */
    const Data_Type& GetData() const;

    //! Get the solution vector
    /*!
     * @return Solution vector
     */
    const FluidVector_Type& GetSolution() const;

    //@}

private:

    //! @name Private Methods
    //@{

    //! Setup the global data of the model.
    /*!
     * In particular, it replaces the default local values with the ones in the global container.
     * If a value is already specified in the data file, do not perform the replacement.
     *
     * @param FileName File name of the specific model.
     */
    void SetupGlobalData( const std::string& FileName );

    //! Setup the exporter and the importer
    /*!
     * @param FileName File name of the specific model.
     */
    void SetupExporterImporter( const std::string& FileName );

    //! Setup the mesh for the fluid problem
    void SetupMesh();

    //! Setup the FE space for pressure and velocity
    void SetupFEspace();

    //! Setup the DOF of the model
    void SetupDOF();

    //! Setup the offset for fluxes boundary conditions
    void SetupBCOffset( const boost::shared_ptr< BC_Type >& BC );

    //! Initialize the solution.
    void InitializeSolution();

    //! Impose the coupling perturbation on the correct BC inside the BCHandler
    void ImposePerturbation();

    //! Reset all the coupling perturbations imposed on the BCHandler
    void ResetPerturbation();

    Real BCFunctionDelta_Zero( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);
    Real BCFunctionDelta_One ( const Real& /*t*/, const Real& /*x*/, const Real& /*y*/, const Real& /*z*/, const ID& /*id*/);

    //@}

    IOFile_PtrType                          M_exporter;
    IOFile_PtrType                          M_importer;

    std::string                             M_fileName; //TODO Temporary - To be removed

    // Fluid problem
    boost::shared_ptr< Fluid_Type >         M_fluid;
    BCInterface_PtrType                     M_BC;
    boost::shared_ptr< BDF_Type >           M_BDF;
    boost::shared_ptr< Data_Type >          M_data;
    boost::shared_ptr< DataMesh >           M_dataMesh;
    boost::shared_ptr< PartitionMesh_Type > M_mesh;
    boost::shared_ptr< EpetraMap >          M_map;
    FluidVector_PtrType                     M_solution;

    // Linear Fluid problem
    BC_PtrType                              M_linearBC;
    bool                                    M_updateLinearModel;

    // FE spaces
    FESpace_PtrType                         M_uFESpace;
    FESpace_PtrType                         M_pFESpace;

    // Lagrange multipliers
    UInt                                    M_lmDOF;

    // Problem coefficients
    Real                                    M_alpha;
    FluidVector_PtrType                     M_beta;
    FluidVector_PtrType                     M_RHS;

    // NS parameters
    UInt                                    M_subiterationsMaximumNumber;
    Real                                    M_tolerance;
    generalizedAitken< FluidVector_Type >   M_generalizedAitken;

    // BC Functions for tangent problem
    BCFunctionBase                          M_BCBaseDelta_Zero;
    BCFunctionBase                          M_BCBaseDelta_One;
};

//! Factory create function
inline MS_PhysicalModel* MS_createFluid3D()
{
    return new MS_Model_Fluid3D();
}

} // Namespace LifeV

#endif /* MS_Model_Fluid3D_H */
