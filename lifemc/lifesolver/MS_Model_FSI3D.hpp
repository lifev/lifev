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
 *  @brief MultiScale Model FSI3D
 *
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  @date 19-04-2010
 *
 */

#ifndef MS_MODEL_FSI3D_H
#define MS_MODEL_FSI3D_H 1

// Mathcard includes
#include <lifemc/lifealg/AztecOOPreconditioner.hpp>
#include <lifemc/lifesolver/BCInterface.hpp>
#include <lifemc/lifesolver/MS_PhysicalModel.hpp>

// LifeV includes
#include <life/lifesolver/FSISolver.hpp>

#include <life/lifefilters/ensight.hpp>
#include <life/lifefilters/noexport.hpp>
#ifdef HAVE_HDF5
    #include <life/lifefilters/hdf5exporter.hpp>
#endif

namespace LifeV {

//! MS_Model_FSI3D - MultiScale model for 3D FSI simulations
/*!
 *  @author Paolo Crosetto
 */
class MS_Model_FSI3D: public virtual MS_PhysicalModel
{
public:

    //! @name Public Types
    //@{

    typedef MS_PhysicalModel                                                               super;
    typedef FSISolver                                                                      FSISolverType;
    typedef boost::shared_ptr< FSISolverType >                                             FSISolverPtrType;
    typedef boost::shared_ptr< LifeV::Exporter<LifeV::RegionMesh3D<LifeV::LinearTetra> > > filter_ptrtype;
    typedef LifeV::Ensight<LifeV::FSIOperator::mesh_type>                                  ensightfilter_type;
    typedef boost::shared_ptr<ensightfilter_type>                                          ensightfilter_ptrtype;
#ifdef HAVE_HDF5
    typedef LifeV::Hdf5exporter<LifeV::FSIOperator::mesh_type>                             hdf5filter_type;
    typedef boost::shared_ptr<hdf5filter_type>                                             hdf5filter_ptrtype;
#endif
    typedef RegionMesh3D<LinearTetra>                                                      mesh_type;
    typedef OseenShapeDerivative< mesh_type >                                              fluid_type;
    typedef VenantKirchhofSolver< mesh_type >                                              solid_type;
    typedef FSISolverType::vector_type                                                     vector_type;
    typedef FSISolverType::vector_ptrtype                                                  vector_ptrtype;
    typedef BCInterface< FSIOperator >                                                     solid_bc_type;

    typedef BCHandler                                                                      BC_Type;
    typedef BCInterface< FSIOperator >                                                     BCInterface_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor
    MS_Model_FSI3D();

    //! Destructor
    ~MS_Model_FSI3D();

    //@}


    //! @name MultiScale PhysicalModel Virtual Methods
    //@{

    //! Setup the data of the model.
    /*!
     * @param FileName Name of data file.
     */
    void SetupData( const std::string& FileName );

    //! Setup the global data of the model.
    /*!
     * @param PhysicalData Global data container.
     */
    void SetupGlobalData( const boost::shared_ptr< MS_PhysicalData >& PhysicalData );

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


    //! @name Get Methods (couplings)
    //@{

    //! Get the BCInterface container of the boundary conditions of the model
    /*!
     * @return BCInterface container
     */
    BCInterface_Type& GetBCInterface();

    //! Get the BCInterface container of the boundary conditions of the linear model
    /*!
     * @return BCInterface container
     */
    BCInterface_Type& GetLinearBCInterface();

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

    //! Get the variation of the flux (on a specific boundary face) using the linear model
    /*!
     * @param Flag flag of the boundary face on which quantity should be computed
     * @param SolveLinearSystem a flag to which determine if the linear system has to be solved
     * @return variation of the flux
     */
    Real GetBoundaryDeltaFlux( const BCFlag& Flag, bool& SolveLinearSystem );

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

private:

    //! @name Private Methods
    //@{

    void SetUpExporter( filter_ptrtype& exporter, const std::string name, const GetPot& dataFile );
    void SetExporterFluid( const filter_ptrtype& exporter );
    void SetExporterSolid( const filter_ptrtype& exporter );
    void updateBCSegregated( const std::string& dataFileName );
    void updateBCMonolithic( const std::string& dataFileName );

    //@}

    FSISolverPtrType                      M_solver;
    filter_ptrtype                        M_exporterFluid;
    filter_ptrtype                        M_exporterSolid;
    vector_ptrtype                        M_solidDisp;
    vector_ptrtype                        M_solidVel;
    vector_ptrtype                        M_velAndPressure;
    vector_ptrtype                        M_fluidDisp;

    boost::shared_ptr< BCInterface_Type >    M_fluidBC;
    boost::shared_ptr< solid_bc_type >       M_solidBC;
    boost::shared_ptr< solid_bc_type >       M_harmonicExtensionBC;
    boost::shared_ptr< BCInterface_Type >    M_linearizedFluidBC;
    boost::shared_ptr< solid_bc_type >       M_linearizedSolidBC;
};

//! Factory create function
inline MS_PhysicalModel* createModelFSI3D()
{
    return new MS_Model_FSI3D();
}

} // Namespace LifeV

#endif /* MS_MODEL_FSI3D_H */
