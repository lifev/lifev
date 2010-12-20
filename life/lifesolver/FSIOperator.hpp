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
//  \include ../../../mathcard/testsuite/test_monolithic/fluidstructure.dox


/*!

    @file
    @brief Pure virtual operator class for FSI solvers

    \include ../../doc/api/bibliography/newton
    \include ../../doc/api/bibliography/fluidstructure

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer  Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010

    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton method implemente in NonlinearRichardson. The fluid
    and structure classes are member of this class and different formulations (e.g. Monolithic \ref CDFQ , segregated
    Newton \ref FM05 , Dirichlet--Neumann \ref DDFQ06 , Robin Neumann \ref BNV08 )

 */



#ifndef FSIOPERATOR_H
#define FSIOPERATOR_H

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"


#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/regionMesh3D.hpp>

#include <life/lifealg/generalizedAitken.hpp>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifefem/bdf_template.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifefilters/HDF5Filter3DMesh.hpp>

#include <life/lifesolver/DataFSI.hpp>
#include <life/lifesolver/OseenShapeDerivative.hpp>
//#include <life/lifesolver/NonLinearVenantKirchhofSolver.hpp>
#include <life/lifesolver/LinearVenantKirchhofSolver.hpp>
#include <life/lifesolver/HarmonicExtensionSolver.hpp>

//#include <life/lifefem/refFE.hpp>
//#include <life/lifesolver/fixedPointBase.hpp>


#define FLUID 1
#define SOLID 0

namespace LifeV
{
/*!

  \include ../../doc/api/bibliography/fluidstructure.dox
  \include ../../doc/api/bibliography/newton.dox
  @class FSIOperator
  @brief Fluid-Structure Interface operator class

  This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
  Jacobian matrix, which make it suited for the generalized Newton algorithm implemented in NonlinearRichardson.hpp. The fluid
  and structure classes are members of this class and different formulations (e.g. Monolithic \ref CDFQ10 , segregated
  Newton \ref{FM05} , Dirichlet--Neumann \ref DDFQ06 , Robin Neumann \ref BNV08 ) are implemented in the derived classes.

  @see
  FSIExactJacobian
  FSIFixedPoint
  FSIMonolithic
*/
class FSIOperator
{

public:

    /** @name Typedefs
     */
    //@{

    typedef RegionMesh3D<LinearTetra>                                               mesh_Type;
#ifdef HAVE_HDF5
    typedef HDF5Filter3DMesh<mesh_Type>                                             meshFilter_Type;
#endif
    typedef OseenShapeDerivative   <mesh_Type>                                      fluid_Type;
    typedef VenantKirchhofSolver   <mesh_Type>                                      solid_Type;
    typedef HarmonicExtensionSolver<mesh_Type>                                      meshMotion_Type;
    typedef OseenShapeDerivative   <mesh_Type>                                      fluidLin_Type;
    typedef VenantKirchhofSolver   <mesh_Type>                                      solidLin_Type;
    typedef boost::shared_ptr<fluid_Type>                                           fluidPtr_Type;
    typedef boost::shared_ptr<solid_Type>                                           solidPtr_Type;
    typedef boost::shared_ptr<meshMotion_Type>                                      meshMotionPtr_Type;
    typedef boost::shared_ptr<fluidLin_Type>                                        fluidLinPtr_Type;
    typedef boost::shared_ptr<solidLin_Type>                                        solidLinPtr_Type;
    typedef fluid_Type::vector_Type/*fluidPtr_Type::vector_Type*/                   vector_Type;
    typedef boost::shared_ptr<vector_Type>                                          vectorPtr_Type;
    typedef fluid_Type::source_Type/*fluidPtr_Type::source_Type*/                   fluidSource_Type;
    typedef solid_Type::source_Type                                                 solidSource_Type;
    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )>           function_Type;
    typedef Real                                                                    ( *bcFunction_Type ) ( const Real&, const Real&,
                                                                                    const Real&, const Real&, const ID& );
    typedef boost::shared_ptr<DofInterface3Dto3D>                                   dofInterface3DPtr_Type;
    typedef boost::shared_ptr<DofInterface3Dto2D>                                   dofInterface2DPtr_Type;
    typedef boost::shared_ptr<BCVectorInterface>                                    bcVectorInterfacePtr_Type;
    typedef fluid_Type::bcHandlerPtr_Type/*fluidPtr_Type::bchandlerPtr_Type*/       fluidBchandlerPtr_Type;
    typedef fluid_Type::bcHandler_Type/*fluidPtr_Type::bchandler_Type*/             fluidBchandler_Type;
    typedef BCHandler                                                               solidBchandler_Type;
    typedef boost::shared_ptr<solidBchandler_Type>                                  solidBchandlerPtr_Type;
    typedef DataFSI                                                                 data_Type;
    typedef boost::shared_ptr<data_Type>                                            dataPtr_Type;
    typedef std::map<ID, ID>::const_iterator                                        iterator_Type;
    typedef singleton<factory<FSIOperator, std::string> >                           FSIFactory_Type;
    typedef Displayer::commPtr_Type/*Displayer::commPtr_Type*/                      commPtr_Type;

     //@}



    /** @name Constructors, Destructor
     */
    //@{

    FSIOperator();

    virtual ~FSIOperator();

    //@}


    /** @name Virtual Public Methods
     */
    //@{

    //! initializes the GetPot data file
    virtual void setDataFile( const GetPot& data );

    //! sets the space discretization parameters
    /*!
      The FE discretization is set accordingly to what specified in the FSIdata member (order of accuracy for the fluid
      pressure, velocity and for the structure).
     */
    virtual void setupFEspace();

    //! partitions the meshes for the fluid and the structure
    /**
       This method is not called when the mesh partition is done offline
     */
    virtual void partitionMeshes();

#ifdef HAVE_HDF5
    //! reads the meshes already partitioned for the fluid and the structure
    /**
       The offline partitioning can avoid the call to partitionMesh and the memory overhead of saving the
       entire unpartitioned mesh. The offline partitioned mesh must be saved in HDF5 format.
       \TODO { This method still does not work for the partitioned algorithms }
     */
    void partitionMeshes( meshFilter_Type& fluidMeshFilter, meshFilter_Type& solidMeshFilter  );
#endif
    //! sets up the correspondences between the fluid and structure degrees of freedom across the interface.
    /**
       This method introduces a non scalable loop, in DOFInterface3Dto3D. It is preferable to avoid it for massively
       parallel computetions, using the offline partitioner. However it is much lighter that the correspondent
       method for partitioned algorithms.
     */
    virtual void setupDOF();

#ifdef HAVE_HDF5
    //!reads from HDF5 file the correspondences between the fluid and structure degrees of freedom across the interface.
    /*!still not implemented for all the FSI formulations*/
    virtual void setupDOF( meshFilter_Type& /*filterMesh*/ ) {}
#endif

    //! setup of the fluid and solid solver classes
    /**
       This method computes the number of fluxes assigned at the boundaries and calls setupFluidSolid(UInt fluxes)
     */
    virtual void setupFluidSolid();

    //! setup of the fluid and solid solver classes
    /**
       Fluid, solid and harmonic extension solvers are instantiated
     */
    virtual void setupFluidSolid( UInt const fluxes );

    //! Setup method
    /*!
    the setup is called for the fluid, structure and harmonic extension solvers
    */
    virtual void setupSystem();

    //!Builds the local matrices
    /**
       The matrix for the harmonic extension, and the constant part of the matrices for the fluid and solid
       solvers are built.
     */
    virtual void buildSystem();

    //!Updates the FSI system
    /**
       The system is updated for the next time iteration
     */
    virtual void updateSystem();

    //!Extrapolates an approximation of the solution
    /**
       Extrapolates the solution for the next time step. This method should be handled by a more general time-advance
       class.
     */
    virtual void couplingVariableExtrap( );

    //! solves the Jacobian system
    /**
       The implementation of this method distinguish the various FSI formulations which derive from this class.
       For this reason it must be pure virtual, snd implemented in the child classes.
       \param muk: unknown solution at the k-th nonlinear iteration
       \param res: residual vector (the right hand side of the Jacobian system)
       \param linearRelTol: tolerance for the nonlinear solver
     */
    virtual void solveJac( vector_Type&       muk,
                           const vector_Type& res,
                           const Real         linearRelTol ) = 0;

    //! Evaluates the nonlinear residual of the FSI system
    /**
       The implementation of this method also depends on the child classes, though it does not characterize them.
       \param res:  residual vector  to be computed
       \param disp: current unknown solution
       \param iter: nonlinear iteration counter. The part of th rhs related to the time discretization is computed only for iter=0
     */
    virtual void evalResidual( vector_Type&        res,
                               const vector_Type&  disp,
                               const UInt          iter ) = 0;

    //!increases the time level
    /**
       the Oseen::shiftSolution() method of the fluid solver is called
     */
    virtual void shiftSolution( );

    //!initializes the solution with functions
    /**
       calls the initialize methods for the subproblems. The mesh velocity is used to compute the convective term in the
       fluid equations
       \param u0: initial fluid velocity
       \param p0: initial fluid pressure
       \param d0: initial solid displacement
       \param w0: initial mesh velocity
     */
    virtual void initialize( FSIOperator::fluidPtr_Type::value_type::function_Type const& u0,
                             FSIOperator::fluidPtr_Type::value_type::function_Type const& p0,
                             FSIOperator::solidPtr_Type::value_type::Function const& d0,
                             FSIOperator::solidPtr_Type::value_type::Function const& w0,
                             FSIOperator::fluidPtr_Type::value_type::function_Type const& /*df0=FSIOperator::solidPtr_Type::value_type::Function()*/)
    {
        Debug( 6220 ) << "FSIOperator:: solid init \n";
        if (this->isSolid())
            solid().initialize(d0, w0, w0);
        Debug( 6220 ) << "FSIOperator:: fluid init \n";
        if (this->isFluid())
            fluid().initialize(u0, p0);
    }

    //! getter for the fluid velocity
    virtual const vectorPtr_Type& un() { return M_un; }

    //! Returns the number of imposed fluxes
    virtual UInt imposeFlux();

    // \todo{kill this method?}
    //virtual void mergeBCHandlers() {}

    //@}

    //!\todo{kill this method}
    //virtual void setFluxBC             (fluidBchandlerPtr_Type /*bc_fluid*/) {}

    //!\todo{kill this method}
    //virtual void setRobinBC             (fluidBchandlerPtr_Type /*bc_solid*/) {}



    //!@name MONOLITHIC Solver Methods - Implemented There
    //{@
    virtual void iterateMesh( const vector_Type& /*disp*/ )     { assert(false); }
    virtual void setupBDF( const vector_Type& /*u0*/ ) { }
    virtual void updateRHS() {}
    virtual void applyBoundaryConditions() {}
    //@}

    //    static VenantKirchhofSolver< FSIOperator::mesh_Type, SolverTrilinos >*    createNonLinearStructure(){ return new NonLinearVenantKirchhofSolver< FSIOperator::mesh_Type, SolverTrilinos >(); }

    //!@name Factory Methods
    //@{
    //! Factory method for the linear elasticity solver
    static VenantKirchhofSolver< FSIOperator::mesh_Type, SolverTrilinos >*    createLinearStructure() { return new LinearVenantKirchhofSolver< FSIOperator::mesh_Type, SolverTrilinos >(); }
    //@}

    //!@name Public Methods
    //@{
    //!Initializes the BDF which should handle the fluid time discretization
    /**
       \todo{a general time advancing class should be used everywhere}
     */
    void initializeBDF( const vector_Type& un );

    //! Creates the Epetra maps for the interface
    /**
       Given a std::map holding the numeration of the interface according to the fluid (first) and the solid (second),
       builds the EpetraMaps for the interface dofs (e.g. M_fluidInterfaceMap). Note that when both the fluid and solid
       meshes are partitioned (e.g. in the monolithic solver) the local dofs are those of the FLUID partition of the
       interface, even when the numeration refers to the solid.
     */
    void createInterfaceMaps(std::map<ID, ID> const& locDofMap);

    //! initializes the fluid solver with vectors
    /**
       \param velAndPressure: initial vector containing the velocity and pressure
       \param displacement: initial vector containing the mesh displacement
     */
    void initializeFluid( const vector_Type& velAndPressure,
                          const vector_Type& displacement );

    //! initializes the solid solver with vectors
    /**
       \param displacement: initial vector containing the structure displacement
       \param velocity: initial vector containing the velocity, used for the initialization of the Newmark scheme
     */
    void initializeSolid( vectorPtr_Type displacement,
                          vectorPtr_Type velocity );

    //!\todo{kill this method}
    //void updateJacobian ( const vector_Type& sol, const int& iter );

    //!moves the mesh using the solution of the harmonic extension equation
    /**
       \param disp displacement of the mesh, must be the difference between the current solution of the HE problem and the one at the previous time step.
     */
    void moveMesh       ( const vector_Type& disp );

//     void solveLinearFluid();
//     void solveLinearSolid();

    //!Method to import an EpetraVector defined on the fluid map (i.e. with the fluid numeration of the dofs) to the interface
    /**
       Note that the output vector will have the solid numeration on the interface! By default in fact the vectors on the
       FSI interface in LifeV are numerated according to the solid.
     */
    void transferFluidOnInterface( const vector_Type& _vec1, vector_Type& _vec2 );

    //works in serial but no yet in parallel
    void transferSolidOnFluid    ( const vector_Type& _vec1, vector_Type& _vec2 );

    //!Method to import an EpetraVector defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       The output vector has the solid numeration of the dofs and is partitioned according to the solid partition. This method is not used in the monolithic solvers.
     */
    void transferSolidOnInterface( const vector_Type& _vec1, vector_Type& _vec2 );
//     void transferInterfaceOnFluid( const vector_Type& _vec1, vector_Type& _vec2 );

    //!Method to import an EpetraVector defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       the output vector have the numeration of the solid, as in transferSolidOnInterface, but is partitioned according to the fluid! This method is not used in the monolithic solvers.
     */
    void transferInterfaceOnSolid( const vector_Type& _vec1, vector_Type& _vec2 );

    //! calls bcManage for a vector
    /**
       calls bcManage for the vector rhs and the BCHandler bch
     */
    void bcManageVectorRHS( const fluidBchandlerPtr_Type& bch, vector_Type& rhs );

    //! Method to set the Robin vector coefficient of the Robin--Neumann coupling scheme (as a constant vector vector)
    void setAlphaf     () { M_Alphaf->epetraVector().PutScalar( M_AlphafCoef ); }
    //! Method to compute the scalar coefficient \f$\alpha\f$ of the Robin--Neumann coupling scheme
    void setAlphafCoef ();
    //! Method calling setAlphaf and setAlphafCoef
    void setStructureToFluidParametres();

    //! method to set the number of fluxes imposed
    void setFluxesNumber                     ( const UInt numLM )
    {
        M_fluxes=numLM;
    }
    //! Setter for the right hand side
    void setRHS                              ( vectorPtr_Type& rhs ) {M_rhs = rhs;}


    //@}



    /** @name  Display Methods
     */
    //@{

    bool isLeader() const;

    //! Getter for the Displayer attribute
    Displayer const& displayer();

    //@}



    /** @name  Get Functions
     */
    //@{

//     vector_Type & displacement()                                        { return *M_lambdaSolid; }
//     vector_Type & displacementOld()                                     { return *M_lambdaSolidOld; }
//     vector_Type & residual();
//     vector_Type const & velocity()                                const { return *M_lambdaDotSolid; }
//     vector_Type & residualFSI();

    const vector_Type& lambdaFluid()                              const { return *M_lambdaFluid; }
    const vector_Type& lambdaSolid()                              const { return *M_lambdaSolid; }
    const vector_Type& lambdaSolidOld()                           const { return *M_lambdaSolidOld; }
    const vector_Type& lambdaDotSolid()                           const { return *M_lambdaDotSolid; }
    const vector_Type& sigmaFluid()                               const { return *M_sigmaFluid; }
    const vector_Type& sigmaSolid()                               const { return *M_sigmaSolid; }

    const vector_Type& lambdaFluidRepeated()                      const { return *M_lambdaFluidRepeated; }
    const vector_Type& lambdaSolidRepeated()                      const { return *M_lambdaSolidRepeated; }
    const vector_Type& lambdaDotSolidRepeated()                   const { return *M_lambdaDotSolidRepeated; }
    const vector_Type& sigmaFluidRepeated()                       const { return *M_sigmaFluidRepeated; }
    const vector_Type& sigmaSolidRepeated()                       const { return *M_sigmaSolidRepeated; }

    //!\todo{remove this method}
    // now residual is Ax-b, but we should decide for b-Ax. In the meantime, we need b-Ax:
    const vector_Type& minusSigmaFluid()                          const { return *M_minusSigmaFluid; }
    //!\todo{remove this method}
    const vector_Type& minusSigmaFluidRepeated()                  const { return *M_minusSigmaFluidRepeated; }

    //!coefficient for the Robin--Neumann coupling scheme
    vector_Type&       Alphaf()                                   const { return *M_Alphaf;}

    commPtr_Type worldComm()                                      const { return M_epetraWorldComm; }

    bool isFluid()                                                const { return M_isFluid; }
    bool isSolid()                                                const { return M_isSolid; }

    bool isLinearFluid()                                          const { return M_linearFluid; }
    bool isLinearSolid()                                          const { return M_linearSolid; }

    int getFluidLeaderId()                                        const { return M_fluidLeader; }
    int getSolidLeaderId()                                        const { return M_solidLeader; }

    //! Getter for the fluid solver
    const fluid_Type& fluid()                                 const { return *M_fluid; }
    //! Getter for the solid solver
    const solid_Type& solid()                                 const { return *M_solid; }
    //! Getter for the harmonic extension solver
    const meshMotion_Type& meshMotion()                       const { return *M_meshMotion; }

    //! Getter-Setter for the fluid solver
    /** \todo{mark as deprecated}*/
    fluidPtr_Type::value_type& fluid()                                     { return *M_fluid; }
    //! Getter-Setter for the solid solver
    /** \todo{mark as deprecated}*/
    solidPtr_Type::value_type& solid()                                     { return *M_solid; }
    //! Getter-Setter for the mesh motion solver
    /** \todo{mark as deprecated}*/
    meshMotionPtr_Type::value_type& meshMotion()                           { return *M_meshMotion; }
//     fluidLinPtr_Type::value_type& fluidLin()                               { return *M_fluidLin; }
//     solidLinPtr_Type::value_type& solidLin()                               { return *M_solidLin; }

    //!getter for the FSI data container
    const data_Type& data()                                       const { return *M_data; }
    //!getter for the fluid data container
    const data_Type::dataFluid_PtrType& dataFluid()               const { return M_data->dataFluid(); }
    //!getter for the solid data container
    const data_Type::dataSolid_PtrType& dataSolid()               const { return M_data->dataSolid(); }

    //!getter for the unpartitioned fluid mesh
    mesh_Type& fluidMesh()                                        const { return *M_fluidMesh; }
    //!getter for the unpartitioned solid mesh
    mesh_Type& solidMesh()                                        const { return *M_solidMesh; }

    // const mesh_Type& fluidMesh()                                  const { return *M_fluidMesh; }
    // const mesh_Type& solidMesh()                                  const { return *M_solidMesh; }

    //!getter for the partitioned fluid mesh
    const partitionMesh< mesh_Type >& fluidMeshPart()             const { return *M_fluidMeshPart; }
    //!getter for the partitioned solid mesh
    const partitionMesh< mesh_Type >& solidMeshPart()             const { return *M_solidMeshPart; }

    //!getter for the fluid velocity FESpace
    const FESpace<mesh_Type, EpetraMap>& uFESpace()               const { return *M_uFESpace; }
    //!getter for the fluid pressure FESpace
    const FESpace<mesh_Type, EpetraMap>& pFESpace()               const { return *M_pFESpace; }
    //!getter for the solid displacement FESpace
    const FESpace<mesh_Type, EpetraMap>& dFESpace()               const { return *M_dFESpace; }
    //!getter for the harmonic extension solution FESpace
    const FESpace<mesh_Type, EpetraMap>& mmFESpace()              const { return *M_mmFESpace; }
    //!getter for the harmonic extension solution
    virtual const vector_Type& meshDisp()                         const { return M_meshMotion->disp(); }
    //!getter for the harmonic extension solution of the previous time step
    const         vector_Type& dispFluidMeshOld()                 const { return *M_dispFluidMeshOld; }
    //!getter for the mesh velocity
    virtual       vector_Type& veloFluidMesh()                          { return *M_veloFluidMesh; }
    //!getter for the mesh velocity increment (used for Newton FSI)
    //! \todo{try to remove this method}
    vector_Type& derVeloFluidMesh()                       { return *M_derVeloFluidMesh; }

    const dofInterface3DPtr_Type& dofFluidToStructure()             const { return M_dofFluidToStructure; }
    const dofInterface3DPtr_Type& dofStructureToSolid()             const { return M_dofStructureToSolid; }
    const dofInterface3DPtr_Type& dofStructureToHarmonicExtension() const { return M_dofStructureToHarmonicExtension; }
    const dofInterface3DPtr_Type& dofHarmonicExtensionToFluid()     const { return M_dofHarmonicExtensionToFluid; }

    boost::shared_ptr<EpetraMap>& fluidInterfaceMap()                   { return M_fluidInterfaceMap; }
    boost::shared_ptr<EpetraMap>& solidInterfaceMap()                   { return M_solidInterfaceMap; }

    //! Getter for the map of the variable used for the coupling
    virtual boost::shared_ptr<EpetraMap>& couplingVariableMap()      { return M_solidInterfaceMap; }

    //! Method to implement Robin boundary conditions on the external wall for the structure
    BCFunctionRobin& bcfRobinOuterWall()                                { return M_bcfRobinOuterWall; }

    bcVectorInterfacePtr_Type bcvStructureDisptoFluid()                 const { return M_bcvStructureDispToFluid; }
    bcVectorInterfacePtr_Type bcvStructureToFluid()                     const { return M_bcvStructureToFluid; }
    bcVectorInterfacePtr_Type bcvSolidLoadToStructure()                 const { return M_bcvSolidLoadToStructure; }
    bcVectorInterfacePtr_Type bcvFluidInterfaceDisp()                   const { return M_bcvFluidInterfaceDisp; }
    bcVectorInterfacePtr_Type bcvHarmonicExtensionVelToFluid()          const { return M_bcvHarmonicExtensionVelToFluid; }
    bcVectorInterfacePtr_Type bcvDerHarmonicExtensionVelToFluid()       const { return M_bcvDerHarmonicExtensionVelToFluid; }
    bcVectorInterfacePtr_Type bcvStructureDispToHarmonicExtension()     const { return M_bcvStructureDispToHarmonicExtension; }
    bcVectorInterfacePtr_Type bcvStructureDispToSolid()                 const { return M_bcvStructureDispToSolid; }
    bcVectorInterfacePtr_Type bcvDerStructureDispToSolid()              const { return M_bcvDerStructureDispToSolid; }
    bcVectorInterfacePtr_Type bcvFluidLoadToStructure()                 const { return M_bcvFluidLoadToStructure; }
    bcVectorInterfacePtr_Type bcvDerFluidLoadToStructure()              const { return M_bcvDerFluidLoadToStructure; }
    bcVectorInterfacePtr_Type bcvDerFluidLoadToFluid()                  const { return M_bcvDerFluidLoadToFluid; }
//     bcVectorInterfacePtr_Type bcvDerReducedFluidLoadToStructure()             { return M_bcvDerReducedFluidLoadToStructure; }
//     bcVectorInterfacePtr_Type bcvDerStructureAccToReducedFluid()              { return M_bcvDerStructureAccToReducedFluid; }

    //! Getter for the BCHandler of the fluid problem
    const fluidBchandlerPtr_Type& BCh_fluid()                       const { return M_BCh_u; }
    //! Getter for the BCHandler of the harmonic extension problem
    const fluidBchandlerPtr_Type& BCh_harmonicExtension()           const { return M_BCh_mesh; }
    //! Getter for the BCHandler of the linearized fluid problem (to be used in Newton for the partitioned FSI)
    const fluidBchandlerPtr_Type& BCh_du()                          const { return M_BCh_du; }
    //! Getter for the BCHandler of the linearized inverse of the fluid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const fluidBchandlerPtr_Type& BCh_du_inv()                      const { return M_BCh_du_inv; }
    //! Getter for the BCHandler of the solid problem
    const solidBchandlerPtr_Type& BCh_solid()                       const { return M_BCh_d; }
    //! Getter for the BCHandler of the linearized solid problem
    const solidBchandlerPtr_Type& BCh_dz()                          const { return M_BCh_dz; }
    //! Getter for the BCHandler of the linearized inverse of the solid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const solidBchandlerPtr_Type& BCh_dz_inv()                      const { return M_BCh_dz_inv; }

    //! gets the solution vector by reference
    virtual const vector_Type& solution()                      const { return *M_lambda; }

    //! gets a pointer to the solution vector by reference
    virtual vectorPtr_Type& solutionPtr()                               { return M_lambda; }

    //! gets the solid displacement by copy
    virtual void getSolidDisp( vector_Type& soliddisp )                 { soliddisp = M_solid->getDisplacement(); }

    //! gets the solid velocity by copy
    virtual void getSolidVel( vector_Type& solidvel )                   { solidvel = M_solid->getVelocity(); }

    //! gets the fluid velocity and pressure by copy
    virtual void getFluidVelAndPres( vector_Type& sol )                 { sol = *M_fluid->solution(); }

    //! Getter for the right hand side
    vectorPtr_Type const& getRHS       ( ) const {return M_rhs;}


    //@}



    /** @name  Set Functions
     */
    //@{
    //! Setter for the local and  world communicators
    /**
       The communicator can be different depending on which type of subdomain we are considering
     */
    void setComm( const commPtr_Type& comm, const commPtr_Type& worldComm );

    //! Setter for the FSI data
    void setData( const dataPtr_Type& data ) { M_data = data; }

    //! Setter for the fluid and geometry problems
    void setFluid( const fluidPtr_Type& fluid, const meshMotionPtr_Type& meshmotion );
    //! Setter for the solid problem
    void setSolid( const solidPtr_Type& solid );

    //!Setter for the "fluid" flag
    void setFluid( const bool& isFluid ) { M_isFluid = isFluid; }
    //!Setter for the "solid" flag
    void setSolid( const bool& isSolid ) { M_isSolid = isSolid; }

    //!Setter for the "linear fluid" flag
    void setLinearFluid( const bool& linFluid ) { M_linearFluid = linFluid; }
    //!Setter for the "linear solid" flag
    void setLinearSolid( const bool& linSolid ) { M_linearSolid = linSolid; }

    void setFluidLeader( const int& fluidLeader ) { M_fluidLeader = fluidLeader; }
    void setSolidLeader( const int& solidLeader ) { M_solidLeader = solidLeader; }

//     void setBC( fluidBchandlerPtr_Type& bc_u, solidBchandlerPtr_Type& bc_d, fluidBchandlerPtr_Type& bc_m );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setFluidBC     ( const fluidBchandlerPtr_Type& bc_fluid );
    //!Setter for the BCHandler of the linearized fluid problem (to be used in segregated Newton FSI)
    void setLinFluidBC          ( const fluidBchandlerPtr_Type& bc_dfluid )     { M_BCh_du     = bc_dfluid; }
    //!Setter for the BCHandler of the inverse linearized fluid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
     */
    void setInvLinFluidBC       ( const fluidBchandlerPtr_Type& bc_dfluid_inv ) { M_BCh_du_inv = bc_dfluid_inv; }
    //!Setter for the BCHandler of the gerometry problem (to be used in segregated Newton FSI)
    void setHarmonicExtensionBC ( const fluidBchandlerPtr_Type& bc_he );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setSolidBC     ( const solidBchandlerPtr_Type& bc_solid );
    //!Setter for the BCHandler of the linearized solid problem (to be used in segregated Newton FSI)
    void setLinSolidBC          ( const solidBchandlerPtr_Type& bc_dsolid )     { M_BCh_dz     = bc_dsolid; }
    //!Setter for the BCHandler of the inverse linearized solid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
    */
    void setInvLinSolidBC       ( const solidBchandlerPtr_Type& bc_dsolid_inv ) { M_BCh_dz_inv = bc_dsolid_inv; }

    //! Setter for the interface displacement (partitioned according to the fluid)
    void setLambdaFluid         ( const vector_Type& lambda );
    //! Setter for the interface displacement (partitioned according to the solid)
    void setLambdaSolid         ( const vector_Type& lambda );

    //! Setter for the solid interface displacement at the previous time step
    /**\todo{see if we can remove these}*/
    void setLambdaSolidOld      ( const vector_Type& lambda );
    //! Setter for the solid interface velocity at the previous time step
    /**\todo{see if we can remove these}*/
    void setLambdaDotSolid      ( const vector_Type& lambda );

    //!Setter for the fluid interface stress
    void setSigmaFluid          ( const vector_Type& sigma );
    //!Setter for the solid interface stress
    void setSigmaSolid          ( const vector_Type& sigma );
    //\todo{try to remove this}
    void setMinusSigmaFluid     ( const vector_Type& sigma );

    //! Setter for the Robin coefficient of the Robin--Neumann coupling scheme (as a BCFunction)
    void setAlphafbcf  ( const bcFunction_Type& alphafbcf );

//     void setMpi     (bool mpi  ){M_mpi      = mpi;}
//     void setFluidMpi(bool fluid){M_isFluidMpi = fluid;}
//     void setSolidMpi(bool solid){M_issolidMpi = solid;}

//     bool mpi(){return M_mpi;}

//     void setRobinOuterWall                   ( function_Type const& dload,
//                                                function_Type const& E ) { M_bcfRobinOuterWall.setFunctions_Robin(dload, E); }
    void setStructureDispToHarmonicExtension ( const vector_Type& disp,  UInt type = 0 );
    void setStructureToFluid                 ( const vector_Type& vel,   UInt type = 0 );
    void setStructureDispToFluid             ( const vector_Type& vel,   UInt type = 0 );
    void setStructureDispToSolid             ( const vector_Type& disp,  UInt type = 0 );
    void setDerStructureDispToSolid          ( const vector_Type& ddisp, UInt type = 0 );
    void setSolidLoadToStructure             ( const vector_Type& load,  UInt type = 0 );
    void setHarmonicExtensionVelToFluid      ( const vector_Type& vel,   UInt type = 0 );
    void setDerHarmonicExtensionVelToFluid   ( const vector_Type& dvel,  UInt type = 0 );
    //void setFluidInterfaceDisp               ( const vector_Type& disp,  UInt type = 0 );
    void setFluidLoadToStructure             ( const vector_Type& load,  UInt type = 0 );
    void setDerFluidLoadToStructure          ( const vector_Type& dload, UInt type = 0 );
    void setDerFluidLoadToFluid              ( const vector_Type& dload, UInt type = 0 );
    void setRobinOuterWall                   ( const function_Type& dload, const function_Type& E);

    //! sets the solution vector by reference
    virtual void setSolutionPtr              ( const vectorPtr_Type& sol )      { M_lambda = sol; }

    //! sets the solution vector by copy
    virtual void setSolution                 ( const vector_Type& solution )    { M_lambda.reset( new vector_Type( solution ) ); }
    //! Initialize the fluid solution given a vector (calls setSolution)
    virtual void initialize                  ( vectorPtr_Type u0)               { setSolution(*u0); }


    //! sets the solution time derivative vector by copy
    //void setSolutionDerivative( vectorPtr_Type& lambdaDot ) { M_lambdaDot = lambdaDot; }

    //! Setter for the time derivative of the interface displacement
    void setSolutionDerivative( const vector_Type& solutionDerivative )         { M_lambdaDot.reset( new vector_Type( solutionDerivative ) ); }

    //!\todo{see if can be marked ad deprected}
    virtual void  setRestarts( bool /*restarts*/ ) { /*M_restarts = restarts;*/ }
    //@}

//     void setDerReducedFluidLoadToStructure   ( vector_Type &dload, UInt type = 0 );
//     void setDerStructureAccToReducedFluid    ( vector_Type &acc,   UInt type = 0 );

//     void setReducedLinFluidBC                (fluidBchandlerPtr_Type bc_dredfluid);
//     void setInvReducedLinFluidBC             (fluidBchandlerPtr_Type bc_invdredfluid);

//     void setBCh_fluid             ( const fluidBchandlerPtr_Type& BCh_fluid)       { M_BCh_u      = BCh_fluid; }
//     void setBCh_HarmonicExtension ( const fluidBchandlerPtr_Type& BCh_mesh)        { M_BCh_mesh   = BCh_mesh; }
//     void setBCh_fluidDer          ( const fluidBchandlerPtr_Type& BCh_fluidDer)    { M_BCh_du     = BCh_fluidDer; }
//     void setBCh_fluidDerInv       ( const fluidBchandlerPtr_Type& BCh_fluidDerInv) { M_BCh_du_inv = BCh_fluidDerInv; }
//     void setBCh_solid             ( const solidBchandlerPtr_Type& BCh_solid)       { M_BCh_d      = BCh_solid; }
//     void setBCh_solidDer          ( const solidBchandlerPtr_Type& BCh_solidDer)    { M_BCh_dz     = BCh_solidDer; }
//     void setBCh_solidDerInv       ( const solidBchandlerPtr_Type& BCh_solidDerInv) { M_BCh_dz_inv = BCh_solidDerInv; }

    //@}


protected:

    //virtual void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);
    //!@name Protected Methods
    //@{
    //!initailize the variables
    /**
       instantiates the pointers which are used in the segregated FSI solvers. Reimplemented in the Monolithic class.
       \param dorder: unused parameter
     */
    virtual void variablesInit( const std::string& dOrder );

    //!Interpolates the mesh motion dofs on the fluid
    /**
       The order of the spatial approximation depends on this method: when the mesh motion approximation is first order
       in space the overall approximation is of the first order even if the fluid is solved with hicher order FEs.
       Calls the interpolateVelocity method
     */
    void transferMeshMotionOnFluid( const vector_Type& _vec1,
                                    vector_Type& _vec2 );

    //! Interpolates mesh motion into velocity
    /**
       Interpolates a vector with the map of the harmonic extension into one with the map of the fluid velocity
     */
    void interpolateVelocity(const vector_Type& _vec1,
                             vector_Type& _vec2);


    //!Interpolates to vectors on the interface
    /**
       The two vectors can have different numeration, for different discretizations this method is not tested.
    */
    void interpolateInterfaceDofs(const FESpace<mesh_Type, EpetraMap>& _fespace1,
                                  const vector_Type&                   _vec1,
                                  const FESpace<mesh_Type, EpetraMap>& _fespace2,
                                  vector_Type&                         _vec2,
                                  dofInterface3DPtr_Type&                _dofInterface);
    //@}

    //!@name Protected Attributes
    //@{
    boost::shared_ptr<mesh_Type>                      M_mesh;

    boost::shared_ptr<FESpace<mesh_Type, EpetraMap> > M_uFESpace;
    boost::shared_ptr<FESpace<mesh_Type, EpetraMap> > M_pFESpace;
    boost::shared_ptr<FESpace<mesh_Type, EpetraMap> > M_dFESpace;
    boost::shared_ptr<FESpace<mesh_Type, EpetraMap> > M_mmFESpace;

    boost::shared_ptr<mesh_Type>                      M_fluidMesh;
    boost::shared_ptr<mesh_Type>                      M_solidMesh;

    boost::shared_ptr<partitionMesh< mesh_Type > >    M_fluidMeshPart;
    boost::shared_ptr<partitionMesh< mesh_Type > >    M_solidMeshPart;

    fluidBchandlerPtr_Type                              M_BCh_u;
    solidBchandlerPtr_Type                              M_BCh_d;
    fluidBchandlerPtr_Type                              M_BCh_mesh;

    // interface operators BCs
    fluidBchandlerPtr_Type                              M_BCh_du;
    fluidBchandlerPtr_Type                              M_BCh_du_inv;

    solidBchandlerPtr_Type                              M_BCh_dz;
    solidBchandlerPtr_Type                              M_BCh_dz_inv;

    fluidBchandlerPtr_Type                              M_BCh_dp;
    fluidBchandlerPtr_Type                              M_BCh_dp_inv;

    fluidPtr_Type                                        M_fluid;
    solidPtr_Type                                        M_solid;
    meshMotionPtr_Type                                   M_meshMotion;

//     fluidLinPtr_Type                                     M_fluidLin;
//     solidLinPtr_Type                                     M_solidLin;

    boost::shared_ptr<BdfT<vector_Type> >             M_bdf;

    GetPot                                            M_dataFile;

    boost::shared_ptr<DataMesh>                       M_dataMeshFluid;
    boost::shared_ptr<DataMesh>                       M_dataMeshSolid;

    dataPtr_Type                                      M_data;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMap;

    //!\todo{kill this attribute}
    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMapOnZero;
    //!\todo{kill this attribute}
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMapOnZero;

    dofInterface3DPtr_Type                              M_dofFluidToStructure; // Needed
//     dofInterface3DPtr_Type                              M_dofSolidToFluid;
    dofInterface3DPtr_Type                              M_dofStructureToFluid; // Needed
    dofInterface3DPtr_Type                              M_dofStructureToSolid; // Needed to construct M_bcvStructureDispToSolid
    dofInterface3DPtr_Type                              M_dofStructureToHarmonicExtension; // Needed to construct interface maps
    dofInterface3DPtr_Type                              M_dofHarmonicExtensionToFluid; // Needed to construct M_bcvStructureToFluid
//     dofInterface3DPtr_Type                              M_dofStructureToReducedFluid;
//     dofInterface3DPtr_Type                              M_dofReducedFluidToStructure;

    dofInterface2DPtr_Type                              M_dofFluid;
    dofInterface2DPtr_Type                              M_dofSolid;
    dofInterface2DPtr_Type                              M_dofFluidInv;
    dofInterface2DPtr_Type                              M_dofSolidInv;

    bcVectorInterfacePtr_Type                               M_bcvFluidInterfaceDisp;
    bcVectorInterfacePtr_Type                               M_bcvFluidLoadToStructure;
    bcVectorInterfacePtr_Type                               M_bcvSolidLoadToStructure;
    bcVectorInterfacePtr_Type                               M_bcvStructureToFluid;
    bcVectorInterfacePtr_Type                               M_bcvStructureDispToFluid;
    bcVectorInterfacePtr_Type                               M_bcvStructureDispToSolid;
    bcVectorInterfacePtr_Type                               M_bcvStructureDispToHarmonicExtension;
    bcVectorInterfacePtr_Type                               M_bcvHarmonicExtensionVelToFluid;
//     bcVectorInterfacePtr_Type                               M_bcvStructureToReducedFluid;
//     bcVectorInterfacePtr_Type                               M_bcvReducedFluidToStructure;

    bcVectorInterfacePtr_Type                               M_bcvDerHarmonicExtensionVelToFluid;
    bcVectorInterfacePtr_Type                               M_bcvDerFluidLoadToStructure;
    bcVectorInterfacePtr_Type                               M_bcvDerFluidLoadToFluid;
    bcVectorInterfacePtr_Type                               M_bcvDerStructureDispToSolid;
    BCFunctionRobin                                   M_bcfRobinOuterWall;

//     bcVectorInterfacePtr_Type                               M_bcvDerReducedFluidLoadToStructure;
//     bcVectorInterfacePtr_Type                               M_bcvDerStructureAccToReducedFluid;

    vectorPtr_Type                    M_lambdaFluid;
    vectorPtr_Type                    M_lambdaFluidRepeated;
    vectorPtr_Type                    M_lambda;
    vectorPtr_Type                    M_lambdaDot;


    vectorPtr_Type                    M_un;
    vectorPtr_Type                    M_rhs;
    vectorPtr_Type                    M_Alphaf;

    Real                                              M_AlphafCoef;
    //\todo{try to set as deprecated}
    Real                                              M_betamedio;

    UInt                                              M_fluxes;

    commPtr_Type                                      M_epetraComm;
    commPtr_Type                                      M_epetraWorldComm;

    //@}
private:

    //!@name Private Methods
    //@{
    //!Private Copy Constructor
    FSIOperator( const FSIOperator& copy){}
    //@}

    //! @name Private Attributes
    //@{
    // displacement on the interface
    vectorPtr_Type                    M_lambdaSolid;
    vectorPtr_Type                    M_lambdaSolidRepeated;

    vectorPtr_Type                    M_lambdaSolidOld;
    vectorPtr_Type                    M_lambdaDotSolid;
    vectorPtr_Type                    M_lambdaDotSolidRepeated;

    vectorPtr_Type                    M_sigmaFluid;
    vectorPtr_Type                    M_sigmaSolid;

    vectorPtr_Type                    M_sigmaFluidRepeated;
    vectorPtr_Type                    M_sigmaSolidRepeated;

    //\todo{try to remove}
    vectorPtr_Type                    M_minusSigmaFluid;
    //\todo{try to remove}
    vectorPtr_Type                    M_minusSigmaFluidRepeated;

    vectorPtr_Type                    M_dispFluidMeshOld;
    vectorPtr_Type                    M_veloFluidMesh;
    vectorPtr_Type                    M_derVeloFluidMesh;

    //\todo{try to set as deprecated}
    bool                                              M_mpi;

    bool                                              M_isFluid;
    bool                                              M_isSolid;

    bool                                              M_linearFluid;
    bool                                              M_linearSolid;

    int                                               M_fluidLeader;
    int                                               M_solidLeader;
    //@}

};

} // Namespace LifeV

#endif /* FSIOPERATOR_H */
