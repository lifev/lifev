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
    @file
    @brief Pure virtual operator class for FSI solvers

    @author Simone Deparis <simone.deparis@epfls.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010
    \include ../../tools/bibliography/fluidstructure.dox
    \include ../../tools/bibliography/newton.dox

    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton method implemente in NonlinearRichardson. The fluid
    and structure classes are member of this class and different formulations (e.g. Monolithic \ref CDFQ10 , segregated
    Newton \ref FM05 , Dirichlet--Neumann \ref DDFQ06 , Robin Neumann \ref BNV08 )

 */
#ifndef FSIOPERATOR_H
#define  FSIOPERATOR_H

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>
#include <life/lifesolver/DataFSI.hpp>
#include <life/lifesolver/OseenShapeDerivative.hpp>
//#include <life/lifesolver/NonLinearVenantKirchhofSolver.hpp>
#include <life/lifesolver/LinearVenantKirchhofSolver.hpp>
#include <life/lifesolver/HarmonicExtensionSolver.hpp>

#include <life/lifemesh/dataMesh.hpp>
#include <life/lifemesh/regionMesh3D.hpp>
#include <life/lifefilters/HDF5Filter3DMesh.hpp>

#include <life/lifealg/generalizedAitken.hpp>

#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifefem/bdf_template.hpp>
//#include <life/lifefem/dof.hpp>
#include <life/lifefem/FESpace.hpp>
//#include <life/lifefem/refFE.hpp>
//#include <life/lifesolver/fixedPointBase.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#define FLUID 1
#define SOLID 0

namespace LifeV
{

/*!
    \class FSIOperator
    \brief Fluid-Structure Interface operator class

    \include ../../tools/bibliography/fluidstructure.dox
    \include ../../tools/bibliography/newton.dox
    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton method implemente in NonlinearRichardson. The fluid
    and structure classes are member of this class and different formulations (e.g. Monolithic \ref CDFQ10 , segregated
    Newton \ref FM05 , Dirichlet--Neumann \ref DDFQ06 , Robin Neumann \ref BNV08 )

*/
class FSIOperator
{
public:

    /** @name Typedefs
     */
    //@{

    //typedef VenantKirchhofSolver   < RegionMesh3D_ALE<LinearTetra> > solid_raw_type;
    //typedef NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid_raw_type;

    typedef RegionMesh3D<LinearTetra>              mesh_type;

#ifdef HAVE_HDF5
    typedef HDF5Filter3DMesh<mesh_type>            mesh_filtertype;
#endif
    typedef OseenShapeDerivative   <mesh_type>     fluid_raw_type;
    //typedef NonLinearVenantKirchhofSolver   <mesh_type>     solid_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solid_raw_type;
    typedef HarmonicExtensionSolver<mesh_type>     meshmotion_raw_type;

    typedef OseenShapeDerivative   <mesh_type>     fluidlin_raw_type;
    //typedef NonLinearVenantKirchhofSolver   <mesh_type>     solidlin_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solidlin_raw_type;

    typedef boost::shared_ptr<fluid_raw_type>      fluid_type;
    typedef boost::shared_ptr<solid_raw_type>      solid_type;
    typedef boost::shared_ptr<meshmotion_raw_type> meshmotion_type;

    typedef boost::shared_ptr<fluidlin_raw_type>   fluidlin_type;
    typedef boost::shared_ptr<solidlin_raw_type>   solidlin_type;

    typedef fluid_raw_type::vector_type            vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;

    typedef fluid_raw_type::source_type            fluid_source_type;
    typedef solid_raw_type::source_Type            solid_source_type;

    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )> function_type;

    typedef Real ( *bc_function_type ) ( const Real&, const Real&,
                                         const Real&, const Real&, const ID& );

    typedef boost::shared_ptr<DofInterface3Dto3D>  dof_interface_type3D;
    typedef boost::shared_ptr<DofInterface3Dto2D>  dof_interface_type2D;
    typedef boost::shared_ptr<BCVectorInterface>   bc_vector_interface;

    typedef fluid_raw_type::bchandler_type         fluid_bchandler_type;
    typedef fluid_raw_type::bchandler_raw_type     fluid_bchandler_raw_type;

    typedef BCHandler                              solid_bchandler_raw_type;
    typedef boost::shared_ptr<solid_bchandler_raw_type>  solid_bchandler_type;

    typedef DataFSI                                data_Type;
    typedef boost::shared_ptr<data_Type>           data_PtrType;

    typedef std::map<ID, ID>::const_iterator                        Iterator;


    typedef singleton<factory<FSIOperator, std::string> >           FSIFactory;

    typedef Displayer::comm_PtrType                                 comm_PtrType;

//     typedef boost::shared_ptr<reducedLinFluid>    quasi_newton_type;

    //@}



    /** @name Constructors, destructor
     */
    //@{

    FSIOperator();

    virtual ~FSIOperator();

    //@}


    /** @name Virtual Methods
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
    void partitionMeshes( mesh_filtertype& fluidMeshFilter, mesh_filtertype& solidMeshFilter  );
#endif
    //! sets up the correspondences between the fluid and structure degrees of freedom across the interface.
    /**
       This method introduces a non scalable loop, in DOFInterface3Dto3D. It is preferable to avoid it for massively
       parallel computetions, using the offline partitioner.
     */
    virtual void setupDOF();

#ifdef HAVE_HDF5
    //!reads from HDF5 file the correspondences between the fluid and structure degrees of freedom across the interface.
    /*!still not implemented for all the FSI formulations*/
    virtual void setupDOF( mesh_filtertype& /*filterMesh*/ ) {}
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
       \todo{replace Real with Real& }
     */
    virtual void solveJac( vector_type&       muk,
                           const vector_type& res,
                           const Real         linearRelTol ) = 0;

    //! Evaluates the nonlinear residual of the FSI system
    /**
       The implementation of this method also depends on the child classes, though it does not characterize them.
       \param res:  residual vector  to be computed
       \param disp: current unknown solution
       \param iter: nonlinear iteration counter. The part of th rhs related to the time discretization is computed only for iter=0
     */
    virtual void evalResidual( vector_type&        res,
                               const vector_type&  disp,
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
    virtual void initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::fluid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::fluid_type::value_type::Function const& /*df0=FSIOperator::solid_type::value_type::Function()*/)
    {
        Debug( 6220 ) << "FSIOperator:: solid init \n";
        if (this->isSolid())
            solid().initialize(d0, w0, w0);
        Debug( 6220 ) << "FSIOperator:: fluid init \n";
        if (this->isFluid())
            fluid().initialize(u0, p0);
    }

    //! getter for the fluid velocity
    virtual const vector_ptrtype& un() { return M_un; }

    //! Returns the number of imposed fluxes
    virtual UInt imposeFlux();

    // \todo{kill this method?}
    virtual void mergeBCHandlers() {}

    //@}

    //!\todo{kill this method}
    virtual void setFluxBC             (fluid_bchandler_type /*bc_fluid*/) {}

    //!\todo{kill this method}
    virtual void setRobinBC             (fluid_bchandler_type /*bc_solid*/) {}



    //!@name MONOLITHIC Solver methods - Implemented there
    //{@
    virtual void iterateMesh( const vector_type& /*disp*/ )     { assert(false); }
    virtual void setupBDF( const vector_type& /*u0*/ ) { }
    virtual void updateRHS() {}
    virtual void applyBoundaryConditions() {}
    //@}

    //    static VenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >*    createNonLinearStructure(){ return new NonLinearVenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >(); }

    //!@name Factory Methods
    //@{
    //! Factory method for the linear elasticity solver
    static VenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >*    createLinearStructure() { return new LinearVenantKirchhofSolver< FSIOperator::mesh_type, SolverTrilinos >(); }
    //@}

    //!@name Public Member Functions
    //@{
    //!Initializes the BDF which should handle the fluid time discretization
    /**
       \todo{a general time advancing class should be used everywhere}
     */
    void initializeBDF( const vector_type& un );

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
    void initializeFluid( const vector_type& velAndPressure,
                          const vector_type& displacement );

    //! initializes the solid solver with vectors
    /**
       \param displacement: initial vector containing the structure displacement
       \param velocity: initial vector containing the velocity, used for the initialization of the Newmark scheme
     */
    void initializeSolid( vector_ptrtype displacement,
                          vector_ptrtype velocity );

    //!\todo{kill this method}
    void updateJacobian ( const vector_type& sol, const int& iter );

    //!moves the mesh using the solution of the harmonic extension equation
    /**
       \param disp displacement of the mesh, must be the difference between the current solution of the HE problem and the one at the previous time step.
     */
    void moveMesh       ( const vector_type& disp );

//     void solveLinearFluid();
//     void solveLinearSolid();

    //!Method to import an EpetraVector defined on the fluid map (i.e. with the fluid numeration of the dofs) to the interface
    /**
       Note that the output vector will have the solid numeration on the interface! By default in fact the vectors on the
       FSI interface in LifeV are numerated according to the solid.
     */
    void transferFluidOnInterface( const vector_type& _vec1, vector_type& _vec2 );

    //works in serial but no yet in parallel
    void transferSolidOnFluid    ( const vector_type& _vec1, vector_type& _vec2 );

    //!Method to import an EpetraVector defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       The output vector has the solid numeration of the dofs and is partitioned according to the solid partition. This method is not used in the monolithic solvers.
     */
    void transferSolidOnInterface( const vector_type& _vec1, vector_type& _vec2 );
//     void transferInterfaceOnFluid( const vector_type& _vec1, vector_type& _vec2 );

    //!Method to import an EpetraVector defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       the output vector have the numeration of the solid, as in transferSolidOnInterface, but is partitioned according to the fluid! This method is not used in the monolithic solvers.
     */
    void transferInterfaceOnSolid( const vector_type& _vec1, vector_type& _vec2 );

    //! calls bcManage for a vector
    /**
       calls bcManage for the vector rhs and the BCHandler bch
     */
    void bcManageVectorRHS( const fluid_bchandler_type& bch, vector_type& rhs );

    //! Method to set the Robin vector coefficient of the Robin--Neumann coupling scheme (as a constant vector vector)
    void setAlphaf     () { M_Alphaf->getEpetraVector().PutScalar( M_AlphafCoef ); }
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
    void setRHS                              ( vector_ptrtype& rhs ) {M_rhs = rhs;}


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

//     vector_type & displacement()                                        { return *M_lambdaSolid; }
//     vector_type & displacementOld()                                     { return *M_lambdaSolidOld; }
//     vector_type & residual();
//     vector_type const & velocity()                                const { return *M_lambdaDotSolid; }
//     vector_type & residualFSI();

    const vector_type& lambdaFluid()                              const { return *M_lambdaFluid; }
    const vector_type& lambdaSolid()                              const { return *M_lambdaSolid; }
    const vector_type& lambdaSolidOld()                           const { return *M_lambdaSolidOld; }
    const vector_type& lambdaDotSolid()                           const { return *M_lambdaDotSolid; }
    const vector_type& sigmaFluid()                               const { return *M_sigmaFluid; }
    const vector_type& sigmaSolid()                               const { return *M_sigmaSolid; }

    const vector_type& lambdaFluidRepeated()                      const { return *M_lambdaFluidRepeated; }
    const vector_type& lambdaSolidRepeated()                      const { return *M_lambdaSolidRepeated; }
    const vector_type& lambdaDotSolidRepeated()                   const { return *M_lambdaDotSolidRepeated; }
    const vector_type& sigmaFluidRepeated()                       const { return *M_sigmaFluidRepeated; }
    const vector_type& sigmaSolidRepeated()                       const { return *M_sigmaSolidRepeated; }

    //!\todo{remove this method}
    // now residual is Ax-b, but we should decide for b-Ax. In the meantime, we need b-Ax:
    const vector_type& minusSigmaFluid()                          const { return *M_minusSigmaFluid; }
    //!\todo{remove this method}
    const vector_type& minusSigmaFluidRepeated()                  const { return *M_minusSigmaFluidRepeated; }

    //!coefficient for the Robin--Neumann coupling scheme
    vector_type&       Alphaf()                                   const { return *M_Alphaf;}

    comm_PtrType worldComm()                                      const { return M_epetraWorldComm; }

    bool isFluid()                                                const { return M_isFluid; }
    bool isSolid()                                                const { return M_isSolid; }

    bool isLinearFluid()                                          const { return M_linearFluid; }
    bool isLinearSolid()                                          const { return M_linearSolid; }

    int getFluidLeaderId()                                        const { return M_fluidLeader; }
    int getSolidLeaderId()                                        const { return M_solidLeader; }

    //! Getter for the fluid solver
    const fluid_raw_type& fluid()                                 const { return *M_fluid; }
    //! Getter for the solid solver
    const solid_raw_type& solid()                                 const { return *M_solid; }
    //! Getter for the harmonic extension solver
    const meshmotion_raw_type& meshMotion()                       const { return *M_meshMotion; }

    //! Getter-Setter for the fluid solver
    /** \todo{mark as deprecated}*/
    fluid_type::value_type& fluid()                                     { return *M_fluid; }
    //! Getter-Setter for the solid solver
    /** \todo{mark as deprecated}*/
    solid_type::value_type& solid()                                     { return *M_solid; }
    //! Getter-Setter for the mesh motion solver
    /** \todo{mark as deprecated}*/
    meshmotion_type::value_type& meshMotion()                           { return *M_meshMotion; }
//     fluidlin_type::value_type& fluidLin()                               { return *M_fluidLin; }
//     solidlin_type::value_type& solidLin()                               { return *M_solidLin; }

    //!getter for the FSI data container
    const data_Type& data()                                       const { return *M_data; }
    //!getter for the fluid data container
    const data_Type::dataFluid_PtrType& dataFluid()               const { return M_data->dataFluid(); }
    //!getter for the solid data container
    const data_Type::dataSolid_PtrType& dataSolid()               const { return M_data->dataSolid(); }

    //!getter for the unpartitioned fluid mesh
    mesh_type& fluidMesh()                                        const { return *M_fluidMesh; }
    //!getter for the unpartitioned solid mesh
    mesh_type& solidMesh()                                        const { return *M_solidMesh; }

    // const mesh_type& fluidMesh()                                  const { return *M_fluidMesh; }
    // const mesh_type& solidMesh()                                  const { return *M_solidMesh; }

    //!getter for the partitioned fluid mesh
    const partitionMesh< mesh_type >& fluidMeshPart()             const { return *M_fluidMeshPart; }
    //!getter for the partitioned solid mesh
    const partitionMesh< mesh_type >& solidMeshPart()             const { return *M_solidMeshPart; }

    //!getter for the fluid velocity FESpace
    const FESpace<mesh_type, EpetraMap>& uFESpace()               const { return *M_uFESpace; }
    //!getter for the fluid pressure FESpace
    const FESpace<mesh_type, EpetraMap>& pFESpace()               const { return *M_pFESpace; }
    //!getter for the solid displacement FESpace
    const FESpace<mesh_type, EpetraMap>& dFESpace()               const { return *M_dFESpace; }
    //!getter for the harmonic extension solution FESpace
    const FESpace<mesh_type, EpetraMap>& mmFESpace()              const { return *M_mmFESpace; }
    //!getter for the harmonic extension solution
    virtual const vector_type& meshDisp()                         const { return M_meshMotion->disp(); }
    //!getter for the harmonic extension solution of the previous time step
    const         vector_type& dispFluidMeshOld()                 const { return *M_dispFluidMeshOld; }
    //!getter for the mesh velocity
    virtual       vector_type& veloFluidMesh()                          { return *M_veloFluidMesh; }
    //!getter for the mesh velocity increment (used for Newton FSI)
    //! \todo{try to remove this method}
    vector_type& derVeloFluidMesh()                       { return *M_derVeloFluidMesh; }

    const dof_interface_type3D& dofFluidToStructure()             const { return M_dofFluidToStructure; }
    const dof_interface_type3D& dofStructureToSolid()             const { return M_dofStructureToSolid; }
    const dof_interface_type3D& dofStructureToHarmonicExtension() const { return M_dofStructureToHarmonicExtension; }
    const dof_interface_type3D& dofHarmonicExtensionToFluid()     const { return M_dofHarmonicExtensionToFluid; }

    boost::shared_ptr<EpetraMap>& fluidInterfaceMap()                   { return M_fluidInterfaceMap; }
    boost::shared_ptr<EpetraMap>& solidInterfaceMap()                   { return M_solidInterfaceMap; }

    //! Getter for the map of the variable used for the coupling
    virtual boost::shared_ptr<EpetraMap>& getCouplingVariableMap()      { return M_solidInterfaceMap; }

    //! Method to implement Robin boundary conditions on the external wall for the structure
    BCFunctionMixte& bcfMixteOuterWall()                                { return M_bcfMixteOuterWall; }

    bc_vector_interface bcvStructureDisptoFluid()                 const { return M_bcvStructureDispToFluid; }
    bc_vector_interface bcvStructureToFluid()                     const { return M_bcvStructureToFluid; }
    bc_vector_interface bcvSolidLoadToStructure()                 const { return M_bcvSolidLoadToStructure; }
    bc_vector_interface bcvFluidInterfaceDisp()                   const { return M_bcvFluidInterfaceDisp; }
    bc_vector_interface bcvHarmonicExtensionVelToFluid()          const { return M_bcvHarmonicExtensionVelToFluid; }
    bc_vector_interface bcvDerHarmonicExtensionVelToFluid()       const { return M_bcvDerHarmonicExtensionVelToFluid; }
    bc_vector_interface bcvStructureDispToHarmonicExtension()     const { return M_bcvStructureDispToHarmonicExtension; }
    bc_vector_interface bcvStructureDispToSolid()                 const { return M_bcvStructureDispToSolid; }
    bc_vector_interface bcvDerStructureDispToSolid()              const { return M_bcvDerStructureDispToSolid; }
    bc_vector_interface bcvFluidLoadToStructure()                 const { return M_bcvFluidLoadToStructure; }
    bc_vector_interface bcvDerFluidLoadToStructure()              const { return M_bcvDerFluidLoadToStructure; }
    bc_vector_interface bcvDerFluidLoadToFluid()                  const { return M_bcvDerFluidLoadToFluid; }
//     bc_vector_interface bcvDerReducedFluidLoadToStructure()             { return M_bcvDerReducedFluidLoadToStructure; }
//     bc_vector_interface bcvDerStructureAccToReducedFluid()              { return M_bcvDerStructureAccToReducedFluid; }

    //! Getter for the BCHandler of the fluid problem
    const fluid_bchandler_type& BCh_fluid()                       const { return M_BCh_u; }
    //! Getter for the BCHandler of the harmonic extension problem
    const fluid_bchandler_type& BCh_harmonicExtension()           const { return M_BCh_mesh; }
    //! Getter for the BCHandler of the linearized fluid problem (to be used in Newton for the partitioned FSI)
    const fluid_bchandler_type& BCh_du()                          const { return M_BCh_du; }
    //! Getter for the BCHandler of the linearized inverse of the fluid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const fluid_bchandler_type& BCh_du_inv()                      const { return M_BCh_du_inv; }
    //! Getter for the BCHandler of the solid problem
    const solid_bchandler_type& BCh_solid()                       const { return M_BCh_d; }
    //! Getter for the BCHandler of the linearized solid problem
    const solid_bchandler_type& BCh_dz()                          const { return M_BCh_dz; }
    //! Getter for the BCHandler of the linearized inverse of the solid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const solid_bchandler_type& BCh_dz_inv()                      const { return M_BCh_dz_inv; }

    //! gets the solution vector by reference
    virtual const vector_type& getSolution()                      const { return *M_lambda; }

    //! gets a pointer to the solution vector by reference
    virtual vector_ptrtype& solutionPtr()                               { return M_lambda; }

    //! gets the solid displacement by copy
    virtual void getSolidDisp( vector_type& soliddisp )                 { soliddisp = M_solid->disp(); }

    //! gets the solid velocity by copy
    virtual void getSolidVel( vector_type& solidvel )                   { solidvel = M_solid->vel(); }

    //! gets the fluid velocity and pressure by copy
    virtual void getFluidVelAndPres( vector_type& sol )                 { sol = *M_fluid->solution(); }

    //! Getter for the right hand side
    vector_ptrtype const& getRHS       ( ) const {return M_rhs;}


//     quasi_newton_type getReducedLinFluid()                              { return M_reducedLinFluid; }
//     UInt reducedFluid()                                                 { return M_reducedFluid; }

    //@}



    /** @name  Set Functions
     */
    //@{
    //! Setter for the local and  world communicators
    /**
       The communicator can be different depending on which type of subdomain we are considering
     */
    void setComm( const comm_PtrType& comm, const comm_PtrType& worldComm );

    //! Setter for the FSI data
    void setData( const data_PtrType& data ) { M_data = data; }

    //! Setter for the fluid and geometry problems
    void setFluid( const fluid_type& fluid, const meshmotion_type& meshmotion );
    //! Setter for the solid problem
    void setSolid( const solid_type& solid );

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

//     void setBC( fluid_bchandler_type& bc_u, solid_bchandler_type& bc_d, fluid_bchandler_type& bc_m );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setFluidBC     ( const fluid_bchandler_type& bc_fluid );
    //!Setter for the BCHandler of the linearized fluid problem (to be used in segregated Newton FSI)
    void setLinFluidBC          ( const fluid_bchandler_type& bc_dfluid )     { M_BCh_du     = bc_dfluid; }
    //!Setter for the BCHandler of the inverse linearized fluid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
     */
    void setInvLinFluidBC       ( const fluid_bchandler_type& bc_dfluid_inv ) { M_BCh_du_inv = bc_dfluid_inv; }
    //!Setter for the BCHandler of the gerometry problem (to be used in segregated Newton FSI)
    void setHarmonicExtensionBC ( const fluid_bchandler_type& bc_he );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setSolidBC     ( const solid_bchandler_type& bc_solid );
    //!Setter for the BCHandler of the linearized solid problem (to be used in segregated Newton FSI)
    void setLinSolidBC          ( const solid_bchandler_type& bc_dsolid )     { M_BCh_dz     = bc_dsolid; }
    //!Setter for the BCHandler of the inverse linearized solid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
    */
    void setInvLinSolidBC       ( const solid_bchandler_type& bc_dsolid_inv ) { M_BCh_dz_inv = bc_dsolid_inv; }

    //! Setter for the interface displacement (partitioned according to the fluid)
    void setLambdaFluid         ( const vector_type& lambda );
    //! Setter for the interface displacement (partitioned according to the solid)
    void setLambdaSolid         ( const vector_type& lambda );

    //! Setter for the solid interface displacement at the previous time step
    /**\todo{see if we can remove these}*/
    void setLambdaSolidOld      ( const vector_type& lambda );
    //! Setter for the solid interface velocity at the previous time step
    /**\todo{see if we can remove these}*/
    void setLambdaDotSolid      ( const vector_type& lambda );

    //!Setter for the fluid interface stress
    void setSigmaFluid          ( const vector_type& sigma );
    //!Setter for the solid interface stress
    void setSigmaSolid          ( const vector_type& sigma );
    //\todo{try to remove this}
    void setMinusSigmaFluid     ( const vector_type& sigma );

    //! Setter for the Robin coefficient of the Robin--Neumann coupling scheme (as a BCFunction)
    void setAlphafbcf  ( const bc_function_type& alphafbcf );

//     void setMpi     (bool mpi  ){M_mpi      = mpi;}
//     void setFluidMpi(bool fluid){M_isFluidMpi = fluid;}
//     void setSolidMpi(bool solid){M_issolidMpi = solid;}

//     bool mpi(){return M_mpi;}

//     void setMixteOuterWall                   ( function_type const& dload,
//                                                function_type const& E ) { M_bcfMixteOuterWall.setFunctions_Mixte(dload, E); }
    void setStructureDispToHarmonicExtension ( const vector_type& disp,  UInt type = 0 );
    void setStructureToFluid                 ( const vector_type& vel,   UInt type = 0 );
    void setStructureDispToFluid             ( const vector_type& vel,   UInt type = 0 );
    void setStructureDispToSolid             ( const vector_type& disp,  UInt type = 0 );
    void setDerStructureDispToSolid          ( const vector_type& ddisp, UInt type = 0 );
    void setSolidLoadToStructure             ( const vector_type& load,  UInt type = 0 );
    void setHarmonicExtensionVelToFluid      ( const vector_type& vel,   UInt type = 0 );
    void setDerHarmonicExtensionVelToFluid   ( const vector_type& dvel,  UInt type = 0 );
    //void setFluidInterfaceDisp               ( const vector_type& disp,  UInt type = 0 );
    void setFluidLoadToStructure             ( const vector_type& load,  UInt type = 0 );
    void setDerFluidLoadToStructure          ( const vector_type& dload, UInt type = 0 );
    void setDerFluidLoadToFluid              ( const vector_type& dload, UInt type = 0 );
    void setMixteOuterWall                   ( const function_type& dload, const function_type& E);

    //! sets the solution vector by reference
    virtual void setSolutionPtr              ( const vector_ptrtype& sol )      { M_lambda = sol; }

    //! sets the solution vector by copy
    virtual void setSolution                 ( const vector_type& solution )    { M_lambda.reset( new vector_type( solution ) ); }
    //! Initialize the fluid solution given a vector (calls setSolution)
    virtual void initialize                  ( vector_ptrtype u0)               { setSolution(*u0); }


    //! sets the solution time derivative vector by copy
    //void setSolutionDerivative( vector_ptrtype& lambdaDot ) { M_lambdaDot = lambdaDot; }

    //! Setter for the time derivative of the interface displacement
    void setSolutionDerivative( const vector_type& solutionDerivative )         { M_lambdaDot.reset( new vector_type( solutionDerivative ) ); }

    //!\todo{see if can be marked ad deprected}
    virtual void  setRestarts( bool restarts ) { /*M_restarts = restarts;*/ }
    //@}

//     void setDerReducedFluidLoadToStructure   ( vector_type &dload, UInt type = 0 );
//     void setDerStructureAccToReducedFluid    ( vector_type &acc,   UInt type = 0 );

//     void setReducedLinFluidBC                (fluid_bchandler_type bc_dredfluid);
//     void setInvReducedLinFluidBC             (fluid_bchandler_type bc_invdredfluid);

//     void setBCh_fluid             ( const fluid_bchandler_type& BCh_fluid)       { M_BCh_u      = BCh_fluid; }
//     void setBCh_HarmonicExtension ( const fluid_bchandler_type& BCh_mesh)        { M_BCh_mesh   = BCh_mesh; }
//     void setBCh_fluidDer          ( const fluid_bchandler_type& BCh_fluidDer)    { M_BCh_du     = BCh_fluidDer; }
//     void setBCh_fluidDerInv       ( const fluid_bchandler_type& BCh_fluidDerInv) { M_BCh_du_inv = BCh_fluidDerInv; }
//     void setBCh_solid             ( const solid_bchandler_type& BCh_solid)       { M_BCh_d      = BCh_solid; }
//     void setBCh_solidDer          ( const solid_bchandler_type& BCh_solidDer)    { M_BCh_dz     = BCh_solidDer; }
//     void setBCh_solidDerInv       ( const solid_bchandler_type& BCh_solidDerInv) { M_BCh_dz_inv = BCh_solidDerInv; }

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
    void transferMeshMotionOnFluid( const vector_type& _vec1,
                                    vector_type& _vec2 );

    //! Interpolates mesh motion into velocity
    /**
       Interpolates a vector with the map of the harmonic extension into one with the map of the fluid velocity
     */
    void interpolateVelocity(const vector_type& _vec1,
                             vector_type& _vec2);


    //!Interpolates to vectors on the interface
    /**
       The two vectors can have different numeration, for different discretizations this method is not tested.
    */
    void interpolateInterfaceDofs(const FESpace<mesh_type, EpetraMap>& _fespace1,
                                  const vector_type&                   _vec1,
                                  const FESpace<mesh_type, EpetraMap>& _fespace2,
                                  vector_type&                         _vec2,
                                  dof_interface_type3D&                _dofInterface);
    //@}

    //!@name Protected Attributes
    //@{
    boost::shared_ptr<mesh_type>                      M_mesh;

    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_uFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_pFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_dFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_mmFESpace;

    boost::shared_ptr<mesh_type>                      M_fluidMesh;
    boost::shared_ptr<mesh_type>                      M_solidMesh;

    boost::shared_ptr<partitionMesh< mesh_type > >    M_fluidMeshPart;
    boost::shared_ptr<partitionMesh< mesh_type > >    M_solidMeshPart;

    fluid_bchandler_type                              M_BCh_u;
    solid_bchandler_type                              M_BCh_d;
    fluid_bchandler_type                              M_BCh_mesh;

    // interface operators BCs
    fluid_bchandler_type                              M_BCh_du;
    fluid_bchandler_type                              M_BCh_du_inv;

    solid_bchandler_type                              M_BCh_dz;
    solid_bchandler_type                              M_BCh_dz_inv;

    fluid_bchandler_type                              M_BCh_dp;
    fluid_bchandler_type                              M_BCh_dp_inv;

    fluid_type                                        M_fluid;
    solid_type                                        M_solid;
    meshmotion_type                                   M_meshMotion;

//     fluidlin_type                                     M_fluidLin;
//     solidlin_type                                     M_solidLin;

    boost::shared_ptr<BdfT<vector_type> >             M_bdf;

    GetPot                                            M_dataFile;

    boost::shared_ptr<DataMesh>                       M_dataMeshFluid;
    boost::shared_ptr<DataMesh>                       M_dataMeshSolid;

    data_PtrType                                      M_data;

//     quasi_newton_type         M_reducedLinFluid;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMap;

    //!\todo{kill this attribute}
    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMapOnZero;
    //!\todo{kill this attribute}
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMapOnZero;

    dof_interface_type3D                              M_dofFluidToStructure; // Needed
//     dof_interface_type3D                              M_dofSolidToFluid;
    dof_interface_type3D                              M_dofStructureToFluid; // Needed
    dof_interface_type3D                              M_dofStructureToSolid; // Needed to construct M_bcvStructureDispToSolid
    dof_interface_type3D                              M_dofStructureToHarmonicExtension; // Needed to construct interface maps
    dof_interface_type3D                              M_dofHarmonicExtensionToFluid; // Needed to construct M_bcvStructureToFluid
//     dof_interface_type3D                              M_dofStructureToReducedFluid;
//     dof_interface_type3D                              M_dofReducedFluidToStructure;

    dof_interface_type2D                              M_dofFluid;
    dof_interface_type2D                              M_dofSolid;
    dof_interface_type2D                              M_dofFluidInv;
    dof_interface_type2D                              M_dofSolidInv;

    bc_vector_interface                               M_bcvFluidInterfaceDisp;
    bc_vector_interface                               M_bcvFluidLoadToStructure;
    bc_vector_interface                               M_bcvSolidLoadToStructure;
    bc_vector_interface                               M_bcvStructureToFluid;
    bc_vector_interface                               M_bcvStructureDispToFluid;
    bc_vector_interface                               M_bcvStructureDispToSolid;
    bc_vector_interface                               M_bcvStructureDispToHarmonicExtension;
    bc_vector_interface                               M_bcvHarmonicExtensionVelToFluid;
//     bc_vector_interface                               M_bcvStructureToReducedFluid;
//     bc_vector_interface                               M_bcvReducedFluidToStructure;

    bc_vector_interface                               M_bcvDerHarmonicExtensionVelToFluid;
    bc_vector_interface                               M_bcvDerFluidLoadToStructure;
    bc_vector_interface                               M_bcvDerFluidLoadToFluid;
    bc_vector_interface                               M_bcvDerStructureDispToSolid;
    BCFunctionMixte                                   M_bcfMixteOuterWall;

//     bc_vector_interface                               M_bcvDerReducedFluidLoadToStructure;
//     bc_vector_interface                               M_bcvDerStructureAccToReducedFluid;

    boost::shared_ptr<vector_type>                    M_lambdaFluid;
    boost::shared_ptr<vector_type>                    M_lambdaFluidRepeated;
    boost::shared_ptr<vector_type>                    M_lambda;
    boost::shared_ptr<vector_type>                    M_lambdaDot;


    boost::shared_ptr<vector_type>                    M_un;
    boost::shared_ptr<vector_type>                    M_rhs;
    boost::shared_ptr<vector_type>                    M_Alphaf;

    Real                                              M_AlphafCoef;
    //\todo{try to set as deprecated}
    Real                                              M_betamedio;
//     boost::shared_ptr<vector_type>                    M_w;
//     boost::shared_ptr<vector_type>                    M_dw;

    UInt                                              M_fluxes;

    comm_PtrType                                      M_epetraComm;
    comm_PtrType                                      M_epetraWorldComm;

    //@}
private:
    //! @name Private Attributes
    //@{
    // displacement on the interface
    boost::shared_ptr<vector_type>                    M_lambdaSolid;
    boost::shared_ptr<vector_type>                    M_lambdaSolidRepeated;

    boost::shared_ptr<vector_type>                    M_lambdaSolidOld;
    boost::shared_ptr<vector_type>                    M_lambdaDotSolid;
    boost::shared_ptr<vector_type>                    M_lambdaDotSolidRepeated;

    boost::shared_ptr<vector_type>                    M_sigmaFluid;
    boost::shared_ptr<vector_type>                    M_sigmaSolid;

    boost::shared_ptr<vector_type>                    M_sigmaFluidRepeated;
    boost::shared_ptr<vector_type>                    M_sigmaSolidRepeated;

    //\todo{try to set as deprecated}
    boost::shared_ptr<vector_type>                    M_minusSigmaFluid;
    //\todo{try to set as deprecated}
    boost::shared_ptr<vector_type>                    M_minusSigmaFluidRepeated;

//     boost::shared_ptr<vector_type>                    M_dispStruct;
//     boost::shared_ptr<vector_type>                    M_dispStructOld;
//     boost::shared_ptr<vector_type>                    M_veloStruct;

    boost::shared_ptr<vector_type>                    M_dispFluidMeshOld;
    boost::shared_ptr<vector_type>                    M_veloFluidMesh;
    boost::shared_ptr<vector_type>                    M_derVeloFluidMesh;

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
