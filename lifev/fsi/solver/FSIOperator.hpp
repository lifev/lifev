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

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer  Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010

    This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
    Jacobian matrix, which make it suited for the generalized Newton algorithm implemented in NonlinearRichardson.hpp. The fluid
    and structure classes are members of this class and different formulations (e.g. Monolithic \cite CrosettoEtAl2009 , segregated
    Newton \cite FernandezMoubachir2005 , Dirichlet--Neumann \cite DeparisDiscacciati2006 , Robin Neumann \cite BadiaNobileVergara2008 ) are implemented in the derived classes.

 */


#ifndef FSIOPERATOR_H
#define FSIOPERATOR_H

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


#include <lifev/core/util/Factory.hpp>
#include <lifev/core/util/FactorySingleton.hpp>

#include <lifev/core/mesh/MeshData.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

#include <lifev/core/algorithm/NonLinearAitken.hpp>

#include <lifev/structure/solver/StructuralOperator.hpp>

#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinearPenalized.hpp>
#include <lifev/structure/solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp>

#ifdef ENABLE_ANISOTROPIC_LAW
#include <lifev/structure/solver/anisotropic/HolzapfelMaterialNonLinear.hpp>
#include <lifev/structure/solver/anisotropic/HolzapfelGeneralizedMaterialNonLinear.hpp>
#endif

#include <lifev/core/fem/DOFInterface3Dto3D.hpp>
#include <lifev/core/fem/DOFInterface3Dto2D.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCFunction.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/eta/fem/ETFESpace.hpp>

#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5Mesh3D.hpp>
#endif

#include <lifev/fsi/solver/FSIData.hpp>
#include <lifev/navier_stokes/solver/OseenSolverShapeDerivative.hpp>
#include <lifev/fsi/solver/HarmonicExtensionSolver.hpp>

namespace LifeV
{
/*!
  @class FSIOperator
  @brief Fluid-Structure Interface operator class

  This is the base class for the FSI solvers in LifeV. It contains the methods to evaluate the residual and compute the
  Jacobian matrix, which make it suited for the generalized Newton algorithm implemented in NonlinearRichardson.hpp. The fluid
  and structure classes are members of this class and different formulations (e.g. Monolithic \cite CrosettoEtAl2009 , segregated
  Newton \cite FernandezMoubachir2005 , Dirichlet--Neumann \cite DeparisDiscacciati2006 , Robin Neumann \cite BadiaNobileVergara2008 ) are implemented in the derived classes.

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

    typedef RegionMesh<LinearTetra>                                               mesh_Type;
#ifdef HAVE_HDF5
    typedef ExporterHDF5Mesh3D<mesh_Type>                                           meshFilter_Type;
#endif

    typedef OseenSolverShapeDerivative   <mesh_Type>                                fluid_Type;
    typedef StructuralOperator  <mesh_Type>                                           solid_Type;
    typedef HarmonicExtensionSolver<mesh_Type>                                      meshMotion_Type;
    typedef OseenSolverShapeDerivative   <mesh_Type>                                fluidLin_Type;
    typedef StructuralOperator  <mesh_Type>                                           solidLin_Type;
    typedef std::shared_ptr<fluid_Type>                                           fluidPtr_Type;
    typedef std::shared_ptr<solid_Type>                                           solidPtr_Type;
    typedef std::shared_ptr<meshMotion_Type>                                      meshMotionPtr_Type;
    typedef std::shared_ptr<fluidLin_Type>                                        fluidLinPtr_Type;
    typedef std::shared_ptr<solidLin_Type>                                        solidLinPtr_Type;
    typedef fluid_Type::vector_Type/*fluidPtr_Type::vector_Type*/                   vector_Type;
    typedef std::shared_ptr<vector_Type>                                          vectorPtr_Type;
    typedef vector_Type                                                             solution_Type;
    typedef std::shared_ptr<solution_Type>                                        solutionPtr_Type;
    typedef fluid_Type::source_Type/*fluidPtr_Type::source_Type*/                   fluidSource_Type;
    typedef solid_Type::source_Type                                                 solidSource_Type;
    typedef std::function < Real ( const Real&, const Real&,
                                     const Real&, const Real&, const ID& ) >           function_Type;
    typedef Real                                                                    ( *bcFunction_Type ) ( const Real&, const Real&,
            const Real&, const Real&, const ID& );
    typedef std::shared_ptr<DOFInterface3Dto3D>                                   dofInterface3DPtr_Type;
    typedef std::shared_ptr<DOFInterface3Dto2D>                                   dofInterface2DPtr_Type;
    typedef std::shared_ptr<BCVectorInterface>                                    bcVectorInterfacePtr_Type;
    typedef fluid_Type::bcHandlerPtr_Type/*fluidPtr_Type::bchandlerPtr_Type*/       fluidBchandlerPtr_Type;
    typedef fluid_Type::bcHandler_Type/*fluidPtr_Type::bchandler_Type*/             fluidBchandler_Type;
    typedef BCHandler                                                               solidBchandler_Type;
    typedef std::shared_ptr<solidBchandler_Type>                                  solidBchandlerPtr_Type;
    typedef FSIData                                                                 data_Type;
    typedef std::shared_ptr<data_Type>                                            dataPtr_Type;
    typedef std::map<ID, ID>::const_iterator                                        iterator_Type;
    typedef FactorySingleton<Factory<FSIOperator, std::string> >                    FSIFactory_Type;
    typedef Displayer::commPtr_Type/*Displayer::commPtr_Type*/                      commPtr_Type;
    typedef GetPot                                                                  dataFile_Type;
    typedef std::shared_ptr<dataFile_Type>                                        dataFilePtr_Type;
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
    virtual void setDataFile ( const dataFile_Type& data );

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
    void partitionMeshes ( meshFilter_Type& fluidMeshFilter, meshFilter_Type& solidMeshFilter  );
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
    virtual void setupDOF ( meshFilter_Type& /*filterMesh*/ ) {}
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
    virtual void setupFluidSolid ( UInt const fluxes );

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
    void couplingVariableExtrap( );

    //! solves the Jacobian system
    /**
       The implementation of this method distinguish the various FSI formulations which derive from this class.
       For this reason it must be pure virtual, snd implemented in the child classes.
       \param muk: unknown solution at the k-th nonlinear iteration
       \param res: residual vector (the right hand side of the Jacobian system)
       \param linearRelTol: tolerance for the nonlinear solver
     */
    virtual void solveJac ( vector_Type&       muk,
                            const vector_Type& res,
                            const Real         linearRelTol ) = 0;

    //! Evaluates the nonlinear residual of the FSI system
    /**
       The implementation of this method also depends on the child classes, though it does not characterize them.
       \param res:  residual vector  to be computed
       \param disp: current unknown solution
       \param iter: nonlinear iteration counter. The part of th rhs related to the time discretization is computed only for iter=0
     */
    virtual void evalResidual ( vector_Type& res, const vector_Type&  disp, const UInt iter ) = 0;

    //! Update the solution after NonLinearRichardson is called.
    /*!
     *  Eventually it can update also some post-processing quantity.
     */
    virtual void updateSolution ( const vector_Type& solution );

    //! Set vectors for restart
    /*!
     *  Set vectors for restart
     */
    virtual void setVectorInStencils ( const vectorPtr_Type& /*vel*/,
                                       const vectorPtr_Type& /*pressure*/,
                                       const vectorPtr_Type& /*solidDisp*/,
                                       //                      const vectorPtr_Type& /*fluidDisp*/,
                                       const UInt /*iter*/) {}

    virtual void setFluidVectorInStencil ( const vectorPtr_Type& /*vel*/, const vectorPtr_Type& /*pressure*/, const UInt /*iter*/) {}

    virtual void setSolidVectorInStencil ( const vectorPtr_Type& /*solidDisp*/, const UInt /*iter*/) {}

    virtual void setALEVectorInStencil ( const vectorPtr_Type& /*fluidDisp*/, const UInt /*iter*/, const bool /*lastVector*/ ) {}

    virtual void finalizeRestart( ) {}

    //! Initializes all the quantities using functions
    /*!
     * calls the initialize methods for the subproblems. The mesh velocity is used to compute the convective term in the fluid equations
     * \param u0: initial fluid velocity
     * \param p0: initial fluid pressure
     * \param d0: initial solid displacement
     * \param w0: initial mesh velocity
     */
    virtual void initialize ( fluid_Type::function_Type const& u0,
                              fluid_Type::function_Type const& p0,
                              solid_Type::function const& d0,
                              solid_Type::function const& w0,
                              fluid_Type::function_Type const& df0 );

    //@}



    //!@name MONOLITHIC Solver Methods - Implemented There
    //{@
    virtual void iterateMesh ( const vector_Type& /*disp*/ )
    {
        assert (false);
    }
    virtual void setupBDF ( const vector_Type& /*u0*/ ) { }
    virtual void updateRHS() {}
    virtual void applyBoundaryConditions() {}
    //@}


    //!@name Factory Methods
    //@{
    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createVenantKirchhoffLinear()
    {
        return new VenantKirchhoffMaterialLinear< FSIOperator::mesh_Type >();
    }

    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createVenantKirchhoffNonLinear()
    {
        return new VenantKirchhoffMaterialNonLinear< FSIOperator::mesh_Type >();
    }
    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createExponentialMaterialNonLinear()
    {
        return new ExponentialMaterialNonLinear< FSIOperator::mesh_Type >();
    }

    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createNeoHookeanMaterialNonLinear()
    {
        return new NeoHookeanMaterialNonLinear< FSIOperator::mesh_Type >();
    }
    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createVenantKirchhoffNonLinearPenalized()
    {
        return new VenantKirchhoffMaterialNonLinearPenalized< FSIOperator::mesh_Type >();
    }

    static StructuralIsotropicConstitutiveLaw< FSIOperator::mesh_Type >*    createSecondOrderExponentialMaterialNonLinear()
    {
        return new SecondOrderExponentialMaterialNonLinear< FSIOperator::mesh_Type >();
    }

#ifdef ENABLE_ANISOTROPIC_LAW
    static StructuralAnisotropicConstitutiveLaw< FSIOperator::mesh_Type >*  createHolzapfelMaterialNonLinear()
    {
        return new HolzapfelMaterialNonLinear<MeshType >();
    }
    static StructuralAnisotropicConstitutiveLaw< FSIOperator::mesh_Type >*  createHolzapfelGeneralizedMaterialNonLinear()
    {
        return new HolzapfelGeneralizedMaterialNonLinear<MeshType >();
    }


#endif


    //@}


    //!@name Public Methods
    //@{

    //!@name Public Methods
    //@{
    //!Initializes the TimeAdvance scheme which should handle the fluid time discretization, solid and move mesh
    /**
       Initialization of the time advancing classes for fluid, structure and geometry problem.
     */
    void initializeTimeAdvance ( const std::vector<vectorPtr_Type>& initialFluidVel, const std::vector<vectorPtr_Type>& initialSolidDisp, const std::vector<vectorPtr_Type>&  initialFluiDisp);

    virtual void initializeMonolithicOperator ( std::vector< vectorPtr_Type> /*u0*/, std::vector< vectorPtr_Type> /*ds0*/, std::vector< vectorPtr_Type> /*df0*/) {}

    //! initializes the fluid solver with vectors
    /**
       \param velAndPressure: initial vector containing the velocity and pressure
       \param displacement: initial vector containing the mesh displacement
     */
    void initializeFluid ( const vector_Type& velAndPressure, const vector_Type& displacement );

    //! initializes the solid solver with vectors
    /**
       \param displacement: initial vector containing the structure displacement
       \param velocity: initial vector containing the velocity, used for the initialization of the TimeAdvanceNewmark scheme
     */
    void initializeSolid ( vectorPtr_Type displacement, vectorPtr_Type /*velocity*/ );

    //!moves the mesh using the solution of the harmonic extension equation
    /**
       \param disp displacement of the mesh, must be the difference between the current solution of the HE problem and the one at the previous time step.
     */
    void moveMesh ( const vector_Type& disp );

    //     void solveLinearFluid();
    //     void solveLinearSolid();

    //! Creates the Epetra maps for the interface
    /**
       Given a std::map holding the numeration of the interface according to the fluid (first) and the solid (second),
       builds the MapEpetras for the interface dofs (e.g. M_fluidInterfaceMap). Note that when both the fluid and solid
       meshes are partitioned (e.g. in the monolithic solver) the local dofs are those of the FLUID partition of the
       interface, even when the numeration refers to the solid.
     */
    void createInterfaceMaps (std::map<ID, ID> const& locDofMap);

    //!Method to import an VectorEpetra defined on the fluid map (i.e. with the fluid numeration of the dofs) to the interface
    /**
       Note that the output vector will have the solid numeration on the interface! By default in fact the vectors on the
       FSI interface in LifeV are numerated according to the solid.
     */
    void transferFluidOnInterface ( const vector_Type& _vec1, vector_Type& _vec2 );

    void transferSolidOnFluid ( const vector_Type& _vec1, vector_Type& _vec2 );

    //!Method to import an VectorEpetra defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       The output vector has the solid numeration of the dofs and is partitioned according to the solid partition. This method is not used in the monolithic solvers.
     */
    void transferSolidOnInterface ( const vector_Type& _vec1, vector_Type& _vec2 );
    //     void transferInterfaceOnFluid( const vector_Type& _vec1, vector_Type& _vec2 );

    //!Method to import an VectorEpetra defined on the solid map (i.e. with the solid numeration of the dofs) to the interface
    /**
       the output vector have the numeration of the solid, as in transferSolidOnInterface, but is partitioned according to the fluid! This method is not used in the monolithic solvers.
     */
    void transferInterfaceOnSolid ( const vector_Type& _vec1, vector_Type& _vec2 );

    //! Update the RHS on the base of the fluid BC
    /*!
     *  This method is used by the Multiscale to update the RHS vector for the nonlinear subiterations.
     *  @param bcHandlerFluid fluid BC handler
     *  @param rhs RHS of the FSI problem
     */
    void bcManageVectorRHS ( const fluidBchandlerPtr_Type& bch, vector_Type& rhs );

    //! Update the RHS on the base of the fluid and solid BC
    /*!
     *  This method is used by the Multiscale to update the RHS vector for the nonlinear subiterations.
     *  @param bcHandlerFluid fluid BC handler
     *  @param bcHandlerSolid solid BC handler
     *  @param rhs RHS of the FSI problem
     */
    void bcManageVectorRHS ( const fluidBchandlerPtr_Type& bcHandlerFluid, const solidBchandlerPtr_Type& bcHandlerSolid, vector_Type& rhs );

    //! Method to set the Robin vector coefficient of the Robin--Neumann coupling scheme (as a constant vector vector)
    void setAlphaf()
    {
        M_alphaF->epetraVector().PutScalar ( M_alphaFCoef );
    }
    //! Method to compute the scalar coefficient \f$\alpha\f$ of the Robin--Neumann coupling scheme
    void setAlphafCoef();
    //! Method calling setAlphaf and setAlphafCoef
    void setStructureToFluidParameters();

    //! Reset the right hand side to zero
    /*!
     *  This method is used in the multiscale framework during subiterations
     */
    void resetRHS()
    {
        *M_rhs = 0;
    }

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

    //! Get the extrapolation of the solution
    /*!
     *  @param extrapolation vector to be filled with the extrapolated solution
     */
    void extrapolation ( vector_Type& extrapolation ) const
    {
        M_fluidTimeAdvance->extrapolation ( extrapolation );
    }

    //     vector_Type & displacement()                                        { return *M_lambdaSolid; }
    //     vector_Type & displacementOld()                                     { return *M_lambdaSolidOld; }
    //     vector_Type & residual();
    //     vector_Type const & velocity()                                const { return *M_lambdaDotSolid; }
    //     vector_Type & residualFSI();

    //! Returns the number of imposed fluxes
    UInt imposedFluxes();

    const vector_Type& lambdaFluid()                              const
    {
        return *M_lambdaFluid;
    }
    const vector_Type& lambdaSolid()                              const
    {
        return *M_lambdaSolid;
    }
    const vector_Type& lambdaSolidOld()                           const
    {
        return *M_lambdaSolidOld;
    }
    const vector_Type& lambdaDotSolid()                           const
    {
        return *M_lambdaDotSolid;
    }
    const vector_Type& sigmaFluid()                               const
    {
        return *M_sigmaFluid;
    }
    const vector_Type& sigmaSolid()                               const
    {
        return *M_sigmaSolid;
    }

    const vector_Type& lambdaFluidRepeated()                      const
    {
        return *M_lambdaFluidRepeated;
    }
    const vector_Type& lambdaSolidRepeated()                      const
    {
        return *M_lambdaSolidRepeated;
    }
    const vector_Type& lambdaDotSolidRepeated()                   const
    {
        return *M_lambdaDotSolidRepeated;
    }
    const vector_Type& sigmaFluidRepeated()                       const
    {
        return *M_sigmaFluidRepeated;
    }
    const vector_Type& sigmaSolidRepeated()                       const
    {
        return *M_sigmaSolidRepeated;
    }

    //!\todo{remove this method}
    // now residual is Ax-b, but we should decide for b-Ax. In the meantime, we need b-Ax:
    const vector_Type& minusSigmaFluid()                          const
    {
        return *M_minusSigmaFluid;
    }
    //!\todo{remove this method}
    const vector_Type& minusSigmaFluidRepeated()                  const
    {
        return *M_minusSigmaFluidRepeated;
    }

    //!coefficient for the Robin--Neumann coupling scheme
    vector_Type&       Alphaf()                                   const
    {
        return *M_alphaF;
    }

    commPtr_Type worldComm()                                      const
    {
        return M_epetraWorldComm;
    }

    bool isFluid()                                                const
    {
        return M_isFluid;
    }
    bool isSolid()                                                const
    {
        return M_isSolid;
    }

    bool isLinearFluid()                                          const
    {
        return M_linearFluid;
    }
    bool isLinearSolid()                                          const
    {
        return M_linearSolid;
    }

    int getFluidLeaderId()                                        const
    {
        return M_fluidLeader;
    }
    int getSolidLeaderId()                                        const
    {
        return M_solidLeader;
    }

    //! Getter for the fluid solver
    const fluid_Type& fluid()                                     const
    {
        return *M_fluid;
    }
    //! Getter for the solid solver
    const solid_Type& solid()                                     const
    {
        return *M_solid;
    }
    //! Getter for the harmonic extension solver
    const meshMotion_Type& meshMotion()                           const
    {
        return *M_meshMotion;
    }

    //! Getter-Setter for the fluid solver
    /** \todo{mark as deprecated}*/
    fluid_Type& fluid()
    {
        return *M_fluid;
    }
    //! Getter-Setter for the solid solver
    /** \todo{mark as deprecated}*/
    solid_Type& solid()
    {
        return *M_solid;
    }
    //! Getter-Setter for the mesh motion solver
    /** \todo{mark as deprecated}*/
    meshMotion_Type& meshMotion()
    {
        return *M_meshMotion;
    }
    //     fluidLin_Type & fluidLin()                               { return *M_fluidLin; }
    //     solidLin_Type & solidLin()                               { return *M_solidLin; }

    //!getter for the FSI data container
    const data_Type& data()                                       const
    {
        return *M_data;
    }
    //!getter for the fluid data container
    const data_Type::dataFluidPtr_Type& dataFluid()               const
    {
        return M_data->dataFluid();
    }
    //!getter for the solid data container
    const data_Type::dataSolidPtr_Type& dataSolid()               const
    {
        return M_data->dataSolid();
    }

    //!getter for the unpartitioned fluid mesh
    mesh_Type& fluidMesh()                                        const
    {
        return *M_fluidMesh;
    }
    //!getter for the unpartitioned solid mesh
    mesh_Type& solidMesh()                                        const
    {
        return *M_solidMesh;
    }

    // const mesh_Type& fluidMesh()                                  const { return *M_fluidMesh; }
    // const mesh_Type& solidMesh()                                  const { return *M_solidMesh; }

    //!getter for the partitioned fluid mesh
    mesh_Type& fluidLocalMesh()
    {
        return *M_fluidLocalMesh;
    }
    //!getter for the partitioned solid mesh
    mesh_Type& solidLocalMesh()
    {
        return *M_solidLocalMesh;
    }

    //!getter for the fluid velocity FESpace
    const FESpace<mesh_Type, MapEpetra>& uFESpace()               const
    {
        return *M_uFESpace;
    }
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > uFESpacePtr() const
    {
        return M_uFESpace;
    }
    //!getter for the fluid pressure FESpace
    const FESpace<mesh_Type, MapEpetra>& pFESpace()               const
    {
        return *M_pFESpace;
    }
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > pFESpacePtr() const
    {
        return M_pFESpace;
    }
    //!getter for the solid displacement FESpace
    const FESpace<mesh_Type, MapEpetra>& dFESpace() const
    {
        return *M_dFESpace;
    }
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > dFESpacePtr() const
    {
        return M_dFESpace;
    }
    //!getter for the solid displacement FESpace
    const ETFESpace<mesh_Type, MapEpetra, 3, 3>& dFESpaceET() const
    {
        return *M_dETFESpace;
    }
    std::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 3> > dFESpaceETPtr() const
    {
        return M_dETFESpace;
    }
    //!getter for the harmonic extension solution FESpace
    const FESpace<mesh_Type, MapEpetra>& mmFESpace()              const
    {
        return *M_mmFESpace;
    }
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > mmFESpacePtr() const
    {
        return M_mmFESpace;
    }
    //!getter for the harmonic extension solution
    const vector_Type& meshDisp()                         const
    {
        return M_ALETimeAdvance->singleElement (0);
    }
    //!getter for the harmonic extension solution of the previous time step
    const         vector_Type& dispFluidMeshOld()                 const
    {
        return *M_dispFluidMeshOld;
    }
    //!getter for the mesh velocity
    virtual       vector_Type& veloFluidMesh()
    {
        return *M_veloFluidMesh;
    }
    //!getter for the mesh velocity increment (used for Newton FSI)
    //! \todo{try to remove this method}
    vector_Type& derVeloFluidMesh()
    {
        return *M_derVeloFluidMesh;
    }

    const dofInterface3DPtr_Type& dofFluidToStructure()             const
    {
        return M_dofFluidToStructure;
    }
    const dofInterface3DPtr_Type& dofStructureToSolid()             const
    {
        return M_dofStructureToSolid;
    }
    const dofInterface3DPtr_Type& dofStructureToHarmonicExtension() const
    {
        return M_dofStructureToHarmonicExtension;
    }
    const dofInterface3DPtr_Type& dofHarmonicExtensionToFluid()     const
    {
        return M_dofHarmonicExtensionToFluid;
    }

    std::shared_ptr<MapEpetra>& fluidInterfaceMap()
    {
        return M_fluidInterfaceMap;
    }
    std::shared_ptr<MapEpetra>& solidInterfaceMap()
    {
        return M_solidInterfaceMap;
    }

    //! Getter for the map of the variable used for the coupling
    virtual std::shared_ptr<MapEpetra>& couplingVariableMap()
    {
        return M_solidInterfaceMap;
    }

    //! Method to implement Robin boundary conditions on the external wall for the structure
    BCFunctionRobin& bcfRobinOuterWall()
    {
        return M_bcfRobinOuterWall;
    }

    bcVectorInterfacePtr_Type bcvStructureDisptoFluid()                 const
    {
        return M_bcvStructureDispToFluid;
    }
    bcVectorInterfacePtr_Type bcvStructureToFluid()                     const
    {
        return M_bcvStructureToFluid;
    }
    bcVectorInterfacePtr_Type bcvSolidLoadToStructure()                 const
    {
        return M_bcvSolidLoadToStructure;
    }
    bcVectorInterfacePtr_Type bcvFluidInterfaceDisp()                   const
    {
        return M_bcvFluidInterfaceDisp;
    }
    bcVectorInterfacePtr_Type bcvHarmonicExtensionVelToFluid()          const
    {
        return M_bcvHarmonicExtensionVelToFluid;
    }
    bcVectorInterfacePtr_Type bcvDerHarmonicExtensionVelToFluid()       const
    {
        return M_bcvDerHarmonicExtensionVelToFluid;
    }
    bcVectorInterfacePtr_Type bcvStructureDispToHarmonicExtension()     const
    {
        return M_bcvStructureDispToHarmonicExtension;
    }
    bcVectorInterfacePtr_Type bcvStructureDispToSolid()                 const
    {
        return M_bcvStructureDispToSolid;
    }
    bcVectorInterfacePtr_Type bcvDerStructureDispToSolid()              const
    {
        return M_bcvDerStructureDispToSolid;
    }
    bcVectorInterfacePtr_Type bcvFluidLoadToStructure()                 const
    {
        return M_bcvFluidLoadToStructure;
    }
    bcVectorInterfacePtr_Type bcvDerFluidLoadToStructure()              const
    {
        return M_bcvDerFluidLoadToStructure;
    }
    bcVectorInterfacePtr_Type bcvDerFluidLoadToFluid()                  const
    {
        return M_bcvDerFluidLoadToFluid;
    }
    //     bcVectorInterfacePtr_Type bcvDerReducedFluidLoadToStructure()             { return M_bcvDerReducedFluidLoadToStructure; }
    //     bcVectorInterfacePtr_Type bcvDerStructureAccToReducedFluid()              { return M_bcvDerStructureAccToReducedFluid; }

    //! Getter for the BCHandler of the fluid problem
    const fluidBchandlerPtr_Type& BCh_fluid()                       const
    {
        return M_BCh_u;
    }
    //! Getter for the BCHandler of the harmonic extension problem
    const fluidBchandlerPtr_Type& BCh_harmonicExtension()           const
    {
        return M_BCh_mesh;
    }
    //! Getter for the BCHandler of the linearized fluid problem (to be used in Newton for the partitioned FSI)
    const fluidBchandlerPtr_Type& BCh_du()                          const
    {
        return M_BCh_du;
    }
    //! Getter for the BCHandler of the linearized inverse of the fluid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const fluidBchandlerPtr_Type& BCh_du_inv()                      const
    {
        return M_BCh_du_inv;
    }
    //! Getter for the BCHandler of the solid problem
    const solidBchandlerPtr_Type& BCh_solid()                       const
    {
        return M_BCh_d;
    }
    //! Getter for the BCHandler of the linearized solid problem
    const solidBchandlerPtr_Type& BCh_dz()                          const
    {
        return M_BCh_dz;
    }
    //! Getter for the BCHandler of the linearized inverse of the solid Steklov Poincare' operator (not used)
    //! \todo{mark as deprecated untill debugged}
    const solidBchandlerPtr_Type& BCh_dz_inv()                      const
    {
        return M_BCh_dz_inv;
    }

    //! Getter for the right hand side
    const vectorPtr_Type& getRHS()                                  const
    {
        return M_rhs;
    }

    const std::shared_ptr<TimeAdvance<vector_Type> > ALETimeAdvance() const
    {
        return  M_ALETimeAdvance;
    }
    const std::shared_ptr<TimeAdvance<vector_Type> > fluidTimeAdvance() const
    {
        return  M_fluidTimeAdvance;
    }
    const std::shared_ptr<TimeAdvance<vector_Type> > solidTimeAdvance() const
    {
        return  M_solidTimeAdvance;
    }

    const string ALETimeAdvanceMethod() const
    {
        return  M_ALETimeAdvanceMethod;
    }
    const string fluidTimeAdvanceMethod() const
    {
        return  M_fluidTimeAdvanceMethod;
    }
    const string solidTimeAdvanceMethod() const
    {
        return  M_solidTimeAdvanceMethod;
    }

    //! gets the solution vector by reference
    virtual const vector_Type& solution()                           const
    {
        return *M_lambda;
    }

    //! gets the solid displacement by copy
    virtual void getSolidDisp ( vector_Type& soliddisp )
    {
        soliddisp = M_solid->displacement();
    }

    //! gets the solid velocity by copy
    virtual void getSolidVel ( vector_Type& solidvel )
    {
        solidvel = M_solidTimeAdvance->firstDerivative();
    }

    //! Export the solid displacement by copying it to an external vector
    /*!
     * @param solidDisplacement vector to be filled with the solid displacement
     */
    virtual void exportSolidDisplacement ( vector_Type& solidDisplacement )
    {
        solidDisplacement = M_solid->displacement();
    }

    //! Export the solid velocity by copying it to an external vector
    /*!
     * @param solidVelocity vector to be filled with the solid velocity
     */
    virtual void exportSolidVelocity ( vector_Type& solidVelocity )
    {
        solidVelocity = M_solidTimeAdvance->firstDerivative();
    }

    //! Export the solid acceleration by copying it to an external vector
    /*!
     * @param solidAcc vector to be filled with the solid acceleration
     */
    virtual void exportSolidAcceleration ( vector_Type& solidAcc )
    {
        solidAcc = M_solidTimeAdvance->secondDerivative();
    }

    //! Export the fluid velocity by copying it to an external vector
    /*!
     * @param fluidVelocity vector to be filled with the fluid velocity
     */
    virtual void exportFluidVelocity ( vector_Type& fluidVelocity )
    {
        fluidVelocity = *M_fluid->solution();
    }

    //! Export the fluid pressure by copying it to an external vector
    /*!
     * @param fluidPressure vector to be filled with the fluid pressure
     */
    virtual void exportFluidPressure ( vector_Type& fluidPressure )
    {
        fluidPressure = *M_fluid->solution();
    }

    //! Export the fluid velocity and pressure by copying it to an external vector
    /*!
     * @param fluidVelocityAndPressure vector to be filled with the fluid velocity and pressure
     */
    virtual void exportFluidVelocityAndPressure ( vector_Type& fluidVelocityAndPressure )
    {
        fluidVelocityAndPressure = *M_fluid->solution();
    }

    //! Export the fluid displacement by copying it to an external vector
    /*!
     * @param fluidDisplacement vector to be filled with the fluid displacement
     */
    virtual void exportFluidDisplacement ( vector_Type& fluidDisplacement )
    {
        fluidDisplacement = M_ALETimeAdvance->singleElement (0);
    }

    //@}



    /** @name  Set Functions
     */
    //@{
    //! Setter for the local and  world communicators
    /**
       The communicator can be different depending on which type of subdomain we are considering
     */
    void setComm ( const commPtr_Type& comm, const commPtr_Type& worldComm );

    //! Setter for the FSI data
    void setData ( const dataPtr_Type& data )
    {
        M_data = data;
    }

    //! Setter for the fluid and geometry problems
    void setFluid ( const fluidPtr_Type& fluid, const meshMotionPtr_Type& meshmotion );
    //! Setter for the solid problem
    void setSolid ( const solidPtr_Type& solid );

    //!Setter for the "fluid" flag
    void setFluid ( const bool& isFluid )
    {
        M_isFluid = isFluid;
    }
    //!Setter for the "solid" flag
    void setSolid ( const bool& isSolid )
    {
        M_isSolid = isSolid;
    }

    //!Setter for the "linear fluid" flag
    void setLinearFluid ( const bool& linFluid )
    {
        M_linearFluid = linFluid;
    }
    //!Setter for the "linear solid" flag
    void setLinearSolid ( const bool& linSolid )
    {
        M_linearSolid = linSolid;
    }

    void setFluidLeader ( const int& fluidLeader )
    {
        M_fluidLeader = fluidLeader;
    }
    void setSolidLeader ( const int& solidLeader )
    {
        M_solidLeader = solidLeader;
    }

    //     void setBC( fluidBchandlerPtr_Type& bc_u, solidBchandlerPtr_Type& bc_d, fluidBchandlerPtr_Type& bc_m );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setFluidBC     ( const fluidBchandlerPtr_Type& bc_fluid );
    //!Setter for the BCHandler of the linearized fluid problem (to be used in segregated Newton FSI)
    void setLinFluidBC          ( const fluidBchandlerPtr_Type& bc_dfluid )
    {
        M_BCh_du     = bc_dfluid;
    }
    //!Setter for the BCHandler of the inverse linearized fluid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
     */
    void setInvLinFluidBC       ( const fluidBchandlerPtr_Type& bc_dfluid_inv )
    {
        M_BCh_du_inv = bc_dfluid_inv;
    }
    //!Setter for the BCHandler of the gerometry problem (to be used in segregated Newton FSI)
    void setHarmonicExtensionBC ( const fluidBchandlerPtr_Type& bc_he );

    //!Setter for the fluid BCHandler
    /**
       \todo{see if this needs to be virtual}
     */
    virtual void setSolidBC     ( const solidBchandlerPtr_Type& bc_solid );
    //!Setter for the BCHandler of the linearized solid problem (to be used in segregated Newton FSI)
    void setLinSolidBC          ( const solidBchandlerPtr_Type& bc_dsolid )
    {
        M_BCh_dz     = bc_dsolid;
    }
    //!Setter for the BCHandler of the inverse linearized solid steklov Poincare' operator (to be used in SP FSI formulation)
    /**
       \todo{mark as deprecated until not debugged}
    */
    void setInvLinSolidBC       ( const solidBchandlerPtr_Type& bc_dsolid_inv )
    {
        M_BCh_dz_inv = bc_dsolid_inv;
    }

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

    //! Setter for the time derivative of the interface displacement
    void setSolutionDerivative ( const vector_Type& solutionDerivative )
    {
        M_lambdaDot.reset ( new vector_Type ( solutionDerivative ) );
    }

    //! Setup of the TimeAdvance classes given the input data file
    void
    setupTimeAdvance ( const dataFile_Type& dataFile );

    //@}


protected:

    //virtual void variablesInit(const ReferenceFE* refFE_struct,const LifeV::QuadratureRule*  bdQr_struct, const LifeV::QuadratureRule* qR_struct);
    //!@name Protected Methods
    //@{
    //!initailize the variables
    /**
       instantiates the pointers which are used in the segregated FSI solvers. Reimplemented in the Monolithic class.
       \param dorder: unused parameter
     */
    virtual void variablesInit ( const std::string& dOrder );

    //!Interpolates the mesh motion dofs on the fluid
    /**
       The order of the spatial approximation depends on this method: when the mesh motion approximation is first order
       in space the overall approximation is of the first order even if the fluid is solved with hicher order FEs.
       Calls the interpolateVelocity method
     */
    void transferMeshMotionOnFluid ( const vector_Type& _vec1, vector_Type& _vec2 );

    //! Interpolates mesh motion into velocity
    /**
       Interpolates a vector with the map of the harmonic extension into one with the map of the fluid velocity
     */
    void interpolateVelocity (const vector_Type& _vec1, vector_Type& _vec2);


    //!Interpolates to vectors on the interface
    /**
       The two vectors can have different numeration, for different discretizations this method is not tested.
    */
    void interpolateInterfaceDofs (const FESpace<mesh_Type, MapEpetra>& _fespace1,
                                   const vector_Type&                   _vec1,
                                   const FESpace<mesh_Type, MapEpetra>& _fespace2,
                                   vector_Type&                         _vec2,
                                   dofInterface3DPtr_Type&              _dofInterface);

    //@}

    //!@name Protected Attributes
    //@{

    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_uFESpace;
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_pFESpace;
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_dFESpace;
    std::shared_ptr<ETFESpace<mesh_Type, MapEpetra, 3, 3> > M_dETFESpace;
    std::shared_ptr<FESpace<mesh_Type, MapEpetra> > M_mmFESpace;

    std::shared_ptr<mesh_Type>                      M_fluidMesh;
    std::shared_ptr<mesh_Type>                      M_solidMesh;

    std::shared_ptr<mesh_Type>                      M_fluidLocalMesh;
    std::shared_ptr<mesh_Type>                      M_solidLocalMesh;

    fluidBchandlerPtr_Type                            M_BCh_u;
    solidBchandlerPtr_Type                            M_BCh_d;
    fluidBchandlerPtr_Type                            M_BCh_mesh;

    // interface operators BCs
    fluidBchandlerPtr_Type                            M_BCh_du;
    fluidBchandlerPtr_Type                            M_BCh_du_inv;

    solidBchandlerPtr_Type                            M_BCh_dz;
    solidBchandlerPtr_Type                            M_BCh_dz_inv;

    fluidBchandlerPtr_Type                            M_BCh_dp;
    fluidBchandlerPtr_Type                            M_BCh_dp_inv;

    fluidPtr_Type                                     M_fluid;
    solidPtr_Type                                     M_solid;
    meshMotionPtr_Type                                M_meshMotion;

    std::string                                       M_fluidTimeAdvanceMethod;
    std::string                                       M_solidTimeAdvanceMethod;
    std::string                                       M_ALETimeAdvanceMethod;

    std::shared_ptr<TimeAdvance<vector_Type> >      M_fluidTimeAdvance;
    std::shared_ptr<TimeAdvance<vector_Type> >      M_fluidMassTimeAdvance;
    std::shared_ptr<TimeAdvance<vector_Type> >      M_solidTimeAdvance;
    std::shared_ptr<TimeAdvance<vector_Type> >      M_ALETimeAdvance;


    dataFile_Type                                     M_dataFile;

    std::shared_ptr<MeshData>                       M_meshDataFluid;
    std::shared_ptr<MeshData>                       M_meshDataSolid;

    dataPtr_Type                                      M_data;

    std::shared_ptr<MapEpetra>                      M_fluidInterfaceMap;
    std::shared_ptr<MapEpetra>                      M_solidInterfaceMap;

    //!\todo{kill this attribute}
    std::shared_ptr<MapEpetra>                      M_fluidInterfaceMapOnZero;
    //!\todo{kill this attribute}
    std::shared_ptr<MapEpetra>                      M_solidInterfaceMapOnZero;

    dofInterface3DPtr_Type                            M_dofFluidToStructure; // Needed
    // dofInterface3DPtr_Type                         M_dofSolidToFluid;
    dofInterface3DPtr_Type                            M_dofStructureToFluid; // Needed
    dofInterface3DPtr_Type                            M_dofStructureToSolid; // Needed to construct M_bcvStructureDispToSolid
    dofInterface3DPtr_Type                            M_dofStructureToHarmonicExtension; // Needed to construct interface maps
    dofInterface3DPtr_Type                            M_dofHarmonicExtensionToFluid; // Needed to construct M_bcvStructureToFluid
    // dofInterface3DPtr_Type                         M_dofStructureToReducedFluid;
    // dofInterface3DPtr_Type                         M_dofReducedFluidToStructure;

    dofInterface2DPtr_Type                            M_dofFluid;
    dofInterface2DPtr_Type                            M_dofSolid;
    dofInterface2DPtr_Type                            M_dofFluidInv;
    dofInterface2DPtr_Type                            M_dofSolidInv;

    bcVectorInterfacePtr_Type                         M_bcvFluidInterfaceDisp;
    bcVectorInterfacePtr_Type                         M_bcvFluidLoadToStructure;
    bcVectorInterfacePtr_Type                         M_bcvSolidLoadToStructure;
    bcVectorInterfacePtr_Type                         M_bcvStructureToFluid;
    bcVectorInterfacePtr_Type                         M_bcvStructureDispToFluid;
    bcVectorInterfacePtr_Type                         M_bcvStructureDispToSolid;
    bcVectorInterfacePtr_Type                         M_bcvStructureDispToHarmonicExtension;
    bcVectorInterfacePtr_Type                         M_bcvHarmonicExtensionVelToFluid;
    // bcVectorInterfacePtr_Type                      M_bcvStructureToReducedFluid;
    // bcVectorInterfacePtr_Type                      M_bcvReducedFluidToStructure;

    bcVectorInterfacePtr_Type                         M_bcvDerHarmonicExtensionVelToFluid;
    bcVectorInterfacePtr_Type                         M_bcvDerFluidLoadToStructure;
    bcVectorInterfacePtr_Type                         M_bcvDerFluidLoadToFluid;
    bcVectorInterfacePtr_Type                         M_bcvDerStructureDispToSolid;
    BCFunctionRobin                                   M_bcfRobinOuterWall;

    // bcVectorInterfacePtr_Type                      M_bcvDerReducedFluidLoadToStructure;
    // bcVectorInterfacePtr_Type                      M_bcvDerStructureAccToReducedFluid;

    vectorPtr_Type                                    M_lambdaFluid;
    vectorPtr_Type                                    M_lambdaFluidRepeated;
    vectorPtr_Type                                    M_lambda;
    vectorPtr_Type                                    M_lambdaDot;


    vectorPtr_Type                                    M_rhs;
    vectorPtr_Type                                    M_alphaF;

    Real                                              M_alphaFCoef;
    //\todo{try to set as deprecated}
    Real                                              M_betaMean;

    commPtr_Type                                      M_epetraComm;
    commPtr_Type                                      M_epetraWorldComm;

    bool                                              M_structureNonLinear;
    //@}
private:

    //!@name Private Methods
    //@{
    //!Private Copy Constructor
    FSIOperator ( const FSIOperator& /*copy*/) {}
    //@}

    //! @name Private Attributes
    //@{
    // displacement on the interface
    vectorPtr_Type                                    M_lambdaSolid;
    vectorPtr_Type                                    M_lambdaSolidRepeated;

    vectorPtr_Type                                    M_lambdaSolidOld;
    vectorPtr_Type                                    M_lambdaDotSolid;
    vectorPtr_Type                                    M_lambdaDotSolidRepeated;

    vectorPtr_Type                                    M_sigmaFluid;
    vectorPtr_Type                                    M_sigmaSolid;

    vectorPtr_Type                                    M_sigmaFluidRepeated;
    vectorPtr_Type                                    M_sigmaSolidRepeated;

    //\todo{try to remove}
    vectorPtr_Type                                    M_minusSigmaFluid;
    //\todo{try to remove}
    vectorPtr_Type                                    M_minusSigmaFluidRepeated;

    vectorPtr_Type                                    M_dispFluidMeshOld;
    vectorPtr_Type                                    M_veloFluidMesh;
    vectorPtr_Type                                    M_derVeloFluidMesh;

    //\todo{try to set as deprecated}
    bool                                              M_mpi;

    bool                                              M_isFluid;
    bool                                              M_isSolid;

    bool                                              M_linearFluid;
    bool                                              M_linearSolid;

    int                                               M_fluidLeader;
    int                                               M_solidLeader;

    std::string                                       M_aleOrder;
    //@}

};

} // Namespace LifeV

#endif /* FSIOPERATOR_H */
