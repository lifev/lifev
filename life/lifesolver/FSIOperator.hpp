/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef TWODIM

#ifndef __FSIOPERATOR_H
#define __FSIOPERATOR_H

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>
#include <life/lifesolver/OseenShapeDerivative.hpp>
//#include <life/lifesolver/Oseen.hpp>
#include <life/lifesolver/VenantKirchhofSolver.hpp>
#include <life/lifesolver/HarmonicExtensionSolver.hpp>

#include <life/lifemesh/regionMesh3D.hpp>

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

namespace LifeV {

enum Preconditioner
{
    NO_PRECONDITIONER = -1,
    NEUMANN_DIRICHLET,
    DIRICHLET_NEUMANN,
    NEUMANN_NEUMANN,
    NEWTON
};

enum DDNPreconditioner
{
    DDN_NO_PRECONDITIONER = -1,
    DDN_NEUMANN_DIRICHLET,
    DDN_DIRICHLET_NEUMANN,
    DDN_NEUMANN_NEUMANN
};

/*!
  \class FSIOperator
  \brief Fluid-Structure Interface operator class
*/
class FSIOperator
{
public:

    /** @name Typedefs
     */
    //@{

//     typedef VenantKirchhofSolver   < RegionMesh3D_ALE<LinearTetra> > solid_raw_type;
//     typedef NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid_raw_type;
    typedef RegionMesh3D<LinearTetra>              mesh_type;

    typedef OseenShapeDerivative   <mesh_type>     fluid_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solid_raw_type;
    typedef HarmonicExtensionSolver<mesh_type>     meshmotion_raw_type;

    typedef OseenShapeDerivative   <mesh_type>     fluidlin_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solidlin_raw_type;

    typedef boost::shared_ptr<fluid_raw_type>      fluid_type;
    typedef boost::shared_ptr<solid_raw_type>      solid_type;
    typedef boost::shared_ptr<meshmotion_raw_type> meshmotion_type;

    typedef boost::shared_ptr<fluidlin_raw_type>   fluidlin_type;
    typedef boost::shared_ptr<solidlin_raw_type>   solidlin_type;

    typedef fluid_raw_type::vector_type            vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;

    typedef fluid_raw_type::source_type            fluid_source_type;
    typedef solid_raw_type::source_type            solid_source_type;

    typedef boost::function<Real ( const Real&, const Real&,
                                   const Real&, const Real&, const ID& )> function_type;

    typedef Real ( *bc_function_type ) ( const Real&, const Real&,
                                         const Real&, const Real&, const ID& );

    typedef boost::shared_ptr<DofInterface3Dto3D>  dof_interface_type3D;
    typedef boost::shared_ptr<DofInterface3Dto2D>  dof_interface_type2D;
    typedef boost::shared_ptr<BCVectorInterface>   bc_vector_interface;

    typedef fluid_raw_type::bchandler_type         fluid_bchandler_type;
    typedef fluid_raw_type::bchandler_raw_type     fluid_bchandler_raw_type;

    typedef solid_raw_type::bchandler_type         solid_bchandler_type;
    typedef solid_raw_type::bchandler_raw_type     solid_bchandler_raw_type;

    typedef FSIOperator::fluid_type::value_type::data_type    data_fluid;
    typedef FSIOperator::solid_type::value_type::data_type    data_solid;

    typedef std::map<ID, ID>::const_iterator                        Iterator;

//     typedef boost::shared_ptr<reducedLinFluid>    quasi_newton_type;

    //@}



    /** @name Constructors, destructor
     */
    //@{

    FSIOperator();

    virtual ~FSIOperator();

    //@}


    /** @name Methods
     */
    //@{

    virtual void setDataFromGetPot( const GetPot& data );

    virtual void setupFEspace();

    virtual void setupDOF();

    virtual void setupFluidSolid();

    virtual void setupSystem();

    virtual void buildSystem();

    virtual void updateSystem( );

    virtual void couplingVariableExtrap( );

    virtual void solveJac( vector_type&       _muk,
                           const vector_type& _res,
                           const Real         _linearRelTol ) = 0;

    virtual void evalResidual( vector_type&        res,
                               const vector_type& _disp,
                               const UInt         _iter ) = 0;

    virtual void shiftSolution();

    virtual void initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::fluid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::fluid_type::value_type::Function const& /*df0=FSIOperator::solid_type::value_type::Function()*/)
    {
        Debug( 6220 ) << "FSIOperator:: solid init \n";
        if (this->isSolid())
            solid().initialize(d0, w0);
        Debug( 6220 ) << "FSIOperator:: fluid init \n";
        if (this->isFluid())
            fluid().initialize(u0, p0);
    }


    void createInterfaceMaps(dof_interface_type3D dofStructureToHarmonicExtension);

    void initializeFluid( const vector_type& velAndPressure,
                          const vector_type& displacement );

    void initializeSolid( vector_ptrtype displacement,
                          vector_ptrtype velocity );

    void updateJacobian ( const vector_type& sol, const int& iter );

    void moveMesh       ( const vector_type& dep );

//     void solveLinearFluid();
//     void solveLinearSolid();

  virtual void setFluxBC             (fluid_bchandler_type /*bc_fluid*/){}

  virtual void setRobinBC             (fluid_bchandler_type /*bc_solid*/){}

    void transferFluidOnInterface( const vector_type& _vec1, vector_type& _vec2 );

    //works in serial but no yet in parallel
    void transferSolidOnFluid    ( const vector_type& _vec1, vector_type& _vec2 );
    void transferSolidOnInterface( const vector_type& _vec1, vector_type& _vec2 );
//     void transferInterfaceOnFluid( const vector_type& _vec1, vector_type& _vec2 );
    void transferInterfaceOnSolid( const vector_type& _vec1, vector_type& _vec2 );

    //! MONOLITHIC Solver methods - Implemented there
    //     virtual boost::shared_ptr<EpetraMap>& monolithicMap()        { assert(false); };
    virtual void iterateMesh( const vector_type& /*disp*/ )     { assert(false); }
    virtual vector_ptrtype const& un()              {assert(false); }
    virtual  void initialize( vector_ptrtype /*u0*/){assert(false); }
    virtual void setupBDF(vector_type const& /*u0*/){}


    //@}



    /** @name  Display Methods
     */
    //@{

    bool isLeader() const;

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

    // now residual is Ax-b, but we should decide for b-Ax. In the meantime, we need b-Ax:
    const vector_type& minusSigmaFluid()                          const { return *M_minusSigmaFluid; }
    const vector_type& minusSigmaFluidRepeated()                  const { return *M_minusSigmaFluidRepeated; }

    string             algorithm()	                              const { return M_algorithm;}
    std::string       method()                                   const { return M_method; }

    vector_type&       Alphaf()                                   const { return *M_Alphaf;}
    Real               time()                                     const { return M_time; }
    const UInt&        nbEval()                                   const { return M_nbEval; }

    Preconditioner     preconditioner()                           const { return M_precond; }
    DDNPreconditioner  DDNpreconditioner()                        const { return M_DDNprecond; }

    Epetra_Comm &       worldComm()                                const { return *M_epetraWorldComm; }

    bool isFluid()                                                const { return M_isFluid; }
    bool isSolid()                                                const { return M_isSolid; }

    bool isLinearFluid()                                          const { return M_linearFluid; }
    bool isLinearSolid()                                          const { return M_linearSolid; }

    int getFluidLeaderId()                                        const { return M_fluidLeader; }
    int getSolidLeaderId()                                        const { return M_solidLeader; }

    const fluid_raw_type& fluid()                                 const { return *M_fluid; }
    const solid_raw_type& solid()                                 const { return *M_solid; }
    const meshmotion_raw_type& meshMotion()                       const { return *M_meshMotion; }

    fluid_type::value_type& fluid()                                     { return *M_fluid; }
    solid_type::value_type& solid()                                     { return *M_solid; }
    meshmotion_type::value_type& meshMotion()                           { return *M_meshMotion; }
//     fluidlin_type::value_type& fluidLin()                               { return *M_fluidLin; }
//     solidlin_type::value_type& solidLin()                               { return *M_solidLin; }

    const data_fluid& dataFluid()                                 const { return *M_dataFluid; }
    const data_solid& dataSolid()                                 const { return *M_dataSolid; }

    const partitionMesh< mesh_type >& fluidMesh()                 const { return *M_fluidMeshPart; }
    const partitionMesh< mesh_type >& solidMesh()                 const { return *M_solidMeshPart; }

    const FESpace<mesh_type, EpetraMap>& uFESpace()               const { return *M_uFESpace; }
    const FESpace<mesh_type, EpetraMap>& pFESpace()               const { return *M_pFESpace; }
    const FESpace<mesh_type, EpetraMap>& dFESpace()               const { return *M_dFESpace; }
    const FESpace<mesh_type, EpetraMap>& mmFESpace()              const { return *M_mmFESpace; }
    virtual const vector_type& meshDisp()                         const { return M_meshMotion->disp(); }
    const         vector_type& dispFluidMeshOld()                 const { return *M_dispFluidMeshOld; }
    virtual       vector_type& veloFluidMesh()                          { return *M_veloFluidMesh; }
                  vector_type& derVeloFluidMesh()                       { return *M_derVeloFluidMesh; }

    //const dof_interface_type3D& dofStructureToFluid()             const { return M_dofSolidToFluid; }
    const dof_interface_type3D& dofFluidToStructure()             const { return M_dofFluidToStructure; }
//     const dof_interface_type3D& dofStructureToFluid()             const { return M_dofStructureToFluid; }
    const dof_interface_type3D& dofStructureToSolid()             const { return M_dofStructureToSolid; }
    const dof_interface_type3D& dofStructureToHarmonicExtension() const { return M_dofStructureToHarmonicExtension; }
    const dof_interface_type3D& dofHarmonicExtensionToFluid()     const { return M_dofHarmonicExtensionToFluid; }

//     std::vector<int>& dofInterfaceFluid()                               { return M_dofInterfaceFluid; }
//     std::vector<int>& dofInterfaceSolid()                               { return M_dofInterfaceSolid; }

            boost::shared_ptr<EpetraMap>& fluidInterfaceMap()              { return M_fluidInterfaceMap; }
            boost::shared_ptr<EpetraMap>& solidInterfaceMap()              { return M_solidInterfaceMap; }
    virtual boost::shared_ptr<EpetraMap>& getCouplingVariableMap()            { return M_solidInterfaceMap; }

    BCFunctionMixte& bcfMixteOuterWall()                                   { return M_bcfMixteOuterWall; }

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

    const fluid_bchandler_type& BCh_fluid()                       const { return M_BCh_u; }
    const fluid_bchandler_type& BCh_harmonicExtension()           const { return M_BCh_mesh; }
    const fluid_bchandler_type& BCh_du()                          const { return M_BCh_du; }
    const fluid_bchandler_type& BCh_du_inv()                      const { return M_BCh_du_inv; }
    const solid_bchandler_type& BCh_solid()                       const { return M_BCh_d; }
    const solid_bchandler_type& BCh_dz()                          const { return M_BCh_dz; }
    const solid_bchandler_type& BCh_dz_inv()                      const { return M_BCh_dz_inv; }
    virtual const fluid_bchandler_type& BCh_flux()                const { }


//     quasi_newton_type getReducedLinFluid()                              { return M_reducedLinFluid; }
//     UInt reducedFluid()                                                 { return M_reducedFluid; }

    //@}



    /** @name  Set Functions
     */
    //@{

    void setComm( const boost::shared_ptr<Epetra_MpiComm>& comm, const boost::shared_ptr<Epetra_MpiComm>& worldComm );

    void setPreconditioner   ( const Preconditioner&    P ) { M_precond    = P; }
    void setDDNPreconditioner( const DDNPreconditioner& P ) { M_DDNprecond = P; }

    void setTime( const Real& time );

    void setFluid( const fluid_type& fluid, const meshmotion_type& meshmotion );
    void setSolid( const solid_type& solid );

    void setFluid( const bool& isFluid ) { M_isFluid = isFluid; }
    void setSolid( const bool& isSolid ) { M_isSolid = isSolid; }

    void setLinearFluid( const bool& linFluid ) { M_linearFluid = linFluid; }
    void setLinearSolid( const bool& linSolid ) { M_linearSolid = linSolid; }

    void setFluidLeader( const int& fluidLeader ) { M_fluidLeader = fluidLeader; }
    void setSolidLeader( const int& solidLeader ) { M_solidLeader = solidLeader; }

//     void setBC( fluid_bchandler_type& bc_u, solid_bchandler_type& bc_d, fluid_bchandler_type& bc_m );

    void setFluidBC             ( const fluid_bchandler_type& bc_fluid );
    void setLinFluidBC          ( const fluid_bchandler_type& bc_dfluid )     { M_BCh_du     = bc_dfluid; }
    void setInvLinFluidBC       ( const fluid_bchandler_type& bc_dfluid_inv ) { M_BCh_du_inv = bc_dfluid_inv; }
    void setHarmonicExtensionBC ( const fluid_bchandler_type& bc_he );

    void setSolidBC             ( const solid_bchandler_type& bc_solid );
    void setLinSolidBC          ( const solid_bchandler_type& bc_dsolid )     { M_BCh_dz     = bc_dsolid; }
    void setInvLinSolidBC       ( const solid_bchandler_type& bc_dsolid_inv ) { M_BCh_dz_inv = bc_dsolid_inv; }

    void setLambdaFluid         ( const vector_type& lambda );
    void setLambdaSolid         ( const vector_type& lambda );

    void setLambdaSolidOld      ( const vector_type& lambda );
    void setLambdaDotSolid      ( const vector_type& lambda );

    void setSigmaFluid          ( const vector_type& sigma );
    void setSigmaSolid          ( const vector_type& sigma );
    void setMinusSigmaFluid     ( const vector_type& sigma );

    void setAlphaf     () { M_Alphaf->getEpetraVector().PutScalar( M_AlphafCoef ); }
    void setAlphafCoef ();
    void setAlphafbcf  ( const bc_function_type& alphafbcf );

//     void setMpi     (bool mpi  ){M_mpi      = mpi;}
//     void setFluidMpi(bool fluid){M_isFluidMpi = fluid;}
//     void setSolidMpi(bool solid){M_issolidMpi = solid;}

//     bool mpi(){return M_mpi;}

//     void setMixteOuterWall                   ( function_type const& dload,
//                                                function_type const& E ) { M_bcfMixteOuterWall.setFunctions_Mixte(dload, E); }

    void setStructureToFluidParametres();

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
    void setMixteOuterWall(function_type const& dload, function_type const& E);

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

    //! gets the solution vector by reference
    virtual void getSolution                  (vector_ptrtype& lambda){lambda=M_lambda;}

    //! sets the solution vector by reference
    void setSolutionPtr                             (const vector_ptrtype& sol)
    {M_lambda = sol;}

    //! sets the solution vector by copy
    void setSolution                                (const vector_type& sol)
    {*M_lambda = sol;}

    //! sets the solution time derivative vector by copy
    virtual  void setSolutionDerivative                                (const vector_ptrtype& lambdaDot){ M_lambdaDot=lambdaDot; }

    //! gets the solid displacement by copy
    virtual  void getSolidDisp                               (vector_type& soliddisp)
    { soliddisp = M_solid->disp();}

    //! gets the solid velocity by copy
    virtual  void getSolidVel                                (vector_type& solidvel)
    { solidvel = M_solid->vel();}

    //! gets the fluid velocity and pressure by copy
    virtual  void getFluidVelAndPres                         (vector_type& sol)
    { sol = M_fluid->solution();}


protected:

    virtual UInt imposeFlux();

    //virtual void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);
    virtual void variablesInit( const std::string& dOrder );

    void transferMeshMotionOnFluid( const vector_type& _vec1,
                                          vector_type& _vec2 );

    void interpolateVelocity(const vector_type& _vec1,
                                   vector_type& _vec2);


    void interpolateInterfaceDofs(const FESpace<mesh_type, EpetraMap>& _fespace1,
                                  const vector_type&                   _vec1,
                                  const FESpace<mesh_type, EpetraMap>& _fespace2,
                                  vector_type&                         _vec2,
                                  dof_interface_type3D&                _dofInterface);

    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_uFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_pFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_dFESpace;
    boost::shared_ptr<FESpace<mesh_type, EpetraMap> > M_mmFESpace;

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
    boost::shared_ptr<data_fluid>                     M_dataFluid;
    boost::shared_ptr<data_solid>                     M_dataSolid;

//     quasi_newton_type         M_reducedLinFluid;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMap;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMapOnZero;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMapOnZero;

//     std::vector<int>                                  M_dofInterfaceFluid;
//     std::vector<int>                                  M_dofInterfaceSolid;

    dof_interface_type3D                              M_dofFluidToStructure;
//     dof_interface_type3D                              M_dofSolidToFluid;
    dof_interface_type3D                              M_dofStructureToFluid;
    dof_interface_type3D                              M_dofStructureToSolid;
    dof_interface_type3D                              M_dofStructureToHarmonicExtension;
    dof_interface_type3D                              M_dofHarmonicExtensionToFluid;
//     dof_interface_type3D                              M_dofStructureToReducedFluid;
//     dof_interface_type3D                              M_dofReducedFluidToStructure;

    int                                               M_fluidInterfaceFlag;
    int                                               M_solidInterfaceFlag;
    int                                               M_structureInterfaceFlag;
    int                                               M_harmonicInterfaceFlag;
    Real                                              M_interfaceTolerance;

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

private:
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

    boost::shared_ptr<vector_type>                    M_minusSigmaFluid;
    boost::shared_ptr<vector_type>                    M_minusSigmaFluidRepeated;

//     boost::shared_ptr<vector_type>                    M_dispStruct;
//     boost::shared_ptr<vector_type>                    M_dispStructOld;
//     boost::shared_ptr<vector_type>                    M_veloStruct;

    boost::shared_ptr<vector_type>                    M_dispFluidMeshOld;
    boost::shared_ptr<vector_type>                    M_veloFluidMesh;
    boost::shared_ptr<vector_type>                    M_derVeloFluidMesh;

protected:

    boost::shared_ptr<vector_type>                    M_un;
    boost::shared_ptr<vector_type>                    M_rhs;
    boost::shared_ptr<vector_type>                    M_Alphaf;

    Real                                              M_AlphafCoef;
    Real                                              M_betamedio;
    Real                                              M_time;
//     boost::shared_ptr<vector_type>                    M_w;
//     boost::shared_ptr<vector_type>                    M_dw;

    UInt                                              M_nbEval;

    boost::shared_ptr<Epetra_Comm>                    M_epetraComm;
    boost::shared_ptr<Epetra_Comm>                    M_epetraWorldComm;
    std::string                                       M_method;

private:

//     UInt                      M_reducedFluid;
    std::string                                       M_algorithm;
    Preconditioner                                    M_precond;
    DDNPreconditioner                                 M_DDNprecond;

    bool                                              M_mpi;
    bool                                              M_isFluid;
    bool                                              M_isSolid;

    bool                                              M_linearFluid;
    bool                                              M_linearSolid;

    int                                               M_fluidLeader;
    int                                               M_solidLeader;
};





// ===================================================
//! MACROS
// ===================================================
typedef boost::shared_ptr<FSIOperator>                 oper_fsi_ptr_mpi;
typedef singleton<factory<FSIOperator,  std::string> > FSIFactory;

/*!
   \def  FOR_EACH_INTERFACE_DOF( Expr )

   \c FOR_EACH_INTERFACE_DOF is an helper macro to ease the
   the computation of quantities like the residual on the interface
   of the Fluid and Structure.

 */
//    UInt iBCf = M_fluid->BCh_fluid().getBCbyName("Wall");



#define FOR_EACH_INTERFACE_DOF( Expr )                              \
{   \
    \
    \
    UInt iBCf = M_harmonicExtension->bcHandler().getBCbyName("Interface"); \
    BCBase const &BC_fluidInterface = M_harmonicExtension->bcHandler()[iBCf]; \
    UInt nDofInterface = BC_fluidInterface.list_size();             \
    UInt nDimF = BC_fluidInterface.numberOfComponents();            \
                                                                    \
    UInt totalDofFluid = M_fluid->velFESpace().dof().numTotalDof();             \
    UInt totalDofSolid = M_solid->dFESpace().dof().numTotalDof();             \
                                                                    \
    for (UInt iBC = 1; iBC <= nDofInterface; ++iBC)                 \
    {                                                               \
        ID IDfluid = BC_fluidInterface(iBC)->id();                  \
                                                                    \
        BCVectorInterface const *BCVInterface =                     \
            dynamic_cast <BCVectorInterface const *>                \
            (BC_fluidInterface.pointerToBCVector());                \
                                                                    \
        assert( BCVInterface != 0 );                                \
                                                                    \
        ID IDsolid = BCVInterface->                                 \
            dofInterface().getInterfaceDof(IDfluid);                \
        for (UInt jDim = 0; jDim < nDimF; ++jDim)                   \
        {                                                           \
            ( Expr );                                               \
        }                                                           \
    }                                                               \
}



#define FOR_EACH_INTERFACE_DOF_SOLID( Expr )                              \
{   \
    \
    \
\
    std::map<ID, ID> solidDofMap = M_dofSolid->locDofMap();\
\
    UInt totalDofSolid = M_solid->dFESpace().dof().numTotalDof();\
\
    for (std::map<ID,ID>::const_iterator it = solidDofMap.begin(); it != solidDofMap.end(); ++it)\
    {\
        ID dofS  = it->first - 1;\
        ID dofFS = it->second - 1;\
\
        for (UInt jDim = 0; jDim < 3; ++jDim)\
        {\
            ( Expr );\
        }\
    }\
}



#define FOR_EACH_INTERFACE_DOF_FLUID( Expr )                              \
{   \
    \
    \
\
    std::map<ID, ID> fluidDofMap = M_dofFluid->locDofMap();\
\
    UInt totalDofFluid = M_fluid->velFESpace().dof().numTotalDof();\
\
    for (std::map<ID,ID>::const_iterator it = fluidDofMap.begin(); it != fluidDofMap.end(); ++it)\
    {\
        ID dofF  = it->first - 1;\
        ID dofFS = it->second - 1;\
        for (UInt jDim = 0; jDim < 3; ++jDim)\
        {\
            ( Expr );\
        }\
    }\
}

} // Namespace LifeV

#endif /* __FSIOPERATOR_H */
#endif /* TWODIM */
