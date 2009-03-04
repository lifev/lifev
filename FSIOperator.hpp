/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politechnico di Milano

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

//! previously called operFS in the CVS repository

#ifndef _OPERFS_H
#define _OPERFS_H

#include <life/lifecore/factory.hpp>
#include <life/lifecore/singleton.hpp>

#include <life/lifefem/dofInterface3Dto3D.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>
#include <lifemc/lifesolver/OseenShapeDerivative.hpp>
//#include <life/lifesolver/Oseen.hpp>
#include <lifemc/lifesolver/VenantKirchhofSolver.hpp>
#include <lifemc/lifesolver/HarmonicExtensionSolver.hpp>
#include <lifemc/lifemesh/regionMesh3D.hpp>
#include <life/lifealg/generalizedAitken.hpp>
#include <life/lifefem/bcHandler.hpp>
#include <life/lifefem/bcFunction.hpp>
#include <life/lifefem/bdf_template.hpp>
//#include <life/lifefem/dof.hpp>
#include <life/lifefem/FESpace.hpp>
//#include <life/lifefem/refFE.hpp>
//#include <life/lifesolver/fixedPointBase.hpp>

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif


#define FLUID 1
#define SOLID 0


namespace LifeV
{

//Fluid-Structure operator Class


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

//class reducedLinFluid;

class FSIOperator {

public:

//     typedef VenantKirchhofSolver   < RegionMesh3D_ALE<LinearTetra> > solid_raw_type;
//     typedef NavierStokesAleSolverPC< RegionMesh3D_ALE<LinearTetra> > fluid_raw_type;
    typedef RegionMesh3D<LinearTetra>              mesh_type;

    typedef OseenShapeDerivative   <mesh_type>     fluid_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solid_raw_type;
    typedef HarmonicExtensionSolver<mesh_type>     meshmotion_raw_type;

    typedef OseenShapeDerivative   <mesh_type>     fluidlin_raw_type;
    typedef VenantKirchhofSolver   <mesh_type>     solidlin_raw_type;

    typedef fluid_raw_type::vector_type            vector_type;
    typedef boost::shared_ptr<vector_type>         vector_ptrtype;

    typedef fluid_raw_type::source_type            fluid_source_type;
    typedef solid_raw_type::source_type            solid_source_type;

    typedef boost::shared_ptr<fluid_raw_type>      fluid_type;
    typedef boost::shared_ptr<solid_raw_type>      solid_type;

    typedef boost::shared_ptr<fluidlin_raw_type>   fluidlin_type;
    typedef boost::shared_ptr<solidlin_raw_type>   solidlin_type;

    typedef boost::shared_ptr<meshmotion_raw_type> meshmotion_type;

    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID
                                   & )> function_type;
    typedef Real ( *bc_function_type ) ( const Real&, const Real&, const Real&,
                                         const Real&, const ID& );
    typedef boost::shared_ptr<DofInterface3Dto3D>  dof_interface_type3D;
    typedef boost::shared_ptr<DofInterface3Dto2D>  dof_interface_type2D;
    typedef boost::shared_ptr<BCVectorInterface>   bc_vector_interface;

    typedef solid_raw_type::bchandler_type         solid_bchandler_type;
    typedef solid_raw_type::bchandler_raw_type     solid_bchandler_raw_type;

    typedef fluid_raw_type::bchandler_type         fluid_bchandler_type;
    typedef fluid_raw_type::bchandler_raw_type     fluid_bchandler_raw_type;

    typedef FSIOperator::fluid_type::value_type::data_type    data_fluid;
    typedef FSIOperator::solid_type::value_type::data_type    data_solid;

//     typedef boost::shared_ptr<reducedLinFluid>    quasi_newton_type;

    // constructors
    FSIOperator():
        M_uFESpace  (),
        M_pFESpace  (),
        M_dFESpace  (),
        M_mmFESpace (),
        M_fluidMeshPart(),
        M_solidMeshPart(),
        M_BCh_u     (new fluid_bchandler_raw_type),
        M_BCh_d     (new solid_bchandler_raw_type),
        M_BCh_mesh  (new fluid_bchandler_raw_type),
        M_BCh_du    (new fluid_bchandler_raw_type),
        M_BCh_du_inv(new fluid_bchandler_raw_type),
        M_BCh_dz    (new solid_bchandler_raw_type),
        M_BCh_dz_inv(new solid_bchandler_raw_type),
        M_BCh_dp    (new BCHandler),
        M_BCh_dp_inv(new BCHandler),
        M_fluid(),
        M_solid(),
//         M_fluidLin(),
//         M_solidLin(),
        M_meshMotion(),
        M_bdf(),
        M_dataFluid(),
        M_dataSolid(),
        M_fluidInterfaceMap(),
        M_solidInterfaceMap(),
        M_fluidInterfaceMapOnZero(),
        M_solidInterfaceMapOnZero(),
        M_dofFluidToStructure                ( new DofInterface3Dto3D ),
//         M_dofFluidToSolid                    ( new DofInterface3Dto3D ),
//         M_dofSolidToFluid                    ( new DofInterface3Dto3D ),
        M_dofStructureToFluid                ( new DofInterface3Dto3D ),
        M_dofStructureToSolid                ( new DofInterface3Dto3D ),
        M_dofStructureToHarmonicExtension    ( new DofInterface3Dto3D ),
        M_dofHarmonicExtensionToFluid        ( new DofInterface3Dto3D ),
//         M_dofStructureToReducedFluid         ( new DofInterface3Dto3D ),
//         M_dofReducedFluidToStructure         ( new DofInterface3Dto3D ),
        M_fluidInterfaceFlag                 ( 1 ),
        M_solidInterfaceFlag                 ( 1 ),
        M_structureInterfaceFlag             ( 1 ),
        M_harmonicInterfaceFlag              ( 1 ),
        M_interfaceTolerance                 ( 0. ),
        M_dofFluid                           ( new DofInterface3Dto2D ),
        M_dofSolid                           ( new DofInterface3Dto2D ),
        M_dofSolidInv                        ( new DofInterface3Dto2D ),
        M_dofFluidInv                        ( new DofInterface3Dto2D ),
// boundary vector interfaces
        M_bcvStructureDispToFluid            ( new  BCVectorInterface ),
        M_bcvFluidInterfaceDisp              ( new  BCVectorInterface ),
        M_bcvFluidLoadToStructure            ( new  BCVectorInterface ),
        M_bcvSolidLoadToStructure            ( new  BCVectorInterface ),
        M_bcvStructureDispToSolid            ( new  BCVectorInterface ),
        M_bcvStructureDispToHarmonicExtension( new  BCVectorInterface ),
        M_bcvHarmonicExtensionVelToFluid     ( new  BCVectorInterface ),
        M_bcvStructureToFluid                ( new  BCVectorInterface ),
//         M_bcvStructureToReducedFluid         ( new  BCVectorInterface ),
//         M_bcvReducedFluidToStructure         ( new  BCVectorInterface ),
        M_bcvDerHarmonicExtensionVelToFluid  ( new  BCVectorInterface ),
        M_bcvDerFluidLoadToStructure         ( new  BCVectorInterface ),
        M_bcvDerFluidLoadToFluid             ( new  BCVectorInterface ),
        M_bcvDerStructureDispToSolid         ( new  BCVectorInterface ),
//         M_bcvDerReducedFluidLoadToStructure  ( new  BCVectorInterface ),
//         M_bcvDerStructureAccToReducedFluid   ( new  BCVectorInterface ),
        M_lambdaFluid(),
        M_lambdaSolid(),
        M_lambdaFluidRepeated(),
        M_lambdaSolidRepeated(),
        M_lambdaDotSolid(),
        M_lambdaSolidOld(),
        M_lambdaDotSolidRepeated(),
        M_sigmaFluid(),
        M_sigmaSolid(),
        M_sigmaFluidRepeated(),
        M_sigmaSolidRepeated(),
        M_dispFluidMeshOld(),
        M_veloFluidMesh(),
        M_minusSigmaFluid(),    //sigmafluid: auxiliary variable for Robin.
        M_minusSigmaFluidRepeated(), //sigmafluid: auxiliary variable for Robin.
        M_nbEval( 0 ),
        M_AlphafCoef(0),
        M_Alphaf(),                   //vector_type, for alphaf robin
        M_epetraComm(),
        M_epetraWorldComm (),
        M_method(),
        M_precond( NO_PRECONDITIONER ),
        M_mpi(true),
        M_isFluid(false),
        M_isSolid(false)
        {}

    // destructor

    virtual ~FSIOperator();

    // virtual member functions

    virtual void   setup();

    virtual void   evalResidual(vector_type&        res,
                                const vector_type& _disp,
                                const int          _iter) = 0;

    virtual void   solveJac(vector_type&       _muk,
                            const vector_type& _res,
                            const double       _linearRelTol) = 0;


    void initializeFluid( const vector_type& velAndPressure,
                          const vector_type& displacement );

    void initializeSolid( const vector_type& displacement,
                          const vector_type& velocity );



    // unknows on the structure mesh

    //    vector_type & displacement()    { return *M_lambdaSolid; }
    //    vector_type & displacementOld() { return *M_lambdaSolidOld; }
    //vector_type & residual()     ;
    //vector_type const & velocity() const { return *M_lambdaDotSolid; }
    //vector_type & residualFSI()  ;

    vector_type const& lambdaFluid()    const {return *M_lambdaFluid;}
    vector_type const& sigmaFluid()     const {return *M_sigmaFluid;}
    vector_type const& lambdaSolid()    const {return *M_lambdaSolid;}
    vector_type const& lambdaSolidOld()    const {return *M_lambdaSolidOld;}
    vector_type const& sigmaSolid()     const {return *M_sigmaSolid;}
    vector_type const& lambdaDotSolid() const {return *M_lambdaDotSolid;}

    vector_type const& lambdaFluidRepeated()    const {return *M_lambdaFluidRepeated;}
    vector_type const& sigmaFluidRepeated()     const {return *M_sigmaFluidRepeated;}
    vector_type const& lambdaSolidRepeated()    const {return *M_lambdaSolidRepeated;}
    vector_type const& sigmaSolidRepeated()     const {return *M_sigmaSolidRepeated;}
    vector_type const& lambdaDotSolidRepeated() const {return *M_lambdaDotSolidRepeated;}
    // now residual is Ax-b, but we should decide for b-Ax. In the meantime, we need b-Ax:
    vector_type const& minusSigmaFluid()  {return *M_minusSigmaFluid;}
    vector_type const& minusSigmaFluidRepeated()  {return *M_minusSigmaFluidRepeated;}

    // member functions

    void updateJacobian (vector_type& sol,
                         int    iter);

    void updateSystem(fluid_source_type &fluidSource, solid_source_type &solidSource);
    void shiftSolution();

    void moveMesh(vector_type const &dep);

    void solveLinearFluid();
    void solveLinearSolid();

    // mutators and setters

    UInt   const & nbEval()      const
        {return M_nbEval;}

    meshmotion_type::value_type& meshMotion() {return *M_meshMotion;}

    fluid_type::value_type&      fluid()      {return *M_fluid;}
    solid_type::value_type&      solid()      {return *M_solid;}

//     fluidlin_type::value_type&   fluidLin()   {return *M_fluidLin;}
//     solidlin_type::value_type&   solidLin()   {return *M_solidLin;}

    void setPreconditioner   ( Preconditioner    _p ) { M_precond = _p; }
    void setDDNPreconditioner( DDNPreconditioner _p ) { M_DDNprecond = _p; }

    Preconditioner       preconditioner() const { return M_precond; }
    DDNPreconditioner    DDNpreconditioner() const { return M_DDNprecond; }

    std::string                    method() const {return M_method; }

    void setTime(const Real &time)
        {
            M_time = time;
            M_dataFluid->setTime(time);
            M_dataSolid->setTime(time);
        }
    Real time() {return M_time;}

    vector_type displacementOnInterface();

    void setFluid( fluid_type const& fluid, meshmotion_type const& meshmotion )
        { M_fluid = fluid; M_meshMotion = meshmotion; M_isFluid = true;}
    void setSolid( solid_type const& solid )
        { M_solid = solid; M_isSolid = true;}

    void setFluid( bool isFluid)
        { M_isFluid = isFluid;}
    void setSolid( bool isSolid )
        { M_isSolid = isSolid;}

    void setLinearFluid(bool linFluid)
        { M_linearFluid = linFluid; }

    void setLinearSolid(bool linSolid)
        { M_linearSolid = linSolid; }


//     void setMpi     (bool mpi  ){M_mpi      = mpi;}
//     void setFluidMpi(bool fluid){M_isFluidMpi = fluid;}
//     void setSolidMpi(bool solid){M_issolidMpi = solid;}

//    bool mpi(){return M_mpi;}
    bool    isFluid() const {return M_isFluid;}
    bool    isSolid() const {return M_isSolid;}

    void    setUpSystem( GetPot const& data_file );
    virtual void buildSystem();

    void    setComm     (   boost::shared_ptr<Epetra_MpiComm> comm,
                            boost::shared_ptr<Epetra_MpiComm> worldComm);

    Epetra_Comm& worldComm(){ return *M_epetraWorldComm; }

    //

    boost::shared_ptr<EpetraMap>& fluidInterfaceMap() {return M_fluidInterfaceMap;}
    boost::shared_ptr<EpetraMap>& solidInterfaceMap() {return M_solidInterfaceMap;}


    void setFluidLeader(int fluidLeader){M_fluidLeader = fluidLeader;}
    void setSolidLeader(int solidLeader){M_solidLeader = solidLeader;}

    bool isLeader() const
    {
        if (isFluid())
        {
            if (M_fluid.get() == 0)
                return (M_epetraComm->MyPID() == 0);
            return M_fluid->isLeader();
        }
        if (M_solid.get() == 0)
            return (M_epetraComm->MyPID() == 0);
        return M_solid->isLeader();
    }

    void leaderPrint   (string const message, double const number) const;
    void leaderPrint   (string const message) const;
    void leaderPrintMax(string const message, double const number) const;

    virtual void setDataFromGetPot( GetPot const& data );

    void setBC( fluid_bchandler_type& bc_u,  solid_bchandler_type& bc_d, fluid_bchandler_type& bc_m );

    void setFluidBC             (fluid_bchandler_type bc_fluid);

    void setLinFluidBC          (fluid_bchandler_type bc_dfluid);

    void setInvLinFluidBC       (fluid_bchandler_type bc_dfluid_inv);

    void setHarmonicExtensionBC (fluid_bchandler_type bc_he);

    void setSolidBC             (solid_bchandler_type bc_solid);

    void setLinSolidBC          (solid_bchandler_type bc_dsolid);

    void setInvLinSolidBC       (solid_bchandler_type bc_dsolid_inv);


    void setLambdaFluid(const vector_type& lambda);
    void setLambdaSolid(const vector_type& lambda);

    void setLambdaSolidOld(const vector_type& lambda);

    void setLambdaDotSolid(const vector_type& lambda);

    void setSigmaFluid(const vector_type& sigma);
    void setSigmaSolid(const vector_type& sigma);
    void setMinusSigmaFluid(const vector_type& sigma);

    void transferFluidOnInterface(const vector_type& _vec1,
                                  vector_type&       _vec2);

    //works in serial but no yet in parallel
    void transferSolidOnFluid(const vector_type& _vec1,
                                  vector_type&       _vec2);

    void transferSolidOnInterface(const vector_type& _vec1,
                                  vector_type&       _vec2);

    void transferInterfaceOnSolid(const vector_type& _vec1,
                                  vector_type&       _vec2);

    void transferInterfaceOnFluid(const vector_type& _vec1,
                                  vector_type&       _vec2);
    //setting type of algorithm
    string algorithm()	{ return M_algorithm;}

    const vector_type& minusSigmaFluid() const { return *M_minusSigmaFluid;}

//  setting Alphaf:

    void setAlphafCoef();

    void setAlphaf() { M_Alphaf->getEpetraVector().PutScalar(M_AlphafCoef); }

    void setAlphafbcf(const bc_function_type &alphafbcf) {
        vector_type vec( M_fluid->velFESpace().map());
        M_fluid->velFESpace().interpolate(alphafbcf, vec, 0.0);
        *M_Alphaf = vec ;
    }

  vector_type& Alphaf() { return *M_Alphaf;}
//@ getters

    const fluid_raw_type& fluid()      const {return *M_fluid;}
    const solid_raw_type& solid()      const {return *M_solid;}
    const meshmotion_raw_type& meshMotion() const {return *M_meshMotion;}


    const data_fluid& dataFluid() const {return *M_dataFluid;}
    const data_solid& dataSolid() const {return *M_dataSolid;}

    const partitionMesh< mesh_type >&    fluidMesh() const {return *M_fluidMeshPart;}
    const partitionMesh< mesh_type >&    solidMesh() const {return *M_solidMeshPart;}

    const FESpace<mesh_type, EpetraMap>& uFESpace()  const {return *M_uFESpace;}
    const FESpace<mesh_type, EpetraMap>& pFESpace()  const {return *M_pFESpace;}
    const FESpace<mesh_type, EpetraMap>& dFESpace()  const {return *M_dFESpace;}
    const FESpace<mesh_type, EpetraMap>& mmFESpace() const {return *M_mmFESpace;}

    const vector_type& veloFluidMesh() const {return *M_veloFluidMesh;}
    vector_type& veloFluidMesh() {return *M_veloFluidMesh;}
    const vector_type& dispFluidMeshOld() const {return *M_dispFluidMeshOld;}

    const vector_type& derVeloFluidMesh() const {return *M_derVeloFluidMesh;}
    vector_type& derVeloFluidMesh() {return *M_derVeloFluidMesh;}

    const dof_interface_type3D& dofStructureToFluid()             const {return M_dofSolidToFluid;}
    const dof_interface_type3D& dofFluidToStructure()             const {return M_dofFluidToStructure;}
  //const dof_interface_type3D& dofStructureToFluid()             const {return M_dofStructureToFluid;}
    const dof_interface_type3D& dofStructureToSolid()             const {return M_dofStructureToSolid;}
    const dof_interface_type3D& dofStructureToHarmonicExtension() const {return M_dofStructureToHarmonicExtension;}
    const dof_interface_type3D& dofHarmonicExtensionToFluid()     const {return M_dofHarmonicExtensionToFluid;}

    std::vector<int>& dofInterfaceFluid() {return M_dofInterfaceFluid;}
    std::vector<int>& dofInterfaceSolid() {return M_dofInterfaceSolid;}


//     void setReducedLinFluidBC   (fluid_bchandler_type bc_dredfluid);

//     void setInvReducedLinFluidBC(fluid_bchandler_type bc_invdredfluid);


    // BC Vector Interface setters and getters

    void setStructureDispToFluid(vector_type const& vel, UInt type = 0);

    bc_vector_interface bcvStructureDisptoFluid()   {return M_bcvStructureDispToFluid;}

    void setStructureToFluid(vector_type const& vel, UInt type = 0);

    bc_vector_interface bcvStructureToFluid()       {return M_bcvStructureToFluid;}

    void setSolidLoadToStructure(vector_type const& load, UInt type = 0);

    bc_vector_interface bcvSolidLoadToStructure()   {return M_bcvSolidLoadToStructure;}

    void setFluidInterfaceDisp      (vector_type const& disp, UInt type = 0);
    bc_vector_interface bcvFluidInterfaceDisp()
        {return M_bcvFluidInterfaceDisp;}

    void setHarmonicExtensionVelToFluid(vector_type const& vel, UInt type = 0);
    bc_vector_interface bcvHarmonicExtensionVelToFluid()
        {return M_bcvHarmonicExtensionVelToFluid;}

    void setDerHarmonicExtensionVelToFluid(vector_type const& dvel, UInt type = 0);
    bc_vector_interface bcvDerHarmonicExtensionVelToFluid()
        {return M_bcvDerHarmonicExtensionVelToFluid;}

    void setStructureDispToHarmonicExtension(vector_type const& disp, UInt type = 0);
    bc_vector_interface bcvStructureDispToHarmonicExtension()
        {return M_bcvStructureDispToHarmonicExtension;}

    void setStructureDispToSolid(vector_type const& disp, UInt type = 0);
    bc_vector_interface bcvStructureDispToSolid()
        {return M_bcvStructureDispToSolid;}

    void setDerStructureDispToSolid(vector_type const& ddisp, UInt type = 0);
    bc_vector_interface bcvDerStructureDispToSolid()
        {return M_bcvDerStructureDispToSolid;}

    void setFluidLoadToStructure(vector_type const& load, UInt type = 0);
    bc_vector_interface bcvFluidLoadToStructure()
        {return M_bcvFluidLoadToStructure;}

    void setDerFluidLoadToStructure(vector_type const& dload, UInt type = 0);
    bc_vector_interface bcvDerFluidLoadToStructure()
        {return M_bcvDerFluidLoadToStructure;}

    void setDerFluidLoadToFluid(vector_type const& dload, UInt type = 0);
    bc_vector_interface bcvDerFluidLoadToFluid()
        {return M_bcvDerFluidLoadToFluid;}

//     void setDerReducedFluidLoadToStructure(vector_type &dload, UInt type = 0);
//     bc_vector_interface bcvDerReducedFluidLoadToStructure()
//         {return M_bcvDerReducedFluidLoadToStructure;}

//     void setDerStructureAccToReducedFluid(vector_type &acc, UInt type = 0);
//     bc_vector_interface bcvDerStructureAccToReducedFluid()
//         {return M_bcvDerStructureAccToReducedFluid;}



    //


//     quasi_newton_type getReducedLinFluid(){return M_reducedLinFluid;}

//     UInt reducedFluid(){return M_reducedFluid;}

    fluid_bchandler_type const& BCh_fluid(){return M_BCh_u;}
    void setBCh_fluid(fluid_bchandler_type BCh_fluid){M_BCh_u = BCh_fluid;}

    solid_bchandler_type const& BCh_solid(){return M_BCh_d;}
    void setBCh_solid(solid_bchandler_type BCh_solid){M_BCh_d = BCh_solid;}

    fluid_bchandler_type const& BCh_harmonicExtension(){return M_BCh_mesh;}
    void setBCh_HarmonicExtension(fluid_bchandler_type BCh_mesh){M_BCh_mesh = BCh_mesh;}

    fluid_bchandler_type const& BCh_du(){return M_BCh_du;}
    void setBCh_fluidDer(fluid_bchandler_type BCh_fluidDer){M_BCh_du = BCh_fluidDer;}

    fluid_bchandler_type const& BCh_du_inv(){return M_BCh_du_inv;}
    void setBCh_fluidDerInv(fluid_bchandler_type BCh_fluidDerInv){M_BCh_du_inv = BCh_fluidDerInv;}

    solid_bchandler_type const& BCh_dz(){return M_BCh_dz;}
    void setBCh_solidDer(solid_bchandler_type BCh_solidDer){M_BCh_dz = BCh_solidDer;}

    solid_bchandler_type const& BCh_dz_inv(){return M_BCh_dz_inv;}
    void setBCh_solidDerInv(solid_bchandler_type BCh_solidDerInv){M_BCh_dz_inv = BCh_solidDerInv;}

    //! relevant only for monolitic solver. re-Implemented there

    virtual const boost::shared_ptr<EpetraMap>& monolithicMap() const {};
    //used in Epetra_fullMonolithic
    //! relevant only for monolitic solver. re-Implemented there
    virtual void updateSystem(const vector_type& /*displacement*/) { assert(false); }

    //! relevant only for monolitic solver. re-Implemented there

    virtual void iterateMesh(const vector_type& disp) {}

protected:

    void transferMeshMotionOnFluid(const vector_type &_vec1,
                                   vector_type       &_vec2);

    void interpolateVelocity(const vector_type& _vec1,
                             vector_type& _vec2);

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

//     fluidlin_type                                     M_fluidLin;
//     solidlin_type                                     M_solidLin;

    meshmotion_type                                   M_meshMotion;

    boost::shared_ptr<BdfT<vector_type> >             M_bdf;

    boost::shared_ptr<data_fluid>                     M_dataFluid;
    boost::shared_ptr<data_solid>                     M_dataSolid;

//     quasi_newton_type         M_reducedLinFluid;

    std::vector<int>                                  M_dofInterfaceFluid;
    std::vector<int>                                  M_dofInterfaceSolid;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMap;

    boost::shared_ptr<EpetraMap>                      M_fluidInterfaceMapOnZero;
    boost::shared_ptr<EpetraMap>                      M_solidInterfaceMapOnZero;

    dof_interface_type3D                              M_dofSolidToFluid;
    dof_interface_type3D                              M_dofFluidToStructure;
    dof_interface_type3D                              M_dofStructureToFluid;
    dof_interface_type3D                              M_dofStructureToSolid;
    dof_interface_type3D                              M_dofStructureToHarmonicExtension;
    dof_interface_type3D                              M_dofHarmonicExtensionToFluid;
//     dof_interface_type3D      M_dofStructureToReducedFluid;
//     dof_interface_type3D      M_dofReducedFluidToStructure;

    int                                              M_fluidInterfaceFlag;
    int                                              M_solidInterfaceFlag;
    int                                              M_structureInterfaceFlag;
    int                                              M_harmonicInterfaceFlag;
    Real                                              M_interfaceTolerance;

    dof_interface_type2D                              M_dofFluid;
    dof_interface_type2D                              M_dofSolid;
    dof_interface_type2D                              M_dofSolidInv;
    dof_interface_type2D                              M_dofFluidInv;

    bc_vector_interface                               M_bcvFluidInterfaceDisp;
    bc_vector_interface                               M_bcvFluidLoadToStructure;
    bc_vector_interface                               M_bcvStructureDispToSolid;
    bc_vector_interface                               M_bcvStructureDispToHarmonicExtension;
    bc_vector_interface                               M_bcvHarmonicExtensionVelToFluid;
    bc_vector_interface                               M_bcvStructureToFluid;
    bc_vector_interface                               M_bcvStructureDispToFluid;
    bc_vector_interface                               M_bcvSolidLoadToStructure;
//     bc_vector_interface       M_bcvStructureToReducedFluid;
//     bc_vector_interface       M_bcvReducedFluidToStructure;

    bc_vector_interface                               M_bcvDerHarmonicExtensionVelToFluid;
    bc_vector_interface                               M_bcvDerFluidLoadToStructure;
    bc_vector_interface                               M_bcvDerFluidLoadToFluid;
    bc_vector_interface                               M_bcvDerStructureDispToSolid;
//     bc_vector_interface       M_bcvDerReducedFluidLoadToStructure;
//     bc_vector_interface       M_bcvDerStructureAccToReducedFluid;


private:
    // displacement on the interface
    boost::shared_ptr<vector_type>                    M_lambdaFluid;
    boost::shared_ptr<vector_type>                    M_lambdaSolid;

    boost::shared_ptr<vector_type>                    M_lambdaFluidRepeated;
    boost::shared_ptr<vector_type>                    M_lambdaSolidRepeated;

    boost::shared_ptr<vector_type>                    M_lambdaDotSolid;
    boost::shared_ptr<vector_type>                    M_lambdaSolidOld;

    boost::shared_ptr<vector_type>                    M_lambdaDotSolidRepeated;

    boost::shared_ptr<vector_type>                    M_sigmaFluid;
    boost::shared_ptr<vector_type>                    M_sigmaSolid;

    boost::shared_ptr<vector_type>                    M_sigmaFluidRepeated;
    boost::shared_ptr<vector_type>                    M_sigmaSolidRepeated;

    boost::shared_ptr<vector_type>                    M_minusSigmaFluid;
    boost::shared_ptr<vector_type>                    M_minusSigmaFluidRepeated;

    //    boost::shared_ptr<vector_type>                    M_dispStruct;
    //    boost::shared_ptr<vector_type>                    M_dispStructOld;
    //    boost::shared_ptr<vector_type>                    M_veloStruct;


    boost::shared_ptr<vector_type>                    M_dispFluidMeshOld;
    boost::shared_ptr<vector_type>                    M_veloFluidMesh;
    boost::shared_ptr<vector_type>                    M_derVeloFluidMesh;

protected:

    std::string const& getMethod() const  {return M_method;}
    bool isLinearFluid() const {return M_linearFluid;}
    bool isLinearSolid() const {return M_linearSolid;}

    boost::shared_ptr<vector_type>                    M_un;
    boost::shared_ptr<vector_type>                    M_rhs;
    boost::shared_ptr<vector_type>                    M_Alphaf;

//     boost::shared_ptr<vector_type>                    M_w;
//     boost::shared_ptr<vector_type>                    M_dw;

    Real                                              M_AlphafCoef;
    Real                                              M_betamedio;

    Real                                              M_time;

    UInt                                              M_nbEval;

    boost::shared_ptr<Epetra_Comm>                    M_epetraComm;
    boost::shared_ptr<Epetra_Comm>                    M_epetraWorldComm;
    std::string                                       M_method;

private:

//     UInt                      M_reducedFluid;
    std::string                                       M_algorithm;
    bool                                              M_monolithic;
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

}

#endif
