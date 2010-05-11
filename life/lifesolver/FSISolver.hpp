/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Christophe Prud'homme <christophe.prudhomme@epfl.ch>
       Date: 2004-11-18

  Copyright (C) 2004 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file FSISolver.hpp
   \author Christophe Prud'homme <christophe.prudhomme@epfl.ch>
   \date 2004-11-18
 */


#ifndef TWODIM

#ifndef __FSISolver_H
#define __FSISolver_H 1

#include <life/lifecore/life.hpp>
#include <life/lifearray/tab.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifealg/newton.hpp>

#include <Epetra_ConfigDefs.h>
#ifdef EPETRA_MPI
	#include <Epetra_MpiComm.h>
#else
	#include <Epetra_SerialComm.h>
#endif

namespace LifeV {

/*!
  \class FSISolver
  \brief solver for Fluid-Structure Interaction

  \c FSISolver uses the FSI operators whose base class is \c FSIOperator to
  solve the FSI problem.  It behaves from an interface point of view
  very much like a NS solver or solid solver.

  The temporal loop is externalized: the member function \c
  iterate(time) must be called at each time steps.  The time data is
  the one associated with the fluid.

  \c FSISolver allows to change the FSI operator using the \c setFSIOperator(name)
  member function that needs the \c name of the new FSI operator to be used.

  \todo Generic fluid and Structure solvers
  \todo Allow delayed initialization

  @author Christophe Prud'homme
  @see FSIOperator
*/
class FSISolver
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FSIOperator::mesh_type									mesh_type;

    typedef FSIOperator::fluid_type::value_type						fluid_value_type;
    typedef FSIOperator::solid_type::value_type						solid_value_type;

    typedef fluid_value_type::Function								fluid_function;
    typedef solid_value_type::Function								solid_function;

    typedef fluid_value_type::source_type							fluid_source_type;
    typedef solid_value_type::source_type							solid_source_type;

    typedef FSIOperator::fluid_bchandler_type						fluid_bchandler_type;
    typedef FSIOperator::solid_bchandler_type						solid_bchandler_type;

    typedef FSIOperator::fluid_bchandler_raw_type					fluid_bchandler_raw_type;
    typedef FSIOperator::solid_bchandler_raw_type					solid_bchandler_raw_type;

    typedef fluid_value_type::data_type								data_fluid;
    typedef solid_value_type::data_type								data_solid;

    typedef FSIOperator::vector_type								vector_type;
    typedef FSIOperator::vector_ptrtype								vector_ptrtype;

    //@}

    /** @name Constructors, destructor
     */
    //@{

    //! default/only constructor for the FSI solver
    /*! \attention the last argument of the constructor gives a
     *  possibility to override the operator that is given by the \c
     *  data_file. An empty string (default) means that it uses the
     *  operator defined in \c data_file
     *  \todo allow to change the FSI operator on the fly
     */
    FSISolver( const std::string& method = "" );

    //! default/only destructor for the FSI solver
    virtual ~FSISolver() {}

    //@}



    /** @name Get Functions
     */
    //@{

    //! get the time step
    Real timeStep() const { return M_oper->dataFluid().dataTime()->getTimeStep(); }

    //! get the final time
    Real timeEnd() const { return M_oper->dataFluid().dataTime()->getEndTime(); }

    //! get the initial time
    Real timeInitial() const { return M_oper->dataFluid().dataTime()->getInitialTime(); }

    //! get the FSI operator
    oper_fsi_ptr_mpi const& FSIOper() const { return M_oper; }
    //oper_fsi_ptr_mpi FSIOper() const { return M_oper; }

    //! get the displacement, which will be on the solidInterfaceMap
    vector_type const& displacement() const
    {
        vector_ptrtype solution;
        M_oper->getSolution(solution);
        return *solution;
    }

    //! get access to the \c bchandler_type for the velocity
//    fluid_bchandler_type& bcHandlerU() { return M_BCh_u; }

    //! get access to the \c bchandler_type for the velocity
//    fluid_bchandler_type const& bcHandlerU() const { return M_BCh_u; }

    //! get access to the \c bchandler_type for the displacement
//    solid_bchandler_type const& bcHandlerD() const { return M_BCh_d; }

    bool isFluid() { return M_oper->isFluid(); }
    bool isSolid() { return M_oper->isSolid(); }

    //@}



    /** @name  Set Functions
     */
    //@{

    //     bool setFluid				(bool fluid) { M_fluid = fluid; M_oper->set }
    //     bool setSolid				(bool solid) { M_solid = solid; }
    void initSol                 ( const vector_type& solInit ) { M_oper->setSolution(solInit); }

    void setSourceTerms          ( const fluid_source_type& fluidSource,
                                   const solid_source_type& solidSource );

    //! set the FSI operator to be used to couple the fluid and the structure
    /*!
     * \param __op FSI operator name
     */
    void setFSIOperator          ( const std::string& __op );

    void setFluidBC              ( const fluid_bchandler_type& bc_fluid );
    void setLinFluidBC           ( const fluid_bchandler_type& bc_dfluid );
    void setInvLinFluidBC        ( const fluid_bchandler_type& bc_dfluid_inv );
    void setHarmonicExtensionBC  ( const fluid_bchandler_type& bc_he );
    void setSolidBC              ( const solid_bchandler_type& bc_solid );
    void setLinSolidBC           ( const solid_bchandler_type& bc_dsolid );
    void setInvLinSolidBC        ( const solid_bchandler_type& bc_dsolid_inv );
    void setFluxBC              (const fluid_bchandler_type& bc_fluid);
    void setRobinBC              (const fluid_bchandler_type& bc_fluid);
//     void setReducedLinFluidBC    ( const fluid_bchandler_type& bc_dredfluid );
//     void setInvReducedLinFluidBC ( const fluid_bchandler_type& bc_dredfluid_inv );

    //@}



    /** @name  Methods
     */
    //@{

    void setDataFromGetPot( const GetPot& dataFile );

    void setup( );

    void initialize( const fluid_function& u0,
                     const fluid_function& p0,
                     const solid_function& d0,
                     const solid_function& w0,
                     const fluid_function& df0=fluid_function() );

    void initialize( const std::string& /*velFName*/,
                     const std::string& /*pressName*/,
                     const std::string& /*velwName*/,
                     const std::string& /*depName*/,
                     const std::string& /*velSName*/,
                     const Real&        /*Tstart = 0.*/);

    virtual void initialize(vector_ptrtype u0=vector_ptrtype(), vector_ptrtype v0=vector_ptrtype());

    void iterate( const Real& time );

    void showMe() {}

    bool isMonolithic(){return M_monolithic;}

private:

    //! forbid default constructor
    FSISolver();

    //! forbid copy constructor
    FSISolver( FSISolver const& );

    //! forbid copy operator
    FSISolver& operator=( FSISolver const& );

    //Private members
    oper_fsi_ptr_mpi							M_oper;

    std::string									M_method;
    bool										M_monolithic;

    bool										M_firstIter;
    UInt										M_maxpf;
    Real										M_defomega;

    Real										M_abstol;
    Real										M_reltol;
    Real										M_etamax;

    Int											M_linesearch;

    boost::shared_ptr<EpetraMap>				M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>				M_solidInterfaceMap;

    boost::shared_ptr<Epetra_MpiComm>			M_epetraComm;
    boost::shared_ptr<Epetra_MpiComm>			M_epetraWorldComm;
    boost::shared_ptr<MPI_Comm>					M_localComm;
    boost::shared_ptr<MPI_Comm>					M_interComm;

    std::ofstream								M_out_iter;
    std::ofstream								M_out_res;

//     data_fluid           M_dataFluid;
//     data_solid           M_dataSolid;

    // be careful here: the BCs must be constructed before the solvers
//     fluid_bchandler_type       M_BCh_u;
//     fluid_bchandler_type       M_BCh_d;
//     fluid_bchandler_type       M_BCh_mesh;

//     fluid_bchandler_type       M_BCh_du;
//     fluid_bchandler_type       M_BCh_dz;


//     FESpace<mesh_type, EpetraMap>* M_uFESpace;
//     FESpace<mesh_type, EpetraMap>* M_pFESpace;
//     FESpace<mesh_type, EpetraMap>* M_dFESpace;
//     FESpace<mesh_type, EpetraMap>* M_mmFESpace;

//     partitionMesh< FSIOperator::mesh_type >*  M_fluidMeshPart;
//     partitionMesh< FSIOperator::mesh_type >*  M_solidMeshPart;

//     FSIOperator::fluid_type        M_fluid;
//     FSIOperator::solid_type        M_solid;
//     FSIOperator::meshmotion_type   M_meshMotion;

//     bool                       M_fluid;
//     bool                       M_solid;

//    Epetra_MpiComm*               M_epetraSolidComm;
};

} // Namespace LifeV
#endif /* __FSISolver_H */
#endif /* TWODIM */
