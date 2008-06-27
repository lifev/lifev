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
#ifndef __FSISolver_H
#define __FSISolver_H 1

#include <life/lifearray/tab.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifealg/nonLinRichardson.hpp>
#include <life/lifealg/newton.hpp>
#include <life/lifesolver/fixedPointBase.hpp>
#include <life/lifesolver/exactJacobianBase.hpp>


#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace LifeV
{
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

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> fct_type;

    typedef FSIOperator::mesh_type mesh_type;

    typedef FSIOperator::fluid_type::value_type::source_type fluid_source_type;
    typedef FSIOperator::solid_type::value_type::source_type solid_source_type;

    typedef FSIOperator::fluid_bchandler_type fluid_bchandler_type;
    typedef FSIOperator::solid_bchandler_type solid_bchandler_type;

    typedef FSIOperator::fluid_bchandler_raw_type fluid_bchandler_raw_type;
    typedef FSIOperator::solid_bchandler_raw_type solid_bchandler_raw_type;

    typedef FSIOperator::fluid_type::value_type::data_type    data_fluid;
    typedef FSIOperator::solid_type::value_type::data_type    data_solid;

    typedef FSIOperator::vector_type                          vector_type;
    //@}

    /** @name Constructors, destructor
     */
    //@{

    /*!
      \brief default/only constructor for the FSI solver

      \attention the last argument of the constructor gives a
      possibility to override the operator that is given by the \c
      data_file. An empty string (default) means that it uses the
      operator defined in \c data_file
      \todo allow to change the FSI operator on the fly
     */

    FSISolver( GetPot const& datafile,
               std::string    __oper = "" );
    //@}

    /** @name Operator overloads
     */
    //@{

    ~FSISolver()
        {
        }
    //@}


    //@}

    /** @name Accessors
     */
    //@{

    //! get the time step
    Real timeStep() const { return M_oper->dataFluid().timestep(); }

    //! get the final time
    Real timeEnd() const { return M_oper->dataFluid().endtime(); }

    //! get the FSI operator
    oper_fsi_ptr_mpi const& operFSI() const { return M_oper; }

    //! get the displacement, which will be on the solidInterfaceMap
    vector_type const& displacement() const { return *M_lambda; }

    //! get access to the \c bchandler_type for the velocity
//    fluid_bchandler_type& bcHandlerU() { return M_BCh_u; }

    //! get access to the \c bchandler_type for the velocity
//    fluid_bchandler_type const& bcHandlerU() const { return M_BCh_u; }

    //! get access to the \c bchandler_type for the displacement
//    solid_bchandler_type const& bcHandlerD() const { return M_BCh_d; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setSourceTerms( fluid_source_type const& __f, solid_source_type const& __s )
        {
            M_oper->fluid().setSourceTerm( __f );
            M_oper->solid().setSourceTerm( __s );
        }

    /*!
      \brief set the FSI operator to be used to couple the fluid and the structure

      The name of the operator has to be passed as an argument

      \param __op FSI operator name
    */
    void setFSIOperator( std::string const& __op );
    oper_fsi_ptr_mpi FSIOper(){ return M_oper; }
    //@}

    void setFluidBC             (fluid_bchandler_type bc_fluid);
    void setLinFluidBC          (fluid_bchandler_type bc_dfluid);
    void setInvLinFluidBC       (fluid_bchandler_type bc_dfluid_inv);
    void setHarmonicExtensionBC (fluid_bchandler_type bc_he);
    void setSolidBC             (solid_bchandler_type bc_solid);
    void setLinSolidBC          (solid_bchandler_type bc_dsolid);
    void setInvLinSolidBC       (solid_bchandler_type bc_dsolid_inv);
//     void setReducedLinFluidBC   (fluid_bchandler_type bc_dredfluid);
//     void setInvReducedLinFluidBC(fluid_bchandler_type bc_dredfluid_inv);

    /** @name  Methods
     */
    //@{

    void initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                     FSIOperator::solid_type::value_type::Function const& p0,
                     FSIOperator::solid_type::value_type::Function const& d0,
                     FSIOperator::solid_type::value_type::Function const& w0 )
        {
            Debug( 6220 ) << "FSISover:: solid init \n";
            if (this->isSolid())
                M_oper->solid().initialize(d0, w0);
            Debug( 6220 ) << "FSISover:: fluid init \n";
            if (this->isFluid())
                M_oper->fluid().initialize(u0, p0);
        }

    void initialize( std::string /*velFName*/,
                     std::string /*pressName*/,
                     std::string /*velwName*/,
                     std::string /*depName*/,
                     std::string /*velSName*/,
                     double      /*Tstart = 0.*/)
        {
//             M_oper->fluid().initialize(velFName, pressName, velwName, Tstart);
//             M_oper->solid().initialize(depName, velSName, Tstart);
        }

    void iterate( Real time );

    void showMe()
        {
//             M_fluid->showMe();
//             M_solid->showMe();
        }

    //@}

    bool isFluid(){return M_oper->isFluid();}
    bool isSolid(){return M_oper->isSolid();}
//     bool setFluid(bool fluid){M_fluid = fluid; M_oper->set}
//     bool setSolid(bool solid){M_solid = solid;}

private:

    //! forbid default constructor
    FSISolver();

    //! forbid copy constructor
    FSISolver( FSISolver const& );

    //! forbid copy operator
    FSISolver& operator=( FSISolver const& );

//     data_fluid           M_dataFluid;
//     data_solid           M_dataSolid;

    // be careful here: the BCs must be constructed before the solvers
//     fluid_bchandler_type       M_BCh_u;
//     fluid_bchandler_type       M_BCh_d;
//     fluid_bchandler_type       M_BCh_mesh;

//     fluid_bchandler_type       M_BCh_du;
//     fluid_bchandler_type       M_BCh_dz;

    oper_fsi_ptr_mpi           M_oper;


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

    boost::shared_ptr<vector_type>               M_lambda;
    boost::shared_ptr<vector_type>               M_lambdaDot;

    /* data */
    bool                       M_firstIter;
    std::string                M_method;
    UInt                       M_maxpf;
    Real                       M_defomega;

    Real                       M_abstol;
    Real                       M_reltol;
    Real                       M_etamax;

    int                        M_linesearch;

    boost::shared_ptr<EpetraMap>      M_fluidInterfaceMap;
    boost::shared_ptr<EpetraMap>      M_solidInterfaceMap;

    boost::shared_ptr<Epetra_MpiComm> M_epetraComm;
    boost::shared_ptr<Epetra_MpiComm> M_epetraWorldComm;
    boost::shared_ptr<MPI_Comm>       M_localComm;
    boost::shared_ptr<MPI_Comm>       M_interComm;

//    Epetra_MpiComm*               M_epetraSolidComm;
    /* streams */
    std::ofstream out_iter;
    std::ofstream out_res;

};

}
#endif /* __FSISolver_H */
