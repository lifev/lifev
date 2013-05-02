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
 \include ../../doc/api/bibliography/newton
 \include ../../doc/api/bibliography/fluidstructure

    @file
    @brief solver class for FSI

    @author Christophe Prud'homme
    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gilles Fourestey <fourestey@cscs.ch>
    @contributor Paolo Crosetto <paolo.crosetto@epfl.ch>
    @maintainer  Paolo Crosetto <paolo.crosetto@epfl.ch>

    @date 10-12-2010

    This class handles the FSI iterations, it is generalized in the context of geometric multiscale applications by the
    class MS_ModelFSI, thus it will become obsolete.
 */



#ifndef __FSISolver_H
#define __FSISolver_H 1

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

#include <lifev/core/LifeV.hpp>

//#include <life/lifealg/newton.hpp>

#include <lifev/fsi/solver/FSIOperator.hpp>

namespace LifeV
{

/*!
  \class FSISolver
  \brief solver for Fluid-Structure Interaction

  This class handles the FSI iterations, it is generalized in the context of geometric multiscale applications by the
  class MS_ModelFSI, thus it will become obsolete.

  \c FSISolver uses the FSI operators whose base class is \c FSI to
  solve the FSI problem.  It behaves from an interface point of view
  very much like a NS solver or solid solver.

  The temporal loop is externalized: the member function \c
  iterate(time) must be called at each time steps.  The time data is
  the one associated with the fluid.

  \c FSISolver allows to change the FSI operator using the \c setFSI(name)
  member function that needs the \c name of the new FSI operator to be used.

  \todo Generic fluid and Structure solvers
  \todo Allow delayed initialization

  @see FSI
*/
class FSISolver
{
public:

    /** @name Typedefs
     */
    //@{

    typedef FSIOperator                                                FSIOper_Type;
    typedef boost::shared_ptr<FSIOper_Type>                            FSIOperPtr_Type;

    typedef FSIOperator::mesh_Type                                     mesh_Type;

    typedef FSIOperator::fluidPtr_Type::element_type                     fluid_Type;
    typedef FSIOperator::solidPtr_Type::element_type                     solid_Type;

    typedef fluid_Type::function_Type                                  fluidFunction_Type;
    typedef solid_Type::function                                       solidFunction_Type;

    typedef fluid_Type::source_Type                                    fluidSource_Type;
    typedef solid_Type::source_Type                                    solidSource_Type;

    typedef FSIOperator::fluidBchandlerPtr_Type                        fluidBchandlerPtr_Type;
    typedef FSIOperator::solidBchandlerPtr_Type                        solidBchandlerPtr_Type;

    typedef FSIOperator::fluidBchandler_Type                           fluidBchandler_Type;
    typedef FSIOperator::solidBchandler_Type                           solidBchandler_Type;

    typedef fluid_Type::data_Type                                      fluidData_Type;
    typedef solid_Type::data_Type                                      solidData_Type;

    typedef FSIOperator::dataPtr_Type                                  dataPtr_Type;

    typedef FSIOperator::vector_Type                                   vector_Type;
    typedef FSIOperator::vectorPtr_Type                                vectorPtr_Type;

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
    FSISolver( );

    //! default/only destructor for the FSI solver
    virtual ~FSISolver() {}

    //@}



    /** @name Get Functions
     */
    //@{

    //! get the FSI operator
    FSIOperPtr_Type const& FSIOper() const
    {
        return M_oper;
    }
    //FSIOperPtr_Type FSIOper() const { return M_oper; }

    //! get the displacement, which will be on the solidInterfaceMap
    LIFEV_DEPRECATED ( const vector_Type& displacement() const )
    {
        if ( M_epetraWorldComm->MyPID() == 0 )
        {
            std::cerr << std::endl << "Warning: FSISolver::displacement() is deprecated!" << std::endl
                      << "         You should not access the solution inside FSISolver or FSIOperator!" << std::endl;
        }

        return M_oper->solution();
    }

    //! get access to the \c bchandler_type for the velocity
    //    fluidBchandlerPtr_Type& bcHandlerU() { return M_BCh_u; }

    //! get access to the \c bchandler_type for the velocity
    //    fluidBchandlerPtr_Type const& bcHandlerU() const { return M_BCh_u; }

    //! get access to the \c bchandler_type for the displacement
    //    solidBchandlerPtr_Type const& bcHandlerD() const { return M_BCh_d; }

    bool isFluid()
    {
        return M_oper->isFluid();
    }
    bool isSolid()
    {
        return M_oper->isSolid();
    }

    //@}



    /** @name  Set Functions
     */
    //@{

    void setSourceTerms          ( const fluidSource_Type& fluidSource,
                                   const solidSource_Type& solidSource );

    //! set the FSI operator to be used to couple the fluid and the structure
    /*!
     * \param __op FSI operator name
     */
    void setFSI          ( );

    void setFluidBC              ( const fluidBchandlerPtr_Type& bc_fluid );
    void setLinFluidBC           ( const fluidBchandlerPtr_Type& bc_dfluid );
    void setInvLinFluidBC        ( const fluidBchandlerPtr_Type& bc_dfluid_inv );
    void setHarmonicExtensionBC  ( const fluidBchandlerPtr_Type& bc_he );
    void setSolidBC              ( const solidBchandlerPtr_Type& bc_solid );
    void setLinSolidBC           ( const solidBchandlerPtr_Type& bc_dsolid );
    void setInvLinSolidBC        ( const solidBchandlerPtr_Type& bc_dsolid_inv );
    //!\todo{kill this method}
    void setFluxBC              (const fluidBchandlerPtr_Type& bc_fluid);
    //!\todo{kill this method}
    void setRobinBC              (const fluidBchandlerPtr_Type& bc_fluid);

    //@}



    /** @name  Methods
     */
    //@{

    void setData ( const dataPtr_Type& data );

    void setup( );


    virtual void initialize (std::vector<vectorPtr_Type> u0 = std::vector<vectorPtr_Type> (0), std::vector<vectorPtr_Type> ds0 = std::vector<vectorPtr_Type> (0), std::vector<vectorPtr_Type> df0 = std::vector<vectorPtr_Type> (0) );

    LIFEV_DEPRECATED ( void iterate() );

    void iterate ( vectorPtr_Type& solution );

    void showMe() {}

private:

    //! forbid copy constructor
    FSISolver ( FSISolver const& );

    //! forbid copy operator
    FSISolver& operator= ( FSISolver const& );

    //Private members
    FSIOperPtr_Type                             M_oper;

    dataPtr_Type                                M_data;

    boost::shared_ptr<MapEpetra>                M_fluidInterfaceMap;
    boost::shared_ptr<MapEpetra>                M_solidInterfaceMap;

    boost::shared_ptr<Epetra_MpiComm>            M_epetraComm;
    boost::shared_ptr<Epetra_MpiComm>            M_epetraWorldComm;
    boost::shared_ptr<MPI_Comm>                    M_localComm;
    boost::shared_ptr<MPI_Comm>                    M_interComm;

    std::ofstream                                M_out_iter;
    std::ofstream                                M_out_res;

    //     data_fluid           M_dataFluid;
    //     data_solid           M_dataSolid;

    // be careful here: the BCs must be constructed before the solvers
    //     fluidBchandlerPtr_Type       M_BCh_u;
    //     fluidBchandlerPtr_Type       M_BCh_d;
    //     fluidBchandlerPtr_Type       M_BCh_mesh;

    //     fluidBchandlerPtr_Type       M_BCh_du;
    //     fluidBchandlerPtr_Type       M_BCh_dz;


    //     FESpace<mesh_Type, MapEpetra>* M_uFESpace;
    //     FESpace<mesh_Type, MapEpetra>* M_pFESpace;
    //     FESpace<mesh_Type, MapEpetra>* M_dFESpace;
    //     FESpace<mesh_Type, MapEpetra>* M_mmFESpace;

    //     partitionMesh< FSI::mesh_Type >*  M_fluidMeshPart;
    //     partitionMesh< FSI::mesh_Type >*  M_solidMeshPart;

    //     FSI::fluid_type        M_fluid;
    //     FSI::solid_type        M_solid;
    //     FSI::meshmotion_type   M_meshMotion;

    //     bool                       M_fluid;
    //     bool                       M_solid;

    //    Epetra_MpiComm*               M_epetraSolidComm;
};

} // Namespace LifeV
#endif /* __FSISolver_H */
