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

#include <tab.hpp>
#include <operFS.hpp>
#include <nonLinRichardson.hpp>
#include <newton.hpp>

namespace LifeV
{
/*!
  \class FSISolver
  \brief solver for Fluid-Structure Interaction

  \c FSISolver uses the FSI operators whose base class is \c operFS to
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
  @see operFS
*/
class FSISolver
{
public:

    /** @name Typedefs
     */
    //@{
    typedef operFS::fluid_type::value_type::source_type fluid_source_type;
    typedef operFS::fluid_type::value_type::source_type solid_source_type;
    typedef operFS::bchandler_type bchandler_type;
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
               bchandler_type& __bcu, bchandler_type& __bcd, bchandler_type& __bchext,
               std::string __oper = "" );

    //@}

    /** @name Operator overloads
     */
    //@{


    //@}


    //@}

    /** @name Accessors
     */
    //@{

    //! get the time step
    Real timeStep() const { return M_fluid->timestep(); }

    //! get the final time
    Real timeEnd() const { return M_fluid->endtime(); }

    //! get the FSI operator
    oper_fsi_ptr const& FSIOperator() const { return M_oper; }

    //! get the displacement
    Vector const& displacement() const { return M_disp; }

    //! get access to the \c bchandler_type for the velocity
    bchandler_type& bcHandlerU() { return M_BCh_u; }

    //! get access to the \c bchandler_type for the velocity
    bchandler_type const& bcHandlerU() const { return M_BCh_u; }

    //! get access to the \c bchandler_type for the displacement
    bchandler_type const& bcHandlerD() const { return M_BCh_d; }

    //@}

    /** @name  Mutators
     */
    //@{

    void setSourceTerms( fluid_source_type const& __f, solid_source_type const& __s )
        {
            M_fluid->setSourceTerm( __f );
            M_solid->setSourceTerm( __s );
        }

    /*!
      \brief set the FSI operator to be used to couple the fluid and the structure

      The name of the operator has to be passed as an argument

      \param __op FSI operator name
    */
    void setFSIOperator( std::string const& __op );

    //@}


    /** @name  Methods
     */
    //@{

    void initialize( operFS::fluid_type::value_type::Function const& __u0,
                     operFS::solid_type::value_type::Function const& __d0,
                     operFS::solid_type::value_type::Function const& __w0 )
        {
            M_oper->fluid().initialize(__u0);
            M_oper->solid().initialize(__d0,__w0);
        }
    void iterate( Real time );

    void showMe()
        {
            M_fluid->showMe();
            M_solid->showMe();
        }

    //@}

private:

    //! forbid default constructor
    FSISolver();

    //! forbid copy constructor
    FSISolver( FSISolver const& );

    //! forbid copy operator
    FSISolver& operator=( FSISolver const& );

private:

    // be careful here: the BCs must be constructed before the solvers
    bchandler_type M_BCh_u;
    bchandler_type M_BCh_d;
    bchandler_type M_BCh_mesh;

    oper_fsi_ptr M_oper;

    operFS::fluid_type M_fluid;
    operFS::solid_type M_solid;



    Vector    M_disp;
    Vector    M_velo;

    /* data */
    std::string M_method;
    UInt M_maxpf;
    Real M_defomega;

    Real M_abstol;
    Real M_reltol;
    Real M_etamax;

    int M_linesearch;

    /* streams */
    std::ofstream out_iter;
    std::ofstream out_res;

};

}
#endif /* __FSISolver_H */
