/* -*- mode: c++ -*-
   This program is part of the LifeV library

   Author(s): Alexandra Moura <moura@mate.polimi.it>
       Date: 2004-11-05

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
/**
   \file flowRate.hpp
   \author Alexandra Moura <moura@mate.polimi.it>
   \date 2004-11-05
 */

#include "NavierStokesSolverPC.hpp"
#include "zeroDModelSolver.hpp"
#include "NavierStokesWithFlux.hpp"


#ifndef _C3D0D
#define _C3D0D

/**

    Coupler for the 3D (NS - non compliant) and 0D (lumped parameter model) models

*/

namespace LifeV
{

  class C3d0d {

  public:

    typedef RegionMesh3D<LinearTetra> Mesh;
    typedef NavierStokesSolverPC< Mesh > NS;
    typedef NavierStokesHandler<Mesh>::source_type source_type;
    typedef Real ( *Function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& )> source_type;

    //Constructors

    C3d0d(GetPot& data_file, BCHandler& BCh_u);

    //Destructor
    virtual ~C3d0d();

    //Virtual Member Functions

    virtual void initialize(const Function& u0, const Function& p0, Real t0, Real dt) = 0;

    virtual void timeAdvance(source_type const& source, Real const& time) = 0;

    virtual void iterate(Real const& time) = 0;

    //Member Functions   
 
    void startWrtMatlab();
    
    void finishWrtMatlab();

    //Mutators 

    //Setters
    
    UInt strategy();

    void setStrategy(int const&);

    boost::shared_ptr<NS> & ns() {return _M_ns;}

    zeroDModelSolver & network() {return _M_network;}

    NavierStokesWithFlux< NS > & nsFlux(){return _M_nsFlux;}

  protected: 

    Real _M_dt;
    Real _M_T;

    BCHandler                  _M_BCh_u;  
    boost::shared_ptr<NS>      _M_ns;
    zeroDModelSolver           _M_network;
    NavierStokesWithFlux< NS > _M_nsFlux;

    std::ofstream _M_outfile;

  private: 
    
    //Coupling Strategy
    UInt _M_strategy;
  };
    
//************************************************************************************
//                                   IMPLEMENTATION                                   
//************************************************************************************


//Constructors
C3d0d::C3d0d(GetPot &data_file, BCHandler &BCh_u)
  :
    _M_dt( data_file("fluid/discretization/timestep",0.001) ),
    _M_T( data_file("fluid/physics/endtime",0.002) ),
    _M_ns( new NS (data_file, feTetraP1bubble, feTetraP1,
                   quadRuleTetra15pt,quadRuleTria3pt,
		   quadRuleTetra5pt, quadRuleTria3pt, BCh_u) ),
    _M_network(data_file),
    _M_nsFlux(_M_ns),
    _M_strategy(data_file("problem/strategy",0))   
{}

//
void 
C3d0d::startWrtMatlab()
{

  _M_outfile.open("res_Q.m", std::ios::app);
  _M_outfile << "xx=[0:" << _M_dt << ":" << _M_T << "]; " << std::endl;
  _M_outfile << "Q = [ " << std::endl;
  _M_outfile.close();

  _M_outfile.open("res_DP.m", std::ios::app);
  _M_outfile << "xx=[0:" << _M_dt << ":" << _M_T << "]; " << std::endl;
  _M_outfile << "DP = [ " << std::endl;
  _M_outfile.close();

}

//
void 
C3d0d::finishWrtMatlab()
{
    _M_outfile.open("res_Q.m", std::ios::app);
    _M_outfile << "    ]; " << std::endl;
    _M_outfile << "figure;" << std::endl;
    _M_outfile << "plot(xx,Q);" << std::endl;
    _M_outfile.close();

    _M_outfile.open("res_DP.m", std::ios::app);
    _M_outfile << "     ]; " << std::endl;
    _M_outfile << "figure;" << std::endl;
    _M_outfile << "plot(xx,DP);" << std::endl;
    _M_outfile.close();

}

//Gets coupling strategy
UInt
C3d0d::strategy()
{
  return _M_strategy;
}

//Sets coupling strategy
void 
C3d0d::setStrategy(int const& strategy)
{
    _M_strategy = strategy;
}

// Destructor
C3d0d::~C3d0d()
{
}

} // namespace LifeV

#endif
