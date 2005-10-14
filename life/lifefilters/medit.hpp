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

/*!
  \file medit.h
  \author M.A. Fernandez
  \date 4/2005
  \version 1.0

  \brief This file contains a Navier-Stokes solver class which implements a
  stabilized implicit scheme with coupled solving.
*/
#ifndef _MEDIT_H_
#define _NEDIT_H_

#include <life/lifecore/life.hpp>
#include <life/lifefilters/medit_wrtrs.hpp>
#include <life/lifefilters/gmv_wrtrs.hpp>
#include <vector>

namespace LifeV
{

  /*!
    \class Medit
    
  

  */
  class Medit {

  public:
    
    template<typename Solver>  Medit(Solver& solver);
    template<typename Solver>  void postProcessALE(Solver& solver);
    template<typename Solver>  void postProcess(Solver& solver);

  private:
    UInt _count;
    UInt _verbose;
    Vector _uMagnitude;
    Real _factor;
  };


  //
  // Implementation
  //



  template<typename Solver> Medit::Medit(Solver& solver):
    _count(0),
    _verbose(solver.verbose()),
    _uMagnitude(solver.uDof().numTotalDof()),
    _factor(solver.factor()) {}


  template <typename Solver> void Medit::postProcessALE(Solver& solver)
  {

    std::ostringstream index;
    std::string name;

    const Vector& u=solver.u();
    const Vector& p=solver.p();
    const Vector& d=solver.d();
     
    UInt dim = solver.uDof().numTotalDof();

    ++_count;
    if ( fmod( float( _count ), float( _verbose ) ) == 0.0 )
    {
        std::cout << "  Medit-  Post-processing \n";
        index << ( _count / _verbose );

        switch ( index.str().size() )
        {
        case 1:
            name = "00" + index.str();
            break;
        case 2:
            name = "0" + index.str();
            break;
        case 3:
            name = index.str();
            break;
        }
	
	// vector magnitude 
	std::vector<Real> unode(3);
	Real sum=0;
	
	for (UInt i=0; i < dim ; ++i) {
	  for (UInt ic = 0; ic < nDimensions; ++ic)
	    unode[ic]= u[ ic* dim + i];
	  for (UInt ic = 0; ic < nDimensions; ++ic)
	    sum += unode[ic]*unode[ic];
	  _uMagnitude[i] = sqrt(sum);
	  sum = 0;
	}  

	// postprocess data file for Medit
        wr_medit_ascii_scalar( "press." + name + ".bb", p.data().begin(), dim );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.data().begin(), dim );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.data().begin() + dim, dim );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.data().begin() + 2 * dim,  dim );
	wr_medit_ascii_scalar( "veloM." + name + ".bb", _uMagnitude.data().begin(), dim );

        wr_medit_ascii( "press." + name + ".mesh", solver.mesh(), d, _factor );
       	system( ("ln -s press." + name + ".mesh veloM." + name + ".mesh" ).data() );
        system( ("ln -s press." + name + ".mesh vel_x." + name + ".mesh" ).data() );
        system( ("ln -s press." + name + ".mesh vel_y." + name + ".mesh" ).data() );
        system( ("ln -s press." + name + ".mesh vel_z." + name + ".mesh" ).data() );
  
	// postprocess data file for GMV
        wr_gmv_ascii( "test." + name + ".inp", solver.mesh(), dim,  u.data().begin(),  p.data().begin() , d, _factor);
    }
  }
  


  template <typename Solver> void Medit::postProcess(Solver& solver)
  {

    std::ostringstream index;
    std::string name;

    Vector& u=solver.u();
    Vector& p=solver.p();
     
    UInt dim = solver.uDof().numTotalDof();

    ++_count;
    if ( fmod( float( _count ), float( _verbose ) ) == 0.0 )
    {
        std::cout << "  Medit-  Post-processing \n";
        index << ( _count / _verbose );

        switch ( index.str().size() )
        {
        case 1:
            name = "00" + index.str();
            break;
        case 2:
            name = "0" + index.str();
            break;
        case 3:
            name = index.str();
            break;
        }
	
	// vector magnitude 
	std::vector<Real> unode(3);
	Real sum=0;
	
	for (UInt i=0; i < dim ; ++i) {
	  for (UInt ic = 0; ic < nDimensions; ++ic)
	    unode[ic]= u[ ic* dim + i];
	  for (UInt ic = 0; ic < nDimensions; ++ic)
	    sum += unode[ic]*unode[ic];
	  _uMagnitude[i] = sqrt(sum);
	  sum = 0;
	}  

	// postprocess data file for Medit
        wr_medit_ascii_scalar( "press." + name + ".bb", p.data().begin(), dim );
        wr_medit_ascii_scalar( "vel_x." + name + ".bb", u.data().begin(), dim );
        wr_medit_ascii_scalar( "vel_y." + name + ".bb", u.data().begin() + dim, dim );
        wr_medit_ascii_scalar( "vel_z." + name + ".bb", u.data().begin() + 2 * dim,  dim );
	wr_medit_ascii_scalar( "veloM." + name + ".bb", _uMagnitude.data().begin(), dim );
	system( ("ln -s " + solver.meshDir()+solver.meshFile() + " press." + name + ".mesh" ).data() );
       	system( ("ln -s " + solver.meshDir()+solver.meshFile() + " veloM." + name + ".mesh" ).data() );
        system( ("ln -s " + solver.meshDir()+solver.meshFile() + " vel_x." + name + ".mesh" ).data() );
        system( ("ln -s " + solver.meshDir()+solver.meshFile() + " vel_y." + name + ".mesh" ).data() );
        system( ("ln -s " + solver.meshDir()+solver.meshFile() + " vel_z." + name + ".mesh" ).data() );
  
	// postprocess data file for GMV
        wr_gmv_ascii( "test." + name + ".inp", solver.mesh(), dim,  u.data().begin(),  p.data().begin());
    }
  }
  


}


#endif
