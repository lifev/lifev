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


#ifndef _FP_HPP
#define _FP_HPP

#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifefilters/ensight.hpp>

namespace LifeV
{

class fixedPoint : public FSIOperator
{
public:

    typedef FSIOperator super;
    typedef super::fluid_type           fluid_type;
    typedef super::solid_type           solid_type;
    typedef super::fluid_bchandler_type bchandler_type;

    typedef fluid_raw_type::vector_type  vector_type;
    // constructors

    fixedPoint();

    // destructor

    ~fixedPoint();

    // member functions

    void evalResidual(vector_type&        _res,
                      const vector_type&  _disp,
                      const int           _iter);

    void solveJac     (vector_type&       _muk,
                       const vector_type& _res,
                       const double       _linearRelTol);

    void setUpBC     ();

    Real   defOmega() {return M_defOmega;}

    void setDataFromGetPot( GetPot const& data );

    void setup();

    //vector_type& displacement()    {return *M_displacement;}
    //vector_type& residual()        {return *M_stress;}
    //vector_type& velocity()        {return *M_velocity;}
    //vector_type& residualFSI()     {return *M_residualFSI;}

private:

    Real                                 M_defOmega;

    // If M_updateEvery == 1, normal fixedPoint algorithm
    // If M_updateEvery  > 1, recompute computational domain every M_updateEvery iterations (transpiration)
    // If M_updateEvery <= 0, recompute computational domain and matrices only at first subiteration (semi-implicit)
    int                                  M_updateEvery;
    generalizedAitken<vector_type, Real> M_aitkFS;

//     boost::shared_ptr<vector_type>       M_displacement;
//     boost::shared_ptr<vector_type>       M_stress;
//     boost::shared_ptr<vector_type>       M_velocity;
//     boost::shared_ptr<vector_type>       M_residualFSI;

    boost::shared_ptr<vector_type>       M_rhsNew;
    boost::shared_ptr<vector_type>       M_beta;

    boost::shared_ptr< Ensight< RegionMesh3D<LinearTetra> > > M_ensightFluid;
    boost::shared_ptr<vector_type>                            M_velAndPressure;

    void eval( const vector_type& disp, int status );

//    FSIOperator* createFP(){ return new fixedPoint(); }
//    static bool              reg;

};


//FSIOperator* createFP(){ return new fixedPoint();}

}

#endif
