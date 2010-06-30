//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
/**
   \file lumpedHeart.hpp
   \author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
   \date 2009-06-03
*/
#ifndef __LUMPEDHEART_HPP
#define __LUMPEDHEART_HPP

#include <life/lifecore/life.hpp>
#include <life/lifecore/GetPot.hpp>
#include <life/lifesolver/FSISolver.hpp>
#include <life/lifesolver/FSIOperator.hpp>
#include <lifemc/lifesolver/BCInterface.hpp>

namespace LifeV
{
class LumpedHeart
{
public:
    typedef BCInterface< FSIOperator >                                                     bc_type;
    typedef FSIOperator::vector_type                                                       vector_type;
    typedef FSIOperator::vector_ptrtype                                                    vector_ptrtype;

    LumpedHeart()
        :
        M_time(0.),
        M_BC(),
        M_ODEscheme(1),
        M_dt(0.),
        M_T_max()  ,
        M_E_max()  ,
        M_V_0()    ,
        M_RV_art() ,
        M_RA_V()   ,
        M_LV_art() ,
        M_LA_V()   ,
        M_PV(),
        M_intFlux(),
        M_Vt_ao()
    {}

    void initParameters      ( FSIOperator&  oper,
                               const std::string&    FileName);

    static Real& outPressure         (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    void renewParameters     ( FSIOperator&  oper , const int& flag, const Real& time, const Real& flux);

    Real fZero               (const Real& t, const Real& x, const Real& y, const Real& z, const ID& i);

    static Real                               M_pressure;

private:

    //! @name Private Methods
    //@{

    Real                                      M_elastance(const Real& t);

    //! Short description of this method
    /*!
        Add more details about the method.
        NOTE: short description is automatically added before this part.
     */


    //@}

    Real                                      M_time;
    boost::shared_ptr< bc_type >              M_BC;
    BdfT<Real>                                M_ODEscheme;
    Real                                      M_dt;
    Real                                      M_T_max  ;
    Real                                      M_E_max  ;
    Real                                      M_V_0    ; //reference ventricular volume
    Real                                      M_RV_art ;
    Real                                      M_RA_V   ;
    Real                                      M_LV_art ;
    Real                                      M_LA_V   ;
    Real                                      M_PV     ; //ventricular pressure
    Real                                      M_intFlux;
    Real                                      M_Vt_ao;
    Real                                      M_Tpb  ;
    Real                                      M_Tpw  ;
};
}

#endif /* __LUMPEDHEART_HPP */

