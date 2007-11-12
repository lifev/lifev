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


#ifndef _EJ_HPP
#define _EJ_HPP

#include <life/lifesolver/FSIOperator.hpp>

namespace LifeV
{


class Epetra_ExactJacobian:
    Epetra_Operator
{

public:

    ~Epetra_ExactJacobian ();

    int 	SetUseTranspose (bool UseTranspose);
    int 	Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    int 	ApplyInverse (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    double 	NormInf () const;
    const char * 	Label () const;
    bool 	UseTranspose () const;
    bool 	HasNormInf () const;
    const Epetra_Comm & 	Comm () const;
    const Epetra_Map & 	OperatorDomainMap () const;
    const Epetra_Map & 	OperatorRangeMap () const;

private:

};



class exactJacobian : public FSIOperator
{
public:

    typedef FSIOperator super;
    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;

    typedef super::fluid_bchandler_type  bchandler_type;

    typedef fluid_raw_type::vector_type  vector_type;

    // default constructor
    exactJacobian():
        super(),
        M_dz     (),
        M_rhs_dz (),
        M_dataJacobian(this)
        {}

    // destructor
    ~exactJacobian();

    // member functions

    void evalResidual(vector_type       &_res,
                      const vector_type &_disp,
                      const int     _iter);
    void solveJac(vector_type       &_muk,
                  const vector_type &_res,
                  const double  _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void setup();

//     vector_type dz() {return *M_dz;}

//     vector_type & displacement()    {return this->solid().disp();}
//     vector_type & residual()        {return M_interfaceStress;}
//     vector_type & velocity()        {return this->solid().w();}
//     vector_type & residualFSI()     {return M_residualFSI;}

    //
    // BCs treatment
    //

    // Solid, Lin. Solid, and inverses
    void setSolidInterfaceDisp         (vector_type& disp  , UInt type = 0);
    void setSolidLinInterfaceDisp      (vector_type& disp  , UInt type = 0);
    void setSolidInvLinInterfaceStress (vector_type& stress, UInt type = 0);
    void setReducedFluidInterfaceAcc   (vector_type& acc   , UInt type = 0);
    void setReducedFluidInvInterfaceAcc(vector_type& acc   , UInt type = 0);

//     void setfluidLoadToStructure(vector_type &acc,
//                                  UInt type);
//    void setSolidInterfaceStress (vector_type &stress, UInt type = 0);

    bc_vector_interface bcvFluidLoadToStructure()       {return M_bcvFluidLoadToStructure;}
//     bc_vector_interface bcvSolidLinInterfaceDisp()    {return M_bcvSolidLinInterfaceDisp;}
// //     bc_vector_interface bcvSolidInvInterfaceDisp()    {return M_bcvSolidInvInterfaceDisp;}
//     bc_vector_interface bcvSolidInvLinInterfaceStress() {return M_bcvSolidInvLinInterfaceStress;}

    // Fluid, Lin. Fluid, and inverses

    void setFluidInterfaceDisp      (vector_type &disp, UInt type = 0);
//    void setFluidLinInterfaceStress (vector_type &disp, UInt type = 0);

//     bc_vector_interface bcvFluidLoadToSolid()       {return M_bcvFluidLoadToSolid;}
//     bc_vector_interface bcvFluidLinInterfaceDisp()    {return M_bcvFluidLinInterfaceDisp;}
// //     bc_vector_interface bcvFluidInvInterfaceDisp()    {return M_bcvFluidInvInterfaceDisp;}
// //     bc_vector_interface bcvFluidInvLinInterfaceDisp() {return M_bcvFluidInvLinInterfaceDisp;}

//     bc_vector_interface bcvReducedFluidInterfaceAcc()    {return M_bcvReducedFluidInterfaceAcc;}
//     bc_vector_interface bcvReducedFluidInvInterfaceAcc() {return M_bcvReducedFluidInvInterfaceAcc;}

    struct dataJacobian
    {

        dataJacobian()
            :
            M_pFS( 0 )
            {}

        dataJacobian(exactJacobian* oper)
            :
            M_pFS(oper){}

        exactJacobian* M_pFS;

    };

private:

     boost::shared_ptr<vector_type>      M_dz;
     boost::shared_ptr<vector_type>      M_rhs_dz;

     boost::shared_ptr<vector_type>      M_residualFSI;

     boost::shared_ptr<vector_type>      M_interfaceDisplacement;
     boost::shared_ptr<vector_type>      M_interfaceStress;

    boost::shared_ptr<vector_type>       M_rhsNew;
    boost::shared_ptr<vector_type>       M_beta;

//     bc_vector_interface     M_bcvFluidLoadToStructure;

    void eval(vector_type& dispNew, vector_type& velo, const vector_type& disp, int status);

    dataJacobian                         M_dataJacobian;
};


void my_matvecJacobianEJ(double *z, double *Jz, AZ_MATRIX* J, int proc_config[]);

Real fzeroEJ(const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID& i);


}

#endif
