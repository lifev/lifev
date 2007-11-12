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


#ifndef _STEKLOV_HPP
#define _STEKLOV_HPP

#include <life/lifesolver/FSIOperator.hpp>
#include <life/lifefem/dofInterface3Dto2D.hpp>

namespace LifeV
{

class steklovPoincare : public FSIOperator
{
public:

    typedef FSIOperator super;
    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;
    typedef super::fluid_bchandler_type bchandler_type;

    typedef super::vector_type vector_type;

    // default constructor
    steklovPoincare();

    // destructor
    ~steklovPoincare();

    // Solvers

    void evalResidual(Vector       &_res,
                      const Vector &_disp,
                      const int     _iter);

    void solveJac(Vector       &_muk,
                  const Vector &_res,
                  const double  _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void solveInvLinearFluid();
    void solveInvLinearSolid();


    // Preconditioners

    Vector DDNprecond(Vector const &_z);

    //setters and getters


    Vector dzSolid()       {return M_dzSolid;}
    Vector dzFluid()       {return M_dzFluid;}
    Real   defOmega()      {return M_defOmega;}

    void setResidualS  ( Vector const& _res ){M_residualS = _res;}
    void setResidualF  ( Vector const& _res ){M_residualF = _res;}

    void setResidualFSI( double const* _res );
    void setResidualFSI( Vector const& _res );

    void setDataFromGetPot( GetPot const& data );

    void computeResidualFSI();
    Vector getResidualFSIOnSolid();
    Vector transferFluidOnInterface(Vector const & _vec);
    Vector transferInterfaceOnFluid(Vector const & _vec);
    Vector transferSolidOnInterface(Vector const & _vec);
    Vector transferInterfaceOnSolid(Vector const & _vec);
    void   getResidualFSIOnSolid(Vector& _vec);
    Vector getSolidInterfaceOnSolid(double const* _vec);
    Vector getSolidInterfaceOnFluid(Vector const& _vec);
    Vector getFluidInterfaceOnSolid(Vector const& _vec);

//     void   transferOnLocalInterface  (Vector const &_vec, dofInterface3Dto2D const &_dofInterface);
//     void   transferFromLocalInterface(Vector       &_vec, dofInterface3Dto2D const &_dofInterface);

    Vector & displacement()    {return M_interfaceDisplacement;}
    Vector & residual()        {return M_interfaceStress;}
    Vector & velocity()        {return M_interfaceVelocity;}
    Vector & residualFSI()     {return M_residualFSI;}

    //
    // BCs treatment
    //

    // Solid, Lin. Solid, and inverses
    void setSolidInterfaceDisp         (Vector& disp  , UInt type = 0);
    void setSolidLinInterfaceDisp      (Vector& disp  , UInt type = 0);
    void setSolidInvLinInterfaceStress (Vector& stress, UInt type = 0);
    void setReducedFluidInterfaceAcc   (Vector& acc   , UInt type = 0);
    void setReducedFluidInvInterfaceAcc(Vector& acc   , UInt type = 0);

//    void setSolidInterfaceStress (Vector &stress, UInt type = 0);

    bc_vector_interface bcvSolidInterfaceDisp()       {return M_bcvSolidInterfaceDisp;}
    bc_vector_interface bcvSolidLinInterfaceDisp()    {return M_bcvSolidLinInterfaceDisp;}
//     bc_vector_interface bcvSolidInvInterfaceDisp()    {return M_bcvSolidInvInterfaceDisp;}
    bc_vector_interface bcvSolidInvLinInterfaceStress() {return M_bcvSolidInvLinInterfaceStress;}

    // Fluid, Lin. Fluid, and inverses

    void setFluidInterfaceDisp      (Vector& disp, UInt type = 0);
    void setFluidLinInterfaceVel    (Vector& vel,  UInt type = 0);
//    void setFluidLinInterfaceStress (Vector &disp, UInt type = 0);

    bc_vector_interface bcvFluidInterfaceDisp()       {return M_bcvFluidInterfaceDisp;}
    bc_vector_interface bcvFluidLinInterfaceVel()     {return M_bcvFluidLinInterfaceVel;}
//     bc_vector_interface bcvFluidInvInterfaceDisp()    {return M_bcvFluidInvInterfaceDisp;}
//     bc_vector_interface bcvFluidInvLinInterfaceDisp() {return M_bcvFluidInvLinInterfaceDisp;}

    bc_vector_interface bcvReducedFluidInterfaceAcc()    {return M_bcvReducedFluidInterfaceAcc;}
    bc_vector_interface bcvReducedFluidInvInterfaceAcc() {return M_bcvReducedFluidInvInterfaceAcc;}

    //
    // Misc.
    //

    void setup();

    struct DataJacobian
    {
        DataJacobian():
            M_pFS(0)
            {}
        DataJacobian(steklovPoincare* oper):
            M_pFS(oper)
            {}

        steklovPoincare* M_pFS;
    };


private:

    vector_type*            M_dzSolid;
    vector_type*            M_dzFluid;

    vector_type*            M_rhs_dz;

    vector_type*            M_residualS;
    vector_type*            M_residualF;
    vector_type*            M_residualFSI;
    //! FS strong residual
    vector_type*            M_strongResidualFSI;

    Real                    M_defOmega;
    Real                    M_defOmegaS;
    Real                    M_defOmegaF;

    generalizedAitken<vector_type, Real> M_aitkFS;

    UInt                    M_interfaceNbreDof;
    vector_type*            M_interfaceDisplacement;
    vector_type*            M_interfaceStress;
    vector_type*            M_interfaceVelocity;

    bc_vector_interface     M_bcvSolidInterfaceDisp;
    bc_vector_interface     M_bcvSolidLinInterfaceDisp;
    bc_vector_interface     M_bcvSolidInterfaceStress;
    bc_vector_interface     M_bcvSolidLinInterfaceStress;
    bc_vector_interface     M_bcvSolidInvLinInterfaceStress;

    bc_vector_interface     M_bcvFluidInterfaceDisp;
    bc_vector_interface     M_bcvFluidLinInterfaceDisp;
    bc_vector_interface     M_bcvFluidLinInterfaceVel;
    bc_vector_interface     M_bcvFluidInterfaceStress;
    bc_vector_interface     M_bcvFluidLinInterfaceStress;

    bc_vector_interface     M_bcvReducedFluidInterfaceAcc;
    bc_vector_interface     M_bcvReducedFluidInvInterfaceAcc;

    Real                    M_linearRelTol;


    void eval           (const  vector_type &disp,
                         int    status,
                         vector_type &dispNew,
                         vector_type &veloStruct);

    void  invSfPrime    (const vector_type &res,
                         double       linear_rel_tol,
                         vector_type       &step);

    void  invSsPrime    (const vector_type &res,
                         double       linear_rel_tol,
                         vector_type       &step);


    void  invSfSsPrime  (const vector_type &res,
                         double       linear_rel_tol,
                         vector_type       &step);

    void computeStrongResidualFSI();


    DataJacobian            M_dataJacobian;
};

Real fzeroSP(const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID& i);

void my_matvecSfSsPrime(double *z,
                        double *Jz,
                        AZ_MATRIX* J,
                        int proc_config[]);
}


#define FOR_EACH_DOF_INTERFACE( Expr )                              \
{   \
 for (std::map<ID,ID>::const_iterator it = FSIDofMap.begin(); it != FSIDofMap.end(); ++it) \
     {                                                                                     \
        ID dofF = it->first;                        \
        ID dofS = it->second;                       \
        ID localS = solidDofMap.find(dofS)->second; \
        ID localF = fluidDofMap.find(dofF)->second; \
        for (UInt jDim = 0; jDim < 3; ++jDim)       \
            {                                       \
               (Expr);                              \
            }                                       \
     }                                              \
}


#endif
