/* -*- mode: c++ -*-
   This program is part of the LifeV library
   Copyright (C) 2001,2002,2003,2004 EPFL, INRIA, Politecnico di Milano

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


#ifndef TWODIM
#ifndef _EJ_HPP
#define _EJ_HPP

#include <life/lifesolver/FSIOperator.hpp>

namespace LifeV
{


class Epetra_ExactJacobian;
//class Epetra_ExactJacobian;




class exactJacobian : public FSIOperator
{
public:

    typedef FSIOperator super;

    typedef super::fluid_type fluid_type;
    typedef super::solid_type solid_type;

    typedef super::fluid_bchandler_type  bchandler_type;

    typedef fluid_raw_type::vector_type  vector_type;
    typedef fluid_raw_type::matrix_type      matrix_type;
    typedef boost::shared_ptr<matrix_type>        matrix_ptrtype;

    // default constructor
    exactJacobian();

    // destructor
    ~exactJacobian();

    // member functions

    void evalResidual(vector_type       &_res,
                      const vector_type &_disp,
                      const UInt          _iter);
    void solveJac(vector_type       &_muk,
                  const vector_type &_res,
                  const Real       _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void setupFEspace();

    void setupFluidSolid();

    void setDataFile( GetPot const& data );

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

    void bcManageVec   ( bchandler_type& bch, vector_type& rhs );

    void setFluidInterfaceDisp      (vector_type &disp, UInt type = 0);

    void registerMyProducts( );
//    void setFluidLinInterfaceStress (vector_type &disp, UInt type = 0);

//     bc_vector_interface bcvFluidLoadToSolid()       {return M_bcvFluidLoadToSolid;}
//     bc_vector_interface bcvFluidLinInterfaceDisp()    {return M_bcvFluidLinInterfaceDisp;}
// //     bc_vector_interface bcvFluidInvInterfaceDisp()    {return M_bcvFluidInvInterfaceDisp;}
// //     bc_vector_interface bcvFluidInvLinInterfaceDisp() {return M_bcvFluidInvLinInterfaceDisp;}

//     bc_vector_interface bcvReducedFluidInterfaceAcc()    {return M_bcvReducedFluidInterfaceAcc;}
//     bc_vector_interface bcvReducedFluidInvInterfaceAcc() {return M_bcvReducedFluidInvInterfaceAcc;}

//     struct dataJacobian
//     {

//         dataJacobian():
//             M_pFS( 0 )
//             {}

//         dataJacobian(exactJacobian* oper):
//             M_pFS(oper)
//             {}

//         exactJacobian* M_pFS;
//     };

    //Displayer const& getDisplayer(){return displayer();}
private:

    UInt imposeFlux( );

    int M_updateEvery;

//     boost::shared_ptr<vector_type>       M_residualFSI;

//     boost::shared_ptr<vector_type>       M_interfaceDisplacement;
//     boost::shared_ptr<vector_type>       M_interfaceStress;

    boost::shared_ptr<vector_type>       M_rhsNew;
    boost::shared_ptr<vector_type>       M_beta;

    generalizedAitken<vector_type> M_aitkFS;


//     bc_vector_interface     M_bcvFluidLoadToStructure;

    void eval(const vector_type& _res, UInt status);

//    dataJacobian                         M_dataJacobian;

    LifeV::SolverTrilinos        M_linearSolver;

    boost::shared_ptr<Epetra_ExactJacobian> M_epetraOper;

    matrix_ptrtype M_matrShapeDer;
    bool M_recomputeShapeDer;

}; // end class exactJacobian


class Epetra_ExactJacobian:
        public Epetra_Operator
{

public:

    typedef exactJacobian::vector_type  vector_type;

    Epetra_ExactJacobian(exactJacobian* ej):
            M_ej               (ej),
            M_operatorDomainMap(*M_ej->solidInterfaceMap()->getMap(Repeated)),
            M_operatorRangeMap (*M_ej->solidInterfaceMap()->getMap(Repeated)),
            M_comm             (M_ej->worldComm())
    {
//             std::cout << ej << std::endl;
//             std::cout << M_ej->fluidInterfaceMap().getEpetra_Map() << std::endl;
//             std::cout << M_ej->solidInterfaceMap().getEpetra_Map() << std::endl;
//             std::cout << "ok" << std::endl;
    };

    virtual ~Epetra_ExactJacobian() {};

    int 	SetUseTranspose (bool  /*UseTranspose*/)
    {std::cout << "********* EJ : transpose not available\n"; return -1;}
    int 	Apply           (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    int 	ApplyInverse    (const Epetra_MultiVector &/*X*/, Epetra_MultiVector &/*Y*/) const
    {std::cout << "********* EJ : inverse not available\n"; return -1;}
    double 	NormInf         () const
    {std::cout << "********* EJ : NormInf not available\n"; return 1.;}
    const char * Label      () const {return "exactJacobian";}
    bool 	UseTranspose    () const {return false;}
    bool 	HasNormInf      () const {return false;}

    const Epetra_Comm&  Comm () const { return *M_comm; }
    const Epetra_Map & 	OperatorDomainMap () const {return M_operatorDomainMap;}
    const Epetra_Map & 	OperatorRangeMap  () const {return M_operatorRangeMap;}

    void setOperator( exactJacobian* ej) {M_ej = ej;}

private:



    exactJacobian*                       M_ej;

    const Epetra_Map                     M_operatorDomainMap;
    const Epetra_Map                     M_operatorRangeMap;

    boost::shared_ptr<Epetra_Comm>       M_comm;

};


Real fzeroEJ(const Real& t,
             const Real& x,
             const Real& y,
             const Real& z,
             const ID& i);


inline FSIOperator* createEJ() { return new exactJacobian(); }

namespace
{
static bool registerEJ = FSIOperator::FSIFactory::instance().registerProduct( "exactJacobian", &createEJ );
}

}

#endif
#endif
