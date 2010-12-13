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
    @file
    @brief DataFSI - File containing a data container for FSI problems

    @author Miguel Fernandez
    @author Gilles Fourestey
    @author Paolo Crosetto <paolo.crosetto@epfl.ch>
    @date 10-06-2010

    @contributor Simone Deparis <simone.deparis@epfl.ch>
    @maintainer Simone Deparis <simone.deparis@epfl.ch>
 */


#ifndef EXACTJACOBIANBASE_HPP
#define EXACTJACOBIANBASE_HPP

#include <life/lifesolver/FSIOperator.hpp>

namespace LifeV
{

class Epetra_ExactJacobian;

//! ExampleClass - Short description of the class
/*!
\include ../../doc/api/bibliography/newton.dox

    @author Name Surname
    @see  \ref{FM05}
    FSIOperator

    This class implements an FSIOperator whose Jacobian is computed exacly, i.e., using
    shape dericatives.

*/

class exactJacobian : public FSIOperator
{
public:

    //! @name Public Types
    //@{
    typedef FSIOperator                     super;

    typedef super::vector_Type              vector_Type;
    typedef super::vectorPtr_Type           vectorPtr_type;

    typedef fluid_Type::matrix_Type         matrix_Type;
    typedef fluid_Type::matrixPtr_Type      matrixPtr_Type;

//! OBSOLETE typedefs
    //typedef super::fluid_Type               fluid_Type;
    //typedef super::solid_Type               solid_Type;

    //typedef super::fluidBchandler_Type      fluidBchandler_Type;

    typedef super::fluid_Type               fluid_type;
    typedef super::solid_Type               solid_type;
    typedef super::vector_Type              vector_type;

    typedef super::fluidBchandler_Type      bchandler_type;

    typedef fluid_Type::matrix_Type         matrix_types;
    typedef fluid_Type::matrixPtr_Type      matrix_ptrtype;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    exactJacobian();

    //! Destructor
    ~exactJacobian();
    //@}

    //! @name Methods
    //@{

    void evalResidual(vector_Type       &_res,
                      const vector_Type &_disp,
                      const UInt          _iter);
    void solveJac(vector_Type       &_muk,
                  const vector_Type &_res,
                  const Real       _linearRelTol);

    void solveLinearFluid();
    void solveLinearSolid();

    void setupFEspace();

    void setupFluidSolid();

    void setDataFile( GetPot const& data );


    //
    // BCs treatment
    //

    // Solid, Lin. Solid, and inverses
    void setSolidInterfaceDisp         (vector_Type& disp  , UInt type = 0);
    void setSolidLinInterfaceDisp      (vector_Type& disp  , UInt type = 0);
    void setSolidInvLinInterfaceStress (vector_Type& stress, UInt type = 0);
    void setReducedFluidInterfaceAcc   (vector_Type& acc   , UInt type = 0);
    void setReducedFluidInvInterfaceAcc(vector_Type& acc   , UInt type = 0);


    bc_vector_interface bcvFluidLoadToStructure()       {return M_bcvFluidLoadToStructure;}

    // Fluid, Lin. Fluid, and inverses

    void bcManageVec   ( super::fluidBchandler_Type& bch, vector_Type& rhs );

    void setFluidInterfaceDisp      (vector_Type &disp, UInt type = 0);

    void registerMyProducts( );

    //! Display general information about the content of the class
    /*!
        List of things displayed in the class
        @param output specify the output format (std::cout by default)
     */
    void showMe( std::ostream& output = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //@}

private:

    //! @name Private Methods
    //@{

    UInt imposeFlux( );

    //@}

    int M_updateEvery;


    vectorPtr_Type       M_rhsNew;
    vectorPtr_Type       M_beta;

    generalizedAitken<vector_Type> M_aitkFS;


//     bc_vector_interface     M_bcvFluidLoadToStructure;

    void eval(const vector_Type& _res, UInt status);

//    dataJacobian                         M_dataJacobian;

    LifeV::SolverTrilinos        M_linearSolver;

    boost::shared_ptr<Epetra_ExactJacobian> M_epetraOper;

    matrixPtr_Type M_matrShapeDer;
    bool M_recomputeShapeDer;

}; // end class exactJacobian


class Epetra_ExactJacobian:
        public Epetra_Operator
{

public:

    typedef exactJacobian::vector_Type  vector_Type;
    // OBSOLETE typedef
    typedef exactJacobian::vector_Type  vector_type;

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

}  // Namespace LifeV

#endif // EXACTJACOBIANBASE_HPP
