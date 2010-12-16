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
    @brief Implementation of an  FSIOperator with Newton algorithm.

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

//! exactJacobian - Implementation of an  FSIOperator with Newton algorithm.
/*!
\include ../../doc/api/bibliography/newton.dox

    @author Miguel Fernandez
    @author Gilles Fourestey
    @author Paolo Crosetto <paolo.crosetto@epfl.ch>
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

    //typedef super::fluid_Type               fluid_Type;
    typedef super::solid_Type               solid_Type;

    //! OBSOLETE typedefs
//     //typedef super::fluidBchandler_Type      fluidBchandler_Type;

//     typedef super::fluid_Type               fluid_type;
//     typedef super::solid_Type               solid_type;
//     typedef super::vector_Type              vector_Type;

//     typedef super::fluidBchandler_Type      bchandler_type;

//     typedef fluid_Type::matrix_Type         matrix_types;
//     typedef fluid_Type::matrixPtr_Type      matrix_ptrtype;

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

    //! solves the Jacobian system
    /**
       The implementation of this method distinguish the various FSI formulations which derive from this class.
       For this reason it must be pure virtual, snd implemented in the child classes.
       \param muk: unknown solution at the k-th nonlinear iteration
       \param res: residual vector (the right hand side of the Jacobian system)
       \param linearRelTol: tolerance for the nonlinear solver
       \todo{replace Real with Real& }
     */
    void solveJac(vector_Type       &_muk,
                  const vector_Type &_res,
                  const Real       _linearRelTol);

    //! Evaluates the nonlinear residual of the FSI system
    /**
       The implementation of this method also depends on the child classes, though it does not characterize them.
       \param res:  residual vector  to be computed
       \param disp: current unknown solution
       \param iter: nonlinear iteration counter. The part of th rhs related to the time discretization is computed only for iter=0
    */
    void evalResidual(vector_Type       &_res,
                      const vector_Type &_disp,
                      const UInt          _iter);

    //! Solves the linear fluid problem
    /** this method is called only by the class Epetra_ExactJacobian
     */
    void solveLinearFluid();

    //! Solves the linear structure problem
    /** this method is called only by the class Epetra_ExactJacobian
     */
    void solveLinearSolid();

    //! sets the space discretization parameters
    /*!
      The FE discretization is set accordingly to what specified in the FSIdata member (order of accuracy for the fluid
      pressure, velocity and for the structure).
     */
    void setupFEspace();

    //! setup of the fluid and solid solver classes
    /**
       This method computes the number of fluxes assigned at the boundaries and calls setupFluidSolid(UInt fluxes)
     */
    void setupFluidSolid();

    //! initializes the GetPot data file
    void setDataFile( GetPot const& data );

    //! should call bcManage for a vector, but the implementation is empty
    void bcManageVec   ( super::fluidBchandler_Type& bch, vector_Type& rhs ) {};

    //! register the product for the factory
    void registerMyProducts( );

    //@}

private:

    //! @name Private Methods
    //@{

    UInt imposeFlux( );

    void eval(const vector_Type& _res, UInt status);


    //@}

//! Epetra_ExactJacobian  This class implements an Epetra_Operator to be passed to AztecOO.
/*!

    @author Gilles Fourestey
    @see
    exactJacobian

    This class relies on exactJacobian to solve the linear jacobian problem

*/

    class Epetra_ExactJacobian:
        public Epetra_Operator
    {

    public:

        typedef exactJacobian::vector_Type  vector_Type;
        typedef Epetra_Map                  map_Type;
        typedef boost::shared_ptr<map_Type> mapPtr_Type;

        // OBSOLETE typedef
//         typedef exactJacobian::vector_Type  vector_Type;

        //! @name Constructor & Destructor
        //@{
        //! Empty Constructor
        Epetra_ExactJacobian();

        //! Destructor
        virtual ~Epetra_ExactJacobian() {};
        //@}

        //! @name Methods
        //@{

        //! sets the exactJacobian pointer and some contents thereof
        void setOperator(exactJacobian* ej);

        //! apply the jacobian to X and returns the result in Y
        int 	Apply           (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

        //! These are the methods necessary to implement Epetra_Operator but that are not used.
        int 	SetUseTranspose (bool  /*UseTranspose*/)
        {std::cout << "********* EJ : transpose not available\n"; return -1;}

        int 	ApplyInverse    (const Epetra_MultiVector &/*X*/, Epetra_MultiVector &/*Y*/) const
        {std::cout << "********* EJ : inverse not available\n"; return -1;}
        double 	NormInf         () const
        {std::cout << "********* EJ : NormInf not available\n"; return 1.;}
        const char * Label      () const {return "exactJacobian";}
        bool 	UseTranspose    () const {return false;}
        bool 	HasNormInf      () const {return false;}

        const Epetra_Comm&  Comm () const { return *M_comm; }
        const Epetra_Map & 	OperatorDomainMap () const {return *M_operatorDomainMap;}
        const Epetra_Map & 	OperatorRangeMap  () const {return *M_operatorRangeMap;}
        //@}


    private:

        exactJacobian*                  M_ej;

        mapPtr_Type                     M_operatorDomainMap;
        mapPtr_Type                     M_operatorRangeMap;

        boost::shared_ptr<Epetra_Comm>  M_comm;

    }; // end of class Epetra_ExactJacobian


    vectorPtr_Type       M_rhsNew;
    vectorPtr_Type       M_beta;

    generalizedAitken<vector_Type> M_aitkFS;

    LifeV::SolverTrilinos  M_linearSolver;
    Epetra_ExactJacobian   M_epetraOper;

    matrixPtr_Type M_matrShapeDer;
    bool           M_recomputeShapeDer;


}; // end class exactJacobian


inline FSIOperator* createEJ() { return new exactJacobian(); }

namespace
{
static bool registerEJ = FSIOperator::FSIFactory_Type::instance().registerProduct( "exactJacobian", &createEJ );
}

}  // Namespace LifeV

#endif // EXACTJACOBIANBASE_HPP
