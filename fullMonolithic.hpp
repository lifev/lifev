/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Paolo Crosetto <crosetto@iacspc70.epfl.ch>
       Date: 2008-09-18

  Copyright (C) 2008

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
   \file fullMonolithic.hpp
   \author crosetto <Paolo Crosetto>
   \date 18/09/2008
 */
#ifndef _FULLMONOLITHIC_HPP
#define _FULLMONOLITHIC_HPP

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV
{
class Epetra_FullMonolithic;

class fullMonolithic : public Monolithic
{
public:

    typedef Monolithic super;

    fullMonolithic();
    ///constructor


    virtual ~fullMonolithic(){}
    ///destructor

    void couplingMatrix(matrix_ptrtype & bigMatrix, bool solidCoupling=true);
    void couplingRhs(vector_ptrtype rhs);

    void   evalResidual(vector_type&        res,
                        const vector_type& _sol,
                        const int          _iter);
    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: fluid domain displacement solution
       \param iter: current nonLinRichardson (block Gauss Seidel for the tangent system) iteration
    */

    /*    void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const double       _linearRelTol);*/
    /**
       solves the Jacobian system
       \param _muk: output, solution at the current block GS step
       \param _res: linear residual
       \param _linearRelTol: not used
    */
    //    void   updateSystem(const vector_type&       solution);

    void setup();
    //    void updateSystem(const vector_type& _sol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */
    void buildSystem ();
    EpetraMap& mapWithoutMesh(){return *M_mapWithoutMesh;}
    //    vector_type& solution(){return *M_un;}
    void updateSystem(const vector_type& displacement);

    void solveJac(vector_type       &_muk,
                  const vector_type &_res,
                  const double       _linearRelTol);
    matrix_ptrtype getMatrixPtr(){return this->M_monolithicMatrix;}
    vector_ptrtype uk(){return M_uk;}
    vector_type& meshVel();
    private:
    boost::shared_ptr<EpetraMap>   M_mapWithoutMesh;
    LifeV::SolverTrilinos        M_linearSolver;
    boost::shared_ptr<Epetra_FullMonolithic> M_epetraOper;
    vector_ptrtype                       M_uk;
    vector_ptrtype                       M_meshVel;
    //    static bool              reg;
    vector_ptrtype                       M_unOld;//****************************
};

class Epetra_FullMonolithic : public Epetra_Operator
{
public :

    typedef Monolithic::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    Epetra_FullMonolithic(fullMonolithic& MOperator/*, Epetra_Map& Mappa*/):
        M_FMOper               (&MOperator),
        M_comm             (&M_FMOper->worldComm()),
        M_operatorDomainMap(/*Mappa*/ *M_FMOper->couplingVariableMap()->getMap(Repeated)),
        M_operatorRangeMap(/*Mappa*/ *M_FMOper->couplingVariableMap()->getMap(Repeated))
    {
        std::cout << "M_FMOper->solidInterfaceMap() = " << M_FMOper->solidInterfaceMap() << std::endl;
        std::cout << "M_FMOper->solidInterfaceMap()->getMap(Repeated) = " << M_FMOper->solidInterfaceMap()->getMap(Repeated) << std::endl;

    };
    virtual ~Epetra_FullMonolithic(){};

    int 	SetUseTranspose (bool  /*UseTranspose*/)
        {std::cout << "********* EJ : transpose not available\n"; return -1;}
    int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
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

    void setOperator( fullMonolithic* fm) {M_FMOper = fm;}

private:

    fullMonolithic* M_FMOper;
    const Epetra_Map M_operatorDomainMap;
    const Epetra_Map M_operatorRangeMap;
    Epetra_Comm* M_comm;
};

#ifdef UNDEFINED
class Epetra_FullMonolithic : public Epetra_Operator
{
public :

    typedef Monolithic::vector_type  vector_type;
    typedef boost::shared_ptr<vector_type> vector_ptrtype;

    Epetra_FullMonolithic(Monolithic* MOperator):
        M_FMOper               (MOperator),
        //        M_operatorDomainMap(*M_FMOper->solidInterfaceMap()->getMap(Repeated)),
        //        M_operatorRangeMap(*M_FMOper->solidInterfaceMap()->getMap(Repeated)),
        M_comm             (&M_FMOper->worldComm())
    {
        std::cout << "M_FMOper->solidInterfaceMap() = " << M_FMOper->solidInterfaceMap() << std::endl;
        std::cout << "M_FMOper->solidInterfaceMap()->getMap(Repeated) = " << M_FMOper->solidInterfaceMap()->getMap(Repeated) << std::endl;

    };
    virtual ~Epetra_FullMonolithic(){};

    int 	SetUseTranspose (bool  /*UseTranspose*/)
        {std::cout << "********* EJ : transpose not available\n"; return -1;}
    int Apply (const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;
    int 	ApplyInverse    (const Epetra_MultiVector &/*X*/, Epetra_MultiVector &/*Y*/) const
        {std::cout << "********* EJ : inverse not available\n"; return -1;}
    double 	NormInf         () const
        {std::cout << "********* EJ : NormInf not available\n"; return 1.;}
    const char * Label      () const {return "exactJacobian";}
    bool 	UseTranspose    () const {return false;}
    bool 	HasNormInf      () const {return false;}
    const Epetra_Comm&  Comm () const { return *M_comm; }
    const Epetra_Map & 	OperatorDomainMap () const {/*return M_operatorDomainMap;*/}
    const Epetra_Map & 	OperatorRangeMap  () const {/*return M_operatorRangeMap;*/}

    void setOperator( Monolithic* fm) {M_FMOper = fm;}

private:

    //    const Epetra_Map M_operatorDomainMap;
    //    const Epetra_Map M_operatorRangeMap;
    Monolithic* M_FMOper;
    Epetra_Comm* M_comm;
};
#endif
}
#endif
