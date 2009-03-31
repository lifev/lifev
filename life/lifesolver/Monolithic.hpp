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

#ifndef _MONOLITHIC_HPP
#define _MONOLITHIC_HPP


#include <life/lifesolver/FSIOperator.hpp>
//#include <Epetra_IntVector.h>

namespace LifeV
{
class Monolithic : public FSIOperator
{
public:


    typedef FSIOperator                            super;
    typedef FSIOperator::fluid_type::value_type::matrix_type   matrix_type;
    typedef boost::shared_ptr<matrix_type>                    matrix_ptrtype;


    // constructors
    Monolithic();

    // destructor

    ~Monolithic();



    void   evalResidual(vector_type&        res,
                        const vector_type& _disp,
                        const int          _iter);

    void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const double       _linearRelTol);

    template <typename SolverType>
    void iterateMonolithic(vector_type& rhs, vector_type& sol, matrix_ptrtype prec, SolverType linearSolver); // used for monolithic

    void updateSystem(const vector_type& displacement);

    void setDataFromGetPot( GetPot const& data );

    void monolithicToInterface(    vector_type& lambdaSolid, const vector_type& disp) ;

    void monolithicToSolid(const vector_type& disp, vector_type& dispSolid);

    void setDispSolid(const vector_type& sol);

    void monolithicToFluid(const vector_type& disp, vector_type& dispFluid);

    void iterateMesh(const vector_type& disp);
    //void coupling(matrix_ptrtype & bigMatrix, vector_ptrtype rhs = vector_ptrtype());
    void couplingMatrix(matrix_ptrtype & bigMatrix);
    void couplingRhs(vector_ptrtype rhs, vector_ptrtype un);

    void buildSystem();
    void addDiagonalEntries(Real& entry, matrix_ptrtype & bigMatrix);
    void setup();

    virtual    boost::shared_ptr<EpetraMap>& couplingVariableMap(){return M_monolithicMap;}

    //void evalResidual(bchandler_raw_type & bcFluid, bchandler_raw_type & bcSolid, const vector_type& sol, vector_type& res/*, matrix_type& bigMatrix*/); // used for monolithic
    void evalResidual(fluid_bchandler_raw_type & bcFluid, solid_bchandler_raw_type & bcSolid, const vector_type& sol, const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling=false); // used for monolithic

    void evalResidual( const vector_type& sol, const vector_ptrtype& rhs,  vector_type& res, bool diagonalScaling=false); // used for monolithic
    void evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, const vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner);

    //void evalResidual(fluid_bchandler_raw_type & bcFluid, solid_bchandler_raw_type & bcSolid, const vector_type& sol, vector_type& res/*, matrix_type& bigMatrix*/); // used for monolithic
    void setBlockPreconditioner(matrix_ptrtype blockPrec){*blockPrec += *M_solid->getMassStiff();}
    void setFullPreconditioner(matrix_ptrtype fullPrec){*fullPrec += *M_monolithicMatrix;}
    bool recomputePrec(){return(!M_reusePrec || M_resetPrec);}
    void setOperator(Epetra_Operator& epetraOperator){M_linearSolver->setOperator(epetraOperator);}
    void setMatrix(){M_linearSolver->setMatrix(*M_monolithicMatrix);}
    //void setOverlappingPreconditioner(matrix_ptrtype fullPrec, short overlap);
    void applyPreconditioner(matrix_ptrtype robinCoupling);
    //    void setAztecooPreconditioner(const GetPot& dataFile, const std::string& section){M_linearSolver->setAztecooPreconditioner( dataFile, section);}
    void updateCoupling(matrix_type couplingMatrix);
    void updateMatrix(matrix_type & bigMatrixStokes);
    void updateSolidSystem(vector_ptrtype& rhsFluidCoupling);
    void setUpSystem( GetPot const& data_file );
    void setUp( const GetPot& dataFile );
    void    diagonalScale(vector_type&, matrix_ptrtype matrFull);
    void shiftSolution(){}

    void couplingVariableExtrap(vector_type& lambda, vector_type& /*lambdaDot*/, bool& /*firstIter*/)
{    leaderPrint("norm( solution ) init = ", lambda.NormInf() );};

protected:

    void resetHeAndFluid(){}
    void solidInit(const RefFE* refFE_struct, const LifeV::QuadRule* bdQr_struct, const LifeV::QuadRule* qR_struct);
    void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);

private:

    int                                               M_updateEvery;
    bool                                              firstIter;
    boost::shared_ptr<EpetraMap>&                     M_monolithicMap;
    UInt                                              M_interface;
    bool                                              M_semiImplicit;
    EpetraMap                                         M_interfaceMap;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
    bool                                              M_isDiagonalBlockPrec;
    boost::shared_ptr<vector_type>                    M_beta;
    matrix_ptrtype                                    M_matrFull;
    boost::shared_ptr<SolverTrilinos>                 M_linearSolver;
    bool                                              M_reusePrec;
    bool                                              M_resetPrec;
    UInt                                              M_maxIterSolver;
    matrix_ptrtype                                    M_couplingMatrix;
    matrix_ptrtype                                    M_monolithicMatrix;
    //    UInt                                              M_maxIterForReuse;

    //    boost::shared_ptr<Epetra_IntVector>               M_numerationInterfaceInt;


};


template <typename SolverType>
void Monolithic::
iterateMonolithic(vector_type& rhs, vector_type& step, matrix_ptrtype prec, SolverType linearSolver)
{
    Chrono chrono;

    M_solid->getDisplayer().leaderPrint("  S-  Solving the system ... \n" );

    M_solid->getDisplayer().leaderPrint("  S-  Updating the boundary conditions ... ");

    chrono.start();

    // boundary conditions applied in the residual evaluation

    //M_matrFull->spy("jacobian");
    //prec->spy("blockPreconditioner");
    chrono.stop();
    M_solid->getDisplayer().leaderPrintMax("done in ", chrono.diff() );

    //M_comm->Barrier();

    M_monolithicMatrix->spy("jacobian");


    /*        for (UInt i = 7980; i < 9240; i++ )//lines to kill
        {
            double* tmp;
            int err;
            int entries;
            err = M_matrFull->getEpetraMatrix().ExtractGlobalRowView(i, entries, tmp);
            if(entries == 0)
                {
                    std::cout<<"ERROR in line " << i << " err " << err << std::endl;
                    //                    break;
                }
            //            std::cout << "nonzero entries for row " << i << " ==> " << entries << std::endl;
            }*/

    //M_comm->Barrier();


    //    M_disp.spy("disp0");
    int numIter = M_linearSolver->solveSystem(M_monolithicMatrix, rhs, step, prec, (M_reusePrec)&&(!M_resetPrec));

    if (numIter < 0)
    {
        chrono.start();

        M_solid->getDisplayer().leaderPrint("  s- Iterative solver failed, numiter = ", -numIter );

        //M_prec->buildPreconditioner(M_solid->getMatrixPtr());

        //double condest = M_prec->Condest();

        //M_linearSolver->setPreconditioner(M_prec);

        if (numIter <= -M_maxIterSolver)
           M_solid->getDisplayer().leaderPrint("  s- ERROR: Iterative solver failed twice.\n");
    }

    //    M_disp += step;

    M_solid->getDisplayer().leaderPrint("  S- system solved.\n ");

    //    M_dispSolid.spy("dispSolid0");
}


}
#endif
