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
/**
 * \file Monolithic.hpp
 * \author Paolo Crosetto
 * \date 13-09-2008
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
                        const vector_type& _sol,
                        const int          _iter);
    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */


    virtual void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const double       _linearRelTol);

    /**
       solves the Jacobian system
       \param _muk: output, solution at the current Newton step
       \param _res: nonlinear residual
       \param _linearRelTol: not used
    */

    virtual void updateSystem(const vector_type& _sol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */

    void setDataFromGetPot( GetPot const& data );

    void monolithicToInterface(    vector_type& _lambdaSolid, const vector_type& _sol) ;

    /**
       restricts a vector with a monolithic map on the solid interface map
       \param _lambdasolid: vector on the solid interface
       \param _disp: monolithic vector
    */

    void monolithicToSolid(const vector_type& _disp, vector_type& _dispSolid);
    /**
       restricts a vector with a monolithic map on the solid map
       \param _dispSolid: vector on the solid domain
       \param _disp: monolithic vector
    */

    void setDispSolid(const vector_type& sol);
    /**
       sets the vector M_solid->dispSolid() to the monolithic solution M_solid->disp() in the solid nodes, to 0 in the fluid nodes
    */

    void monolithicToFluid(const vector_type& _disp, vector_type& _dispFluid);
    /**
        restricts a vector with a monolithic map on the fluid map
       \param _dispFluid: vector on the fluid domain
       \param _disp: monolithic vector
    */

    void iterateMesh(const vector_type& disp);
    /**
       iterates the mesh
       \param disp: monolithic solution
    */

    void buildSystem();
    /**
      builds the constant part of the monolithic matrix
    */
    virtual void setup();
    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map
    */

    virtual    boost::shared_ptr<EpetraMap>& couplingVariableMap(){return M_monolithicMap;}
    /** returns the monolithic map*/

    void shapeDerivatives(vector_ptrtype rhsShapeDerivatives, vector_ptrtype meshDeltaDisp, const vector_type& sol);
    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    //! getters
    vector_type getRhsShapeDerivatives(){return *M_rhsShapeDerivatives;}
    //    const boost::shared_ptr<EpetraMap>& monolithicMap() {return M_monolithicMap;}


    const UInt dimInterface() const {return nDimensions*M_interface ;}
    void setUpSystem( GetPot const& data_file );
    void setUp( const GetPot& dataFile );

protected:

template <typename SolverType>
void iterateMonolithic(vector_type& rhs, vector_type& step, matrix_ptrtype prec, SolverType linearSolver);
    void addDiagonalEntries(Real entry, matrix_ptrtype matrix, EpetraMap& Map, UInt offset=0, bool replacs=false);
    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
   */

        void addDiagonalEntries(Real entry, matrix_ptrtype matrix, std::map<ID, ID> const& Map, UInt offset=0, bool replacs=false);
    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
   */

    virtual void couplingMatrix(matrix_ptrtype & bigMatrix, bool couplingSolid = true);
    /**
      \small  to the monolithic matrix.
      \param bigMatrix: matrix to be modified
    */
    void couplingRhs( vector_ptrtype rhs , vector_ptrtype un );
    /**
       \small adds the part due to coupling to the rhs
       \param rhs: right hand side
       \param un: current solution
    */
    void robinCoupling(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse = false*/);
    /**
       \small adds the matrix [0,0,0,0;0,0,0,alphaf;0,0,0,0;0,alphas,0,0] to IdentityMatrix
       \param IdentityMatrix: a matrix, if it is an identity the method builds the Robin-Robin linear combination of boundary conditions
       \param alphaf: coefficient
       \param alphas: coefficient
       \param inverse if true the inverse of the linear combination is computed. the output matrix is [];
    */

    void robinCouplingInv(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas/*, bool inverse*/); // not working with non-matching grids
    //not used now

    void zeroBlock( matrix_ptrtype matrixPtr,vector_type& colNumeration , const std::map<ID, ID>& map, UInt rowOffset=0, UInt colOffset=0);


    void evalResidual(fluid_bchandler_raw_type & bcFluid, solid_bchandler_raw_type & bcSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling=false); // used for monolithic

    void evalResidual( const vector_type& sol, const vector_ptrtype& rhs,  vector_type& res, bool diagonalScaling=false); // used for monolithic
    void evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner);

    void setBlockPreconditioner(matrix_ptrtype blockPrec){
*blockPrec += *M_solid->getMassStiff();
 bcManageMatrix( *blockPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1.,
                  dataSolid().time() );
 bcManageMatrix( *blockPrec, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u, M_uFESpace->feBd(), 1.,
              dataSolid().time() );

}
    void setFullPreconditioner(matrix_ptrtype fullPrec){
        *fullPrec += *M_monolithicMatrix;
        fullPrec->GlobalAssemble();
//     bcManageMatrix( *fullPrec, *M_dFESpace->mesh(), M_dFESpace->dof(), *M_BCh_d, M_dFESpace->feBd(), 1.,
//                   dataSolid().time() );
//     bcManageMatrix( *fullPrec, *M_uFESpace->mesh(), M_uFESpace->dof(), *M_BCh_u, M_uFESpace->feBd(), 1.,
//               dataSolid().time() );
}
    bool recomputePrec(){return(!M_reusePrec || M_resetPrec);}
    void setOperator(Epetra_Operator& epetraOperator){M_linearSolver->setOperator(epetraOperator);}
    void setMatrix(){M_linearSolver->setMatrix(*M_monolithicMatrix);}

    void applyPreconditioner(matrix_ptrtype robinCoupling);
    void updateCoupling(matrix_type couplingMatrix);
    void updateMatrix(matrix_type & bigMatrixStokes);
    void updateSolidSystem(vector_ptrtype& rhsFluidCoupling);
    void    diagonalScale(vector_type&, matrix_ptrtype matrFull);

    void shiftSolution(){}

    void couplingVariableExtrap(vector_ptrtype& lambda, vector_ptrtype& /*lambdaDot*/, bool& /*firstIter*/)
{    leaderPrint("norm( solution ) init = ", lambda->NormInf() );};

protected:

    void resetHeAndFluid(){}
    void solidInit(const RefFE* refFE_struct, const LifeV::QuadRule* bdQr_struct, const LifeV::QuadRule* qR_struct);
    void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);


    boost::shared_ptr<EpetraMap>                      M_monolithicMap;
    matrix_ptrtype                                    M_couplingMatrix;
    UInt                                              M_interface;
    EpetraMap                                         M_interfaceMap;///the solid interface map
    boost::shared_ptr<vector_type>                    M_beta;
    matrix_ptrtype                                    M_monolithicMatrix;
    matrix_ptrtype                                    M_bigPrecPtr;
    int                                               M_DDBlockPrec;
    boost::shared_ptr<vector_type>                    M_rhsFull;

private:
    int                                               M_updateEvery;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
    boost::shared_ptr<vector_type>                    M_rhsShapeDerivatives;
    boost::shared_ptr<vector_type>                    M_rhsNew;
    bool                                              M_fullMonolithic;
    Real                                              M_entry;
    matrix_ptrtype                                    M_robinCoupling;
    //    matrix_ptrtype                                    M_robinCouplingInv;
    Real                                              M_alphaf;
    Real                                              M_alphas;
    bool                                              M_diagonalScale;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
    bool                                              M_reusePrec;
    bool                                              M_resetPrec;
    UInt                                              M_maxIterSolver;
    boost::shared_ptr<SolverTrilinos>                 M_linearSolver;
    //    UInt                                              M_maxIterForReuse;
// =======
//     bool                                              M_isDiagonalBlockPrec;
//     boost::shared_ptr<vector_type>                    M_beta;
//     matrix_ptrtype                                    M_matrFull;
//     boost::shared_ptr<SolverTrilinos>                 M_linearSolver;
//     bool                                              M_reusePrec;
//     bool                                              M_resetPrec;
//     UInt                                              M_maxIterSolver;
//     matrix_ptrtype                                    M_couplingMatrix;
//     //    UInt                                              M_maxIterForReuse;
// >>>>>>> 1.3

    //    boost::shared_ptr<Epetra_IntVector>               M_numerationInterfaceInt;


};


template <typename SolverType>
void Monolithic::
iterateMonolithic(vector_type& rhs, vector_type& step, matrix_ptrtype prec, SolverType linearSolver)
{
    Chrono chrono;

    M_solid->getDisplayer().leaderPrint("    Solving the system ... \n" );

    M_solid->getDisplayer().leaderPrint("    Updating the boundary conditions ... ");

    chrono.start();

    // boundary conditions applied in the residual evaluation

    //M_matrFull->spy("jacobian");
    //prec->spy("blockPreconditioner");
    chrono.stop();
    M_solid->getDisplayer().leaderPrintMax("done in ", chrono.diff() );

    //M_comm->Barrier();

    //M_monolithicMatrix->spy("jacobian");


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
    int numIter = M_linearSolver->solveSystem(/*M_monolithicMatrix,*/ rhs, step, prec, (M_reusePrec)&&(!M_resetPrec));

    if (numIter < 0)
    {
        chrono.start();

        M_solid->getDisplayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

        //M_prec->buildPreconditioner(M_solid->getMatrixPtr());

        //double condest = M_prec->Condest();

        //M_linearSolver->setPreconditioner(M_prec);

        if (numIter <= -M_maxIterSolver)
           M_solid->getDisplayer().leaderPrint("   ERROR: Iterative solver failed twice.\n");
    }

    //    M_disp += step;

    M_solid->getDisplayer().leaderPrint("   system solved.\n ");

    //    M_dispSolid.spy("dispSolid0");
}

}
#endif
