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
    typedef SolverTrilinos                                      solver_type;


    // constructors
    Monolithic();

    // destructor
    ~Monolithic();



    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */
    void   evalResidual(vector_type&        res,
                        const vector_type& _sol,
                        const UInt          _iter);



    /**
       solves the Jacobian system
       \param _muk: output, solution at the current Newton step
       \param _res: nonlinear residual
       \param _linearRelTol: not used
    */
    virtual void   solveJac(vector_type&       _muk,
                    const vector_type& _res,
                    const Real       _linearRelTol);


    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */
    virtual void updateSystem(const vector_type& _sol);

    /** sets the parameters from the data file*/
    void setDataFromGetPot( GetPot const& data );



    /**
       restricts a vector with a monolithic map on the solid interface map
       \param _lambdasolid: vector on the solid interface
       \param _disp: monolithic vector
    */
    void monolithicToInterface(    vector_type& _lambdaSolid, const vector_type& _sol) ;

    /**
       restricts a vector with a monolithic map on the solid map
       \param _dispSolid: vector on the solid domain
       \param _disp: monolithic vector
    */
    void monolithicToSolid(const vector_type& _disp, vector_type& _dispSolid);

    /**
       sets the vector M_solid->dispSolid() to the monolithic solution M_solid->disp() in the solid nodes, to 0 in the fluid nodes
    */
    void setDispSolid(const vector_type& sol);

    /**
        restricts a vector with a monolithic map on the fluid map
       \param _dispFluid: vector on the fluid domain
       \param _disp: monolithic vector
    */
    void monolithicToFluid(const vector_type& _disp, vector_type& _dispFluid);

    /**
       iterates the mesh
       \param disp: monolithic solution
    */
    void iterateMesh(const vector_type& disp);

    /**
      builds the constant part of the monolithic matrix
    */
    void buildSystem();

    /**
       create FEspace
    */
    void setupFEspace();

    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map
    */
    virtual void setupFluidSolid();

    /** returns the monolithic map*/
    virtual    boost::shared_ptr<EpetraMap>& couplingVariableMap(){return M_monolithicMap;}

#if OBSOLETE
    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    void shapeDerivatives(vector_ptrtype rhsShapeDerivatives, vector_ptrtype meshDeltaDisp, const vector_type& sol);
#endif

    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    void shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol);
    //! getters
#if OBSOLETE
    /** get the shape derivatives vector*/
    vector_type getRhsShapeDerivatives(){return *M_rhsShapeDerivatives;}
#endif
    //    const boost::shared_ptr<EpetraMap>& monolithicMap() {return M_monolithicMap;}

    //!get the total dimension of the FS interface
    const UInt dimInterface() const {return nDimensions*M_interface ;}

    void setupSystem( );
    void setUp( const GetPot& dataFile );

protected:

template <typename SolverType>
/**
\small solves the monolithic system, once a solver, a preconditioner and a rhs have been defined.
 */
void iterateMonolithic(vector_type& rhs, vector_type& step, matrix_ptrtype prec, SolverType linearSolver);
    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
       \param offset: offset from which starts the insertion
       \param replace=false: flag stating wether the diagonal elements already exist or have to be inserted in the map
   */
    void addDiagonalEntries(Real entry, matrix_ptrtype matrix, EpetraMap& Map, UInt offset=0, bool replacs=false);

    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: standard map specifying the location where to add diagonal entries
       \param offset: offset from which starts the insertion
       \param replace=false: flag stating wether the diagonal elements already exist or have to be inserted in the map
   */
        void addDiagonalEntries(Real entry, matrix_ptrtype matrix, std::map<ID, ID> const& Map, UInt offset=0, bool replace=false);

    /**
      \small builds the coupling matrix.
      \param bigMatrix: the coupling matrix to be built
      \param couplingSolid: if false the solid coupling part is neglected in the matrix
    */
    virtual void couplingMatrix(matrix_ptrtype & bigMatrix, bool couplingSolid = true);

    /**
       \small adds the part due to coupling to the rhs
       \param rhs: right hand side
       \param un: current solution
    */
    void couplingRhs( vector_ptrtype rhs , vector_ptrtype un );

    /**
       \small builds the matrix [1,0,0,0;0,1,0,alphaf;0,0,1,0;0,alphas,0,1] to IdentityMatrix. This matrix is intended to be a preconditioner and has not been tested yet.
       \param IdentityMatrix: a matrix, if it is an identity the method builds the Robin-Robin linear combination of boundary conditions
       \param alphaf: coefficient
       \param alphas: coefficient
       \param inverse if true the inverse of the linear combination is computed. the output matrix is [];
    */
    void robinCoupling(matrix_ptrtype & IdentityMatrix, Real& alphaf, Real& alphas);


    /**
       \small Method to replace a block in the matrix with a zero block. It is inefficient and works only serially. It is implemented only for testing purposes.
    */
    void zeroBlock( matrix_ptrtype matrixPtr,vector_type& colNumeration , const std::map<ID, ID>& map, UInt rowOffset=0, UInt colOffset=0);

    /**\small evaluates the nonlinear residual
     \param bcFluid: fluid BCHandler
     \param bcSolid: solid BCHandler
     \param sol: solution vector
     \param rhs: right-hand side
     \param res: the output residual
     \param diagonalScaling: flag stating wether to perform diagonal scaling
*/
    void evalResidual(fluid_bchandler_raw_type & bcFluid, solid_bchandler_raw_type & bcSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling=false);

    /**\small evaluates the linear residual
     \param sol: solution vector
     \param rhs: right-hand side
     \param res: the output residual
     \param diagonalScaling: flag stating wether to perform diagonal scaling
*/
    void evalResidual( const vector_type& sol, const vector_ptrtype& rhs,  vector_type& res, bool diagonalScaling=false);

    /**\small evaluates the nonlinear residual
     \param bcFluid: fluid BCHandler
     \param bcSolid: solid BCHandler
     \param sol: solution vector
     \param rhs: right-hand side
     \param res: the output residual
     \param diagonalScaling: flag stating wether to perform diagonal scaling
     \param prec: preconditioner P, the final residual will be r=PAx-Pb
*/
    void evalResidual( fluid_bchandler_raw_type& bchFluid, solid_bchandler_raw_type& bchSolid, const vector_type& sol, vector_ptrtype& rhs, vector_type& res, bool diagonalScaling, matrix_ptrtype preconditioner);

    //!\small says if the preconditioner will be recomputed
    bool recomputePrec(){return(!M_reusePrec || M_resetPrec);}

#if OBSOLETE
    void setOperator(Epetra_Operator& epetraOperator){M_linearSolver->setOperator(epetraOperator);}
#endif

    //!\small inserts the input matrix in the linear solver
    void setMatrix(matrix_type& matr){M_linearSolver->setMatrix(matr);}

    //!\small left-multiply both the monolithic matrix and the rhs times a preconditioner matrix
    void applyPreconditioner(matrix_ptrtype robinCoupling, vector_ptrtype& rhs);
    //!\small left-multiply the monolithic matrix, the rhs and the preconditioner times a preconditioner matrix
    void applyPreconditioner( matrix_ptrtype robinCoupling, matrix_ptrtype& prec );
    //void setAztecooPreconditioner(const GetPot& dataFile, const std::string& section){M_linearSolver->setAztecooPreconditioner( dataFile, section);}

    //!\small adds the constant part to the monolithic matrix
    void updateMatrix(matrix_type & bigMatrixStokes);

    //!\small updates the rhs of the solid block.
    void updateSolidSystem(vector_ptrtype& rhsFluidCoupling);

    /**\small scales matrix and rhs
    \param rhs: the output rhs
    \param matrFull: the output matrix*/
    void    diagonalScale(vector_type& rhs, matrix_ptrtype matrFull);

    void shiftSolution(){}

    void couplingVariableExtrap(vector_ptrtype& lambda, vector_ptrtype& /*lambdaDot*/, bool& /*firstIter*/)
{    displayer().leaderPrint("norm( solution ) init = ", lambda->NormInf() );}

virtual  const vector_type& meshDisp()const{return this->M_meshMotion->disp();}
virtual  const vector_type&  veloFluidMesh() const;
virtual vector_type& veloFluidMesh();

protected:

	//void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);
	void variablesInit(const std::string& dOrder);

    boost::shared_ptr<EpetraMap>                      M_monolithicMap;
    matrix_ptrtype                                    M_couplingMatrix;
    UInt                                              M_interface;
    EpetraMap                                         M_interfaceMap;///the solid interface map
    boost::shared_ptr<vector_type>                    M_beta;
    matrix_ptrtype                                    M_monolithicMatrix;
    matrix_ptrtype                                    M_bigPrecPtr;
    UInt                                              M_DDBlockPrec;
    boost::shared_ptr<vector_type>                    M_rhsFull;

private:
    int                                               M_updateEvery;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
#if OBSOLETE
    boost::shared_ptr<vector_type>                    M_rhsShapeDerivatives;
#endif
    boost::shared_ptr<vector_type>                    M_rhsNew;
    bool                                              M_fullMonolithic;

//! perturbation value added to a diagonal zero block in the preconditioner
    Real                                              M_entry;
protected:
    matrix_ptrtype                                    M_robinCoupling;
    //    matrix_ptrtype                                    M_robinCouplingInv;
//! coefficient of the Robin preconditioner
    Real                                              M_alphaf;
//! coefficient of the Robin preconditioner
    Real                                              M_alphas;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
private:
    bool                                              M_diagonalScale;
    bool                                              M_reusePrec;
    bool                                              M_resetPrec;
    UInt                                              M_maxIterSolver;
    boost::shared_ptr<solver_type>                    M_linearSolver;
    matrix_ptrtype                                    M_fluidBlock;
    matrix_ptrtype                                    M_solidBlock;
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
    //M_monolithicMatrix->spy("j");
    //prec->spy("p");

    chrono.stop();
    M_solid->getDisplayer().leaderPrintMax("done in ", chrono.diff() );

    M_linearSolver->setMatrix(*M_monolithicMatrix);
    int numIter = M_linearSolver->solveSystem( rhs, step, prec, (M_reusePrec)&&(!M_resetPrec));

    if (numIter < 0)
    {
        chrono.start();

        M_solid->getDisplayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

        if (numIter <= -M_maxIterSolver)
           M_solid->getDisplayer().leaderPrint("   ERROR: Iterative solver failed twice.\n");
    }

    M_solid->getDisplayer().leaderPrint("   system solved.\n ");
}
}
#endif
