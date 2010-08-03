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
   \file monolithicGI.hpp
   \author crosetto <Paolo Crosetto>
   \date 18/09/2008
*/
#ifndef _MONOLITHICGI_HPP
#define _MONOLITHICGI_HPP

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV
{
#ifdef OBSOLETE
class Epetra_FullMonolithic;
#endif

/**
   Monolithic Geomitry-Implicit solver
 * Class handling the nonlinear monolithic solver for FSI problems. The (exact or inexact)
 * Newton algorithm is used to solve the nonlinearity.
 * The block structure of the jacobian matrix is
 *\f$\left(\begin{array}{ccc}
 C&B&S\\
 D&N&0\\
 0&E&H
 \end{array}\right)\f$
 * where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, \f$H\f$ is the harmonic extension block,
 * while the extra
 * diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 * is obtained by means of an augmented formulation.
 * Different possible preconditioners are implemented.


 Important parameters to set properly in the data file:
 - useShapeDerivatives: MUST be false, because in the GE approach the geometry is explicit;
 - domainVelImplicit: MUST be false, because in the GE approach the geometry is explicit;
 - convectiveTermDer: false if the convective term is linearized (\f$u^{n+1}\nabla(u^n-w^n)\f$),
 otherwise it can be either true (if we use the Newton method to solve the convective term nonlinearity) or false
 (fixed-point method). For the GCE must be false;
 - semiImplicit:  if true only one iteration of the nonlinear solver is performed. Otherwise
 the nonlinear iterations continue up to the specified tolerance. Set it to true for the GCE;
 - method: can be either monolithicGE, monolithicGI if the geometry is treated respectively explicitly or implicitly,
 or exactJacobians, fixedPoint for partitioned strategies;
 - blockOper: specifies the matrix type to be used for the linear system: if AdditiveSchwarz, the matrix is the standard
 ine for GE; if AdditiveSchwarzRN the coupling blocks are of Robin type instead of Dirichlet and Neumann. The parameters
 for the Robin coupling are alphaf and alphas in the data file. NOTE: this method has currently been tested only for
 alphas=0.
 - DDBlockPrec: specifies the possible preconditioners to use. Can be: AdditiveSchwarz, ComposedDN, ComposedDN2,
 ComposedNN, ComposedDNND.

 */
class MonolithicGI : public Monolithic
{
public:

    typedef Monolithic                                         super;
    typedef EpetraPreconditioner                               prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>                   prec_type;

    ///constructor
    MonolithicGI();


    ///destructor
    virtual ~MonolithicGI(){}

    /**
       constructs the matrix handling the coupling and sums it to matrix
       \param matrix: output matrix
       \param coupling:   flag handling the coupling of the different blocks. Chosing properly it's value, ranging from 0 to 31, one can decide which couplings to keep and which to neglect (it works like the chmod command in bash)
    */
    void                        setUp( const GetPot& dataFile );

    //! initializes the fluid and mesh problems, creates the map of the global matrix
    void                       setupFluidSolid();

    /**
       updates the meshmotion, advances of a time step
       \param displacement: solution
    */
    void                        updateSystem();

    //! temporary methods (to be re-implemented whithin the BlockInterface framework):
    //@{
    /**
       builds the coupling
    */
    void                      couplingMatrix(matrix_ptrtype & matrix, const std::map<ID, ID>& locDofMap,
                                             int coupling=31);

    /**
       builds the coupling
    */
    void                      couplingMatrix2(matrix_ptrtype & matrix, const std::map<ID, ID>& locDofMap,
                                             int coupling=31);


    //!sets the block preconditioner
    int                        setupBlockPrec(vector_type& rhs);

    //@}

    //! builds the constant part of the fluid-structure-mesh motion matrix
    void                        buildSystem ();

    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: fluid domain displacement solution
       \param iter: current nonLinRichardson (block Gauss Seidel for the tangent system) iteration
    */
    void                      evalResidual(vector_type&        res,
                                           const vector_type& _sol,
                                           const UInt          _iter);




    /**
       solves the Jacobian system
       \param _muk: output, solution at the current block GS step
       \param _res: linear residual
       \param _linearRelTol: not used
    */
    void                        solveJac(vector_type       &_muk,
                                         const vector_type &_res,
                                         const Real       _linearRelTol);

    //! @getters
    //!{
    //! getter for the map of fluid-structure-interface (without the mesh motion)
    const EpetraMap&            mapWithoutMesh() const {return *M_mapWithoutMesh;}

    //! getter for the global matrix of the system
    const matrix_ptrtype        getMatrixPtr() const {return this->M_monolithicMatrix->getMatrix();}

    //! getter for the current iteration solution
    const vector_ptrtype        uk()      const {return M_uk;}

    //! getter for the domain displacement at the previous time step (to correctly visualize the solition of both GE
    //! and GI)
    const vector_type&          meshDisp()const
    {
        return this->M_meshMotion->dispOld();
    }

//     //! getter for the current iteration solution
//     void getSolution                  (vector_ptrtype& sol){sol = M_uk;}
    //}

    //! initialize the system with functions
    void                        initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                                            FSIOperator::solid_type::value_type::Function const& p0,
                                            FSIOperator::solid_type::value_type::Function const& d0,
                                            FSIOperator::solid_type::value_type::Function const& w0,
                                            FSIOperator::solid_type::value_type::Function const& df0 );

private:

    //! creates the
    void createOperator( std::string& operType )
    {
        M_monolithicMatrix.reset(BlockMatrix::Factory::instance().createObject( operType ));
    }

    void                        initialize( vector_type const& u0, vector_type const& p0, vector_type const& d0, vector_type const& df0);

    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    void                        shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol,  bool fullImplicit, bool convectiveTerm);

    void applyPreconditioner( const matrix_ptrtype prec, matrix_ptrtype& oper );

    void updateMatrix(matrix_type & bigMatrixStokes);

    template <typename SolverType, typename PrecOperatorPtr>
    void iterateFullMonolithic(const vector_type& rhs, vector_type& step, PrecOperatorPtr prec, SolverType linearSolver);

    void assembleFluidBlock(UInt iter);

    void assembleMeshBlock(UInt iter);



    /**
       \small builds the matrix [1,0,0,0;0,1,0,alphaf;0,0,1,0;0,alphas,0,1] to IdentityMatrix. This matrix is intended to be a preconditioner and has not been tested yet.
       \param IdentityMatrix: a matrix, if it is an identity the method builds the Robin-Robin linear combination of boundary conditions
       \param alphaf: coefficient
       \param alphas: coefficient
       \param inverse if true the inverse of the linear combination is computed. the output matrix is [];
    */
    void robinCoupling(matrix_ptrtype& IdentityMatrix, const Real& alphaf, const Real& alphas, int coupling = 4);

    void setMatrix(matrix_type& matr){M_linearSolver->setMatrix(matr);}

    boost::shared_ptr<EpetraMap>         M_mapWithoutMesh;
    vector_ptrtype                       M_uk;
    vector_ptrtype                       M_meshVel;
    bool                                 M_domainVelImplicit;
    bool                                 M_convectiveTermDer;
    IfpackComposedPrec::operator_type    M_solidOper;
    IfpackComposedPrec::operator_type    M_fluidOper;
    IfpackComposedPrec::operator_type    M_meshOper;
    UInt                                 M_interface;
    prec_type                            M_compPrecPtr;
    matrix_ptrtype                       M_precMatrPtr;
    matrix_ptrtype                       M_meshBlock;
    boost::shared_ptr<EpetraMap>         M_monolithicInterfaceMap;
    static bool reg;
};

template <typename SolverType, typename PrecOperatorPtr>
void MonolithicGI::
iterateFullMonolithic(const vector_type& rhs, vector_type& step, PrecOperatorPtr prec, SolverType linearSolver)
{
    M_solid->getDisplayer().leaderPrint("  M-  Preconditioner type:                     ", M_data->DDBlockPreconditioner(), "\n" );
    Chrono chrono;

    M_monolithicMatrix->getMatrix()->GlobalAssemble();
    //necessary if we did not imposed Dirichlet b.c.
    M_linearSolver->setMatrix(*M_monolithicMatrix->getMatrix());

    M_linearSolver->setReusePreconditioner( (M_reusePrec) && (!M_resetPrec) );
    int numIter = M_linearSolver->solveSystem( rhs, step, boost::static_pointer_cast< Epetra_Operator >(prec) );
    //int numIter = M_linearSolver->solveSystem( *const_cast<const vector_type*>(&rhs), step,  );

    if (numIter < 0)
        {
            chrono.start();

            M_solid->getDisplayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

            if (numIter <= -M_maxIterSolver)
                M_solid->getDisplayer().leaderPrint("   ERROR: Iterative solver failed.\n");
        }

    M_solid->getDisplayer().leaderPrint("  M-  Jacobian system solved                   \n");
}


}
#endif

