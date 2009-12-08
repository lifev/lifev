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
 * Class handling the monolithic solver for FSI problems. The block structure of the matrix is
 *$\left(\begin{array}{cc}
 *C&B
 *D&N
 *\end{array}\right)$
 */
#ifndef _MONOLITHIC_HPP
#define _MONOLITHIC_HPP


#include <life/lifesolver/FSIOperator.hpp>

#include <lifemc/lifealg/IfpackComposedPrec.hpp>
#include <lifemc/lifealg/ComposedPreconditioner.hpp>
//#include <Epetra_IntVector.h>

namespace LifeV
{
/**
 * Class handling the monolithic solver for FSI problems. The block structure of the matrix is
 *\f$\left(\begin{array}{cc}
 C&B\\
 D&N
 \end{array}\right)\f$
 * where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, while the extra
 * diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 * is obtained by means of an augmented formulation.
 * Different possible preconditioners are implemented.
 */

class WRONG_PREC_EXCEPTION;

class Monolithic : public FSIOperator
{
public:


    typedef FSIOperator                            super;
    typedef FSIOperator::fluid_type::value_type::matrix_type   matrix_type;
    //typedef IfpackComposedPrec                               prec_raw_type;
    typedef EpetraPreconditioner                               prec_raw_type;
    typedef boost::shared_ptr<matrix_type>                     matrix_ptrtype;
    typedef boost::shared_ptr<prec_raw_type>                   prec_type;
    typedef SolverTrilinos                                     solver_type;


    // constructors

//!@name constructors, destructor
//!@{
    Monolithic();
    ~Monolithic();
//!@}

//!@name setup methods
//!@{

    /**
       create FEspace
    */
    void setupFEspace();

    /**
       sets the interface map between the fluid and solid meshes
    */
    void setupDOF( void );

    //!@}

    //!@name Methods
    //!@{
    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */
    void   evalResidual(vector_type&        res,
                        const vector_type& _sol,
                        const UInt          _iter);




    /** sets the parameters from the data file*/
    void setDataFromGetPot( GetPot const& data );



    /**
       restricts a vector with a monolithic map on the solid interface map
       \param _lambdasolid: vector on the solid interface
       \param _disp: monolithic vector
    */
    void monolithicToInterface(    vector_type& _lambdaSolid, const vector_type& _sol) ;

    /**
       sets the vector M_solid->dispSolid() to the monolithic solution M_solid->disp() in the solid nodes, to 0 in the fluid nodes
    */
    void setDispSolid(const vector_type& sol);


    /**
       restricts a vector with a monolithic map on another map that must
       have a sequential numbering (not the interface map)
       \param _dispFluid: vector on the fluid domain
       \param _disp: monolithic vector
    */
    void monolithicToX(const vector_type& _disp, vector_type& _dispFluid, EpetraMap& map, UInt offset=(UInt)0);

    /**
       iterates the mesh
       \param disp: monolithic solution
    */
    void iterateMesh(const vector_type& disp);

    /**
       builds the constant part of the monolithic matrix
    */
    void buildSystem();


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
    void shapeDerivatives(matrix_ptrtype sdMatrix, const vector_type& sol,  bool fullImplicit, bool convectiveTerm);
    //! getters

    void setupSystem( );
    virtual void setUp( const GetPot& dataFile );

    /**
       \small adds a constant scalar entry to a diagonal block of a matrix
       \param entry: entry
       \param matrix: matrix
       \param Map: map specifying the location where to add diagonal entries
       \param offset: offset from which starts the insertion
       \param replace=false: flag stating wether the diagonal elements already exist or have to be inserted in the map
    */
    void addDiagonalEntries(Real entry, matrix_ptrtype matrix, const EpetraMap& Map, UInt offset=0, bool replacs=false);



    /**
       \small initialize the current solution vector. Note: this is not sufficient for the correct initialization
       of bdf!
    */
    void initialize( vector_ptrtype u0)
    {M_un=u0;}


#ifdef HAVE_TRILINOS_ANASAZI
    /**
       \small Computes the maximum singular value of the preconditioned system \f$P^-1A\f$ where \f$P\f$ is an
       instance of ComposedPreconditioner and \f$A\f$ is the system matrix.
    */
    void computeMaxSingularValue();
#endif
    /**
       \small Computes the wall shear stress. Some issues when working in parallel still has to be fixed.
    */
    vector_ptrtype computeWS();

    /**
       \small Enables the computation of the wall shear stress on the specified boundary.
       \param flag : flag specifying the boundary of interest.
    */
    void enableWssComputation(EntityFlag flag);

    //!@}
    //!@name Getters
    //!@{

    //! Returns the wall stress
    //!@note still not fixed in parallel
    vector_ptrtype getWS( )
    {
        return M_wss;
    }

    //! returns a boost shared pointer to the preconditioner
    prec_raw_type & getPrec(){return *M_precPtr;}

#ifdef OBSOLETE
    /** get the shape derivatives vector*/
    vector_type getRhsShapeDerivatives(){return *M_rhsShapeDerivatives;}
#endif
    //    const boost::shared_ptr<EpetraMap>& monolithicMap() {return M_monolithicMap;}

    //!get the total dimension of the FS interface
    const UInt getDimInterface() const {return nDimensions*M_interface ;}

    //! Returns the solution at the previous time step
    vector_ptrtype const& un(){return M_un;}

    //! Returns true if CE of FI methods are used, false otherwise (GCE)
    bool const isFullMonolithic(){return M_fullMonolithic;}

    /** returns the monolithic map*/
    virtual    boost::shared_ptr<EpetraMap>& getCouplingVariableMap(){return M_monolithicMap;}

    //!@}
    //!@name Virtual methods
    //!@{
    /**
       assigns each mesh partition to the corresponding processor, builds the monolithic map
    */
    virtual void setupFluidSolid();

    /**
       solves the Jacobian system
       \param _muk: output, solution at the current Newton step
       \param _res: nonlinear residual
       \param _linearRelTol: not used

       \small The preconditioner type is usually an algebraic additive Schwarz. The following values
       assigned to the field DDBlockPrec in the data file correspond to different variants:

       - DDBlockPrec = 0 is AAS on a the system matrix. Can give problem in parallel, we suggest to use DDBlockPrec = 2 instead
       - DDBlockPrec = 2 is AAS on a the system matrix, where the system is left-premultiplied times a preconditioner
       that mix the interface boundary conditions
       - DDBlockPrec = 3 only for testing purposes
       - DDBlockPrec = 4 only for testing purposes
       - DDBlockPrec = 5 no preconditioner is set (if you use the native Aztecoo preconditioners set this option)

       Only for the Monolithic Geometry-Convective Explicit:
       - DDBlockPrec = 1 is AAS on a Dirichlet-Neumann preconditioner
       - DDBlockPrec = 7 is AAS on a Dirichlet-Neumann preconditioner using the ComposedPreconditioner strategy
       - DDBlockPrec = 8 is AAS on an alternative Dirichlet-Neumann preconditioner using the ComposedPreconditioner strategy

       Only for the fullMonolithic Convective Explicit:
       - DDBlockPrec = 6 is AAS on the quasi-newton matrix
       - DDBlockPrec = 9 is AAS on the quasi-newton matrix obtained with the ComposedPreconditioner strategy
       - DDBlockPrec = 10 is AAS on an alternative matrix obtained with the ComposedPreconditioner strategy
       - DDBlockPrec = 11 is AAS on the quasi-newton matrix obtained with the ComposedPreconditioner strategy, composing
       3 preconditioners
       - DDBlockPrec = 12 is AAS on an alternative matrix obtained with the ComposedPreconditioner strategy, composing
       3 preconditioners
    */
    virtual void   solveJac(vector_type&       _muk,
                            const vector_type& _res,
                            const Real       _linearRelTol);

    /**
       updates the meshmotion, advances of a time step
       \param _sol: solution
    */
    virtual void updateSystem(const vector_type& _sol);


    /**
       \small initialize with functions
    */
    virtual void initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                             FSIOperator::solid_type::value_type::Function const& p0,
                             FSIOperator::solid_type::value_type::Function const& d0,
                             FSIOperator::solid_type::value_type::Function const& w0,
                             FSIOperator::solid_type::value_type::Function const& df0);


    virtual void initializeMesh(vector_ptrtype fluid_dispOld);
    //!@}

protected:

    //!@name Protected methods
    //!@{
    //!@name Templated methods
    //!@{
    template <typename SolverType, typename PrecOperatorPtr>
    /**
       \small solves the monolithic system, once a solver, a preconditioner and a rhs have been defined.
    */
    void iterateMonolithic(const vector_type& rhs, vector_type& step, PrecOperatorPtr prec, SolverType linearSolver);
    //!@}

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
    virtual void couplingMatrix(matrix_ptrtype & bigMatrix, int couplingSolid = 31);

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
       \small Method to replace a block in the matrix with a zero block. It is inefficient and works only serially. It is implemented only for testing different preconditioners.
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

    void couplingVariableExtrap(vector_ptrtype& /*lambda*/, vector_ptrtype& /*lambdaDot*/, bool& /*firstIter*/)
    { }
    //void solidInit(const RefFE* refFE_struct, const LifeV::QuadRule* bdQr_struct, const LifeV::QuadRule* qR_struct);
    void solidInit(std::string const& dOrder);

	//void variablesInit(const RefFE* refFE_struct,const LifeV::QuadRule*  bdQr_struct, const LifeV::QuadRule* qR_struct);

    void variablesInit(std::string const& dOrder);

    //!@}
    //!@name protected setters
    //!@{
#ifdef OBSOLETE
    void setOperator(Epetra_Operator& epetraOperator){M_linearSolver->setOperator(epetraOperator);}
#endif

    //!\small inserts the input matrix in the linear solver
    void setMatrix(matrix_type& matr){M_linearSolver->setMatrix(matr);}

    void setFluxBC             (fluid_bchandler_type bc_fluid);

    void setRobinBC             (fluid_bchandler_type bc_solid);
    //!@}

    //!@name protected getters
    //!@{
    const fluid_bchandler_type& BCh_flux()                      const { return M_BCh_flux; }
    //!@}


    //!@name Protected attributes
    //@{
    boost::shared_ptr<EpetraMap>                      M_monolithicMap;
    matrix_ptrtype                                    M_couplingMatrix;
    UInt                                              M_interface;
    EpetraMap                                         M_interfaceMap;///the solid interface map
    EpetraMap                                         M_monolithicInterfaceMap;
    boost::shared_ptr<vector_type>                    M_beta;
    matrix_ptrtype                                    M_monolithicMatrix;
    matrix_ptrtype                                    M_precMatrPtr;
    prec_type                                         M_precPtr;
    UInt                                              M_DDBlockPrec;
    boost::shared_ptr<vector_type>                    M_rhsFull;

    fluid_bchandler_type                              M_BCh_flux;
    solid_bchandler_type                              M_BCh_Robin;
    UInt                                              M_fluxes;
    solid_bchandler_type                              M_BChWSS;
    matrix_ptrtype                                    M_bdMass;
    BCFunctionMixte                                   M_bcfWss;
    matrix_ptrtype                                    M_robinCoupling;
    //    matrix_ptrtype                                    M_robinCouplingInv;
    //! coefficient of the Robin preconditioner
    Real                                              M_alphaf;
    //! coefficient of the Robin preconditioner
    Real                                              M_alphas;
    UInt                                              M_offset;
    UInt                                              M_solidAndFluidDim;
    IfpackComposedPrec::operator_type                 M_solidOper;
    IfpackComposedPrec::operator_type                 M_fluidOper;
    IfpackComposedPrec::operator_type                 M_meshOper;
    matrix_ptrtype                                    M_fluidBlock;
    matrix_ptrtype                                    M_solidBlock;
    matrix_ptrtype                                    M_solidBlockPrec;
    matrix_ptrtype                                    M_meshBlock;
    boost::shared_ptr<solver_type>                    M_linearSolver;
    vector_ptrtype                                    M_wss;
    //!}
private:
    //!@name Private methods
    //!@{
    void initialize( vector_type const& u0, vector_type const& p0, vector_type const& d0);
    //!@}
    //!@name Private attributes
    //!@{
    boost::shared_ptr<ComposedPreconditioner<Epetra_Operator> > M_PAAP;
    //UInt                                              M_updateEvery;
    boost::shared_ptr<vector_type>                    M_numerationInterface;
#ifdef OBSOLETE
    boost::shared_ptr<vector_type>                    M_rhsShapeDerivatives;
#endif
    boost::shared_ptr<vector_type>                    M_rhsNew;
    bool                                              M_fullMonolithic;
    //! perturbation value added to a diagonal zero block in the preconditioner
    Real                                              M_entry;
    bool                                              M_diagonalScale;
    bool                                              M_reusePrec;
    bool                                              M_resetPrec;
    UInt                                              M_maxIterSolver;
    //!@}
};


template <typename SolverType, typename PrecOperatorPtr>
void Monolithic::
iterateMonolithic(const vector_type& rhs, vector_type& step, PrecOperatorPtr prec, SolverType linearSolver)
{
    M_solid->getDisplayer().leaderPrint("    preconditioner type : ", M_DDBlockPrec );
    Chrono chrono;

    M_solid->getDisplayer().leaderPrint("    Solving the system ... \n" );

    M_solid->getDisplayer().leaderPrint("    Updating the boundary conditions ... ");

    M_linearSolver->setMatrix(*M_monolithicMatrix);


    M_linearSolver->setReusePreconditioner( (M_reusePrec) && (!M_resetPrec) );
    int numIter = M_linearSolver->solveSystem( rhs, step, prec );

    if (numIter < 0)
        {
            chrono.start();

            M_solid->getDisplayer().leaderPrint("   Iterative solver failed, numiter = ", -numIter );

            if (numIter <= -M_maxIterSolver)
                M_solid->getDisplayer().leaderPrint("   ERROR: Iterative solver failed.\n");
        }

    M_solid->getDisplayer().leaderPrint("   system solved.\n ");
}

class WRONG_PREC_EXCEPTION{
public:
    WRONG_PREC_EXCEPTION(){}
    virtual ~WRONG_PREC_EXCEPTION(){}
};


}
#endif
