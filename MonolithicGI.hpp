/* -*- mode: c++ -*- */
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

/**
   \file monolithicGI.hpp
   @breif Monolithic Geometry--Implicit FSI Solver
   \author crosetto <Paolo Crosetto>
   \date 18 Sep 2008

   \include ../../testsuite/test_monolithic/fluidstructure.dox
   This file implements the Monolithic Geometry--Implicit solver, see \ref CDFQ for details

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
 - useShapeDerivatives: if true the shape derivatives block is added to the Jacobian matrix;
 - domainVelImplicit: if true the domain velocity w in the convective term is considered an unknown (at the time n+1);
 - convectiveTermDer: false if the convective term is linearized (\f$u^{n+1}\nabla(u^n-w^n)\f$),
 otherwise it can be either true (if we use the Newton method to solve the convective term nonlinearity) or false
 (fixed-point method);
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

    //!@name Constructor and Destructor
    //@{
    //! Empty Constructor
    MonolithicGI();
    //! Destructor
    virtual ~MonolithicGI(){}
    //@}

    //!@name Public Methods
    //@{
    /**
       constructs the matrix handling the coupling and sums it to matrix
       \param matrix: output matrix
       \param coupling:   flag handling the coupling of the different blocks. Chosing properly it's value, ranging from 0 to 31, one can decide which couplings to keep and which to neglect (it works like the chmod command in bash)
    */
    void                        setUp( const GetPot& dataFile );

    //! initializes the fluid and mesh problems, creates the map of the global matrix
    void                       setupFluidSolid( UInt const fluxes );

    /**
       updates the meshmotion, advances of a time step
       \param displacement: solution
    */
    void                        updateSystem();


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



    void applyBoundaryConditions();


    /**
       solves the Jacobian system
       \param _muk: output, solution at the current block GS step
       \param _res: linear residual
       \param _linearRelTol: not used
    */
    void                        solveJac(vector_type       &_muk,
                                         const vector_type &_res,
                                         const Real       _linearRelTol);
    //! initialize the system with functions
    void                        initialize( FSIOperator::fluid_type::value_type::Function const& u0,
                                            FSIOperator::solid_type::value_type::Function const& p0,
                                            FSIOperator::solid_type::value_type::Function const& d0,
                                            FSIOperator::solid_type::value_type::Function const& w0,
                                            FSIOperator::solid_type::value_type::Function const& df0 );

    //@}
    //!@name Get Methods
    //@{

    //! getter for the map of fluid-structure-interface (without the mesh motion)
    const EpetraMap&            mapWithoutMesh() const {return *M_mapWithoutMesh;}

    //! getter for the global matrix of the system
    const matrixPtr_Type        getMatrixPtr() const {return this->M_monolithicMatrix->getMatrix();}

    //! getter for the current iteration solution
    const vector_ptrtype  uk()  const      {return M_uk;}

    //! getter for the domain displacement at the previous time step (to correctly visualize the solition of both GE
    //! and GI)
    const vector_type&          meshDisp()const
    {
        return this->M_meshMotion->dispOld();
    }

    //! get the solution.
    const vector_type& getSolution() const { return *M_uk; }

    //! get the solution.
    vector_ptrtype& solutionPtr() { return M_uk; }

    //@}
    //!@name Set Methods
    //@{
    //! set the solution
    void setSolution( const vector_type& solution ) { M_uk.reset( new vector_type( solution ) ); }

    void setSolutionPtr                     ( const vector_ptrtype& sol) { M_uk = sol; }
    //@}

protected:

    //!@name Protected Methods
    //@{
    //!sets the block preconditioner
    int                        setupBlockPrec( );
    //@}

private:

    //! @name Private Methods
    //@{

    //! Factory method for the system matrix, of type MonolithicBlockBase
    void createOperator( std::string& operType )
    {
        M_monolithicMatrix.reset(BlockMatrix::Factory::instance().createObject( operType ));
    }

    /**
       calculates the terms due to the shape derivatives on the rhs of the monolithic system rhsShapeDerivatives
       given the mesh increment deltaDisp.
       \param rhsShapeDerivatives: output. Shape derivative terms.
       \param meshDeltaDisp: input. Mesh displacement increment.
    */
    void shapeDerivatives(matrixPtr_Type sdMatrix, const vector_type& sol,  bool fullImplicit, bool convectiveTerm);

    //! assembles the mesh motion matrix.
    /*!In Particular it diagonalize the part of the matrix corresponding to the
      Dirichlet condition expressing the coupling
      \param iter: current iteration: used as flag to distinguish the first nonlinear iteration from the others
     */
    void assembleMeshBlock(UInt iter);

    //@}
    //!@name Private Members
    //@{

    boost::shared_ptr<EpetraMap>         M_mapWithoutMesh;
    vector_ptrtype                       M_uk;
    bool                                 M_domainVelImplicit;
    bool                                 M_convectiveTermDer;
    UInt                                 M_interface;
    matrixPtr_Type                       M_meshBlock;
    matrixPtr_Type                       M_shapeDerivativesBlock;
    matrixPtr_Type                       M_solidDerBlock;
    //std::vector<fluid_bchandler_type>    M_BChsLin;
    static bool                          reg;
    //@}
};

}
#endif

