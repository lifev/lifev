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

/*!
    @file
    @brief Monolithic Geometry--Explicit FSI Solver

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 26 Jul 2010

    \include ../../testsuite/test_monolithic/fluidstructure.dox
    This file implements the Monolithic Geometry--Explicit solver, see \ref CDFQ for details
 */

#ifndef MONOLITHICGE_H
#define MONOLITHICGE_H 1

#include <lifemc/lifesolver/BlockMatrix.hpp>
#include <lifemc/lifesolver/BlockMatrixRN.hpp>
#include <lifemc/lifesolver/ComposedDN.hpp>
#include <lifemc/lifesolver/ComposedNN.hpp>
#include <lifemc/lifesolver/ComposedDNND.hpp>

#include <life/lifecore/life.hpp>

#include <lifemc/lifesolver/Monolithic.hpp>

namespace LifeV
{

//! monolithicGE - Monolithic Geometry-Explicit solver
/*!
  \include ../../testsuite/test_monolithic/fluidstructure.dox
  @author Paolo Crosetto
  @see \ref CDFQ


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
class MonolithicGE : public Monolithic
{
public:

    typedef Monolithic super_Type;
    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    MonolithicGE():
        super_Type()
    {}

    //! Destructor
    ~MonolithicGE()
    {}

    //@}



    //! @name Public Methods
    //@{

    //! Setup method for the subfroblem
    /**
       Sets up the fluid, solid and harmonic extension finite element spaces and initializes most of the variables
       used in the solver
     */
    void setupFluidSolid( UInt const fluxes );

    //! setup of the dofs
    /** calls super_Type::setupDof and instantiate the boundary condition vector needed to couple fluid--structure and harmonic extention*/
    void setupDOF();

    //! setUp from data file
    /**
       calls the setup for the fluid, solid and mesh motion problems
     */
    void setupSystem();

    //!Updates the system for the next time step
    /**
       Calls the updateSystem of the mother class and updates the solid displacement in the solid problem
     */
    void updateSystem();


    /**
       evaluates the residual b-Ax
       \param res: output
       \param _sol: monolithic solution
       \param iter: current nonLinRichardson (Newton) iteration
    */
    void   evalResidual(vector_Type&        res,
                        const vector_Type& _sol,
                        const UInt          _iter);

    /**
      iterates the mesh
      \param disp: monolithic solution
    */
    void iterateMesh( const vector_Type& disp );

    //@}
    //!@name Getter Methods
    //@{

    //! Gets the solution
    const vector_Type& getSolution() const { return *M_un; }

    //! Gets the solution ptr
    vectorPtr_Type& solutionPtr() { return M_un; }

    //! Sets the solution
    void setSolution( const vector_Type& solution ) { M_un.reset( new vector_Type( solution ) ); }
    //@}

    //!@name Setter Methods
    //@{
    //! Sets the solution ptr
    void setSolutionPtr( const vectorPtr_Type& sol) { M_un = sol; }

    //! Applies the bounsary conditions to the matrix
    void applyBoundaryConditions();
    //@}

    //! Factory method
    static FSIOperator* createM(){ return new MonolithicGE(); }

private:
    //!@name Private Methods
    //@{
    void createOperator( std::string& operType )
    {
        M_monolithicMatrix.reset(BlockMatrix::Factory::instance().createObject( operType ));
    }
    //@}

    //!@name Private Members
    static bool reg;
    //@}
};

} // Namespace LifeV

#endif /* MONOLITHICGCE_H */
