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

    This file implements the Monolithic Geometry--Explicit solver, see \cite CrosettoEtAl2009 for details
 */

#ifndef MONOLITHICGE_H
#define MONOLITHICGE_H 1

#include <lifev/core/LifeV.hpp>

#include <lifev/fsi/solver/MonolithicBlockMatrix.hpp>
#include <lifev/fsi/solver/MonolithicBlockMatrixRN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedNN.hpp>
#include <lifev/fsi/solver/MonolithicBlockComposedDNND.hpp>

#include <lifev/fsi/solver/FSIMonolithic.hpp>

namespace LifeV
{

//! FSIMonolithicGE - FSIMonolithic Geometry-Explicit solver
/*!
  @author Paolo Crosetto
  @see \cite CrosettoEtAl2009

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
 - DDBlockPrec: specifies the possible preconditioners to use. Can be: AdditiveSchwarz, MonolithicBlockComposedDN, MonolithicBlockComposedDN2,
 MonolithicBlockComposedNN, MonolithicBlockComposedDNND.
 */
typedef FactorySingleton<Factory<FSIOperator, std::string> >                    FSIFactory_Type;
class FSIMonolithicGE : public FSIMonolithic
{
public:

    typedef FSIMonolithic super_Type;

    //! @name Constructor & Destructor
    //@{

    //! Empty Constructor
    FSIMonolithicGE() :
        super_Type()
    {}

    //! Destructor
    ~FSIMonolithicGE()
    {}

    //@}


    //! @name Public Methods
    //@{

    //! Setup method for the subfroblem
    /**
       Sets up the fluid, solid and harmonic extension finite element spaces and initializes most of the variables
       used in the solver
     */
    void setupFluidSolid ( UInt const fluxes );

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

    //! Set vectors for restart
    /*!
     *  Set vectors for restart
     */
    void setALEVectorInStencil (const vectorPtr_Type& fluidDisp, const UInt iter, const bool /*lastVector*/);

    /**
       evaluates the residual Ax-b
       \param res: output
       \param _sol: monolithic solution
       \param iter: current NonLinearRichardson (Newton) iteration
    */
    void evalResidual ( vector_Type& res, const vector_Type& sol, const UInt iter );

    /**
      iterates the mesh
      \param disp: monolithic solution
    */
    void iterateMesh ( const vector_Type& disp );

    //! Applies the bounsary conditions to the matrix
    void applyBoundaryConditions();

    void updateSolution ( const vector_Type& solution );

    //@}


    //! @name Get Methods
    //@{

    //! Gets the solution
    LIFEV_DEPRECATED ( const vector_Type& solution() const )
    {
        if ( M_epetraWorldComm->MyPID() == 0 )
        {
            std::cerr << std::endl << "Warning: FSIMonolithic::solution() is deprecated!" << std::endl
                      << "         You should not access the solution inside FSIOperator or FSIMonolithic!" << std::endl;
        }

        return  M_fluidTimeAdvance->singleElement (0);
    }

    //@}


    //! Factory method
    static FSIOperator* instantiate()
    {
        return new FSIMonolithicGE();
    }

private:

    //!@name Private Methods
    //@{

    void createOperator ( std::string& operType )
    {
        M_monolithicMatrix.reset (MonolithicBlockMatrix::Factory_Type::instance().createObject ( operType ) );
    }

    //@}


public:

    static bool S_register;
};

//! Factory create function
inline FSIMonolithic* createFSIMonolithicGE()
{
    return new FSIMonolithicGE();
}

} // Namespace LifeV

#endif /* MONOLITHICGCE_H */
