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

/*!
 *  @file
 *  @brief File containing the Monolithic Geometry--Implicit FSI Solver
 *
 *  @date 18-09-2008
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *
 *  @maintainer Paolo Crosetto <paolo.crosetto@epfl.ch>
 */

#ifndef _MONOLITHICGI_HPP
#define _MONOLITHICGI_HPP

#include <lifev/fsi/solver/FSIMonolithic.hpp>
#include <lifev/core/array/MatrixBlockMonolithicEpetra.hpp>
#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>


namespace LifeV
{
#ifdef OBSOLETE
class Epetra_FullMonolithic;
#endif

typedef FactorySingleton< Factory< FSIOperator, std::string > > FSIFactory_Type;

//! FSIMonolithicGI Geometry-Implicit solver
/*!
 *  @author Paolo Crosetto <paolo.crosetto@epfl.ch>
 *  Class handling the nonlinear monolithic solver for FSI problems. The (exact or inexact)
 *  Newton algorithm is used to solve the nonlinearity.
 *  The block structure of the jacobian matrix is
 *  \f$\left(\begin{array}{ccc} C&B&S\\ D&N&0\\ 0&E&H \end{array}\right)\f$
 *  where \f$N\f$ represents the solid block, \f$C\f$ the fluid block, \f$H\f$ is the harmonic extension block,
 *  while the extra
 *  diagonal blocks represent the coupling. The implementation of the stress continuity coupling condition
 *  is obtained by means of an augmented formulation.
 *  Different possible preconditioners are implemented.
 *
 *  Important parameters to set properly in the data file:
 *  - useShapeDerivatives: if true the shape derivatives block is added to the Jacobian matrix;
 *  - domainVelImplicit: if true the domain velocity w in the convective term is considered an unknown (at the time n+1);
 *  - convectiveTermDer: false if the convective term is linearized (\f$u^{n+1}\nabla(u^n-w^n)\f$),
 *  otherwise it can be either true (if we use the Newton method to solve the convective term nonlinearity) or false
 *  (fixed-point method);
 *  - semiImplicit:  if true only one iteration of the nonlinear solver is performed. Otherwise
 *  the nonlinear iterations continue up to the specified tolerance. Set it to true for the GCE;
 *  - method: can be either monolithicGE, monolithicGI if the geometry is treated respectively explicitly or implicitly,
 *  or exactJacobians, fixedPoint for partitioned strategies;
 *  - blockOper: specifies the matrix type to be used for the linear system: if AdditiveSchwarz, the matrix is the standard
 *  ine for GE; if AdditiveSchwarzRN the coupling blocks are of Robin type instead of Dirichlet and Neumann. The parameters
 *  for the Robin coupling are alphaf and alphas in the data file. NOTE: this method has currently been tested only for
 *  alphas=0.
 *  - DDBlockPrec: specifies the possible preconditioners to use. Can be: AdditiveSchwarz, MonolithicBlockComposedDN, MonolithicBlockComposedDN2,
 *  MonolithicBlockComposedNN, MonolithicBlockComposedDNND.
 */

class FSIMonolithicGI : public FSIMonolithic
{
public:

    typedef FSIMonolithic super_Type;
    typedef Preconditioner prec_Type;
    typedef boost::shared_ptr< prec_Type > prec_type;

    //!@name Constructor and Destructor
    //@{

    //! Empty Constructor
    FSIMonolithicGI();

    //! Destructor
    virtual ~FSIMonolithicGI() {}

    //@}

    //!@name Public Methods
    //@{

    /**
     Sets the parameters read from data file
     */
    void setup ( const GetPot& dataFile );

    //! initializes the fluid and mesh problems, creates the map of the global matrix
    void setupFluidSolid ( UInt const fluxes );

    //! builds the constant part of the fluid-structure-mesh motion matrix
    void buildSystem();

    /**
     evaluates the residual b-Ax
     \param res: output
     \param _sol: fluid domain displacement solution
     \param iter: current NonLinearRichardson (block Gauss Seidel for the tangent system) iteration
     */
    void evalResidual ( vector_Type& res, const vector_Type& sol, const UInt iter );

    //!Apply the boundary conditions to each block composing the monolithic problem
    /**
     Sets the vectors of: boundary conditions, FESpaces, couplings, offsets, and sets the blocks in the composed operator
     which constitutes the monolithic problem. Then calls the applyBoundaryConditions of the MonolithicBlockMatrix operator, passing
     also the right hand side.
     */
    void applyBoundaryConditions();

    void updateSolution ( const vector_Type& solution )
    {
        super_Type::updateSolution ( solution );

        //The size of the vectors for the ALE is = dimension of the ALE problem
        //To do the shift right we first need to extract the fluid displacement
        //And then push it into the ALE timeAdvance class.
        vectorPtr_Type displacementToSave ( new vector_Type (M_mmFESpace->map() ) );
        UInt offset ( M_solidAndFluidDim + nDimensions * M_interface );
        displacementToSave->subset (solution, offset);

        //This updateRHSFirstDerivative has to be done before the shiftRight
        //In fact it updates the right hand side of the velocity using the
        //previous times. The method velocity() uses it and then, the compuation
        //of the velocity is done using the current time and the previous times.
        M_ALETimeAdvance->updateRHSFirstDerivative ( M_data->dataFluid()->dataTime()->timeStep() );
        M_ALETimeAdvance->shiftRight ( *displacementToSave );
    }

    //! Set vectors for restart
    /*!
     *  Set vectors for restart
     */
    void setALEVectorInStencil (const vectorPtr_Type& fluidDisp,
                                const UInt iter,
                                const bool lastVector);

    //!@name Get Methods
    //@{

    //! get the current solution vector.
    LIFEV_DEPRECATED ( const vector_Type& solution() const )
    {
        if ( M_epetraWorldComm->MyPID() == 0 )
        {
            std::cerr << std::endl << "Warning: FSIMonolithic::solution() is deprecated!" << std::endl
                      << "         You should not access the solution inside FSIOperator or FSIMonolithic!" << std::endl;
        }

        return  M_fluidTimeAdvance->singleElement (0);
    }

    //! getter for the map of fluid-structure-interface (without the mesh motion)
    const MapEpetra& mapWithoutMesh() const
    {
        return *M_mapWithoutMesh;
    }

    //! getter for the global matrix of the system
    const matrixPtr_Type matrixPtr() const
    {
        return M_monolithicMatrix->matrix();
    }

    static bool S_register;

    //@}

protected:

    //!@name Protected Methods
    //@{

    //! set the block preconditioner
    void setupBlockPrec();

    //@}

private:

    //! @name Private Methods
    //@{

    //! Factory method for the system matrix, of type MonolithicBlockBase
    void createOperator ( std::string& operType )
    {
        M_monolithicMatrix.reset ( MonolithicBlockMatrix::Factory_Type::instance().createObject ( operType ) );
        M_monolithicMatrix.reset ( MonolithicBlockMatrix::Factory_Type::instance().createObject ( operType ) );
    }

    /**
     calculates the terms due to the shape derivatives given the mesh increment deltaDisp. The shape derivative block is assembled in a matrix
     (not in a right hand side representing the matrix-vector multiplication)
     \param sdMatrix: output. Shape derivatives block to be summed to the Jacobian matrix.
     */
    void shapeDerivatives ( FSIOperator::fluidPtr_Type::value_type::matrixPtr_Type sdMatrix );

    //! assembles the mesh motion matrix.
    /*!In Particular it diagonalize the part of the matrix corresponding to the
     Dirichlet condition expressing the coupling
     \param iter: current iteration: used as flag to distinguish the first nonlinear iteration from the others
     */
    void assembleMeshBlock ( UInt iter );

    //@}


    //!@name Private Members
    //@{

    boost::shared_ptr<MapEpetra>         M_mapWithoutMesh;
    //This vector is used in the shapeDerivatives method since a
    //copy of the solution at the current iteration k is necessary
    vectorPtr_Type                       M_uk;
    UInt                                 M_interface;
    matrixPtr_Type                       M_meshBlock;
    FSIOperator::fluidPtr_Type::value_type::matrixPtr_Type M_shapeDerivativesBlock;
    matrixPtr_Type                       M_solidDerBlock;
    //std::vector<fluidBchandlerPtr_Type>    M_BChsLin;
    //@}

    //! Factory method
    static FSIOperator* instantiate()
    {
        return new FSIMonolithicGI();
    }

};

//! Factory create function
inline FSIMonolithic* createFSIMonolithicGI()
{
    return new FSIMonolithicGI();
}

}
#endif
