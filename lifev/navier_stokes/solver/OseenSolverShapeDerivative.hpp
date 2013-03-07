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
    @brief This file contains an Oseen equation solver class with shape derivative
           for fluid structure interaction problem

    @author MiGlobaluel A. Fernandez <miGlobaluel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @date 09-06-2003

    @author G. Fourestey
    @date 00-02-2007

    @contributor Zhen Wang <zhen.wang@emory.edu>

 */


#ifndef OSEENSOLVERSHAPEDERIVATIVE_H
#define OSEENSOLVERSHAPEDERIVATIVE_H 1

#include <lifev/navier_stokes/solver/OseenSolver.hpp>

namespace LifeV
{

//! @class OseenSolverShapeDerivative
/*!

    @brief This class contains an Oseen equation solver class with shape derivative
           for fluid structure interaction problem

    @author MiGlobaluel A. Fernandez <miGlobaluel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @date 09-06-2003
    @author G. Fourestey
    @date 02-2007

    @contributor Zhen Wang <zhen.wang@emory.edu>

 */

template< typename MeshType, typename SolverType = LifeV::SolverAztecOO >
class OseenSolverShapeDerivative:
    public OseenSolver< MeshType, SolverType >
{

public:

    //! @name Public Types
    //@{

    typedef MeshType                                          mesh_Type;
    typedef SolverType                                        linearSolver_Type;
    typedef OseenSolver< mesh_Type, linearSolver_Type >       oseenSolver_Type;
    typedef typename oseenSolver_Type::vector_Type            vector_Type;
    typedef typename oseenSolver_Type::matrix_Type            matrix_Type;
    typedef typename oseenSolver_Type::matrixPtr_Type         matrixPtr_Type;
    typedef typename oseenSolver_Type::data_Type              data_Type;
    typedef typename oseenSolver_Type::preconditioner_Type    preconditioner_Type;
    typedef typename oseenSolver_Type::preconditionerPtr_Type preconditionerPtr_type;
    typedef typename oseenSolver_Type::bcHandler_Type         bcHandler_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    OseenSolverShapeDerivative();

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */
    OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                                 FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                                 FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                                 boost::shared_ptr<Epetra_Comm>& communicator,
                                 const Int                       lagrangeMultiplier = 0);

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param monolithicMap MapEpetra class
        @param offset
     */
    OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                                 FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                                 FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                                 boost::shared_ptr<Epetra_Comm>& communicator,
                                 const MapEpetra                 monolithicMap,
                                 const UInt                      offset = 0 );

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param mmFESpace FE space
        @param communicator MPI communicator
        @param monolithicMap MapEpetra class
        @param offset
     */
    OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                                 FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                                 FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                                 FESpace<mesh_Type, MapEpetra>&  mmFESpace,
                                 boost::shared_ptr<Epetra_Comm>& communicator,
                                 const MapEpetra                 monolithicMap,
                                 const UInt                      offset = 0 );

    //! Virtual destructor
    virtual ~OseenSolverShapeDerivative();

    //@}


    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    void setUp ( const GetPot& dataFile );

    //!
    /*!
        @param bcHandler BC handler
     */
    void solveLinearSystem ( bcHandler_Type& bcHandler );

    //! Update linear system.
    /*!
        @param matrixNoBC Fluid matrix withoud BC
        @param alpha alpha
        @param un Beta
        @param uk Fluid solution
        @param disp mesh_Type deltaX
        @param w mesh_Type Velocity
        @param dw mesh_Type deltaVelocity
        @param sourceVector RHS (usually 0 )
     */
    void updateLinearSystem ( const matrix_Type& matrixNoBC,
                              Real&              alpha,
                              const vector_Type& un,
                              const vector_Type& uk,
                              const vector_Type& disp,
                              const vector_Type& w,
                              const vector_Type& dw,
                              const vector_Type& sourceVector);

    //! Update shape derivatives.
    /*!
        @param matrixNoBC Fluid matrix withoud BC
        @param alpha alpha
        @param un Beta
        @param uk Fluid solution
        @param disp mesh_Type deltaX
        @param w mesh_Type Velocity
        @param offset
        @param dFESpace
        @param wImplicit
        @param convectiveTermDerivative
     */
    void updateShapeDerivatives ( matrix_Type&                   matrixNoBC,
                                  Real&                          alpha,
                                  const vector_Type&             un,
                                  const vector_Type&             uk,
                                  //const vector_Type&           disp,
                                  const vector_Type&             w,
                                  UInt                           offset,
                                  FESpace<mesh_Type, MapEpetra>& dFESpace,
                                  bool                           wImplicit = true,
                                  bool                           convectiveTermDerivative = false);

    //@}


    //! @name Set Methods
    //@{

    //! Set
    /*!
        @param rightHandSide
     */
    void updateLinearRightHandSideNoBC ( const vector_Type& rightHandSide)
    {
        M_linearRightHandSideNoBC = rightHandSide;
    }

    //@}

    //! @name Get Methods
    //@{

    //! Return
    /*!
        @return M_linearRightHandSideNoBC
     */
    vector_Type& linearRightHandSideNoBC()
    {
        return M_linearRightHandSideNoBC;
    }

    const vector_Type& linearRightHandSideNoBC() const
    {
        return M_linearRightHandSideNoBC;
    }

    //! Get the solution of the Shape Derivative problem.
    /*!
        @return vector containing the solution of the Shape Derivative problem.
     */
    const vector_Type& linearSolution() const
    {
        return M_linearSolution;
    }

    //! Return
    /*!
        @return stabilization
     */
    bool stabilization()
    {
        return this->M_stabilization;
    }

    const bool& stabilization() const
    {
        return this->M_stabilization;
    }

    //! Compute the derivative of the flow rate on a boundary face
    /*!
     *  @param flag boundary flag
     *  @return derivative of flow rate
     */
    Real linearFlux ( const markerID_Type& flag );
    LIFEV_DEPRECATED ( Real getLinearFlux ( const markerID_Type& flag ) );

    //! Compute the derivative of the pressure on a boundary face
    /*!
     *  @param flag boundary flag
     *  @return derivative of pressure
     */
    Real linearPressure ( const markerID_Type& flag );
    LIFEV_DEPRECATED ( Real getLinearPressure ( const markerID_Type& flag ) );

    //! Compute the derivative of a Lagrange multiplier (which correspond to the value of the derivative of the mean normal stress on a boundary face)
    /*!
     *  @param flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     *  @param BC BChandler containing the boundary conditions of the problem.
     *  @return derivative of the Lagrange multiplier
     */
    Real linearLagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler );
    LIFEV_DEPRECATED ( Real getLinearLagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler ) );

    //! Compute the derivative of the kinetic normal stress (i.e., the normal stress due to the kinetic energy) on a boundary face with given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag boundary flag
     *  @return derivative of the kinetic normal stress
     */
    Real linearKineticNormalStress ( const markerID_Type& flag );

    //! Compute the derivative of the kinetic normal stress (i.e., the normal stress due to the kinetic energy) on a boundary face with a given flag and a given solution
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag boundary flag
     *  @param solution problem solution
     *  @param linearSolution linear problem solution
     *  @return derivative of the kinetic normal stress
     */
    Real linearKineticNormalStress ( const markerID_Type& flag, const vector_Type& solution, const vector_Type& linearSolution );

    //! Compute the derivative of the mean normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     *  @param BC BChandler containing the boundary conditions of the problem.
     *  @return derivative of the mean normal stress
     */
    Real linearMeanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler );

    //! Compute the derivative of the mean normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     *  @param BC BChandler containing the boundary conditions of the problem.
     *  @param linearSolution linear problem solution
     *  @return derivative of the mean normal stress
     */
    Real linearMeanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& linearSolution );

    //! Compute the derivative of the mean total normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     *  @param BC BChandler containing the boundary conditions of the problem.
     *  @return derivative of the mean total normal stress
     */
    Real linearMeanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler );

    //! Compute the derivative of the mean total normal stress on a boundary face with a given flag
    /*!
     *  @see \cite BlancoMalossi2012 \cite Malossi-Thesis
     *
     *  @param flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
     *  @param BC BChandler containing the boundary conditions of the problem.
     *  @param solution problem solution
     *  @param linearSolution linear problem solution
     *  @return derivative of the mean total normal stress
     */
    Real linearMeanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution, const vector_Type& linearSolution  );

    //@}


private:

    //@{

    //! Empty copy constructor
    OseenSolverShapeDerivative ( const OseenSolverShapeDerivative& oseenShapeDerivative );

    //@}

    vector_Type                       M_linearRightHandSideNoBC;
    vector_Type                       M_linearRightHandSideFull;

    vector_Type                       M_linearSolution;

    linearSolver_Type                 M_linearLinSolver;
    preconditionerPtr_type            M_linearPreconditioner;


    VectorElemental                   M_elementVectorVelocity; // Elementary right hand side for the linearized velocity
    VectorElemental                   M_elementVectorPressure; // Elementary right hand side for the linearized pressure
    //    boost::shared_ptr<MatrixElemental>                   M_elementMatrixVelocity;
    //    boost::shared_ptr<MatrixElemental>                   M_elementMatrixConvective;
    //    boost::shared_ptr<MatrixElemental>                   M_elementMatrixPressure;    // Elementary displacement for right hand side
    VectorElemental                   M_elementMeshVelocity;    // Elementary mesh velocity
    VectorElemental                   M_elementVelocity;   // Elementary velocity
    VectorElemental                   M_elementPressure;   // Elementary pressure
    VectorElemental                   M_elementConvectionVelocity;    // Elementary convection velocity
    VectorElemental                   M_elementDisplacement;    // Elementary displacement for right hand side
    VectorElemental                   M_elementVelocityRightHandSide;   // Elementary mesh velocity for right hand side
    VectorElemental                   M_u_loc;
    bool                              M_reuseLinearPreconditioner;
    FESpace<mesh_Type, MapEpetra>*    M_mmFESpace;
};



// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
OseenSolverShapeDerivative<MeshType, SolverType>::
OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                             FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                             FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                             boost::shared_ptr<Epetra_Comm>& communicator,
                             const Int                       lagrangeMultiplier) :
    oseenSolver_Type ( dataType,
                       velocityFESpace,
                       pressureFESpace,
                       communicator,
                       lagrangeMultiplier ),
    M_linearRightHandSideNoBC     ( this->getMap() ),
    M_linearRightHandSideFull     ( this->getMap() ),
    M_linearSolution         ( this->getMap() ),
    M_linearLinSolver ( communicator ),
    M_linearPreconditioner        ( ),
    M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVectorPressure ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    //    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbFEDof(), nDimensions ),
    M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocity         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementPressure         ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementDisplacement          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_u_loc          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_reuseLinearPreconditioner   ( true )
{

}

template<typename MeshType, typename SolverType>
OseenSolverShapeDerivative<MeshType, SolverType>::
OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                             FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                             FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                             boost::shared_ptr<Epetra_Comm>& communicator,
                             const MapEpetra                 monolithicMap,
                             const UInt                      offset ) :
    oseenSolver_Type ( dataType,
                       velocityFESpace,
                       pressureFESpace,
                       communicator,
                       monolithicMap,
                       offset),
    M_linearRightHandSideNoBC     ( this->getMap() ),
    M_linearRightHandSideFull     ( this->getMap() ),
    M_linearSolution         ( this->getMap() ),
    M_linearLinSolver ( communicator ),
    M_linearPreconditioner        ( ),
    M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVectorPressure ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    //    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbFEDof(), nDimensions ),
    M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocity         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementPressure         ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementDisplacement          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_u_loc          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_reuseLinearPreconditioner   ( true )
{

}

template<typename MeshType, typename SolverType>
OseenSolverShapeDerivative<MeshType, SolverType>::
OseenSolverShapeDerivative ( boost::shared_ptr<data_Type>    dataType,
                             FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                             FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                             FESpace<mesh_Type, MapEpetra>&  mmFESpace,
                             boost::shared_ptr<Epetra_Comm>& communicator,
                             const MapEpetra                 monolithicMap,
                             const UInt                      offset) :
    oseenSolver_Type  (dataType,
                       velocityFESpace,
                       pressureFESpace,
                       communicator,
                       monolithicMap,
                       offset),
    M_linearRightHandSideNoBC     ( this->getMap() ),
    M_linearRightHandSideFull     ( this->getMap() ),
    M_linearSolution         ( this->getMap() ),
    M_linearLinSolver( ),
    M_linearPreconditioner        ( ),
    M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVectorPressure ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    //    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbFEDof(), nDimensions ),
    M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocity         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementPressure         ( this->M_pressureFESpace.fe().nbFEDof(), 1 ),
    M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementDisplacement          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_u_loc          ( this->M_velocityFESpace.fe().nbFEDof(), nDimensions ),
    M_reuseLinearPreconditioner   ( true ),
    M_mmFESpace      ( &mmFESpace )
{

}

template<typename MeshType, typename SolverType>
OseenSolverShapeDerivative<MeshType, SolverType>::
~OseenSolverShapeDerivative()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearFlux ( const markerID_Type& flag )
{
    return this->flux ( flag, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::getLinearFlux ( const markerID_Type& flag )
{
    if ( this->M_displayer->isLeader() )
    {
        std::cerr << "Warning: getLinearFlux is deprecated!" << std::endl
                  << "         You should use linearFlux instead!" << std::endl;
    }

    return linearFlux ( flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearPressure ( const markerID_Type& flag )
{
    return this->pressure ( flag, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::getLinearPressure ( const markerID_Type& flag )
{
    if ( this->M_displayer->isLeader() )
    {
        std::cerr << "Warning: getLinearPressure is deprecated!" << std::endl
                  << "         You should use linearPressure instead!" << std::endl;
    }

    return linearPressure ( flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearLagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return lagrangeMultiplier ( flag, bcHandler, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::getLinearLagrangeMultiplier ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    if ( this->M_displayer->isLeader() )
    {
        std::cerr << "Warning: getLinearLagrangeMultiplier is deprecated!" << std::endl
                  << "         You should use linearLagrangeMultiplier instead!" << std::endl;
    }

    return linearLagrangeMultiplier ( flag, bcHandler );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearKineticNormalStress ( const markerID_Type& flag )
{
    return this->linearKineticNormalStress ( flag, *this->M_solution, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearKineticNormalStress ( const markerID_Type& flag, const vector_Type& solution, const vector_Type& linearSolution )
{
    vector_Type velocityAndPressure ( solution, Repeated );
    vector_Type velocity ( this->M_velocityFESpace.map(), Repeated );
    velocity.subset ( velocityAndPressure );

    vector_Type linearVelocityAndPressure ( linearSolution, Repeated );
    vector_Type linearVelocity ( this->M_velocityFESpace.map(), Repeated );
    linearVelocity.subset ( linearVelocityAndPressure );

    return this->M_postProcessing->kineticNormalStressDerivative ( velocity, linearVelocity, this->M_oseenData->density(), flag );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearMeanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return this->linearMeanNormalStress ( flag, bcHandler, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearMeanNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& linearSolution )
{
    return this->meanNormalStress ( flag, bcHandler, linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearMeanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler )
{
    return this->meanNormalStress ( flag, bcHandler, M_linearSolution ) - linearKineticNormalStress ( flag, *this->M_solution, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenSolverShapeDerivative<MeshType, SolverType>::linearMeanTotalNormalStress ( const markerID_Type& flag, bcHandler_Type& bcHandler, const vector_Type& solution, const vector_Type& linearSolution )
{
    return this->meanNormalStress ( flag, bcHandler, linearSolution ) - linearKineticNormalStress ( flag, solution, linearSolution );
}

template<typename MeshType, typename SolverType>
void OseenSolverShapeDerivative<MeshType, SolverType>::setUp ( const GetPot& dataFile )
{
    // M_linearLinSolver.setDataFromGetPot( dataFile, "lin_fluid/solver" );
    // M_linearLinSolver.setAztecooPreconditioner( dataFile, "lin_fluid/solver" );

    oseenSolver_Type::setUp ( dataFile );

    M_reuseLinearPreconditioner = dataFile ( "lin_fluid/prec/reuse", true );

    // std::string preconditionerType = dataFile( "lin_fluid/prec/prectype", "Ifpack" );

    // M_linearPreconditioner = prec_type( PRECFactory::instance().createObject( preconditionerType ) ); //linPrec is not used
    // M_linearPreconditioner->setDataFromGetPot( dataFile, "lin_fluid/prec" );

}


template<typename MeshType, typename SolverType>
void OseenSolverShapeDerivative<MeshType, SolverType>::solveLinearSystem ( bcHandler_Type& bcHandler )
{
    this->M_Displayer.leaderPrint ( " LF-  Finalizing the matrix and vectors ...    " );

    LifeChrono chrono;
    chrono.start();

    // matrix and vector assembling communication
    this->M_matrixNoBC->globalAssemble();

    M_linearRightHandSideNoBC.globalAssemble();

    matrixPtr_Type matrixFull ( new matrix_Type ( this->M_localMap, this->M_matrixNoBC->meanNumEntries() ) );

    this->updateStabilization ( *matrixFull );
    this->getFluidMatrix ( *matrixFull );

    vector_Type    rightHandSideFull ( M_linearRightHandSideNoBC );

    chrono.stop();
    this->M_Displayer.leaderPrintMax ( "done in " , chrono.diff() );

    // boundary conditions update
    this->M_Displayer.leaderPrint ( " LF-  Applying boundary conditions ...         " );
    chrono.start();

    this->applyBoundaryConditions ( *matrixFull, rightHandSideFull, bcHandler );

    chrono.stop();
    this->M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );

    // solving the system

    // using the same preconditioner as for the non linear problem (the matrix changes only in the
    // boundary terms).
    matrixFull->globalAssemble();
    this->M_linearSolver->setMatrix ( *matrixFull );
    this->M_linearSolver->setReusePreconditioner ( M_reuseLinearPreconditioner );
    boost::shared_ptr<MatrixEpetra<Real> > staticCast = boost::static_pointer_cast<MatrixEpetra<Real> > (matrixFull);
    this->M_linearSolver->solveSystem ( rightHandSideFull, M_linearSolution, staticCast );

    *this->M_residual  = M_linearRightHandSideNoBC;
    *this->M_residual -= *this->M_matrixNoBC * this->M_linearSolution;

    //     if(S_verbose)
    //         {
    //             this->M_Displayer.leaderPrintMax( "NormInf Residual Lin = " , this->M_residual.NormInf());
    //             this->M_Displayer.leaderPrintMax( "NormInf Solution Lin = " , this->M_linearSolution.NormInf());
    //         }
} // solveLinearSystem

template<typename MeshType, typename SolverType>
void
OseenSolverShapeDerivative<MeshType, SolverType>::updateLinearSystem ( const matrix_Type& /*matrixNoBC*/,
                                                                       Real&              /*alpha*/,
                                                                       const vector_Type& un,
                                                                       const vector_Type& uk,
                                                                       const vector_Type& disp,
                                                                       const vector_Type& w,
                                                                       const vector_Type& dw,
                                                                       const vector_Type& sourceVector )
{
    this->M_Displayer.leaderPrint ( " LF-  Updating the right hand side ...         " );
    LifeChrono chrono;
    chrono.start();

    Int numVelocityComponent = nDimensions;
    // sourceVector is usually zero
    M_linearRightHandSideNoBC = sourceVector;

    if ( this->M_oseenData->useShapeDerivatives() )
    {
        // right hand side for the linearized ale system
        // Loop on elements
        vector_Type unRepeated   ( un  , Repeated );
        vector_Type ukRepeated   ( uk  , Repeated );
        vector_Type dispRepeated ( disp, Repeated );
        vector_Type wRepeated    ( w   , Repeated );
        vector_Type dwRepeated   ( dw  , Repeated );

        //     std::cout << wRepeated.NormInf() << std::endl;
        //     std::cout << dwRepeated.NormInf() << std::endl;
        //     std::cout << dispRepeated.NormInf() << std::endl;

        vector_Type linearRightHandSideNoBC ( M_linearRightHandSideNoBC.map(), Repeated );

        for ( UInt i = 0; i < this->M_velocityFESpace.mesh()->numVolumes(); i++ )
        {

            this->M_pressureFESpace.fe().update ( this->M_pressureFESpace.mesh()->volumeList ( i ) );
            this->M_velocityFESpace.fe().updateFirstDerivQuadPt ( this->M_velocityFESpace.mesh()->volumeList ( i ) );

            // initialization of elementary vectors
            M_elementVectorVelocity.zero();
            M_elementVectorPressure.zero();

            for ( UInt iNode = 0 ; iNode < this->M_velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = this->M_velocityFESpace.fe().patternFirst ( iNode ); // iLocal = iNode

                for ( Int iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = this->M_velocityFESpace.dof().localToGlobalMap ( i, iLocal ) + iComponent * this->dimVelocity();

                    // u^n - w^iNode local
                    M_elementConvectionVelocity.vec( )  [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = unRepeated (iGlobal)
                            - wRepeated ( iGlobal );
                    // w^iNode local
                    M_elementMeshVelocity.vec( )  [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = wRepeated ( iGlobal );
                    // u^iNode local
                    M_elementVelocity.vec( ) [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = ukRepeated ( iGlobal );
                    // d local
                    M_elementDisplacement.vec( )  [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = dispRepeated ( iGlobal );
                    // dw local
                    M_elementVelocityRightHandSide.vec( ) [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = dwRepeated ( iGlobal );
                    // un local
                    M_u_loc.vec()   [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = unRepeated ( iGlobal );
                }
            }
            /*
              std::cout << M_elementConvectionVelocity.vec() << std::endl;
              std::cout << M_elementMeshVelocity.vec() << std::endl;
              std::cout << M_elementVelocity.vec() << std::endl;
              std::cout << M_elementDisplacement.vec() << std::endl;
              std::cout << M_elementVelocityRightHandSide.vec() << std::endl;
              std::cout << M_u_loc.vec() << std::endl;
            */
            for ( UInt iNode = 0 ; iNode < this->M_pressureFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = this->M_pressureFESpace.fe().patternFirst ( iNode ); // iLocal = iNode
                UInt iGlobal   = this->M_pressureFESpace.dof().localToGlobalMap ( i, iLocal ) + numVelocityComponent * this->dimVelocity();
                M_elementPressure[ iLocal ] = ukRepeated[ iGlobal ];  // p^iNode local

            }

            // Elementary vectors

            //commented the code to print out the elementary data. Useful for debugging.

            //  - \rho ( \grad( u^n-w^iNode ):[I\div d - (\grad d)^T] u^iNode + ( u^n-w^iNode )^T[I\div d - (\grad d)^T] (\grad u^iNode)^T , v  )
            source_mass1 ( - this->M_oseenData->density(),
                           M_elementVelocity,
                           M_elementMeshVelocity,
                           M_elementConvectionVelocity,
                           M_elementDisplacement,
                           M_elementVectorVelocity,
                           this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_mass1 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  + \rho * ( \grad u^iNode dw, v  )
            source_mass2 ( this->M_oseenData->density(),
                           M_elementVelocity,
                           M_elementVelocityRightHandSide,
                           M_elementVectorVelocity,
                           this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_mass2 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^iNode , v  )
            source_mass3 ( - 0.5 * this->M_oseenData->density(),
                           M_u_loc,
                           M_elementVelocity,
                           M_elementDisplacement,
                           M_elementVectorVelocity,
                           this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_mass3 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  - ( [-p^iNode I + 2*mu e(u^iNode)] [I\div d - (\grad d)^T] , \grad v  )
            source_stress ( - 1.0,
                            this->M_oseenData->viscosity(),
                            M_elementVelocity,
                            M_elementPressure,
                            M_elementDisplacement,
                            M_elementVectorVelocity,
                            this->M_velocityFESpace.fe(),
                            this->M_pressureFESpace.fe() );
            /*
            std::cout << "source_stress -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            // + \mu ( \grad u^iNode \grad d + [\grad d]^T[\grad u^iNode]^T : \grad v )
            source_stress2 ( this->M_oseenData->viscosity(),
                             M_elementVelocity,
                             M_elementDisplacement,
                             M_elementVectorVelocity,
                             this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_stress2 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  + ( (\grad u^iNode):[I\div d - (\grad d)^T] , q  )
            source_press ( -1.0,
                           M_elementVelocity,
                           M_elementDisplacement,
                           M_elementVectorPressure,
                           this->M_velocityFESpace.fe(),
                           this->M_pressureFESpace.fe() );
            /*
            std::cout << "source_press -> norm_inf(M_elementVectorVelocity)"  << std::endl;
            M_elementVectorPressure.showMe(std::cout);
            */
            //
            // Assembling
            //
            /*
              std::cout << "debut ====================" << std::endl;
              M_elementVectorPressure.showMe(std::cout);
              M_elementVectorVelocity.showMe(std::cout);
              std::cout << "fin   ====================" << std::endl;
            */
            // assembling pressure right hand side
            assembleVector ( linearRightHandSideNoBC,
                             M_elementVectorPressure,
                             this->M_pressureFESpace.fe(),
                             this->M_pressureFESpace.dof(),
                             0,
                             numVelocityComponent * this->dimVelocity() );

            // loop on velocity components
            for ( Int iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
            {
                // assembling velocity right hand side
                assembleVector ( linearRightHandSideNoBC,
                                 M_elementVectorVelocity,
                                 this->M_velocityFESpace.fe(),
                                 this->M_velocityFESpace.dof(),
                                 iComponent,
                                 iComponent * this->dimVelocity() );
            }
        }

        linearRightHandSideNoBC.globalAssemble();
        M_linearRightHandSideNoBC += linearRightHandSideNoBC;
        //       M_linearRightHandSideNoBC *= -1.;
        //       if( S_verbose )
        //           this->M_Displayer.leaderPrint( "norm( M_linearRightHandSideNoBC)  = " , M_linearRightHandSideNoBC.NormInf() );

    }

    chrono.stop();
    this->M_Displayer.leaderPrintMax ( "done in ", chrono.diff() );
}


//#if UNDEF
template<typename MeshType, typename SolverType>
void
OseenSolverShapeDerivative<MeshType, SolverType>::
updateShapeDerivatives ( matrix_Type&                   matrix,
                         Real&                          alpha,
                         const vector_Type&             un,
                         const vector_Type&             uk,
                         //const vector_Type&           disp,
                         const vector_Type&             w,
                         UInt                           offset,
                         FESpace<mesh_Type, MapEpetra>& mmFESpace,
                         bool                           wImplicit,
                         bool                           convectiveTermDerivative )
{
    LifeChrono chrono;

    UInt numVelocityComponent = nDimensions;

    //    M_linearRightHandSideNoBC = sourceVector;//which is usually zero

    if ( this->M_oseenData->useShapeDerivatives() )
    {
        this->M_Displayer.leaderPrint ( " LF-  Updating shape derivative blocks ...     " );

        //
        // RIGHT HAND SIDE FOR THE LINEARIZED ALE SYSTEM
        //
        chrono.start();

        // Loop on elements

        vector_Type unRepeated  ( un  , Repeated );
        vector_Type ukRepeated  ( uk  , Repeated );
        // vector_Type dispRepeated( disp, Repeated );
        vector_Type wRepeated   ( w   , Repeated );
        // vector_Type dwRepeated  ( dw  , Repeated );

        //     std::cout << wRepeated.NormInf() << std::endl;
        //     std::cout << dwRepeated.NormInf() << std::endl;
        //     std::cout << dispRepeated.NormInf() << std::endl;

        //            vector_Type rhsLinNoBC( M_linearRightHandSideNoBC.map(), Repeated);

        for ( UInt i = 0; i < this->M_velocityFESpace.mesh()->numVolumes(); i++ )
        {

            this->M_pressureFESpace.fe().update ( this->M_pressureFESpace.mesh()->volumeList ( i ) );
            this->M_velocityFESpace.fe().updateFirstDerivQuadPt ( this->M_velocityFESpace.mesh()->volumeList ( i ) );
            this->M_pressureFESpace.fe().updateFirstDerivQuadPt ( this->M_velocityFESpace.mesh()->volumeList ( i ) );

            // just to provide the id number in the assem_mat_mixed
            //this->M_pressureFESpace.fe().updateFirstDeriv( this->M_velocityFESpace.mesh()->volumeList( i ) );
            //as updateFirstDer
            //this->M_velocityFESpace.fe().updateFirstDeriv( this->M_velocityFESpace.mesh()->volumeList( i ) );
            mmFESpace.fe().updateFirstDerivQuadPt ( mmFESpace.mesh()->volumeList ( i ) );

            // initialization of elementary vectors
            boost::shared_ptr<MatrixElemental> elementMatrixPressure ( new MatrixElemental ( this->M_pressureFESpace.fe().nbFEDof(),
                                                                       1,
                                                                       0,
                                                                       mmFESpace.fe().nbFEDof(),
                                                                       0,
                                                                       nDimensions ) );
            boost::shared_ptr<MatrixElemental> elementMatrixVelocity ( new MatrixElemental ( this->M_velocityFESpace.fe().nbFEDof(),
                                                                       nDimensions,
                                                                       0,
                                                                       this->M_velocityFESpace.fe().nbFEDof(),
                                                                       0,
                                                                       nDimensions ) );
            boost::shared_ptr<MatrixElemental> elementMatrixConvective;

            if ( convectiveTermDerivative )
            {
                elementMatrixConvective.reset( new MatrixElemental( this->M_velocityFESpace.fe().nbFEDof(),
                                                            nDimensions,
                                                            0,
                                                            this->M_velocityFESpace.fe().nbFEDof(),
                                                            0,
                                                            nDimensions ) );
                elementMatrixConvective->zero();
            }

            elementMatrixPressure->zero();
            elementMatrixVelocity->zero();

            for ( UInt iNode = 0 ; iNode < this->M_velocityFESpace.fe().nbFEDof() ; iNode++ )
            {
                UInt iLocal = this->M_velocityFESpace.fe().patternFirst ( iNode ); // iLocal = iNode

                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = this->M_velocityFESpace.dof().localToGlobalMap ( i, iLocal ) + iComponent * this->dimVelocity();

                    // if(!wImplicit)
                    // u^n - w^iNode local
                    M_elementConvectionVelocity.vec() [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = unRepeated (iGlobal)
                            - wRepeated ( iGlobal );
                    // else
                    // u^n - w^iNode local
                    // M_elementConvectionVelocity.vec() [ iLocal + iComponent*this->M_velocityFESpace.fe().nbFEDof() ] = ukRepeated(iGlobal)
                    // - wRepeated(iGlobal);
                    // w^iNode local
                    M_elementMeshVelocity.vec( )  [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = wRepeated ( iGlobal );
                    // u^iNode local
                    M_elementVelocity.vec( ) [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = ukRepeated ( iGlobal );
                    // dw local
                    //M_elementDisplacement.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbFEDof() ] = dispRepeated( iGlobal );
                    // dw local
                    //M_elementVelocityRightHandSide.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbFEDof() ] = dwRepeated( iGlobal );
                    // un local
                    M_u_loc.vec()   [ iLocal + iComponent * this->M_velocityFESpace.fe().nbFEDof() ] = unRepeated ( iGlobal );
                }
            }
            /*
            std::cout << M_elementConvectionVelocity.vec() << std::endl;
            std::cout << M_elementMeshVelocity.vec() << std::endl;
            std::cout << M_elementVelocity.vec() << std::endl;
            std::cout << M_elementDisplacement.vec() << std::endl;
            std::cout << M_elementVelocityRightHandSide.vec() << std::endl;
            std::cout << M_u_loc.vec() << std::endl;
            */
            for ( UInt iNode = 0 ; iNode < this->M_pressureFESpace.fe().nbFEDof() ; iNode++ )
            {
                // iLocal = iNode
                UInt iLocal = this->M_pressureFESpace.fe().patternFirst ( iNode );
                UInt iGlobal = this->M_pressureFESpace.dof().localToGlobalMap ( i, iLocal ) + numVelocityComponent * this->dimVelocity();
                // p^iNode local
                M_elementPressure[ iLocal ] = ukRepeated[ iGlobal ];
            }
	    
	    shape_terms( //M_elementDisplacement,
			this->M_oseenData->density(),
			this->M_oseenData->viscosity(),
			M_u_loc,
			M_elementVelocity,
			M_elementMeshVelocity,
			M_elementConvectionVelocity,
			M_elementPressure,
			*elementMatrixVelocity,
			mmFESpace.fe(),
			this->M_pressureFESpace.fe(),
			(ID) mmFESpace.fe().nbFEDof(),
			*elementMatrixPressure,
			0,
			wImplicit,
			alpha//,
			//elementMatrixConvective
			 );

            //elementMatrixVelocity->showMe(std::cout);

            /*
            source_mass2( this->M_oseenData->density(),
                          M_elementVelocity,
                          *M_elementMatrixConvective,
                          this->M_velocityFESpace.fe(),
                          alpha );
            */

            source_press( 1.0,
                          M_elementVelocity,
                          *elementMatrixPressure,
                          mmFESpace.fe(),
                          this->M_pressureFESpace.fe(),
                          (ID) mmFESpace.fe().nbFEDof() );

            //derivative of the convective term
            if ( convectiveTermDerivative )
                mass_gradu ( this->M_oseenData->density(),
                             M_elementVelocity,
                             *elementMatrixConvective,
                             this->M_velocityFESpace.fe() );
            /*
              std::cout << "source_press -> norm_inf( M_elementVectorVelocity )"  << std::endl;
            M_elementVectorPressure.showMe( std::cout );
            */
            //
            // Assembling
            //
            /*
            std::cout << "debut ====================" << std::endl;
            M_elementVectorPressure.showMe( std::cout );
            M_elementVectorVelocity.showMe( std::cout );
            std::cout << "fin   ====================" << std::endl;
            */
            UInt const velocityTotalDof ( this->M_velocityFESpace.dof().numTotalDof() );
            UInt const meshTotalDof ( mmFESpace.dof().numTotalDof() );
            for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
            {
                for ( UInt jComponent = 0; jComponent < numVelocityComponent; ++jComponent )
                {
                    assembleMatrix( matrix,
                                    *elementMatrixVelocity,
                                    this->M_velocityFESpace.fe(),
                                    mmFESpace.fe(),
                                    this->M_velocityFESpace.dof(),
                                    mmFESpace.dof(),
                                    iComponent,
                                    jComponent,
                                    iComponent * velocityTotalDof,
                                    offset + jComponent * meshTotalDof );

                    //assembling the derivative of the convective term
                    if ( convectiveTermDerivative )
                        assembleMatrix ( matrix,
                                         *elementMatrixConvective,
                                         this->M_velocityFESpace.fe(),
                                         this->M_velocityFESpace.fe(),
                                         this->M_velocityFESpace.dof(),
                                         this->M_velocityFESpace.dof(),
                                         iComponent,
                                         jComponent,
                                         iComponent * velocityTotalDof,
                                         jComponent * velocityTotalDof );
                }

	    assembleMatrix( matrix,
			    *elementMatrixPressure,
			    this->M_pressureFESpace.fe(),
			    mmFESpace.fe(),
			    this->M_pressureFESpace.dof(),
			    mmFESpace.dof(),
			    (UInt) 0,
			    iComponent,
			    (UInt) numVelocityComponent * velocityTotalDof,
			    offset + iComponent * meshTotalDof );
            }
        }
    }

    chrono.stop();
    this->M_Displayer.leaderPrintMax ("done in ", chrono.diff() );
}

} // namespace LifeV

#endif // OSEENSOLVERSHAPEDERIVATIVE_H
