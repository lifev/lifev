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


#ifndef OSEENSHAPEDERIVATIVE_H
#define OSEENSHAPEDERIVATIVE_H 1

#include <life/lifesolver/Oseen.hpp>

namespace LifeV
{

//! @class OseenShapeDerivative
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

template< typename MeshType, typename SolverType = LifeV::SolverTrilinos >
class OseenShapeDerivative:
        public Oseen< MeshType, SolverType >
{

public:

    //! @name Public Types
    //@{

    ///////// old: to be removed
    typedef Oseen< MeshType, SolverType >                     oseenSolver_Type;
    typedef typename oseenSolver_Type::vector_type            vector_type;
    typedef typename oseenSolver_Type::matrix_type            matrix_type;
    typedef typename oseenSolver_Type::matrix_ptrtype         matrix_ptrtype;
    typedef typename oseenSolver_Type::data_type              data_type;
    typedef typename oseenSolver_Type::prec_type              prec_type;
    typedef typename oseenSolver_Type::prec_raw_type          prec_raw_type;
    typedef typename oseenSolver_Type::bchandler_raw_type     bchandler_raw_type;
    /////////

    typedef MeshType                                          mesh_Type;
    typedef SolverType                                        linearSolver_Type;
    // typedef Oseen< mesh_Type, linearSolver_Type >             oseenSolver_Type;
    typedef typename oseenSolver_Type::vector_Type            vector_Type;
    typedef typename oseenSolver_Type::matrix_Type            matrix_Type;
    typedef typename oseenSolver_Type::matrixPtr_Type         matrixPtr_Type;
    typedef typename oseenSolver_Type::data_Type              data_Type;
    typedef typename oseenSolver_Type::preconditionerPtr_Type preconditionerPtr_type;
    typedef typename oseenSolver_Type::preconditioner_Type    preconditioner_Type;
    typedef typename oseenSolver_Type::bcHandler_Type         bcHandler_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    OseenShapeDerivative();

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */
    OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                          FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                          FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                          boost::shared_ptr<Epetra_Comm>& communicator,
                          const Int                       lagrangeMultiplier = 0);

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param monolithicMap EpetraMap class
        @param offset
     */
    OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                          FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                          FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                          boost::shared_ptr<Epetra_Comm>& communicator,
                          const EpetraMap                 monolithicMap,
                          const UInt                      offset = 0 );

    //! Constructor
    /*!
        @param dataType DataNavierStokes class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param mmFESpace FE space
        @param communicator MPI communicator
        @param monolithicMap EpetraMap class
        @param offset
     */
    OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                          FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                          FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                          FESpace<mesh_Type, EpetraMap>&  mmFESpace,
                          boost::shared_ptr<Epetra_Comm>& communicator,
                          const EpetraMap                 monolithicMap,
                          const UInt                      offset = 0 );

    //! Virtual destructor
    virtual ~OseenShapeDerivative();

    //@}


    //! @name Methods
    //@{

    //! Set up data from GetPot
    /*!
        @param dataFile GetPot object
     */
    void setUp( const GetPot& dataFile );

    //!
    /*!
        @param bcHandler BC handler
     */
    void iterateLin( bcHandler_Type& bcHandler );

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
    void updateLinearSystem( const matrix_Type& matrixNoBC,
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
    void updateShapeDerivatives( matrix_Type&                   matrixNoBC,
                                 Real&                          alpha,
                                 const vector_Type&             un,
                                 const vector_Type&             uk,
                                 //const vector_Type&           disp,
                                 const vector_Type&             w,
                                 UInt                           offset,
                                 FESpace<mesh_Type, EpetraMap>& dFESpace,
                                 bool                           wImplicit = true,
                                 bool                           convectiveTermDerivative = false);

    //@}


    //! @name Set Methods
    //@{

    //! Set
    /*!
        @param rightHandSide
     */
    void updateRhsLinNoBC( const vector_Type& rightHandSide)
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
    vector_Type& rhsLinNoBC()
    {
        return M_linearRightHandSideNoBC;
    }

    const vector_Type& rhsLinNoBC() const
    {
        return M_linearRightHandSideNoBC;
    }

    //! Get the solution of the Shape Derivative problem.
    /*!
        @return vector containing the solution of the Shape Derivative problem.
     */
    const vector_Type& LinearSolution() const
    {
        return M_linearSolution;
    }

    //! Return
    /*!
        @return stab
     */
    bool stab()
    {
        return this->M_stabilization;
    }

    const bool& stab() const
    {
        return this->M_stabilization;
    }

    //! Return
    /*!
        @return linear flux
     */
    Real GetLinearFlux ( const EntityFlag& flag );

    //! Return
    /*!
        @return linear pressure
     */
    Real GetLinearPressure( const EntityFlag& flag );

    //! Get the Lagrange multiplier related to a flux imposed on a given part of the boundary.
    /*!
        @param Flag flag of the boundary face associated with the flux and the Lagrange multiplier we want.
        @param BC BChandler containing the boundary conditions of the problem.
        @return Lagrange multiplier
     */
    Real LinearLagrangeMultiplier( const EntityFlag& flag, bcHandler_Type& bcHandler );

    //@}


private:

    vector_Type               M_linearRightHandSideNoBC;
    vector_Type               M_linearRightHandSideFull;

    vector_Type               M_linearSolution;

    linearSolver_Type         M_linearLinSolver;
    prec_type                 M_linearPreconditioner;


    ElemVec                   M_elementVectorVelocity; // Elementary right hand side for the linearized pressure
    ElemVec                   M_elementVectorPressure; // Elementary right hand side for the linearized pressure
    //    boost::shared_ptr<ElemMat>                   M_elementMatrixVelocity;
    //    boost::shared_ptr<ElemMat>                   M_elementMatrixConvective;
    //    boost::shared_ptr<ElemMat>                   M_elementMatrixPressure;    // Elementary displacement for right hand side
    ElemVec                   M_elementMeshVelocity;    // Elementary mesh velocity
    ElemVec                   M_elementVelocity;   // Elementary velocity
    ElemVec                   M_elementPressure;   // Elementary pressure
    ElemVec                   M_elementConvectionVelocity;    // Elementary convection velocity
    ElemVec                   M_elementDisplacement;    // Elementary displacement for right hand side
    ElemVec                   M_elementVelocityRightHandSide;   // Elementary mesh velocity for right hand side
    ElemVec                   M_u_loc;
    bool                      M_reuseLinearPreconditioner;
    FESpace<mesh_Type, EpetraMap>* M_mmFESpace;
};



// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
OseenShapeDerivative<MeshType, SolverType>::
OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                      FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                      FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                      boost::shared_ptr<Epetra_Comm>& communicator,
                      const Int                       lagrangeMultiplier):
        oseenSolver_Type ( dataType,
                           velocityFESpace,
                           pressureFESpace,
                           communicator,
                           lagrangeMultiplier ),
        M_linearRightHandSideNoBC     ( this->getMap() ),
        M_linearRightHandSideFull     ( this->getMap() ),
        M_linearSolution         ( this->getMap() ),
        M_linearLinSolver( communicator ),
        M_linearPreconditioner        ( ),
        M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVectorPressure ( this->M_pressureFESpace.fe().nbNode, 1 ),
//    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbNode, nDimensions ),
        M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocity         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementPressure         ( this->M_pressureFESpace.fe().nbNode, 1 ),
        M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementDisplacement          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_u_loc          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_reuseLinearPreconditioner   ( true )
{

}

template<typename MeshType, typename SolverType>
OseenShapeDerivative<MeshType, SolverType>::
OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                      FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                      FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                      boost::shared_ptr<Epetra_Comm>& communicator,
                      const EpetraMap                 monolithicMap,
                      const UInt                      offset ):
        oseenSolver_Type ( dataType,
                           velocityFESpace,
                           pressureFESpace,
                           communicator,
                           monolithicMap,
                           offset),
        M_linearRightHandSideNoBC     ( this->getMap()),
        M_linearRightHandSideFull     ( this->getMap()),
        M_linearSolution         ( this->getMap()),
        M_linearLinSolver( communicator ),
        M_linearPreconditioner        ( ),
        M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVectorPressure ( this->M_pressureFESpace.fe().nbNode, 1 ),
//    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbNode, nDimensions ),
        M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocity         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementPressure         ( this->M_pressureFESpace.fe().nbNode, 1 ),
        M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementDisplacement          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_u_loc          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_reuseLinearPreconditioner   ( true )
{

}

template<typename MeshType, typename SolverType>
OseenShapeDerivative<MeshType, SolverType>::
OseenShapeDerivative( boost::shared_ptr<data_Type>    dataType,
                      FESpace<mesh_Type, EpetraMap>&  velocityFESpace,
                      FESpace<mesh_Type, EpetraMap>&  pressureFESpace,
                      FESpace<mesh_Type, EpetraMap>&  mmFESpace,
                      boost::shared_ptr<Epetra_Comm>& communicator,
                      const EpetraMap                 monolithicMap,
                      const UInt                      offset):
        oseenSolver_Type  (dataType,
                           velocityFESpace,
                           pressureFESpace,
                           communicator,
                           monolithicMap,
                           offset),
        M_linearRightHandSideNoBC     ( this->getMap()),
        M_linearRightHandSideFull     ( this->getMap()),
        M_linearSolution         ( this->getMap()),
        M_linearLinSolver( ),
        M_linearPreconditioner        ( ),
        M_elementVectorVelocity ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVectorPressure ( this->M_pressureFESpace.fe().nbNode, 1 ),
//    M_elementVectorPressure   ( this->M_pressureFESpace.fe().nbNode, nDimensions ),
        M_elementMeshVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocity         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementPressure         ( this->M_pressureFESpace.fe().nbNode, 1 ),
        M_elementConvectionVelocity          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementDisplacement          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_elementVelocityRightHandSide         ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_u_loc          ( this->M_velocityFESpace.fe().nbNode, nDimensions ),
        M_reuseLinearPreconditioner   ( true ),
        M_mmFESpace      ( &mmFESpace )
{

}

template<typename MeshType, typename SolverType>
OseenShapeDerivative<MeshType, SolverType>::
~OseenShapeDerivative()
{

}


// ===================================================
// Methods
// ===================================================

template<typename MeshType, typename SolverType>
Real
OseenShapeDerivative<MeshType, SolverType>::GetLinearFlux( const EntityFlag& flag )
{
    return flux( flag, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenShapeDerivative<MeshType, SolverType>::GetLinearPressure( const EntityFlag& flag )
{
    return pressure( flag, M_linearSolution );
}

template<typename MeshType, typename SolverType>
Real
OseenShapeDerivative<MeshType, SolverType>::LinearLagrangeMultiplier( const EntityFlag& flag, bcHandler_Type& bcHandler )
{
    return LagrangeMultiplier( flag, bcHandler, M_linearSolution );
}


template<typename MeshType, typename SolverType>
void OseenShapeDerivative<MeshType, SolverType>::setUp( const GetPot& dataFile )
{
    // M_linearLinSolver.setDataFromGetPot( dataFile, "lin_fluid/solver" );
    // M_linearLinSolver.setAztecooPreconditioner( dataFile, "lin_fluid/solver" );

    oseenSolver_Type::setUp( dataFile );

    M_reuseLinearPreconditioner = dataFile( "lin_fluid/prec/reuse", true );

    // std::string precType = dataFile( "lin_fluid/prec/prectype", "Ifpack" );

    // M_linearPreconditioner = prec_type( PRECFactory::instance().createObject( precType ) ); //linPrec is not used
    // M_linearPreconditioner->setDataFromGetPot( dataFile, "lin_fluid/prec" );

}


template<typename MeshType, typename SolverType>
void OseenShapeDerivative<MeshType, SolverType>::iterateLin( bcHandler_Type& bcHandler )
{
    this->M_Displayer.leaderPrint( " LF-  Finalizing the matrix and vectors ...    " );

    Chrono chrono;
    chrono.start();

    // matrix and vector assembling communication
    this->M_matrixNoBC->GlobalAssemble();

    M_linearRightHandSideNoBC.GlobalAssemble();

    matrix_ptrtype matrixFull( new matrix_Type( this->M_localMap, this->M_matrixNoBC->getMeanNumEntries() ) );

    updateStab( *matrixFull );
    getFluidMatrix( *matrixFull );

    vector_Type    rightHandSideFull( M_linearRightHandSideNoBC );

    chrono.stop();
    this->M_Displayer.leaderPrintMax( "done in " , chrono.diff() );

    // boundary conditions update
    this->M_Displayer.leaderPrint( " LF-  Applying boundary conditions ...         " );
    chrono.start();

    applyBoundaryConditions( *matrixFull, rightHandSideFull, bcHandler );

    chrono.stop();
    this->M_Displayer.leaderPrintMax( "done in ", chrono.diff() );

    // solving the system

    // using the same preconditioner as for the non linear problem (the matrix changes only in the
    // boundary terms).
    matrixFull->GlobalAssemble();
    this->M_linearSolver.setMatrix( *matrixFull );
    this->M_linearSolver.setReusePreconditioner( M_reuseLinearPreconditioner );
    this->M_linearSolver.solveSystem( rightHandSideFull, M_linearSolution, matrixFull );

    this->M_residual  = M_linearRightHandSideNoBC;
    this->M_residual -= *this->M_matrixNoBC*this->M_linearSolution;

//     if(S_verbose)
//         {
//             this->M_Displayer.leaderPrintMax( "NormInf Residual Lin = " , this->M_residual.NormInf());
//             this->M_Displayer.leaderPrintMax( "NormInf Solution Lin = " , this->M_linearSolution.NormInf());
//         }
} // iterateLin

template<typename MeshType, typename SolverType>
void
OseenShapeDerivative<MeshType, SolverType>::updateLinearSystem( const matrix_Type& matrixNoBC,
                                                                  Real&              alpha,
                                                                  const vector_Type& un,
                                                                  const vector_Type& uk,
                                                                  const vector_Type& disp,
                                                                  const vector_Type& w,
                                                                  const vector_Type& dw,
                                                                  const vector_Type& sourceVector )
{
    this->M_Displayer.leaderPrint( " LF-  Updating the right hand side ...         " );
    Chrono chrono;
    chrono.start();

    Int numVelocityComponent = nDimensions;
    // sourceVector is usually zero
    M_linearRightHandSideNoBC = sourceVector;

    if ( this->M_dataNavierStokes->useShapeDerivatives() )
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

        vector_Type rhsLinNoBC( M_linearRightHandSideNoBC.getMap(), Repeated );

        for ( UInt i = 1; i <= this->M_velocityFESpace.mesh()->numVolumes(); i++ )
        {

            this->M_pressureFESpace.fe().update( this->M_pressureFESpace.mesh()->volumeList( i ) );
            this->M_velocityFESpace.fe().updateFirstDerivQuadPt( this->M_velocityFESpace.mesh()->volumeList( i ) );

            // initialization of elementary vectors
            M_elementVectorVelocity.zero();
            M_elementVectorPressure.zero();

            for ( UInt iNode = 0 ; iNode < this->M_velocityFESpace.fe().nbNode ; iNode++ )
            {
                UInt iLocal = this->M_velocityFESpace.fe().patternFirst( iNode ); // iLocal = iNode

                for ( Int iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = this->M_velocityFESpace.dof().localToGlobal( i, iLocal + 1 ) + iComponent * this->dim_u();

                    // u^n - w^iNode local
                    M_elementConvectionVelocity.vec( )  [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = unRepeated(iGlobal)
                            - wRepeated( iGlobal );
                    // w^iNode local
                    M_elementMeshVelocity.vec( )  [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = wRepeated( iGlobal );
                    // u^iNode local
                    M_elementVelocity.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = ukRepeated( iGlobal );
                    // d local
                    M_elementDisplacement.vec( )  [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = dispRepeated( iGlobal );
                    // dw local
                    M_elementVelocityRightHandSide.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = dwRepeated( iGlobal );
                    // un local
                    M_u_loc.vec()   [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = unRepeated( iGlobal );
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
            for ( UInt iNode = 0 ; iNode < this->M_pressureFESpace.fe().nbNode ; iNode++ )
            {
                UInt iLocal = this->M_pressureFESpace.fe().patternFirst( iNode ); // iLocal = iNode
                UInt iGlobal   = this->M_pressureFESpace.dof().localToGlobal( i, iLocal + 1 ) + numVelocityComponent*this->dim_u();
                M_elementPressure[ iLocal ] = ukRepeated[ iGlobal ];  // p^iNode local

            }

            // Elementary vectors

            //commented the code to print out the elementary data. Useful for debugging.

            //  - \rho ( \grad( u^n-w^iNode ):[I\div d - (\grad d)^T] u^iNode + ( u^n-w^iNode )^T[I\div d - (\grad d)^T] (\grad u^iNode)^T , v  )
            source_mass1( - this->M_dataNavierStokes->density(),
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
            source_mass2( this->M_dataNavierStokes->density(),
                          M_elementVelocity,
                          M_elementVelocityRightHandSide,
                          M_elementVectorVelocity,
                          this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_mass2 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  - \rho/2 ( \nabla u^n:[2 * I\div d - (\grad d)^T]  u^iNode , v  )
            source_mass3( - 0.5*this->M_dataNavierStokes->density(),
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
            source_stress( - 1.0,
                           this->M_dataNavierStokes->viscosity(),
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
            source_stress2( this->M_dataNavierStokes->viscosity(),
                            M_elementVelocity,
                            M_elementDisplacement,
                            M_elementVectorVelocity,
                            this->M_velocityFESpace.fe() );
            /*
            std::cout << "source_stress2 -> norm_inf(M_elementVectorVelocity)" << std::endl;
            M_elementVectorVelocity.showMe(std::cout);
            */
            //  + ( (\grad u^iNode):[I\div d - (\grad d)^T] , q  )
            source_press( -1.0,
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
            assembleVector( rhsLinNoBC,
                            M_elementVectorPressure,
                            this->M_pressureFESpace.fe(),
                            this->M_pressureFESpace.dof(),
                            0,
                            numVelocityComponent*this->dim_u() );

            // loop on velocity components
            for ( Int iComponent = 0; iComponent < numVelocityComponent; iComponent++ )
            {
                // assembling velocity right hand side
                assembleVector( rhsLinNoBC,
                                M_elementVectorVelocity,
                                this->M_velocityFESpace.fe(),
                                this->M_velocityFESpace.dof(),
                                iComponent,
                                iComponent*this->dim_u() );
            }
        }

        rhsLinNoBC.GlobalAssemble();
        M_linearRightHandSideNoBC += rhsLinNoBC;
//       M_linearRightHandSideNoBC *= -1.;
//       if( S_verbose )
//           this->M_Displayer.leaderPrint( "norm( M_linearRightHandSideNoBC)  = " , M_linearRightHandSideNoBC.NormInf() );

    }

    chrono.stop();
    this->M_Displayer.leaderPrintMax( "done in ", chrono.diff() );
}


//#if UNDEF
template<typename MeshType, typename SolverType>
void
OseenShapeDerivative<MeshType, SolverType>::
updateShapeDerivatives( matrix_Type&                   matrix,
                        Real&                          alpha,
                        const vector_Type&             un,
                        const vector_Type&             uk,
                        //const vector_Type&           disp,
                        const vector_Type&             w,
                        UInt                           offset,
                        FESpace<mesh_Type, EpetraMap>& mmFESpace,
                        bool                           wImplicit,
                        bool                           convectiveTermDerivative )
{
    Chrono chrono;

    UInt numVelocityComponent = nDimensions;

    //    M_linearRightHandSideNoBC = sourceVector;//which is usually zero

    if ( this->M_dataNavierStokes->useShapeDerivatives() )
    {
        this->M_Displayer.leaderPrint( " LF-  Updating shape derivative blocks ...     " );

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

//            vector_Type rhsLinNoBC( M_linearRightHandSideNoBC.getMap(), Repeated);

        for ( UInt i = 1; i <= this->M_velocityFESpace.mesh()->numVolumes(); i++ )
        {

            this->M_pressureFESpace.fe().update( this->M_pressureFESpace.mesh()->volumeList( i ) );
            this->M_velocityFESpace.fe().updateFirstDerivQuadPt( this->M_velocityFESpace.mesh()->volumeList( i ) );
            this->M_pressureFESpace.fe().updateFirstDerivQuadPt( this->M_velocityFESpace.mesh()->volumeList( i ) );

            // just to provide the id number in the assem_mat_mixed
            //this->M_pressureFESpace.fe().updateFirstDeriv( this->M_velocityFESpace.mesh()->volumeList( i ) );
            //as updateFirstDer
            //this->M_velocityFESpace.fe().updateFirstDeriv( this->M_velocityFESpace.mesh()->volumeList( i ) );
            mmFESpace.fe().updateFirstDerivQuadPt( mmFESpace.mesh()->volumeList( i ) );

            // initialization of elementary vectors
            boost::shared_ptr<ElemMat> elementMatrixPressure ( new ElemMat( this->M_pressureFESpace.fe().nbNode,
                                                                            1,
                                                                            0,
                                                                            mmFESpace.fe().nbNode,
                                                                            0,
                                                                            nDimensions ) );
            boost::shared_ptr<ElemMat> elementMatrixVelocity ( new ElemMat( this->M_velocityFESpace.fe().nbNode,
                                                                            nDimensions,
                                                                            0,
                                                                            this->M_velocityFESpace.fe().nbNode,
                                                                            0,
                                                                            nDimensions ) );
            boost::shared_ptr<ElemMat> elementMatrixConvective;

            if ( convectiveTermDerivative )
            {
                elementMatrixConvective.reset( new ElemMat( this->M_velocityFESpace.fe().nbNode,
                                                            nDimensions,
                                                            0,
                                                            mmFESpace.fe().nbNode,
                                                            0,
                                                            nDimensions ) );
                elementMatrixConvective->zero();
            }

            elementMatrixPressure->zero();
            elementMatrixVelocity->zero();

            for ( UInt iNode = 0 ; iNode < this->M_velocityFESpace.fe().nbNode ; iNode++ )
            {
                UInt iLocal = this->M_velocityFESpace.fe().patternFirst( iNode ); // iLocal = iNode

                for ( UInt iComponent = 0; iComponent < numVelocityComponent; ++iComponent )
                {
                    UInt iGlobal = this->M_velocityFESpace.dof().localToGlobal( i, iLocal + 1 ) + iComponent * this->dim_u();

                    // if(!wImplicit)
                    // u^n - w^iNode local
                    M_elementConvectionVelocity.vec() [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = unRepeated(iGlobal)
                            - wRepeated( iGlobal );
                    // else
                    // u^n - w^iNode local
                    // M_elementConvectionVelocity.vec() [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = ukRepeated(iGlobal)
                    - wRepeated(iGlobal);
                    // w^iNode local
                    M_elementMeshVelocity.vec( )  [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = wRepeated( iGlobal );
                    // u^iNode local
                    M_elementVelocity.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = ukRepeated( iGlobal );
                    // dw local
                    //M_elementDisplacement.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = dispRepeated( iGlobal );
                    // dw local
                    //M_elementVelocityRightHandSide.vec( ) [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = dwRepeated( iGlobal );
                    // un local
                    M_u_loc.vec()   [ iLocal + iComponent*this->M_velocityFESpace.fe().nbNode ] = unRepeated( iGlobal );
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
            for ( UInt iNode = 0 ; iNode < this->M_pressureFESpace.fe().nbNode ; iNode++ )
            {
                // iLocal = iNode
                UInt iLocal = this->M_pressureFESpace.fe().patternFirst( iNode );
                UInt iGlobal = this->M_pressureFESpace.dof().localToGlobal( i, iLocal + 1 ) + numVelocityComponent*this->dim_u();
                // p^iNode local
                M_elementPressure[ iLocal ] = ukRepeated[ iGlobal ];
            }


            shape_terms( //M_elementDisplacement,
                this->M_dataNavierStokes->density(),
                this->M_dataNavierStokes->viscosity(),
                M_u_loc,
                M_elementVelocity,
                M_elementMeshVelocity,
                M_elementConvectionVelocity,
                M_elementPressure,
                *elementMatrixVelocity,
                this->M_velocityFESpace.fe(),
                this->M_pressureFESpace.fe(),
                (ID) mmFESpace.fe().nbNode,
                *elementMatrixPressure,
                0,
                wImplicit,
                alpha//,
                //elementMatrixConvective
            );

            //elementMatrixVelocity->showMe(std::cout);

            /*
            source_mass2( this->M_dataNavierStokes->density(),
                          M_elementVelocity,
                          *M_elementMatrixConvective,
                          this->M_velocityFESpace.fe(),
                          alpha );
            */

            source_press( 1.0,
                          M_elementVelocity,
                          *elementMatrixPressure,
                          this->M_velocityFESpace.fe(),
                          this->M_pressureFESpace.fe(),
                          (ID) mmFESpace.fe().nbNode );

            //derivative of the convective term
            if ( convectiveTermDerivative )
                mass_gradu( this->M_dataNavierStokes->density(),
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
                                    offset + jComponent * velocityTotalDof );

                    //assembling the derivative of the convective term
                    if ( convectiveTermDerivative )
                        assembleMatrix( matrix,
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
                                offset + iComponent * velocityTotalDof );
            }
        }

    }
    chrono.stop();
    this->M_Displayer.leaderPrintMax("done in ", chrono.diff() );
}

} // namespace LifeV

#endif // OSSENSHAPEDERIVATIVE_H
