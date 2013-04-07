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
 *  @file
 *  @brief This file contains solvers for different materials.
 *  @warning: This is the most important issue related with this class.
 *  At the moment, the BC are applied on the matrix and on rhsNoBc for VK models
 *  but for NH and EXP they are applied on the residual directly. This does
 *  not work for nonhomogeneus Dirichlet conditions!!
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
*/

#ifndef _STRUCTURALOPERATOR_H_
#define _STRUCTURALOPERATOR_H_ 1

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/core/util/LifeChrono.hpp>
#include <lifev/core/util/Displayer.hpp>

#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/structure/fem/AssemblyElementalStructure.hpp>
#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/BCManage.hpp>
#include <lifev/core/fem/FESpace.hpp>

#include <lifev/core/mesh/MeshEntityContainer.hpp>

#include <lifev/core/algorithm/NonLinearRichardson.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

#ifdef COMPUTATION_JACOBIAN
#include <Epetra_SerialDenseMatrix.h>
#endif

//Linear Solver includes
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>


//ET includes
#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/eta/expression/Integrate.hpp>

// Time Advance includes
#include <lifev/core/fem/TimeAdvance.hpp>
#include <lifev/core/fem/TimeAdvanceNewmark.hpp>
#include <lifev/core/fem/TimeAdvanceBDF.hpp>

namespace LifeV
{

using namespace ExpressionAssembly;

template < typename MeshEntityType,
         typename ComparisonPolicyType = boost::function2 < bool,
         const UInt,
         const UInt > >
class MarkerSelector
{
public:
    typedef MeshEntityType       meshEntity_Type;
    typedef ComparisonPolicyType comparisonPolicy_Type;

    MarkerSelector ( const UInt materialFlagReference,
                     comparisonPolicy_Type const& policy = std::equal_to<UInt>() )
        : M_reference ( materialFlagReference ),
          M_policy ( policy ) {}

    bool operator() ( const meshEntity_Type& entity ) const
    {
        //Extract the flag from the mesh entity
        UInt flagChecked = entity.markerID();

        return M_policy ( flagChecked, M_reference );
    }

private:
    const UInt M_reference;
    const comparisonPolicy_Type M_policy;

}; // Marker selector

/*!
  \class StructuralSolver
  \brief
  This class solves the linear elastodynamics equations for different kinds of materials
  (St. Venant-Kirchoff materials right now)

*/
template <typename Mesh>
class StructuralOperator
{
public:

    //!@name Type definitions
    //@{
    typedef Real ( *function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );

    typedef boost::function<Real ( Real const&, Real const&, Real const&, Real const&, ID const& ) > source_Type;

    typedef StructuralConstitutiveLaw<Mesh>               material_Type;
    typedef boost::shared_ptr<material_Type>              materialPtr_Type;

    typedef BCHandler                                     bcHandlerRaw_Type;
    typedef boost::shared_ptr<bcHandlerRaw_Type>          bcHandler_Type;

    typedef LinearSolver                                  solver_Type;

    typedef typename solver_Type::matrix_Type             matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                matrixPtr_Type;
    typedef typename solver_Type::vector_Type             vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;
    typedef vector_Type                                   solution_Type;
    typedef boost::shared_ptr<solution_Type>              solutionPtr_Type;

    typedef StructuralConstitutiveLawData                 data_Type;

    typedef RegionMesh<LinearTetra >                      mesh_Type;
    typedef std::vector< mesh_Type::element_Type* >       vectorVolumes_Type;
    typedef std::vector< UInt >                           vectorIndexes_Type;

    typedef std::map< UInt, vectorVolumes_Type>           mapMarkerVolumes_Type;
    typedef std::map< UInt, vectorIndexes_Type>           mapMarkerIndexes_Type;
    typedef boost::shared_ptr<mapMarkerVolumes_Type>      mapMarkerVolumesPtr_Type;
    typedef boost::shared_ptr<mapMarkerIndexes_Type>      mapMarkerIndexesPtr_Type;
    typedef mapMarkerVolumes_Type::const_iterator         mapIterator_Type;

    typedef typename mesh_Type::element_Type              meshEntity_Type;

    typedef typename boost::function2<bool, const UInt, const UInt> comparisonPolicy_Type;

    typedef MarkerSelector<meshEntity_Type, comparisonPolicy_Type> markerSelector_Type;
    typedef boost::scoped_ptr<markerSelector_Type>          markerSelectorPtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >  ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>                      ETFESpacePtr_Type;

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >          FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;

    //Preconditioners typedef
    typedef LifeV::Preconditioner                   basePrec_Type;
    typedef boost::shared_ptr<basePrec_Type>        basePrecPtr_Type;
    typedef LifeV::PreconditionerIfpack             precIfpack_Type;
    typedef boost::shared_ptr<precIfpack_Type>      precIfpackPtr_Type;
    typedef LifeV::PreconditionerML                 precML_Type;
    typedef boost::shared_ptr<precML_Type>          precMLPtr_Type;

    // Time advance
    typedef TimeAdvance< vector_Type >                                  timeAdvance_Type;
    typedef boost::shared_ptr< timeAdvance_Type >                       timeAdvancePtr_Type;

    //@}


#ifdef COMPUTATION_JACOBIAN
    typedef Epetra_SerialDenseMatrix                     matrixSerialDense_Type;
    typedef boost::shared_ptr<matrixSerialDense_Type>    matrixSerialDensePtr_Type;
    typedef std::vector<LifeV::Real>                     vectorInvariants_Type;
    typedef boost::shared_ptr<vectorInvariants_Type>     vectorInvariantsPtr_Type;
#endif
    //@}

    //! @name Constructor & Destructor
    //@{

    StructuralOperator();

    virtual ~StructuralOperator() {};

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class Venantkirchhof
    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param BCh boundary conditions for the displacement
      \param comm the Epetra Comunicator
    */
    void setup ( boost::shared_ptr<data_Type>  data,
                 const FESpacePtr_Type&        dFESpace,
                 const ETFESpacePtr_Type&      dETFESpace,
                 bcHandler_Type&       BCh,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param comm the Epetra Comunicator
    */
    void setup ( boost::shared_ptr<data_Type> data,
                 const FESpacePtr_Type&       dFESpace,
                 const ETFESpacePtr_Type&     dETFESpace,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param comm the comunicator parameter
      \param monolithicMap the MapEpetra
      \param offset the offset parameter
    */
    void setup ( boost::shared_ptr<data_Type> data,
                 const FESpacePtr_Type&       dFESpace,
                 const ETFESpacePtr_Type&     dETFESpace,
                 boost::shared_ptr<Epetra_Comm>&     comm,
                 const boost::shared_ptr<const MapEpetra>&       monolithicMap,
                 UInt       offset = 0
               );


    //! Updates the system at the end of each time step
    void updateSystem ( void );

    //! Updates the system at the end of each time step when the matrix is passed from outside
    /*!
      \param stiff stiffness matrix provided from outside
    */
    void updateSystem ( matrixPtr_Type& mat_stiff);

    //! Updates the system at the end of each time step given a source term
    /*!
      \param source volumic source
      \param time present time
    */
    void updateSystem ( source_Type const& source );


    //! Updates the system at the end of each time step given a source term
    /*!
      \param source volumic source
      \param time present time
    */
    void updateSourceTerm ( source_Type const& source );

    //! Updates the rhs at the start of each time step
    /*!
      \param rhs: solid  right hand side
      !*/
    void setRightHandSide (const vector_Type& rightHandSide)
    {
        *M_rhsNoBC = rightHandSide;
    };

    //! Comuptes the right hand side in the updateSystem methods
    void computeRHSNoBC ( void );

    //! Compute the mass matrix and it calls the method to build the linear part of the stiffness matrix of the material class
    void buildSystem ( const Real coefficient );

    //void buildSystem(matrix_Type & bigMatrixStokes, const Real& timeAdvanceCoefficient, const Real& factor); // used for monolithic

    //! Compute the mass matrix and the linear part of the stiffness matrix
    /*!
      \param matrix the matrix containing the mass matrix and the linear part of he stiffness matrix
      \param rescale factor for FSI problems
    */
    void computeMassMatrix ( const Real factor = 1.);

    //! Solve the non-linear system
    /*!
      \param bch BCHander object containing the boundary conditions
    */
    void iterate ( const bcHandler_Type& bch );

    //! Solve the linearized problem. Used in FSI segregated in ExactJacobian
    /*!
      \param bch BCHander object containing the boundary conditions
    */
    void iterateLin ( bcHandler_Type& bch );


    //! Output
    /*!
      \param c output file
    */
    void showMe ( std::ostream& c = std::cout ) const;

    //! Update the Jacobian Matrix at each iteration of the nonLinearRichardson method
    /*!
      \param solution the current solution at each iteration of the nonLinearRichardson method
      \param jacobian the Jacobian matrix that must be updated
    */
    void updateJacobian ( const vector_Type& solution, matrixPtr_Type& jacobian );

    //! Update the Jacobian Matrix at each iteration of the nonLinearRichardson method
    /*! Note: this method is used in FSIExactJacobian
      \param solution the current solution at each iteration of the nonLinearRichardson method
    */
    void updateJacobian (const vector_Type& solution );

    //! Solves the tangent problem for newton iterations
    /*!
      \param step the vector containing the solution of the sistem J*step=-Res
      \param res the vector conteining the residual
      \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
      eta is determined by the modified Eisenstat-Walker formula
    */
    void solveJac ( vector_Type&       step,
                    const vector_Type& residual,
                    Real&            linear_rel_tol ) ;
    //    void solveJac( const Vector& res, Real& linear_rel_tol, Vector &step);

    //! Solves the tangent problem with custom BC
    /*!
      \param step the vector containing the solution of the sistem J*step=-Res
      \param res the vector conteining the residual
      \param lin_res_tol linear_rel_tol send for the relative tolerance to the linear solver is therefore eta.
      eta is determined by the modified Eisenstat-Walker formula
      \param BCd BCHandler object containing the boundary condition
    */
    void solveJacobian ( vector_Type&       step,
                         const vector_Type& residual,
                         Real&            linear_rel_tol,
                         bcHandler_Type&    BCd ) ;


    //! Evaluates residual for newton interations
    /*!
      \param res residal vector that is update every time the method is called
      \param sol solution vector from which the residual is computed
      \param iter iteration of the nonLinearRichardson method
    */
    void evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter);

    //! Evaluates residual of the displacement for FSI problems
    /*!
      \param sol, the current displacement of he sturcture
    */
    void evalResidualDisplacement ( const vector_Type& solution );

    //! Evaluates residual of the displacement in the Linearized problem of ExactJcobian. FSI problems
    /*!
      \param sol, the current displacement of he sturcture
    */
    void evalResidualDisplacementLin ( const vector_Type& solution );

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 space function describing the initial velocity
      \param a0 space function describing the initial acceleration
    */
    void initialize ( const function& d0 );

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 empty vector
      \param a0 empty vector
    */
    void initialize ( vectorPtr_Type d0 );

    //! Computes the velocity and acceleration vector at the n-th time step
    //void updateVelAndAcceleration();

    //! Reduce the complete solution to the solution on the pocessor with rank 0
    /*!
      \param disp displacement solution
      \param vel velocity solution
    */
    void reduceSolution ( Vector& displacement, Vector& velocity );

    //! Multiply the mass matrix and the linear stiffness matrix by the rescaleFactor
    //  void rescaleMatrices(); // used for monolithic

    /**
       in the linear case the solid matrix is constant, thus it does not need to be recomputed.
    */

    //! Update (in the case of nonlinear material) the solid matrix
    /*!
      \param stiff stiffness matrix
      \param sol the current solution
      \param factor the rescaleFactor
    */
    void computeMatrix ( matrixPtr_Type& stiff, const vector_Type& sol, Real const& factor );


#ifdef COMPUTATION_JACOBIAN
    //! compute the value of the determinant of F in all the volumes of the mesh
    /*!
      \param displacement the solution at a certain time
      \return the vector with the values for J
    */
    void jacobianDistribution ( vectorPtr_Type displacement, vector_Type& jacobianDistribution );
#endif


#ifdef COLORING_MESH
    //! compute the value of the determinant of F in all the volumes of the mesh
    /*!
      \param displacement the solution at a certain time
      \return the vector with the values for J
    */
    void colorMesh ( vector_Type& meshColors );
#endif

    //void updateMatrix(matrix_Type & bigMatrixStokes);// used for monolithic
    //void updateCoupling(matrix_Type couplingMatrix);// used for monolithic

    //@}

    //! @name Set Methods
    //@{

    //!Setters
    //! Set the BCHandler object
    void setBC (const bcHandler_Type& BCd)
    {
        M_BCh = BCd;
    }

    //! Set the source object
    void setSourceTerm ( source_Type const& s )
    {
        M_source = s;
    }

    // //! Set the preconditioner
    // void resetPrec(bool reset = true) { if (reset) M_linearSolver.precReset(); }

    // //! Set the displacement
    // virtual void setDisp(const vector_Type& disp) {*M_disp = disp;} // used for monolithic

    //! Set the recur parameter
    void setRecur (UInt recur)
    {
        M_recur = recur;
    }

    //! Set the data fields with the Getpot data file for preconditioners and solver
    void setDataFromGetPot ( const GetPot& dataFile );

    void setTimeAdvance( const timeAdvancePtr_Type& timeAdvancePtr )
    {
        M_timeAdvance = timeAdvancePtr;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getters
    //! Get the Epetramap
    MapEpetra   const& map()       const
    {
        return *M_localMap;
    }

    //! Get the Displayer object
    Displayer   const& displayer() const
    {
        return *M_Displayer;
    }

    boost::shared_ptr<const Displayer>   const& displayerPtr() const
    {
        return M_Displayer;
    }

    //! Get the matrix containing the mass mtrix and the linear part of the stiffness matrix
    //matrixPtr_Type const MassStiff() const {return M_massStiff; }

    //! Get the mass matrix
    matrixPtr_Type const massMatrix() const
    {
        return M_massMatrix;
    }

    //! Get the FESpace object
    FESpace_Type& dispFESpace()
    {
        return *M_dispFESpace;
    }

    //! Get the ETFESpace object
    ETFESpace_Type& dispETFESpace()
    {
        return *M_dispETFESpace;
    }

    //! Get the bCHandler object
    bcHandler_Type const& bcHandler() const
    {
        return M_BCh;
    }

    //! Get the residual
    vector_Type& residual()
    {
        return *M_residual_d;
    }

    //! Get the source term
    source_Type const& sourceTerm() const
    {
        return M_source;
    }

    //! Get the displacement
    vector_Type& displacement()
    {
        return *M_disp;
    }

    vector_Type& displacementPtr()
    {
        return M_disp;
    }

    //! Get the right hand sde without BC
    vectorPtr_Type& rhsWithoutBC()
    {
        return M_rhsNoBC;
    }

    //! Get the right hand. The member rhsCopy is used for Debug purposes!
    vector_Type& rhsCopy()
    {
        return *M_rhsCopy;
    }
    vector_Type& residualCopy()
    {
        return *M_residualCopy;
    }

    //! Get the comunicator object
    boost::shared_ptr<Epetra_Comm> const& comunicator() const
    {
        return M_Displayer->comm();
    }

    //! Get the rescaleFactor
    Real rescaleFactor()
    {
        return M_rescaleFactor;
    }

    /*! Get the offset parameter. It is taken into account when the boundary conditions
      are applied and the matrices are assembled.
    */
    const UInt& offset() const
    {
        return M_offset;
    }

    /*! Get the offset parameter. It is taken into account when the boundary conditions
      are applied and the matrices are assembled.
    */
    const materialPtr_Type& material() const
    {
        return M_material;
    }

    /**
       Do nothing in the linear case: the matrix remains constant. Otherwise substitute the matrix with an updated one
    */
    //! Get the Solid Matrix
    void solidMatrix ( matrixPtr_Type& /*matrix*/ )
    {
    }

    // Physic constant
    //! Get the thickness
    Real thickness() const
    {
        return M_data->thickness();
    }

    //! Get the Young modulus
    Real young ( UInt material = 1)            const
    {
        return M_data->young ( material );
    }

    //! Get the Poisson coefficient
    Real poisson ( UInt material = 1 )          const
    {
        return M_data->poisson ( material );
    }

    //! Get the density
    Real rho()       const
    {
        return M_data->rho();
    }

    //! Get the data container
    const boost::shared_ptr<data_Type>& data() const
    {
        return M_data;
    }

    void apply ( const vector_Type& sol, vector_Type& res) const;

    //! Get the density
    mapMarkerVolumesPtr_Type mapMarkersVolumes() const
    {
        return M_mapMarkersVolumes;
    }

    //! Get the density
    mapMarkerIndexesPtr_Type mapMarkersIndexes() const
    {
        return M_mapMarkersIndexes;
    }

    const timeAdvancePtr_Type& timeAdvancePtr() const
    {
        return M_timeAdvance;
    }

    //@}

protected:

    //! Apply boundary condition
    /*!
      \param matrix the matrix of the system
      \param rhs the right hand side of the system
      \param BCh BCHandler object
      \param offset the offset parameter
    */
    void applyBoundaryConditions (matrix_Type& matrix,
                                  vector_Type& rhs,
                                  bcHandler_Type& BCh,
                                  UInt         offset = 0);


    UInt dim() const
    {
        return M_dispFESpace->dim();
    }


    //! construct the map between the markers and the volumes
    /*!
      \param VOID
      \return VOID
    */
    void setupMapMarkersVolumes ( void );

    //!Protected Members

#ifdef COMPUTATION_JACOBIAN
    //! constructPatchAreaVector: This method build the patch area vector used in the reconstruction process
    /*!
      \param NONE
    */
    void constructPatchAreaVector ( vector_Type& patchArea, const vector_Type& solution );


    //! reconstructElementaryVector: This method applies a reconstruction procedure on the elvec that is passed
    /*!
      \param elvecTens VectorElemental over which the reconstruction is applied
    */
    void reconstructElementaryVector ( VectorElemental& elVecSigma, vector_Type& patchArea, UInt nVol );
#endif


    boost::shared_ptr<data_Type>         M_data;

    FESpacePtr_Type                      M_dispFESpace;

    ETFESpacePtr_Type                    M_dispETFESpace;

    boost::shared_ptr<const Displayer>   M_Displayer;

    Int                                  M_me;

    //! data for solving tangent problem with aztec + preconditioner
    boost::shared_ptr<solver_Type>       M_linearSolver;
    basePrecPtr_Type                     M_preconditioner;

    //! Elementary matrices and vectors
    boost::shared_ptr<MatrixElemental>   M_elmatM;

    //! linearized velocity
    vectorPtr_Type                       M_disp;

    //! right  hand  side displacement
    vectorPtr_Type                       M_rhs;
    vectorPtr_Type                       M_rhsCopy;
    vectorPtr_Type                       M_residualCopy;

    //! right  hand  side
    vectorPtr_Type                       M_rhsNoBC;

    //! right  hand  side
    //boost::shared_ptr<vector_Type>       M_f;

    //! residual
    boost::shared_ptr<vector_Type>       M_residual_d;

    //! files for lists of iterations and residuals per timestep
    std::ofstream                        M_out_iter;
    std::ofstream                        M_out_res;

    //! BCHandler object
    bcHandler_Type                       M_BCh;

    //! Map Epetra
    boost::shared_ptr<const MapEpetra>   M_localMap;

    //! Matrix M: mass
    matrixPtr_Type                       M_massMatrix;

    //! Matrix Temp: Temporary matrix to compute residuals or rhs
    matrixPtr_Type                       M_systemMatrix;


    //! Jacobian Matrix: Matrix to store the jacobian of the newton method
    matrixPtr_Type                       M_jacobian;


    //! level of recursion for Aztec (has a sens with FSI coupling)
    UInt                                 M_recur;

    source_Type                          M_source;

    UInt                                 M_offset;
    Real                                 M_rescaleFactor;
    //  Real                                 M_zeta;
    //  Real                                 M_theta;

    //! Material class
    materialPtr_Type                     M_material;

    //! Map between markers and volumes on the mesh
    mapMarkerVolumesPtr_Type             M_mapMarkersVolumes;

    //! Map between markers and volumes on the mesh
    mapMarkerIndexesPtr_Type             M_mapMarkersIndexes;

#ifdef COMPUTATION_JACOBIAN
    //! Elementary matrix for the tensor F
    matrixSerialDensePtr_Type            M_deformationF;
    vectorInvariants_Type                M_invariants;
#endif

    timeAdvancePtr_Type                  M_timeAdvance;
};

//====================================
// Constructor
//=====================================

template <typename Mesh>
StructuralOperator<Mesh>::StructuralOperator( ) :
    M_data                       ( ),
    M_dispFESpace                    ( ),
    M_dispETFESpace                  ( ),
    M_Displayer                  ( ),
    M_me                         ( 0 ),
    M_linearSolver               ( ),
    M_preconditioner             ( ),
    M_elmatM                     ( ),
    M_disp                       ( ),
    M_rhsNoBC                    ( ),
    M_rhsCopy                    ( ),
    M_residualCopy               ( ),
    M_residual_d                 ( ),
    M_out_iter                   ( ),
    M_out_res                    ( ),
    M_BCh                        ( ),
    M_localMap                   ( ),
    M_massMatrix                 ( ),
    M_systemMatrix               ( ),
    M_jacobian                   ( ),
    M_recur                      ( ),
    M_source                     ( ),
    M_offset                     ( 0 ),
    M_rescaleFactor              ( 1. ),
    M_material                   ( ),
#ifdef COMPUTATION_JACOBIAN
    M_deformationF               ( ),
    M_invariants                 ( ),
#endif
    M_mapMarkersVolumes          ( ),
    M_mapMarkersIndexes          ( ),
    M_timeAdvance                ( )
{

    //    M_Displayer->leaderPrint("I am in the constructor for the solver");
}

template <typename Mesh>
void
StructuralOperator<Mesh>::setup (boost::shared_ptr<data_Type>          data,
                                 const FESpacePtr_Type& dFESpace,
                                 const ETFESpacePtr_Type& dETFESpace,
                                 bcHandler_Type&                    BCh,
                                 boost::shared_ptr<Epetra_Comm>&   comm)
{
    setup (data, dFESpace, dETFESpace, comm);
    M_BCh = BCh;
}

template <typename Mesh>
void
StructuralOperator<Mesh>::setup (boost::shared_ptr<data_Type>        data,
                                 const FESpacePtr_Type& dFESpace,
                                 const ETFESpacePtr_Type& dETFESpace,
                                 boost::shared_ptr<Epetra_Comm>&     comm)
{
    setup ( data, dFESpace, dETFESpace, comm, dFESpace->mapPtr(), (UInt) 0 );

    M_rhs.reset                        ( new vector_Type (*M_localMap) );
    M_rhsCopy.reset                    ( new vector_Type (*M_localMap) );
    M_residualCopy.reset               ( new vector_Type (*M_localMap) );
    M_rhsNoBC.reset                    ( new vector_Type (*M_localMap) );
    M_linearSolver.reset               ( new LinearSolver ( comm ) );
    M_disp.reset                       ( new vector_Type (*M_localMap) );
}

template <typename Mesh>
void
StructuralOperator<Mesh>::setup (boost::shared_ptr<data_Type>        data,
                                 const FESpacePtr_Type& dFESpace,
                                 const ETFESpacePtr_Type& dETFESpace,
                                 boost::shared_ptr<Epetra_Comm>&     comm,
                                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                 UInt                                offset)
{
    M_data                            = data;
    M_dispFESpace                     = dFESpace;
    M_dispETFESpace                   = dETFESpace;
    M_Displayer.reset                 (new Displayer (comm) );
    M_me                              = comm->MyPID();
    M_elmatM.reset                    ( new MatrixElemental ( M_dispFESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_localMap                        = monolithicMap;
    M_massMatrix.reset                (new matrix_Type (*M_localMap) );
    M_systemMatrix.reset              (new matrix_Type (*M_localMap) );
    M_jacobian.reset                  (new matrix_Type (*M_localMap) );

    M_offset                          = offset;

    M_material.reset ( material_Type::StructureMaterialFactory::instance().createObject ( M_data->solidType() ) );
    M_material->setup ( dFESpace, dETFESpace, M_localMap, M_offset, M_data, M_Displayer );

    if ( M_data->verbose() )
    {
        M_out_iter.open ( "out_iter_solid" );
        M_out_res.open ( "out_res_solid" );
    }
    M_mapMarkersVolumes.reset ( new mapMarkerVolumes_Type() );
    M_mapMarkersIndexes.reset ( new mapMarkerIndexes_Type() );
    //this->setupMapMarkersVolumes();
}


template <typename Mesh>
void StructuralOperator<Mesh>::setupMapMarkersVolumes ( void )
{

    LifeChrono time;

    this->M_Displayer->leaderPrint (" S-  Starting the time:  \n");

    time.start();


    this->M_Displayer->leaderPrint (" S-  Building the map between volumesMarkers <--> volumes \n");

    //We first loop over the vector of the material_flags
    for (  UInt i (0); i < M_data->vectorFlags().size(); i++ )
    {

        //Create the functor to extract volumes
        markerSelectorPtr_Type ref ( new markerSelector_Type (M_data->vectorFlags() [i]) );

        //Number of volumes with the current marker
        UInt numExtractedVolumes = this->M_dispFESpace->mesh()->elementList().countAccordingToPredicate ( *ref );

        this->M_Displayer->leaderPrint (" Current marker: ", M_data->vectorFlags() [i]);
        this->M_Displayer->leaderPrint (" \n");
        this->M_Displayer->leaderPrint (" Number of volumes:", numExtractedVolumes);
        this->M_Displayer->leaderPrint (" \n");

        //Vector large enough to contain the number of volumes with the current marker
        vectorVolumes_Type extractedVolumes ( numExtractedVolumes );
        vectorIndexes_Type extractedIndexes;

        //Extracting the volumes and the corresponding position
        extractedVolumes = this->M_dispFESpace->mesh()->elementList().extractAccordingToPredicateNonConstElement ( *ref, extractedIndexes );

        //Insert the correspondande Marker <--> List of Volumes inside the map
        M_mapMarkersVolumes->insert ( pair<UInt, vectorVolumes_Type> (M_data->vectorFlags() [i], extractedVolumes) ) ;
        M_mapMarkersIndexes->insert ( pair<UInt, vectorIndexes_Type> (M_data->vectorFlags() [i], extractedIndexes) ) ;

        // for( UInt i(0); i<extractedIndexes.size(); i++ )
        //     std::cout << "Element: " << extractedIndexes[i] << std::endl;

        //Cleaning the vector
        extractedVolumes.clear();
        extractedIndexes.clear();

    }

    time.stop();

    this->M_Displayer->leaderPrint (" S-  Time to build the map:", time.diff() );
}

template <typename Mesh>
void StructuralOperator<Mesh>::updateSystem ( void )
{
    updateSystem (M_systemMatrix);
}

template <typename Mesh>
void StructuralOperator<Mesh>::updateSystem ( matrixPtr_Type& mat_stiff)
{
    M_Displayer->leaderPrint (" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix
    M_material->computeStiffness (*M_disp, M_rescaleFactor, M_data, M_mapMarkersVolumes, M_mapMarkersIndexes, M_Displayer);

    if ( M_data->solidType() == "linearVenantKirchhoff" )
    {
        *mat_stiff += *M_material->stiffMatrix();
        mat_stiff->globalAssemble();
    }


    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );

}

template <typename Mesh>
void StructuralOperator<Mesh>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs (vector_Type (*M_localMap) );

    VectorElemental M_elvec (M_dispFESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= M_dispFESpace->mesh()->numVolumes(); ++i )
    {

        M_dispFESpace->fe().updateFirstDerivQuadPt ( M_dispFESpace->mesh()->volumeList ( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            //compute_vec( source, M_elvec, M_dispFESpace->fe(),  M_data->dataTime()->time(), ic ); // compute local vector
            assembleVector ( *rhs, M_elvec, M_dispFESpace->fe(), M_dispFESpace->dof(), ic, ic * M_dispFESpace->fieldDim() ); // assemble local vector into global one
        }
    }
    M_rhsNoBC += rhs;
}

template <typename Mesh>
void StructuralOperator<Mesh>::buildSystem ( const Real coefficient )
{
    M_Displayer->leaderPrint ("  S-  Computing constant matrices ...          ");
    LifeChrono chrono;
    chrono.start();

    computeMassMatrix ( coefficient );
    M_material->computeLinearStiff (M_data, M_mapMarkersVolumes, M_mapMarkersIndexes);

    chrono.stop();
    M_Displayer->leaderPrintMax ( "done in ", chrono.diff() );

}

template <typename Mesh>
void
StructuralOperator<Mesh>::computeMassMatrix ( const Real factor)
{
    using namespace ExpressionAssembly;

    //UInt totalDof = M_dispFESpace->dof().numTotalDof();

    //! Number of displacement components
    /*UInt nc = nDimensions;*/
    const Real factorMassMatrix = factor * M_data->rho();

    //Assembling using the Expression Template
    //At the moment, that the ET is not the standard, when the integration is
    //performed the quadrature rule set in FESpace is used.
    //Otherwhise, one can set his preferred quadrature rule for the integrals.

    integrate ( elements (M_dispETFESpace->mesh() ),
                M_dispFESpace->qr(),
                M_dispETFESpace,
                M_dispETFESpace,
                value (factorMassMatrix) *  dot ( phi_i , phi_j ) ) >> M_massMatrix;

    M_massMatrix->globalAssemble();

    //M_massMatrix->spy("massMatrixStructure.m");

    //*massStiff *= factor; //M_data.dataTime()->timeStep() * M_rescaleFactor;
}

template <typename Mesh>
void
StructuralOperator<Mesh>::iterate ( const bcHandler_Type& bch )
{
    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer->leaderPrint ("  S-  Solving the system ... \n");

    M_BCh = bch;

    Real abstol  = 1.e-7;
    Real reltol  = 1.e-7;
    UInt maxiter = 200;
    Real etamax  = 1e-7;
    Int NonLinearLineSearch = 0;

    Real time = M_data->dataTime()->time();

    Int status = 0;

    if ( M_data->verbose() )
    {
        status = NonLinearRichardson ( *M_disp, *this, abstol, reltol, maxiter, etamax, NonLinearLineSearch, 0, 2, M_out_res, M_data->dataTime()->time() );
    }
    else
    {
        status = NonLinearRichardson ( *M_disp, *this, abstol, reltol, maxiter, etamax, NonLinearLineSearch );
    }


    if ( status == 1 )
    {
        std::ostringstream ex;
        ex << "StructuralOperator::iterate() Inners nonLinearRichardson iterations failed to converge\n";
        throw std::logic_error ( ex.str() );
    }
    else // if status == 0 NonLinearrRichardson converges
    {
        // std::cout << std::endl;

        // std::cout <<" Number of inner iterations       : " << maxiter <<  std::endl;

        // std::cout <<" We are at the time step          : "  << M_data->dataTime()->time() << std::endl;
        if ( M_data->verbose() )
        {
            M_out_iter << time << " " << maxiter << std::endl;
        }
    }

    // std::cout << "iterate: d norm       = " << M_disp->norm2() << std::endl;

    //These two lines mut be checked fo FSI. With the linear solver, they have a totally
    //different expression. For structural problems it is not used.
    evalResidualDisplacement (*M_disp); //\todo related to FSI. Should be caled from outside
}

template <typename Mesh>
void
StructuralOperator<Mesh>::iterateLin ( bcHandler_Type& bch )
{
    vector_Type rhsFull (M_rhsNoBC->map() );
    Real zero (0.);
    if ( !bch->bcUpdateDone() )
    {
        bch->bcUpdate ( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
    }
    bcManageVector ( rhsFull, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *bch, M_dispFESpace->feBd(), M_data->dataTime()->time(), 1.0 );
    solveJacobian (*M_disp, rhsFull, zero, bch);
    evalResidualDisplacementLin (*M_disp);
}


template <typename Mesh>
void
StructuralOperator<Mesh>::showMe ( std::ostream& c  ) const
{
    c << "\n*** StructuralOperator::showMe method" << std::endl;

    M_data->showMe ( c );

}

template <typename Mesh>
void StructuralOperator<Mesh>::computeMatrix ( matrixPtr_Type& stiff, const vector_Type& sol,  Real const& /*factor*/)
{
    M_Displayer->leaderPrint ( " Computing residual ... \t\t\t");

    LifeChrono chrono;
    chrono.start();

    //! It is right to do globalAssemble() inside the M_material class
    M_material->computeStiffness ( sol, 1., M_data, M_mapMarkersVolumes, M_mapMarkersIndexes, M_Displayer);

    if ( M_data->solidType() == "linearVenantKirchhoff" )
    {
        *stiff = *M_material->stiffMatrix();
        *stiff += *M_massMatrix;
        stiff->globalAssemble();
    }

    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}

#ifdef COMPUTATION_JACOBIAN

template <typename Mesh>
void StructuralOperator<Mesh>::jacobianDistribution ( vectorPtr_Type displacement, vector_Type& jacobianDistribution )
{
    M_Displayer->leaderPrint ( " Computing the jacobian for all the volumes ... \t\t\t");

    //Initialization of the deformationF matrix
    M_deformationF.reset  ( new matrixSerialDense_Type ( M_dispFESpace->fieldDim(), M_dispFESpace->fieldDim() ) );

    vector_Type vectorJacobian ( jacobianDistribution );
    vectorJacobian *= 0.0;

    LifeChrono chrono;
    chrono.start();

    //construct a vector to store the are
    vector_Type patchArea (*displacement, Unique, Add);
    patchArea *= 0.0;

    constructPatchAreaVector ( patchArea, *displacement );

    //Before assembling the reconstruction process is done
    vector_Type patchAreaR (patchArea, Repeated);


    //Loop over the volumes to compute J = det(F)
    //Inside the loop, the determinant is store in the appropriate positions
    vector_Type dRep (*displacement, Repeated);
    UInt totalDof = M_dispFESpace->dof().numTotalDof();
    VectorElemental dk_loc ( M_dispFESpace->fe().nbFEDof(), this->M_dispFESpace->fieldDim() );
    VectorElemental elVecDet ( M_dispFESpace->fe().nbFEDof(), this->M_dispFESpace->fieldDim() );

    //Building fake quadrature rules to compute the deformation gradient F at the nodes
    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < M_dispFESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += M_dispFESpace->qr().weight (iq);
    }

    Real wQuad (refElemArea / M_dispFESpace->refFE().nbDof() );

    //Setting the quadrature Points = DOFs of the element and weight = 1
    std::vector<GeoVector> coords = M_dispFESpace->refFE().refCoor();
    std::vector<Real> weights (M_dispFESpace->fe().nbFEDof(), wQuad);
    fakeQuadratureRule.setDimensionShape ( shapeDimension (M_dispFESpace->refFE().shape() ), M_dispFESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //Set the new quadrature rule
    M_dispFESpace->setQuadRule (fakeQuadratureRule);


    //Loop over the volumes. No markerIDs are necessary on this loop for J = det(F)
    for ( UInt i = 0; i < M_dispFESpace->mesh()->numVolumes(); ++i )
    {

        //Vectors for the deformation tensor
        std::vector<matrixSerialDense_Type> vectorDeformationF (M_dispFESpace->fe().nbFEDof(), *M_deformationF);
        M_invariants.resize   ( M_dispFESpace->fieldDim() + 1 );


        M_dispFESpace->fe().updateFirstDerivQuadPt ( M_dispFESpace->mesh()->volumeList ( i ) );

        UInt eleID = M_dispFESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) M_dispFESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = M_dispFESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_dispFESpace->fieldDim(); ++iComp )
            {
                UInt ig = M_dispFESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * M_dispFESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * M_dispFESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //computing the tensor F
        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, M_dispFESpace->fe() );


        //Cycle over the nDofs/element
        for ( UInt nDOF = 0; nDOF < ( UInt ) M_dispFESpace->fe().nbFEDof(); nDOF++ )
        {
            UInt  iloc = M_dispFESpace->fe().patternFirst ( nDOF );

            //computing the determinant
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor ( M_invariants, vectorDeformationF[nDOF] );

            //Assembling elemental vector
            (elVecDet) [ iloc ] = M_invariants[3];
            (elVecDet) [ iloc + M_dispFESpace->fe().nbFEDof() ] = 0.0;
            (elVecDet) [ iloc + 2 * M_dispFESpace->fe().nbFEDof() ] = 0.0;

        }

        //multiplying it for the patch area
        reconstructElementaryVector ( elVecDet, patchAreaR, i );

        //assembling it into the global vector
        for ( UInt ic = 0; ic < this->M_dispFESpace->fieldDim(); ++ic )
        {
            assembleVector (vectorJacobian, elVecDet, M_dispFESpace->fe(), M_dispFESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }

        vectorDeformationF.clear();
        M_invariants.clear();

    }

    vectorJacobian.globalAssemble();

    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );

    jacobianDistribution = vectorJacobian;
}


template <typename Mesh>
void StructuralOperator<Mesh >::constructPatchAreaVector ( vector_Type& patchArea,
                                                           const vector_Type& solution )
{

    vector_Type patchAreaR (solution, Repeated);
    patchAreaR *= 0.0;

    Real refElemArea (0); //area of reference element
    UInt totalDof = M_dispFESpace->dof().numTotalDof();
    //compute the area of reference element
    for (UInt iq = 0; iq < M_dispFESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += M_dispFESpace->qr().weight (iq);
    }

    // Define a special quadrature rule for the interpolation
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_dispFESpace->refFE().shape() ), M_dispFESpace->refFE().shape() );
    Real wQuad (refElemArea / M_dispFESpace->refFE().nbDof() );

    for (UInt i (0); i < M_dispFESpace->refFE().nbDof(); ++i) //nbRefCoor
    {
        interpQuad.addPoint (QuadraturePoint (M_dispFESpace->refFE().xi (i), M_dispFESpace->refFE().eta (i), M_dispFESpace->refFE().zeta (i), wQuad) );
    }

    UInt totalNumberVolumes (M_dispFESpace->mesh()->numVolumes() );
    UInt numberLocalDof (M_dispFESpace->dof().numLocalDof() );

    CurrentFE interpCFE (M_dispFESpace->refFE(), getGeometricMap (* (M_dispFESpace->mesh() ) ), interpQuad);

    // Loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberVolumes; iterElement++)
    {
        interpCFE.update (M_dispFESpace->mesh()->volumeList ( iterElement ), UPDATE_WDET );

        for (UInt iterDof (0); iterDof < numberLocalDof; iterDof++)
        {
            for (UInt iDim (0); iDim < M_dispFESpace->fieldDim(); ++iDim)
            {
                ID globalDofID (M_dispFESpace->dof().localToGlobalMap (iterElement, iterDof) + iDim * totalDof);
                patchAreaR[globalDofID] += interpCFE.measure();
            }
        }
    }

    vector_Type final (patchAreaR, Unique, Add);

    patchArea.add (final);

}

template <typename Mesh>
void
StructuralOperator<Mesh >::reconstructElementaryVector ( VectorElemental& elVecDet,
                                                         vector_Type& patchArea,
                                                         UInt nVol )
{
    //UpdateElement Infos
    //M_dispFESpace->fe().updateFirstDerivQuadPt( M_dispFESpace->mesh()->volumeList( nVol ) );

    Real measure = M_dispFESpace->fe().measure();
    UInt eleID = M_dispFESpace->fe().currentLocalId();

    for (UInt iDof = 0; iDof < M_dispFESpace->fe().nbFEDof(); iDof++)
    {
        UInt  iloc = M_dispFESpace->fe().patternFirst ( iDof );

        for ( UInt icoor = 0;  icoor < M_dispFESpace->fieldDim(); icoor++ )
        {
            ID globalDofID (M_dispFESpace->dof().localToGlobalMap (eleID, iDof) + icoor * M_dispFESpace->dof().numTotalDof() );

            elVecDet[iloc + icoor * M_dispFESpace->fe().nbFEDof()] *= ( measure / patchArea[globalDofID] );
        }

    }
}

#endif


#ifdef COLORING_MESH
template <typename Mesh>
void StructuralOperator<Mesh>::colorMesh ( vector_Type& meshColors )
{
    UInt totalDof = this->M_dispFESpace->dof().numTotalDof();

    mapIterator_Type it;

    for ( it = (*M_mapMarkersVolumes).begin(); it != (*M_mapMarkersVolumes).end(); it++ )
    {

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        for ( UInt j (0); j < it->second.size(); j++ )
        {
            this->M_dispFESpace->fe().updateFirstDerivQuadPt ( * (it->second[j]) );

            UInt eleID = this->M_dispFESpace->fe().currentLocalId();

            for ( UInt iNode = 0; iNode < ( UInt ) this->M_dispFESpace->fe().nbFEDof(); iNode++ )
            {
                UInt  iloc = this->M_dispFESpace->fe().patternFirst ( iNode );

                //Extract the global ID of the x-component of the field
                UInt globalIDofDOF = this->M_dispFESpace->dof().localToGlobalMap ( eleID, iloc );

                if ( meshColors.blockMap().LID (globalIDofDOF) != -1 ) // The Global ID is on the calling processors
                {
                    Int LIDid = meshColors.blockMap().LID ( globalIDofDOF );
                    Int GIDid = meshColors.blockMap().GID ( LIDid );
                    meshColors[ GIDid ] = marker;

                }

            }
        }

    }
}

#endif

template <typename Mesh>
void
StructuralOperator<Mesh>::evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter)
{

    //This method call the M_material computeStiffness
    computeMatrix (M_systemMatrix, solution, 1.);


    M_Displayer->leaderPrint ("    S- Updating the boundary conditions ... \t");
    LifeChrono chrono;

    if ( !M_BCh->bcUpdateDone() )
    {
        M_BCh->bcUpdate ( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    if ( M_data->solidType() == "linearVenantKirchhoff" )
    {
        chrono.start();

        matrix_Type matrixFull (*M_systemMatrix);

        if (iter == 0)
        {
            *M_rhs = *M_rhsNoBC;

            bcManageVector ( *M_rhs, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_BCh, M_dispFESpace->feBd(),  M_data->dataTime()->time(), 1.0 );

            //To export for check
            M_rhsCopy = M_rhs;

            // std::string nameFile="residualAfterBC";
            // M_rhs->spy(nameFile);
            // int n;
            // std::cin >> n;
        }

        bcManageMatrix ( matrixFull, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_BCh, M_dispFESpace->feBd(), 1.0 );

        residual  = matrixFull * solution;
        residual -= *M_rhs;
        chrono.stop();
        M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
    }
    else //NH and Exp and SVK VK-Penalized
    {
        chrono.start();
        *M_rhs = *M_rhsNoBC;
        residual = *M_massMatrix * solution;
        residual += *M_material->stiffVector();
        vector_Type solRep (solution, Repeated);
        bcManageResidual ( residual, *M_rhs, solRep, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_BCh, M_dispFESpace->feBd(), M_data->dataTime()->time(), 1.0 );
        residual -= *M_rhs;
        chrono.stop();
        M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
    }

    if ( iter == 0 )
    {
        *M_residualCopy = residual;
        M_rhsCopy = M_rhs;
    }
}

template <typename Mesh>
void
StructuralOperator<Mesh>::evalResidualDisplacement ( const vector_Type& solution )
{

    M_Displayer->leaderPrint ("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    if ( M_data->solidType() == "linearVenantKirchhoff" )
    {
        M_residual_d.reset (new vector_Type ( *M_systemMatrix * solution ) );
        *M_residual_d -= *M_rhsNoBC;
    }
    else //Cases: NH, Exp, VK-Penalized
    {
        M_residual_d.reset (new vector_Type ( *M_material->stiffVector() ) );
        *M_residual_d -= *M_rhsNoBC;
    }
    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}


template <typename Mesh>
void
StructuralOperator<Mesh>::evalResidualDisplacementLin ( const vector_Type& solution )
{

    M_Displayer->leaderPrint ("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    //This definition of residual_d is similar to the one of iterateLin in VenantKirchhoffSolver
    M_residual_d.reset (new vector_Type ( (*M_jacobian) *solution) );

    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}


template <typename Mesh>
void
StructuralOperator<Mesh>::initialize ( vectorPtr_Type disp )
{
    *M_disp = *disp;
}

template <typename Mesh>
void
StructuralOperator<Mesh>::initialize ( const function& d0 )
{
    M_FESpace->interpolate ( static_cast<typename FESpace<Mesh, MapEpetra>::function_Type> ( d0 ), *M_disp, 0.0);
    //M_FESpace->interpolate(w0, *M_vel , 0.0);
    // M_FESpace->interpolate(a0, *M_acc , 0.0);
}

}

template<typename Mesh>
void
StructuralOperator<Mesh>::reduceSolution ( Vector& displacement, Vector& velocity )
{
    vector_Type disp (*M_disp, 0);
    //vector_Type vel(*M_vel , 0);

    if ( comunicator()->MyPID() == 0 )
    {
        for ( UInt iDof = 0; iDof < nDimensions * dim(); ++iDof )
        {
            disp[ iDof ] = displacement[ iDof + 1 ];
        }
    }
}


template <typename Mesh>
void
StructuralOperator<Mesh>::setDataFromGetPot ( const GetPot& dataFile )
{

    M_Displayer->leaderPrint ( "Setting up Preconditioner... \n" );
    //Setting up the preconditioner
    const std::string preconditionerType = dataFile ( "solid/prec/prectype", "Ifpack" );
    const std::string xmlFileName = dataFile ( "solid/prec/xmlName", "xmlParameters.xml" );
    basePrecPtr_Type precPtr; //Abstract class for preconditioners

    if (  ! ( preconditionerType.compare ("Ifpack") ) ) //The preconditioner if Ifpack
    {
        precIfpack_Type* precRawPtr;
        precRawPtr = new precIfpack_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "solid/prec" );

        //Initializing the preconditioner
        M_preconditioner.reset ( precRawPtr );
    }
    else
    {
        precML_Type* precRawPtr;
        precRawPtr = new precML_Type;
        precRawPtr->setDataFromGetPot ( dataFile, "solid/prec" );

        //Initializing the preconditioner
        M_preconditioner.reset ( precRawPtr );
    }


    M_Displayer->leaderPrint ( "Setting up LinearSolver... \n" );

    Teuchos::RCP< Teuchos::ParameterList > paramList = Teuchos::rcp ( new Teuchos::ParameterList );
    paramList = Teuchos::getParametersFromXmlFile ( xmlFileName );


    M_linearSolver->setParameters ( *paramList );
    M_linearSolver->setPreconditioner ( M_preconditioner );
    M_rescaleFactor = dataFile ( "solid/rescaleFactor", 0. );
}


//Method UpdateJacobian
template <typename Mesh>
void StructuralOperator<Mesh>::updateJacobian ( const vector_Type& sol, matrixPtr_Type& jacobian  )
{
    M_Displayer->leaderPrint ("  S-  Solid: Updating JACOBIAN... ");

    LifeChrono chrono;
    chrono.start();

    M_material->updateJacobianMatrix (sol, M_data, M_mapMarkersVolumes, M_mapMarkersIndexes,  M_Displayer);

    jacobian.reset (new matrix_Type (*M_localMap) );
    *jacobian += * (M_material->jacobian() );

    *jacobian += *M_massMatrix;

    jacobian->globalAssemble();

    chrono.stop();
    M_Displayer->leaderPrintMax ("   ... done in ", chrono.diff() );

}

//Method UpdateJacobian
template <typename Mesh>
void StructuralOperator<Mesh>::updateJacobian ( const vector_Type& sol)
{
    updateJacobian (sol, M_jacobian);
}

//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh>
void StructuralOperator<Mesh>::
solveJac ( vector_Type& step, const vector_Type& res, Real& linear_rel_tol)
{
    updateJacobian ( *M_disp, M_jacobian );
    solveJacobian (step,  res, linear_rel_tol, M_BCh);
}


//Method SolveJacobian
template <typename Mesh>
void StructuralOperator<Mesh>::
solveJacobian (vector_Type&           step,
               const vector_Type&     res,
               Real&                 /*linear_rel_tol*/,
               bcHandler_Type&       /* BCh*/)
{
    LifeChrono chrono;

    //Creating pointers to vectors for linear solver
    vectorPtr_Type pointerToRes ( new vector_Type (res) );
    vectorPtr_Type pointerToStep ( new vector_Type (*M_localMap) );

    //Initializing the pointer
    *pointerToStep *= 0.0;

    matrixPtr_Type matrFull (new matrix_Type (*M_localMap) );
    *matrFull += *M_jacobian;

    M_Displayer->leaderPrint ("\tS'-  Solving the linear system ... \n");

    M_Displayer->leaderPrint ("\tS'-  Applying boundary conditions      ... ");


    if ( !M_BCh->bcUpdateDone() )
    {
        M_BCh->bcUpdate ( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
    }
    bcManageMatrix ( *matrFull, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *M_BCh, M_dispFESpace->feBd(), 1.0 );

    M_Displayer->leaderPrintMax ( "done in ", chrono.diff() );

    M_Displayer->leaderPrint ("\tS'-  Solving system                    ... \n");
    chrono.start();

    //Setting up the quantities
    M_linearSolver->setOperator ( matrFull );
    M_linearSolver->setRightHandSide ( pointerToRes );

    //Solving the system
    M_linearSolver->solve ( pointerToStep );

    step = *pointerToStep;

    chrono.stop();
}

template<typename Mesh>
void StructuralOperator<Mesh>::apply ( const vector_Type& sol, vector_Type& res) const
{
    M_material->apply (sol, res, M_mapMarkersVolumes, M_mapMarkersIndexes);
    res += (*M_massMatrix) * sol;
}

template<typename Mesh>
void
StructuralOperator<Mesh>::applyBoundaryConditions ( matrix_Type&        matrix,
                                                    vector_Type&        rhs,
                                                    bcHandler_Type&     BCh,
                                                    UInt                offset)
{
    // BC manage for the velocity
    if (offset)
    {
        BCh->setOffset (offset);
    }
    if ( !BCh->bcUpdateDone() )
    {
        BCh->bcUpdate ( *M_dispFESpace->mesh(), M_dispFESpace->feBd(), M_dispFESpace->dof() );
    }

    // vector_Type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rhsFull (rhs, Unique); // bcManages now manages the also repeated parts

    bcManage ( matrix, rhsFull, *M_dispFESpace->mesh(), M_dispFESpace->dof(), *BCh, M_dispFESpace->feBd(), 1., M_data->dataTime()->time() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

}



}
#endif
