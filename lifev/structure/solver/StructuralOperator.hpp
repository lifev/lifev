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


#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>


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

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/NonLinearRichardson.hpp>

#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>


namespace LifeV
{

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
template < typename Mesh,
         typename SolverType = LifeV::SolverAztecOO >

class StructuralOperator
{
public:

    //!@name Type definitions
    //@{
    typedef Real ( *function ) ( const Real&, const Real&, const Real&, const Real&, const ID& );
    typedef boost::function<Real ( const Real&, const Real&, const Real&, const Real&, const ID& ) > source_Type;

    typedef StructuralConstitutiveLaw<Mesh>               material_Type;
    typedef boost::shared_ptr<material_Type>              materialPtr_Type;

    typedef BCHandler                                     bcHandlerRaw_Type;
    typedef boost::shared_ptr<bcHandlerRaw_Type>          bcHandler_Type;

    typedef SolverType                                    solver_Type;

    typedef typename solver_Type::matrix_type             matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                matrixPtr_Type;
    typedef typename solver_Type::vector_type             vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;

    typedef typename SolverType::prec_raw_type            precRaw_Type;
    typedef typename SolverType::prec_type                prec_Type;

    typedef StructuralConstitutiveLawData                 data_Type;

    typedef RegionMesh<LinearTetra >                      mesh_Type;
    typedef std::vector< mesh_Type::element_Type const*> vectorVolumes_Type;

    typedef std::map< UInt, vectorVolumes_Type>           mapMarkerVolumes_Type;
    typedef boost::shared_ptr<mapMarkerVolumes_Type>      mapMarkerVolumesPtr_Type;

    typedef typename mesh_Type::element_Type                        meshEntity_Type;

    typedef typename boost::function2<bool, const UInt, const UInt> comparisonPolicy_Type;

    typedef MarkerSelector<meshEntity_Type, comparisonPolicy_Type> markerSelector_Type;
    typedef boost::scoped_ptr<markerSelector_Type>          markerSelectorPtr_Type;
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
    void setup ( boost::shared_ptr<data_Type> data,
                 const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   FESpace,
                 bcHandler_Type&       BCh,
                 boost::shared_ptr<Epetra_Comm>&     comm
               );

    /*!
      \param data_file GetPot data file
      \param refFE reference FE for the displacement
      \param comm the Epetra Comunicator
    */
    void setup ( boost::shared_ptr<data_Type> data,
                 const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   FESpace,
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
                 const boost::shared_ptr< FESpace<Mesh, MapEpetra> >&   dFESpace,
                 boost::shared_ptr<Epetra_Comm>&     comm,
                 const boost::shared_ptr<const MapEpetra>&       monolithicMap,
                 UInt       offset = 0
                                     //boost::shared_ptr<FESpace<Mesh, MapEpetra> >   uFESpace=0
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
    void iterate ( bcHandler_Type& bch );

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
    /*!
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

    void evalConstraintTensor();

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 space function describing the initial velocity
      \param a0 space function describing the initial acceleration
    */
    void initialize ( const function& d0, const function& w0, const function& a0 );

    //! Sets the initial velocity
    /*!
      \param w0 space function describing the initial velocity
    */
    void initializeVel ( const vector_Type& w0);

    //! Sets the initial displacement, velocity, acceleration
    /*!
      \param d0 space function describing the initial displacement
      \param w0 empty vector
      \param a0 empty vector
    */
    void initialize ( vectorPtr_Type d0,  vectorPtr_Type w0 = vectorPtr_Type(),  vectorPtr_Type a0 = vectorPtr_Type() );

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

    //void updateMatrix(matrix_Type & bigMatrixStokes);// used for monolithic
    //void updateCoupling(matrix_Type couplingMatrix);// used for monolithic

    //@}

    //! @name Set Methods
    //@{

    //!Setters
    //! Set the BCHandler object
    void setBC (bcHandler_Type& BCd)
    {
        M_BCh = BCd;
    }

    //! Set the source object
    void setSourceTerm ( source_Type const& s )
    {
        M_source = s;
    }

    //! Set the preconditioner
    void resetPrec (bool reset = true)
    {
        if (reset)
        {
            M_linearSolver.precReset();
        }
    }

    // //! Set the displacement
    // virtual void setDisp(const vector_Type& disp) {*M_disp = disp;} // used for monolithic

    //! Set the recur parameter
    void setRecur (UInt recur)
    {
        M_recur = recur;
    }

    //! Set the data fields with the Getpot data file
    void setDataFromGetPot ( const GetPot& dataFile );

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
    FESpace<Mesh, MapEpetra>& dFESpace()
    {
        return M_FESpace;
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
    //! Get the velocity
    //vector_Type& velocity()         { return *M_vel; }

    //! Get the velocity
    //vector_Type& acceleration()         { return *M_acc; }

    //! Get the right hand sde without BC
    vectorPtr_Type& rhsWithoutBC()
    {
        return M_rhsNoBC;
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
    boost::shared_ptr<data_Type> data()       const
    {
        return M_data;
    }

    void apply ( const vector_Type& sol, vector_Type& res) const;


    //! Get the density
    mapMarkerVolumesPtr_Type mapMarkersVolumes() const
    {
        return M_mapMarkersVolumes;
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
        return M_FESpace->dim();
    }


    //! construct the map between the markers and the volumes
    /*!
      \param VOID
      \return VOID
    */
    void setupMapMarkersVolumes ( void );

    //!Protected Members

    boost::shared_ptr<data_Type>         M_data;

    boost::shared_ptr<FESpace<Mesh, MapEpetra> >      M_FESpace;

    boost::shared_ptr<const Displayer>         M_Displayer;

    Int                                  M_me;

    //! data for solving tangent problem with aztec
    boost::shared_ptr<solver_Type>       M_linearSolver;

    //! Elementary matrices and vectors
    boost::shared_ptr<MatrixElemental>   M_elmatM;

    //! linearized velocity
    vectorPtr_Type                       M_disp;
    //vectorPtr_Type                       M_vel;
    //vectorPtr_Type                       M_acc;

    //! right  hand  side displacement
    vectorPtr_Type                       M_rhs;

    //! right  hand  side velocity
    //  vectorPtr_Type                       M_rhsW;

    //! right  hand  side velocity
    //vectorPtr_Type                       M_rhsA;

    //! right  hand  side
    vectorPtr_Type                       M_rhsNoBC;

    //! right  hand  side
    //boost::shared_ptr<vector_Type>       M_f;

    //! residual
    boost::shared_ptr<vector_Type>       M_residual_d;

    //! Components of the Constraint Tensor
    vectorPtr_Type                       M_sxx;
    vectorPtr_Type                       M_syy;
    vectorPtr_Type                       M_szz;

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

};

//====================================
// Constructor
//=====================================

template <typename Mesh, typename SolverType>
StructuralOperator<Mesh, SolverType>::StructuralOperator( ) :
    M_data                       ( ),
    M_FESpace                    ( ),
    M_Displayer                  ( ),
    M_me                         ( 0 ),
    M_linearSolver               ( ),
    M_elmatM                     ( ),
    M_disp                       ( ),
    M_rhsNoBC                    ( ),
    M_residual_d                 ( ),
    M_sxx                        (/*M_localMap*/),//useless
    M_syy                        (/*M_localMap*/),//useless
    M_szz                        (/*M_localMap*/),//useless
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
    M_mapMarkersVolumes          ( )
{
    //    M_Displayer->leaderPrint("I am in the constructor for the solver");
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::setup (boost::shared_ptr<data_Type>          data,
                                             const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                             bcHandler_Type&                    BCh,
                                             boost::shared_ptr<Epetra_Comm>&   comm)
{
    setup (data, dFESpace, comm);
    M_BCh = BCh;
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::setup (boost::shared_ptr<data_Type>        data,
                                             const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                             boost::shared_ptr<Epetra_Comm>&     comm)
{
    setup ( data, dFESpace, comm, dFESpace->mapPtr(), (UInt) 0 );

    M_rhs.reset                        ( new vector_Type (*M_localMap) );
    M_rhsNoBC.reset                    ( new vector_Type (*M_localMap) );
    M_sxx.reset                        ( new vector_Type (*M_localMap) );
    M_syy.reset                        ( new vector_Type (*M_localMap) );
    M_szz.reset                        ( new vector_Type (*M_localMap) );
    M_linearSolver.reset               ( new SolverType ( comm ) );
    M_disp.reset                       ( new vector_Type (*M_localMap) );
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::setup (boost::shared_ptr<data_Type>        data,
                                             const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                             boost::shared_ptr<Epetra_Comm>&     comm,
                                             const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                             UInt                                offset)
{
    M_data                            = data;
    M_FESpace                         = dFESpace;
    M_Displayer.reset                 (new Displayer (comm) );
    M_me                              = comm->MyPID();
    M_elmatM.reset                    ( new MatrixElemental ( M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    M_localMap                        = monolithicMap;
    M_massMatrix.reset                (new matrix_Type (*M_localMap) );
    M_systemMatrix.reset              (new matrix_Type (*M_localMap) );
    M_jacobian.reset                  (new matrix_Type (*M_localMap) );

    M_offset                          = offset;

    M_material.reset ( material_Type::StructureMaterialFactory::instance().createObject ( M_data->solidType() ) );
    M_material->setup ( dFESpace, M_localMap, M_offset, M_data, M_Displayer );

    if ( M_data->verbose() )
    {
        M_out_iter.open ( "out_iter_solid" );
        M_out_res.open ( "out_res_solid" );
    }
    M_mapMarkersVolumes.reset ( new mapMarkerVolumes_Type() );
    this->setupMapMarkersVolumes();
}


template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::setupMapMarkersVolumes ( void )
{

    this->M_Displayer->leaderPrint (" S-  Building the map between volumesMarkers <--> volumes \n");

    //We first loop over the vector of the material_flags
    for (  UInt i (0); i < M_data->vectorFlags().size(); i++ )
    {

        //Create the functor to extract volumes
        markerSelectorPtr_Type ref ( new markerSelector_Type (M_data->vectorFlags() [i]) );

        //Number of volumes with the current marker
        UInt numExtractedVolumes = this->M_FESpace->mesh()->elementList().countAccordingToPredicate ( *ref );

        this->M_Displayer->leaderPrint (" Current marker: ", M_data->vectorFlags() [i]);
        this->M_Displayer->leaderPrint (" \n");
        this->M_Displayer->leaderPrint (" Number of volumes:", numExtractedVolumes);
        this->M_Displayer->leaderPrint (" \n");

        //Vector large enough to contain the number of volumes with the current marker
        vectorVolumes_Type extractedVolumes ( numExtractedVolumes );

        //Extracting the volumes
        extractedVolumes = this->M_FESpace->mesh()->elementList().extractAccordingToPredicate ( *ref );

        //Insert the correspondande Marker <--> List of Volumes inside the map
        M_mapMarkersVolumes->insert ( pair<UInt, vectorVolumes_Type> (M_data->vectorFlags() [i], extractedVolumes) ) ;

        //Cleaning the vector
        extractedVolumes.clear();

    }
}




template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateSystem ( void )
{
    updateSystem (M_systemMatrix);
}

template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateSystem ( matrixPtr_Type& mat_stiff)
{
    M_Displayer->leaderPrint (" S-  Updating mass term on right hand side... ");

    LifeChrono chrono;
    chrono.start();

    //Compute the new Stiffness Matrix
    M_material->computeStiffness (*M_disp, M_rescaleFactor, M_data, M_mapMarkersVolumes, M_Displayer);

    if ( M_data->solidType() == "linearVenantKirchhoff" || M_data->solidType() == "nonLinearVenantKirchhoff" )
    {
        *mat_stiff += *M_material->stiffMatrix();
        mat_stiff->globalAssemble();
    }


    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );

}

template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateSourceTerm ( source_Type const& source )
{
    vector_Type rhs ( *M_localMap );

    VectorElemental M_elvec (M_FESpace->fe().nbFEDof(), nDimensions);
    UInt nc = nDimensions;

    // loop on volumes: assembling source term
    for ( UInt i = 1; i <= M_FESpace->mesh()->numVolumes(); ++i )
    {

        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( i ) );

        M_elvec.zero();

        for ( UInt ic = 0; ic < nc; ++ic )
        {
            compute_vec ( source, M_elvec, M_FESpace->fe(),  M_data->dataTime()->time(), ic ); // compute local vector
            assembleVector ( *rhs, M_elvec, M_FESpace->fe(), M_FESpace->dof(), ic, ic * M_FESpace->getDim() ); // assemble local vector into global one
        }
    }
    M_rhsNoBC += rhs;
}

template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::buildSystem ( const Real coefficient )
{
    M_Displayer->leaderPrint ("  S-  Computing constant matrices ...          ");
    LifeChrono chrono;
    chrono.start();

    computeMassMatrix ( coefficient );
    M_material->computeLinearStiff (M_data, M_mapMarkersVolumes);

    chrono.stop();
    M_Displayer->leaderPrintMax ( "done in ", chrono.diff() );

}


// template <typename Mesh, typename SolverType>
// void  StructuralOperator<Mesh, SolverType>::buildSystem(matrix_Type & bigMatrixStokes, const Real& timeAdvanceCoefficient, const Real& factor)
// {}
// ;


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::computeMassMatrix ( const Real factor)
{
    UInt totalDof = M_FESpace->dof().numTotalDof();

    //! Number of displacement components
    UInt nc = nDimensions;
    const Real factorMassMatrix = factor * M_data->rho();

    //! Elementary computation and matrix assembling
    //! Loop on elements
    for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); i++ )
    {

        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( i ) );

        M_elmatM->zero();

        // mass
        // The method mass is implemented in AssemblyElemental.cpp
        mass ( factorMassMatrix , *M_elmatM, M_FESpace->fe(), 0, 0, nDimensions );

        //! assembling
        for ( UInt ic = 0; ic < nc; ic++ )
        {
            //mass
            assembleMatrix ( *M_massMatrix, *M_elmatM, M_FESpace->fe(), M_FESpace->dof(), ic, ic, M_offset +  ic * totalDof, M_offset +  ic * totalDof);
        }
    }

    //getComunicator()->Barrier();

    M_massMatrix->globalAssemble();

    //*massStiff *= factor; //M_data.dataTime()->timeStep() * M_rescaleFactor;
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::iterate ( bcHandler_Type& bch )
{
    LifeChrono chrono;

    // matrix and vector assembling communication
    M_Displayer->leaderPrint ("  S-  Solving the system ... \n");

    M_BCh = bch;

    Real abstol  = 1.e-7;
    Real reltol  = 1.e-7;
    UInt maxiter = 50;
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

    // updateVelAndAcceleration();

    // std::cout << "iterate: d norm       = " << M_disp->norm2() << std::endl;

    //These two lines mut be checked fo FSI. With the linear solver, they have a totally
    //different expression. For structural problems it is not used.
    evalResidualDisplacement (*M_disp); //\todo related to FSI. Should be caled from outside
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::iterateLin ( bcHandler_Type& bch )
{
    vector_Type rhsFull (M_rhsNoBC->map() );
    Real zero (0.);
    if ( !bch->bcUpdateDone() )
    {
        bch->bcUpdate ( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );
    }
    bcManageVector ( rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *bch, M_FESpace->feBd(), M_data->dataTime()->time(), 1.0 );
    solveJacobian (*M_disp, rhsFull, zero, bch);
    evalResidualDisplacementLin (*M_disp);
}


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::showMe ( std::ostream& c  ) const
{
    c << "\n*** StructuralOperator::showMe method" << std::endl;

    M_data->showMe ( c );

}

template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::computeMatrix ( matrixPtr_Type& stiff, const vector_Type& sol,  Real const& /*factor*/)
{
    M_Displayer->leaderPrint ( " Computing residual ... \t\t\t");

    LifeChrono chrono;
    chrono.start();

    //! It is right to do globalAssemble() inside the M_material class
    M_material->computeStiffness ( sol, 1., M_data, M_mapMarkersVolumes, M_Displayer);

    if ( M_data->solidType() == "linearVenantKirchhoff" || M_data->solidType() == "nonLinearVenantKirchhoff" )
    {
        *stiff = *M_material->stiffMatrix();
        *stiff += *M_massMatrix;
        stiff->globalAssemble();
    }

    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::evalResidual ( vector_Type& residual, const vector_Type& solution, Int iter)
{

    //This method call the M_material computeStiffness
    computeMatrix (M_systemMatrix, solution, 1.);


    M_Displayer->leaderPrint ("    S- Updating the boundary conditions ... \t");
    LifeChrono chrono;

    if ( !M_BCh->bcUpdateDone() )
    {
        M_BCh->bcUpdate ( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );
    }

    // ignoring non-local entries, Otherwise they are summed up lately
    if ( M_data->solidType() == "linearVenantKirchhoff" || M_data->solidType() == "nonLinearVenantKirchhoff" )
    {
        chrono.start();

        matrix_Type matrixFull (*M_systemMatrix);

        if (iter == 0)
        {
            *M_rhs = *M_rhsNoBC;
            bcManageVector ( *M_rhs, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(),  M_data->dataTime()->time(), 1.0 );
        }

        bcManageMatrix ( matrixFull, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), 1.0 );

        residual  = matrixFull * solution;
        residual -= *M_rhs;
        chrono.stop();
        M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
    }
    else //NH and Exp
    {
        chrono.start();
        *M_rhs = *M_rhsNoBC;
        residual = *M_massMatrix * solution;
        residual += *M_material->stiffVector();
        vector_Type solRep (solution, Repeated);
        bcManageResidual ( residual, *M_rhs, solRep, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), M_data->dataTime()->time(), 1.0 );
        residual -= *M_rhs;
        chrono.stop();
        M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
    }
}

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::evalResidualDisplacement ( const vector_Type& solution )
{

    M_Displayer->leaderPrint ("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    if ( M_data->solidType() == "linearVenantKirchhoff" || M_data->solidType() == "nonLinearVenantKirchhoff" )
    {
        M_residual_d.reset (new vector_Type ( *M_systemMatrix * solution ) );
        *M_residual_d -= *M_rhsNoBC;
    }
    else
    {
        M_residual_d.reset (new vector_Type ( *M_material->stiffVector() ) );
        *M_residual_d -= *M_rhsNoBC;
    }
    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::evalResidualDisplacementLin ( const vector_Type& solution )
{

    M_Displayer->leaderPrint ("    S- Computing the residual displacement for the structure..... \t");
    LifeChrono chrono;
    chrono.start();

    //This definition of residual_d is similar to the one of iterateLin in VenantKirchhoffSolver
    M_residual_d.reset (new vector_Type ( (*M_jacobian) *solution) );

    chrono.stop();
    M_Displayer->leaderPrintMax ("done in ", chrono.diff() );
}



template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::evalConstraintTensor()
{
    vector_Type count (*M_localMap);

    *M_sxx *= 0.;
    *M_syy *= 0.;
    *M_szz *= 0.;

    for ( UInt ielem = 0; ielem < M_FESpace->mesh()->numVolumes(); ielem++ )
    {
        //UInt elem = M_FESpace->mesh()->volumeList( ielem ).id();
        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( ielem ) );

        //int    marker = M_FESpace->mesh()->volumeList( ielem ).marker();
        Real s      = 0;
        Real volume = M_FESpace->fe().detJac (0);

        for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
            {
                Int i    = M_FESpace->fe().patternFirst (k);
                Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

                s += (2 * M_data->mu (1) + M_data->lambda (1) ) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 0 , ig ) *
                     (*M_disp) [idof + 0 * M_FESpace->dim()];

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 1 , ig ) *
                     (*M_disp) [idof + 1 * M_FESpace->dim()];

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 2 , ig ) *
                     (*M_disp) [idof + 2 * M_FESpace->dim()];

                count[idof]++;
            }
        }

        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
        {
            Int i    = M_FESpace->fe().patternFirst (k);
            Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

            (*M_sxx) [idof] += s / M_FESpace->fe().detJac (0);
        }

        s = 0;

        for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
            {
                Int i    = M_FESpace->fe().patternFirst (k);
                Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 0 , ig ) *
                     (*M_disp) [idof + 0 * M_FESpace->dim()];

                s += (2 * M_data->mu (1) + M_data->lambda (1) ) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 1 , ig ) *
                     (*M_disp) [idof + 1 * M_FESpace->dim()];

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 2 , ig ) *
                     (*M_disp) [idof + 2 * M_FESpace->dim()];
                //                         M_sxx[idof] += s;

                //                M_syy[idof] += s/volume;
            }
        }

        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
        {
            Int i    = M_FESpace->fe().patternFirst (k);
            Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

            (*M_syy) [idof] += s / volume;
        }


        s = 0;

        for ( Int ig = 0; ig < M_FESpace->fe().nbQuadPt; ++ig )
        {
            for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
            {
                Int i    = M_FESpace->fe().patternFirst (k);
                Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 0 , ig ) *
                     (*M_disp) [idof + 0 * M_FESpace->dim()];

                s += M_data->lambda (1) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 1 , ig ) *
                     (*M_disp) [idof + 1 * M_FESpace->dim()];

                s += (2 * M_data->mu (1) + M_data->lambda (1) ) *
                     M_FESpace->fe().weightDet ( ig ) *
                     M_FESpace->fe().phiDer ( k, 2 , ig ) *
                     (*M_disp) [idof + 2 * M_FESpace->dim()];

                //                         M_sxx[idof] += s;
            }
        }

        for ( Int k = 0; k < M_FESpace->fe().nbFEDof(); ++k )
        {
            Int i    = M_FESpace->fe().patternFirst (k);
            Int idof = M_FESpace->dof().localToGlobal (M_FESpace->fe().currentLocalId(), i + 1);

            (*M_szz) [idof] += s / M_FESpace->fe().detJac (0);
        }

    }

    for (UInt ii = 1; ii <= M_FESpace->dim(); ++ii)
    {
        (*M_sxx) [ii] /= count[ii];
        (*M_syy) [ii] /= count[ii];
        (*M_szz) [ii] /= count[ii];
    }
}


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::initialize ( vectorPtr_Type disp, vectorPtr_Type /*vel*/, vectorPtr_Type /*acc*/)
{
    *M_disp = *disp;
    //  if (vel.get())
    //    initializeVel(*vel);
}
/*
  template <typename Mesh, typename SolverType>
  void
  StructuralOperator<Mesh, SolverType>::initializeVel( const vector_Type& vel)
  {
  *M_vel = vel;
  }
*/

template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::initialize ( const function& d0, const function& w0, const function& a0 )
{
    M_FESpace->interpolate ( static_cast<typename FESpace<Mesh, MapEpetra>::function_Type> ( d0 ), *M_disp, 0.0);
    //M_FESpace->interpolate(w0, *M_vel , 0.0);
    // M_FESpace->interpolate(a0, *M_acc , 0.0);
}

/*
//Matteo this method isn't necessary timeAdvance compute the accelerate and velocity
template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateVelAndAcceleration()
{
Real DeltaT = M_data->dataTime()->timeStep();

*M_acc = (2.0 /( M_zeta * pow(DeltaT,2) ))  * (*M_disp)  - *M_rhsA;
*M_vel = *M_rhsW + M_theta * DeltaT * (*M_acc) ;
}
*/

template<typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::reduceSolution ( Vector& displacement, Vector& velocity )
{
    vector_Type disp (*M_disp, 0);
    //vector_Type vel(*M_vel , 0);

    if ( comunicator()->MyPID() == 0 )
    {
        for ( UInt iDof = 0; iDof < nDimensions * dim(); ++iDof )
        {
            disp[ iDof ] = displacement[ iDof + 1 ];
            //vel [ iDof ] = velocity    [ iDof + 1 ];
        }
    }
}


template <typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::setDataFromGetPot ( const GetPot& dataFile )
{
    M_linearSolver->setDataFromGetPot ( dataFile, "solid/solver" );
    M_linearSolver->setupPreconditioner (dataFile, "solid/prec");
    M_rescaleFactor = dataFile ( "solid/rescaleFactor", 0. );
    // This is done in DataClass
    // UInt marker = M_FESpace->mesh()->volumeList( 1 ).marker();
    // if (!M_data->young(marker))
    //     M_data->setYoung(dataFile( "solid/physics/young", 0. ), marker);
    // if (!M_data->poisson(marker))
    //     M_data->setPoisson(dataFile( "solid/physics/poisson", 0. ), marker);
}


//Method UpdateJacobian
template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateJacobian ( const vector_Type& sol, matrixPtr_Type& jacobian  )
{
    M_Displayer->leaderPrint ("  S-  Solid: Updating JACOBIAN... ");

    LifeChrono chrono;
    chrono.start();

    M_material->updateJacobianMatrix (sol, M_data, M_mapMarkersVolumes, M_Displayer);

    jacobian.reset (new matrix_Type (*M_localMap) );
    *jacobian += *M_material->jacobian();
    *jacobian += *M_massMatrix;

    chrono.stop();
    M_Displayer->leaderPrintMax ("   ... done in ", chrono.diff() );

}

//Method UpdateJacobian
template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::updateJacobian ( const vector_Type& sol)
{
    updateJacobian (sol, M_jacobian);
}

//solveJac( const Vector& res, Real& linear_rel_tol, Vector &step)
template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::
solveJac ( vector_Type& step, const vector_Type& res, Real& linear_rel_tol)
{
    updateJacobian ( *M_disp, M_jacobian );
    solveJacobian (step,  res, linear_rel_tol, M_BCh);
}


//Method SolveJacobian
template <typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::
solveJacobian (vector_Type&           step,
               const vector_Type&     res,
               Real&                /*linear_rel_tol*/,
               bcHandler_Type&       /* BCh*/)
{
    LifeChrono chrono;

    M_jacobian->globalAssemble();
    matrixPtr_Type matrFull (new matrix_Type (*M_localMap) );
    *matrFull += *M_jacobian;

    M_Displayer->leaderPrint ("\tS'-  Solving the linear system ... \n");

    M_Displayer->leaderPrint ("\tS'-  Applying boundary conditions      ... ");


    if ( !M_BCh->bcUpdateDone() )
    {
        M_BCh->bcUpdate ( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );
    }
    bcManageMatrix ( *matrFull, *M_FESpace->mesh(), M_FESpace->dof(), *M_BCh, M_FESpace->feBd(), 1.0 );

    M_Displayer->leaderPrintMax ( "done in ", chrono.diff() );

    M_Displayer->leaderPrint ("\tS'-  Solving system                    ... \n");
    chrono.start();

    M_linearSolver->setMatrix (*matrFull);

    M_linearSolver->solveSystem ( res, step, matrFull );

    //     matrFull->spy("J");
    //     M_material->stiffMatrix()->spy("S");
    //     M_systemMatrix->spy("M");
    chrono.stop();
}

template<typename Mesh, typename SolverType>
void StructuralOperator<Mesh, SolverType>::apply ( const vector_Type& sol, vector_Type& res) const
{
    M_material->apply (sol, res, M_mapMarkersVolumes);
    res += (*M_massMatrix) * sol;
}

template<typename Mesh, typename SolverType>
void
StructuralOperator<Mesh, SolverType>::applyBoundaryConditions ( matrix_Type&        matrix,
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
        BCh->bcUpdate ( *M_FESpace->mesh(), M_FESpace->feBd(), M_FESpace->dof() );
    }

    // vector_Type rhsFull(rhs, Repeated, Zero); // ignoring non-local entries, Otherwise they are summed up lately
    vector_Type rhsFull (rhs, Unique); // bcManages now manages the also repeated parts

    bcManage ( matrix, rhsFull, *M_FESpace->mesh(), M_FESpace->dof(), *BCh, M_FESpace->feBd(), 1., M_data->dataTime()->time() );

    // matrix should be GlobalAssembled by  bcManage

    rhs = rhsFull;

}



}
#endif
