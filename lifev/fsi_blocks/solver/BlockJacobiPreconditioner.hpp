//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/interpolation/Interpolation.hpp>

#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>

#ifndef _BlockJacobiPreconditioner_H_
#define _BlockJacobiPreconditioner_H_

namespace LifeV{
namespace Operators
{

class BlockJacobiPreconditioner: public LinearOperatorAlgebra
{
public:
    //! @name Public Types
    //@{

    typedef  LinearOperatorAlgebra                     super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  std::shared_ptr<matrix_Type>            matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  std::shared_ptr<matrixEpetra_Type>      matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  std::shared_ptr<lumpedMatrix_Type>      lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  std::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  std::shared_ptr<mapEpetra_Type>         mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  std::shared_ptr<VectorEpetra_Type>      VectorEpetraPtr_Type;
    typedef  std::shared_ptr<Interpolation>          interpolationPtr_Type;

    typedef FESpace< RegionMesh<LinearTetra> , mapEpetra_Type > FESpace_Type;
    typedef std::shared_ptr<FESpace_Type> 					FESpacePtr_Type;

    typedef std::shared_ptr<BCHandler> BCHandlerPtr_Type;

    //@}

    //! @name Constructors
    //@{
    //! Empty constructor
    BlockJacobiPreconditioner();
    //@}
    virtual ~BlockJacobiPreconditioner();

    //! @name SetUp
    //@{
    //! SetUp
    /*!
     */

    //! @name Set Methods
    //@{

    //! set the communicator
    void setComm ( const commPtr_Type & comm ) { M_comm = comm; };
    //! \warning Transpose of this operator is not supported
    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}
    //! set the domain map
    void setDomainMap(const std::shared_ptr<BlockEpetra_Map> & domainMap){M_operatorDomainMap = domainMap;}
    //! set the range map
    void setRangeMap(const std::shared_ptr<BlockEpetra_Map> & rangeMap){M_operatorRangeMap = rangeMap;}
    //@}


    //! @name
    //@{
    //! \warning No method \c Apply defined for this operator. It return an error code.
    int Apply(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};
    //! Returns the High Order Yosida approximation of the inverse pressure Schur Complement applied to \c X.
    int ApplyInverse(const vector_Type &X, vector_Type &Y) const;
    //! \warning Infinity norm not defined for this operator
    double NormInf() const {return -1.0;}
    //@}

    // @name Attribute access functions
    //@{
    //! Return a character string describing the operator
    const char * Label() const {return M_label.c_str();}
    //! Return the current UseTranspose setting \warning Not Supported Yet.
    bool UseTranspose() const {return M_useTranspose;}
    //! Return false.
    bool HasNormInf() const {return false;}
    //! return a reference to the Epetra_Comm communicator associated with this operator
    const comm_Type & Comm() const {return *M_comm;}
    //! Returns the Epetra_Map object associated with the domain of this operator
    const map_Type & OperatorDomainMap() const {return *(M_monolithicMap->map(Unique));}
    //! Returns the Epetra_Map object associated with the range of this operator
    const map_Type & OperatorRangeMap() const {return *(M_monolithicMap->map(Unique));}
    //@}

    // @name Set the blocks
    //@{

    //! Set the structure block
    void setStructureBlock ( const matrixEpetraPtr_Type & S );

    //! Set the geometry block
    void setGeometryBlock ( const matrixEpetraPtr_Type & G );

    //! Set the fluid blocks
    void setFluidBlocks ( const matrixEpetraPtr_Type & F,
    					  const matrixEpetraPtr_Type & Btranspose,
    					  const matrixEpetraPtr_Type & B);

    //! Set the fluid blocks
    void setFluidBlocks(const matrixEpetraPtr_Type & F,
    					const matrixEpetraPtr_Type & Btranspose,
    		  	  	  	const matrixEpetraPtr_Type & B,
    		  	  	  	const matrixEpetraPtr_Type & D);

    //! Set the coupling blocks
    void setCouplingBlocks ( const matrixEpetraPtr_Type & C1transpose,
			  	  	   	     const matrixEpetraPtr_Type & C2transpose,
			  	  	   	     const matrixEpetraPtr_Type & C2,
			  	  	   	     const matrixEpetraPtr_Type & C1,
		   	   	   	   	  	 const matrixEpetraPtr_Type & C3);

    //! Set the shape derivatives
    void setShapeDerivativesBlocks( const matrixEpetraPtr_Type & ShapeVelocity,
    								const matrixEpetraPtr_Type & ShapePressure);

    //@}

    // @name Update approximations of the block preconditioners
    //@{

    //! Update the approximation of the structure momentum
    void updateApproximatedStructureMomentumOperator();

    //! Update the approximation of the the geometry
    void updateApproximatedGeometryOperator();

    //! Update the approximation of the the geometry
    void updateApproximatedFluidOperator();

    //@}

    //! Show information about the class
    void showMe();

    // @name Set the parameters
    //@{

    //! Interface to set the parameters of each block
    void setOptions(const Teuchos::ParameterList& solversOptions);

    //! Set the monolithic map
    void setMonolithicMap(const mapEpetraPtr_Type& monolithicMap);

    //! Set the use of shape derivatives
    void setUseShapeDerivatives(const bool & useShapeDerivatives);

    //! Set the preconditioner type
    void setFluidPreconditioner(const std::string& type);

    //! Set the blocks needed by the PCD preconditioner
    void setPCDBlocks(const matrixEpetraPtr_Type & Fp, const matrixEpetraPtr_Type & Mp, const matrixEpetraPtr_Type & Mu);

    const char * preconditionerTypeFluid() const { return M_FluidPrec->Label(); }

    void setSubiterateFluidDirichlet(const bool & subiterateFluidDirichlet);

    //! Copy the pointer of the interpolation objects
    void setCouplingOperators_nonconforming( interpolationPtr_Type fluidToStructure, interpolationPtr_Type structureToFluid, mapEpetraPtr_Type lagrangeMap);

    //! Copy the pointer of the fluid velocity fespace
    void setVelocityFESpace( const FESpacePtr_Type& fluid_vel_FESpace) { M_velocityFESpace = fluid_vel_FESpace; };

    //! Copy the pointer of the fluid velocity fespace
    void setBC( const BCHandlerPtr_Type& bc) { M_myBC = bc; };

    //! Copy the value of the timestep
    void setTimeStep( Real dt) { M_timeStep = dt; };

    void setGamma (Real gamma){ M_gamma = gamma;};

    void setBeta (Real beta){ M_beta = beta;};

    void setBDFcoeff (Real coef){ M_bdfCoef = coef; M_useBDFStructure = true; };

    //@}

private:

    // @name Set the parameters
    //@{

    //! Set the list of the structure momentum
    void setStructureMomentumOptions(const parameterListPtr_Type & _oList);

    //! Set the list of the geometry
    void setGeometryOptions(const parameterListPtr_Type & _oList);

    //! Set the list of the fluid momentum
    void setFluidMomentumOptions(const parameterListPtr_Type & _oList);

    //! Set the list of the shur complement of the fluid
    void setSchurOptions(const parameterListPtr_Type & _oList);

    //! Set the list of the shur complement of the fluid
    void setPressureMassOptions(const parameterListPtr_Type & _oList);

    //@}

    //! Create the domain and the range maps
    void setMaps();

    std::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    std::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    // @name Approximation of the blocks
    //@{
    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedStructureMomentumOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedGeometryOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedFluidMomentumOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedSchurComplementOperator;

    std::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedSchurComplementCouplingOperator;
    //@}

    // @name Parameter list of each block
    //@{

    //! Parameters for the structure
    parameterListPtr_Type M_structureMomentumOptions;

    //! Parameters for the geometry
    parameterListPtr_Type M_geometryOptions;

    //! Parameters for the fluid momentum
    parameterListPtr_Type M_fluidMomentumOptions;

    //! Parameters for the shur complent of the fluid
    parameterListPtr_Type M_schurOptions;

    //! Parameters for the shur complent of the couplig
    parameterListPtr_Type M_schurCouplingOptions;

    //! Parameters for the pressure mass of the PCD
    parameterListPtr_Type M_pressureMassOptions;
    //@}

    // @name Parameter list of each block
    //@{

    //! Structure block
    matrixEpetraPtr_Type M_S;

    //! Geometry block
    matrixEpetraPtr_Type M_G;

    //! Fluid blocks
    matrixEpetraPtr_Type M_F;
    matrixEpetraPtr_Type M_Btranspose;
    matrixEpetraPtr_Type M_B;
    matrixEpetraPtr_Type M_D;

    std::shared_ptr<Operators::NavierStokesPreconditionerOperator> M_FluidPrec;
    // Epetra Operator needed to solve the linear system
    std::shared_ptr<Operators::InvertibleOperator> M_invOper;
    std::shared_ptr<Operators::NavierStokesOperator> M_oper;

    //! Coupling blocks
    matrixEpetraPtr_Type M_C1transpose;
    matrixEpetraPtr_Type M_C2transpose;
    matrixEpetraPtr_Type M_C2;
    matrixEpetraPtr_Type M_C1;
    matrixEpetraPtr_Type M_C3;
    matrixEpetraPtr_Type M_shapeVelocity;
    matrixEpetraPtr_Type M_shapePressure;
    //@}

    mapEpetraPtr_Type M_monolithicMap;

    std::shared_ptr<Epetra_Vector> M_invD;

    matrixEpetraPtr_Type M_schurComplement;
    matrixEpetraPtr_Type M_schurComplementCoupling;

    //! Label
    const std::string M_label;

    bool M_shapeDerivatives;

    // Offsets
    Real M_fluidVelocity;
    Real M_fluid;
    Real M_structure;
    Real M_lambda;

    //! Vectors needed for the applyInverse - input vectors associated to each part
    std::shared_ptr<VectorEpetra_Type > M_X_velocity;
    std::shared_ptr<VectorEpetra_Type > M_X_pressure;
    std::shared_ptr<VectorEpetra_Type > M_X_displacement;
    std::shared_ptr<VectorEpetra_Type > M_X_lambda;
    std::shared_ptr<VectorEpetra_Type > M_X_geometry;

    //! Vectors needed for the applyInverse - output vectors associated to each part
    std::shared_ptr<VectorEpetra_Type > M_Y_velocity;
    std::shared_ptr<VectorEpetra_Type > M_Y_pressure;
    std::shared_ptr<VectorEpetra_Type > M_Y_displacement;
    std::shared_ptr<VectorEpetra_Type > M_Y_lambda;
    std::shared_ptr<VectorEpetra_Type > M_Y_geometry;

    //! PCD blocks
    matrixEpetraPtr_Type M_Fp;
    matrixEpetraPtr_Type M_Mp;
    matrixEpetraPtr_Type M_Mu;

    bool M_subiterateFluidDirichlet;

    bool M_useStabilization;

    interpolationPtr_Type M_FluidToStructureInterpolant;
    interpolationPtr_Type M_StructureToFluidInterpolant;
    std::shared_ptr<MapEpetra> M_lagrangeMap;

    bool M_nonconforming;

    FESpacePtr_Type M_velocityFESpace;
    std::shared_ptr<BCHandler> M_myBC;

    Real M_timeStep;

    Real M_gamma;

    Real M_beta;

    Real M_bdfCoef;
    bool M_useBDFStructure;
};

} /* end namespace Operators */
} //end namespace
#endif
