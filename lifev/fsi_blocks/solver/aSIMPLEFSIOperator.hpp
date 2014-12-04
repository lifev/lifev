//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/LinearOperator.hpp>
#include <lifev/navier_stokes/solver/aSIMPLEOperator.hpp>
#include <lifev/navier_stokes/solver/NavierStokesOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#ifndef _aSIMPLEFSIOPERATOR_H_
#define _aSIMPLEFSIOPERATOR_H_

namespace LifeV{
namespace Operators
{

class aSIMPLEFSIOperator: public LinearOperator
{
public:
    //! @name Public Types
    //@{

    typedef  LinearOperator                            super;
    typedef  Epetra_CrsMatrix                          matrix_Type;
    typedef  boost::shared_ptr<matrix_Type>            matrixPtr_Type;
    typedef  MatrixEpetra<Real>                        matrixEpetra_Type;
    typedef  boost::shared_ptr<matrixEpetra_Type>      matrixEpetraPtr_Type;
    typedef  Epetra_Vector                             lumpedMatrix_Type;
    typedef  boost::shared_ptr<lumpedMatrix_Type>      lumpedMatrixPtr_Type;
    typedef  super::comm_Type                          comm_Type;
    typedef  super::commPtr_Type                       commPtr_Type;
    typedef  boost::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                 mapEpetra_Type;
    typedef  boost::shared_ptr<mapEpetra_Type>         mapEpetraPtr_Type;
    typedef  VectorEpetra                              VectorEpetra_Type;
    typedef  boost::shared_ptr<VectorEpetra_Type>      VectorEpetraPtr_Type;

    //@}

    //! @name Constructors
    //@{
    //! Empty constructor
    aSIMPLEFSIOperator();
    //@}
    virtual ~aSIMPLEFSIOperator();

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
    void setDomainMap(const boost::shared_ptr<BlockEpetra_Map> & domainMap){M_operatorDomainMap = domainMap;}
    //! set the range map
    void setRangeMap(const boost::shared_ptr<BlockEpetra_Map> & rangeMap){M_operatorRangeMap = rangeMap;}
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
    const map_Type & OperatorDomainMap() const {return *(M_operatorDomainMap->monolithicMap());}
    //! Returns the Epetra_Map object associated with the range of this operator
    const map_Type & OperatorRangeMap() const {return *(M_operatorRangeMap->monolithicMap());}
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

    //@}

    //! Create the domain and the range maps
    void setMaps();

    boost::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    boost::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    // @name Approximation of the blocks
    //@{
    boost::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedStructureMomentumOperator;

    boost::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedGeometryOperator;

    boost::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedFluidMomentumOperator;

    boost::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedSchurComplementOperator;

    boost::shared_ptr<Operators::ApproximatedInvertibleRowMatrix> M_approximatedSchurComplementCouplingOperator;
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

    boost::shared_ptr<Operators::aSIMPLEOperator> M_FluidDirichlet;
    // Epetra Operator needed to solve the linear system
    boost::shared_ptr<Operators::InvertibleOperator> M_invOper;
    boost::shared_ptr<Operators::NavierStokesOperator> M_oper;

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

    boost::shared_ptr<Epetra_Vector> M_invD;

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
    boost::shared_ptr<VectorEpetra_Type > M_X_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_X_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_X_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_X_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_X_geometry;

    //! Vectors needed for the applyInverse - output vectors associated to each part
    boost::shared_ptr<VectorEpetra_Type > M_Y_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_Y_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_Y_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_Y_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_Y_geometry;


};

} /* end namespace Operators */
} //end namespace
#endif
