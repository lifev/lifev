//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/LinearOperator.hpp>

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
			  	  	   	     const matrixEpetraPtr_Type & C1);

    //@}

    // @name Update approximations of the block preconditioners
    //@{

    //! Update the approximation of the structure momentum
    void updateApproximatedStructureMomentumOperator();

    //! Update the approximation of the the geometry
    void updateApproximatedGeometryOperator();

    //! Update the approximation of the fluid momentum
    void updateApproximatedFluidMomentumOperator();

    //! Update the shur complement associated to the fluid
    void updateApproximatedSchurComplementOperator();

    //! Update the shur complement associated to the coupling
    void updateApproximatedSchurComplementCouplingOperator();
    //@}

    //! Show information about the class
    void showMe();

    // @name Set the parameters
    //@{

    //! Interface to set the parameters of each block
    void setOptions(const Teuchos::ParameterList& solversOptions);
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

    //! Set the list of the shur complement of the coupling
    void setSchurCouplingOptions(const parameterListPtr_Type & _oList);
    //@}

    //! Create the domain and the range maps
    void setMaps();

    //! create the matrix B*diag(F)^-1*Btranspose
    void buildShurComplement();

    boost::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    boost::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    // @name Approximation of the blocks
    //@{
    Operators::ApproximatedInvertibleRowMatrix * M_approximatedStructureMomentumOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedGeometryOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedFluidMomentumOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedSchurComplementOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedSchurComplementCouplingOperator;
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

    //! Coupling blocks
    matrixEpetraPtr_Type M_C1transpose;
    matrixEpetraPtr_Type M_C2transpose;
    matrixEpetraPtr_Type M_C2;
    matrixEpetraPtr_Type M_C1;
    //@}

    mapEpetraPtr_Type M_monolithicMap;

    //! Label
    const std::string M_label;
};

} /* end namespace Operators */
} //end namespace
#endif
