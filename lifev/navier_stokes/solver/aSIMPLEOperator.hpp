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

#ifndef _aSIMPLEOPERATOR_H_
#define _aSIMPLEOPERATOR_H_

namespace LifeV{
namespace Operators
{

class aSIMPLEOperator: public LinearOperator
{
public:
    //! @name Public Types
    //@{
    typedef boost::numeric::ublas::vector< boost::shared_ptr<vector_Type> > Zdata;
    typedef boost::numeric::ublas::matrix< boost::shared_ptr<vector_Type> > ZZdata;
    
    typedef LinearOperator                            super;
    typedef  Epetra_CrsMatrix                         matrix_Type;
    typedef  boost::shared_ptr<matrix_Type>           matrixPtr_Type;
    typedef  MatrixEpetra<Real>                       matrixEpetra_Type;
    typedef  boost::shared_ptr<matrixEpetra_Type>     matrixEpetraPtr_Type;
    typedef  Epetra_Vector                            lumpedMatrix_Type;
    typedef  boost::shared_ptr<lumpedMatrix_Type>     lumpedMatrixPtr_Type;
    typedef  super::comm_Type                         comm_Type;
    typedef  super::commPtr_Type                      commPtr_Type;
    typedef  boost::shared_ptr<Teuchos::ParameterList> parameterListPtr_Type;
    typedef  MapEpetra                                mapEpetra_Type;
    typedef  boost::shared_ptr<mapEpetra_Type>        mapEpetraPtr_Type;
    typedef  VectorEpetra                             VectorEpetra_Type;
    typedef  boost::shared_ptr<VectorEpetra_Type>     VectorEpetraPtr_Type;
    //@}
    
    //! @name Constructors
    //@{
    //! Empty constructor
    aSIMPLEOperator();
    //@}
    virtual ~aSIMPLEOperator();
    
    //! @name SetUp
    //@{
    //! SetUp
    /*!
     */
    void setUp(const matrixEpetraPtr_Type & F,
               const matrixEpetraPtr_Type & B,
               const matrixEpetraPtr_Type & Btranspose,
               const commPtr_Type & comm);
    
    //! @name Set Methods
    //@{
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
    
    void updateApproximatedMomentumOperator();
    
    void updateApproximatedSchurComplementOperator();
    
    void setMomentumOptions(const parameterListPtr_Type & _oList);
    
    void setSchurOptions(const parameterListPtr_Type & _oList);

    //! Show information about the class
    void showMe();
    
    void setOptions(const Teuchos::ParameterList& solversOptions);
    
private:
    
    //! Create the domain and the range maps
    void setMaps();
    
    //! create the matrix B*diag(F)^-1*Btranspose
    void buildShurComplement();
    
    boost::shared_ptr<BlockEpetra_Map> M_operatorDomainMap;
    //! Range Map
    boost::shared_ptr<BlockEpetra_Map> M_operatorRangeMap;
    
    matrixEpetraPtr_Type M_F;

    matrixEpetraPtr_Type M_B;

    matrixEpetraPtr_Type M_Btranspose;
    
    matrixEpetraPtr_Type M_schurComplement;
    
    //! Communicator
    commPtr_Type M_comm;
    
    bool M_useTranspose;
    
    Operators::ApproximatedInvertibleRowMatrix * M_approximatedMomentumOperator;
    
    Operators::ApproximatedInvertibleRowMatrix * M_approximatedSchurComplementOperator;
    
    parameterListPtr_Type M_momentumOptions;
    
    parameterListPtr_Type M_schurOptions;
    
    mapEpetraPtr_Type M_monolithicMap;
    
    boost::shared_ptr<Epetra_Vector> M_invD;
    
    //! Label
    const std::string M_label;
};
    
} /* end namespace Operators */
} //end namespace
#endif