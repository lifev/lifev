//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/navier_stokes/solver/NavierStokesPreconditionerOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#ifndef _aPCDOPERATOR_H_
#define _aPCDOPERATOR_H_1

namespace LifeV{
namespace Operators
{

class aPCDOperator: public NavierStokesPreconditionerOperator
{
public:
    //! @name Public Types
    //@{

    typedef  LinearOperator                            super;
    typedef  Epetra_MultiVector                        vector_Type;
    typedef  boost::shared_ptr<vector_Type>            vectorPtr_Type;
    typedef  Epetra_Map                                map_Type;
    typedef  boost::shared_ptr<map_Type> 			   mapPtr_Type;
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
    aPCDOperator();
    //@}
    virtual ~aPCDOperator();

    //! @name SetUp
    //@{
    //! SetUp
    /*!
     */
    void setUp(const matrixEpetraPtr_Type & F,
               const matrixEpetraPtr_Type & B,
               const matrixEpetraPtr_Type & Btranspose,
               const matrixEpetraPtr_Type & Fp,
               const matrixEpetraPtr_Type & Mp,
               const matrixEpetraPtr_Type & Mu);

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

    //! Returns the High Order Yosida approximation of the inverse pressure Schur Complement applied to \c (Xu, Xp).
    int ApplyInverse( VectorEpetra_Type const& X_velocity,
                      VectorEpetra_Type const& X_pressure,
                      VectorEpetra_Type & Y_velocity,
                      VectorEpetra_Type & Y_pressure) const;

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

    void updateApproximatedPressureMassOperator();

    void setMomentumOptions(const parameterListPtr_Type & _oList);

    void setPressureMassOptions(const parameterListPtr_Type & _oList);

    void setSchurOptions(const parameterListPtr_Type & _oList);

    //! Show information about the class
    void showMe();

    void setOptions(const Teuchos::ParameterList& solversOptions);

    matrixEpetraPtr_Type const& F() const { return M_F; }

    matrixEpetraPtr_Type const& B() const { return M_B; }

    matrixEpetraPtr_Type const& Btranspose() const { return M_Btranspose; }

    matrixEpetraPtr_Type const& Fp() const { return M_Fp; }

    matrixEpetraPtr_Type const& Mp() const { return M_Mp; }

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

    matrixEpetraPtr_Type M_Fp;

    matrixEpetraPtr_Type M_Mp;

    matrixEpetraPtr_Type M_Mu;

    matrixEpetraPtr_Type M_schurComplement;

    //! Communicator
    commPtr_Type M_comm;

    bool M_useTranspose;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedMomentumOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedSchurComplementOperator;

    Operators::ApproximatedInvertibleRowMatrix * M_approximatedPressureMassOperator;

    parameterListPtr_Type M_momentumOptions;

    parameterListPtr_Type M_schurOptions;

    parameterListPtr_Type M_pressureMassOptions;

    mapEpetraPtr_Type M_monolithicMap;

    boost::shared_ptr<Epetra_Vector> M_invDiagMassVelocity;

    //! Label
    const std::string M_label;

    //! Vectors needed for the apply inverse
    boost::shared_ptr<VectorEpetra_Type> M_Zu;
    boost::shared_ptr<VectorEpetra_Type> M_Zp;

    boost::shared_ptr<VectorEpetra_Type> M_X_velocity;
    boost::shared_ptr<VectorEpetra_Type> M_X_pressure;
    boost::shared_ptr<VectorEpetra_Type> M_Y_velocity;
    boost::shared_ptr<VectorEpetra_Type> M_Y_pressure;

};

//! Factory create function
inline NavierStokesPreconditionerOperator * create_aPCD()
{
    return new aPCDOperator ();
}
namespace
{
static bool S_register_aPCD = NSPreconditionerFactory::instance().registerProduct ( "PCD", &create_aPCD );
}

} /* end namespace Operators */
} //end namespace
#endif
