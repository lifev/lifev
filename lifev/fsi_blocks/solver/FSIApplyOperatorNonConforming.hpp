//@HEADER

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/BlockEpetra_Map.hpp>
#include <lifev/core/linear_algebra/LinearOperatorAlgebra.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesPreconditionerOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aSIMPLEOperator.hpp>
#include <lifev/navier_stokes_blocks/solver/aPCDOperator.hpp>

#include <lifev/navier_stokes_blocks/solver/NavierStokesOperator.hpp>

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/linear_algebra/ApproximatedInvertibleRowMatrix.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <lifev/core/interpolation/RBFInterpolation.hpp>
#include <lifev/core/interpolation/RBFlocallyRescaledVectorial.hpp>

#ifndef _FSIAPPLYOPERATORNONCONFORMING_H_
#define _FSIAPPLYOPERATORNONCONFORMING_H_

namespace LifeV
{
namespace Operators
{
class FSIApplyOperatorNonConforming: public LinearOperatorAlgebra
{
public:
    //! @name Public Types
    //@{

    typedef  LinearOperatorAlgebra                     super;
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
    typedef  boost::shared_ptr<RBFInterpolation<RegionMesh<LinearTetra> > > interpolationPtr_Type;
    //@}

    //! @name Constructors
    //@{
    //! Empty constructor
    FSIApplyOperatorNonConforming();
    //@}
    virtual ~FSIApplyOperatorNonConforming();

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
    //! Returns the Application of the Jacobian applied to \c X.
    int Apply(const vector_Type &X, vector_Type &Y) const;
    //! \warning No method \c Apply defined for this operator. It return an error code.
    int ApplyInverse(const vector_Type &/*X*/, vector_Type &/*Y*/) const {return -1;};
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

    // @name Set
    //@{

    //! Set the monolithic map
    void setMonolithicMap(const mapEpetraPtr_Type& monolithicMap);

    //@}

    // @name Update approximations of the block preconditioners
    //@{


    //@}

    //! Show information about the class
    void showMe();

    // @name Set the parameters
    //@{

    //@}

private:

    // @name Set the parameters
    //@{


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

    //@}

    // @name Parameter list of each block
    //@{

    //@}

    // @name Parameter list of each block
    //@{

    //! Structure block

    //@}

    mapEpetraPtr_Type M_monolithicMap;


    //! Label
    const std::string M_label;

    //! Vectors needed for the Apply - input vectors associated to each part of the residual
    boost::shared_ptr<VectorEpetra_Type > M_X_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_X_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_X_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_X_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_X_geometry;

    //! Vectors needed for the Apply - output vectors associated to the application of the
    // Jacobian to each part of the residual
    boost::shared_ptr<VectorEpetra_Type > M_Y_velocity;
    boost::shared_ptr<VectorEpetra_Type > M_Y_pressure;
    boost::shared_ptr<VectorEpetra_Type > M_Y_displacement;
    boost::shared_ptr<VectorEpetra_Type > M_Y_lambda;
    boost::shared_ptr<VectorEpetra_Type > M_Y_geometry;

};

} /* end namespace Operators */
} //end namespace
#endif
