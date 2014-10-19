//@HEADER

/*!
* @file aSIMPLE.hpp
* @brief A class for high-order-splitting on Navier-Stokes.
* @author Umberto Villa <uvilla@emory.edu>
* @date 12-12-2009
*
* This file contains an implementation of the High Order Yosida preconditioner for the unsteady Navier-Stokes
* equations as described in A. Veneziani, U. Villa <i> ALADINS: an ALgebraic splitting time ADaptive solver for the Incompressible Navier-Stokes equations. </i>
*/

#include <boost/numeric/ublas/matrix.hpp>

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/core/linear_algebra/LinearOperator.hpp>

#ifndef _ASIMPLE_H_
#define _ASIMPLE_H_

namespace LifeV{
namespace Operators
{

class aSIMPLE: public LinearOperator
{
public:
    //! @name Public Types
    //@{
    typedef boost::numeric::ublas::vector< boost::shared_ptr<vector_Type> > Zdata;
    typedef boost::numeric::ublas::matrix< boost::shared_ptr<vector_Type> > ZZdata;
    
    typedef  Epetra_CrsMatrix                     matrix_Type;
    typedef  boost::shared_ptr<matrix_Type>       matrixPtr_Type;
    typedef  Epetra_Vector                        lumpedMatrix_Type;
    typedef  boost::shared_ptr<lumpedMatrix_Type> lumpedMatrixPtr_Type;
    //@}
    
    //! @name Constructors
    //@{
    //! Empty constructor
    aSIMPLE();
    //@}
    virtual ~aSIMPLE();
    
    //! @name SetUp
    //@{
    
    //! SetUp
    /*!
     @param Bt pointer to the pressure gradient matrix
     @param ml pointer to the lumped velocity mass matrix
     @param A  pointer to the stiffness + convection velocity matrix
     @param DL discrete laplacian object
     @param p  order of the Yosida Pressure Correction ( error = dt^(p+2) )
     */
    void SetUp(const matrixPtr_Type & Bt,
               const lumpedMatrixPtr_Type & ml,
               const matrixPtr_Type & A,
               const operatorPtr_Type & DL,
               int p);
    //@}
    
    //! @name Set Methods
    //@{
    //  //! Set the communicator
    //  void setComm(const boost::shared_ptr<Epetra_Comm> & comm) {M_comm = comm;}
    //! Set the constant @f$ \sigma = \frac{\alpha_0}{\Delta t}@f$,  @f$\alpha_0@f$ being the leading bdf coefficient.
    void setSigma(Real sigma);
    //! set the stiffness matrix
    void setA(const matrixPtr_Type & A);
    //! Set the order of the splitting
    void setSplittingOrder(int p);
    //! \warning Transpose of this operator is not supported
    int SetUseTranspose(bool UseTranspose){M_useTranspose = UseTranspose; return 0;}
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
    const comm_Type & Comm() const {return M_DL->Comm();}
    //! Returns the Epetra_Map object associated with the domain of this operator
    const map_Type & OperatorDomainMap() const {return M_Bt->OperatorDomainMap();}
    //! Returns the Epetra_Map object associated with the range of this operator
    const map_Type & OperatorRangeMap() const {return M_Bt->OperatorDomainMap();}
    //@}
    
    
    //! Show information about the class
    void showMe();
    
private:
    
    //! Solve the Approximate Shur-Complement S = 1/sigma*BHB'
    void solveS(const vector_Type & X, vector_Type & Y) const;
    
    //! Label
    const std::string M_label;
    //! Discrete Laplacian
    operatorPtr_Type M_DL;
    //! stiff(u) + convection(u)
    matrixPtr_Type M_A;
    //! Pressure divergence
    matrixPtr_Type M_Bt;
    //! velocity lumped mass "matrix"
    lumpedMatrixPtr_Type     M_ml;
    //! inverse velocity lumped mass "matrix"
    lumpedMatrixPtr_Type     M_h;
    //! \alpha_0/delta_t from bdf
    Real M_sigma;
    //! number of corrections
    UInt M_p;
    
    bool M_useTranspose;
};
    
} /* end namespace Operators */
} //end namespace
#endif
