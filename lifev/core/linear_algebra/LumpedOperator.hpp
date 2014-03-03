/*
 * LumpedOperator.hpp
 *
 * An Operator for Mass lumping:
 * \hat{M}_{i,i} = M_{i,i}* \frac{\sum_{i,j} M_{i,j}}{\sum_i{M_{i,i}}
 *
 *  Created on: Sep 17, 2010
 *      Author: uvilla
 */

#ifndef LUMPED_OPERATOR_HPP_
#define LUMPED_OPERATOR_HPP_

#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>

#include <lifev/operator/linear_algebra/LinearOperator.hpp>

namespace LifeV
{
namespace Operators
{
//! @class LumpedOperator
/*!
 * @brief A class for Mass Lumping
 * The lumped mass matrix is such that
 * (1) The total mass is conserved
 * (2) \hat{M}_{i,i} is proportional to M_{i,i}.
 *
 * The formula for the lumping is
 * \hat{M}_{i,i} = M_{i,i}* \frac{\sum_{i,j} M_{i,j}}{\sum_i{M_{i,i}}
 *
 * In the case of P1 F.E. such lumping is in agreement with the
 * one obtained using the Trapezoidal quadrature rule.
 */

class LumpedOperator: public LinearOperator
{
public:

    typedef Epetra_CrsMatrix rowMatrix_Type;
    typedef boost::shared_ptr<rowMatrix_Type> rowMatrixPtr_Type;
    typedef Epetra_Vector lumpedMatrix_Type;
    typedef boost::shared_ptr<lumpedMatrix_Type> lumpedMatrixPtr_Type;

    //! Constructor
    LumpedOperator();
    //! SetUp
    void setUp(const rowMatrixPtr_Type & _matrix);
    //! set the scaling constant when applying the operator
    void setAlpha(const Real & _alpha){ M_alpha = _alpha;}
    //! compute the lumpded operator
    int compute();

    //! if true apply the traspose of the operator
    int SetUseTranspose(bool UseTranspose) {M_useTranspose = UseTranspose; return 0;};
    //! Apply the Operator
    int Apply(const vector_Type& X, vector_Type& Y) const;
    //! Apply the inverse of the operator
    int ApplyInverse(const vector_Type& X, vector_Type& Y) const;
    //! return the infinity norm of the operator
    double NormInf() const {return M_normInf;}

    //! return the name of the operator
    const char * Label() const{return M_name.c_str();}
    //! return true if we are using the transpose
    bool UseTranspose() const {return M_useTranspose;}
    //! return true if the operator has the norm inf.
    bool HasNormInf() const {return true;}
    //! return the communicator
    const comm_Type & Comm() const {return M_matrix->Comm();}
    //! return the domain map
    const map_Type & OperatorDomainMap() const{return M_matrix->OperatorDomainMap();}
    //! return the range map
    const map_Type & OperatorRangeMap() const {return M_matrix->OperatorRangeMap();}
    //! return the lumped mass in vectorial form.
    const lumpedMatrix_Type & getLumpedVector(){return *M_lumped;}
    //! return the lumped mass in vectorial form (pointer).
    const lumpedMatrixPtr_Type & getLumpedVector_ptr(){return M_lumped;}

private:
    //! the name of the operator
    std::string M_name;
    //! the matrix to be lumped
    rowMatrixPtr_Type M_matrix;
    //! the lumped matrix in vectorial form
    lumpedMatrixPtr_Type M_lumped;
    //! the infinity norm of the operator
    Real M_normInf;
    //! whenever to use the transpose
    bool M_useTranspose;
    //! the scalar coefficient for scaling the operator
    Real M_alpha;
};
} /*end namespace Operators*/
} /*end namespace */
#endif /* LUMPEDOPERATOR_HPP_ */
