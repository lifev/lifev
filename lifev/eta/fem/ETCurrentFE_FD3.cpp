#include <lifev/eta/fem/ETCurrentFE.hpp>


namespace LifeV
{

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<1, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0];
}

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<2, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] = M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                             - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][0][1];
}

// Full specialization for the computation of the determinant
template<>
void
ETCurrentFE<3, 3>::
updateDetJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its determinant");

#ifndef NDEBUG
    M_isDetJacobianUpdated = true;
#endif

    M_detJacobian[iQuadPt] =
        M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
        + M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
        + M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]
        - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]
        - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0];
}


template<>
void
ETCurrentFE<1, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    M_tInverseJacobian[iQuadPt][0][0] = 1.0 / M_jacobian[iQuadPt][0][0];
}


template<>
void
ETCurrentFE<2, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{

    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] =  M_jacobian[iQuadPt][0][0] / det;
    M_tInverseJacobian[iQuadPt][1][0] = -M_jacobian[iQuadPt][1][0] / det;
    M_tInverseJacobian[iQuadPt][0][1] = -M_jacobian[iQuadPt][0][1] / det;
    M_tInverseJacobian[iQuadPt][1][1] =  M_jacobian[iQuadPt][1][1] / det;
}

template<>
void
ETCurrentFE<3, 3>::
updateInverseJacobian (const UInt& iQuadPt)
{
    ASSERT (M_isJacobianUpdated, "Jacobian must be updated to compute its inverse");
    ASSERT (M_isDetJacobianUpdated, "The determinant of the jacobian must be updated to compute its inverse");

#ifndef NDEBUG
    M_isInverseJacobianUpdated = true;
#endif

    Real det = M_detJacobian[iQuadPt];

    M_tInverseJacobian[iQuadPt][0][0] = ( M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][0][1] = ( M_jacobian[iQuadPt][1][2] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][0][2] = ( M_jacobian[iQuadPt][1][0] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][1][1] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][0] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][2]) / det;

    M_tInverseJacobian[iQuadPt][1][1] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][2][0]) / det;

    M_tInverseJacobian[iQuadPt][1][2] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][2][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][2][1]) / det;

    M_tInverseJacobian[iQuadPt][2][0] = ( M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][2]
                                          - M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][1]) / det;

    M_tInverseJacobian[iQuadPt][2][1] = ( M_jacobian[iQuadPt][0][2] * M_jacobian[iQuadPt][1][0]
                                          - M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][2]) / det;

    M_tInverseJacobian[iQuadPt][2][2] = ( M_jacobian[iQuadPt][0][0] * M_jacobian[iQuadPt][1][1]
                                          - M_jacobian[iQuadPt][0][1] * M_jacobian[iQuadPt][1][0]) / det;
}


} //Namespace LifeV
