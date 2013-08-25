//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
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
 *  @brief This file contains the definition for the St. Venant Kirchhoff linear material
 *
 *  @version 1.0
 *  @date 29-07-2010
 *  @author Paolo Tricerri,
 *
 *  @maintainer  Paolo Tricerri      <paolo.tricerri@epfl.ch>
 */

#ifndef _SECONDORDEREXPONENTIALMATERIAL_H_
#define _SECONDORDEREXPONENTIALMATERIAL_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <lifev/structure/solver/isotropic/StructuralIsotropicConstitutiveLaw.hpp>

namespace LifeV
{
template <typename MeshType>
class SecondOrderExponentialMaterialNonLinear : public StructuralIsotropicConstitutiveLaw<MeshType>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralIsotropicConstitutiveLaw<MeshType>       super;

    typedef typename super::data_Type                data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;

    typedef typename super::vectorVolumes_Type       vectorVolumes_Type;
    typedef boost::shared_ptr<vectorVolumes_Type>    vectorVolumesPtr_Type;

    typedef typename super::mapMarkerIndexesPtr_Type mapMarkerIndexesPtr_Type;
    typedef typename super::mapMarkerIndexes_Type    mapMarkerIndexes_Type;
    typedef typename mapMarkerIndexes_Type::const_iterator mapIteratorIndex_Type;

    typedef std::vector<UInt>                        vectorIndexes_Type;
    typedef boost::shared_ptr<vectorIndexes_Type>    vectorIndexesPtr_Type;

    typedef typename super::FESpacePtr_Type          FESpacePtr_Type;
    typedef typename super::ETFESpacePtr_Type        ETFESpacePtr_Type;

    //Vector for vector parameters
    typedef typename super::vectorsParameters_Type       vectorsParameters_Type;
    typedef typename super::vectorsParametersPtr_Type    vectorsParametersPtr_Type;

    typedef MatrixSmall<3, 3>                          matrixSmall_Type;

    // Typedefs for expression definitions
    typedef typename super::tensorF_Type               tensorF_Type;
    typedef typename super::determinantF_Type          determinantF_Type;
    typedef typename super::tensorC_Type               tensorC_Type;
    typedef typename super::minusT_Type                minusT_Type;
    typedef typename super::traceTensor_Type           traceTensor_Type;
    //@}



    //! @name Constructor &  Destructor
    //@{

    SecondOrderExponentialMaterialNonLinear();

    virtual  ~SecondOrderExponentialMaterialNonLinear();

    //@}



    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralIsotropicConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup ( const FESpacePtr_Type&                       dFESpace,
                 const ETFESpacePtr_Type&                     dETFESpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial);


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/, const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                const displayerPtr_Type& displayer );


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                        const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                        const displayerPtr_Type& displayer );


    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& disp, Real factor, const dataPtr_Type& dataMaterial, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                            const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/, const displayerPtr_Type& displayer );


    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    // void computeVector( const vector_Type& sol,
    //                     Real factor,
    //                     const dataPtr_Type& dataMaterial,
    //                     const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
    //                     const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
    //                     const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F), the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralIsotropicConstitutiveLaw::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc );

    //! Computes the deformation Gradient F, the cofactor of F Cof(F), the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    //void computeStress( const vector_Type& sol );

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain );

    //! Compute the First Piola Kirchhoff Tensor
    /*!
       \param firstPiola Epetra_SerialDenseMatrix that has to be filled
       \param tensorF Epetra_SerialDenseMatrix the deformation gradient
       \param cofactorF Epetra_SerialDenseMatrix cofactor of F
       \param invariants std::vector with the invariants of C and the detF
       \param material UInt number to get the material parameteres form the VenantElasticData class
    */
    void computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
                                                 const Epetra_SerialDenseMatrix& tensorF,
                                                 const Epetra_SerialDenseMatrix& cofactorF,
                                                 const std::vector<Real>& invariants,
                                                 const UInt marker);

    //! Compute the First Piola Kirchhoff Tensor
    /*!
       \param disp the displacement field from which we compute the fisrt piola-Kirchhoff tensor
       \param sigma_1 the first column of the Cauchy stress tensor
       \param sigma_2 the second column of the Cauchy stress tensor
       \param sigma_3 the third column of the Cauchy stress tensor
    */
    void computeCauchyStressTensor ( const vectorPtr_Type disp,
				     const QuadratureRule& evalQuad,
				     vectorPtr_Type sigma_1,
				     vectorPtr_Type sigma_2,
				     vectorPtr_Type sigma_3);

    //@}

    //! @name Get Methods
    //@{

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const
    {
        return super::M_jacobian;
    }


    //! Get the stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        return M_stiff;
    }

    void apply ( const vector_Type& sol, vector_Type& res,
                 const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                 const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                 const displayerPtr_Type displayer) ;

    //@}



protected:

    //! construct the vectors for the parameters
    /*!
      \param VOID
      \return VOID
    */
    void setupVectorsParameters ( void );

    //! Vector: stiffness non-linear
    vectorPtr_Type                         M_stiff;

    //Create the indentity for F
    matrixSmall_Type                      M_identity;
};





template <typename MeshType>
SecondOrderExponentialMaterialNonLinear<MeshType>::SecondOrderExponentialMaterialNonLinear() :
    super       ( ),
    M_stiff     ( ),
    M_identity ( )
{
}





template <typename MeshType>
SecondOrderExponentialMaterialNonLinear<MeshType>::~SecondOrderExponentialMaterialNonLinear()
{}





template <typename MeshType>
void
SecondOrderExponentialMaterialNonLinear<MeshType>::setup (const FESpacePtr_Type&                       dFESpace,
                                                          const ETFESpacePtr_Type&                     dETFESpace,
                                                          const boost::shared_ptr<const MapEpetra>&    monolithicMap,
                                                          const UInt                                   offset,
                                                          const dataPtr_Type& dataMaterial)
{

    this->M_dataMaterial  = dataMaterial;
    //    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_dispFESpace                     = dFESpace;
    this->M_dispETFESpace                = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;

    M_stiff.reset                     ( new vector_Type (*this->M_localMap) );

    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;

    // The 3 is because the law uses three parameters (alpha, gamma, bulk).
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 3 ) );

    this->setupVectorsParameters();
}


template <typename MeshType>
void
SecondOrderExponentialMaterialNonLinear<MeshType>::setupVectorsParameters ( void )
{
    // Paolo Tricerri: February, 20th
    // In each class, the name of the parameters has to inserted in the law
    // TODO: move the saving of the material parameters to more abstract objects
    //       such that in the class of the material we do not need to call explicitly
    //       the name of the parameter.

    // Number of volume on the local part of the mesh
    UInt nbElements = this->M_dispFESpace->mesh()->numVolumes();

    // Parameter alpha
    // 1. resize the vector in the first element of the vector.
    (* (this->M_vectorsParameters) ) [0].resize ( nbElements );

    // Parameter gamma
    (* (this->M_vectorsParameters) ) [1].resize ( nbElements );

    // Parameter bulk
    (* (this->M_vectorsParameters) ) [2].resize ( nbElements );

    for (UInt i (0); i < nbElements; i++ )
    {
        // Extracting the marker
        UInt markerID = this->M_dispFESpace->mesh()->element ( i ).markerID();

        Real alpha = this->M_dataMaterial->alpha ( markerID );
        Real gamma = this->M_dataMaterial->gamma ( markerID );
        Real bulk = this->M_dataMaterial->bulk ( markerID );

        ( (* (this->M_vectorsParameters) ) [0]) [ i ] = alpha;
        ( (* (this->M_vectorsParameters) ) [1]) [ i ] = gamma;
        ( (* (this->M_vectorsParameters) ) [2]) [ i ] = bulk;
    }
}


template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/, const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    //! Empty method for second order exponential material
}

template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                               const dataPtr_Type&      dataMaterial,
                                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                               const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n************************************************\n ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes,  displayer);
    displayer->leaderPrint (" \n************************************************\n ");
}





template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type&         jacobian,
        const vector_Type&     disp,
        const dataPtr_Type&     dataMaterial,
        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
        const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
        const displayerPtr_Type& displayer )
{

    using namespace ExpressionAssembly;

    *jacobian *= 0.0;

    displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Second Order Exponential)");

    //! Nonlinear part of jacobian
    //! loop on volumes (i)

    // mapIterator_Type it;
    // //mapIteratorIndex_Type itIndex;

    // vectorVolumesPtr_Type pointerListOfVolumes;
    // vectorIndexesPtr_Type pointerListOfIndexes;

    // for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    // {

    //     //Given the marker pointed by the iterator, let's extract the material parameters
    //     UInt marker = it->first;

    //     // Debug
    //     // UInt markerIndex = itIndex->first;
    //     // ASSERT( marker == markerIndex, "The list of volumes is referring to a marker that is not the same as the marker of index!!!");

    //     pointerListOfVolumes.reset( new vectorVolumes_Type(it->second) );
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker] ) );

    //     Real bulk = dataMaterial->bulk(marker);
    //     Real alpha = dataMaterial->alpha(marker);
    //     Real gamma = dataMaterial->gamma(marker);

    // Definition of F
    tensorF_Type F = ExpressionDefinitions::deformationGradient( this->M_dispETFESpace,  disp, this->M_offset, this->M_identity );

    // Definition of J
    determinantF_Type J = ExpressionDefinitions::determinantF( F );

    // Definition of tensor C
    tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

    // Definition of F^-T
    minusT_Type  F_T = ExpressionDefinitions::minusT( F );

    // Definition of tr( C )
    traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

    //! VOLUMETRIC PART
    //! 1. Stiffness matrix: int { 1/2 * bulk * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [2] ) * ( value (2.0) *pow (J, 2.0) - J + value (1.0) ) * dot ( F_T, grad (phi_j) ) * dot ( F_T, grad (phi_i) )
              ) >> jacobian;

    //! 2. Stiffness matrix: int { 1/2 * bulk * ( 1/J- 1 - log(J)/J^2 ) * ( CofF1 [\nabla \delta]^t CofF ) : \nabla v }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( - 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [2] ) * ( pow (J, 2.0) - J + log (J) ) * dot ( F_T * transpose (grad (phi_j) ) * F_T,  grad (phi_i) )
              ) >> jacobian;


    //! ISOCHORIC PART
    //! 1. Stiffness matrix : int { (- 4/3 * alpha * J^(-2/3) * exp( gamma*( Ic_iso - 3)^2 )*
    //!    * ( Ic_isoK + ( Ic_isoK - 3 ) * ( 2 * gamma * ( Ic_isoK - 3) +1 ) ) ) * ( F^-T : \nabla \delta ) ( F : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( -4.0 / 3.0 ) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow ( J, -2.0 / 3.0 ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 )* I_C - value (3.0) ) ) *
                ( pow( J, -2.0/3.0 ) * I_C + ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( value ( 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) + value (1.0) ) ) *
                ( dot ( F_T, grad (phi_j) ) * dot ( F, grad (phi_i) ) )
              ) >> jacobian;


    //! 2. Stiffness matrix : int { 4 * alpha * J^(-4/3) * exp( gamma*( Ic_iso - 3)^2 ) * ( 1.0 + 2 * gamma * (Ic_isoK - 3 )^2 ) *
    //!                           ( F : \nabla \delta ) ( F : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 4.0 ) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow ( J, (-4.0 / 3.0) ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) *
                ( value (1.0) + value ( 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) *
                ( dot ( F, grad (phi_j) ) * dot ( F, grad (phi_i) ) )
              ) >> jacobian;

    //! 3. Stiffness matrix :
    //!int { ( 4.0/9.0 *  alpha * Ic_isoK *  exp( gamma*( Ic_iso - 3)^2 ) * ( Ic_isoK + 2 * gamma * (Ic_isoK - 3)^2 * Ic_isoK + (Ic_isoK - 3) ) ) *
    //!     ( F^-T : \nabla \delta ) ( F^-T : \nabla \v )}
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (4.0 / 9.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) *  pow( J, -2.0/3.0 ) * I_C *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0) * I_C - value (3.0) ) ) *
                ( ( pow( J, -2.0/3.0) * I_C - value (3.0) ) * ( value (1.0) + value (2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * pow(J, -2.0/3.0) * I_C * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) + pow( J, -2.0/3.0 ) * I_C ) *
                ( dot ( F_T, grad (phi_j) ) * dot ( F_T, grad (phi_i) ) )
              ) >> jacobian;


    //! 4. Stiffness matrix :
    //! int { (-4.0/3.0 *  alpha * J^(-2/3) * exp( gamma*( Ic_iso - 3)*( Ic_iso - 3) ) * ( Ic_isoK + 2*gamma*(Ic_isok - 3)Ic + 1) ) *
    //!       ( F : \nabla \delta ) ( F^-T : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (-4.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow ( J, (-2.0 / 3.0) ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) *
                ( pow( J, -2.0/3.0 ) * I_C  + ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( value ( 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * pow( J, -2.0/3.0 ) * I_C + value (1.0) ) ) *
                ( dot ( F, grad (phi_j) ) * dot ( F_T, grad (phi_i) ) )
              ) >> jacobian;

    //! 5. Stiffness matrix : int { (2 * alpha * J^(-2/3) * exp( gamma*( Ic_iso - 3)*( Ic_iso - 3))*(Ic_iso - 3) ) * (\nabla \delta: \nabla \v)}
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow ( J, (-2.0 / 3.0) ) * ( pow(J, -2.0/3.0) * I_C - value (3.0) ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow(J, -2.0/3.0) * I_C - value (3.0) ) * ( pow(J, -2.0/3.0) * I_C - value (3.0) ) ) *
                dot ( grad (phi_j), grad (phi_i) )
              ) >> jacobian;

    //! 6. Stiffness matrix : int { 2.0/3.0 * alpha * Ic_iso * ( Ic_iso - 3) * exp(gamma*( Ic_iso - 3)*( Ic_iso - 3)) *
    //!                       (CofF [\nabla \delta]^t CofF ) : \nabla \v  }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (2.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow( J, -2.0/3.0 ) * I_C * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) *
                dot ( F_T * transpose (grad (phi_j) ) * F_T , grad (phi_i) )
              ) >> jacobian;


    //    }

    jacobian->globalAssemble();
}





template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::computeStiffness ( const vector_Type& disp,
                                                                           Real /*factor*/,
                                                                           const dataPtr_Type& dataMaterial,
                                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                           const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                           const displayerPtr_Type& displayer )
{
    using namespace ExpressionAssembly;

    this->M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

    displayer->leaderPrint (" \n*********************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the Second Order Exponential nonlinear stiffness vector ");
    displayer->leaderPrint (" \n*********************************\n  ");


    // mapIterator_Type it;
    // //mapIteratorIndex_Type itIndex;

    // vectorVolumesPtr_Type pointerListOfVolumes;
    // vectorIndexesPtr_Type pointerListOfIndexes;


    // for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++)
    // {

    //     //Given the marker pointed by the iterator, let's extract the material parameters
    //     UInt marker = it->first;

    //     // Debug
    //     // UInt markerIndex = itIndex->first;
    //     // ASSERT( marker == markerIndex, "The list of volumes is referring to a marker that is not the same as the marker of index!!!");

    //     pointerListOfVolumes.reset( new vectorVolumes_Type(it->second) );
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker] ) );


    //     Real bulk = dataMaterial->bulk(marker);
    //     Real alpha = dataMaterial->alpha(marker);
    //     Real gamma = dataMaterial->gamma(marker);

    // Definition of F
    tensorF_Type F = ExpressionDefinitions::deformationGradient( this->M_dispETFESpace,  disp, this->M_offset, this->M_identity );

    // Definition of J
    determinantF_Type J = ExpressionDefinitions::determinantF( F );

    // Definition of tensor C
    tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

    // Definition of F^-T
    minusT_Type  F_T = ExpressionDefinitions::minusT( F );

    // Definition of tr( C )
    traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

    //! Stiffness for non-linear terms of the Neo-Hookean model
    /*!
      The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
    */
    //! Volumetric part
    /*!
      Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
    */
    //! Stiffness for non-linear terms of the Neo-Hookean model
    //! Volumetric part : int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [2] ) * ( pow ( J , 2.0) - J + log (J) ) * dot ( F_T, grad (phi_i) )
              ) >> M_stiff;

    //! Isochoric part
    /*!
      Source term P1iso_Exp: int { 2 * alpha * ( Ic1_iso - 3 ) * exp(gamma *(  Ic1_iso -3 )^2) *
      ( J1^(-2/3)* (F1 : \nabla v) - 1/3 * (Ic1_iso / J1) * (CofF1 : \nabla v) ) }
    */
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                value ( 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * ( pow( J, -2.0/3.0) * I_C - value (3.0) ) *
                exp ( parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) * ( pow( J, -2.0/3.0 ) * I_C - value (3.0) ) ) *
                pow ( J, (-2.0 / 3.0) ) * dot ( F - value (1.0 / 3.0) * pow( J, -2.0/3.0 ) * I_C * F_T , grad (phi_i) )
              ) >> M_stiff;



    //    }

    this->M_stiff->globalAssemble();
}





template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::computeKinematicsVariables ( const VectorElemental& dk_loc )
{

    // Real s;

    // //! loop on quadrature points (ig)
    // for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
    // {
    //     //! loop on space coordinates (icoor)
    //     for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
    //     {
    //         //! loop  on space coordinates (jcoor)
    //         for ( UInt jcoor = 0; jcoor < nDimensions; jcoor++ )
    //         {
    //             s = 0.0;
    //             for ( UInt i = 0; i < this->M_dispFESpace->fe().nbFEDof(); i++ )
    //             {
    //                 //! \grad u^k at a quadrature point
    //                 s += this->M_dispFESpace->fe().phiDer( i, jcoor, ig ) * dk_loc[ i + icoor * this->M_dispFESpace->fe().nbFEDof() ];
    //             }
    //             //! gradient of displacement
    //             (*M_Fk)[ icoor ][ jcoor ][ig ] = s;
    //         }
    //     }
    // }

    // //! loop on quadrature points (ig)
    // for ( UInt ig = 0; ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
    // {
    //     //! loop on space coordinates (icoor)
    //     for ( UInt  icoor = 0;icoor < nDimensions; icoor++ )
    //     {
    //         //! deformation gradient Fk
    //         (*M_Fk)[ icoor ][ icoor ][ ig ] +=  1.0;
    //     }
    // }

    // Real a,b,c,d,e,f,g,h,i;

    // for( UInt ig=0; ig< this->M_dispFESpace->fe().nbQuadPt(); ig++ )
    // {
    //     a = (*M_Fk)[ 0 ][ 0 ][ ig ];
    //     b = (*M_Fk)[ 0 ][ 1 ][ ig ];
    //     c = (*M_Fk)[ 0 ][ 2 ][ ig ];
    //     d = (*M_Fk)[ 1 ][ 0 ][ ig ];
    //     e = (*M_Fk)[ 1 ][ 1 ][ ig ];
    //     f = (*M_Fk)[ 1 ][ 2 ][ ig ];
    //     g = (*M_Fk)[ 2 ][ 0 ][ ig ];
    //     h = (*M_Fk)[ 2 ][ 1 ][ ig ];
    //     i = (*M_Fk)[ 2 ][ 2 ][ ig ];

    //     //! determinant of deformation gradient Fk
    //     (*M_Jack)[ig] = a*( e*i - f*h ) - b*( d*i - f*g ) + c*( d*h - e*g );

    //     ASSERT_PRE((*M_Jack)[ig] > 0, "Negative Jacobian. Error!" );

    //     (*M_CofFk)[ 0 ][ 0 ][ ig ] =   ( e*i - f*h );
    //     (*M_CofFk)[ 0 ][ 1 ][ ig ] = - ( d*i - g*f );
    //     (*M_CofFk)[ 0 ][ 2 ][ ig ] =   ( d*h - e*g );
    //     (*M_CofFk)[ 1 ][ 0 ][ ig ] = - ( b*i - c*h );
    //     (*M_CofFk)[ 1 ][ 1 ][ ig ] =   ( a*i - c*g );
    //     (*M_CofFk)[ 1 ][ 2 ][ ig ] = - ( a*h - g*b );
    //     (*M_CofFk)[ 2 ][ 0 ][ ig ] =   ( b*f - c*e );
    //     (*M_CofFk)[ 2 ][ 1 ][ ig ] = - ( a*f - c*d );
    //     (*M_CofFk)[ 2 ][ 2 ][ ig ] =   ( a*e - d*b );
    // }

    // //! loop on quadrature points
    // for ( UInt ig = 0;ig < this->M_dispFESpace->fe().nbQuadPt(); ig++ )
    // {
    //     s = 0.0;
    //     for ( UInt i = 0; i < nDimensions; i++)
    //     {
    //         for ( UInt j = 0; j < nDimensions; j++)
    //         {
    //             //! trace of  C1 = (F1k^t F1k)
    //             s +=  (*M_Fk)[ i ][ j ][ ig ] * (*M_Fk)[ i ][ j ][ ig ];
    //         }
    //     }
    //     (*M_trCk)[ ig ] = s;
    // }

    // for ( UInt ig = 0; ig <  this->M_dispFESpace->fe().nbQuadPt(); ig++ )
    // {
    //     //! trace of deviatoric C
    //     (*M_trCisok)[ ig ] =  std::pow((*M_Jack)[ ig ], -2./3.) * (*M_trCk)[ ig ];
    // }
}





template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::showMe ( std::string const& fileNameStiff,
                                                                 std::string const& fileNameJacobian )
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}


template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::apply ( const vector_Type& sol, vector_Type& res, const mapMarkerVolumesPtr_Type mapsMarkerVolumes,  const mapMarkerIndexesPtr_Type mapsMarkerIndexes,const displayerPtr_Type displayer )
{
    computeStiffness (sol, 0, this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    res += *M_stiff;
}


template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
        const Epetra_SerialDenseMatrix& tensorF,
        const Epetra_SerialDenseMatrix& cofactorF,
        const std::vector<Real>& invariants,
        const UInt marker)
{
    //Get the material parameters
    Real alpha    = this->M_dataMaterial->alpha (marker);
    Real gamma    = this->M_dataMaterial->gamma (marker);
    Real bulk     = this->M_dataMaterial->bulk (marker);


    //Computing the first term \alphaJ^{-2/3}[F-(1/3)tr(C)F^{-T}]exp(\gamma(tr(Ciso) - 3)
    Epetra_SerialDenseMatrix firstTerm (tensorF);
    Epetra_SerialDenseMatrix copyCofactorF (cofactorF);

    Real scale (0.0);
    scale = -invariants[0] / 3.0;
    copyCofactorF.Scale ( scale );
    firstTerm += copyCofactorF;

    //Computation trace of the isochoric C
    Real trCiso (0.0);
    trCiso = std::pow (invariants[3], - (2.0 / 3.0) ) * invariants[0];

    Real coef ( 0.0 );
    coef = 2.0 * alpha * std::pow (invariants[3], - (2.0 / 3.0) ) * ( trCiso - 3.0 ) * std::exp ( gamma * ( trCiso - 3 ) * ( trCiso - 3 )  );
    firstTerm.Scale ( coef );

    //Computing the second term (volumetric part) J*(bulk/2)(J-1+(1/J)*ln(J))F^{-T}
    Epetra_SerialDenseMatrix secondTerm (cofactorF);
    Real secCoef (0);
    secCoef = invariants[3] * (bulk / 2.0) * (invariants[3] - 1 + (1.0 / invariants[3]) * std::log (invariants[3]) );

    secondTerm.Scale ( secCoef );

    firstPiola += firstTerm;
    firstPiola += secondTerm;

}

template <typename MeshType>
void SecondOrderExponentialMaterialNonLinear<MeshType>::computeCauchyStressTensor ( const vectorPtr_Type disp,
										    const QuadratureRule& evalQuad,
										    vectorPtr_Type sigma_1,
										    vectorPtr_Type sigma_2,
										    vectorPtr_Type sigma_3) 

{
  ASSERT( 2 < 0, "This method has to be implemented for the 2nd Exponential law");
}

template <typename MeshType>
inline StructuralIsotropicConstitutiveLaw<MeshType>* createSecondOrderExponentialMaterialNonLinear()
{
    return new SecondOrderExponentialMaterialNonLinear<MeshType >();
}

namespace
{
static bool registerSOEXP = StructuralIsotropicConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureIsotropicMaterialFactory::instance().registerProduct ( "secondOrderExponential", &createSecondOrderExponentialMaterialNonLinear<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __SECONDORDEREXPONENTIALMATERIAL_H */
