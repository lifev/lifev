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
 *  @author Gianmarco Mengaldo
 *
 *  @maintainer  Paolo Tricerri      <paolo.tricerri@epfl.ch>
 *  @contributor Gianmarco Mengaldo  <gianmarco.mengaldo@gmail.com>
 */

#ifndef _NEOHOOKEANMATERIAL_H_
#define _NEOHOOKEANMATERIAL_H_


#include <lifev/structure/solver/isotropic/StructuralIsotropicConstitutiveLaw.hpp>

namespace LifeV
{

template <typename MeshType>
class NeoHookeanMaterialNonLinear :
    public StructuralIsotropicConstitutiveLaw<MeshType>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralIsotropicConstitutiveLaw<MeshType>          super;

    typedef StructuralConstitutiveLawData            data_Type;

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

    typedef std::vector<UInt>                        vectorIndexes_Type;
    typedef boost::shared_ptr<vectorIndexes_Type>    vectorIndexesPtr_Type;

    typedef typename super::mapMarkerIndexesPtr_Type mapMarkerIndexesPtr_Type;
    typedef typename super::mapMarkerIndexes_Type    mapMarkerIndexes_Type;
    typedef typename mapMarkerIndexes_Type::const_iterator mapIteratorIndex_Type;

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

    NeoHookeanMaterialNonLinear();

    virtual  ~NeoHookeanMaterialNonLinear();

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralIsotropicConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup ( const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type& dETFESpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial);


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                              const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                const displayerPtr_Type& displayer);


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                        const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                        const displayerPtr_Type& displayer);

    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& disp, Real factor, const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                            const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                            const displayerPtr_Type& displayer );

    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in
    //! StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem
    //! since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    // void computeVector( const vector_Type& sol,
    //                     Real factor,
    //                     const dataPtr_Type& dataMaterial,
    //                     const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
    //                     const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F),
    //! the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralIsotropicConstitutiveLaw::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc ) {}

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain);


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
                 const displayerPtr_Type displayer);

    //@}



protected:

    //! construct the vectors for the parameters
    /*!
      \param VOID
      \return VOID
    */
    void setupVectorsParameters ( void );

    //! Vector: stiffness non-linear
    vectorPtr_Type                      M_stiff;

    //Create the indentity for F
    matrixSmall_Type                      M_identity;

};

template <typename MeshType>
NeoHookeanMaterialNonLinear<MeshType>::NeoHookeanMaterialNonLinear() :
    super           ( ),
    M_stiff         ( ),
    M_identity        ( )
{
}





template <typename MeshType>
NeoHookeanMaterialNonLinear<MeshType>::~NeoHookeanMaterialNonLinear()
{}





template <typename MeshType>
void
NeoHookeanMaterialNonLinear<MeshType>::setup ( const FESpacePtr_Type&                      dFESpace,
                                               const ETFESpacePtr_Type&                    dETFESpace,
                                               const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                                               const UInt                                  offset,
                                               const dataPtr_Type& dataMaterial)
{
    this->M_dataMaterial                    = dataMaterial;
    this->M_dispFESpace                     = dFESpace;
    this->M_dispETFESpace                     = dETFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;
    M_stiff.reset                   ( new vector_Type (*this->M_localMap) );

    M_identity (0, 0) = 1.0;
    M_identity (0, 1) = 0.0;
    M_identity (0, 2) = 0.0;
    M_identity (1, 0) = 0.0;
    M_identity (1, 1) = 1.0;
    M_identity (1, 2) = 0.0;
    M_identity (2, 0) = 0.0;
    M_identity (2, 1) = 0.0;
    M_identity (2, 2) = 1.0;

    // The 2 is because the law uses three parameters (mu, bulk).
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );

    this->setupVectorsParameters();
}

template <typename MeshType>
void
NeoHookeanMaterialNonLinear<MeshType>::setupVectorsParameters ( void )
{
    // Paolo Tricerri: February, 20th
    // In each class, the name of the parameters has to inserted in the law
    // TODO: move the saving of the material parameters to more abstract objects
    //       such that in the class of the material we do not need to call explicitly
    //       the name of the parameter.

    // Number of volume on the local part of the mesh
    UInt nbElements = this->M_dispFESpace->mesh()->numVolumes();

    // Parameter mu
    // 1. resize the vector in the first element of the vector.
    (* (this->M_vectorsParameters) ) [0].resize ( nbElements );

    // Parameter bulk
    (* (this->M_vectorsParameters) ) [1].resize ( nbElements );

    for (UInt i (0); i < nbElements; i++ )
    {
        // Extracting the marker
        UInt markerID = this->M_dispFESpace->mesh()->element ( i ).markerID();

        Real mu = this->M_dataMaterial->mu ( markerID );
        Real bulk = this->M_dataMaterial->bulk ( markerID );

        ( (* (this->M_vectorsParameters) ) [0]) [ i ] = mu;
        ( (* (this->M_vectorsParameters) ) [1]) [ i ] = bulk;
    }
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/,
                                                                const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    //! Empty method for neo-hookean material
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                   const dataPtr_Type&      dataMaterial,
                                                                   const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                   const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                   const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    displayer->leaderPrint (" \n*********************************\n  ");
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type&       jacobian,
                                                                           const vector_Type&    disp,
                                                                           const dataPtr_Type&   dataMaterial,
                                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                           const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                           const displayerPtr_Type&  displayer )
{

    using namespace ExpressionAssembly;

    displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Neo-Hookean)");

    * (jacobian) *= 0.0;

    // mapIterator_Type it;
    // //mapIteratorIndex_Type itIndex;

    // vectorVolumesPtr_Type pointerListOfVolumes;
    // vectorIndexesPtr_Type pointerListOfIndexes;

    // for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++)
    //  {

    //     //Given the marker pointed by the iterator, let's extract the material parameters
    //     UInt marker = it->first;

    //     // Debug
    //     // UInt markerIndex = itIndex->first;
    //     // ASSERT( marker == markerIndex, "The list of volumes is referring to a marker that is not the same as the marker of index!!!");

    //     pointerListOfVolumes.reset( new vectorVolumes_Type(it->second) );
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker]) );

    //     Real mu     = dataMaterial->mu(marker);
    //     Real bulk   = dataMaterial->bulk(marker);

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

    //Assembling Volumetric Part
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( value (2.0) *pow (J, 2.0) - J + value (1.0) ) * dot ( F_T, grad (phi_j) ) * dot ( F_T, grad (phi_i) )
              ) >> jacobian;

    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( - 1.0 / 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow (J, 2.0) - J + log (J) ) * dot ( F_T * transpose (grad (phi_j) ) * F_T,  grad (phi_i) )
              ) >> jacobian;


    //! ISOCHORIC PART
    //! 1. Stiffness matrix : int { -2/3 * mu * J^(-2/3) *( F^-T : \nabla \delta ) ( F : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (-2.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, - (2.0 / 3.0) )  * dot ( F_T , grad (phi_j) ) * dot ( F , grad (phi_i) )
              ) >> jacobian;


    //! 2. Stiffness matrix : int { 2/9 * mu * ( Ic_iso )( F^-T : \nabla \delta ) ( F^-T : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (2.0 / 9.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) *
                pow( J, (-2.0/3.0) ) * I_C  * dot ( F_T , grad (phi_j) ) * dot ( F_T , grad (phi_i) )
              ) >> jacobian;

    //! 3. Stiffness matrix : int { mu * J^(-2/3) (\nabla \delta : \nabla \v)}
    integrate ( elements ( this->M_dispETFESpace->mesh() ) ,
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, - (2.0 / 3.0) )  * dot ( grad (phi_j), grad (phi_i) )
              ) >> jacobian;

    //! 4. Stiffness matrix : int { -2/3 * mu * J^(-2/3) ( F : \nabla \delta ) ( F^-T : \nabla \v ) }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (-2.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) *
                pow (J, - (2.0 / 3.0) )  * dot ( F , grad (phi_j) ) * dot ( F_T , grad (phi_i) )
              ) >> jacobian;



    //! 5. Stiffness matrix : int { 1/3 * mu * Ic_iso * (F^-T [\nabla \delta]^t F^-T ) : \nabla \v }
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value (1.0 / 3.0) * parameter ( (* (this->M_vectorsParameters) ) [0] ) *
                pow( J, -2.0/3.0) * I_C  * dot ( ( F_T * transpose (grad (phi_j) ) * F_T ), grad (phi_i) )
              ) >> jacobian;

    //    }

    jacobian->globalAssemble();
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::apply ( const vector_Type& sol, vector_Type& res,
                                                    const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                    const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                    const displayerPtr_Type displayer)
{
    computeStiffness (sol, 0., this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    res += *M_stiff;
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeStiffness ( const vector_Type&       disp,
                                                               Real                     /*factor*/,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                               const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                               const displayerPtr_Type& displayer )
{
    using namespace ExpressionAssembly;

    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n******************************************************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the Neo-Hookean nonlinear stiffness vector"     );
    displayer->leaderPrint (" \n******************************************************************\n  ");

    M_stiff.reset (new vector_Type (*this->M_localMap) );
    * (M_stiff) *= 0.0;

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
    //     pointerListOfIndexes.reset( new vectorIndexes_Type( (*mapsMarkerIndexes)[marker]) );

    //     Real mu     = dataMaterial->mu(marker);
    //     Real bulk   = dataMaterial->bulk(marker);

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

    //Computation of the volumetric part
    //
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow ( J , 2.0) - J + log (J) ) * dot (  F_T, grad (phi_i) )
              ) >> M_stiff;

    //Computation of the isochoric part
    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, -2.0 / 3.0) * (dot ( F - value (1.0 / 3.0) * I_C * F_T, grad (phi_i) ) )
              ) >> M_stiff;

    //    }
    this->M_stiff->globalAssemble();
}

template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::showMe ( std::string const& fileNameStiff,
                                                     std::string const& fileNameJacobian)
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}

template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
        const Epetra_SerialDenseMatrix& tensorF,
        const Epetra_SerialDenseMatrix& cofactorF,
        const std::vector<Real>& invariants,
        const UInt marker)
{

    //Get the material parameters
    Real mu       = this->M_dataMaterial->mu (marker);
    Real bulk     = this->M_dataMaterial->bulk (marker);

    //Computing the first term \muJ^{-2/3}[F-(1/3)tr(C)F^{-T}]
    Epetra_SerialDenseMatrix firstTerm (tensorF);
    Epetra_SerialDenseMatrix copyCofactorF (cofactorF);
    Real scale ( 0.0 );
    scale = -1 * (1.0 / 3.0) * invariants[0];
    copyCofactorF.Scale ( scale );
    firstTerm += copyCofactorF;

    Real coef ( 0.0 );
    coef = mu * std::pow (invariants[3], -2.0 / 3.0);
    firstTerm.Scale ( coef );

    //Computing the second term (volumetric part) J*(bulk/2)(J-1+(1/J)*ln(J))F^{-T}
    Epetra_SerialDenseMatrix secondTerm (cofactorF);
    Real sCoef (0);
    sCoef = invariants[3] * (bulk / 2.0) * (invariants[3] - 1 + (1 / invariants[3]) * std::log (invariants[3]) );
    secondTerm.Scale ( sCoef );

    firstPiola += firstTerm;
    firstPiola += secondTerm;
}


template <typename MeshType>
void NeoHookeanMaterialNonLinear<MeshType>::computeCauchyStressTensor ( const vectorPtr_Type disp,
									const QuadratureRule& evalQuad,
									vectorPtr_Type sigma_1,
									vectorPtr_Type sigma_2,
									vectorPtr_Type sigma_3)

{

  using namespace ExpressionAssembly;

  // Definition that are needed
  // Definition of F
  tensorF_Type F = ExpressionDefinitions::deformationGradient( this->M_dispETFESpace,  *disp, this->M_offset, this->M_identity );

  // Definition of J
  determinantF_Type J = ExpressionDefinitions::determinantF( F );

  //Definition of C
  tensorC_Type C = ExpressionDefinitions::tensorC( transpose(F), F );

  // Definition of F^-T
  minusT_Type  F_T = ExpressionDefinitions::minusT( F );

  // Definition I_C
  traceTensor_Type I_C = ExpressionDefinitions::traceTensor( C );

  evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
		evalQuad,
		this->M_dispETFESpace,
		meas_K *  dot ( vectorFromMatrix( ( 1 / J )*
						  ( value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow ( J , 2.0) - J + log (J) ) * F_T + // volumetric
						    parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, -2.0 / 3.0) * ( F - value (1.0 / 3.0) * I_C * F_T )
						    ) *
						  transpose(F), 0 ), phi_i)
		) >> sigma_1;
  sigma_1->globalAssemble();


  evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
		evalQuad,
		this->M_dispETFESpace,
		meas_K *  dot ( vectorFromMatrix( ( 1 / J )*
						  ( value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow ( J , 2.0) - J + log (J) ) * F_T + // volumetric
						    parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, -2.0 / 3.0) * ( F - value (1.0 / 3.0) * I_C * F_T ) ) *
						  transpose( F ) ,  1 ), phi_i)
		) >> sigma_2;
  sigma_2->globalAssemble();

  evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
		evalQuad,
		this->M_dispETFESpace,
		meas_K *  dot ( vectorFromMatrix( ( 1 / J )*
						  ( value (1.0 / 2.0) * parameter ( (* (this->M_vectorsParameters) ) [1] ) * ( pow ( J , 2.0) - J + log (J) ) * F_T + // volumetric
						    parameter ( (* (this->M_vectorsParameters) ) [0] ) * pow (J, -2.0 / 3.0) * ( F - value (1.0 / 3.0) * I_C * F_T )
						    ) *
						  transpose( F ) , 2 ), phi_i)
		) >> sigma_3;
  sigma_3->globalAssemble();

}


template <typename MeshType>
inline StructuralIsotropicConstitutiveLaw<MeshType>* createNeoHookeanMaterialNonLinear()
{
    return new NeoHookeanMaterialNonLinear<MeshType >();
}
namespace
{
static bool registerNH = StructuralIsotropicConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureIsotropicMaterialFactory::instance().registerProduct ( "neoHookean", &createNeoHookeanMaterialNonLinear<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __NEOHOOKENANMATERIAL_H */
