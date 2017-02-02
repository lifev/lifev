//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
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
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _LINVENANTKIRCHHOFFMATERIAL_H_
#define _LINVENANTKIRCHHOFFMATERIAL_H_

#include <lifev/structure/solver/isotropic/StructuralIsotropicConstitutiveLaw.hpp>


namespace LifeV
{
template <typename MeshType>
class VenantKirchhoffMaterialLinear :
    public StructuralIsotropicConstitutiveLaw<MeshType>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralIsotropicConstitutiveLaw<MeshType>                 super;

    typedef typename super::data_Type                         data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;

    typedef typename super::vectorVolumes_Type       vectorVolumes_Type;
    typedef std::shared_ptr<vectorVolumes_Type>    vectorVolumesPtr_Type;

    typedef std::vector<UInt>                             vectorIndexes_Type;
    typedef std::shared_ptr<vectorIndexes_Type>    vectorIndexesPtr_Type;
    typedef std::map< UInt, vectorIndexes_Type>           mapMarkerIndexes_Type;
    typedef std::shared_ptr<mapMarkerIndexes_Type>      mapMarkerIndexesPtr_Type;
    typedef typename mapMarkerIndexes_Type::const_iterator mapIteratorIndex_Type;

    typedef typename super::FESpacePtr_Type          FESpacePtr_Type;
    typedef typename super::ETFESpacePtr_Type        ETFESpacePtr_Type;

    //Vector for vector parameters
    typedef typename super::vectorsParameters_Type       vectorsParameters_Type;
    typedef typename super::vectorsParametersPtr_Type    vectorsParametersPtr_Type;

    typedef MatrixSmall<3, 3>                         matrixSmall_Type;

    typedef typename super::tensorF_Type               tensorF_Type;

    //@}

    //! @name Constructor &  Destructor
    //@{

    VenantKirchhoffMaterialLinear();

    virtual  ~VenantKirchhoffMaterialLinear();

    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralIsotropicConstitutiveLaw
    /*!
    \param dFespace: the FiniteElement Space
    \param monolithicMap: the MapEpetra
    \param offset: the offset parameter used assembling the matrices
    */
    void setup (const FESpacePtr_Type& dFESpace,
                const ETFESpacePtr_Type& ETFESpace,
                const std::shared_ptr<const MapEpetra>&  monolithicMap,
                const UInt offset,const dataPtr_Type& dataMaterial);

    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/,
                              const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                              const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/  );

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
    void updateNonLinearJacobianTerms (  matrixPtr_Type& jacobian,
                                         const vector_Type& /*disp*/,
                                         const dataPtr_Type& /*dataMaterial*/,
                                         const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                         const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                         const displayerPtr_Type& /*displayer*/);

    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the
        material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */

    void computeStiffness ( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                            const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                            const displayerPtr_Type& displayer );

    void computeKinematicsVariables ( const VectorElemental& /*dk_loc*/ ) {}

    //! ShowMe method of the class (saved on a file the two matrices)
    void showMe ( std::string const& fileNameStiff,
                  std::string const& fileNameJacobian);

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

    //! Get the linear part of the matrix
    matrixPtr_Type const linearStiff() const
    {
        return M_linearStiff;
    }

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const
    {
        return M_stiff;
    }

    //! Get the Stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        vectorPtr_Type zero ( new vector_Type() );
        return zero;
    }

    void apply ( const vector_Type& sol, vector_Type& res,
                 const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                 const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                 const displayerPtr_Type /*displayer*/)
    {
        res += *M_stiff * sol;
    }

    //@}

protected:

    //! construct the vectors for the parameters
    /*!
      \param VOID
      \return VOID
    */
    void setupVectorsParameters ( void);

    //! Protected members

    //! Matrix Kl: stiffness linear
    matrixPtr_Type                                 M_linearStiff;

    //! Matrix Kl: stiffness linear
    matrixPtr_Type                                 M_stiff;

};

template <typename MeshType>
VenantKirchhoffMaterialLinear<MeshType>::VenantKirchhoffMaterialLinear() :
    super             ( ),
    M_linearStiff                ( ),
    M_stiff                      ( )
{
}

template <typename MeshType>
VenantKirchhoffMaterialLinear<MeshType>::~VenantKirchhoffMaterialLinear()
{}


template <typename MeshType>
void
VenantKirchhoffMaterialLinear<MeshType>::setup (const FESpacePtr_Type& dFESpace,
                                                const ETFESpacePtr_Type& dETFESpace,
                                                const std::shared_ptr<const MapEpetra>&  monolithicMap,
                                                const UInt offset,const dataPtr_Type& dataMaterial)
{
    this->M_dataMaterial                  = dataMaterial;
    this->M_dispFESpace                   = dFESpace;
    this->M_dispETFESpace                 = dETFESpace;
    this->M_localMap                      = monolithicMap;
    this->M_linearStiff.reset             (new matrix_Type (*this->M_localMap) );
    this->M_offset                        = offset;

    // The 2 is because the law uses two parameters.
    // another way would be to set up the number of constitutive parameters of the law
    // in the data file to get the right size. Note the comment below.
    this->M_vectorsParameters.reset ( new vectorsParameters_Type ( 2 ) );

    this->setupVectorsParameters( );

}


template <typename MeshType>
void
VenantKirchhoffMaterialLinear<MeshType>::setupVectorsParameters ( void )
{
    // Paolo Tricerri: February, 20th
    // In each class, the name of the parameters has to inserted in the law
    // TODO: move the saving of the material parameters to more abstract objects
    //       such that in the class of the material we do not need to call explicitly
    //       the name of the parameter.

    // Number of volume on the local part of the mesh
    UInt nbElements = this->M_dispFESpace->mesh()->numVolumes();

    // Parameter lambda
    // 1. resize the vector in the first element of the vector.
    (* (this->M_vectorsParameters) ) [0].resize ( nbElements );

    // Parameter mu
    (* (this->M_vectorsParameters) ) [1].resize ( nbElements );

    for (UInt i (0); i < nbElements; i++ )
    {
        // Extracting the marker
        UInt markerID = this->M_dispFESpace->mesh()->element ( i ).markerID();

        Real lambda = this->M_dataMaterial->lambda ( markerID );
        Real mu = this->M_dataMaterial->mu ( markerID );

        ( (* (this->M_vectorsParameters) ) [0]) [ i ] = lambda;
        ( (* (this->M_vectorsParameters) ) [1]) [ i ] = mu;
    }
}


template <typename MeshType>
void VenantKirchhoffMaterialLinear<MeshType>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/,
                                                                  const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                  const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/)
{
    using namespace ExpressionAssembly;

    * (this->M_linearStiff) *= 0.0;

    //Compute the linear part of the Stiffness Matrix.
    //In the case of Linear Material it is the Stiffness Matrix.
    //In the case of NonLinear Materials it must be added of the non linear part.

    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                parameter ( (* (this->M_vectorsParameters) ) [0]) * div ( phi_i ) * div ( phi_j )
              ) >> M_linearStiff;

    integrate ( elements ( this->M_dispETFESpace->mesh() ),
                this->M_dispFESpace->qr(),
                this->M_dispETFESpace,
                this->M_dispETFESpace,
                value ( 2.0 ) * parameter ( (* (this->M_vectorsParameters) ) [1]) * dot ( sym (grad (phi_j) ) , grad (phi_i) )
              ) >> M_linearStiff;

    this->M_linearStiff->globalAssemble();

    //Initialization of the pointer M_stiff to what is pointed by M_linearStiff
    this->M_stiff = this->M_linearStiff;
    this->M_jacobian = this->M_linearStiff;
}


template <typename MeshType>
void VenantKirchhoffMaterialLinear<MeshType>::updateJacobianMatrix (const vector_Type& disp,
                                                                    const dataPtr_Type& dataMaterial,
                                                                    const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                    const mapMarkerIndexesPtr_Type mapsMarkerIndexes,
                                                                    const displayerPtr_Type& displayer)
{
    //displayer->leaderPrint(" \n*********************************\n  ");
    displayer->leaderPrint ("  S-  Updating the Jacobian Matrix (constant, Linear Elastic)\n");
    //displayer->leaderPrint(" \n*********************************\n  ");

    //displayer->leaderPrint(" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, this->M_dataMaterial, mapsMarkerVolumes, mapsMarkerIndexes, displayer);
    //displayer->leaderPrint(" \n*********************************\n  ");

}

template <typename MeshType>
void VenantKirchhoffMaterialLinear<MeshType>::updateNonLinearJacobianTerms ( matrixPtr_Type& /*jacobian*/,
                                                                             const  vector_Type& /*disp*/,
                                                                             const dataPtr_Type& /*dataMaterial*/,
                                                                             const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                             const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                                                             const displayerPtr_Type& /*displayer*/ )
{
    //    this->M_stiff->globalAssemble();
    //  displayer->leaderPrint("  S- Doing nothing - Updating non linear terms in Jacobian Matrix (constant, Linear Elastic)\n");
}

template <typename MeshType>
void VenantKirchhoffMaterialLinear<MeshType>::computeStiffness ( const vector_Type& /*disp*/,
                                                                 Real /*factor*/,
                                                                 const dataPtr_Type& /*dataMaterial*/,
                                                                 const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                 const mapMarkerIndexesPtr_Type /*mapsMarkerIndexes*/,
                                                                 const displayerPtr_Type& displayer )

{
    displayer->leaderPrint (" \n*********************************\n  ");
    displayer->leaderPrint ("  S- Using the the Stiffness Matrix (constant, Linear Elastic)");
    displayer->leaderPrint (" \n*********************************\n  ");
}


template <typename MeshType>
void
VenantKirchhoffMaterialLinear<MeshType>::showMe ( std::string const& fileNameStiff,
                                                  std::string const& fileNameJacobian
                                                )
{
    //This string is to save the linear part
    std::string fileNamelinearStiff =  fileNameStiff;
    fileNamelinearStiff += "linear";

    this->M_linearStiff->spy (fileNamelinearStiff);
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}

template <typename MeshType>
void
VenantKirchhoffMaterialLinear<MeshType>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
										 const Epetra_SerialDenseMatrix& tensorF,
										 const Epetra_SerialDenseMatrix& cofactorF,
										 const std::vector<Real>& invariants,
										 const UInt marker)
{

    //Get the material parameters
    Real lambda   = this->M_dataMaterial->lambda (marker);
    Real mu       = this->M_dataMaterial->mu (marker);

    Epetra_SerialDenseMatrix copyF (tensorF);
    Epetra_SerialDenseMatrix identity (nDimensions, nDimensions);

    copyF (0,0) += -1.0;
    copyF (1,1) += -1.0;
    copyF (2,2) += -1.0;

    identity (0,0) = 1.0;
    identity (1,1) = 1.0;
    identity (2,2) = 1.0;

    Real divergenceU ( copyF (0, 0) + copyF (1, 1) + copyF (2, 2) ); //DivU = tr(copyF)
    Real coefIdentity (0.0);
    coefIdentity = divergenceU * lambda;
    identity.Scale ( coefIdentity );

    Epetra_SerialDenseMatrix transposed (copyF);
    transposed.SetUseTranspose (true);

    Epetra_SerialDenseMatrix secondTerm (nDimensions, nDimensions);

    secondTerm = copyF;
    secondTerm += transposed;
    secondTerm.Scale (mu);

    firstPiola = identity;
    firstPiola += secondTerm;

    // Here we multiply the Piola tensor for J and F^-T because in the WallTensionEstimator class
    // the transformation from the Piola to the Cauchy stress tensor is always done while for this model
    // the two tensors coincide
    firstPiola.Scale( invariants[3] );

    Epetra_SerialDenseMatrix tensorP (nDimensions, nDimensions);

    tensorP.Multiply ('N', 'N', 1.0, firstPiola, cofactorF, 0.0); //see Epetra_SerialDenseMatrix

    firstPiola.Scale(0.0);
    firstPiola += tensorP;
}

template <typename MeshType>
void VenantKirchhoffMaterialLinear<MeshType>::computeCauchyStressTensor ( const vectorPtr_Type disp,
                                                                          const QuadratureRule& evalQuad,
                                                                          vectorPtr_Type sigma_1,
                                                                          vectorPtr_Type sigma_2,
                                                                          vectorPtr_Type sigma_3)

{

    using namespace ExpressionAssembly;

    MatrixSmall<3, 3> identity;

    // Defining the div of u
    identity (0, 0) = 1.0; identity (0, 1) = 0.0; identity (0, 2) = 0.0;
    identity (1, 0) = 0.0; identity (1, 1) = 1.0; identity (1, 2) = 0.0;
    identity (2, 0) = 0.0; identity (2, 1) = 0.0; identity (2, 2) = 1.0;

    ExpressionTrace<ExpressionInterpolateGradient<MeshType, MapEpetra, 3, 3> > Div( grad( this->M_dispETFESpace,  *disp, this->M_offset ) );

    evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
                  evalQuad,
                  this->M_dispETFESpace,
                  meas_K *  dot ( vectorFromMatrix(  parameter ( (* (this->M_vectorsParameters) ) [0]) * Div * identity +
                                                     parameter ( (* (this->M_vectorsParameters) ) [1]) *
                                                     ( grad ( this->M_dispETFESpace,  *disp, this->M_offset )  +
                                                       transpose ( grad ( this->M_dispETFESpace,  *disp, this->M_offset ) )
                                                       )
                                                     ,  0 ), phi_i)
                  ) >> sigma_1;
    sigma_1->globalAssemble();



    evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
                  evalQuad,
                  this->M_dispETFESpace,
                  meas_K *  dot ( vectorFromMatrix(  parameter ( (* (this->M_vectorsParameters) ) [0]) * Div * identity +
                                                     parameter ( (* (this->M_vectorsParameters) ) [1]) *
                                                     ( grad ( this->M_dispETFESpace,  *disp, this->M_offset )  +
                                                       transpose ( grad ( this->M_dispETFESpace,  *disp, this->M_offset ) )
                                                       )
                                                     ,  1 ), phi_i)
                  ) >> sigma_2;
    sigma_2->globalAssemble();

    evaluateNode( elements ( this->M_dispETFESpace->mesh() ),
                  evalQuad,
                  this->M_dispETFESpace,
                  meas_K *  dot ( vectorFromMatrix(  parameter ( (* (this->M_vectorsParameters) ) [0]) * Div * identity  +
                                                     parameter ( (* (this->M_vectorsParameters) ) [1]) *
                                                     ( grad ( this->M_dispETFESpace,  *disp, this->M_offset )  +
                                                       transpose ( grad ( this->M_dispETFESpace,  *disp, this->M_offset ) )
                                                      )
                                                     ,  2 ), phi_i)
                  ) >> sigma_3;
    sigma_3->globalAssemble();

}


template <typename MeshType>
inline StructuralIsotropicConstitutiveLaw<MeshType>* createVenantKirchhoffLinear()
{
    return new VenantKirchhoffMaterialLinear<MeshType >();
}
namespace
{
static bool registerVKL = StructuralIsotropicConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureIsotropicMaterialFactory::instance().registerProduct ( "linearVenantKirchhoff", &createVenantKirchhoffLinear<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __LINVENANTKIRCHHOFF_H */
