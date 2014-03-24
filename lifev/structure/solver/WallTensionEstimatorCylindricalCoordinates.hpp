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
 *  @brief This file contains a class which can be used to evaluate the wall tension in the arterial wall.
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 *
 */

#ifndef WALLTENSIONCYLINDRICAL_H
#define WALLTENSIONCYLINDRICAL_H 1

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <boost/scoped_ptr.hpp>
#include <boost/multi_array.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

//Trilinos include
#include <Epetra_SerialDenseMatrix.h>
#include <iostream>
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// LifeV core includes
#include <lifev/core/array/MatrixElemental.hpp>
#include <lifev/core/array/VectorElemental.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/fem/Assembly.hpp>
#include <lifev/core/fem/AssemblyElemental.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/LifeV.hpp>
#include <lifev/core/util/Displayer.hpp>
#include <lifev/core/array/MapEpetra.hpp>

#include <lifev/core/filter/ExporterEnsight.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#endif
#include <lifev/core/filter/ExporterEmpty.hpp>

// Structure module include
#include <lifev/structure/fem/AssemblyElementalStructure.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>

//Mother class
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>
#include <lifev/structure/solver/WallTensionEstimator.hpp>

namespace LifeV
{
/*!
  \class WallTensionEstimator
  \brief
  This class lets to compute the wall tensions inside the arterial wall using the cylindrical coordinates
  r, \theta, z. The tensorial operations
  that are needed to compute the stress tensor are defined in AssemblyElementalStructure. When a new
  type of analysis wants to be performed new methods can be added
*/

template <typename Mesh>
class WallTensionEstimatorCylindricalCoordinates : public WallTensionEstimator<Mesh>
{
public:

    //!@name Type definitions
    //@{

    // Data classes
    typedef WallTensionEstimator<Mesh>                    super;

    // Communicator
    typedef typename super::comm_Type                     comm_Type;
    typedef typename super::commPtr_Type                  commPtr_Type;

    // FE space
    typedef typename super::feSpace_Type                  feSpace_Type;
    typedef typename super::feSpacePtr_Type               feSpacePtr_Type;

    typedef typename super::feSpaceET_Type                feSpaceET_Type;
    typedef typename super::feSpaceETPtr_Type             feSpaceETPtr_Type;

    typedef StructuralConstitutiveLawData                 data_Type;
    typedef WallTensionEstimatorData                      analysisData_Type;
    typedef typename boost::shared_ptr<data_Type>         dataPtr_Type;
    typedef typename boost::shared_ptr<analysisData_Type> analysisDataPtr_Type;

    //Matrices 3x3 and std::vector for the invariants
    typedef Epetra_SerialDenseMatrix                      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                matrixPtr_Type;
    typedef std::vector<LifeV::Real>                      vector_Type;
    typedef boost::shared_ptr<vector_Type>                vectorPtr_Type;

    // These two are to handle the vector displacement read from hdf5
    typedef VectorEpetra                                  solutionVect_Type;
    typedef boost::shared_ptr<VectorEpetra>               solutionVectPtr_Type;

    // Displayer and Exporter classes
    typedef typename boost::shared_ptr<const Displayer>   displayerPtr_Type;
    typedef typename boost::shared_ptr< Exporter<Mesh> >  exporterPtr_Type;

    // Materials
    typedef StructuralConstitutiveLaw<Mesh>               material_Type;
    typedef boost::shared_ptr<material_Type>              materialPtr_Type;

    //@}


    //! @name Constructor &  Deconstructor
    //@{

    WallTensionEstimatorCylindricalCoordinates();

    virtual ~WallTensionEstimatorCylindricalCoordinates() {};

    //@}


    //!@name Methods
    //@{

    //! Setup the created object of the class WallTensionEstimatorCylindricalCoordinates
    /*!
      \param dataMaterial: the class containing the VenantKirchhoffElasticData
      \param tensionData: the class containing the WallTensionEstimatorData
      \param dFESpace: the FiniteElement Space
      \param displayer: the displayer object
    */
    void setup ( const dataPtr_Type& dataMaterial,
                 const analysisDataPtr_Type& tensionData,
                 const feSpacePtr_Type& FESpace,
                 const feSpaceETPtr_Type& ETFESpace,
                 const commPtr_Type& comm,
                 UInt marker);



    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE FESpace<Mesh, MapEpetra>& copyFESpace
    */
    void analyzeTensions();

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryDisplacementCylindrical();

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryEigenvaluesCylindrical();

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryCauchyStressesCylindrical();

    //@}


    //! @name Set Methods
    //@{

    //@}


    //! @name Get Methods
    //@{

    //@}

protected:

    //! @name Protected methods
    //@{
    //! moveToCylindricalCoordinates: This methods brings the gradient of the displacement field computed with respect to
    //! the classical x,y,z reference system to the cylindrical one r, \theta, zeta
    /*!
      \param deformationF: the gradient of the displacement with respect to the normal coordinates x,y,z
      \param deformationCylindricalF: the tensor F with respect to the second reference system
    */

    void moveToCylindricalCoordinates (matrix_Type& gradientDispl,
                                       UInt iloc,
                                       matrix_Type& deformationCylindricalF);

    //! constructGlobalStressVector: This method construct the vectors \sigma_{.,i} for i=x,y,z to have for each DOF the tensor \sigma
    /*!
      \param NONE
    */
    void constructGlobalStressVector();

    //@}


    //! @name Protected members
    //@{

    //! Elementary matrix for the tensor F
    matrixPtr_Type M_deformationCylindricalF;

    //Construct the matrix of change of variables given the position of the DOF
    matrixPtr_Type M_changeOfVariableMatrix;

    //@}
};

//=====================================
// Constructor
//=====================================

template <typename Mesh>
WallTensionEstimatorCylindricalCoordinates<Mesh>::WallTensionEstimatorCylindricalCoordinates( ) :
    super( ),
    M_deformationCylindricalF( ),
    M_changeOfVariableMatrix( )
{

}



//====================================
// Public Methods
//===================================
template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::setup ( const dataPtr_Type& dataMaterial,
                                                           const analysisDataPtr_Type& tensionData,
                                                           const feSpacePtr_Type& FESpace,
                                                           const feSpaceETPtr_Type& ETFESpace,
                                                           const commPtr_Type& comm,
                                                           UInt marker )
{
    super::setup (dataMaterial, tensionData, FESpace, ETFESpace, comm, marker);
    M_deformationCylindricalF.reset  ( new matrix_Type ( this->M_FESpace->fieldDim(), this->M_FESpace->fieldDim() ) );
    M_changeOfVariableMatrix.reset   ( new matrix_Type (this->M_FESpace->fieldDim(), this->M_FESpace->fieldDim() ) );
}

template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensions ( void )
{
    //Initialize the global vector to zero
    * (this->M_globalEigenvalues) *= 0.0;
    if ( !this->M_analysisData->recoveryVariable().compare ("displacement") )
    {
        analyzeTensionsRecoveryDisplacementCylindrical();
    }
    else if ( !this->M_analysisData->recoveryVariable().compare ("eigenvalues") )
    {
        analyzeTensionsRecoveryEigenvaluesCylindrical();
    }
    else
    {
        analyzeTensionsRecoveryCauchyStressesCylindrical();
    }
}

template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryDisplacementCylindrical ( void )
{

    LifeChrono chrono;

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the displacement..., ", this->M_dataMaterial->solidType() );
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    solutionVectPtr_Type grDisplX ( new solutionVect_Type (* (this->M_FESpace->mapPtr() ) ) );
    solutionVectPtr_Type grDisplY ( new solutionVect_Type (* (this->M_FESpace->mapPtr() ) ) );
    solutionVectPtr_Type grDisplZ ( new solutionVect_Type (* (this->M_FESpace->mapPtr() ) ) );

    //Compute the deformation gradient tensor F of the displ field
    super::computeDisplacementGradient ( grDisplX, grDisplY, grDisplZ);

    solutionVect_Type grXRep (*grDisplX, Repeated);
    solutionVect_Type grYRep (*grDisplY, Repeated);
    solutionVect_Type grZRep (*grDisplZ, Repeated);
    solutionVect_Type dRep (* (this->M_displacement), Repeated);

    //For each of the DOF, the Cauchy tensor is computed.
    //Therefore the tensor C,P, \sigma are computed for each DOF

    chrono.start();

    for ( UInt i (0); i < this->M_FESpace->mesh()->numVolumes(); ++i )
    {

        //Setup quantities
        this->M_FESpace->fe().updateFirstDerivQuadPt ( this->M_FESpace->mesh()->volumeList ( i ) );
        UInt eleID = this->M_FESpace->fe().currentLocalId();
        this->M_marker = this->M_FESpace->mesh()->volumeList ( i ).markerID();

        //store the local \grad(u)
        matrix_Type gradientDispl ( this->M_FESpace->fieldDim(), this->M_FESpace->fieldDim() );
        gradientDispl.Scale ( 0.0 );

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * this->M_FESpace->dim() + this->M_offset;

                gradientDispl (iComp, 0) = grXRep[ig];
                gradientDispl (iComp, 1) = grYRep[ig];
                gradientDispl (iComp, 2) = grZRep[ig];
            }

            //Reinitialization of matrices and arrays
            ( * (this->M_cofactorF) ).Scale (0.0);
            ( * (this->M_firstPiola) ).Scale (0.0);
            ( * (this->M_sigma) ).Scale (0.0);

            //Moving to cylindrical coordinates
            moveToCylindricalCoordinates ( gradientDispl, iloc, *M_deformationCylindricalF );

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (this->M_invariants, *M_deformationCylindricalF, * (this->M_cofactorF) );

            //Compute the first Piola-Kirchhoff tensor
            this->M_material->computeLocalFirstPiolaKirchhoffTensor (* (this->M_firstPiola), *M_deformationCylindricalF, * (this->M_cofactorF), this->M_invariants, this->M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (* (this->M_sigma), * (this->M_firstPiola), this->M_invariants[3], *M_deformationCylindricalF);

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (* (this->M_sigma), this->M_eigenvaluesR, this->M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( UInt j (0); j < this->M_eigenvaluesI.size(); ++j )
            {
                sum += std::abs ( this->M_eigenvaluesI[j] );
            }

            ASSERT ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( this->M_eigenvaluesR.begin(), this->M_eigenvaluesR.end() );

            std::cout << "Saving eigenvalues" << i << std::endl;

            //Save the eigenvalues in the global vector
            for ( UInt icoor = 0; icoor < this->M_FESpace->fieldDim(); ++icoor )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + icoor * this->M_FESpace->dim() + this->M_offset;
                (* (this->M_globalEigenvalues) ) (ig) = this->M_eigenvaluesR[icoor];
                // Int LIDid = this->M_displacement->blockMap().LID(ig);
                // if( this->M_globalEigenvalues->blockMap().LID(ig) != -1  )
                // {
                //     Int GIDid = this->M_globalEigenvalues->blockMap().GID(LIDid);

                // }
            }
        }
    }

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );

}


template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryEigenvaluesCylindrical ( void )
{

    LifeChrono chrono;

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the tensions..., ", this->M_dataMaterial->solidType() );
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    solutionVect_Type patchArea (* (this->M_displacement), Unique, Add);
    patchArea *= 0.0;

    super::constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < this->M_FESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += this->M_FESpace->qr().weight (iq);
    }

    Real wQuad (refElemArea / this->M_FESpace->refFE().nbDof() );

    //Setting the quadrature Points = DOFs of the element and weight = 1
    std::vector<GeoVector> coords = this->M_FESpace->refFE().refCoor();
    std::vector<Real> weights (this->M_FESpace->fe().nbFEDof(), wQuad);
    fakeQuadratureRule.setDimensionShape ( shapeDimension (this->M_FESpace->refFE().shape() ), this->M_FESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //Set the new quadrature rule
    this->M_FESpace->setQuadRule (fakeQuadratureRule);

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (this->M_FESpace->fe().nbFEDof(), * (this->M_deformationF) );
    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (* (this->M_displacement), Repeated);

    VectorElemental elVecTens (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    chrono.start();

    //Loop on each volume
    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); ++i )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt ( this->M_FESpace->mesh()->volumeList ( i ) );
        elVecTens.zero();

        this->M_marker = this->M_FESpace->mesh()->volumeList ( i ).markerID();

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * this->M_FESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * this->M_FESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradientWithoutIdentity ( dk_loc, vectorDeformationF, this->M_FESpace->fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF = 0; nDOF < ( UInt ) this->M_FESpace->fe().nbFEDof(); nDOF++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( nDOF );
            vector_Type localDisplacement (this->M_FESpace->fieldDim(), 0.0);

            for ( UInt coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                localDisplacement[coor] = iloc + coor * this->M_FESpace->fe().nbFEDof();
            }

            this->M_sigma->Scale (0.0);
            this->M_firstPiola->Scale (0.0);
            this->M_cofactorF->Scale (0.0);
            M_deformationCylindricalF->Scale (0.0);

            moveToCylindricalCoordinates (vectorDeformationF[nDOF], iloc, *M_deformationCylindricalF);

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (this->M_invariants, *M_deformationCylindricalF, * (this->M_cofactorF) );

            //Compute the first Piola-Kirchhoff tensor
            this->M_material->computeLocalFirstPiolaKirchhoffTensor (* (this->M_firstPiola), *M_deformationCylindricalF, * (this->M_cofactorF), this->M_invariants, this->M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (* (this->M_sigma), * (this->M_firstPiola), this->M_invariants[3], *M_deformationCylindricalF);

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (* (this->M_sigma), this->M_eigenvaluesR, this->M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( int i = 0; i < this->M_eigenvaluesI.size(); i++ )
            {
                sum += std::abs (this->M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( this->M_eigenvaluesR.begin(), this->M_eigenvaluesR.end() );

            //Assembling the local vector
            for ( int coor = 0; coor < this->M_eigenvaluesR.size(); coor++ )
            {
                elVecTens[iloc + coor * this->M_FESpace->fe().nbFEDof()] = this->M_eigenvaluesR[coor];
            }
        }

        super::reconstructElementaryVector ( elVecTens, patchAreaR, *this->M_FESpace );

        //Assembling the local into global vector
        for ( UInt ic = 0; ic < this->M_FESpace->fieldDim(); ++ic )
        {
            assembleVector (* (this->M_globalEigenvalues), elVecTens, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }
    }

    this->M_globalEigenvalues->globalAssemble();

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );
}

template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryCauchyStressesCylindrical ( void )
{

    LifeChrono chrono;

    chrono.start();
    UInt dim = this->M_FESpace->dim();

    constructGlobalStressVector();


    for ( UInt iDOF = 0; iDOF < ( UInt ) this->M_FESpace->dof().numTotalDof(); iDOF++ )
    {

        if ( this->M_displacement->blockMap().LID (iDOF) != -1 ) // The Global ID is on the calling processors
        {

            (* (this->M_sigma) ).Scale (0.0);

            //Extracting the gradient of U on the current DOF
            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                Int LIDid = this->M_displacement->blockMap().LID (iDOF + iComp * dim + this->M_offset);
                Int GIDid = this->M_displacement->blockMap().GID (LIDid);
                (* (this->M_sigma) ) (iComp, 0) = (*this->M_sigmaX) (GIDid); // (d_xX,d_yX,d_zX)
                (* (this->M_sigma) ) (iComp, 1) = (*this->M_sigmaY) (GIDid); // (d_xY,d_yY,d_zY)
                (* (this->M_sigma) ) (iComp, 2) = (*this->M_sigmaZ) (GIDid); // (d_xZ,d_yZ,d_zZ)
            }

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (* (this->M_sigma), this->M_eigenvaluesR, this->M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( int i = 0; i < this->M_eigenvaluesI.size(); i++ )
            {
                sum += std::abs (this->M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( this->M_eigenvaluesR.begin(), this->M_eigenvaluesR.end() );

            //Save the eigenvalues in the global vector
            for ( UInt icoor = 0; icoor < this->M_FESpace->fieldDim(); ++icoor )
            {
                Int LIDid = this->M_displacement->blockMap().LID (iDOF + icoor * dim + this->M_offset);
                Int GIDid = this->M_displacement->blockMap().GID (LIDid);
                (* (this->M_globalEigenvalues) ) (GIDid) = this->M_eigenvaluesR[icoor];
            }

        }
    }

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );

}

template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::constructGlobalStressVector()
{

    //Creating the local stress tensors
    VectorElemental elVecSigmaX (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );
    VectorElemental elVecSigmaY (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );
    VectorElemental elVecSigmaZ (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    LifeChrono chrono;

    //Constructing the patch area vector for reconstruction purposes
    solutionVect_Type patchArea (* (this->M_displacement), Unique, Add);
    patchArea *= 0.0;

    super::constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < this->M_FESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += this->M_FESpace->qr().weight (iq);
    }

    Real wQuad (refElemArea / this->M_FESpace->refFE().nbDof() );

    //Setting the quadrature Points = DOFs of the element and weight = 1
    std::vector<GeoVector> coords = this->M_FESpace->refFE().refCoor();
    std::vector<Real> weights (this->M_FESpace->fe().nbFEDof(), wQuad);
    fakeQuadratureRule.setDimensionShape ( shapeDimension (this->M_FESpace->refFE().shape() ), this->M_FESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //Set the new quadrature rule
    this->M_FESpace->setQuadRule (fakeQuadratureRule);

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the Cauchy stresses..., ", this->M_dataMaterial->solidType() );
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (this->M_FESpace->fe().nbFEDof(), * (this->M_deformationF) );
    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (* (this->M_displacement), Repeated);

    chrono.start();

    //Loop on each volume
    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); ++i )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt ( this->M_FESpace->mesh()->volumeList ( i ) );

        elVecSigmaX.zero();
        elVecSigmaY.zero();
        elVecSigmaZ.zero();

        this->M_marker = this->M_FESpace->mesh()->volumeList ( i ).markerID();

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * this->M_FESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * this->M_FESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradientWithoutIdentity ( dk_loc, vectorDeformationF, this->M_FESpace->fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF = 0; nDOF < ( UInt ) this->M_FESpace->fe().nbFEDof(); nDOF++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst ( nDOF );

            vector_Type localDisplacement (this->M_FESpace->fieldDim(), 0.0);

            for ( UInt coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                localDisplacement[coor] = iloc + coor * this->M_FESpace->fe().nbFEDof();
            }

            this->M_sigma->Scale (0.0);
            this->M_firstPiola->Scale (0.0);
            this->M_cofactorF->Scale (0.0);
            this->M_deformationCylindricalF->Scale (0.0);

            moveToCylindricalCoordinates (vectorDeformationF[nDOF], iloc, *M_deformationCylindricalF);

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (this->M_invariants, *M_deformationCylindricalF, * (this->M_cofactorF) );

            //Compute the first Piola-Kirchhoff tensor
            this->M_material->computeLocalFirstPiolaKirchhoffTensor (* (this->M_firstPiola), *M_deformationCylindricalF, * (this->M_cofactorF), this->M_invariants, this->M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (* (this->M_sigma), * (this->M_firstPiola), this->M_invariants[3], *M_deformationCylindricalF);

            //Assembling the local vectors for local tensions Component X
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaX) [iloc + coor * this->M_FESpace->fe().nbFEDof()] = (* (this->M_sigma) ) (coor, 0);
            }

            //Assembling the local vectors for local tensions Component Y
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaY) [iloc + coor * this->M_FESpace->fe().nbFEDof()] = (* (this->M_sigma) ) (coor, 1);
            }

            //Assembling the local vectors for local tensions Component Z
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaZ) [iloc + coor * this->M_FESpace->fe().nbFEDof()] = (* (this->M_sigma) ) (coor, 2);
            }

        }

        super::reconstructElementaryVector ( elVecSigmaX, patchAreaR, *this->M_FESpace );
        super::reconstructElementaryVector ( elVecSigmaY, patchAreaR, *this->M_FESpace );
        super::reconstructElementaryVector ( elVecSigmaZ, patchAreaR, *this->M_FESpace );

        //Assembling the three elemental vector in the three global
        for ( UInt ic = 0; ic < this->M_FESpace->fieldDim(); ++ic )
        {
            assembleVector (*this->M_sigmaX, elVecSigmaX, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
            assembleVector (*this->M_sigmaY, elVecSigmaY, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
            assembleVector (*this->M_sigmaZ, elVecSigmaZ, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }
    }


    this->M_sigmaX->globalAssemble();
    this->M_sigmaY->globalAssemble();
    this->M_sigmaZ->globalAssemble();
}

template <typename Mesh>
void
WallTensionEstimatorCylindricalCoordinates<Mesh >::moveToCylindricalCoordinates ( matrix_Type& gradientDispl,
        UInt iloc,
        matrix_Type& deformationCylindricalF)
{
    //Reinitialize
    (*M_changeOfVariableMatrix).Scale ( 0.0 );
    (*M_deformationCylindricalF).Scale ( 0.0 );

    //Extracting the coordinates of the iloc DOF on the reference element
    Real xi ( this->M_FESpace->refFE().xi (iloc) );
    Real eta ( this->M_FESpace->refFE().eta (iloc) );
    Real zeta ( this->M_FESpace->refFE().zeta (iloc) );
    //Trasforming them in the local coordinates
    Real x (0);
    Real y (0);
    Real z (0);
    this->M_FESpace->fe().coorMap (x, y, z, xi, eta, zeta);

    // Defining the new variables
    //Real radius= std::sqrt( x*x + y*y  );
    Real theta = std::atan ( y / x );

    //Filling the change of variable Matrix and its derivative with respect to theta
    (*M_changeOfVariableMatrix) (0, 0) =   std::cos (theta);
    (*M_changeOfVariableMatrix) (1, 0) = - std::sin (theta);
    (*M_changeOfVariableMatrix) (2, 0) =   0.0;
    (*M_changeOfVariableMatrix) (0, 1) =   std::sin (theta);
    (*M_changeOfVariableMatrix) (1, 1) =   std::cos (theta);
    (*M_changeOfVariableMatrix) (2, 1) =   0.0;
    (*M_changeOfVariableMatrix) (0, 2) =   0.0;
    (*M_changeOfVariableMatrix) (1, 2) =   0.0;
    (*M_changeOfVariableMatrix) (2, 2) =   1.0;

    //Move to cylindrical gradient at the DOF
    deformationCylindricalF.Multiply ('N', 'T', 1.0, gradientDispl, *M_changeOfVariableMatrix, 0.0 );

    //Add the identity
    for ( Int icoor (0); icoor < this->M_FESpace->fieldDim(); icoor++ )
    {
        deformationCylindricalF ( icoor, icoor ) += 1.0;
    }

}

}

#endif /*WALLTENSIONCYLINDRICAL_H 1*/
