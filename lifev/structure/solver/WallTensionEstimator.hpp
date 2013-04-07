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
 */

#ifndef WALLTENSION_H_
#define WALLTENSION_H_ 1

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

// STL classes
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>

// Boost classes
#include <boost/scoped_ptr.hpp>
#include <boost/multi_array.hpp>

// Trilinos classes
#include <Epetra_SerialDenseMatrix.h>

// Tell the compiler to restore the warning previously silented
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

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

#include <lifev/eta/fem/ETFESpace.hpp>

// Structure module include
#include <lifev/structure/fem/AssemblyElementalStructure.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLawData.hpp>
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

//Materials
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinearPenalized.hpp>
#include <lifev/structure/solver/isotropic/SecondOrderExponentialMaterialNonLinear.hpp>

namespace LifeV
{
/*!
  \class WallTensionEstimator
  \brief
  This class lets to compute the wall tensions inside the arterial wall. The tensorial operations
  that are needed to compute the stress tensor are defined in AssemblyElementalStructure. When a new
  type of analysis wants to be performed new methods can be added
*/

template <typename Mesh>
class WallTensionEstimator
{
public:

    //!@name Type definitions
    //@{

    // Communicator
    typedef Displayer::comm_Type                          comm_Type;
    typedef Displayer::commPtr_Type                       commPtr_Type;

    // FE space
    typedef FESpace < Mesh, MapEpetra >                   feSpace_Type;
    typedef boost::shared_ptr < feSpace_Type >            feSpacePtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 > feSpaceET_Type;
    typedef boost::shared_ptr<feSpaceET_Type>                     feSpaceETPtr_Type;

    // Data classes
    typedef StructuralConstitutiveLawData                 data_Type;
    typedef typename boost::shared_ptr<data_Type>         dataPtr_Type;
    typedef WallTensionEstimatorData                      analysisData_Type;
    typedef typename boost::shared_ptr<analysisData_Type> analysisDataPtr_Type;

    //Matrices 3x3 and std::vector for the invariants
    typedef Epetra_SerialDenseMatrix                      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>                matrixPtr_Type;
    typedef std::vector< Real >                           vector_Type;
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

    WallTensionEstimator();

    virtual ~WallTensionEstimator() {};

    //@}


    //!@name Methods
    //@{

    //! Setup the created object of the class WallTensionEstimator
    /*!
      \param dataMaterial: the class containing the VenantKirchhoffElasticData
      \param tensionData: the class containing the WallTensionEstimatorData
      \param dFESpace: the FiniteElement Space
      \param displayer: the displayer object
      \param marker volume marker
    */
    void setup ( const dataPtr_Type& dataMaterial,
                 const analysisDataPtr_Type& tensionData,
                 const feSpacePtr_Type& feSpace,
                 const feSpaceETPtr_Type& feSpaceET,
                 const commPtr_Type& comm,
                 UInt marker);

    //! Setup the created object of the class WallTensionEstimator
    /*!
      \param dataMaterial: the class containing the VenantKirchhoffElasticData
      \param dFESpace: the FiniteElement Space
      \param displayer: the displayer object
      \param marker volume marker
    */
    void setup ( const dataPtr_Type& dataMaterial,
                 const feSpacePtr_Type& feSpace,
                 const feSpaceETPtr_Type& feSpaceET,
                 const commPtr_Type& comm,
                 UInt marker);

    //! This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    virtual void analyzeTensions();

    //! This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    void analyzeTensionsRecoveryDisplacement();

    //! This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    void analyzeTensionsRecoveryEigenvalues();

    //! This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    void analyzeTensionsRecoveryCauchyStresses();

    //! This method computes the Von Mises stress. It uses the displacement vector that has to be set
    void analyzeTensionsRecoveryVonMisesStress ();

    //@}


    //! @name Set Methods
    //@{

    //! Set the displacement vector
    void setDisplacement ( const solutionVect_Type& displacement )
    {
        *M_displacement = displacement;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Get the Epetramap
    MapEpetra const& map() const
    {
        return *M_FESpace->mapPtr();
    }

    //! Get the FESpace object
    LIFEV_DEPRECATED ( const feSpace_Type& dFESpace() const
    {
        return feSpace();
    }
                     )

    //! Get the pointer to the FESpace object
    LIFEV_DEPRECATED ( const feSpacePtr_Type& dFESpacePtr() const
    {
        return feSpacePtr();
    }
                     )

    //! Get the FESpace object
    const feSpace_Type& feSpace() const
    {
        return *M_FESpace;
    }

    //! Get the pointer to the FESpace object
    const feSpacePtr_Type& feSpacePtr() const
    {
        return M_FESpace;
    }

    //! Get the scalar FESpace object
    const feSpace_Type& feSpaceScalar() const
    {
        return *M_FESpaceScalar;
    }

    //! Get the pointer to the scalar FESpace object
    const feSpacePtr_Type& feSpaceScalarPtr() const
    {
        return M_FESpaceScalar;
    }

    //! Get the displacement solution
    solutionVect_Type displacement() const
    {
        return *M_displacement;
    }

    //! Get the displacement solution
    solutionVect_Type gradientX() const
    {
        return *M_gradientX;
    }

    solutionVect_Type gradientY() const
    {
        return *M_gradientY;
    }

    solutionVect_Type gradientZ() const
    {
        return *M_gradientZ;
    }

    solutionVect_Type sigmaX() const
    {
        return *M_sigmaX;
    }

    solutionVect_Type sigmaY() const
    {
        return *M_sigmaY;
    }

    solutionVect_Type sigmaZ() const
    {
        return *M_sigmaZ;
    }

    //! Export the XX component of the stress by copying it to an external vector
    /*!
     * @param stressXX vector to be filled with the XX component of the stress
     */
    void exportStressXX ( solutionVect_Type& stressXX )
    {
        stressXX.subset ( *M_sigmaX, static_cast<UInt> ( 0 ) );
    }

    //! Export the XY component of the stress by copying it to an external vector
    /*!
     * @param stressXY vector to be filled with the XY component of the stress
     */
    void exportStressXY ( solutionVect_Type& stressXY )
    {
        stressXY.subset ( *M_sigmaX, static_cast<UInt> ( 1 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the XZ component of the stress by copying it to an external vector
    /*!
     * @param stressXZ vector to be filled with the XZ component of the stress
     */
    void exportStressXZ ( solutionVect_Type& stressXZ )
    {
        stressXZ.subset ( *M_sigmaX, static_cast<UInt> ( 2 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the YX component of the stress by copying it to an external vector
    /*!
     * @param stressYX vector to be filled with the YX component of the stress
     */
    void exportStressYX ( solutionVect_Type& stressYX )
    {
        stressYX.subset ( *M_sigmaY, static_cast<UInt> ( 0 ) );
    }

    //! Export the YY component of the stress by copying it to an external vector
    /*!
     * @param stressYY vector to be filled with the YY component of the stress
     */
    void exportStressYY ( solutionVect_Type& stressYY )
    {
        stressYY.subset ( *M_sigmaY, static_cast<UInt> ( 1 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the YZ component of the stress by copying it to an external vector
    /*!
     * @param stressYZ vector to be filled with the YZ component of the stress
     */
    void exportStressYZ ( solutionVect_Type& stressYZ )
    {
        stressYZ.subset ( *M_sigmaY, static_cast<UInt> ( 2 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the ZX component of the stress by copying it to an external vector
    /*!
     * @param stressZX vector to be filled with the ZX component of the stress
     */
    void exportStressZX ( solutionVect_Type& stressZX )
    {
        stressZX.subset ( *M_sigmaZ, static_cast<UInt> ( 0 ) );
    }

    //! Export the ZY component of the stress by copying it to an external vector
    /*!
     * @param stressZY vector to be filled with the ZY component of the stress
     */
    void exportStressZY ( solutionVect_Type& stressZY )
    {
        stressZY.subset ( *M_sigmaZ, static_cast<UInt> ( 1 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the ZZ component of the stress by copying it to an external vector
    /*!
     * @param stressZZ vector to be filled with the ZZ component of the stress
     */
    void exportStressZZ ( solutionVect_Type& stressZZ )
    {
        stressZZ.subset ( *M_sigmaZ, static_cast<UInt> ( 2 * M_FESpace->dof().numTotalDof() ) );
    }

    //! Export the ZZ component of the stress by copying it to an external vector
    /*!
     * @param stressZZ vector to be filled with the ZZ component of the stress
     */
    void exportStressVonMises ( solutionVect_Type& stressVonMises )
    {
        stressVonMises = *M_sigmaVonMises;
    }

    //! Get the global vector for the eigenvalues
    solutionVect_Type principalStresses() const
    {
        return *M_globalEigenvalues;
    }

    //@}


protected:


    //! @name Protected methods
    //@{

    //! computeDeformation: This method computes the tensor F given the displacement on the element.
    /*!
      \param NONE
    */
    void computeDisplacementGradient ( solutionVectPtr_Type grDisplX,
                                       solutionVectPtr_Type grDisplY,
                                       solutionVectPtr_Type grDisplZ );

    //! constructGlobalStressVector: This method construct the vectors \sigma_{.,i} for i=x,y,z to have for each DOF the tensor \sigma
    /*!
      \param NONE
    */
    void constructGlobalStressVector ();

    //! constructPatchAreaVector: This method build the patch area vector used in the reconstruction process
    /*!
      \param NONE
    */
    void constructPatchAreaVector ( solutionVect_Type& patchArea );

    //! reconstructElementaryVector: This method applies a reconstruction procedure on the elvec that is passed
    /*!
      \param elvecTens VectorElemental over which the reconstruction is applied
    */
    void reconstructElementaryVector ( VectorElemental& elVecSigma, const solutionVect_Type& patchArea, const feSpace_Type& feSpace );

    //@}


    //! @name Protected members
    //@{

    //! Vectorial FE space
    feSpacePtr_Type                                M_FESpace;

    //! Scalar FE space
    feSpacePtr_Type                                M_FESpaceScalar;

    //! Elementary matrix for the tensor F
    matrixPtr_Type                                 M_deformationF;

    //! Elementary matrix for the tensor F
    matrixPtr_Type                                 M_cofactorF;
    //! Elementary matrix for the tensor P
    matrixPtr_Type                                 M_firstPiola;

    //! Elementary matrix for the tensor \sigma (Cauchy tensor on the current config)
    matrixPtr_Type                                 M_sigma;

    //! Vector of the invariants of C and detF (length = 4)
    vector_Type                                    M_invariants;

    //! Vector of the eigenvalues of \sigma on the DOF (length = 3)
    vector_Type                                    M_eigenvaluesR;
    vector_Type                                    M_eigenvaluesI;

    //! Vector for the displacement field
    solutionVectPtr_Type                           M_displacement;

    //! Vector for the gradient along X of the displacement field
    solutionVectPtr_Type                           M_gradientX;

    //! Vector for the gradient along Y of the displacement field
    solutionVectPtr_Type                           M_gradientY;

    //! Vector for the gradient along Z of the displacement field
    solutionVectPtr_Type                           M_gradientZ;

    //! Vector for the X component of the stress tensor
    solutionVectPtr_Type                           M_sigmaX;

    //! Vector for the Y component of the stress tensor
    solutionVectPtr_Type                           M_sigmaY;

    //! Vector for the Z component of the stress tensor
    solutionVectPtr_Type                           M_sigmaZ;

    //! Vector for the Von Mises stress
    solutionVectPtr_Type                           M_sigmaVonMises;

    //! Vector for the eigenvalues of the Cauchy stress tensor
    solutionVectPtr_Type                           M_globalEigenvalues;

    //! The Offset parameter
    UInt                                           M_offset;

    //! The volume marker
    UInt                                           M_marker;

    //Class for material parameter
    dataPtr_Type                                   M_dataMaterial;

    //Class for analysis parameter
    analysisDataPtr_Type                           M_analysisData;

    //Displayer
    displayerPtr_Type                              M_displayer;

    //! Material class
    materialPtr_Type                               M_material;

    //@}
};

//=====================================
// Constructor
//=====================================
template <typename Mesh>
WallTensionEstimator<Mesh>::WallTensionEstimator( ) :
    M_FESpace                    ( ),
    M_FESpaceScalar              ( ),
    M_deformationF               ( ),
    M_cofactorF                  ( ),
    M_firstPiola                 ( ),
    M_sigma                      ( ),
    M_invariants                 ( ),
    M_eigenvaluesR               ( ),
    M_eigenvaluesI               ( ),
    M_displacement               ( ),
    M_gradientX                  ( ),
    M_gradientY                  ( ),
    M_gradientZ                  ( ),
    M_globalEigenvalues          ( ),
    M_sigmaX                     ( ),
    M_sigmaY                     ( ),
    M_sigmaZ                     ( ),
    M_sigmaVonMises              ( ),
    M_offset                     ( 0 ),
    M_marker                     ( 1 ),
    M_dataMaterial               ( ),
    M_analysisData               ( ),
    M_displayer                  ( ),
    M_material                   ( )
{
}



//====================================
// Public Methods
//===================================
template <typename Mesh>
void
WallTensionEstimator<Mesh >::setup ( const dataPtr_Type& dataMaterial,
                                     const analysisDataPtr_Type& tensionData,
                                     const feSpacePtr_Type& feSpace,
                                     const feSpaceETPtr_Type& feSpaceET,
                                     const commPtr_Type& comm,
                                     UInt marker)
{
    M_analysisData = tensionData;
    setup ( dataMaterial, feSpace, feSpaceET, comm, marker );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::setup ( const dataPtr_Type& dataMaterial,
                                     const feSpacePtr_Type& feSpace,
                                     const feSpaceETPtr_Type& feSpaceET,
                                     const commPtr_Type& comm,
                                     UInt marker)
{
    // Data classes & Volumes markers
    M_dataMaterial = dataMaterial;
    M_marker       = marker;

    // FESpace
    M_FESpace      = feSpace;

    // Create a scalar FESpace
    M_FESpaceScalar.reset ( new feSpace_Type ( M_FESpace->mesh(),
                                               M_FESpace->refFE(),
                                               M_FESpace->qr(),
                                               M_FESpace->bdQr(),
                                               1,
                                               comm ) );

    // Displayer
    M_displayer.reset     ( new Displayer ( comm ) );

    // Vector and Tensors
    M_sigma.reset         ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );

    M_deformationF.reset  ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_cofactorF.reset     ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_firstPiola.reset    ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );

    M_displacement.reset  ( new solutionVect_Type (*M_FESpace->mapPtr() ) );

    M_gradientX.reset     ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    M_gradientY.reset     ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    M_gradientZ.reset     ( new solutionVect_Type (*M_FESpace->mapPtr() ) );

    M_sigmaX.reset        ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    M_sigmaY.reset        ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    M_sigmaZ.reset        ( new solutionVect_Type (*M_FESpace->mapPtr() ) );

    M_sigmaVonMises.reset ( new solutionVect_Type (*M_FESpaceScalar->mapPtr() ) );

    M_globalEigenvalues.reset ( new solutionVect_Type (*M_FESpace->mapPtr() ) );

    M_invariants.resize   ( M_FESpace->fieldDim() + 1 );
    M_eigenvaluesR.resize ( M_FESpace->fieldDim() );
    M_eigenvaluesI.resize ( M_FESpace->fieldDim() );

    // Materials
    M_material.reset( new material_Type() );
    M_material->setup ( feSpace, feSpaceET, M_FESpace->mapPtr(), M_offset, M_dataMaterial, M_displayer );

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensions( )
{
    //TODO: Maybe here a case for the VonMises should be added
    *M_globalEigenvalues *= 0.0;
    if ( !M_analysisData->recoveryVariable().compare ("displacement") )
    {
        analyzeTensionsRecoveryDisplacement();
    }
    else if ( !M_analysisData->recoveryVariable().compare ("eigenvalues") )
    {
        analyzeTensionsRecoveryEigenvalues();
    }
    else
    {
        analyzeTensionsRecoveryCauchyStresses();
    }
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryDisplacement ( void )
{

    solutionVectPtr_Type grDisplX ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    solutionVectPtr_Type grDisplY ( new solutionVect_Type (*M_FESpace->mapPtr() ) );
    solutionVectPtr_Type grDisplZ ( new solutionVect_Type (*M_FESpace->mapPtr() ) );

    //Compute the deformation gradient tensor F of the displ field
    computeDisplacementGradient ( grDisplX, grDisplY, grDisplZ );

    M_gradientX = grDisplX;
    M_gradientY = grDisplY;
    M_gradientZ = grDisplZ;

    //Initially used for debug
    // M_displayer->leaderPrint(" \n*********************************\n  ");
    // M_displayer->leaderPrint("   Norm of the gradient with respect to x ", grDisplX->norm2() );
    // M_displayer->leaderPrint(" \n*********************************\n  ");
    // M_displayer->leaderPrint("   Norm of the gradient with respect to y ", grDisplY->norm2() );
    // M_displayer->leaderPrint(" \n*********************************\n  ");
    // M_displayer->leaderPrint("   Norm of the gradient with respect to z ", grDisplZ->norm2() );

    //For each of the DOF, the Cauchy tensor is computed.
    //Therefore the tensor C,P, \sigma are computed for each DOF
    UInt dim = M_FESpace->dim();

    LifeChrono chrono;

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the displacement..., ");
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    chrono.start();

    for ( UInt iDOF = 0; iDOF < ( UInt ) M_FESpace->dof().numTotalDof(); ++iDOF )
    {

        if ( M_displacement->blockMap().LID (iDOF) != -1 ) // The Global ID is on the calling processors
        {
            std::vector<Real> dX (3, 0.0);
            std::vector<Real> dY (3, 0.0);
            std::vector<Real> dZ (3, 0.0);

            //Reinitialization of matrices and arrays
            (*M_deformationF).Scale (0.0);
            (*M_cofactorF).Scale (0.0);
            (*M_firstPiola).Scale (0.0);
            (*M_sigma).Scale (0.0);

            //Extracting the gradient of U on the current DOF
            for ( UInt iComp = 0; iComp < M_FESpace->fieldDim(); ++iComp )
            {
                Int LIDid = M_displacement->blockMap().LID (iDOF + iComp * dim + M_offset);
                Int GIDid = M_displacement->blockMap().GID (LIDid);
                dX[iComp] = (*grDisplX) (GIDid); // (d_xX,d_yX,d_zX)
                dY[iComp] = (*grDisplY) (GIDid); // (d_xY,d_yY,d_zY)
                dZ[iComp] = (*grDisplZ) (GIDid); // (d_xZ,d_yZ,d_zZ)
            }

            //Fill the matrix F
            for ( UInt icoor (0); icoor < M_FESpace->fieldDim(); ++icoor )
            {
                (*M_deformationF) (icoor, 0) = dX[icoor];
                (*M_deformationF) (icoor, 1) = dY[icoor];
                (*M_deformationF) (icoor, 2) = dZ[icoor];

                (*M_deformationF) (icoor, icoor) += 1.0;
            }

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (M_invariants, *M_deformationF, *M_cofactorF);

            // Real sumI(0);
            // for( UInt i(0); i < M_invariants.size(); i++ )
            //     sumI += M_invariants[i];

            //Compute the first Piola-Kirchhoff tensor
            M_material->computeLocalFirstPiolaKirchhoffTensor (*M_firstPiola, *M_deformationF, *M_cofactorF, M_invariants, M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (*M_sigma, *M_firstPiola, M_invariants[3], *M_deformationF);

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( UInt i (0); i < M_eigenvaluesI.size(); ++i )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( M_eigenvaluesR.begin(), M_eigenvaluesR.end() );

            //Save the eigenvalues in the global vector
            for ( UInt icoor (0); icoor < M_FESpace->fieldDim(); ++icoor )
            {
                Int LIDid = M_displacement->blockMap().LID (iDOF + icoor * dim + M_offset);
                Int GIDid = M_displacement->blockMap().GID (LIDid);
                (*M_globalEigenvalues) (GIDid) = M_eigenvaluesR[icoor];
            }

        }
    }

    chrono.stop();
    M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryEigenvalues ( void )
{
    LifeChrono chrono;

    solutionVect_Type patchArea (*M_displacement, Unique, Add);
    patchArea *= 0.0;

    constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for ( UInt iq (0); iq < M_FESpace->qr().nbQuadPt(); ++iq )
    {
        refElemArea += M_FESpace->qr().weight (iq);
    }

    Real wQuad (refElemArea / M_FESpace->refFE().nbDof() );

    //Setting the quadrature Points = DOFs of the element and weight = 1
    std::vector<GeoVector> coords = M_FESpace->refFE().refCoor();
    std::vector<Real> weights (M_FESpace->fe().nbFEDof(), wQuad);
    fakeQuadratureRule.setDimensionShape ( shapeDimension (M_FESpace->refFE().shape() ), M_FESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //Set the new quadrature rule
    M_FESpace->setQuadRule (fakeQuadratureRule);

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the tensions..., " );
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    UInt totalDof = M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (M_FESpace->fe().nbFEDof(), M_FESpace->fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (M_FESpace->fe().nbFEDof(), *M_deformationF);
    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (*M_displacement, Repeated);

    VectorElemental elVecTens (M_FESpace->fe().nbFEDof(), M_FESpace->fieldDim() );

    chrono.start();

    //Loop on each volume
    for ( UInt i (0); i < M_FESpace->mesh()->numVolumes(); ++i )
    {
        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( i ) );
        elVecTens.zero();

        M_marker = M_FESpace->mesh()->volumeList ( i ).markerID();

        UInt eleID = M_FESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode (0); iNode < ( UInt ) M_FESpace->fe().nbFEDof(); ++iNode )
        {
            UInt  iloc = M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * M_FESpace->dim() + M_offset;
                dk_loc[iloc + iComp * M_FESpace->fe().nbFEDof() ] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, M_FESpace->fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF (0); nDOF < ( UInt ) M_FESpace->fe().nbFEDof(); ++nDOF )
        {
            UInt  iloc = M_FESpace->fe().patternFirst ( nDOF );

            M_sigma->Scale (0.0);
            M_firstPiola->Scale (0.0);
            M_cofactorF->Scale (0.0);

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (M_invariants, vectorDeformationF[nDOF], *M_cofactorF);

            //Compute the first Piola-Kirchhoff tensor
            M_material->computeLocalFirstPiolaKirchhoffTensor (*M_firstPiola, vectorDeformationF[nDOF], *M_cofactorF, M_invariants, M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (*M_sigma, *M_firstPiola, M_invariants[3], vectorDeformationF[nDOF]);

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( UInt i (0); i < M_eigenvaluesI.size(); ++i )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( M_eigenvaluesR.begin(), M_eigenvaluesR.end() );

            //Assembling the local vector
            for ( UInt coor (0); coor < M_eigenvaluesR.size(); ++coor )
            {
                elVecTens[iloc + coor * M_FESpace->fe().nbFEDof() ] = M_eigenvaluesR[coor];
            }
        }

        reconstructElementaryVector ( elVecTens, patchAreaR, *M_FESpace );

        //Assembling the local into global vector
        for ( UInt ic (0); ic < M_FESpace->fieldDim(); ++ic )
        {
            assembleVector (*M_globalEigenvalues, elVecTens, M_FESpace->fe(), M_FESpace->dof(), ic, M_offset +  ic * totalDof );
        }
    }

    M_globalEigenvalues->globalAssemble();

    chrono.stop();
    M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryCauchyStresses ( void )
{
    //Chrono
    LifeChrono chrono;
    chrono.start();

    //Construct stress tensor
    constructGlobalStressVector ();

    for ( UInt iDOF (0); iDOF < ( UInt ) M_FESpace->dof().numTotalDof(); ++iDOF )
    {

        if ( M_displacement->blockMap().LID (iDOF) != -1 ) // The Global ID is on the calling processors
        {

            (*M_sigma).Scale (0.0);

            //Extracting the gradient of U on the current DOF
            for ( UInt iComp = 0; iComp < M_FESpace->fieldDim(); ++iComp )
            {
                Int LIDid = M_displacement->blockMap().LID (iDOF + iComp * M_FESpace->dim() + M_offset);
                Int GIDid = M_displacement->blockMap().GID (LIDid);
                (*M_sigma) (iComp, 0) = (*M_sigmaX) (GIDid); // (d_xX,d_yX,d_zX)
                (*M_sigma) (iComp, 1) = (*M_sigmaY) (GIDid); // (d_xY,d_yY,d_zY)
                (*M_sigma) (iComp, 2) = (*M_sigmaZ) (GIDid); // (d_xZ,d_yZ,d_zZ)
            }

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( UInt i (0); i < M_eigenvaluesI.size(); ++i )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            std::sort ( M_eigenvaluesR.begin(), M_eigenvaluesR.end() );

            //Save the eigenvalues in the global vector
            for ( UInt icoor = 0; icoor < M_FESpace->fieldDim(); ++icoor )
            {
                Int LIDid = M_displacement->blockMap().LID (iDOF + icoor * M_FESpace->dim() + M_offset);
                Int GIDid = M_displacement->blockMap().GID (LIDid);
                (*M_globalEigenvalues) (GIDid) = M_eigenvaluesR[icoor];
            }
        }
    }

    chrono.stop();
    M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryVonMisesStress ()
{
    //Chrono
    LifeChrono chrono;
    chrono.start();

    //Construct stress tensor
    constructGlobalStressVector ();

    //Compute Von Mises stress
    *M_sigmaVonMises *= 0.;

    //Initialize vectors
    solutionVect_Type stressComponent1 ( *M_FESpaceScalar->mapPtr() );
    solutionVect_Type stressComponent2 ( *M_FESpaceScalar->mapPtr() );
    solutionVect_Type stressComponent3 ( *M_FESpaceScalar->mapPtr() );

    //Off-diagonal elements
    exportStressXY ( stressComponent1 );
    exportStressXZ ( stressComponent2 );
    exportStressYZ ( stressComponent3 );
    *M_sigmaVonMises += stressComponent1 * stressComponent1
                        +  stressComponent2 * stressComponent2
                        +  stressComponent3 * stressComponent3;
    *M_sigmaVonMises *= 6.;

    //Diagonal elements
    exportStressXX ( stressComponent1 );
    exportStressYY ( stressComponent2 );
    exportStressZZ ( stressComponent3 );
    *M_sigmaVonMises += ( stressComponent1 - stressComponent2 ) * ( stressComponent1 - stressComponent2 )
                        +  ( stressComponent1 - stressComponent3 ) * ( stressComponent1 - stressComponent3 )
                        +  ( stressComponent2 - stressComponent3 ) * ( stressComponent2 - stressComponent3 );

    //Final operations
    *M_sigmaVonMises *= 0.5;
    M_sigmaVonMises->sqrt();

    //Chrono
    chrono.stop();
    M_displayer->leaderPrint ("  S-  Von Mises stress computed in:             ", chrono.globalDiff ( *M_displayer->comm() ), " s\n" );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::computeDisplacementGradient ( solutionVectPtr_Type grDisplX,
                                                           solutionVectPtr_Type grDisplY,
                                                           solutionVectPtr_Type grDisplZ )
{
    //The map of the displacement field is not transformed in a Repeated map
    //because it is done inside the gradientRecovery method

    //Compute the gradient along X of the displacement field
    *grDisplX = M_FESpace->gradientRecovery (*M_displacement, 0);

    //Compute the gradient along Y of the displacement field
    *grDisplY = M_FESpace->gradientRecovery (*M_displacement, 1);

    //Compute the gradient along Z of the displacement field
    *grDisplZ = M_FESpace->gradientRecovery (*M_displacement, 2);
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::constructGlobalStressVector ()
{
    //Chrono
    LifeChrono chrono;
    chrono.start();

    //Constructing the patch area vector for reconstruction purposes
    solutionVect_Type patchArea (*M_displacement, Unique, Add);
    patchArea *= 0.0;

    constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    //Compute the area of reference element
    Real refElemArea (0);
    for (UInt iq (0); iq < M_FESpace->qr().nbQuadPt(); ++iq)
    {
        refElemArea += M_FESpace->qr().weight (iq);
    }

    //Setting the quadrature Points = DOFs of the element and weight = 1
    Real wQuad (refElemArea / M_FESpace->refFE().nbDof() );
    std::vector<Real> weights (M_FESpace->fe().nbFEDof(), wQuad);
    std::vector<GeoVector> coords = M_FESpace->refFE().refCoor();

    QuadratureRule fakeQuadratureRule;
    fakeQuadratureRule.setDimensionShape ( shapeDimension (M_FESpace->refFE().shape() ), M_FESpace->refFE().shape() );
    fakeQuadratureRule.setPoints (coords, weights);

    //Creating a copy of the FESpace
    feSpace_Type fakeFESpace ( M_FESpace->mesh(), M_FESpace->refFE(), M_FESpace->qr(), M_FESpace->bdQr(), 3, M_FESpace->map().commPtr() );

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the Cauchy stresses..., ");
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    //Set the new quadrature rule
    fakeFESpace.setQuadRule (fakeQuadratureRule);

    //Preliminary variables
    UInt totalDof = fakeFESpace.dof().numTotalDof();
    VectorElemental dk_loc (fakeFESpace.fe().nbFEDof(), fakeFESpace.fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (fakeFESpace.fe().nbFEDof(), *M_deformationF);

    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (*M_displacement, Repeated);

    //Creating the local stress tensors
    VectorElemental elVecSigmaX (fakeFESpace.fe().nbFEDof(), fakeFESpace.fieldDim() );
    VectorElemental elVecSigmaY (fakeFESpace.fe().nbFEDof(), fakeFESpace.fieldDim() );
    VectorElemental elVecSigmaZ (fakeFESpace.fe().nbFEDof(), fakeFESpace.fieldDim() );

    //Loop on each volume
    for ( UInt i (0); i < fakeFESpace.mesh()->numVolumes(); ++i )
    {
        //fakeFESpace.fe().update ( fakeFESpace.mesh()->volumeList ( i ), UPDATE_DPHI | UPDATE_WDET );
        fakeFESpace.fe().updateFirstDerivQuadPt ( fakeFESpace.mesh()->volumeList ( i ) );

        elVecSigmaX.zero();
        elVecSigmaY.zero();
        elVecSigmaZ.zero();

        M_marker = fakeFESpace.mesh()->volumeList ( i ).markerID();

        UInt eleID = fakeFESpace.fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode (0); iNode < ( UInt ) fakeFESpace.fe().nbFEDof(); ++iNode )
        {
            UInt  iloc = fakeFESpace.fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < fakeFESpace.fieldDim(); ++iComp )
            {
                UInt ig = fakeFESpace.dof().localToGlobalMap ( eleID, iloc ) + iComp * fakeFESpace.dim() + M_offset;
                dk_loc[iloc + iComp * fakeFESpace.fe().nbFEDof() ] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, fakeFESpace.fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF (0); nDOF < ( UInt ) fakeFESpace.fe().nbFEDof(); ++nDOF )
        {
            UInt  iloc = fakeFESpace.fe().patternFirst ( nDOF );

            M_sigma->Scale (0.0);
            M_firstPiola->Scale (0.0);
            M_cofactorF->Scale (0.0);

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (M_invariants, vectorDeformationF[nDOF], *M_cofactorF);

            //Compute the first Piola-Kirchhoff tensor
            M_material->computeLocalFirstPiolaKirchhoffTensor (*M_firstPiola, vectorDeformationF[nDOF], *M_cofactorF, M_invariants, M_marker);

            //Compute the Cauchy tensor
            AssemblyElementalStructure::computeCauchyStressTensor (*M_sigma, *M_firstPiola, M_invariants[3], vectorDeformationF[nDOF]);

            //Assembling the local vectors for local tensions Component X, Y, Z
            for ( UInt coor (0); coor < fakeFESpace.fieldDim(); ++coor )
            {
                (elVecSigmaX) [iloc + coor * fakeFESpace.fe().nbFEDof() ] = (*M_sigma) (coor, 0);
                (elVecSigmaY) [iloc + coor * fakeFESpace.fe().nbFEDof() ] = (*M_sigma) (coor, 1);
                (elVecSigmaZ) [iloc + coor * fakeFESpace.fe().nbFEDof() ] = (*M_sigma) (coor, 2);
            }
        }

        reconstructElementaryVector ( elVecSigmaX, patchAreaR, fakeFESpace );
        reconstructElementaryVector ( elVecSigmaY, patchAreaR, fakeFESpace );
        reconstructElementaryVector ( elVecSigmaZ, patchAreaR, fakeFESpace );

        //Assembling the three elemental vector in the three global
        for ( UInt ic = 0; ic < fakeFESpace.fieldDim(); ++ic )
        {
            assembleVector (*M_sigmaX, elVecSigmaX, fakeFESpace.fe(), fakeFESpace.dof(), ic, M_offset +  ic * totalDof );
            assembleVector (*M_sigmaY, elVecSigmaY, fakeFESpace.fe(), fakeFESpace.dof(), ic, M_offset +  ic * totalDof );
            assembleVector (*M_sigmaZ, elVecSigmaZ, fakeFESpace.fe(), fakeFESpace.dof(), ic, M_offset +  ic * totalDof );
        }
    }

    M_sigmaX->globalAssemble();
    M_sigmaY->globalAssemble();
    M_sigmaZ->globalAssemble();

    //Chrono
    chrono.stop();
    M_displayer->leaderPrint ("  S-  Cauchy stresses recovered in:             ", chrono.globalDiff ( *M_displayer->comm() ), " s\n" );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::constructPatchAreaVector ( solutionVect_Type& patchArea )
{
    solutionVect_Type patchAreaR (*M_displacement, Repeated);
    patchAreaR *= 0.0;

    Real refElemArea (0); //area of reference element
    UInt totalDof = M_FESpace->dof().numTotalDof();
    //compute the area of reference element
    for (UInt iq (0); iq < M_FESpace->qr().nbQuadPt(); ++iq)
    {
        refElemArea += M_FESpace->qr().weight (iq);
    }

    // Define a special quadrature rule for the interpolation
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_FESpace->refFE().shape() ), M_FESpace->refFE().shape() );
    Real wQuad (refElemArea / M_FESpace->refFE().nbDof() );

    for (UInt i (0); i < M_FESpace->refFE().nbDof(); ++i ) //nbRefCoor
    {
        interpQuad.addPoint (QuadraturePoint (M_FESpace->refFE().xi (i), M_FESpace->refFE().eta (i), M_FESpace->refFE().zeta (i), wQuad) );
    }

    UInt totalNumberVolumes (M_FESpace->mesh()->numVolumes() );
    UInt numberLocalDof (M_FESpace->dof().numLocalDof() );

    CurrentFE interpCFE (M_FESpace->refFE(), getGeometricMap (* (M_FESpace->mesh() ) ), interpQuad);

    // Loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberVolumes; ++iterElement)
    {
        interpCFE.update (M_FESpace->mesh()->volumeList ( iterElement ), UPDATE_WDET );

        for (UInt iterDof (0); iterDof < numberLocalDof; ++iterDof)
        {
            for (UInt iDim (0); iDim < M_FESpace->fieldDim(); ++iDim)
            {
                ID globalDofID (M_FESpace->dof().localToGlobalMap (iterElement, iterDof) + iDim * totalDof);
                patchAreaR[globalDofID] += interpCFE.measure();
            }
        }
    }

    solutionVect_Type final (patchAreaR, Unique, Add);

    patchArea.add (final);
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::reconstructElementaryVector ( VectorElemental& elVecSigma,
                                                           const solutionVect_Type& patchArea,
                                                           const feSpace_Type& feSpace )
{
    Real measure = feSpace.fe().measure();
    UInt eleID   = feSpace.fe().currentLocalId();

    for (UInt iDof = 0; iDof < feSpace.fe().nbFEDof(); ++iDof)
    {
        UInt  iloc = feSpace.fe().patternFirst ( iDof );

        for ( UInt icoor = 0;  icoor < feSpace.fieldDim(); ++icoor )
        {
            ID globalDofID (feSpace.dof().localToGlobalMap (eleID, iDof) + icoor * feSpace.dof().numTotalDof() );

            elVecSigma[iloc + icoor * feSpace.fe().nbFEDof() ] *= ( measure / patchArea[globalDofID] );
        }
    }
}

}

#endif /*WALLTENSION_H 1*/
