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

#ifndef _WALLTENSION_H_
#define _WALLTENSION_H_ 1

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
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>
#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

//Materials
#include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialLinear.hpp>
// #include <lifev/structure/solver/isotropic/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/isotropic/ExponentialMaterialNonLinear.hpp>
// #include <lifev/structure/solver/isotropic/NeoHookeanMaterialNonLinear.hpp>

#include <lifev/eta/fem/ETFESpace.hpp>

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

    // Data classes
    typedef StructuralConstitutiveLawData                 data_Type;
    typedef WallTensionEstimatorData                      analysisData_Type;
    typedef typename boost::shared_ptr<data_Type>         dataPtr_Type;
    typedef typename boost::shared_ptr<analysisData_Type> analysisDataPtr_Type;

    //Matrices 3x3 and std::vector for the invariants
    typedef Epetra_SerialDenseMatrix                     matrix_Type;
    typedef boost::shared_ptr<matrix_Type>               matrixPtr_Type;
    typedef std::vector<LifeV::Real>                     vector_Type;
    typedef boost::shared_ptr<vector_Type>               vectorPtr_Type;

    // These two are to handle the vector displacement read from hdf5
    typedef VectorEpetra                                 solutionVect_Type;
    typedef boost::shared_ptr<VectorEpetra>              solutionVectPtr_Type;

    // Displayer and Exporter classes
    typedef typename boost::shared_ptr<const Displayer>   displayerPtr_Type;
    typedef typename boost::shared_ptr< Exporter<Mesh> >  exporterPtr_Type;

    // Materials
    typedef StructuralConstitutiveLaw<Mesh>                material_Type;
    typedef boost::shared_ptr<material_Type>               materialPtr_Type;

    typedef ETFESpace< RegionMesh<LinearTetra>, MapEpetra, 3, 3 >  ETFESpace_Type;
    typedef boost::shared_ptr<ETFESpace_Type>                      ETFESpacePtr_Type;

    typedef FESpace< RegionMesh<LinearTetra>, MapEpetra >          FESpace_Type;
    typedef boost::shared_ptr<FESpace_Type>                        FESpacePtr_Type;

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
    */
    void setup ( const dataPtr_Type& dataMaterial,
                 const analysisDataPtr_Type& tensionData,
                 const FESpacePtr_Type& dFESpace,
                 const ETFESpacePtr_Type&      dETFESpace,
                 boost::shared_ptr<Epetra_Comm>&     comm,
                 UInt marker);



    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE FESpace<Mesh, MapEpetra>& copyFESpace
    */
    virtual void analyzeTensions( );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryDisplacement ( void );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryEigenvalues ( void );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryCauchyStresses ( void );

    //! @name Set Methods
    //@{

    //! Set the displacement vector
    void setDisplacement (solutionVect_Type& displVect)
    {
        *M_displ = displVect;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getters
    //! Get the Epetramap
    MapEpetra   const& map()     const
    {
        return *M_localMap;
    }

    //! Get the FESpace object
    FESpace<Mesh, MapEpetra>& dFESpace()
    {
        return *M_FESpace;
    }

    //! Get the pointer to the FESpace object
    boost::shared_ptr<FESpace<Mesh, MapEpetra> > dFESpacePtr()
    {
        return M_FESpace;
    }

    //! Get the displacement solution
    solutionVect_Type displacement()
    {
        return *M_displ;
    }

    //Getters initially used for debug to check the gradient of the displacement
    // //! Get the displacement solution
    solutionVect_Type gradientX()
    {
        return *M_displX;
    }
    solutionVect_Type gradientY()
    {
        return *M_displY;
    }
    solutionVect_Type gradientZ()
    {
        return *M_displZ;
    }

    solutionVect_Type sigmaX()
    {
        return *M_sigmaX;
    }
    solutionVect_Type sigmaY()
    {
        return *M_sigmaY;
    }
    solutionVect_Type sigmaZ()
    {
        return *M_sigmaZ;
    }

    //! Get the global vector for the eigenvalues
    solutionVect_Type principalStresses()
    {
        return *M_globalEigen;
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
                                       solutionVectPtr_Type grDisplZ);

    //! reconstructElementaryVector: This method applies a reconstruction procedure on the elvec that is passed
    /*!
      \param elvecTens VectorElemental over which the reconstruction is applied
    */
    void reconstructElementaryVector ( VectorElemental& elVecSigma, solutionVect_Type& patchArea, UInt nVol );

    //! constructPatchAreaVector: This method build the patch area vector used in the reconstruction process
    /*!
      \param NONE
    */
    void constructPatchAreaVector ( solutionVect_Type& patchArea );

    //! orderEigenvalues it puts in an increasing order the eigenvalues
    /*!
      \param std::vector of the real part of the eigenvalues that has been found
    */
    void orderEigenvalues ( std::vector<Real>& eigenvaluesR );

    //! constructGlobalStressVector: This method construct the vectors \sigma_{.,i} for i=x,y,z to have for each DOF the tensor \sigma
    /*!
      \param NONE
    */
    void constructGlobalStressVector ( solutionVect_Type& sigmaX, solutionVect_Type& sigmaY, solutionVect_Type& sigmaZ );


    //@}

    //! @name Protected members
    //@{

    boost::shared_ptr<FESpace<Mesh, MapEpetra> >   M_FESpace;

    boost::shared_ptr<const MapEpetra>             M_localMap;

    //! Elementary matrix for the tensor F
    matrixPtr_Type                                 M_deformationF;

    //! Elementary matrix for the tensor F
    matrixPtr_Type                                 M_cofactorF;
    //! Elementary matrix for the tensor P
    matrixPtr_Type                                 M_firstPiola;

    //! Elementary matrix for the tensor \sigma (Cauchy tensor on the current config)
    matrixPtr_Type                                 M_sigma;

    //! Vector of the invariants of C and detF (length = 4)
    vector_Type                                 M_invariants;

    //! Vector of the eigenvalues of \sigma on the DOF (length = 3)
    vector_Type                                 M_eigenvaluesR;
    vector_Type                                 M_eigenvaluesI;

    //! Vector for the displacement field
    solutionVectPtr_Type                            M_displ;


    // Vector initially used for debug
    //! Elementary vector for the tensions on the element
    boost::scoped_ptr<VectorElemental>             M_elVecTens;
    // //! Vector for the gradient along X of the displacement field
    solutionVectPtr_Type                            M_displX;
    //! Vector for the gradient along Y of the displacement field
    solutionVectPtr_Type                            M_displY;
    //! Vector for the gradient along Z of the displacement field
    solutionVectPtr_Type                            M_displZ;


    //! Vector for the gradient along X of the displacement field
    solutionVectPtr_Type                            M_sigmaX;
    //! Vector for the gradient along Y of the displacement field
    solutionVectPtr_Type                            M_sigmaY;
    //! Vector for the gradient along Z of the displacement field
    solutionVectPtr_Type                            M_sigmaZ;


    //! Vector for the eigenvalues of the Cauchy stress tensor
    solutionVectPtr_Type                            M_globalEigen;

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
    M_localMap                   ( ),
    M_offset                     ( 0 ),
    M_marker                     ( 1 ),
    M_dataMaterial               ( ),
    M_analysisData               ( ),
    M_displayer                  ( ),
    M_sigma                      ( ),
    M_deformationF               ( ),
    M_cofactorF                  ( ),
    M_firstPiola                 ( ),
    M_invariants                 ( ),
    M_eigenvaluesR               ( ),
    M_eigenvaluesI               ( ),
    M_displ                      ( ),
    M_globalEigen                ( ),
    M_displX                     ( ),
    M_displY                     ( ),
    M_displZ                     ( ),
    M_sigmaX                     ( ),
    M_sigmaY                     ( ),
    M_sigmaZ                     ( ),
    M_elVecTens                  ( ),
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
                                     const FESpacePtr_Type&        dFESpace,
                                     const ETFESpacePtr_Type&      dETFESpace,
                                     boost::shared_ptr<Epetra_Comm>&     comm,
                                     UInt marker)

{

    // Data classes & Volumes markers
    M_dataMaterial = dataMaterial;
    M_analysisData = tensionData;
    M_marker = marker;

    // FESpace and EpetraMap
    M_FESpace      = dFESpace;
    M_localMap     = dFESpace->mapPtr();

    // Displayer
    M_displayer.reset    (new Displayer (comm) );

    // Vector and Tensors
    M_sigma.reset         ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_deformationF.reset  ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_cofactorF.reset     ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_firstPiola.reset    ( new matrix_Type ( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
    M_displ.reset         ( new solutionVect_Type (*M_localMap) );
    M_displX.reset        ( new solutionVect_Type (*M_localMap) );
    M_displY.reset        ( new solutionVect_Type (*M_localMap) );
    M_displZ.reset        ( new solutionVect_Type (*M_localMap) );

    M_sigmaX.reset        ( new solutionVect_Type (*M_localMap) );
    M_sigmaY.reset        ( new solutionVect_Type (*M_localMap) );
    M_sigmaZ.reset        ( new solutionVect_Type (*M_localMap) );
    M_elVecTens.reset     ( new VectorElemental (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() ) );

    M_globalEigen.reset   ( new solutionVect_Type (*M_localMap) );
    M_invariants.resize   ( M_FESpace->fieldDim() + 1 );
    M_eigenvaluesR.resize ( M_FESpace->fieldDim() );
    M_eigenvaluesI.resize ( M_FESpace->fieldDim() );

    // Materials
    M_material.reset( new material_Type() );
    M_material->setup ( dFESpace, dETFESpace, M_localMap, M_offset, M_dataMaterial, M_displayer );
}



template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensions( )
{

    *M_globalEigen *= 0.0;
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

    solutionVectPtr_Type grDisplX ( new solutionVect_Type (*M_localMap) );
    solutionVectPtr_Type grDisplY ( new solutionVect_Type (*M_localMap) );
    solutionVectPtr_Type grDisplZ ( new solutionVect_Type (*M_localMap) );

    //Compute the deformation gradient tensor F of the displ field
    computeDisplacementGradient ( grDisplX, grDisplY, grDisplZ);

    M_displX = grDisplX;
    M_displY = grDisplY;
    M_displZ = grDisplZ;

    //Initially used for debug
    // this->M_displayer->leaderPrint(" \n*********************************\n  ");
    // this->M_displayer->leaderPrint("   Norm of the gradient with respect to x ", grDisplX->norm2() );
    // this->M_displayer->leaderPrint(" \n*********************************\n  ");
    // this->M_displayer->leaderPrint("   Norm of the gradient with respect to y ", grDisplY->norm2() );
    // this->M_displayer->leaderPrint(" \n*********************************\n  ");
    // this->M_displayer->leaderPrint("   Norm of the gradient with respect to z ", grDisplZ->norm2() );

    //For each of the DOF, the Cauchy tensor is computed.
    //Therefore the tensor C,P, \sigma are computed for each DOF
    UInt dim = M_FESpace->dim();

    LifeChrono chrono;

    this->M_displayer->leaderPrint (" \n*********************************\n  ");
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the displacement..., ");
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    chrono.start();

    for ( UInt iDOF = 0; iDOF < ( UInt ) this->M_FESpace->dof().numTotalDof(); iDOF++ )
    {

        if ( M_displ->blockMap().LID (iDOF) != -1 ) // The Global ID is on the calling processors
        {
            std::vector<LifeV::Real> dX (3, 0.0);
            std::vector<LifeV::Real> dY (3, 0.0);
            std::vector<LifeV::Real> dZ (3, 0.0);

            //Reinitialization of matrices and arrays
            (*M_deformationF).Scale (0.0);
            (*M_cofactorF).Scale (0.0);
            (*M_firstPiola).Scale (0.0);
            (*M_sigma).Scale (0.0);

            //Extracting the gradient of U on the current DOF
            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                Int LIDid = M_displ->blockMap().LID (iDOF + iComp * dim + M_offset);
                Int GIDid = M_displ->blockMap().GID (LIDid);
                dX[iComp] = (*grDisplX) (GIDid); // (d_xX,d_yX,d_zX)
                dY[iComp] = (*grDisplY) (GIDid); // (d_xY,d_yY,d_zY)
                dZ[iComp] = (*grDisplZ) (GIDid); // (d_xZ,d_yZ,d_zZ)
            }

            //Fill the matrix F
            for ( UInt icoor = 0; icoor < M_FESpace->fieldDim(); icoor++ )
            {
                (*M_deformationF) (icoor, 0) = dX[icoor];
                (*M_deformationF) (icoor, 1) = dY[icoor];
                (*M_deformationF) (icoor, 2) = dZ[icoor];

                (*M_deformationF) (icoor, icoor) += 1.0;
            }

            //Compute the rightCauchyC tensor
            AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor (M_invariants, *M_deformationF, *M_cofactorF);

            // LifeV::Real sumI(0);
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
            for ( int i = 0; i < M_eigenvaluesI.size(); i++ )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            orderEigenvalues ( M_eigenvaluesR );

            //Save the eigenvalues in the global vector
            for ( UInt icoor = 0; icoor < this->M_FESpace->fieldDim(); ++icoor )
            {
                Int LIDid = M_displ->blockMap().LID (iDOF + icoor * dim + M_offset);
                Int GIDid = M_displ->blockMap().GID (LIDid);
                (*M_globalEigen) (GIDid) = M_eigenvaluesR[icoor];
            }

        }
    }

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::computeDisplacementGradient ( solutionVectPtr_Type grDisplX,
                                                           solutionVectPtr_Type grDisplY,
                                                           solutionVectPtr_Type grDisplZ)
{
    //The map of the displacement field is not transformed in a Repeated map
    //because it is done inside the gradientRecovery method

    //Compute the gradient along X of the displacement field
    *grDisplX = M_FESpace->gradientRecovery (*M_displ, 0);

    //Compute the gradient along Y of the displacement field
    *grDisplY = M_FESpace->gradientRecovery (*M_displ, 1);

    //Compute the gradient along Z of the displacement field
    *grDisplZ = M_FESpace->gradientRecovery (*M_displ, 2);

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryEigenvalues ( void )
{

    LifeChrono chrono;

    solutionVect_Type patchArea (*M_displ, Unique, Add);
    patchArea *= 0.0;

    constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < M_FESpace->qr().nbQuadPt(); iq++)
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
    VectorElemental dk_loc (M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (M_FESpace->fe().nbFEDof(), *M_deformationF);
    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (*M_displ, Repeated);

    VectorElemental elVecTens (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    chrono.start();

    //Loop on each volume
    for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); ++i )
    {
        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( i ) );
        elVecTens.zero();

        M_marker = M_FESpace->mesh()->volumeList ( i ).markerID();

        UInt eleID = M_FESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * M_FESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * M_FESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, M_FESpace->fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF = 0; nDOF < ( UInt ) M_FESpace->fe().nbFEDof(); nDOF++ )
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
            for ( int i = 0; i < M_eigenvaluesI.size(); i++ )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            orderEigenvalues ( M_eigenvaluesR );

            //Assembling the local vector
            for ( int coor = 0; coor < M_eigenvaluesR.size(); coor++ )
            {
                elVecTens[iloc + coor * M_FESpace->fe().nbFEDof()] = M_eigenvaluesR[coor];
            }
        }

        reconstructElementaryVector ( elVecTens, patchAreaR, i );

        //Assembling the local into global vector
        for ( UInt ic = 0; ic < this->M_FESpace->fieldDim(); ++ic )
        {
            assembleVector (*M_globalEigen, elVecTens, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }
    }

    M_globalEigen->globalAssemble();

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::constructPatchAreaVector ( solutionVect_Type& patchArea )
{

    solutionVect_Type patchAreaR (*M_displ, Repeated);
    patchAreaR *= 0.0;

    Real refElemArea (0); //area of reference element
    UInt totalDof = M_FESpace->dof().numTotalDof();
    //compute the area of reference element
    for (UInt iq = 0; iq < M_FESpace->qr().nbQuadPt(); iq++)
    {
        refElemArea += M_FESpace->qr().weight (iq);
    }

    // Define a special quadrature rule for the interpolation
    QuadratureRule interpQuad;
    interpQuad.setDimensionShape (shapeDimension (M_FESpace->refFE().shape() ), M_FESpace->refFE().shape() );
    Real wQuad (refElemArea / M_FESpace->refFE().nbDof() );

    for (UInt i (0); i < M_FESpace->refFE().nbDof(); ++i) //nbRefCoor
    {
        interpQuad.addPoint (QuadraturePoint (M_FESpace->refFE().xi (i), M_FESpace->refFE().eta (i), M_FESpace->refFE().zeta (i), wQuad) );
    }

    UInt totalNumberVolumes (M_FESpace->mesh()->numVolumes() );
    UInt numberLocalDof (M_FESpace->dof().numLocalDof() );

    CurrentFE interpCFE (M_FESpace->refFE(), getGeometricMap (* (M_FESpace->mesh() ) ), interpQuad);

    // Loop over the cells
    for (UInt iterElement (0); iterElement < totalNumberVolumes; iterElement++)
    {
        interpCFE.update (M_FESpace->mesh()->volumeList ( iterElement ), UPDATE_WDET );

        for (UInt iterDof (0); iterDof < numberLocalDof; iterDof++)
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
void WallTensionEstimator<Mesh >::orderEigenvalues ( std::vector<Real>& eigenvaluesR )
{
    //The number of elements is 3. Thefore, a simple bubble scheme is implemented.

    int ordered = 0;

    do
    {
        ordered = 0;
        for ( UInt i (0); i < eigenvaluesR.size() - 1; i++ )
        {
            if ( eigenvaluesR[i] > eigenvaluesR[i + 1] )
            {
                Real tmp = eigenvaluesR[i];
                eigenvaluesR[i] = eigenvaluesR[i + 1];
                eigenvaluesR[i + 1] = tmp;
                ordered = 1;
            }
        }
    }
    while ( ordered );

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::analyzeTensionsRecoveryCauchyStresses ( void )
{

    LifeChrono chrono;

    chrono.start();
    UInt dim = M_FESpace->dim();

    //Construction of the global tensionsVector
    solutionVectPtr_Type sigmaX ( new solutionVect_Type (*M_localMap) );
    solutionVectPtr_Type sigmaY ( new solutionVect_Type (*M_localMap) );
    solutionVectPtr_Type sigmaZ ( new solutionVect_Type (*M_localMap) );

    constructGlobalStressVector (*sigmaX, *sigmaY, *sigmaZ);


    for ( UInt iDOF = 0; iDOF < ( UInt ) this->M_FESpace->dof().numTotalDof(); iDOF++ )
    {

        if ( M_displ->blockMap().LID (iDOF) != -1 ) // The Global ID is on the calling processors
        {

            (*M_sigma).Scale (0.0);

            //Extracting the gradient of U on the current DOF
            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                Int LIDid = M_displ->blockMap().LID (iDOF + iComp * dim + M_offset);
                Int GIDid = M_displ->blockMap().GID (LIDid);
                (*M_sigma) (iComp, 0) = (*sigmaX) (GIDid); // (d_xX,d_yX,d_zX)
                (*M_sigma) (iComp, 1) = (*sigmaY) (GIDid); // (d_xY,d_yY,d_zY)
                (*M_sigma) (iComp, 2) = (*sigmaZ) (GIDid); // (d_xZ,d_yZ,d_zZ)
            }

            //Compute the eigenvalue
            AssemblyElementalStructure::computeEigenvalues (*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

            //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
            //Check on the imaginary part of eigen values given by the Lapack method
            Real sum (0);
            for ( int i = 0; i < M_eigenvaluesI.size(); i++ )
            {
                sum += std::abs (M_eigenvaluesI[i]);
            }
            ASSERT_PRE ( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

            orderEigenvalues ( M_eigenvaluesR );

            //Save the eigenvalues in the global vector
            for ( UInt icoor = 0; icoor < this->M_FESpace->fieldDim(); ++icoor )
            {
                Int LIDid = M_displ->blockMap().LID (iDOF + icoor * dim + M_offset);
                Int GIDid = M_displ->blockMap().GID (LIDid);
                (*M_globalEigen) (GIDid) = M_eigenvaluesR[icoor];
            }

        }
    }

    chrono.stop();
    this->M_displayer->leaderPrint ("Analysis done in: ", chrono.diff() );

}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::constructGlobalStressVector ( solutionVect_Type& sigmaX, solutionVect_Type& sigmaY, solutionVect_Type& sigmaZ )
{

    //Creating the local stress tensors
    VectorElemental elVecSigmaX (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );
    VectorElemental elVecSigmaY (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );
    VectorElemental elVecSigmaZ (this->M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    LifeChrono chrono;

    //Constructing the patch area vector for reconstruction purposes
    solutionVect_Type patchArea (*M_displ, Unique, Add);
    patchArea *= 0.0;

    constructPatchAreaVector ( patchArea );

    //Before assembling the reconstruction process is done
    solutionVect_Type patchAreaR (patchArea, Repeated);

    QuadratureRule fakeQuadratureRule;

    Real refElemArea (0); //area of reference element
    //compute the area of reference element
    for (UInt iq = 0; iq < M_FESpace->qr().nbQuadPt(); iq++)
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
    this->M_displayer->leaderPrint ("   Performing the analysis recovering the Cauchy stresses..., ");
    this->M_displayer->leaderPrint (" \n*********************************\n  ");

    UInt totalDof = M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (M_FESpace->fe().nbFEDof(), this->M_FESpace->fieldDim() );

    //Vectors for the deformation tensor
    std::vector<matrix_Type> vectorDeformationF (M_FESpace->fe().nbFEDof(), *M_deformationF);
    //Copying the displacement field into a vector with repeated map for parallel computations
    solutionVect_Type dRep (*M_displ, Repeated);

    chrono.start();

    //Loop on each volume
    for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); ++i )
    {
        M_FESpace->fe().updateFirstDerivQuadPt ( M_FESpace->mesh()->volumeList ( i ) );

        elVecSigmaX.zero();
        elVecSigmaY.zero();
        elVecSigmaZ.zero();

        M_marker = M_FESpace->mesh()->volumeList ( i ).markerID();

        UInt eleID = M_FESpace->fe().currentLocalId();

        //Extracting the local displacement
        for ( UInt iNode = 0; iNode < ( UInt ) M_FESpace->fe().nbFEDof(); iNode++ )
        {
            UInt  iloc = M_FESpace->fe().patternFirst ( iNode );

            for ( UInt iComp = 0; iComp < this->M_FESpace->fieldDim(); ++iComp )
            {
                UInt ig = M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * M_FESpace->dim() + this->M_offset;
                dk_loc[iloc + iComp * M_FESpace->fe().nbFEDof()] = dRep[ig];
            }
        }

        //Compute the element tensor F
        AssemblyElementalStructure::computeLocalDeformationGradient ( dk_loc, vectorDeformationF, M_FESpace->fe() );

        //Compute the local vector of the principal stresses
        for ( UInt nDOF = 0; nDOF < ( UInt ) M_FESpace->fe().nbFEDof(); nDOF++ )
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

            //Assembling the local vectors for local tensions Component X
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaX) [iloc + coor * M_FESpace->fe().nbFEDof()] = (*M_sigma) (coor, 0);
            }

            //Assembling the local vectors for local tensions Component Y
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaY) [iloc + coor * M_FESpace->fe().nbFEDof()] = (*M_sigma) (coor, 1);
            }

            //Assembling the local vectors for local tensions Component Z
            for ( int coor = 0; coor < this->M_FESpace->fieldDim(); coor++ )
            {
                (elVecSigmaZ) [iloc + coor * M_FESpace->fe().nbFEDof()] = (*M_sigma) (coor, 2);
            }

        }

        reconstructElementaryVector ( elVecSigmaX, patchAreaR, i );
        reconstructElementaryVector ( elVecSigmaY, patchAreaR, i );
        reconstructElementaryVector ( elVecSigmaZ, patchAreaR, i );

        //Assembling the three elemental vector in the three global
        for ( UInt ic = 0; ic < this->M_FESpace->fieldDim(); ++ic )
        {
            assembleVector (sigmaX, elVecSigmaX, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
            assembleVector (sigmaY, elVecSigmaY, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
            assembleVector (sigmaZ, elVecSigmaZ, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
        }
    }


    sigmaX.globalAssemble();
    sigmaY.globalAssemble();
    sigmaZ.globalAssemble();
}

template <typename Mesh>
void
WallTensionEstimator<Mesh >::reconstructElementaryVector ( VectorElemental& elVecSigma,
                                                           solutionVect_Type& patchArea,
                                                           UInt nVol )
{
    //UpdateElement Infos
    //M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( nVol ) );

    Real measure = M_FESpace->fe().measure();
    UInt eleID = M_FESpace->fe().currentLocalId();

    for (UInt iDof = 0; iDof < M_FESpace->fe().nbFEDof(); iDof++)
    {
        UInt  iloc = M_FESpace->fe().patternFirst ( iDof );

        for ( UInt icoor = 0;  icoor < M_FESpace->fieldDim(); icoor++ )
        {
            ID globalDofID (M_FESpace->dof().localToGlobalMap (eleID, iDof) + icoor * M_FESpace->dof().numTotalDof() );

            elVecSigma[iloc + icoor * M_FESpace->fe().nbFEDof()] *= ( measure / patchArea[globalDofID] );
        }

    }
}

}

#endif /*_WALLTENSION_H_ 1*/
