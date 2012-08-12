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
#include <lifev/structure/solver/VenantKirchhoffElasticData.hpp>
#include <lifev/structure/solver/WallTensionEstimatorData.hpp>
#include <lifev/structure/solver/StructuralMaterial.hpp>

//Materials
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialNonLinear.hpp>
#include <lifev/structure/solver/ExponentialMaterialNonLinear.hpp>
#include <lifev/structure/solver/NeoHookeanMaterialNonLinear.hpp>
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
class WallTensionEstimatorCylindricalCoordinates : public WallTensionEstimator
{
public:

//!@name Type definitions
//@{

    // Data classes
    typedef WallTensionEstimator                          super;


    typedef VenantKirchhoffElasticData                    data_Type;
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
    typedef StructuralMaterial<Mesh>                       material_Type;
    typedef boost::shared_ptr<material_Type>               materialPtr_Type;

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
    void setup( const dataPtr_Type& dataMaterial,
    		const analysisDataPtr_Type& tensionData,
    	        const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
    		boost::shared_ptr<Epetra_Comm>&     comm,
    		UInt marker);



    // //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    // /*!
    //   \param NONE FESpace<Mesh, MapEpetra>& copyFESpace 
    // */
    // void analyzeTensions( );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryDisplacementCylindrical( void );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryEigenvaluesCylindrical( void );

    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensionsRecoveryCauchyStressesCylindrical( void );

//! @name Set Methods
//@{
 
//@}


//! @name Get Methods
//@{

//@}

protected:

//! @name Protected methods
//@{

//@}

//! @name Protected members
//@{
        
    //! Elementary matrix for the tensor F
    matrixPtr_Type                                 M_deformationCylindricalF;

//@}



};

//=====================================
// Constructor
//=====================================

template <typename Mesh>
WallTensionEstimatorCylindricalCoordinates<Mesh>::WallTensionEstimatorCylindricalCoordinates( ):
  super( ),
  M_deformationCylindricalF( )
{
  
}



//====================================
// Public Methods
//===================================
template <typename Mesh>
void 
WallTensionEstimatorCylindricalCoordinates<Mesh >::setup( const dataPtr_Type& dataMaterial,
				    const analysisDataPtr_Type& tensionData,
				    const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
				    boost::shared_ptr<Epetra_Comm>&     comm,
				    UInt marker)
  
{
  super::setup(dataMaterial, tensionData, dFESpace, comm, marker);
  M_deformationCylindricalF.reset  ( new matrix_Type( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
}

template <typename Mesh>
void 
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryDisplacementCylindrical( void )
{
  
  solutionVectPtr_Type grDisplX( new solutionVect_Type(*M_localMap) );
  solutionVectPtr_Type grDisplY( new solutionVect_Type(*M_localMap) );
  solutionVectPtr_Type grDisplZ( new solutionVect_Type(*M_localMap) );

  //Compute the deformation gradient tensor F of the displ field
  computeDisplacementGradient( grDisplX, grDisplY, grDisplZ);

  //For each of the DOF, the Cauchy tensor is computed. 
  //Therefore the tensor C,P, \sigma are computed for each DOF
  UInt dim = M_FESpace->dim();  

  LifeChrono chrono;

  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Performing the analysis recovering the displacement..., ", M_dataMaterial->solidType() );
  this->M_displayer->leaderPrint(" \n*********************************\n  ");

  chrono.start();

  for ( UInt iDOF = 0; iDOF <( UInt ) this->M_FESpace->dof().numTotalDof(); iDOF++ )
    {      
      
      if ( M_displ->blockMap().LID(iDOF) != -1 ) // The Global ID is on the calling processors
	{
	  std::vector<LifeV::Real> dX(3,0.0); 
	  std::vector<LifeV::Real> dY(3,0.0); 
	  std::vector<LifeV::Real> dZ(3,0.0);

	  //Reinitialization of matrices and arrays
	  (*M_deformationF).Scale(0.0);
	  (*M_cofactorF).Scale(0.0);
	  (*M_firstPiola).Scale(0.0);
	  (*M_sigma).Scale(0.0);

	  //Extracting the gradient of U on the current DOF
	  for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
	    {		    
	      Int LIDid = M_displ->blockMap().LID(iDOF + iComp * dim + M_offset); 
	      Int GIDid = M_displ->blockMap().GID(LIDid); 
	      dX[iComp] = (*grDisplX)(GIDid); // (d_xX,d_yX,d_zX)
	      dY[iComp] = (*grDisplY)(GIDid); // (d_xY,d_yY,d_zY)
	      dZ[iComp] = (*grDisplZ)(GIDid); // (d_xZ,d_yZ,d_zZ)
	    }
	  	      
	  //Fill the matrix F
	  for( UInt icoor = 0; icoor < M_FESpace->fieldDim(); icoor++ )
	    {
	      (*M_deformationF)(icoor,0)=dX[icoor];
	      (*M_deformationF)(icoor,1)=dY[icoor];
	      (*M_deformationF)(icoor,2)=dZ[icoor];
	  
	      (*M_deformationF)(icoor,icoor) += 1.0;
	    }
	  
	  //Compute the rightCauchyC tensor
	  AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor(M_invariants, *M_deformationF, *M_cofactorF);	 
      
	  LifeV::Real sumI(0);
	  for( UInt i(0); i < M_invariants.size(); i++ )
	    sumI += M_invariants[i];
      
	  //Compute the first Piola-Kirchhoff tensor
	  M_material->computeLocalFirstPiolaKirchhoffTensor(*M_firstPiola, *M_deformationF, *M_cofactorF, M_invariants, M_marker);

	  //Compute the Cauchy tensor
	  AssemblyElementalStructure::computeCauchyStressTensor(*M_sigma, *M_firstPiola, M_invariants[3], *M_deformationF);

	  //Compute the eigenvalue
	  AssemblyElementalStructure::computeEigenvalues(*M_sigma, M_eigenvaluesR, M_eigenvaluesI);
      
	  //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
	  //Check on the imaginary part of eigen values given by the Lapack method 
	  Real sum(0);
	  for( int i=0; i < M_eigenvaluesI.size(); i++ )
	    sum += std::abs(M_eigenvaluesI[i]);
	  ASSERT_PRE( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );
	  
	  orderEigenvalues( M_eigenvaluesR );

	  //Save the eigenvalues in the global vector
	  for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	    {
	      Int LIDid = M_displ->blockMap().LID(iDOF + icoor * dim + M_offset); 
	      Int GIDid = M_displ->blockMap().GID(LIDid); 
	      (*M_globalEigen)(GIDid) = M_eigenvaluesR[icoor];
	    }
	
	}
    }

  chrono.stop();
  this->M_displayer->leaderPrint("Analysis done in: ", chrono.diff());

}


template <typename Mesh>
void 
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryEigenvaluesCylindrical( void )
{

  LifeChrono chrono;

  solutionVect_Type patchArea(*M_displ,Unique,Add);
  patchArea *= 0.0;

  constructPatchAreaVector( patchArea );
  
  //Before assembling the reconstruction process is done
  solutionVect_Type patchAreaR(patchArea,Repeated);

  QuadratureRule fakeQuadratureRule;

  Real refElemArea(0); //area of reference element
  //compute the area of reference element
  for(UInt iq=0; iq< M_FESpace->qr().nbQuadPt(); iq++)
    refElemArea += M_FESpace->qr().weight(iq);

  Real wQuad(refElemArea/M_FESpace->refFE().nbDof());

  //Setting the quadrature Points = DOFs of the element and weight = 1
  std::vector<GeoVector> coords = M_FESpace->refFE().refCoor();
  std::vector<Real> weights(M_FESpace->fe().nbFEDof(), wQuad);
  fakeQuadratureRule.setDimensionShape ( shapeDimension(M_FESpace->refFE().shape()), M_FESpace->refFE().shape() );
  fakeQuadratureRule.setPoints(coords,weights);

  //Set the new quadrature rule
  M_FESpace->setQuadRule(fakeQuadratureRule);

  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Performing the analysis recovering the tensions..., ", M_dataMaterial->solidType() );
  this->M_displayer->leaderPrint(" \n*********************************\n  ");

  UInt totalDof = M_FESpace->dof().numTotalDof();
  VectorElemental dk_loc(M_FESpace->fe().nbFEDof(), nDimensions);

  //Vectors for the deformation tensor
  std::vector<matrix_Type> vectorDeformationF(M_FESpace->fe().nbFEDof(),*M_deformationF);
  //Copying the displacement field into a vector with repeated map for parallel computations
  solutionVect_Type dRep(*M_displ, Repeated);

  VectorElemental elVecTens(this->M_FESpace->fe().nbFEDof(), nDimensions);

  chrono.start();

  //Loop on each volume
  for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); ++i )
    {
      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );
      elVecTens.zero();

      M_marker = M_FESpace->mesh()->volumeList( i ).marker();

      UInt eleID = M_FESpace->fe().currentLocalId();

      //Extracting the local displacement
      for ( UInt iNode = 0; iNode < ( UInt ) M_FESpace->fe().nbFEDof(); iNode++ )
	{
	  UInt  iloc = M_FESpace->fe().patternFirst( iNode );

	  for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
	    {
	      UInt ig = M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*M_FESpace->dim() + this->M_offset;
	      dk_loc[iloc + iComp*M_FESpace->fe().nbFEDof()] = dRep[ig];
	    }
	}

      //Compute the element tensor F
      AssemblyElementalStructure::computeLocalDeformationGradient( dk_loc, vectorDeformationF, M_FESpace->fe() );

      //Compute the local vector of the principal stresses
      for( UInt nDOF=0; nDOF < ( UInt ) M_FESpace->fe().nbFEDof(); nDOF++ )
	{
	  UInt  iloc = M_FESpace->fe().patternFirst( nDOF );

	  M_sigma->Scale(0.0);
	  M_firstPiola->Scale(0.0);
	  M_cofactorF->Scale(0.0);

	  //Compute the rightCauchyC tensor
	  AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor(M_invariants, vectorDeformationF[nDOF], *M_cofactorF);	  

	  //Compute the first Piola-Kirchhoff tensor
	  M_material->computeLocalFirstPiolaKirchhoffTensor(*M_firstPiola, vectorDeformationF[nDOF], *M_cofactorF, M_invariants, M_marker);

	  //Compute the Cauchy tensor
	  AssemblyElementalStructure::computeCauchyStressTensor(*M_sigma, *M_firstPiola, M_invariants[3], vectorDeformationF[nDOF]);
	 
	  //Compute the eigenvalue
	  AssemblyElementalStructure::computeEigenvalues(*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

	  //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
	  //Check on the imaginary part of eigen values given by the Lapack method 
	  Real sum(0);
	  for( int i=0; i < M_eigenvaluesI.size(); i++ )
	    sum += std::abs(M_eigenvaluesI[i]);
	  ASSERT_PRE( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

	  orderEigenvalues( M_eigenvaluesR );

	  //Assembling the local vector
	  for ( int coor=0; coor < M_eigenvaluesR.size(); coor++ )
	    {
	      elVecTens[iloc + coor*M_FESpace->fe().nbFEDof()] = M_eigenvaluesR[coor];
	    }
	}

      reconstructElementaryVector( elVecTens, patchAreaR, i );

      //Assembling the local into global vector
      for ( UInt ic = 0; ic < nDimensions; ++ic )
	  assembleVector(*M_globalEigen, elVecTens, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic*totalDof );
    }
  
  M_globalEigen->globalAssemble();

  chrono.stop();
  this->M_displayer->leaderPrint("Analysis done in: ", chrono.diff());
}

template <typename Mesh>
void 
WallTensionEstimatorCylindricalCoordinates<Mesh >::analyzeTensionsRecoveryCauchyStressesCylindrical( void )
{

  LifeChrono chrono;

  chrono.start();
  UInt dim = M_FESpace->dim();  
  
  //Construction of the global tensionsVector
  solutionVectPtr_Type sigmaX( new solutionVect_Type(*M_localMap) );
  solutionVectPtr_Type sigmaY( new solutionVect_Type(*M_localMap) );
  solutionVectPtr_Type sigmaZ( new solutionVect_Type(*M_localMap) );
  
  constructGlobalStressVector(*sigmaX,*sigmaY,*sigmaZ);


  for ( UInt iDOF = 0; iDOF <( UInt ) this->M_FESpace->dof().numTotalDof(); iDOF++ )
    {      
      
      if ( M_displ->blockMap().LID(iDOF) != -1 ) // The Global ID is on the calling processors
	{

	  (*M_sigma).Scale(0.0);

	  //Extracting the gradient of U on the current DOF
	  for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
	    {		    
	      Int LIDid = M_displ->blockMap().LID(iDOF + iComp * dim + M_offset); 
	      Int GIDid = M_displ->blockMap().GID(LIDid); 
	      (*M_sigma)(iComp,0) = (*sigmaX)(GIDid); // (d_xX,d_yX,d_zX)
	      (*M_sigma)(iComp,1) = (*sigmaY)(GIDid); // (d_xY,d_yY,d_zY)
	      (*M_sigma)(iComp,2) = (*sigmaZ)(GIDid); // (d_xZ,d_yZ,d_zZ)
	    }

	  //Compute the eigenvalue
	  AssemblyElementalStructure::computeEigenvalues(*M_sigma, M_eigenvaluesR, M_eigenvaluesI);

	  //The Cauchy tensor is symmetric and therefore, the eigenvalues are real
	  //Check on the imaginary part of eigen values given by the Lapack method 
	  Real sum(0);
	  for( int i=0; i < M_eigenvaluesI.size(); i++ )
	    sum += std::abs(M_eigenvaluesI[i]);
	  ASSERT_PRE( sum < 1e-6 , "The eigenvalues of the Cauchy stress tensors have to be real!" );

	  orderEigenvalues( M_eigenvaluesR );

	  //Save the eigenvalues in the global vector
	  for( UInt icoor = 0; icoor < nDimensions; ++icoor )
	    {
	      Int LIDid = M_displ->blockMap().LID(iDOF + icoor * dim + M_offset); 
	      Int GIDid = M_displ->blockMap().GID(LIDid); 
	      (*M_globalEigen)(GIDid) = M_eigenvaluesR[icoor];
	    }
	
	}
    }

  chrono.stop();
  this->M_displayer->leaderPrint("Analysis done in: ", chrono.diff());
	 
}

template <typename Mesh>
void 
WallTensionEstimatorCylindricalCoordinates<Mesh >::constructGlobalStressVector( solutionVect_Type& sigmaX, solutionVect_Type& sigmaY, solutionVect_Type& sigmaZ )
{

  //Creating the local stress tensors
  VectorElemental elVecSigmaX(this->M_FESpace->fe().nbFEDof(), nDimensions);
  VectorElemental elVecSigmaY(this->M_FESpace->fe().nbFEDof(), nDimensions);
  VectorElemental elVecSigmaZ(this->M_FESpace->fe().nbFEDof(), nDimensions);

  LifeChrono chrono;

  //Constructing the patch area vector for reconstruction purposes
  solutionVect_Type patchArea(*M_displ,Unique,Add);
  patchArea *= 0.0;

  constructPatchAreaVector( patchArea );
  
  //Before assembling the reconstruction process is done
  solutionVect_Type patchAreaR(patchArea,Repeated);

  QuadratureRule fakeQuadratureRule;

  Real refElemArea(0); //area of reference element
  //compute the area of reference element
  for(UInt iq=0; iq< M_FESpace->qr().nbQuadPt(); iq++)
    refElemArea += M_FESpace->qr().weight(iq);

  Real wQuad(refElemArea/M_FESpace->refFE().nbDof());

  //Setting the quadrature Points = DOFs of the element and weight = 1
  std::vector<GeoVector> coords = M_FESpace->refFE().refCoor();
  std::vector<Real> weights(M_FESpace->fe().nbFEDof(), wQuad);
  fakeQuadratureRule.setDimensionShape ( shapeDimension(M_FESpace->refFE().shape()), M_FESpace->refFE().shape() );
  fakeQuadratureRule.setPoints(coords,weights);

  //Set the new quadrature rule
  M_FESpace->setQuadRule(fakeQuadratureRule);

  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Performing the analysis recovering the Cauchy stresses..., ", M_dataMaterial->solidType() );
  this->M_displayer->leaderPrint(" \n*********************************\n  ");

  UInt totalDof = M_FESpace->dof().numTotalDof();
  VectorElemental dk_loc(M_FESpace->fe().nbFEDof(), nDimensions);

  //Vectors for the deformation tensor
  std::vector<matrix_Type> vectorDeformationF(M_FESpace->fe().nbFEDof(),*M_deformationF);
  //Copying the displacement field into a vector with repeated map for parallel computations
  solutionVect_Type dRep(*M_displ, Repeated);

  chrono.start();

  //Loop on each volume
  for ( UInt i = 0; i < M_FESpace->mesh()->numVolumes(); ++i )
    {
      M_FESpace->fe().updateFirstDerivQuadPt( M_FESpace->mesh()->volumeList( i ) );
      
      elVecSigmaX.zero();
      elVecSigmaY.zero();
      elVecSigmaZ.zero();

      M_marker = M_FESpace->mesh()->volumeList( i ).marker();

      UInt eleID = M_FESpace->fe().currentLocalId();

      //Extracting the local displacement
      for ( UInt iNode = 0; iNode < ( UInt ) M_FESpace->fe().nbFEDof(); iNode++ )
	{
	  UInt  iloc = M_FESpace->fe().patternFirst( iNode );

	  for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
	    {
	      UInt ig = M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*M_FESpace->dim() + this->M_offset;
	      dk_loc[iloc + iComp*M_FESpace->fe().nbFEDof()] = dRep[ig];
	    }
	}

      //Compute the element tensor F
      AssemblyElementalStructure::computeLocalDeformationGradient( dk_loc, vectorDeformationF, M_FESpace->fe() );

      //Compute the local vector of the principal stresses
      for( UInt nDOF=0; nDOF < ( UInt ) M_FESpace->fe().nbFEDof(); nDOF++ )
	{
	  UInt  iloc = M_FESpace->fe().patternFirst( nDOF );

	  M_sigma->Scale(0.0);
	  M_firstPiola->Scale(0.0);
	  M_cofactorF->Scale(0.0);

	  //Compute the rightCauchyC tensor
	  AssemblyElementalStructure::computeInvariantsRightCauchyGreenTensor(M_invariants, vectorDeformationF[nDOF], *M_cofactorF);	  

	  //Compute the first Piola-Kirchhoff tensor
	  M_material->computeLocalFirstPiolaKirchhoffTensor(*M_firstPiola, vectorDeformationF[nDOF], *M_cofactorF, M_invariants, M_marker);

	  //Compute the Cauchy tensor
	  AssemblyElementalStructure::computeCauchyStressTensor(*M_sigma, *M_firstPiola, M_invariants[3], vectorDeformationF[nDOF]);

	  //Assembling the local vectors for local tensions Component X
	  for ( int coor=0; coor < nDimensions; coor++ )
	      (elVecSigmaX)[iloc + coor*M_FESpace->fe().nbFEDof()] = (*M_sigma)(coor,0);

	  //Assembling the local vectors for local tensions Component Y
	  for ( int coor=0; coor < nDimensions; coor++ )
	      (elVecSigmaY)[iloc + coor*M_FESpace->fe().nbFEDof()] = (*M_sigma)(coor,1);

	  //Assembling the local vectors for local tensions Component Z
	  for ( int coor=0; coor < nDimensions; coor++ )
	      (elVecSigmaZ)[iloc + coor*M_FESpace->fe().nbFEDof()] = (*M_sigma)(coor,2);

	}

      reconstructElementaryVector( elVecSigmaX, patchAreaR, i );
      reconstructElementaryVector( elVecSigmaY, patchAreaR, i );
      reconstructElementaryVector( elVecSigmaZ, patchAreaR, i );

      //Assembling the three elemental vector in the three global
      for ( UInt ic = 0; ic < nDimensions; ++ic )
	{
	  assembleVector(sigmaX, elVecSigmaX, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic*totalDof );
	  assembleVector(sigmaY, elVecSigmaY, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic*totalDof );
	  assembleVector(sigmaZ, elVecSigmaZ, M_FESpace->fe(), M_FESpace->dof(), ic, this->M_offset +  ic*totalDof );
	}
    }


  sigmaX.globalAssemble();
  sigmaY.globalAssemble();
  sigmaZ.globalAssemble();
}    

}

#endif /*_WALLTENSION_H_ 1*/
