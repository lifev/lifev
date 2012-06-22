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
    void setup( const dataPtr_Type& dataMaterial,
		const analysisDataPtr_Type& tensionData,
	        const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
		boost::shared_ptr<Epetra_Comm>&     comm,
		UInt marker);


    //! Analyze Tensions: This method computes the Cauchy stress tensor and its principal values. It uses the displacement vector that has to be set
    /*!
      \param NONE
    */
    void analyzeTensions( void );

    //! Clean Matrices: It is to clean the term that are less than 10-7 to set them to zero.
    /*!
      \param NONE
    */
    void cleanMatrices( void );


//! @name Set Methods
//@{

     //! Set the displacement vector
     void setDisplacement(solutionVect_Type& displVect) {*M_displ = displVect;}

//@}


//! @name Get Methods
//@{

  //! Getters
  //! Get the Epetramap
  MapEpetra   const& map()     const { return *M_localMap; }

  //! Get the FESpace object
  FESpace<Mesh, MapEpetra>& dFESpace()  {return *M_FESpace;}

  //! Get the pointer to the FESpace object
  boost::shared_ptr<FESpace<Mesh, MapEpetra> > dFESpacePtr()  {return M_FESpace;}
  
  //! Get the displacement solution
  solutionVect_Type displacement() {return *M_displ;}

  //! Get the displacement solution
  solutionVect_Type gradientX() {return *M_displX;}
  solutionVect_Type gradientY() {return *M_displY;}
  solutionVect_Type gradientZ() {return *M_displZ;}

  //! Get the global vector for the eigenvalues
  solutionVect_Type principalStresses() {return *M_globalEigen;}
//@}

protected:

//! @name Protected methods
//@{
    //! computeDeformation: This method computes the tensor F given the displacement on the element.
    /*!
      \param M_deformationF the local 3x3 tensor F to be filled
      \param dk_loc the local displacement field
    */
   void computeDisplacementGradient( void );    
    
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

    //! Vector for the gradient along X of the displacement field
    solutionVectPtr_Type                            M_displX;
    //! Vector for the gradient along Y of the displacement field
    solutionVectPtr_Type                            M_displY;
    //! Vector for the gradient along Z of the displacement field
    solutionVectPtr_Type                            M_displZ;

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
WallTensionEstimator<Mesh>::WallTensionEstimator( ):
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
    M_material                   ( )
{
  
}



//====================================
// Public Methods
//===================================
template <typename Mesh>
void 
WallTensionEstimator<Mesh >::setup( const dataPtr_Type& dataMaterial,
				    const analysisDataPtr_Type& tensionData,
				    const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
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
  M_displayer.reset    (new Displayer(comm));

  // Vector and Tensors
  M_sigma.reset         (new matrix_Type( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
  M_deformationF.reset  (new matrix_Type( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
  M_cofactorF.reset     (new matrix_Type( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
  M_firstPiola.reset    (new matrix_Type( M_FESpace->fieldDim(), M_FESpace->fieldDim() ) );
  M_displ.reset         (new solutionVect_Type(*M_localMap) );
  M_displX.reset        (new solutionVect_Type(*M_localMap) );
  M_displY.reset        (new solutionVect_Type(*M_localMap) );
  M_displZ.reset        (new solutionVect_Type(*M_localMap) );
  M_globalEigen.reset   (new solutionVect_Type(*M_localMap) );
  M_invariants.resize   ( M_FESpace->fieldDim() + 1 );
  M_eigenvaluesR.resize ( M_FESpace->fieldDim() );
  M_eigenvaluesI.resize ( M_FESpace->fieldDim() );

  // Materials
  M_material.reset( material_Type::StructureMaterialFactory::instance().createObject( M_dataMaterial->solidType() ) );
  M_material->setup( dFESpace,M_localMap,M_offset, M_dataMaterial, M_displayer );
}





template <typename Mesh>
void 
WallTensionEstimator<Mesh >::analyzeTensions( void )
{
  //Compute the deformation gradient tensor F of the displ field
  computeDisplacementGradient();

  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Norm of the gradient with respect to x ", M_displX->norm2() );
  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Norm of the gradient with respect to y ", M_displY->norm2() );
  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Norm of the gradient with respect to z ", M_displZ->norm2() );

  std::string gradX = "gradientX";
  std::string gradY = "gradientY";
  std::string gradZ = "gradientZ";

  M_displX->spy(gradX);
  M_displY->spy(gradY);
  M_displZ->spy(gradZ);
  
  //For each of the DOF, the Cauchy tensor is computed. 
  //Therefore the tensor C,P, \sigma are computed for each DOF
  UInt dim = M_FESpace->dim();  

  LifeChrono chrono;

  this->M_displayer->leaderPrint(" \n*********************************\n  ");
  this->M_displayer->leaderPrint("   Performing the analysis..., EXP ");
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
	      dX[iComp] = (*M_displX)(GIDid); // (d_xX,d_yX,d_zX)
	      dY[iComp] = (*M_displY)(GIDid); // (d_xY,d_yY,d_zY)
	      dZ[iComp] = (*M_displZ)(GIDid); // (d_xZ,d_yZ,d_zZ)
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

	  cleanMatrices();
      
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
WallTensionEstimator<Mesh >::computeDisplacementGradient( void )
{
  //The map of the displacement field is not transformed in a Repeated map
  //because it is done inside the gradientRecovery method

  //Compute the gradient along X of the displacement field
  *M_displX = M_FESpace->gradientRecovery(*M_displ, 0);

  //Compute the gradient along Y of the displacement field
  *M_displY = M_FESpace->gradientRecovery(*M_displ, 1);

  //Compute the gradient along Z of the displacement field
  *M_displZ = M_FESpace->gradientRecovery(*M_displ, 2);

}  

template <typename Mesh>
void 
WallTensionEstimator<Mesh >::cleanMatrices( void )
{

  UInt N(M_deformationF->N());
  UInt M(M_deformationF->M());

  //Clearning the deformationGradientF
  for ( UInt i(0); i < N; i++ )
      for ( UInt j(0); j < M; j++ )
	if( (i-j) ) (*M_deformationF)(i,j) = 0;
  
  //Cleaning the cofactorF
  for ( UInt i(0); i < N; i++ )
      for ( UInt j(0); j < M; j++ )
	if( (i-j) ) (*M_cofactorF)(i,j) = 0;

}


}
#endif /*_WALLTENSION_H_ 1*/
