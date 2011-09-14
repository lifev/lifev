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

#ifndef _EXPONENTIALMATERIAL_H_
#define _EXPONENTIALMATERIAL_H_

#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/NeoHookeanMaterialNonLinear.hpp>

namespace LifeV
{
template <typename Mesh>
class ExponentialMaterialNonLinear : public StructuralMaterial<Mesh>
{
//!@name Type definitions
//@{

    public:
    typedef StructuralMaterial<Mesh>                 super;
    
    typedef typename super::data_Type                data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename boost::shared_ptr<data_Type>    dataPtr_Type;
    typedef typename boost::scoped_ptr<Displayer>    displayerPtr_Type;
  
    //! Definition of local tensors
    typedef KN<Real> 				     KN_Type;
    typedef boost::shared_ptr<KN_Type>		     KNPtr_Type;

    typedef KNM<Real> 				     KNM_Type;
    typedef boost::shared_ptr<KNM_Type>		     KNMPtr_Type;

    typedef KNMK<Real> 				     KNMK_Type;
    typedef boost::shared_ptr<KNMK_Type>	     KNMKPtr_Type;

//@}



//! @name Constructor &  Destructor
//@{

    ExponentialMaterialNonLinear();
	
    virtual  ~ExponentialMaterialNonLinear();

//@}



//!@name Methods
//@{

    //! Setup the created object of the class StructuralMaterial
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
	        const boost::shared_ptr<const MapEpetra>&  monolithicMap,
	        const UInt offset );


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff( dataPtr_Type& /*dataMaterial*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix( const vector_Type& disp,
                               const dataPtr_Type& dataMaterial,
                               const displayerPtr_Type& displayer );


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms( matrixPtr_Type& jacobian,
                                       const vector_Type& /*disp*/,
                                       const dataPtr_Type& /*dataMaterial*/,
                                       const displayerPtr_Type& /*displayer*/ );


    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in 
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );


    //! Computes the new Stiffness vector for Neo-Hookean and Exponential materials in StructuralSolver given a certain displacement field. 
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeVector( const vector_Type& sol,
			Real factor,
			const dataPtr_Type& dataMaterial,
			const displayerPtr_Type& displayer );


    //! Computes the deformation gradient F, the cofactor matrix Cof(F), the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralMaterial::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables( const VectorElemental& dk_loc );


    //! Computes the deformation Gradient F, the cofactor of F Cof(F), the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    void computeStress( const vector_Type& sol );


    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe( std::string const& fileNameVectStiff,
		 std::string const& fileNameJacobain );

//@}

//! @name Get Methods
//@{

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const { }

  
    //! Get the stiffness vector
    vectorPtr_Type const stiffVector() const {return M_stiff; }

//@}



protected: 

    //! Local stress vector
    boost::scoped_ptr<VectorElemental>   	      	M_elvecK;

    //! Elementary matrices
    boost::scoped_ptr<MatrixElemental>             	M_elmatK;

    //! Vector: stiffness non-linear
    vectorPtr_Type		     			M_stiff;

    //! First Piola-Kirchhoff stress tensor
    vectorPtr_Type	      		      		M_FirstPiolaKStress;

    //! Local tensors initialization
    KNMPtr_Type						M_FirstPiolaKStressEle;
    KNMKPtr_Type					M_Fk;
    KNMKPtr_Type					M_CofFk;
    KNMKPtr_Type					M_Identity;
    KNPtr_Type						M_Jack;
    KNPtr_Type						M_trCisok;
    KNPtr_Type						M_trCk;

};





template <typename Mesh>
ExponentialMaterialNonLinear<Mesh>::ExponentialMaterialNonLinear():
    super			( ),
    M_elvecK 			( ),
    M_elmatK                    ( ),
    M_stiff	     	        ( ),
    M_FirstPiolaKStress		( )
{
}





template <typename Mesh>
ExponentialMaterialNonLinear<Mesh>::~ExponentialMaterialNonLinear()
{}





template <typename Mesh>
void
ExponentialMaterialNonLinear<Mesh>::setup( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                           const boost::shared_ptr<const MapEpetra>&            monolithicMap,
                                           const UInt                                           offset )
{
    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_FESpace                     = dFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;

    M_stiff.reset                 	( new vector_Type(*this->M_localMap) );

    M_FirstPiolaKStress.reset		( new vector_Type(*this->M_localMap) );
    M_elvecK.reset			( new VectorElemental (this->M_FESpace->fe().nbFEDof(), nDimensions) );
    this->M_elmatK.reset                ( new MatrixElemental( this->M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );

    //! Local tensors initilization
    M_FirstPiolaKStressEle.reset	( new KNM_Type( nDimensions, nDimensions ) );
    M_Fk.reset			        ( new KNMK_Type( nDimensions, nDimensions,dFESpace->fe().nbQuadPt() ) );
    M_CofFk.reset		   	( new KNMK_Type( nDimensions, nDimensions,dFESpace->fe().nbQuadPt() ) );		
    M_Identity.reset 	         	( new KNMK_Type( nDimensions, nDimensions,dFESpace->fe().nbQuadPt() ) );	
    M_Jack.reset			( new KN_Type( dFESpace->fe().nbQuadPt() ) );
    M_trCisok.reset			( new KN_Type( dFESpace->fe().nbQuadPt() ) );
    M_trCk.reset			( new KN_Type( dFESpace->fe().nbQuadPt() ) );

}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::computeLinearStiff(dataPtr_Type& /*dataMaterial*/)
{
    //! Empty method for neo-hookean material
}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::updateJacobianMatrix( const vector_Type&       disp,
                                                               const dataPtr_Type&      dataMaterial,
                                                               const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset(new matrix_Type(*this->M_localMap));

    displayer->leaderPrint(" ************************************************ ");
    updateNonLinearJacobianTerms(this->M_jacobian, disp, dataMaterial, displayer);
    displayer->leaderPrint(" ************************************************ ");

	
    this->M_jacobian->globalAssemble();
}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::updateNonLinearJacobianTerms( matrixPtr_Type& 		jacobian,
                                                                       const vector_Type& 	disp,
                                                                       const dataPtr_Type& 	dataMaterial,
                                                                       const displayerPtr_Type& displayer )
{
    displayer->leaderPrint("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (Exponential)");
    std::cout << std::endl;

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc(this->M_FESpace->fe().nbFEDof(), nDimensions);

    vector_Type dRep(disp, Repeated);

    //! Number of displacement components
    UInt nc = nDimensions;

    //! Nonlinear part of jacobian
    //! loop on volumes (i)
    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); ++i )
    {
	this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );
       	this->M_elmatK->zero();

	UInt marker = this->M_FESpace->mesh()->volumeList( i ).marker();
	Real mu     	= dataMaterial->mu(marker);
	Real bulk   	= dataMaterial->bulk(marker);
	Real alpha     	= dataMaterial->alpha(marker);
	Real gamma   	= dataMaterial->gamma(marker);

       	UInt eleID = this->M_FESpace->fe().currentLocalId();

       	for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
       	{
       	    	UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );

       	    	for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
       		{
       	        	UInt ig = this->M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*this->M_FESpace->dim() + this->M_offset;
       	        	dk_loc[iloc + iComp*this->M_FESpace->fe().nbFEDof()] = dRep[ig];
       		}
       	}

	//! Computes F, Cof(F), J = det(F), Tr(C)
	computeKinematicsVariables( dk_loc );

	//! Stiffness for non-linear terms of the Exponential model
	/*! 
	 The results of the integrals are stored at each step into elmatK, until to build K matrix of the bilinear form
	*/	

	//! VOLUMETRIC PART
	//! 1. Stiffness matrix: int { 1/2 * bulk * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
	AssemblyElementalStructure::stiff_Jac_Pvol_1term( bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );
	
	//! 2. Stiffness matrix: int { 1/2 * bulk * ( 1/J- 1 - log(J)/J^2 ) * ( CofF1 [\nabla \delta]^t CofF ) : \nabla v }
	AssemblyElementalStructure::stiff_Jac_Pvol_2term( bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );  		
    
    	//! ISOCHORIC PART
	//! 1. Stiffness matrix : int { - 2/3 * zeta * alpha * J^(-5/3) * exp( gamma*( Ic_iso - 3) )* ( 1. + coefExp * Ic_iso )
	//!                      *( CofF : \nabla \delta ) ( F : \nabla \v ) }
  	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_1term( alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );
	    
   	//! 2. Stiffness matrix : int { 2 * zeta * alpha * gamma * J^(-4/3) * exp( gamma*( Ic_iso - 3) ) *
	//!			 ( F : \nabla \delta ) ( F : \nabla \v ) }
   	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_2term( alpha * gamma, gamma, (*M_Fk), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );
	    
   	//! 3. Stiffness matrix : int { 2.0/9.0 * zeta * alpha * J^-2 * Ic_iso * exp( gamma*( Ic_iso - 3) )* 
	//!                       ( 1. + gamma * Ic_iso )( CofF : \nabla \delta ) ( CofF : \nabla \v )}
 	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_3term( alpha, gamma, (*M_CofFk), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );
	    
  	//! 4. Stiffness matrix : int { -2.0/3.0 * zeta * alpha * J^(-5/3) * exp( gamma*( Ic_iso - 3) )
	//!                    * ( 1. + gamma * Ic_iso )( F : \nabla \delta ) ( CofF : \nabla \v ) }
  	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_4term( alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );
	    
  	//! 5. Stiffness matrix : int { zeta * alpha * J^(-2/3) * exp( gamma*( Ic_iso - 3)) (\nabla \delta: \nabla \v)}
   	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_5term( alpha, gamma, (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );

  	//! 6. Stiffness matrix : int { 1.0/3.0 * zeta * alpha * J^(-2) * Ic_iso *  exp(gamma*( Ic_iso - 3)) *  
	//!                       (CofF [\nabla \delta]^t CofF ) : \nabla \v  }
   	AssemblyElementalStructure::stiff_Jac_P1iso_Exp_6term( alpha, gamma, (*M_CofFk), (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );

        //! assembling
        for ( UInt ic = 0; ic < nc; ++ic )
	{			
            	for ( UInt jc = 0; jc < nc; jc++ )
	   	{
                	assembleMatrix( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof );
		}
	}
    }
}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::computeStiffness( const vector_Type& sol,
						           Real /*factor*/,
						           const dataPtr_Type& dataMaterial,
						           const displayerPtr_Type& displayer )
{
    this->M_stiff.reset(new vector_Type(*this->M_localMap));

    std::cout << std::endl;
    std::cout << "***********************************************************************************" << std::endl;
    displayer->leaderPrint(" Non-Linear S-  Computing the Exponential nonlinear stiffness vector \n");
    std::cout << "***********************************************************************************" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << "\n The Euclidean norm of the solution is:" << sol.norm2() << std::endl;
    std::cout << "----------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    UInt dim = this->M_FESpace->dim();

    VectorElemental dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );
    vector_Type disp(sol);

    vector_Type dRep(disp, Repeated);

    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); i++ )
    {
	this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

	UInt marker = this->M_FESpace->mesh()->volumeList( i ).marker();

	//Getting the proper coefficients
	Real mu     	= dataMaterial->mu(marker);
	Real bulk   	= dataMaterial->bulk(marker);
	Real alpha     	= dataMaterial->alpha(marker);
	Real gamma   	= dataMaterial->gamma(marker);

       	UInt eleID = this->M_FESpace->fe().currentLocalId();

	for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
	{
		UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );

       		for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
		{
               		UInt ig = this->M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*dim + this->M_offset;
               		dk_loc[ iloc + iComp*this->M_FESpace->fe().nbFEDof() ] = dRep[ig];
		}
	}

        this->M_elvecK->zero(); 

     	this->computeKinematicsVariables( dk_loc );

     	//! Stiffness for non-linear terms of the Neo-Hookean model
     	/*! 
     	The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
     	*/	 		
     	//! Volumetric part
     	/*! 
	Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) } 
	*/
	AssemblyElementalStructure::source_Pvol( bulk, (*M_CofFk), (*M_Jack), *this->M_elvecK,  this->M_FESpace->fe() );

	//! Isochoric part
	/*!
 	Source term P1iso_Exp: int { alpha * exp(gamma *(  Ic1_iso -3 )) * 
				( J1^(-2/3)* (F1 : \nabla v) - 1/3 * (Ic1_iso / J1) * (CofF1 : \nabla v) ) }
	*/ 
       	AssemblyElementalStructure::source_P1iso_Exp( alpha, gamma, (*M_CofFk), (*M_Fk), (*M_Jack), (*M_trCisok), *this->M_elvecK, this->M_FESpace->fe() );

     	for ( UInt ic = 0; ic < nDimensions; ++ic )
     	{            
		/*! 
		M_elvecK is assemble into *vec_stiff vector that is recall
		from updateSystem(matrix_ptrtype& mat_stiff, vector_ptr_type& vec_stiff) 
        	*/
		assembleVector( *this->M_stiff, *this->M_elvecK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic*totalDof );
     	}
     }

     this->M_stiff->globalAssemble();
}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::computeKinematicsVariables( const VectorElemental& dk_loc )
{

    Real s;

    //! loop on quadrature points (ig)
    for ( int ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
	//! loop on space coordinates (icoor)
	for ( int icoor = 0; icoor < nDimensions; icoor++ )
	{
		//! loop  on space coordinates (jcoor)
		for ( int jcoor = 0; jcoor < nDimensions; jcoor++ )
		{
			s = 0.0;
			for ( int i = 0; i < this->M_FESpace->fe().nbFEDof(); i++ )
			{
				//! \grad u^k at a quadrature point
				s += this->M_FESpace->fe().phiDer( i, jcoor, ig ) * dk_loc[ i + icoor * this->M_FESpace->fe().nbFEDof() ]; 
			}
			//! gradient of displacement
			(*M_Fk)( icoor , jcoor ,ig ) = s;  
		}
	}
    }

    //! loop on quadrature points (ig)
    for ( int ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
	//! loop on space coordinates (icoor)
	for ( int  icoor = 0;icoor < nDimensions; icoor++ )
	{
                //! deformation gradient Fk
		(*M_Fk)( icoor , icoor , ig ) +=  1.0;
	}
    }
    
    Real a,b,c,d,e,f,g,h,i;
  
    for( int ig=0; ig< this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
	a = (*M_Fk)( 0 , 0 , ig );
	b = (*M_Fk)( 0 , 1 , ig );
	c = (*M_Fk)( 0 , 2 , ig );
	d = (*M_Fk)( 1 , 0 , ig );
	e = (*M_Fk)( 1 , 1 , ig );
	f = (*M_Fk)( 1 , 2 , ig );
	g = (*M_Fk)( 2 , 0 , ig );
	h = (*M_Fk)( 2 , 1 , ig );
	i = (*M_Fk)( 2 , 2 , ig );
      
	//! determinant of deformation gradient Fk  
	(*M_Jack)(ig) = a*( e*i - f*h ) - b*( d*i - f*g ) + c*( d*h - e*g );

	(*M_CofFk)( 0 , 0 , ig ) =   ( e*i - f*h );
	(*M_CofFk)( 0 , 1 , ig ) = - ( d*i - g*f );
	(*M_CofFk)( 0 , 2 , ig ) =   ( d*h - e*g );
	(*M_CofFk)( 1 , 0 , ig ) = - ( b*i - c*h );
	(*M_CofFk)( 1 , 1 , ig ) =   ( a*i - c*g );
	(*M_CofFk)( 1 , 2 , ig ) = - ( a*h - g*b );
	(*M_CofFk)( 2 , 0 , ig ) =   ( b*f - c*e );
	(*M_CofFk)( 2 , 1 , ig ) = - ( a*f - c*d );
	(*M_CofFk)( 2 , 2 , ig ) =   ( a*e - d*b );
    }

    for ( int ig = 0; ig <  this->M_FESpace->fe().nbQuadPt()  ;ig++ )
    {
	if ((*M_Jack)(ig) < 0)
	{
		std::cout <<"************************************************************"<<  std::endl;
		std::cout <<"negative jacobian !!!!!! ERROR in computeKinematicsVariables"<<  std::endl;
		std::cout <<"************************************************************"<<  std::endl;
	}
    }
    
    //! loop on quadrature points
    for ( int ig = 0;ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
	s = 0.0;
	for ( int i = 0; i < nDimensions; i++)
	{
		for ( int j = 0; j < nDimensions; j++)
		{
		        //! trace of  C1 = (F1k^t F1k)
			s +=  (*M_Fk)( i , j , ig ) * (*M_Fk)( i , j , ig ); 
		}
	}
	(*M_trCk)( ig ) = s;
    }
    
    for ( int ig = 0; ig <  this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        //! trace of deviatoric C
	(*M_trCisok)( ig ) =  pow((*M_Jack)( ig ), -2./3.) * (*M_trCk)( ig );
    }
}





template <typename Mesh>
void ExponentialMaterialNonLinear<Mesh>::showMe( std::string const& fileNameStiff,
						 std::string const& fileNameJacobian )
{
    this->M_stiff->spy(fileNameStiff);
    this->M_jacobian->spy(fileNameJacobian);
}





template <typename Mesh>
inline StructuralMaterial<Mesh>* createExponentialMaterialNonLinear() { return new ExponentialMaterialNonLinear<Mesh >(); }
namespace
{
static bool registerEXP = StructuralMaterial<LifeV::RegionMesh3D<LinearTetra> >::StructureMaterialFactory::instance().registerProduct( "exponential", &createExponentialMaterialNonLinear<LifeV::RegionMesh3D<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __EXPONENTIALMATERIAL_H */
