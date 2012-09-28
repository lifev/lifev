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
    @file
    @brief This file contains an Oseen equation solver class with shape derivative
           for fluid structure interaction problem

    @author MiGlobaluel A. Fernandez <miGlobaluel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @date 09-06-2003

    @author G. Fourestey
    @date 00-02-2007

    @contributor Zhen Wang <zhen.wang@emory.edu>

 */


#ifndef OSEENSOLVERBOUNDARYDERIVATIVE_H
#define OSEENSOLVERBOUNDARYDERIVATIVE_H 1


#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/navier_stokes/solver/OseenSolver.hpp>
#include <lifev/core/mesh/ElementShapes.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/fem/AssemblyElementalBoundary.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

//! @class OseenSolverBoundaryDerivative
/*!

    @brief This class contains an Oseen equation solver class with shape derivative
           for fluid structure interaction problem

    @author MiGlobaluel A. Fernandez <miGlobaluel.fernandez@inria.fr>
            Christoph Winkelmann <christoph.winkelmann@epfl.ch>
    @date 09-06-2003
    @author G. Fourestey
    @date 02-2007

    @contributor Zhen Wang <zhen.wang@emory.edu>

 */


  template< typename MeshType, typename SolverType > // = LifeV::SolverAztecOO >
class OseenSolverBoundaryDerivative:
        public OseenSolver< MeshType, SolverType >
{

public:

    //! @name Public Types
    //@{

    typedef MeshType                                          mesh_Type;
    typedef SolverType                                        linearSolver_Type;
    typedef OseenSolver< mesh_Type, linearSolver_Type >       oseenSolver_Type;
    typedef typename oseenSolver_Type::vector_Type            vector_Type;
    typedef typename oseenSolver_Type::vectorPtr_Type         vectorPtr_Type;
    typedef typename oseenSolver_Type::matrix_Type            matrix_Type;
    typedef typename oseenSolver_Type::matrixPtr_Type         matrixPtr_Type;
    typedef typename oseenSolver_Type::data_Type              data_Type;
    typedef typename oseenSolver_Type::preconditioner_Type    preconditioner_Type;
    typedef typename oseenSolver_Type::preconditionerPtr_Type preconditionerPtr_type;
    typedef typename oseenSolver_Type::bcHandler_Type         bcHandler_Type;

    //@}

    //! @name Constructors & Destructor
    //@{

    //! Empty constructor
    OseenSolverBoundaryDerivative();

    //! Constructor
    /*!
        @param dataType OseenData class
        @param velocityFESpace Velocity FE space
        @param pressureFESpace Pressure FE space
        @param communicator MPI communicator
        @param lagrangeMultiplier Lagrange multiplier
     */
    OseenSolverBoundaryDerivative( boost::shared_ptr<data_Type>    dataType,
				   FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
				   FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
				   boost::shared_ptr<Epetra_Comm>& communicator,
				   const UInt boundaryFlag,
				   const Int                       lagrangeMultiplier = 0);

    //! Virtual destructor
  virtual ~OseenSolverBoundaryDerivative() {}

    //@}

    //! Build linear system.
  virtual void buildSystem( const Real coeffStrain, const Real coeffDiv, const Real coeffTransp );

    //! Setup Stabilization Fluxes
  void setupStabilizationFluxes( std::vector<UInt> fluxesFlags, const Real coeffStab );

    //! Update system    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
     */
  virtual void updateSystem( const Real         alpha,
                               const vector_Type& betaVector,
			     const vector_Type& sourceVector )
  { oseenSolver_Type::updateSystem( alpha, betaVector, sourceVector ); }

    /*!
        @param alpha
        @param betaVector
        @param sourceVector
     */
    virtual void updateSystem( const Real         alpha,
                               const vector_Type& betaVector,
                               const vector_Type& sourceVector,
                               const vector_Type& sourceLBVector );

    //! Update system
    /*!
        @param alpha
        @param betaVector
        @param sourceVector
        @param matrix
        @param un
     */
    virtual void updateSystem( const Real         alpha,
                               const vector_Type& betaVector,
                               const vector_Type& sourceVector,
                               const vector_Type& sourceLBVector,
                               matrixPtr_Type     matrix,
                               vectorPtr_Type     un );



    void addTranspiration(const matrixPtr_Type& matrix,
                                   const vector_Type& beta,
			           const Real& coef,
			           const BCBase& boundaryCondition );

    void addStabilizationFluxes(const matrixPtr_Type& matrix,
                                   const vector_Type& beta,
			           const Real& coef,
			           const BCBase& boundaryCondition );


private:


    typedef CurrentFE                                    currentFE_type;
    typedef CurrentBoundaryFE                            currentBdFE_type;  
    typedef boost::shared_ptr<currentFE_type>            currentFE_ptrType;
    typedef boost::shared_ptr<currentBdFE_type>          currentBdFE_ptrType;

    typedef MatrixElemental                                      localMatrix_type;
    typedef boost::shared_ptr<localMatrix_type>          localMatrix_ptrType;

    typedef VectorElemental                                      localVector_type;
    typedef boost::shared_ptr<localVector_type>          localVector_ptrType;

    //Boundary for the laplace beltrami
    UInt   M_boundaryFlagLB;

    //Boundary Flag for the fluxes
    std::vector<UInt>  M_fluxesFlags;
    bool M_useStabilizationFluxes;
    Real M_stabilizationFluxesCoeff;

    currentBdFE_ptrType M_IPFaceCFE;

    currentFE_ptrType M_IPQuad1CFE;
    currentFE_ptrType M_IPQuad2CFE;

    currentFE_ptrType M_IP1CFE;
    currentFE_ptrType M_IP2CFE;

    currentFE_ptrType M_IPBetaCFE;

    // local matrix for the Galerkin stencil entries
    localMatrix_ptrType M_localIPGalerkin_11;
    localMatrix_ptrType M_localIPGalerkin_22;

    // local matrix for the Extended stencil entries
    localMatrix_ptrType M_localIPExtended_12;
    localMatrix_ptrType M_localIPExtended_21;

    // CurrentBdFE for the Transpiration
    currentBdFE_ptrType M_transpirationCBdFE;

    // CurrentFE for the Transpiration
    currentFE_ptrType M_transpirationCFE;

    // Local matrix for the Transpiration
    localMatrix_ptrType M_localTranspiration;

    Real M_transpirationCoeff;

    // Chronos   
    LifeChrono M_transpirationAssemblyChrono;

    BCHandler M_bcHLB;

    // CurrentBdFE for the Laplace-Beltrami
    currentBdFE_ptrType M_laplaceBeltramiCBdFE;

    // CurrentFE for the Laplace-Beltrami
    currentFE_ptrType M_laplaceBeltramiCFE;

    // Local matrix for the Laplace-Beltrami
    localMatrix_ptrType M_localLaplaceBeltrami;

    // Chronos   
    LifeChrono M_laplaceBeltramiAssemblyChrono;

    //! LaplaceBeltrami matrix
    matrixPtr_Type                 M_matrixLB;

    //!
    bool M_isTranspirationOn;

};    

// ===================================================
// Constructors & Destructor
// ===================================================

template<typename MeshType, typename SolverType>
OseenSolverBoundaryDerivative<MeshType, SolverType>::
OseenSolverBoundaryDerivative( boost::shared_ptr<data_Type>    dataType,
                      FESpace<mesh_Type, MapEpetra>&  velocityFESpace,
                      FESpace<mesh_Type, MapEpetra>&  pressureFESpace,
                      boost::shared_ptr<Epetra_Comm>& communicator,		      
		      const UInt                   boundaryFlag,
		      const Int                       lagrangeMultiplier ):
        oseenSolver_Type ( dataType,
                           velocityFESpace,
                           pressureFESpace,
                           communicator,
                           lagrangeMultiplier ),
	//M_laplaceBeltramiCFE(),
        //M_laplaceBeltramiCBdFE(),
        //M_localLaplaceBeltrami(),
        //M_laplaceBeltramiAssemblyChrono(),
	//M_transpirationCFE(),
        //M_transpirationCBdFE(),
        //M_localTranspiration(),
	//M_isTranspirationOn(),
        //M_transpirationAssemblyChrono(),
	M_boundaryFlagLB(boundaryFlag),
        M_useStabilizationFluxes(false),
        M_stabilizationFluxesCoeff(0),
	//M_fluxesFlags(),
	//M_bcHLB(),
	//M_transpirationCoeff(),
	M_matrixLB()
{
    
    M_laplaceBeltramiCFE.reset(new currentFE_type(velocityFESpace.refFE(),
						      velocityFESpace.fe().geoMap(),
						  velocityFESpace.qr()));
    M_laplaceBeltramiCBdFE.reset(new currentBdFE_type(velocityFESpace.refFE().boundaryFE(),
							  velocityFESpace.fe().geoMap().boundaryMap(), 
						      velocityFESpace.bdQr()));
    M_localLaplaceBeltrami.reset(new localMatrix_type(velocityFESpace.fe().nbFEDof(),
							  velocityFESpace.fieldDim(),
						      velocityFESpace.fieldDim()));
    
    M_transpirationCFE.reset(new currentFE_type(velocityFESpace.refFE(),
						      velocityFESpace.fe().geoMap(),
						  velocityFESpace.qr()));
    M_transpirationCBdFE.reset(new currentBdFE_type(velocityFESpace.refFE().boundaryFE(),
							  velocityFESpace.fe().geoMap().boundaryMap(), 
						      velocityFESpace.bdQr()));
    M_localTranspiration.reset(new localMatrix_type(velocityFESpace.fe().nbFEDof(),
							  velocityFESpace.fieldDim(),
						      velocityFESpace.fieldDim()));


}


template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
setupStabilizationFluxes( std::vector<UInt> fluxesFlags , const Real coeffStab )
{

  //M_useStabilizationFluxes=true;
    M_fluxesFlags = fluxesFlags;
    M_stabilizationFluxesCoeff=coeffStab;

}

template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
buildSystem(const Real coeffStrain, const Real coeffDiv , const Real coeffTransp )
{
    
    OseenSolver<MeshType, SolverType>::buildSystem();

    M_transpirationCoeff = coeffTransp;

    //BCHandler bcH_LB;

    BCFunctionBase uNull;
    M_bcHLB.addBC( "LB",  M_boundaryFlagLB,     Natural,   Full,     uNull, 3 );

    M_bcHLB.bcUpdate(*(this->M_velocityFESpace.mesh()),this->M_velocityFESpace.feBd(),this->M_velocityFESpace.dof());

    //M_boundaryMembrane = bcH_LB[0];

    M_matrixLB.reset( new matrix_Type( this->M_localMap ) );

    M_laplaceBeltramiAssemblyChrono.start();

    // Some constants
    const UInt nbElements(M_bcHLB[0].list_size());    
    const UInt fieldDim(this->M_velocityFESpace.fieldDim());
    const UInt nbTotalDof(this->M_velocityFESpace.dof().numTotalDof());
    const UInt nbQuadPt(this->M_velocityFESpace.bdQr().nbQuadPt());


    const BCIdentifierNatural* pId;
    ID ibF, ibE, ibFE;

    //std::cout<<"Numero di nodi sulla faccia: "<< M_laplaceBeltramiCBdFE->nbNode()<<std::endl;

    if ( M_laplaceBeltramiCBdFE->nbNode() == 3 /*P1*/) {

        typedef RegionMesh<LinearTetra>::elementShape_Type elementShape_Type;

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {

        // Pointer to the i-th identifier in the list
        pId = static_cast< const BCIdentifierNatural* >( M_bcHLB[0][ iterElement ] );
	
	// Number of the current boundary face
	ibF = pId->id();
	ibE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementIdentity();
	ibFE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

	// Updating the laplace Beltrami current BdFE and FE
	M_laplaceBeltramiCBdFE->updateMeasNormalQuadPt( this->M_velocityFESpace.mesh()->boundaryFace( ibF ) );
	//M_laplaceBeltramiCFE->updateFirstDeriv( this->M_velocityFESpace.mesh()->element( ibE ) );
	M_laplaceBeltramiCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_ONLY_CELL_NODES );

        QuadratureRule faceQR1("custom quad 1",TETRA,3,0,0);
        for (UInt iQuad(0); iQuad< nbQuadPt; ++iQuad) // Here we do not use UInt because of KNM, but we should
        {
            Real x(0.0),y(0.0),z(0.0);
            M_laplaceBeltramiCFE->coorBackMap( M_laplaceBeltramiCBdFE->quadPt(iQuad,0),
                                       M_laplaceBeltramiCBdFE->quadPt(iQuad,1),
                                       M_laplaceBeltramiCBdFE->quadPt(iQuad,2),
                                       x,y,z);
            QuadraturePoint newPoint(x,y,z,this->M_velocityFESpace.bdQr().weight(iQuad));
            faceQR1.addPoint(newPoint);
        }
        M_laplaceBeltramiCFE->setQuadRule(faceQR1);
	M_laplaceBeltramiCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_DPHI | UPDATE_WDET );

	// Clean the local matrix
        M_localLaplaceBeltrami->zero();
	
        // local stiffness

	AssemblyElementalBoundary<elementShape_Type>::instance().stiffStrainBoundary(*M_localLaplaceBeltrami,
								    *M_laplaceBeltramiCFE,
								    *M_laplaceBeltramiCBdFE, 
								    coeffStrain ,ibFE,fieldDim);
	//AssemblyElementalBoundary<elementShape_Type>::instance().stiffStrainNormalBoundary(*M_localLaplaceBeltrami,
	//							    *M_laplaceBeltramiCFE,
	//							    *M_laplaceBeltramiCBdFE, 
	//							    -(1.0/6.0)*coeffStrain ,ibFE,fieldDim);
	AssemblyElementalBoundary<elementShape_Type>::instance().divDivBoundary(*M_localLaplaceBeltrami,*M_laplaceBeltramiCFE,*M_laplaceBeltramiCBdFE, coeffDiv ,ibFE,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            assembleMatrix( *M_matrixLB,
                            *M_localLaplaceBeltrami,
                            *M_laplaceBeltramiCFE,
                            *M_laplaceBeltramiCFE,
                            this->M_velocityFESpace.dof(),
                            this->M_velocityFESpace.dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*nbTotalDof, iFieldDim*nbTotalDof );
        }
    }

    }else {

    typedef RegionMesh<QuadraticTetra>::elementShape_Type elementShape_Type;

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {

        // Pointer to the i-th identifier in the list
        pId = static_cast< const BCIdentifierNatural* >( M_bcHLB[0][ iterElement ] );
	
	// Number of the current boundary face
	ibF = pId->id();
	ibE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementIdentity();
	ibFE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

	// Updating the laplace Beltrami current BdFE and FE
	M_laplaceBeltramiCBdFE->updateMeasNormalQuadPt( this->M_velocityFESpace.mesh()->boundaryFace( ibF ) );
	//M_laplaceBeltramiCFE->updateFirstDeriv( this->M_velocityFESpace.mesh()->element( ibE ) );
	M_laplaceBeltramiCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_ONLY_CELL_NODES );

        QuadratureRule faceQR1("custom quad 1",TETRA,3,0,0);
        for (UInt iQuad(0); iQuad< nbQuadPt; ++iQuad) // Here we do not use UInt because of KNM, but we should
        {
            Real x(0.0),y(0.0),z(0.0);
            M_laplaceBeltramiCFE->coorBackMap( M_laplaceBeltramiCBdFE->quadPt(iQuad,0),
                                       M_laplaceBeltramiCBdFE->quadPt(iQuad,1),
                                       M_laplaceBeltramiCBdFE->quadPt(iQuad,2),
                                       x,y,z);
            QuadraturePoint newPoint(x,y,z,this->M_velocityFESpace.bdQr().weight(iQuad));
            faceQR1.addPoint(newPoint);
        }
        M_laplaceBeltramiCFE->setQuadRule(faceQR1);
	M_laplaceBeltramiCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_DPHI | UPDATE_WDET );

	// Clean the local matrix
        M_localLaplaceBeltrami->zero();
	
        // local stiffness

	AssemblyElementalBoundary<elementShape_Type>::instance().stiffStrainBoundary(*M_localLaplaceBeltrami,
								    *M_laplaceBeltramiCFE,
								    *M_laplaceBeltramiCBdFE, 
								    coeffStrain ,ibFE,fieldDim);
	//AssemblyElementalBoundary<elementShape_Type>::instance().stiffStrainNormalBoundary(*M_localLaplaceBeltrami,
	//							    *M_laplaceBeltramiCFE,
	//							    *M_laplaceBeltramiCBdFE, 
	//							    -(1.0/6.0)*coeffStrain ,ibFE,fieldDim);
	AssemblyElementalBoundary<elementShape_Type>::instance().divDivBoundary(*M_localLaplaceBeltrami,
										    *M_laplaceBeltramiCFE,
										    *M_laplaceBeltramiCBdFE, 
										    coeffDiv ,ibFE,fieldDim);

        // Assembly
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            assembleMatrix( *M_matrixLB,
                            *M_localLaplaceBeltrami,
                            *M_laplaceBeltramiCFE,
                            *M_laplaceBeltramiCFE,
                            this->M_velocityFESpace.dof(),
                            this->M_velocityFESpace.dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*nbTotalDof, iFieldDim*nbTotalDof );
        }
    }
    
    }

    M_matrixLB->globalAssemble();

    //M_matrixLB->spy("MatrixLB");

    M_laplaceBeltramiAssemblyChrono.stop(); 
}


template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
updateSystem( const Real         alpha,
              const vector_Type& betaVector,
              const vector_Type& sourceVector, 
	      const vector_Type& sourceLBVector )
{
    if ( this->M_matrixNoBC.get() )
        this->M_matrixNoBC.reset( new matrix_Type( this->M_localMap, this->M_matrixNoBC->meanNumEntries() ) );
    else
        this->M_matrixNoBC.reset( new matrix_Type( this->M_localMap ) );

    updateSystem( alpha, betaVector, sourceVector, sourceLBVector, this->M_matrixNoBC, this->M_un );
    if ( alpha != 0. )
        this->M_matrixNoBC->globalAssemble();

    //this->M_matrixNoBC->spy("MatrixLBNoBC");

}

template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
updateSystem( const Real         alpha,
              const vector_Type& betaVector,
              const vector_Type& sourceVector,
	      const vector_Type& sourceLBVector,
              matrixPtr_Type     matrixNoBC,
              vectorPtr_Type     un )
{
    OseenSolver<MeshType, SolverType>::updateSystem( alpha, betaVector, sourceVector, matrixNoBC, un );

    *matrixNoBC += *M_matrixLB ;

    //addTranspiration( matrixNoBC, -sourceLBVector , M_transpirationCoeff, M_bcHLB[0] );

    this->M_rightHandSideNoBC +=  *M_matrixLB * sourceLBVector;

}


template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
addTranspiration(const matrixPtr_Type& matrix,
		 const vector_Type& beta,
		 const Real& coef,
		 const BCBase& boundaryCondition)
{
    if (beta.mapType() != Repeated)
    {
        addTranspiration(matrix,vector_Type(beta,Repeated),coef,boundaryCondition);
        return;
    }

    //typedef LifeV::RegionMesh3D<lifeV::LinearTetra>::ElementShape geoShape_Type;

    const UInt nbElements(boundaryCondition.list_size());    
    const UInt fieldDim(this->M_velocityFESpace.fieldDim());
    const UInt nbTotalDof(this->M_velocityFESpace.dof().numTotalDof());
    const UInt nbQuadPt(this->M_velocityFESpace.bdQr().nbQuadPt());
    const UInt nbFEDof(this->M_velocityFESpace.feBd().nbNode());

    // Temporaries
    Real localValue(0.0);
    std::vector<Real> betaN(nbQuadPt,0.0);
    Real hFace2(0.0);

    const BCIdentifierNatural* pId;
    ID ibF, ibE, ibFE, iDofE, jDofE;

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {

        // Pointer to the i-th identifier in the list
        pId = static_cast< const BCIdentifierNatural* >( boundaryCondition[ iterElement ] );
	
	// Number of the current boundary face
	ibF = pId->id();
	ibE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementIdentity();
	ibFE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementPosition(); // local id of the face in its adjacent element

	// Updating the laplace Beltrami current BdFE and FE
	M_transpirationCBdFE->updateMeasNormalQuadPt( this->M_velocityFESpace.mesh()->boundaryFace( ibF ) );
	//M_laplaceBeltramiCFE->updateFirstDeriv( this->M_velocityFESpace.mesh()->element( ibE ) );
	M_transpirationCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_ONLY_CELL_NODES );

        QuadratureRule faceQR1("custom quad 1",TETRA,3,0,0);
        for (int iQuad(0); iQuad< nbQuadPt; ++iQuad) // Here we do not use UInt because of KNM, but we should
        {
            Real x(0.0),y(0.0),z(0.0);
            M_transpirationCFE->coorBackMap( M_transpirationCBdFE->quadPt(iQuad,0),
                                       M_transpirationCBdFE->quadPt(iQuad,1),
                                       M_transpirationCBdFE->quadPt(iQuad,2),
                                       x,y,z);
            QuadraturePoint newPoint(x,y,z,this->M_velocityFESpace.bdQr().weight(iQuad));
            faceQR1.addPoint(newPoint);
        }
        M_transpirationCFE->setQuadRule(faceQR1);
	M_transpirationCFE->update( this->M_velocityFESpace.mesh()->element( ibE ),UPDATE_DPHI | UPDATE_WDET);

        for (UInt iQuadPt(0); iQuadPt<nbQuadPt; ++iQuadPt)
        {
            betaN[iQuadPt]=0.0;
            for (UInt iDof(0); iDof< nbFEDof; ++iDof)
            {
	        iDofE = RegionMesh<LinearTetra>::elementShape_Type::faceToPoint( ibFE , iDof ); // local vertex number (in element)

                for (UInt iDim(0); iDim<3; ++iDim)
                {
                    betaN[iQuadPt] += beta[this->M_velocityFESpace.dof().localToGlobalMap(ibE,iDofE) 
					   + nbTotalDof*iDim]
		                      * M_transpirationCBdFE->phi(iDof,iQuadPt);
		}
            }
            betaN[iQuadPt] = std::fabs(betaN[iQuadPt]);
        }

        // Now we can start the assembly. There are 4 parts, depending on
        // which sides we consider ( side1 with side1, side1 with side2,...)
        // We have also to be carefull about the signs, because we want the jumps
        // Therefore, there is a "-" when considering two different sides.

        // Zero out the local matrices
        M_localTranspiration->zero();
 

        // Loop on the components
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            // Extract the views
            localMatrix_type::matrix_view viewTranspiration = M_localTranspiration->block(iFieldDim,iFieldDim);

            for (UInt iDofFace(0); iDofFace< nbFEDof; ++iDofFace)
            {
	      iDofE = RegionMesh<LinearTetra>::elementShape_Type::faceToPoint( ibFE, iDofFace ); // local vertex number (in element)

                for (UInt jDofFace(0); jDofFace<nbFEDof; ++jDofFace)
                {
                    localValue = 0.0;
		
		    jDofE = RegionMesh<LinearTetra>::elementShape_Type::faceToPoint( ibFE, jDofFace ); // local vertex number (in element)

                    for (UInt iQuadPt(0); iQuadPt< nbQuadPt; ++iQuadPt)
                    {
                        for (UInt iDim(0); iDim<3; ++iDim)
                        {
                            localValue += betaN[iQuadPt]
                                             * M_transpirationCFE->dphi(iDofE,iDim,iQuadPt)
			                     * M_transpirationCBdFE->phi(jDofFace,iQuadPt)
                                             * M_transpirationCBdFE->weightMeas(iQuadPt);                            
                        }
                    }

                    // Here we put the values in the local matrices
                    viewTranspiration(iDofE,jDofE) += coef*localValue;

                }
            }
        }

        // Global Assembly
        for (UInt iFieldDim(0); iFieldDim<fieldDim; ++iFieldDim)
        {
            assembleMatrix( *matrix,
                            *M_localTranspiration,
                            *M_transpirationCFE,
                            *M_transpirationCFE,
                            this->M_velocityFESpace.dof(),
                            this->M_velocityFESpace.dof(),
                            iFieldDim, iFieldDim,
                            iFieldDim*this->M_velocityFESpace.dof().numTotalDof(), iFieldDim*this->M_velocityFESpace.dof().numTotalDof() );
        }
    }

}



  /*template<typename MeshType, typename SolverType>
void
OseenSolverBoundaryDerivative<MeshType, SolverType>::
addStabilizationFluxes(const matrixPtr_Type& matrix,
		 const vector_Type& beta,
		 const Real& coef,
		 const BCBase& boundaryCondition)
{
    if (beta.mapType() != Repeated)
    {
        addStabilizationFluxes(matrix,vector_Type(beta,Repeated),coef,boundaryCondition);
        return;
    }

    const UInt nbBoundaryFaces(boundaryCondition.list_size());    
    const UInt nbComponents(this->M_velocityFESpace.fieldDim());
    const UInt betaTotalDof(this->M_velocityFESpace.dof().numTotalDof());
    const UInt nbQuadPt(this->M_velocityFESpace.bdQr().nbQuadPt());
    const UInt nbFEDof(this->M_velocityFESpace.feBd().nbNode());
    const UInt nbQuadPt(M_fespace->bdQr().nbQuadPt());
    const UInt nbLocalDof(M_fespace->fe().nbFEDof());
    const UInt nbLocalBetaDof(M_betaFESpace->fe().nbFEDof());

    // Temporaries
    Real localValue(0.0);
    std::vector<Real> betaN(nbQuadPt,0.0);
    Real hFace2(0.0);

    const BCIdentifierNatural* pId;
    ID ibF, ibE, ibFE, iDofE, jDofE;

    // Loop over the elements
    for (UInt iterElement(0); iterElement < nbElements; ++iterElement)
    {

        // Pointer to the i-th identifier in the list
        pId = static_cast< const BCIdentifierNatural* >( boundaryCondition[ iterElement ] );
	
	// Number of the current boundary face
	ibF = pId->id();
	ibE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementIdentity();
	ibFE = this->M_velocityFESpace.mesh()->boundaryFace( ibF ).firstAdjacentElementPosition(); // local id of the face in its adjacent element
    }

}*/
}

#endif
