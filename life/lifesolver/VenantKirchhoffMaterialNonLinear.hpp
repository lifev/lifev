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
 *  @brief This file contains the definition of the St. Venant Kirchhoff material.
 *         ATTENTION: At the moment, to get \mu and \lambda it is used the default value of the marker.
 *                    This must be changed using the same approach used to built up the linear stiffness matrix!
 *                     So: getMu(marker) and the same for lambda!!!
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _NONLINVENANTKIRCHHOFFMATERIAL_H_
#define _NONLINVENANTKIRCHHOFFMATERIAL_H_

#include <life/lifesolver/StructuralMaterial.hpp>
#include <life/lifesolver/VenantKirchhoffMaterialLinear.hpp>

namespace LifeV
{
template <typename Mesh>
class VenantKirchhoffMaterialNonLinear :
        public VenantKirchhoffMaterialLinear<Mesh>
{
 //!@name Type definitions
 //@{

  public:
    typedef VenantKirchhoffMaterialLinear<Mesh>      super;
  
    typedef VenantKirchhoffElasticData               data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename boost::shared_ptr<data_Type>    dataPtr_Type;
    typedef typename boost::scoped_ptr<Displayer>    displayerPtr_Type;
  
 //@}

 //! @name Constructor &  Deconstructor
 //@{

    VenantKirchhoffMaterialNonLinear();

    virtual  ~VenantKirchhoffMaterialNonLinear();
 //@}

 //!@name Methods
 //@{

  //! Setup the created object of the class StructuralMaterial
  /*!
    \param dFespace: the FiniteElement Space
    \param monolithicMap: the MapEpetra
    \param offset: the offset parameter used assembling the matrices
  */
  void setup(const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
	     const boost::shared_ptr<const MapEpetra>&  monolithicMap,
	     const UInt offset
	     );

  //! Compute the linear part Stiffness matrix in StructuralSolver::buildSystem()
  /*!
    \param dataMaterial the class with Material properties data
  */
  void computeLinearStiffMatrix( dataPtr_Type& dataMaterial );

  //! Updates the Jacobian matrix
  /*!
    \param disp: solution at the k-th iteration of NonLinearRichardson Method
    \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
    \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
  */
  void updateJacobianMatrix( const vector_Type& disp,
			     const dataPtr_Type& dataMaterial,
			     const displayerPtr_Type& displayer);
  //! Updates the nonlinear terms in the Jacobian matrix
  /*!
    \param stiff: stiffness matrix provided from outside
    \param disp: solution at the k-th iteration of NonLinearRichardson Method
    \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
    \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
  */
  void updateNonLinearJacobianTerms(  matrixPtr_Type& stiff,
				      const vector_Type& disp,
				      const dataPtr_Type& dataMaterial,
				      const displayerPtr_Type& displayer);
  
    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in 
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );

    //! Computes the nonlinear part of Stiffness matrix in StructuralSolver given a certain displacement field. This function is used both in StructuralSolver::evalResidual and in 
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same. This is virtual and not pure virtual since in the linear St. Venant-Kirchhoff law it is not needed.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeNonLinearMatrix( matrixPtr_Type& stiff, const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );

    //! Missing Documentation !!!
    void computeKinematicsVariables( const VectorElemental& /*dk_loc*/ ){}

  //@}

};

template <typename Mesh>
VenantKirchhoffMaterialNonLinear<Mesh>::VenantKirchhoffMaterialNonLinear():
    super()
{
}

template <typename Mesh>
VenantKirchhoffMaterialNonLinear<Mesh>::~VenantKirchhoffMaterialNonLinear()
{}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::setup(const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                                   const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                                   const UInt offset
               )
{
    super::setup(dFESpace,monolithicMap,offset);
    this->M_stiff.reset               (new matrix_Type(*this->M_localMap));
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeLinearStiffMatrix(dataPtr_Type& dataMaterial)
{
  super::computeLinearStiffMatrix(dataMaterial);
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::updateJacobianMatrix(const vector_Type& disp,
								  const dataPtr_Type& dataMaterial,
								  const displayerPtr_Type& displayer)
{
    this->M_jacobian.reset(new matrix_Type(*this->M_localMap));

    *this->M_jacobian += *this->M_linearStiff;
    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
    updateNonLinearJacobianTerms(this->M_jacobian,disp,dataMaterial,displayer);
    std::cout << "*********************************" << std::endl;
    std::cout << std::endl;
    this->M_jacobian->globalAssemble();
}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::updateNonLinearJacobianTerms( matrixPtr_Type& jacobian,
                                                                           const  vector_Type& disp,
                                                                           const dataPtr_Type& dataMaterial,
                                                                           const displayerPtr_Type& displayer )
{
      displayer->leaderPrint("   NonLin S-  Updating non linear terms in the Jacobian Matrix (in updateJacobian)");
      std::cout << std::endl;

      UInt totalDof = this->M_FESpace->dof().numTotalDof();
      VectorElemental dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

      vector_Type dRep(disp, Repeated);

      //! Number of displacement components
      UInt nc = nDimensions;


      for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); ++i )
	{
	  this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

	  this->M_elmatK->zero();

	  UInt marker = this->M_FESpace->mesh()->volumeList( i ).marker();

	  Real mu = dataMaterial->mu(marker);
	  Real lambda = dataMaterial->lambda(marker);
	  
	  UInt eleID = this->M_FESpace->fe().currentLocalId();
	  UInt dim = this->M_FESpace->dim();

	  for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
	    {
	      UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
	      for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
		{
		  UInt ig = this->M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*dim + this->M_offset;
		  dk_loc[iloc + iComp*this->M_FESpace->fe().nbFEDof()] = dRep[ig]; // BASEINDEX + 1
		}
	    }

	  //  3):  \lambda * ( \tr { [\grad d^k]^T \grad \delta d }, \div v  )
	  this->M_assembler->stiff_derdiv( lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );
	  //  4):  \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )
	  this->M_assembler->stiff_dergrad( mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  // the sum of these terms is the Jacobian of the divgrad term
	  // 5):  \lambda * ( (\div u_k) \grad \delta u : \grad v  )
	  this->M_assembler->stiff_divgrad(  lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //  \lambda * ( (\div u) \grad u_k : \grad v  )
	  this->M_assembler->stiff_divgrad_2(  lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  // the sum of these terms is the Jacobian of the gradgrad term
	  // 6): 1/2 * \lambda * ( \grad u_k : \grad  u_k) *( \grad \delta u : \grad v  )
	  this->M_assembler->stiff_gradgrad(   0.5 * lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //\lambda * ( \grad u_k : \grad \delta u) *( \grad u_k : \grad v  )
	  this->M_assembler->stiff_gradgrad_2(  lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  // the sum of these terms is he jacobian of the stiff_dergrad_gradbis term
	  // 7A) : \mu *  ( \grad u^k \grad \delta u : \grad v  )
	  this->M_assembler->stiff_dergrad_gradbis(  mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //  \mu *  ( \grad \delta u \grad u^k : \grad v  )
	  this->M_assembler->stiff_dergrad_gradbis_2(  mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //  the sum of these terms is he jacobian of the stiff_dergrad_gradbis_Tr term
	  // 7B) :  \mu *  ( \grad u^k [\grad \delta u]^T : \grad v  )
	  this->M_assembler->stiff_dergrad_gradbis_Tr(  mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  // \mu *  ( \grad \delta u [\grad u^k]^T : \grad v  )
	  this->M_assembler->stiff_dergrad_gradbis_Tr_2( mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //   the sum of these terms is he jacobian of the stiff_gradgradTr_gradbis term
	  // 8) :   \mu * (  \grad d^k [\grad d^k]^T \grad \delta d : \grad v  )
	  this->M_assembler->stiff_gradgradTr_gradbis(  mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //  \mu * (  \grad d^k [\grad \delta d]^T \grad d^k : \grad v  )
	  this->M_assembler->stiff_gradgradTr_gradbis_2( mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  //  \mu * (  \grad \delta u [\grad u^k]^T \grad u^k : \grad v  )
	  this->M_assembler->stiff_gradgradTr_gradbis_3(  mu , dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

	  // assembling
	  for ( UInt ic = 0; ic < nc; ++ic )
            for ( UInt jc = 0; jc < nc; jc++ )
	      assembleMatrix( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset + jc*totalDof );
	}



}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeStiffness( const vector_Type& disp,
							       Real factor,
							       const dataPtr_Type& dataMaterial,
							       const displayerPtr_Type& displayer )
{

    this->M_stiff.reset(new matrix_Type(*this->M_localMap));

    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
    computeNonLinearMatrix(this->M_stiff,disp,factor,dataMaterial,displayer);
    std::cout << "*********************************" << std::endl;
    std::cout << std::endl;

    *this->M_stiff += *this->M_linearStiff;

    this->M_stiff->globalAssemble();
}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeNonLinearMatrix(matrixPtr_Type& stiff,
								    const vector_Type& sol,
								    Real /*factor*/,
								    const dataPtr_Type& dataMaterial,
								    const displayerPtr_Type& displayer)
{
    displayer->leaderPrint("   NonLin S-  Updating non linear terms in the Stiffness Matrix");
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

	Real mu = dataMaterial->mu(marker);
	Real lambda = dataMaterial->lambda(marker);

        UInt eleID = this->M_FESpace->fe().currentLocalId();

        for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
        {
            UInt  iloc = this->M_FESpace->fe().patternFirst( iNode );
            for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
            {
                UInt ig = this->M_FESpace->dof().localToGlobalMap( eleID, iloc ) + iComp*dim + this->M_offset;
                dk_loc[ iloc + iComp*this->M_FESpace->fe().nbFEDof() ] = dRep[ig]; // BASEINDEX + 1
            }
        }

        this->M_elmatK->zero();

        // non-linear terms of the stiffness matrix

        // 3) 1/2 * \lambda  ( \tr { [\grad d^k]^T \grad d }, \div v  )
        this->M_assembler->stiff_derdiv( 0.5 * lambda , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //4)  \mu *( [\grad d^k]^T \grad d : \grad v  )
        this->M_assembler->stiff_dergradbis( mu , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        //  5): \lambda * (div u_k) \grad d : \grad v
        this->M_assembler->stiff_divgrad( lambda, dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 6) 1/2  * \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
        this->M_assembler->stiff_gradgrad( 0.5 * lambda , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
        this->M_assembler->stiff_dergrad_gradbis( mu , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );
        // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
        this->M_assembler->stiff_dergrad_gradbis_Tr( mu , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        // 8) // \mu *  (  \grad d^k [\grad d^k]^T \grad d : \grad v  )
        this->M_assembler->stiff_gradgradTr_gradbis( mu , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

        for ( UInt ic = 0; ic < nDimensions; ++ic )
        {
            // stiff is the nonlinear matrix of the bilinear form
            for ( UInt jc = 0; jc < nDimensions; jc++ )

                assembleMatrix( *stiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);
        }

    }

}

template <typename Mesh>
inline StructuralMaterial<Mesh>* createVenantKirchhoffNonLinear() { return new VenantKirchhoffMaterialNonLinear<Mesh >(); }
namespace
{
static bool registerVKNL = StructuralMaterial<LifeV::RegionMesh3D<LinearTetra> >::StructureMaterialFactory::instance().registerProduct( "nonlinearVenantKirchhoff", &createVenantKirchhoffNonLinear<LifeV::RegionMesh3D<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __NONLINVENANTKIRCHHOFF_H */
