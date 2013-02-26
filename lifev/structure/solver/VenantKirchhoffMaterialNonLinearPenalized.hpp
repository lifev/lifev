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
 *
 *  @maintainer  Paolo Tricerri      <paolo.tricerri@epfl.ch>
 */

#ifndef _VENANTKIRCHHOFFPENALIZED_H_
#define _VENANTKIRCHHOFFPENALIZED_H_

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"


#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>

namespace LifeV
{
template <typename Mesh>
class VenantKirchhoffMaterialNonLinearPenalized : public StructuralConstitutiveLaw<Mesh>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralConstitutiveLaw<Mesh>                 super;

    typedef typename super::data_Type                data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;
    //@}



    //! @name Constructor &  Destructor
    //@{

    VenantKirchhoffMaterialNonLinearPenalized();

    virtual  ~VenantKirchhoffMaterialNonLinearPenalized();

    //@}



    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup ( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                 const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                 const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer );


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/ );


    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const displayerPtr_Type& displayer );


    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms ( matrixPtr_Type& jacobian,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                        const displayerPtr_Type& displayer );


    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness ( const vector_Type& sol,
                            Real factor,
                            const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                            const displayerPtr_Type& displayer );

    //! Computes the new Stiffness vector for Neo-Hookean and VK-Penalized materials in StructuralSolver given a certain displacement field.
    //! This function is used both in StructuralSolver::evalResidual and in StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeVector ( const vector_Type& sol,
                         Real factor,
                         const dataPtr_Type& dataMaterial,
                         const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                         const displayerPtr_Type& displayer );

    //! Computes the deformation gradient F, the cofactor matrix Cof(F), the determinant of F (J = det(F)), the trace of right Cauchy-Green tensor tr(C)
    //! This function is used in StructuralConstitutiveLaw::computeStiffness
    /*!
      \param dk_loc: the elemental displacement
    */
    void computeKinematicsVariables ( const VectorElemental& dk_loc );

    //! Computes the deformation Gradient F, the cofactor of F Cof(F), the determinant of F J = det(F), the trace of C Tr(C).
    /*!
      \param dk_loc: local displacement vector
    */
    //void computeStress( const vector_Type& sol );

    //! ShowMe method of the class (saved on a file the stiffness vector and the jacobian)
    void showMe ( std::string const& fileNameVectStiff,
                  std::string const& fileNameJacobain );

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

    void apply ( const vector_Type& sol, vector_Type& res, const mapMarkerVolumesPtr_Type mapsMarkerVolumes) ;

    //@}



protected:

    //! Local stress vector
    boost::scoped_ptr<VectorElemental>                 M_elvecK;

    //! Elementary matrices
    boost::scoped_ptr<MatrixElemental>                 M_elmatK;

    //! Vector: stiffness non-linear
    vectorPtr_Type                         M_stiff;

    //! First Piola-Kirchhoff stress tensor
    vectorPtr_Type                                M_FirstPiolaKStress;

    //! Local tensors initialization
    boost::shared_ptr<boost::multi_array<Real, 3> > M_Fk;
    boost::shared_ptr<boost::multi_array<Real, 3> > M_FkMinusTransposed;
    boost::shared_ptr<boost::multi_array<Real, 3> > M_CofFk;
    boost::shared_ptr<boost::multi_array<Real, 3> > M_Ck;
    boost::shared_ptr<boost::multi_array<Real, 3> > M_FkCk;

    boost::shared_ptr<std::vector<Real> > M_Jack;
    boost::shared_ptr<std::vector<Real> > M_trCisok;
    boost::shared_ptr<std::vector<Real> > M_trCk;
    boost::shared_ptr<std::vector<Real> > M_trCkSquared;

};





template <typename Mesh>
VenantKirchhoffMaterialNonLinearPenalized<Mesh>::VenantKirchhoffMaterialNonLinearPenalized() :
    super                  ( ),
    M_elvecK               ( ),
    M_elmatK               ( ),
    M_stiff                ( ),
    M_FirstPiolaKStress    ( )
{
}





template <typename Mesh>
VenantKirchhoffMaterialNonLinearPenalized<Mesh>::~VenantKirchhoffMaterialNonLinearPenalized()
{}





template <typename Mesh>
void
VenantKirchhoffMaterialNonLinearPenalized<Mesh>::setup ( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                                         const boost::shared_ptr<const MapEpetra>&            monolithicMap,
                                                         const UInt  offset,
                                                         const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer  )
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;
    //    std::cout<<"I am setting up the Material"<<std::endl;

    this->M_FESpace                     = dFESpace;
    this->M_localMap                    = monolithicMap;
    this->M_offset                      = offset;

    M_stiff.reset                     ( new vector_Type (*this->M_localMap) );

    M_FirstPiolaKStress.reset        ( new vector_Type (*this->M_localMap) );
    M_elvecK.reset            ( new VectorElemental (this->M_FESpace->fe().nbFEDof(), nDimensions) );
    this->M_elmatK.reset                ( new MatrixElemental ( this->M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );

    //! Local tensors initilization
    M_Fk.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_CofFk.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_FkMinusTransposed.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_Ck.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );
    M_FkCk.reset ( new boost::multi_array<Real, 3> (boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );

    M_Jack.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_trCisok.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_trCk.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );
    M_trCkSquared.reset ( new std::vector<Real> (dFESpace->fe().nbQuadPt(), 0.0) );

}





template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::computeLinearStiff (dataPtr_Type& /*dataMaterial*/, const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/)
{
    //! Empty method for this law
}





template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::updateJacobianMatrix ( const vector_Type&       disp,
                                                                             const dataPtr_Type&      dataMaterial,
                                                                             const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                             const displayerPtr_Type& displayer )
{
    this->M_jacobian.reset (new matrix_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n************************************************\n ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, dataMaterial, mapsMarkerVolumes, displayer);
    displayer->leaderPrint (" \n************************************************\n ");


    this->M_jacobian->globalAssemble();
}





template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::updateNonLinearJacobianTerms ( matrixPtr_Type&        jacobian,
        const vector_Type&     disp,
        const dataPtr_Type&    dataMaterial,
        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
        const displayerPtr_Type& displayer )
{

    *jacobian *= 0.0;

    displayer->leaderPrint ("   Non-Linear S-  updating non linear terms in the Jacobian Matrix (VK-Penalized)");

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc (this->M_FESpace->fe().nbFEDof(), nDimensions);

    vector_Type dRep (disp, Repeated);

    //! Number of displacement components
    UInt nc = nDimensions;

    //! Nonlinear part of jacobian
    //! loop on volumes (i)

    mapIterator_Type it;

    for ( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    {

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        Real bulk = dataMaterial->bulk (marker);
        Real mu = dataMaterial->mu (marker);
        Real lambda = dataMaterial->lambda (marker);

        for ( UInt j (0); j < it->second.size(); j++ )
        {
            this->M_FESpace->fe().updateFirstDerivQuadPt ( * (it->second[j]) );

            UInt eleID = this->M_FESpace->fe().currentLocalId();

            for ( UInt iNode = 0; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof(); iNode++ )
            {
                UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * this->M_FESpace->dim() + this->M_offset;
                    dk_loc[iloc + iComp * this->M_FESpace->fe().nbFEDof()] = dRep[ig];
                }
            }

            this->M_elmatK->zero();

            //! Computes F, Cof(F), J = det(F), Tr(C)
            this->computeKinematicsVariables ( dk_loc );

            //! Stiffness for non-linear terms of the VK-Penalized model
            /*!
              The results of the integrals are stored at each step into elmatK, until to build K matrix of the bilinear form
            */

            //! VOLUMETRIC PART
            //! 1. Stiffness matrix: int { 1/2 * bulk * ( 2 - 1/J + 1/J^2 ) * ( CofF : \nabla \delta ) (CofF : \nabla v) }
            AssemblyElementalStructure::stiff_Jac_Pvol_1term ( 0.5 *  bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

            //! 2. Stiffness matrix: int { 1/2 * bulk * ( 1/J- 1 - log(J)/J^2 ) * ( CofF1 [\nabla \delta]^t CofF ) : \nabla v }
            AssemblyElementalStructure::stiff_Jac_Pvol_2term ( 0.5 *  bulk, (*M_CofFk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

            //! ISOCHORIC PART
            //! 0. Stiffness matrix : int { -(2.0/3.0) * Jk^(-2.0/3.0) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * F^-T:\nabla \delta ) * ( F - (1.0/3.0) * Ic_k * F^-T ): \nabla \v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_0term ( lambda, mu, (*M_FkMinusTransposed), (*M_Fk), (*M_Jack), (*M_trCk), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );


            // //! 1. Stiffness matrix : int { J^(-2/3) * (lambda / 2) * ( (-2/3) * Ic_k * J^(-2/3) * F^-T:\nabla \delta ) * ( F : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_1term ( ( 0.5 * lambda ), (*M_FkMinusTransposed), (*M_Fk), (*M_Jack), (*M_trCk),  *this->M_elmatK, this->M_FESpace->fe() );

            // //! 2. Stiffness matrix : int { J^(-2/3) * (lambda / 2) * ( ( 2/9 ) * J^(-2/3) * Ic_k^2 ) * ( F^-T : \nabla \delta ) ( F^-T : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_2term ( ( 0.5 * lambda ), (*M_FkMinusTransposed), (*M_Jack), (*M_trCk), *this->M_elmatK, this->M_FESpace->fe() );


            // //! 3. Stiffness matrix:int { J^(-2/3) * (lambda / 2) * ( 2 * J^(-2/3) * F : \nabla \delta ) * ( F : \nabla v) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_3term ( ( 0.5 * lambda  ), (*M_Fk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 4. Stiffness matrix:int{ J^(-2/3) * (lambda / 2) * ( -2.0/3.0 * J^(-2/3) * Ic_k ) * ( F : \nabla \delta ) * ( F^-T : \nabla \v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_4term ( ( 0.5 * lambda ), (*M_FkMinusTransposed), (*M_Fk), (*M_Jack), (*M_trCk), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 5. Stiffness matrix : int { J^(-2/3) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * \nabla \delta : \nabla v }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_5term ( lambda, mu, (*M_Jack), (*M_trCisok), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 6. Stiffness matrix : int { J^(-2.0/3.0) * ( (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * ( (-2/3) * ( F :\nabla \delta ) ) * ( F^-T : \nabla v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_6term ( lambda, mu, (*M_Jack), (*M_trCisok), (*M_Fk), (*M_FkMinusTransposed), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 7. Stiffness matrix : int { ( J^(-2.0/3.0) * (lambda/2) * Ic_isok - ( (3/2)*lambda + mu ) ) * ( (1/3) * Ic_k * ( F^-T \nabla \delta^T F-T ) : \nabla v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_7term ( lambda, mu, (*M_FkMinusTransposed), (*M_trCisok), (*M_trCk), (*M_Jack), *this->M_elmatK, this->M_FESpace->fe() );

            //! 8. Stiffness matrix : int { ( -4.0/3.0) * ( mu * J^(-4/3) ) * ( F^-T: \grad \delta ) * ( F C ) : \nabla v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_8term ( mu, (*M_Jack), (*M_FkMinusTransposed), (*M_FkCk), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 9. Stiffness matrix : int { ( 4.0/9.0) * ( mu * J^(-4/3) ) Ic_kSquared * (F^-T : \grad \delta ) * F^-T : \nabla \v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_9term ( mu, (*M_Jack), (*M_trCkSquared), (*M_FkMinusTransposed), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 10. Stiffness matrix : int { ( mu * J^(-4/3) ) * ( \nabla \delta * C ) : \nabla v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_10term ( mu, (*M_Jack), (*M_Ck), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 11. Stiffness matrix : int { ( mu * J^(-4/3) ) * (F [\nabla \delta]^T F ) : \nabla \v  }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_11term ( mu, (*M_Jack), (*M_Fk), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 12. Stiffness matrix : int  { ( mu * J^(-4/3) ) * (F * F^T * [\nabla \delta] ) : \nabla \v }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_12term ( mu, (*M_Jack), (*M_Fk), *this->M_elmatK, this->M_FESpace->fe() );

            // //! 13. Stiffness matrix : int {  ( mu * J^(-4/3) ) * ( (1/3) *  Ic_SquaredK * ( F^-T [\nabla \delta ]^T F^-T) : \nabla v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_13term ( mu, (*M_Jack), (*M_trCkSquared), (*M_FkMinusTransposed), *this->M_elmatK, this->M_FESpace->fe() );

            //! 14. Stiffness matrix : int {  ( mu * J^(-4/3) ) * ( (-4.0/3.0) * ( FkCk : \nabla \delta ) ) * F^-T : \nabla v ) }
            AssemblyElementalStructure::stiff_Jac_P1iso_VKPenalized_14term ( mu, (*M_Jack), (*M_FkCk), (*M_FkMinusTransposed), *this->M_elmatK, this->M_FESpace->fe() );

            //! assembling
            for ( UInt ic = 0; ic < nc; ++ic )
            {
                for ( UInt jc = 0; jc < nc; jc++ )
                {
                    assembleMatrix ( *jacobian, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, jc, this->M_offset +  ic * totalDof, this->M_offset +  jc * totalDof );
                }
            }
        }
    }
}





template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::computeStiffness ( const vector_Type& sol,
                                                                         Real /*factor*/,
                                                                         const dataPtr_Type& dataMaterial,
                                                                         const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                         const displayerPtr_Type& displayer )
{
    this->M_stiff.reset (new vector_Type (*this->M_localMap) );

    displayer->leaderPrint (" \n*********************************\n  ");
    displayer->leaderPrint (" Non-Linear S-  Computing the VK-Penalized nonlinear stiffness vector ");
    displayer->leaderPrint (" \n*********************************\n  ");

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    UInt dim = this->M_FESpace->dim();

    VectorElemental dk_loc ( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type dRep (sol, Repeated);

    mapIterator_Type it;

    for ( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    {

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        Real bulk = dataMaterial->bulk (marker);
        Real mu = dataMaterial->mu (marker);
        Real lambda = dataMaterial->lambda (marker);

        for ( UInt j (0); j < it->second.size(); j++ )
        {
            this->M_FESpace->fe().updateFirstDerivQuadPt ( * (it->second[j]) );

            UInt eleID = this->M_FESpace->fe().currentLocalId();

            for ( UInt iNode = 0 ; iNode < ( UInt ) this->M_FESpace->fe().nbFEDof() ; iNode++ )
            {
                UInt  iloc = this->M_FESpace->fe().patternFirst ( iNode );

                for ( UInt iComp = 0; iComp < nDimensions; ++iComp )
                {
                    UInt ig = this->M_FESpace->dof().localToGlobalMap ( eleID, iloc ) + iComp * dim + this->M_offset;
                    dk_loc[ iloc + iComp * this->M_FESpace->fe().nbFEDof() ] = dRep[ig];
                }
            }

            this->M_elvecK->zero();

            this->computeKinematicsVariables ( dk_loc );

            //! Stiffness for non-linear terms of the Neo-Hookean model
            /*!
              The results of the integrals are stored at each step into elvecK, until to build K matrix of the bilinear form
            */
            //! Volumetric part
            /*!
              Source term Pvol: int { bulk /2* (J1^2 - J1  + log(J1) ) * 1/J1 * (CofF1 : \nabla v) }
            */
            AssemblyElementalStructure::source_Pvol ( 0.5 * bulk, (*M_CofFk), (*M_Jack), *this->M_elvecK,  this->M_FESpace->fe() );

            //! Isochoric part
            /*!
              Source term P1iso_VKPenalized: int { J^(-2.0/3.0) * ( \frac{lambda}{2}*Ic_isoK -  \frac{3}{2}*lambda - mu )*
              ( (F1 : \nabla v) - 1/3 * (Ic) * (F1^-T : \nabla v) ) }
            */
            AssemblyElementalStructure::source_P1iso_VKPenalized ( lambda, mu, (*M_FkMinusTransposed), (*M_Fk), (*M_trCisok), (*M_trCk), (*M_Jack), *this->M_elvecK, this->M_FESpace->fe() );

            //! Isochoric part second term
            /*!
              Source term P2iso_VKPenalized: int { ( mu * Jk^(-4.0/3.0) )*
              ( (F*C : \nabla v) - 1/3 * (Ic_squared) * (F1^-T : \nabla v) ) }
            */
            AssemblyElementalStructure::source_P2iso_VKPenalized ( mu, (*M_FkMinusTransposed), (*M_FkCk), (*M_trCkSquared), (*M_Jack), *this->M_elvecK, this->M_FESpace->fe() );

            for ( UInt ic = 0; ic < nDimensions; ++ic )
            {
                /*!
                  M_elvecK is assemble into *vec_stiff vector that is recall
                  from updateSystem(matrix_ptrtype& mat_stiff, vector_ptr_type& vec_stiff)
                */
                assembleVector ( *this->M_stiff, *this->M_elvecK, this->M_FESpace->fe(), this->M_FESpace->dof(), ic, this->M_offset +  ic * totalDof );
            }
        }
    }

    this->M_stiff->globalAssemble();
}





template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::computeKinematicsVariables ( const VectorElemental& dk_loc )
{

    Real s;
    Real sum;

    //! loop on quadrature points (ig)
    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        //! loop on space coordinates (icoor)
        for ( UInt icoor = 0; icoor < nDimensions; icoor++ )
        {
            //! loop  on space coordinates (jcoor)
            for ( UInt jcoor = 0; jcoor < nDimensions; jcoor++ )
            {
                s = 0.0;
                for ( UInt i = 0; i < this->M_FESpace->fe().nbFEDof(); i++ )
                {
                    //! \grad u^k at a quadrature point
                    s += this->M_FESpace->fe().phiDer ( i, jcoor, ig ) * dk_loc[ i + icoor * this->M_FESpace->fe().nbFEDof() ];
                }
                //! gradient of displacement
                (*M_Fk) [ icoor ][ jcoor ][ig ] = s;
            }
        }
    }

    //! loop on quadrature points (ig)
    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        //! loop on space coordinates (icoor)
        for ( UInt  icoor = 0; icoor < nDimensions; icoor++ )
        {
            //! deformation gradient Fk
            (*M_Fk) [ icoor ][ icoor ][ ig ] +=  1.0;
        }
    }

    Real a, b, c, d, e, f, g, h, i;

    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        a = (*M_Fk) [ 0 ][ 0 ][ ig ];
        b = (*M_Fk) [ 0 ][ 1 ][ ig ];
        c = (*M_Fk) [ 0 ][ 2 ][ ig ];
        d = (*M_Fk) [ 1 ][ 0 ][ ig ];
        e = (*M_Fk) [ 1 ][ 1 ][ ig ];
        f = (*M_Fk) [ 1 ][ 2 ][ ig ];
        g = (*M_Fk) [ 2 ][ 0 ][ ig ];
        h = (*M_Fk) [ 2 ][ 1 ][ ig ];
        i = (*M_Fk) [ 2 ][ 2 ][ ig ];

        //! determinant of deformation gradient Fk
        (*M_Jack) [ig] = a * ( e * i - f * h ) - b * ( d * i - f * g ) + c * ( d * h - e * g );

        ASSERT_PRE ( (*M_Jack) [ig] > 0, "Negative Jacobian. Error!" );

        (*M_FkMinusTransposed) [ 0 ][ 0 ][ ig ] =   ( e * i - f * h ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 0 ][ 1 ][ ig ] = - ( d * i - g * f ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 0 ][ 2 ][ ig ] =   ( d * h - e * g ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 1 ][ 0 ][ ig ] = - ( b * i - c * h ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 1 ][ 1 ][ ig ] =   ( a * i - c * g ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 1 ][ 2 ][ ig ] = - ( a * h - g * b ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 2 ][ 0 ][ ig ] =   ( b * f - c * e ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 2 ][ 1 ][ ig ] = - ( a * f - c * d ) / ( (*M_Jack) [ig] );
        (*M_FkMinusTransposed) [ 2 ][ 2 ][ ig ] =   ( a * e - d * b ) / ( (*M_Jack) [ig] );

        (*M_CofFk) [ 0 ][ 0 ][ ig ] =   ( e * i - f * h );
        (*M_CofFk) [ 0 ][ 1 ][ ig ] = - ( d * i - g * f );
        (*M_CofFk) [ 0 ][ 2 ][ ig ] =   ( d * h - e * g );
        (*M_CofFk) [ 1 ][ 0 ][ ig ] = - ( b * i - c * h );
        (*M_CofFk) [ 1 ][ 1 ][ ig ] =   ( a * i - c * g );
        (*M_CofFk) [ 1 ][ 2 ][ ig ] = - ( a * h - g * b );
        (*M_CofFk) [ 2 ][ 0 ][ ig ] =   ( b * f - c * e );
        (*M_CofFk) [ 2 ][ 1 ][ ig ] = - ( a * f - c * d );
        (*M_CofFk) [ 2 ][ 2 ][ ig ] =   ( a * e - d * b );
    }

    //! Compute the tensor C
    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        for ( UInt i = 0; i < nDimensions; i++)
        {
            for ( UInt j = 0; j < nDimensions; j++)
            {
                s = 0.0;
                for ( UInt p = 0; p < nDimensions; p++ )
                {
                    s += (*M_Fk) [ p ][ i ][ ig ] * (*M_Fk) [ p ][ j ][ ig ];
                }
                (*M_Ck) [ i ][ j ][ ig ] = s;
            }
        }
    }

    //! compute the trace of C and the trace of C squared
    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        s = 0.0;
        sum = 0.0;

        for ( UInt i = 0; i < nDimensions; i++)
        {
            for ( UInt j = 0; j < nDimensions; j++)
            {
                //! trace of  C1 = (F1k^t F1k) = \sum_{i,j} F[i][j]*F[i][j]
                s +=  (*M_Fk) [ i ][ j ][ ig ] * (*M_Fk) [ i ][ j ][ ig ];
                //! trace of  C^2 = (Ck^t Ck) = \sum_{i,j} C[i][j]*C[i][j]
                sum +=  (*M_Ck) [ i ][ j ][ ig ] * (*M_Ck) [ i ][ j ][ ig ];
            }
        }

        (*M_trCk) [ ig ] = s;
        (*M_trCkSquared) [ ig ] = sum;
    }

    for ( UInt ig = 0; ig <  this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        //! trace of deviatoric C
        (*M_trCisok) [ ig ] =  pow ( (*M_Jack) [ ig ], -2.0 / 3.0) * (*M_trCk) [ ig ];
    }

    //! Compute the product FC
    for ( UInt ig = 0; ig < this->M_FESpace->fe().nbQuadPt(); ig++ )
    {
        for ( UInt i = 0; i < nDimensions; i++)
        {
            for ( UInt j = 0; j < nDimensions; j++)
            {
                s = 0.0;
                for ( UInt p = 0; p < nDimensions; p++ )
                {
                    s += (*M_Fk) [ i ][ p ][ ig ] * (*M_Ck) [ p ][ j ][ ig ];
                }
                (*M_FkCk) [ i ][ j ][ ig ] = s;
            }
        }
    }

}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::showMe ( std::string const& fileNameStiff,
                                                               std::string const& fileNameJacobian )
{
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::apply ( const vector_Type& sol, vector_Type& res, const mapMarkerVolumesPtr_Type mapsMarkerVolumes )
{
    computeStiffness (sol, 0, this->M_dataMaterial, mapsMarkerVolumes, this->M_displayer);
    res += *M_stiff;
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinearPenalized<Mesh>::computeLocalFirstPiolaKirchhoffTensor ( Epetra_SerialDenseMatrix& firstPiola,
        const Epetra_SerialDenseMatrix& tensorF,
        const Epetra_SerialDenseMatrix& cofactorF,
        const std::vector<Real>& invariants,
        const UInt marker)
{

    //Get the material parameters
    Real alpha    = this->M_dataMaterial->alpha (marker);
    Real gamma    = this->M_dataMaterial->gamma (marker);
    Real bulk   = this->M_dataMaterial->bulk (marker);


    //Computing the first term \alphaJ^{-2/3}[F-(1/3)tr(C)F^{-T}]exp(\gamma(tr(Ciso) - 3)
    Epetra_SerialDenseMatrix firstTerm (tensorF);
    Epetra_SerialDenseMatrix copyCofactorF (cofactorF);

    Real scale (0.0);
    scale = -invariants[0] / 3.0;
    copyCofactorF.Scale ( scale );
    firstTerm += copyCofactorF;

    //Computation trace of the isochoric C
    Real trCiso (0.0);
    trCiso = std::pow (invariants[3], - (2.0 / 3.0) ) * invariants[0];

    Real coef ( 0.0 );
    coef = alpha * std::pow (invariants[3], - (2.0 / 3.0) ) * std::exp ( gamma * ( trCiso - 3 ) );
    firstTerm.Scale ( coef );

    //Computing the second term (volumetric part) J*(bulk/2)(J-1+(1/J)*ln(J))F^{-T}
    Epetra_SerialDenseMatrix secondTerm (cofactorF);
    Real secCoef (0);
    secCoef = invariants[3] * (bulk / 2.0) * (invariants[3] - 1 + (1.0 / invariants[3]) * std::log (invariants[3]) );

    secondTerm.Scale ( secCoef );

    firstPiola += firstTerm;
    firstPiola += secondTerm;

}

template <typename Mesh>
inline StructuralConstitutiveLaw<Mesh>* createVenantKirchhoffMaterialNonLinearPenalized()
{
    return new VenantKirchhoffMaterialNonLinearPenalized<Mesh >();
}

namespace
{
static bool registerVKP = StructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ( "nonLinearVenantKirchhoffPenalized", &createVenantKirchhoffMaterialNonLinearPenalized<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __VENANTKIRCHHOFFPENALIZED_H */
