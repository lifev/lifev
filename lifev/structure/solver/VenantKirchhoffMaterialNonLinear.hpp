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

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>
#include <lifev/structure/solver/VenantKirchhoffMaterialLinear.hpp>

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

    typedef StructuralConstitutiveLawData            data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;
    typedef typename super::vectorPtr_Type              vectorPtr_Type;
    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapIterator_Type mapIterator_Type;

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
               const UInt offset, const dataPtr_Type& dataMaterial,
               const displayerPtr_Type& displayer
               );

    //! Compute the linear part Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiffMatrix( dataPtr_Type& dataMaterial,
                                   const mapMarkerVolumesPtr_Type mapsMarkerVolumes );

    //! Updates the Jacobian matrix
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix( const vector_Type& disp,
                               const dataPtr_Type& dataMaterial,
                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                               const displayerPtr_Type& displayer);
    //! Updates the nonlinear terms in the Jacobian matrix
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms(  matrixPtr_Type& stiff,
                                        const vector_Type& disp,
                                        const dataPtr_Type& dataMaterial,
                                        const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                        const displayerPtr_Type& displayer);


    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeStiffness( const vector_Type& sol, Real factor,
                           const dataPtr_Type& dataMaterial,
                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                           const displayerPtr_Type& displayer );

    //! Computes the nonlinear part of Stiffness matrix in StructuralSolver given a certain
    //! displacement field. This function is used both in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    //! This is virtual and not pure virtual since in the linear St. Venant-Kirchhoff law it is not needed.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void computeNonLinearMatrix( matrixPtr_Type& stiff, const vector_Type& sol, Real factor,
                                 const dataPtr_Type& dataMaterial, const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                 const displayerPtr_Type& displayer );

    void computeKinematicsVariables( const VectorElemental& /*dk_loc*/ ){}

    //@}

protected:
    //KNMKPtr_Type					M_gradientLocalDisplacement;
    boost::shared_ptr<boost::multi_array<Real, 3> >	M_gradientLocalDisplacement;

};

template <typename Mesh>
VenantKirchhoffMaterialNonLinear<Mesh>::VenantKirchhoffMaterialNonLinear():
    super(),
    M_gradientLocalDisplacement()
{
}

template <typename Mesh>
VenantKirchhoffMaterialNonLinear<Mesh>::~VenantKirchhoffMaterialNonLinear()
{}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::setup(const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                                   const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                                   const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer
                                                   )
{
    super::setup(dFESpace,monolithicMap,offset, dataMaterial, displayer);
    this->M_stiff.reset               (new matrix_Type(*this->M_localMap));
    //M_gradientLocalDisplacement.reset ( new KNMK_Type( nDimensions, nDimensions, dFESpace->fe().nbQuadPt() ) );
    M_gradientLocalDisplacement.reset ( new boost::multi_array<Real, 3>(boost::extents[nDimensions][nDimensions][dFESpace->fe().nbQuadPt()]) );

    //    M_gradientLocalDisplacement.
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeLinearStiffMatrix(dataPtr_Type& dataMaterial,
                                                                      const mapMarkerVolumesPtr_Type mapsMarkerVolumes)
{
    super::computeLinearStiffMatrix(dataMaterial, mapsMarkerVolumes);
}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::updateJacobianMatrix(const vector_Type& disp,
                                                                  const dataPtr_Type& dataMaterial,
                                                                  const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                  const displayerPtr_Type& displayer)
{
    this->M_jacobian.reset(new matrix_Type(*this->M_localMap));

    *this->M_jacobian += *this->M_linearStiff;

    displayer->leaderPrint(" \n*********************************\n  ");
    updateNonLinearJacobianTerms(this->M_jacobian,disp,dataMaterial,mapsMarkerVolumes,displayer);
    displayer->leaderPrint(" \n*********************************\n  ");
    this->M_jacobian->globalAssemble();
}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::updateNonLinearJacobianTerms( matrixPtr_Type& jacobian,
                                                                           const  vector_Type& disp,
                                                                           const dataPtr_Type& dataMaterial,
                                                                           const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                           const displayerPtr_Type& displayer )
{
    displayer->leaderPrint("   NonLin S-  Updating non linear terms in the Jacobian Matrix (in updateJacobian)");
    // std::cout << std::endl;

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    VectorElemental dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );

    vector_Type dRep(disp, Repeated);

    //! Number of displacement components
    UInt nc = nDimensions;

    mapIterator_Type it;

    for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
	{

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        Real mu = dataMaterial->mu(marker);
        Real lambda = dataMaterial->lambda(marker);

        for ( UInt j(0); j < it->second.size(); j++ )
	    {
            this->M_FESpace->fe().updateFirstDerivQuadPt( *(it->second[j]) );

            this->M_elmatK->zero();


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

            //Reset the local gradient of the Displacement
            //this->M_gradientLocalDisplacement.reset (new KNMK_Type(nDimensions , nDimensions,this->M_FESpace->fe().nbQuadPt() ) );
            //*M_gradientLocalDisplacement *= 0.0;
            AssemblyElementalStructure::computeGradientLocalDisplacement(*M_gradientLocalDisplacement, dk_loc, this->M_FESpace->fe());

            //  3):  \lambda * ( \tr { [\grad d^k]^T \grad \delta d }, \div v  )
            AssemblyElementalStructure::stiff_derdiv( lambda, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );
            //  4):  \mu * ( [\grad \delta d]^T \grad d^k + [\grad d^k]^T \grad \delta d : \grad v  )
            AssemblyElementalStructure::stiff_dergrad( mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            // the sum of these terms is the Jacobian of the divgrad term
            // 5):  \lambda * ( (\div u_k) \grad \delta u : \grad v  )
            AssemblyElementalStructure::stiff_divgrad(  lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

            //  \lambda * ( (\div u) \grad u_k : \grad v  )
            AssemblyElementalStructure::stiff_divgrad_2(  lambda, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            // the sum of these terms is the Jacobian of the gradgrad term
            // 6): 1/2 * \lambda * ( \grad u_k : \grad  u_k) *( \grad \delta u : \grad v  )
            AssemblyElementalStructure::stiff_gradgrad(   0.5 * lambda, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

            //\lambda * ( \grad u_k : \grad \delta u) *( \grad u_k : \grad v  )
            AssemblyElementalStructure::stiff_gradgrad_2(  lambda, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            // the sum of these terms is he jacobian of the stiff_dergrad_gradbis term
            // 7A) : \mu *  ( \grad u^k \grad \delta u : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis(  mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            //  \mu *  ( \grad \delta u \grad u^k : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis_2(  mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            //  the sum of these terms is he jacobian of the stiff_dergrad_gradbis_Tr term
            // 7B) :  \mu *  ( \grad u^k [\grad \delta u]^T : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis_Tr(  mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            // \mu *  ( \grad \delta u [\grad u^k]^T : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis_Tr_2( mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            //   the sum of these terms is he jacobian of the stiff_gradgradTr_gradbis term
            // 8) :   \mu * (  \grad d^k [\grad d^k]^T \grad \delta d : \grad v  )
            AssemblyElementalStructure::stiff_gradgradTr_gradbis(  mu, dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

            //  \mu * (  \grad d^k [\grad \delta d]^T \grad d^k : \grad v  )
            AssemblyElementalStructure::stiff_gradgradTr_gradbis_2( mu, *M_gradientLocalDisplacement, *this->M_elmatK, this->M_FESpace->fe() );

            //  \mu * (  \grad \delta u [\grad u^k]^T \grad u^k : \grad v  )
            AssemblyElementalStructure::stiff_gradgradTr_gradbis_3(  mu , dk_loc, *this->M_elmatK, this->M_FESpace->fe() );

            // assembling
            for ( UInt ic = 0; ic < nc; ++ic )
                for ( UInt jc = 0; jc < nc; jc++ )
                    assembleMatrix( *jacobian,
                                    *this->M_elmatK,
                                    this->M_FESpace->fe(),
                                    this->M_FESpace->dof(),
                                    ic, jc,
                                    this->M_offset +  ic*totalDof, this->M_offset + jc*totalDof );
	    }
	}

}


template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeStiffness( const vector_Type& disp,
                                                               Real factor,
                                                               const dataPtr_Type& dataMaterial,
                                                               const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                               const displayerPtr_Type& displayer )
{

    this->M_stiff.reset(new matrix_Type(*this->M_localMap));

    displayer->leaderPrint(" \n*********************************\n  ");
    computeNonLinearMatrix(this->M_stiff,disp,factor,dataMaterial,mapsMarkerVolumes,displayer);
    displayer->leaderPrint(" \n*********************************\n  ");

    *this->M_stiff += *this->M_linearStiff;

    this->M_stiff->globalAssemble();
}

template <typename Mesh>
void VenantKirchhoffMaterialNonLinear<Mesh>::computeNonLinearMatrix(matrixPtr_Type& stiff,
                                                                    const vector_Type& sol,
                                                                    Real /*factor*/,
                                                                    const dataPtr_Type& dataMaterial,
                                                                    const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                    const displayerPtr_Type& displayer)
{
    displayer->leaderPrint("   NonLin S-  Updating non linear terms in the Stiffness Matrix (Venant-Kirchhoff)");

    UInt totalDof   = this->M_FESpace->dof().numTotalDof();
    UInt dim = this->M_FESpace->dim();

    VectorElemental dk_loc( this->M_FESpace->fe().nbFEDof(), nDimensions );
    vector_Type disp(sol);

    vector_Type dRep(disp, Repeated);


    mapIterator_Type it;

    for( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    {

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        Real mu = dataMaterial->mu(marker);
        Real lambda = dataMaterial->lambda(marker);

        for ( UInt j(0); j < it->second.size(); j++ )
        {
            this->M_FESpace->fe().updateFirstDerivQuadPt( *(it->second[j]) );


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
            //*M_gradientLocalDisplacement *= 0.0;
            AssemblyElementalStructure::computeGradientLocalDisplacement(*M_gradientLocalDisplacement, dk_loc, this->M_FESpace->fe());

            // 3) 1/2 * \lambda  ( \tr { [\grad d^k]^T \grad d }, \div v  )
            AssemblyElementalStructure::stiff_derdiv( 0.5 * lambda , *M_gradientLocalDisplacement, *this->M_elmatK,  this->M_FESpace->fe() );

            //4)  \mu *( [\grad d^k]^T \grad d : \grad v  )
            AssemblyElementalStructure::stiff_dergradbis( mu , *M_gradientLocalDisplacement, *this->M_elmatK,  this->M_FESpace->fe() );

            //  5): \lambda * (div u_k) \grad d : \grad v
            AssemblyElementalStructure::stiff_divgrad( lambda, dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

            // 6) 1/2  * \lambda * ( \grad u_k : \grad u_k) *( \grad u : \grad v  )
            AssemblyElementalStructure::stiff_gradgrad( 0.5 * lambda , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

            // 7A) \mu *  ( \grad d^k \grad d : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis( mu , *M_gradientLocalDisplacement, *this->M_elmatK,  this->M_FESpace->fe() );
            // 7B) \mu *  ( \grad d^k [\grad d]^T : \grad v  )
            AssemblyElementalStructure::stiff_dergrad_gradbis_Tr( mu , *M_gradientLocalDisplacement, *this->M_elmatK,  this->M_FESpace->fe() );

            // 8) // \mu *  (  \grad d^k [\grad d^k]^T \grad d : \grad v  )
            AssemblyElementalStructure::stiff_gradgradTr_gradbis( mu , dk_loc, *this->M_elmatK,  this->M_FESpace->fe() );

            for ( UInt ic = 0; ic < nDimensions; ++ic )
            {
                // stiff is the nonlinear matrix of the bilinear form
                for ( UInt jc = 0; jc < nDimensions; jc++ )

                    assembleMatrix( *stiff,
                                    *this->M_elmatK,
                                    this->M_FESpace->fe(),
                                    this->M_FESpace->dof(),
                                    ic, jc,
                                    this->M_offset +  ic*totalDof, this->M_offset +  jc*totalDof);
            }

        }
    }
}

template <typename Mesh>
inline StructuralConstitutiveLaw<Mesh>* createVenantKirchhoffNonLinear() { return new VenantKirchhoffMaterialNonLinear<Mesh >(); }
namespace
{
static bool registerVKNL = StructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct( "nonLinearVenantKirchhoff", &createVenantKirchhoffNonLinear<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __NONLINVENANTKIRCHHOFF_H */
