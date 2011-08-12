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
 *  @brief This file contains the definition for the St. Venant Kirchhoff linear material
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _LINVENANTKIRCHHOFFMATERIAL_H_
#define _LINVENANTKIRCHHOFFMATERIAL_H_

#include <life/lifesolver/StructuralMaterial.hpp>

namespace LifeV
{
template <typename Mesh>
class VenantKirchhoffMaterialLinear :
        public StructuralMaterial<Mesh>
{
 //!@name Type definitions
 //@{

  public:
    typedef StructuralMaterial<Mesh>                 super;

    typedef typename  super::data_Type                         data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;


 //@}

 //! @name Constructor &  Destructor
 //@{

    VenantKirchhoffMaterialLinear();

    virtual  ~VenantKirchhoffMaterialLinear();
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


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiffMatrix( dataPtr_Type& dataMaterial );

    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix( const vector_Type& disp,
                               const dataPtr_Type& dataMaterial,
                               const displayerPtr_Type& displayer);

    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianMatrix( matrixPtr_Type& stiff,
                                        const vector_Type& /*disp*/,
                                        const dataPtr_Type& /*dataMaterial*/,
                                        const displayerPtr_Type& /*displayer*/);


    //! Computes the new Stiffness matrix in StructuralSolver given a certain displacement field. This function is used both in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    virtual void computeMatrix( const vector_Type& /*sol*/,
                                Real /*factor*/,
                                const dataPtr_Type& /*dataMaterial*/,
                                const displayerPtr_Type& displayer);

    //@}

};

template <typename Mesh>
VenantKirchhoffMaterialLinear<Mesh>::VenantKirchhoffMaterialLinear():
    super()
{
}

template <typename Mesh>
VenantKirchhoffMaterialLinear<Mesh>::~VenantKirchhoffMaterialLinear()
{}


template <typename Mesh>
void
VenantKirchhoffMaterialLinear<Mesh>::setup(const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                           const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                           const UInt offset
				)
{
  std::cout<<"I am setting up the Material"<<std::endl;

  this->M_FESpace                       = dFESpace;
  this->M_elmatK.reset                        (new MatrixElemental( this->M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
  this->M_localMap                      = monolithicMap;
  this->M_linearStiff.reset             (new matrix_Type(*this->M_localMap));
  this->M_offset                        = offset;
}

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::computeLinearStiffMatrix(dataPtr_Type& dataMaterial)
{
    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    // Number of displacement components
    UInt nc = nDimensions;

    //Compute the linear part of the Stiffness Matrix.
    //In the case of Linear Material it is the Stiffness Matrix.
    //In the case of NonLinear Materials it must be added of the non linear part.

    for ( UInt i = 0; i < this->M_FESpace->mesh()->numVolumes(); i++ )
    {
        this->M_FESpace->fe().updateFirstDerivQuadPt( this->M_FESpace->mesh()->volumeList( i ) );

        this->M_elmatK->zero();

        UInt marker = this->M_FESpace->mesh()->volumeList( i ).marker();

	Real mu = dataMaterial->mu(marker);
	Real lambda = dataMaterial->lambda(marker);

    stiff_strain( mu, *this->M_elmatK, this->M_FESpace->fe() );// here in the previous version was 1. (instead of 2.)
    stiff_div   ( lambda, *this->M_elmatK, this->M_FESpace->fe() );// here in the previous version was 0.5 (instead of 1.)

        // assembling
        for ( UInt ic = 0; ic < nc; ic++ )
        {
            for ( UInt jc = 0; jc < nc; jc++ )
            {
                assembleMatrix( *this->M_linearStiff, *this->M_elmatK, this->M_FESpace->fe(), this->M_FESpace->fe(), this->M_FESpace->dof(), this->M_FESpace->dof(),  ic,  jc,  this->M_offset +ic*totalDof, this->M_offset + jc*totalDof );

            }
        }

    }

    this->M_linearStiff->globalAssemble();

    //Initialization of the pointer M_stiff to what is pointed by M_linearStiff
    this->M_stiff = this->M_linearStiff;
}


template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::updateJacobianMatrix(const vector_Type& disp,
                                                               const dataPtr_Type& dataMaterial,
                                                               const displayerPtr_Type& displayer)
{

    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
    displayer->leaderPrint("   Linear S-  Using the Stiffness Matrix (constant) in UpdateJacobian");
    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;

    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
    updateNonLinearJacobianMatrix(this->M_stiff,disp,dataMaterial,displayer);
    std::cout << "*********************************" << std::endl;
    std::cout << std::endl;
}

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::updateNonLinearJacobianMatrix( matrixPtr_Type& /*stiff*/,
                                                                         const  vector_Type& /*disp*/,
                                                                         const dataPtr_Type& /*dataMaterial*/,
                                                                         const displayerPtr_Type& displayer )
    {
      displayer->leaderPrint("   Linear S-  Doing nothing (updating non linear terms in the Jacobian Matrix (in updateJacobian)");
      std::cout << std::endl;
    }

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::computeMatrix(const vector_Type& /*disp*/,
							Real /*factor*/,
							const dataPtr_Type& /*dataMaterial*/,
							const displayerPtr_Type& displayer)
{
    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
    displayer->leaderPrint("   Linear S-  Using the linear part of the Stiffness Matrix\n");
    std::cout << std::endl;
    std::cout << "*********************************" << std::endl;
}


template <typename Mesh>
inline StructuralMaterial<Mesh>* createVenantKirchhoffLinear() { return new VenantKirchhoffMaterialLinear<Mesh >(); }
namespace
{
static bool registerVKL = StructuralMaterial<LifeV::RegionMesh3D<LinearTetra> >::StructureMaterialFactory::instance().registerProduct( "linearVenantKirchhoff", &createVenantKirchhoffLinear<LifeV::RegionMesh3D<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __LINVENANTKIRCHHOFF_H */
