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

#include <lifev/structure/solver/StructuralConstitutiveLaw.hpp>


namespace LifeV
{
template <typename Mesh>
class VenantKirchhoffMaterialLinear :
    public StructuralConstitutiveLaw<Mesh>
{
    //!@name Type definitions
    //@{

public:
    typedef StructuralConstitutiveLaw<Mesh>                 super;

    typedef typename super::data_Type                         data_Type;

    typedef typename super::vector_Type              vector_Type;
    typedef typename super::matrix_Type              matrix_Type;

    typedef typename super::matrixPtr_Type           matrixPtr_Type;
    typedef typename super::dataPtr_Type             dataPtr_Type;
    typedef typename super::displayerPtr_Type        displayerPtr_Type;
    typedef typename super::vectorPtr_Type           vectorPtr_Type;

    typedef typename super::mapMarkerVolumesPtr_Type mapMarkerVolumesPtr_Type;
    typedef typename super::mapMarkerVolumes_Type mapMarkerVolumes_Type;
    typedef typename mapMarkerVolumes_Type::const_iterator mapIterator_Type;


    //@}

    //! @name Constructor &  Destructor
    //@{

    VenantKirchhoffMaterialLinear();

    virtual  ~VenantKirchhoffMaterialLinear();
    //@}

    //!@name Methods
    //@{

    //! Setup the created object of the class StructuralConstitutiveLaw
    /*!
      \param dFespace: the FiniteElement Space
      \param monolithicMap: the MapEpetra
      \param offset: the offset parameter used assembling the matrices
    */
    void setup (const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer
               );


    //! Compute the Stiffness matrix in StructuralSolver::buildSystem()
    /*!
      \param dataMaterial the class with Material properties data
    */
    void computeLinearStiff ( dataPtr_Type& dataMaterial, const mapMarkerVolumesPtr_Type mapsMarkerVolumes );

    //! Updates the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateJacobianMatrix ( const vector_Type& disp,
                                const dataPtr_Type& dataMaterial,
                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                const displayerPtr_Type& displayer);

    //! Updates the nonlinear terms in the Jacobian matrix in StructualSolver::updateJacobian
    /*!
      \param stiff: stiffness matrix provided from outside
      \param disp: solution at the k-th iteration of NonLinearRichardson Method
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get
                           the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    void updateNonLinearJacobianTerms (  matrixPtr_Type& jacobian,
                                         const vector_Type& /*disp*/,
                                         const dataPtr_Type& /*dataMaterial*/,
                                         const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                         const displayerPtr_Type& /*displayer*/);

    //! Interface method to compute the new Stiffness matrix in StructuralSolver::evalResidual and in
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the
        material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */

    void computeStiffness ( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial,
                            const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                            const displayerPtr_Type& displayer );

    void computeKinematicsVariables ( const VectorElemental& /*dk_loc*/ ) {}

    //! ShowMe method of the class (saved on a file the two matrices)
    void showMe ( std::string const& fileNameStiff,
                  std::string const& fileNameJacobian);

    //@}

    //! @name Get Methods
    //@{

    //! Get the linear part of the matrix
    matrixPtr_Type const linearStiff() const
    {
        return M_linearStiff;
    }

    //! Get the Stiffness matrix
    matrixPtr_Type const stiffMatrix() const
    {
        return M_stiff;
    }

    //! Get the Stiffness vector
    vectorPtr_Type const stiffVector() const
    {
        vectorPtr_Type zero ( new vector_Type() );
        return zero;
    }

    void apply ( const vector_Type& sol, vector_Type& res,
                 const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/)
    {
        res += *M_stiff * sol;
    }

    //@}

protected:
    //! Protected members

    //! Elementary matrices
    boost::scoped_ptr<MatrixElemental>             M_elmatK;

    //! Matrix Kl: stiffness linear
    matrixPtr_Type                                 M_linearStiff;

    //! Matrix Kl: stiffness linear
    matrixPtr_Type                                 M_stiff;

};

template <typename Mesh>
VenantKirchhoffMaterialLinear<Mesh>::VenantKirchhoffMaterialLinear() :
    super             ( ),
    M_elmatK                     ( ),
    M_linearStiff                ( ),
    M_stiff                      ( )
{
}

template <typename Mesh>
VenantKirchhoffMaterialLinear<Mesh>::~VenantKirchhoffMaterialLinear()
{}


template <typename Mesh>
void
VenantKirchhoffMaterialLinear<Mesh>::setup (const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                                            const boost::shared_ptr<const MapEpetra>&  monolithicMap,
                                            const UInt offset, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer
                                           )
{
    this->M_displayer = displayer;
    this->M_dataMaterial  = dataMaterial;

    //    std::cout<<"I am setting up the Material "<<std::endl;

    this->M_FESpace                       = dFESpace;
    this->M_elmatK.reset                  (new MatrixElemental ( this->M_FESpace->fe().nbFEDof(), nDimensions, nDimensions ) );
    this->M_localMap                      = monolithicMap;
    this->M_linearStiff.reset             (new matrix_Type (*this->M_localMap) );
    this->M_offset                        = offset;
}

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::computeLinearStiff (dataPtr_Type& dataMaterial,
                                                              const mapMarkerVolumesPtr_Type mapsMarkerVolumes)
{
    //  std::cout<<"compute LinearStiff Matrix start\n";

    UInt totalDof = this->M_FESpace->dof().numTotalDof();
    // Number of displacement components
    UInt nc = nDimensions;

    //Compute the linear part of the Stiffness Matrix.
    //In the case of Linear Material it is the Stiffness Matrix.
    //In the case of NonLinear Materials it must be added of the non linear part.

    mapIterator_Type it;

    for ( it = (*mapsMarkerVolumes).begin(); it != (*mapsMarkerVolumes).end(); it++ )
    {

        //Given the marker pointed by the iterator, let's extract the material parameters
        UInt marker = it->first;

        Real mu = dataMaterial->mu (marker);
        Real lambda = dataMaterial->lambda (marker);

        //Given the parameters I loop over the volumes with that marker
        for ( UInt j (0); j < it->second.size(); j++ )
        {
            this->M_FESpace->fe().updateFirstDerivQuadPt ( * (it->second[j]) );

            this->M_elmatK->zero();

            //These methods are implemented in AssemblyElemental.cpp
            //They have been kept in AssemblyElemental in order to avoid repetitions
            stiff_strain ( 2 * mu, *this->M_elmatK, this->M_FESpace->fe() ); // here in the previous version was 1. (instead of 2.)
            stiff_div   ( lambda, *this->M_elmatK, this->M_FESpace->fe() );// here in the previous version was 0.5 (instead of 1.)

            //this->M_elmatK->showMe();

            // assembling
            for ( UInt ic = 0; ic < nc; ic++ )
            {
                for ( UInt jc = 0; jc < nc; jc++ )
                {
                    assembleMatrix ( *this->M_linearStiff,
                                     *this->M_elmatK,
                                     this->M_FESpace->fe(),
                                     this->M_FESpace->fe(),
                                     this->M_FESpace->dof(),
                                     this->M_FESpace->dof(),
                                     ic,  jc,
                                     this->M_offset + ic * totalDof, this->M_offset + jc * totalDof );

                }
            }


        }

    }

    this->M_linearStiff->globalAssemble();

    //Initialization of the pointer M_stiff to what is pointed by M_linearStiff
    this->M_stiff = this->M_linearStiff;
    //   std::cout<<"compute LinearStiff Matrix end\n";
    this->M_jacobian = this->M_linearStiff;
}


template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::updateJacobianMatrix (const vector_Type& disp,
                                                                const dataPtr_Type& dataMaterial,
                                                                const mapMarkerVolumesPtr_Type mapsMarkerVolumes,
                                                                const displayerPtr_Type& displayer)
{
    //displayer->leaderPrint(" \n*********************************\n  ");
    displayer->leaderPrint ("  S-  Updating the Jacobian Matrix (constant, Linear Elastic)\n");
    //displayer->leaderPrint(" \n*********************************\n  ");

    //displayer->leaderPrint(" \n*********************************\n  ");
    updateNonLinearJacobianTerms (this->M_jacobian, disp, dataMaterial, mapsMarkerVolumes, displayer);
    //displayer->leaderPrint(" \n*********************************\n  ");

}

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::updateNonLinearJacobianTerms ( matrixPtr_Type& /*jacobian*/,
                                                                         const  vector_Type& /*disp*/,
                                                                         const dataPtr_Type& /*dataMaterial*/,
                                                                         const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                                         const displayerPtr_Type& /*displayer*/ )
{
    //    this->M_stiff->globalAssemble();
    //  displayer->leaderPrint("  S- Doing nothing - Updating non linear terms in Jacobian Matrix (constant, Linear Elastic)\n");
}

template <typename Mesh>
void VenantKirchhoffMaterialLinear<Mesh>::computeStiffness ( const vector_Type& /*disp*/,
                                                             Real /*factor*/,
                                                             const dataPtr_Type& /*dataMaterial*/,
                                                             const mapMarkerVolumesPtr_Type /*mapsMarkerVolumes*/,
                                                             const displayerPtr_Type& /*displayer*/ )
{
    //displayer->leaderPrint(" \n*********************************\n  ");
    //displayer->leaderPrint("  S- Using the the Stiffness Matrix (constant, Linear Elastic)");
    //displayer->leaderPrint(" \n*********************************\n  ");
}


template <typename Mesh>
void
VenantKirchhoffMaterialLinear<Mesh>::showMe ( std::string const& fileNameStiff,
                                              std::string const& fileNameJacobian
                                            )
{
    //This string is to save the linear part
    std::string fileNamelinearStiff =  fileNameStiff;
    fileNamelinearStiff += "linear";

    this->M_linearStiff->spy (fileNamelinearStiff);
    this->M_stiff->spy (fileNameStiff);
    this->M_jacobian->spy (fileNameJacobian);
}

template <typename Mesh>
inline StructuralConstitutiveLaw<Mesh>* createVenantKirchhoffLinear()
{
    return new VenantKirchhoffMaterialLinear<Mesh >();
}
namespace
{
static bool registerVKL = StructuralConstitutiveLaw<LifeV::RegionMesh<LinearTetra> >::StructureMaterialFactory::instance().registerProduct ( "linearVenantKirchhoff", &createVenantKirchhoffLinear<LifeV::RegionMesh<LinearTetra> > );
}

} //Namespace LifeV

#endif /* __LINVENANTKIRCHHOFF_H */
