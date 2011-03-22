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
 *  @brief This file contains an abstract class to implement different kinds of materials for structural dynamic (St. Venant-Kirchhoff materials right now )
 *
 *  @version 1.0
 *  @date 01-01-2010
 *  @author Paolo Tricerri
 *
 *  @maintainer  Paolo Tricerri <paolo.tricerri@epfl.ch>
 */

#ifndef _STRUCTURALMATERIAL_H_
#define _STRUCTURALMATERIAL_H_ 1

#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <boost/scoped_ptr.hpp>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>

#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <life/lifearray/MatrixElemental.hpp>
#include <life/lifearray/VectorElemental.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>

#include <life/lifefem/AssemblyElemental.hpp>
#include <life/lifefem/Assembly.hpp>
#include <life/lifefem/FESpace.hpp>

#include <life/lifecore/Displayer.hpp>
#include <life/lifecore/Factory.hpp>
#include <life/lifecore/FactorySingleton.hpp>

#include <life/lifealg/SolverAztecOO.hpp>

#include <life/lifesolver/VenantKirchhoffElasticData.hpp>


namespace LifeV
{


/*!
  \class StructuralMaterial
  \brief
  This class is an abstract class to define different type of models for the arterial wall.
  This class has just pure virtual methods. They are implemented in the specific class for one material
*/

template <typename Mesh>
class StructuralMaterial
{
public:

  //!@name Type definitions
  //@{

  typedef VenantKirchhoffElasticData             data_Type;

  typedef typename LifeV::SolverAztecOO          solver_Type;

  typedef typename solver_Type::matrix_type      matrix_Type;
  typedef boost::shared_ptr<matrix_Type>         matrixPtr_Type;
  typedef typename solver_Type::vector_type      vector_Type;
  typedef boost::shared_ptr<vector_Type>         vectorPtr_Type;

  typedef typename boost::shared_ptr<data_Type>    dataPtr_Type;
  typedef typename boost::scoped_ptr<Displayer>    displayerPtr_Type;

  typedef FactorySingleton<Factory<StructuralMaterial<Mesh>,std::string> >  StructureMaterialFactory;
  //@}


  //! @name Constructor &  Deconstructor
  //@{

  StructuralMaterial();

  virtual ~StructuralMaterial() {}

  //@}

  //!@name Methods
  //@{

  //! Setup the created object of the class StructuralMaterial
  /*!
    \param dFespace: the FiniteElement Space
    \param monolithicMap: the MapEpetra
    \param offset: the offset parameter used assembling the matrices
  */
    virtual void setup( const boost::shared_ptr< FESpace<Mesh, MapEpetra> >& dFESpace,
                        const boost::shared_ptr<const MapEpetra>&   monolithicMap,
                        const UInt offset
                        ) = 0;

  //! Computes the linear part of the stiffness matrix StructuralSolver::buildSystem
  /*!
    \param dataMaterial the class with Material properties data
  */
  virtual  void computeLinearStiffMatrix( dataPtr_Type& dataMaterial ) = 0;

  //! Updates the Jacobian matrix in StructuralSolver::updateJacobian
  /*!
    \param disp: solution at the k-th iteration of NonLinearRichardson Method
    \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
    \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
  */
    virtual  void updateJacobianMatrix( const vector_Type& disp, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer ) = 0;

  //! Updates the nonlinear terms in the Jacobian matrix in StructuralSolver::updateJacobian
  /*!
    \param stiff: stiffness matrix provided from outside
    \param disp: solution at the k-th iteration of NonLinearRichardson Method
    \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
    \param displayer: a pointer to the Dysplaier member in the StructuralSolver class    
  */
  virtual  void updateNonLinearJacobianMatrix( matrixPtr_Type& /*stiff*/, const vector_Type& /*disp*/, const dataPtr_Type& /*dataMaterial*/, const displayerPtr_Type& /*displayer*/ ) = 0;

    //! Computes the new Stiffness matrix in StructuralSolver given a certain displacement field. This function is used both in StructuralSolver::evalResidual and in 
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
    virtual  void computeMatrix( const vector_Type& sol, Real factor, const dataPtr_Type& dataMaterial, const displayerPtr_Type& displayer ) = 0;

    //! Computes the nonlinear part of Stiffness matrix in StructuralSolver given a certain displacement field. This function is used both in StructuralSolver::evalResidual and in 
    //! StructuralSolver::updateSystem since the matrix is the expression of the matrix is the same. This is virtual and not pure virtual since in the linear St. Venant-Kirchhoff law it is not needed.
    /*!
      \param sol:  the solution vector
      \param factor: scaling factor used in FSI
      \param dataMaterial: a pointer to the dataType member in StructuralSolver class to get the material coefficients (e.g. Young modulus, Poisson ratio..)
      \param displayer: a pointer to the Dysplaier member in the StructuralSolver class
    */
  virtual  void computeNonLinearMatrix( matrixPtr_Type& /*stiff*/, const vector_Type& /*sol*/, Real /*factor*/, const dataPtr_Type& /*dataMaterial*/, const displayerPtr_Type& /*displayer*/ ){};


  //! Output of the class
  /*!
    \param fileNamelinearStiff the filename where to apply the spy method for the linear part of the Stiffness matrix
    \param fileNameStiff the filename where to apply the spy method for the Stiffness matrix
  */
  void showMe( std::string const& fileNamelinearStiff, std::string const& fileNameStiff );

  //! @name Set Methods
  //@{

    //No set Methods

  //@}


  //! @name Get Methods
  //@{

  //! Getters
  //! Get the Epetramap
  MapEpetra   const& map()     const { return *M_localMap; }

  //! Get the Stiffness matrix
  matrixPtr_Type const stiff()    const {return M_stiff; }

  //! Get the Stiffness matrix
  matrixPtr_Type const linearStiff()    const {return M_linearStiff; }

  //! Get the FESpace object
  FESpace<Mesh, MapEpetra>& dFESpace()  {return M_FESpace;}

  //@}

protected:

  //!Protected Members

  boost::shared_ptr<FESpace<Mesh, MapEpetra> >   M_FESpace;

  boost::shared_ptr<const MapEpetra>             M_localMap;

  //! Elementary matrices
  boost::scoped_ptr<MatrixElemental>             M_elmatK;

  //! Matrix Knl: stiffness (linear + nonlinear)
  matrixPtr_Type                                 M_stiff;

  //! Matrix Kl: stiffness linear
  matrixPtr_Type                                 M_linearStiff;

  //! The Offset parameter
  UInt                                           M_offset;

};

//====================================
// Constructor
//=====================================

template <typename Mesh>
StructuralMaterial<Mesh>::StructuralMaterial( ):
  M_FESpace                    ( ),
  M_localMap                   ( ),
  M_elmatK                     ( ),
  M_stiff                      ( ),
  M_linearStiff                ( ),
  M_offset                     ( 0 )
{
  std::cout << "I am in the constructor of StructuralMaterial" << std::endl;
}

template <typename Mesh>
void
StructuralMaterial<Mesh>::showMe( std::string const& fileNamelinearStiff,
				  std::string const& fileNameStiff
				)
{
  this->M_linearStiff(fileNamelinearStiff);
  this->M_stiff(fileNameStiff);
}

}

#endif /*_STRUCTURALMATERIAL_H*/
