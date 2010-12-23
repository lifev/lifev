/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2010-11-29

  Copyright (C) 2010 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file PreconditionerPCD.hpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2010-11-29
 */


#ifndef _PRECONDITIONERPCD_HPP_
#define _PRECONDITIONERPCD_HPP_

#include <boost/shared_ptr.hpp>

#include <life/lifefilters/GetPot.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>
#include <lifemc/lifealg/ComposedPreconditioner.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <life/lifesolver/ADRAssembler.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>

namespace LifeV {

//! PreconditionerPCD
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerPCD class provides the PCD block preconditioner
 */
class PreconditionerPCD:
        public ComposedPreconditioner
{
public:

    /** @name Public Types
     */
    //@{
    typedef RegionMesh3D<LinearTetra>               mesh_type;
    typedef MapEpetra                               map_type;
    typedef MatrixBlock<Real>                       matrix_type;
    typedef MatrixEpetra<Real>                      parent_matrix_type;
    typedef Epetra_FECrsMatrix                      src_matrix_type;
    typedef VectorEpetra                            vector_type;
    typedef boost::shared_ptr<vector_type>          vector_ptr;

    typedef Preconditioner                          super;

    typedef ComposedOperator<Preconditioner>        prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>        prec_type;

    typedef super::operator_raw_type                operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type>    operator_type;

    typedef boost::shared_ptr<FESpace<mesh_type,map_type> >  FESpace_ptr;

    typedef Teuchos::ParameterList                  list_type;
    //@}


    //! @name Constructors, destructor
    //@{
    //! default constructor.
    PreconditionerPCD();

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner(operator_type& A);

    //! default destructor
    ~PreconditionerPCD();

    //@}

    //! @name  Methods
    //@{
    void createList( list_type&         list,
                     const GetPot&      dataFile,
                     const std::string& section,
                     const std::string& subSection );

    static void createPCDList( list_type&         list,
                               const GetPot&      dataFile,
                               const std::string& section,
                               const std::string& subSection = "PCD" );

    //! Return an estimation of the conditionement number of the preconditioner
    double      Condest ();

    //! Update the vector beta of the convective term in Fp
    /*!
        This method updates the value of beta.
        @param beta New vector beta to be used to built the convective term
     */
    void updateBeta(const vector_type& beta);


    //! Build the preconditioner
    int         buildPreconditioner(operator_type& A);
    //@}

    //! @name  Get Methods
    //@{
    int         numBlocksRows() const;
    int         numBlocksCols() const;

    //! Return the name of the preconditioner to be used in the factory
    std::string precType(){return M_precType;}
    //@}

    //! @name  Set Methods
    //@{
    //! Setter using GetPot
    /*!
        This method use GetPot to load data from a file and then set
        the preconditioner.
        @param dataFile is a GetPot dataFile
        @param section is the section containing the data
     */
    void setDataFromGetPot ( const GetPot&      dataFile,
                            const std::string& section );

    //! Setter for the FESpace
    /*!
        This method set the pointer for the FESpaces needed
        for the construction of the operators Ap, Fp and Mp.
        @param uFESpace Boost::shared_ptr on the FESpace for the velocity
        @param pFESpace Boost::shared_ptr on the FESpace for the pressure
     */
    void setFESpace(FESpace_ptr uFESpace,FESpace_ptr pFESpace);

    //! Setter for the timestep
    /*!
        This method set the timestep used to compute Fp.
        @param timestep Timestep used to compute the solution of the Navier-Stokes equations
     */
    void setTimestep(const Real& timestep);

    //! Setter for the viscosity
    /*!
        This method set the viscosity used to compute Fp.
        @param viscosity Viscosity used to compute the solution of the Navier-Stokes equations
     */
    void setViscosity(const Real& viscosity);

    //! Setter for the densitz
    /*!
        This method set the density used to compute Fp.
        @param density Density used to compute the solution of the Navier-Stokes equations
     */
    void setDensity(const Real& density);

    //@}

protected:

    std::string M_precType;
    int         M_velocityBlockSize;
    int         M_pressureBlockSize;
    FESpace_ptr M_uFESpace;
    FESpace_ptr M_pFESpace;

    Real        M_timestep;
    Real        M_viscosity;
    Real        M_density;
    vector_ptr  M_beta;

    ADRAssembler<mesh_type,matrix_type,vector_type> M_adrPressureAssembler;
    ADRAssembler<mesh_type,matrix_type,vector_type> M_adrVelocityAssembler;

};

inline Preconditioner* createPCD(){ return new PreconditionerPCD(); }
namespace
{
	static bool registerPCD = PRECFactory::instance().registerProduct( "PCD", &createPCD );
}

} // namespace LifeV

#endif
