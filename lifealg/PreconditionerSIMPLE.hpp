/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
       Date: 2011-01-24

  Copyright (C) 2011 EPFL

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
   \file PreconditionerSIMPLE.hpp
   \author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
   \date 2011-01-24
 */


#ifndef _PRECONDITIONERSIMPLE_HPP_
#define _PRECONDITIONERSIMPLE_HPP_

#include <boost/shared_ptr.hpp>

#include <life/lifefilters/GetPot.hpp>
#include <life/lifearray/MatrixEpetra.hpp>
#include <life/lifearray/VectorEpetra.hpp>
#include <lifemc/lifealg/PreconditionerComposition.hpp>
#include <lifemc/lifearray/MatrixBlock.hpp>
#include <life/lifesolver/ADRAssembler.hpp>
#include <life/lifefem/FESpace.hpp>
#include <life/lifemesh/RegionMesh3D.hpp>

namespace LifeV {

//! PreconditionerSIMPLE
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerSIMPLE class provides the SIMPLE block preconditioner
 */
class PreconditionerSIMPLE:
        public PreconditionerComposition
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
    typedef boost::shared_ptr<super>                super_PtrType;

    typedef ComposedOperator<Preconditioner>        prec_raw_type;
    typedef boost::shared_ptr<prec_raw_type>        prec_type;

    typedef super::operator_raw_type                operator_raw_type;
    typedef boost::shared_ptr<operator_raw_type>    operator_type;

    typedef boost::shared_ptr<FESpace<mesh_type,map_type> >  FESpace_ptr;

    typedef Teuchos::ParameterList                  list_Type;
    //@}


    //! @name Constructors, destructor
    //@{
    //! default constructor.
    PreconditionerSIMPLE(const  boost::shared_ptr<Epetra_Comm>& comm= boost::shared_ptr<Epetra_Comm>());

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner(operator_type& A);

    //! default destructor
    ~PreconditionerSIMPLE();

    //@}

    //! @name  Methods
    //@{
    void createParametersList( list_Type&         list,
                               const GetPot&      dataFile,
                               const std::string& section,
                               const std::string& subSection );

    static void createSIMPLEList( list_Type&         list,
                                  const GetPot&      dataFile,
                                  const std::string& section,
                                  const std::string& subSection = "SIMPLE" );

    //! Return an estimation of the conditionement number of the preconditioner
    double      condest ();

    //! Build the preconditioner
    int buildPreconditioner(operator_type& A);

    //@}

    //! @name  Get Methods
    //@{
    int numBlocksRows() const;
    int numBlocksColumns() const;
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

    //! Setter for the damping factor
    /*!
        This method set damping factor used to build the preconditioner
        @param dampingFactor Damping factor
    */
    void setDampingFactor(const Real& dampingFactor);

    //@}

protected:

    int         M_velocityBlockSize;
    int         M_pressureBlockSize;

    Real        M_dampingFactor;

    string      M_SIMPLEType;

    // todo: Remove the member dataFile (bad programmation)
    GetPot      M_dataFile;
    string      M_fluidPrec;
    string      M_fluidDataSection;
    string      M_schurPrec;
    string      M_schurDataSection;

private:
    PreconditionerSIMPLE(const PreconditionerSIMPLE& P):
        PreconditionerComposition(P.M_comm){}
    PreconditionerSIMPLE(const boost::shared_ptr<PreconditionerSIMPLE>& /*P*/){}

};

inline Preconditioner* createSIMPLE(){ return new PreconditionerSIMPLE(); }
namespace
{
	static bool registerSIMPLE = PRECFactory::instance().registerProduct( "SIMPLE", &createSIMPLE );
}

} // namespace LifeV

#endif
