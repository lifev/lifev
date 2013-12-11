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
    @brief PreconditionerYosida

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 13-10-2011
 */

#ifndef PRECONDITIONERYOSIDA_HPP
#define PRECONDITIONERYOSIDA_HPP 1

#include <boost/shared_ptr.hpp>

#include <Teuchos_ParameterList.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/algorithm/PreconditionerComposition.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

//! PreconditionerYosida
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerYosida is inspired by the Yosida method
 */
class PreconditionerYosida:
    public PreconditionerComposition
{
public:

    /** @name Public Types
     */
    //@{
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    typedef MapEpetra                               map_Type;
    typedef MatrixEpetraStructured<Real>            matrixBlock_Type;
    typedef MatrixEpetraStructuredView<Real>        matrixBlockView_Type;
    typedef MatrixEpetra<Real>                      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;
    typedef VectorEpetra                            vector_Type;
    typedef boost::shared_ptr<vector_Type>          vectorPtr_Type;

    typedef Preconditioner                          super_Type;
    typedef boost::shared_ptr<super_Type>           superPtr_Type;

    typedef ComposedOperator<Preconditioner>        preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>  preconditionerPtr_Type;

    typedef boost::shared_ptr<FESpace<mesh_Type, map_Type> >  FESpacePtr_Type;

    typedef Teuchos::ParameterList                  list_Type;
    //@}


    //! @name Constructors, destructor
    //@{
    //! default constructor.
    PreconditionerYosida ( const  boost::shared_ptr<Epetra_Comm>& comm = boost::shared_ptr<Epetra_Comm>() );

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner( matrixPtr_Type& A );

    //! default destructor
    virtual ~PreconditionerYosida();

    //@}

    //! @name  Methods
    //@{
    void createParametersList ( list_Type&         list,
                                const GetPot&      dataFile,
                                const std::string& section,
                                const std::string& subsection = "Yosida" );

    //! Return an estimation of the conditionement number of the preconditioner
    double condest ();

    //! Build the preconditioner
    int buildPreconditioner ( matrixPtr_Type& A );

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

    //! Method to setup the solver using Teuchos::ParameterList
    /*!
        @param list Teuchos::ParameterList object
     */
    virtual void setParameters ( Teuchos::ParameterList& list );

    //! Setter for the FESpace
    /*!
        This method set the pointer for the FESpaces needed
        for the construction of the operators Ap, Fp and Mp.
        @param uFESpace Boost::shared_ptr on the FESpace for the velocity
        @param pFESpace Boost::shared_ptr on the FESpace for the pressure
     */
    void setFESpace ( FESpacePtr_Type uFESpace, FESpacePtr_Type pFESpace );

    //! Setter for the timestep
    /*!
        This method set the timestep used to compute the approximate Schur complement.
        @param timestep Timestep used to compute the solution of the Navier-Stokes equations
     */
    void setTimestep ( const Real& timestep );

    //@}

protected:

    int             M_velocityBlockSize;
    int             M_pressureBlockSize;

    FESpacePtr_Type M_uFESpace;
    FESpacePtr_Type M_pFESpace;
    Real            M_timestep;

    ADRAssembler<mesh_Type, matrixBlock_Type, vector_Type> M_adrVelocityAssembler;

    // todo: Remove the member dataFile (bad programmation)
    GetPot          M_dataFile;
    string          M_fluidPrec;
    string          M_fluidDataSection;
    string          M_schurPrec;
    string          M_schurDataSection;

private:
    PreconditionerYosida ( const PreconditionerYosida& P ) :
        PreconditionerComposition ( P.M_comm ) {}
    PreconditionerYosida ( const boost::shared_ptr<PreconditionerYosida>& /*P*/ ) {}

};

inline Preconditioner* createYosida()
{
    return new PreconditionerYosida();
}
namespace
{
static bool registerYosida = PRECFactory::instance().registerProduct ( "Yosida", &createYosida );
}

} // namespace LifeV

#endif /* PRECONDITIONERYOSIDA_HPP */
