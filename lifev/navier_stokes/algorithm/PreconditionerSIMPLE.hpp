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
    @brief PreconditionerSIMPLE

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 24-01-2011
 */

#ifndef PRECONDITIONERSIMPLE_HPP
#define PRECONDITIONERSIMPLE_HPP 1

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
    //! default constructor
#ifdef HAVE_MPI
    PreconditionerSIMPLE ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerSIMPLE ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner( matrixPtr_Type& A );

    //! default destructor
    virtual ~PreconditionerSIMPLE();

    //@}

    //! @name  Methods
    //@{
    void createParametersList ( list_Type&         list,
                                const GetPot&      dataFile,
                                const std::string& section,
                                const std::string& subsection = "SIMPLE" );

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

    //! Setter for the damping factor
    /*!
        This method set damping factor used to build the preconditioner
        @param dampingFactor Damping factor
    */
    void setDampingFactor ( const Real& dampingFactor );

    //@}

protected:

    int             M_velocityBlockSize;
    int             M_pressureBlockSize;
    FESpacePtr_Type M_uFESpace;
    FESpacePtr_Type M_pFESpace;

    Real            M_dampingFactor;

    string          M_SIMPLEType;

    // todo: Remove the member dataFile (bad programmation)
    GetPot          M_dataFile;
    string          M_fluidPrec;
    string          M_fluidDataSection;
    string          M_schurPrec;
    string          M_schurDataSection;

private:
    PreconditionerSIMPLE ( const PreconditionerSIMPLE& P ) :
        PreconditionerComposition ( P.M_comm ) {}
    PreconditionerSIMPLE ( const boost::shared_ptr<PreconditionerSIMPLE>& /*P*/ ) {}

};

inline Preconditioner* createSIMPLE()
{
    return new PreconditionerSIMPLE();
}
namespace
{
static bool registerSIMPLE = PRECFactory::instance().registerProduct ( "SIMPLE", &createSIMPLE );
}

} // namespace LifeV

#endif /* PRECONDITIONERSIMPLE_HPP */
