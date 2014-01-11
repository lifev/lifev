//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2010 EPFL, Politecnico di Milano, INRIA

 This library is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as
 published by the Free Software Foundation; either version 2.1 of the
 License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA

************************************************************************
*/
//@HEADER

/*!
    @file
    @brief This file contains the PreconditionerLSC class.

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @date 08-11-2010
 */

#ifndef PRECONDITIONERLSC_HPP
#define PRECONDITIONERLSC_HPP 1

#include <lifev/core/LifeV.hpp>

#ifdef LIFEV_HAVE_TEKO

#include <boost/shared_ptr.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/algorithm/PreconditionerTeko.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

// Teuchos includes
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

// Teko-Package includes
#include <Teko_Utilities.hpp>
#include <Teko_InverseFactory.hpp>
#include <Teko_InverseLibrary.hpp>
#include <Teko_BlockPreconditionerFactory.hpp>
#include <Teko_InvLSCStrategy.hpp>
#include <Teko_LSCPreconditionerFactory.hpp>

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

namespace LifeV
{

//! PreconditionerLSC
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerLSC class provides the LSC block preconditioner
 *  available in the Teko package of Trilinos
 */
class PreconditionerLSC:
    public PreconditionerTeko
{
public:

    /** @name Public Types
     */
    //@{
    typedef Preconditioner                          super_Type;

    typedef Teko::Epetra::EpetraBlockPreconditioner preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>  preconditionerPtr_Type;
    typedef RegionMesh<LinearTetra>                 mesh_Type;
    typedef MapEpetra                               map_Type;
    typedef boost::shared_ptr<FESpace<mesh_Type, map_Type> >  FESpacePtr_Type;
    typedef MatrixEpetra<Real>                      matrix_Type;
    typedef boost::shared_ptr<matrix_Type>          matrixPtr_Type;

    typedef Teuchos::ParameterList                  list_Type;
    //@}


    /** @name Constructors, destructor
     */
    //@{
    //! default constructor.
#ifdef HAVE_MPI
    PreconditionerLSC ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerLSC ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner( matrixPtr_Type& A );

    //! default destructor
    virtual ~PreconditionerLSC();

    //@}

    /** @name  Methods
     */

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

    void setFESpace ( FESpacePtr_Type uFESpace,
                      FESpacePtr_Type pFESpace );

    void createParametersList ( list_Type&         list,
                                const GetPot&      dataFile,
                                const std::string& section,
                                const std::string& subSection = "LSC" );

    //! Return an estimation of the conditionement number of the preconditioner
    double condest ();

    //! Return the name of the preconditioner to be used in the factory
    std::string preconditionerType()
    {
        return M_precType;
    }

    //! Build the preconditioner
    int buildPreconditioner ( matrixPtr_Type& A );

    int numBlocksRows() const;
    int numBlocksCols() const;

protected:

    std::string M_precType;
    int         M_velocityBlockSize;
    int         M_pressureBlockSize;
    boost::shared_ptr<Epetra_Comm> M_comm;

};

inline Preconditioner* createLSC()
{
    return new PreconditionerLSC();
}
namespace
{
static bool registerLSC = PRECFactory::instance().registerProduct ( "LSC", &createLSC );
}

} // namespace LifeV

#endif // LIFEV_HAVE_TEKO

#endif /* PRECONDITIONERLSC_HPP */
