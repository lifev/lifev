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
    @brief PreconditionerPCD

    @author Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 29-11-2010
 */

#ifndef PRECONDITIONERPCD_HPP
#define PRECONDITIONERPCD_HPP 1

#include <boost/shared_ptr.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

#include <Teuchos_ParameterList.hpp>

// Tell the compiler to ignore specific kind of warnings:
#pragma GCC diagnostic warning "-Wunused-variable"
#pragma GCC diagnostic warning "-Wunused-parameter"

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/algorithm/PreconditionerComposition.hpp>
#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/solver/ADRAssembler.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include <lifev/core/fem/BCBase.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>

namespace LifeV
{

//! PreconditionerPCD
/*!
 *  @author Gwenol Grandperrin
 *
 *  The PreconditionerPCD class provides the PCD block preconditioner
 */
class PreconditionerPCD:
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
    typedef boost::shared_ptr<matrix_Type>          matrixPtr_type;
    typedef VectorEpetra                            vector_Type;
    typedef boost::shared_ptr<vector_Type>          vectorPtr_Type;

    typedef Preconditioner                          super_Type;
    typedef boost::shared_ptr<super_Type>           superPtr_Type;

    typedef ComposedOperator<Preconditioner>        preconditioner_Type;
    typedef boost::shared_ptr<preconditioner_Type>  preconditionerPtr_Type;

    typedef boost::shared_ptr<FESpace<mesh_Type, map_Type> >  FESpacePtr_Type;
    typedef boost::shared_ptr<BCHandler>            BCHandlerPtr_Type;

    typedef Teuchos::ParameterList                  list_Type;
    //@}


    //! @name Constructors, destructor
    //@{
    //! default constructor
#ifdef HAVE_MPI
    PreconditionerPCD ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerPCD ( boost::shared_ptr<Epetra_Comm> comm = boost::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! constructor from matrix A.
    //! @param A EpetraMatrix<double> matrix upon which construct the preconditioner
    //    IfpackPreconditioner( matrixPtr_type& A );

    //! default destructor
    virtual ~PreconditionerPCD();

    //@}

    //! @name  Methods
    //@{
    void createParametersList ( list_Type&         list,
                                const GetPot&      dataFile,
                                const std::string& section,
                                const std::string& subsection = "PCD" );

    //! Return an estimation of the conditionement number of the preconditioner
    double condest ();

    //! Update the vector beta of the convective term in Fp
    /*!
        This method updates the value of beta.
        @param beta New vector beta to be used to built the convective term
     */
    void updateBeta ( const vector_Type& beta );


    //! Build the preconditioner
    int buildPreconditioner ( matrixPtr_type& A );

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

    //! Setter for the BCHandler
    /*!
        This method set the pointer for the FESpaces needed
        for the construction of the operators Ap, Fp and Mp.
        @param bchPtr pointer to a BCHandler boject
    */
    void setBCHandler ( BCHandlerPtr_Type bchPtr );

    //! Setter for the timestep
    /*!
        This method set the timestep used to compute Fp.
        @param timestep Timestep used to compute the solution of the Navier-Stokes equations
     */
    void setTimestep ( const Real& timestep );

    //! Setter for the viscosity
    /*!
        This method set the viscosity used to compute Fp.
        @param viscosity Viscosity used to compute the solution of the Navier-Stokes equations
     */
    void setViscosity ( const Real& viscosity );

    //! Setter for the density
    /*!
        This method set the density used to compute Fp.
        @param density Density used to compute the solution of the Navier-Stokes equations
     */
    void setDensity ( const Real& density );

    //! Setter to know if we used B or -B in the discretization of the Navier-Stokes equations
    /*!
        @param useMinusDivergence is true if -B has been used.
     */
    void setUseMinusDivergence ( const bool& useMinusDivergence );

    //@}

protected:

    int         M_velocityBlockSize;
    int         M_pressureBlockSize;
    FESpacePtr_Type M_uFESpace;
    FESpacePtr_Type M_pFESpace;

    Real        M_timestep;
    Real        M_viscosity;
    Real        M_density;
    vectorPtr_Type  M_beta;

    ADRAssembler<mesh_Type, matrixBlock_Type, vector_Type> M_adrPressureAssembler;
    ADRAssembler<mesh_Type, matrixBlock_Type, vector_Type> M_adrVelocityAssembler;

    // todo: Remove the member dataFile (bad programmation)
    GetPot      M_dataFile;
    BCHandlerPtr_Type M_bcHandlerPtr;
    std::string M_fluidPrec;
    std::string M_fluidPrecDataSection;
    std::string M_pressureLaplacianPrec;
    std::string M_pressureLaplacianPrecDataSection;
    std::string M_pressureMassPrec;
    std::string M_pressureMassPrecDataSection;
    std::string M_pressureLaplacianOperator;
    std::string M_pressureMassOperator;
    bool        M_setApBoundaryConditions;
    bool        M_setFpBoundaryConditions;
    bool        M_setMpBoundaryConditions;
    bool        M_fullFactorization;
    bool        M_schurOperatorReverseOrder;
    bool        M_useStiffStrain;
    bool        M_enableTransient;
    Real        M_divergenceCoeff;
    bool        M_recomputeNormalVectors;

    std::string M_inflowBoundaryType;
    std::string M_outflowBoundaryType;
    std::string M_characteristicBoundaryType;

    vectorPtr_Type M_normalVectors;

private:
    PreconditionerPCD ( const PreconditionerPCD& P ) :
        PreconditionerComposition ( P.M_comm ) {}
    PreconditionerPCD ( const boost::shared_ptr<PreconditionerPCD>& /*P*/ ) {}

    void computeNormalVectors();

    vectorPtr_Type computeRobinCoefficient();

    static Real fZero ( const Real& /* t */,
                        const Real& /* x */,
                        const Real& /* y */,
                        const Real& /* z */,
                        const ID& /* i */ )
    {
        return 0.0;
    }

    void setBCByBoundaryType ( matrixPtr_type Ap, UInt ApOffset,
                               matrixPtr_type Fp, UInt FpOffset,
                               matrixPtr_type Mp, UInt MpOffset );

};

inline Preconditioner* createPCD()
{
    return new PreconditionerPCD();
}
namespace
{
static bool registerPCD = PRECFactory::instance().registerProduct ( "PCD", &createPCD );
}

} // namespace LifeV

#endif /* PRECONDITIONERPCD_HPP */
