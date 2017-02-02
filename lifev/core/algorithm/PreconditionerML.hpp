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
    @brief ML preconditioner

    @author Simone Deparis <simone.deparis@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 09-11-0006
 */

#ifndef _MLPRECONDITIONER_HPP_
#define _MLPRECONDITIONER_HPP_

#include <ml_MultiLevelPreconditioner.h>

#include <lifev/core/LifeV.hpp>

#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>

namespace LifeV
{

//! PreconditionerML - Class of multilevels preconditioner
/*!
  @author Simone Deparis   <simone.deparis@epfl.ch>
  @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
  @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
*/
class PreconditionerML:
    public Preconditioner
{
public:

    //! @name Public Types
    //@{

    typedef Preconditioner                       super;

    typedef ML_Epetra::MultiLevelPreconditioner  prec_raw_type;
    typedef std::shared_ptr<prec_raw_type>     prec_type;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;

    //@}


    //! @name Constructors & Destructor
    //@{
    //! Empty constructor.
#ifdef HAVE_MPI
    PreconditionerML ( std::shared_ptr<Epetra_Comm> comm = std::shared_ptr<Epetra_Comm> ( new Epetra_MpiComm ( MPI_COMM_WORLD ) ) );
#else
    PreconditionerML ( std::shared_ptr<Epetra_Comm> comm = std::shared_ptr<Epetra_Comm> ( new Epetra_SerialComm ) );
#endif

    //! destructor.
    virtual ~PreconditionerML();

    //! Constructor from a matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
    */
    PreconditionerML ( operator_type& matrix );

    //@}

    //! @name Methods
    //@{

    //! Build a preconditioner based on the given matrix
    /*!
      @param matrix Matrix upon which construct the preconditioner
     */
    Int buildPreconditioner ( operator_type& matrix );

    //! Reset the preconditioner
    void resetPreconditioner();

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    virtual void createParametersList ( list_Type& list,
                                        const GetPot& dataFile,
                                        const std::string& section,
                                        const std::string& subSection )
    {
        createMLList ( list, dataFile, section, subSection, M_comm->MyPID() == 0 );
    }

    //! Create the list of parameters of the preconditioner
    /*!
      @param list A Parameter list to be filled
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
      @param subSection The subsection in "dataFile" where to find data about the preconditioner
     */
    static void createMLList ( list_Type& list,
                               const GetPot& dataFile,
                               const std::string& section,
                               const std::string& subSection = "ML",
                               const bool& verbose = true );

    //! Apply the inverse of the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int ApplyInverse ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->ApplyInverse ( vector1, vector2 );
    }

    //! Apply the preconditioner on vector1 and store the result in vector2
    /*!
      @param vector1 Vector to which we apply the preconditioner
      @param vector2 Vector to the store the result
     */
    virtual Int Apply ( const Epetra_MultiVector& vector1, Epetra_MultiVector& vector2 ) const
    {
        return M_preconditioner->Apply ( vector1, vector2 );
    }

    //! Show informations about the preconditioner
    virtual void showMe ( std::ostream& output = std::cout ) const;

    //@}


    //! @name Set Methods
    //@{

    //! Set the data of the preconditioner using a GetPot object
    /*!
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
     */
    void setDataFromGetPot ( const GetPot&      dataFile,
                             const std::string& section );

    //! Set the matrix to be used transposed (or not)
    /*!
      @param useTranspose If true the preconditioner is transposed
     */
    Int SetUseTranspose ( bool useTranspose = false )
    {
        return M_preconditioner->SetUseTranspose (useTranspose);
    }

    //! Set the coordinate to be used for the visualization of the aggregates
    /*!
      According to the Trilinos::ML manual one should setup the following variables in ML:
      <ul>
        <li> viz: enable
        <li> viz: output format
        <li> x-coordinates
        <li> y-coordinates
        <li> z-coordinates
      </ul>
      This is the reason why one has to provide these informations. In order to get an output of the aggregation one has to set
      the variable /ML/visualization/enable to true in the datafile.
      @see M.W. Gee, C.M. Siefert, J.J. Hu, R.S. Tuminaro and M.G. Sala. ML 5.0 Smoothed Aggregation User's Guide. Sandia National Laboratories, 2006
      @param xCoord Shared pointer on a vector of the x coordinates of the vertices of the mesh
      @param yCoord Shared pointer on a vector of the y coordinates of the vertices of the mesh
      @param zCoord Shared pointer on a vector of the z coordinates of the vertices of the mesh
     */
    void setVerticesCoordinates (std::shared_ptr<std::vector<Real> > xCoord,
                                 std::shared_ptr<std::vector<Real> > yCoord,
                                 std::shared_ptr<std::vector<Real> > zCoord);

    //@}


    //! @name Get Methods
    //@{

    //! Return An estimation of the condition number of the preconditioner
    Real condest ();

    //! Return a raw pointer on the preconditioner
    super::prec_raw_type* preconditioner();

    //! Return a shared pointer on the preconditioner
    super::prec_type preconditionerPtr()
    {
        return M_preconditioner;
    }

    //! Return the type of preconditioner
    std::string preconditionerType()
    {
        return M_precType;
    }

    //! Return true if the preconditioner is transposed
    bool UseTranspose()
    {
        return M_preconditioner->UseTranspose();
    }

    //! Return the Range map of the operator
    const Epetra_Map& OperatorRangeMap() const
    {
        return M_preconditioner->OperatorRangeMap();
    }

    //! Return the Domain map of the operator
    const Epetra_Map& OperatorDomainMap() const
    {
        return M_preconditioner->OperatorDomainMap();
    }

    //@}

protected:

    list_Type  M_IfpackSubList;
    std::shared_ptr<Epetra_Comm> M_comm;

private:

    operator_type           M_operator;

    prec_type               M_preconditioner;

    bool                    M_analyze;

    bool                    M_visualizationDataAvailable;
    std::shared_ptr<std::vector<Real> > M_xCoord;
    std::shared_ptr<std::vector<Real> > M_yCoord;
    std::shared_ptr<std::vector<Real> > M_zCoord;

};


inline Preconditioner* createML()
{
    return new PreconditionerML();
}
namespace
{
static bool registerML = PRECFactory::instance().registerProduct ( "ML", &createML );
}


} // namespace LifeV

#endif
