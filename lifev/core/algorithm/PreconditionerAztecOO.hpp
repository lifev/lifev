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
    @brief AztecOO preconditioner

    @author Cristiano Malossi <cristiano.malossi@epfl.ch>
    @contributor Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>
    @maintainer Gwenol Grandperrin <gwenol.grandperrin@epfl.ch>

    @date 17-11-2009
 */

#ifndef AZTECOOPRECONDITIONER_HPP
#define AZTECOOPRECONDITIONER_HPP 1

#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/Preconditioner.hpp>

namespace LifeV
{

//! PreconditionerAztecOO - The implementation of Preconditioner for AztecOO preconditioners
/*!
  @author Cristiano Malossi <cristinao.malossi@epfl.ch>

  This class provides the interface for using AztecOO preconditioners with SolverAztecOO.
 */
class PreconditionerAztecOO : public Preconditioner
{
public:

    //! @name Public Types
    //@{

    typedef Preconditioner     super;

    typedef super::operator_raw_type operator_raw_type;
    typedef super::operator_type     operator_type;

    typedef SolverAztecOO           Solver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    PreconditionerAztecOO();

    //! Destructor
    virtual ~PreconditionerAztecOO() {}

    //@}


    //! @name Methods
    //@{

    //! Build the preconditioner
    /*!
      @param A the base matrix for computing the preconditioner
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
    void createParametersList ( list_Type& list, const GetPot& dataFile, const std::string& section, const std::string& subSection );

    //! Show informations about the preconditioner
    virtual void showMe ( std::ostream& output = std::cout ) const;

    //! Compute the condition number of the preconditioner
    /*!
      @return Condition number of the preconditioner
     */
    Real condest();

    //@}


    //! @name Set Methods
    //@{

    //! Set the data of the preconditioner using a GetPot object
    /*!
      @param dataFile A GetPot object containing the data about the preconditioner
      @param section The section in "dataFile" where to find data about the preconditioner
     */
    void setDataFromGetPot ( const GetPot& dataFile, const std::string& section );

    //! Set the external solver (AztecOO)
    /*!
      @param solver reference to the AztecOO solver
     */
    void setSolver ( SolverAztecOO& solver )
    {
        M_solver = &solver;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Return the pointer to the preconditioner
    /*!
      @return always zero because no external precondtioner is used here!
     */
    super::prec_raw_type* preconditioner();

    //! Return the shared pointer to the preconditioner
    /*!
      @Deprecated
     */
    super::prec_type  preconditionerPtr();

    //! Return the name of the preconditioner
    /*!
      @return "AztecOO"
     */
    std::string preconditionerType()
    {
        return "AztecOO";
    }

    //@}

private:

    Solver_Type*           M_solver;
};

inline Preconditioner* createPreconditionerAztecOO()
{
    return new PreconditionerAztecOO();
}

namespace
{
static bool registerAztecOO = PRECFactory::instance().registerProduct ( "AztecOO", &createPreconditionerAztecOO );
}

} // namespace LifeV

#endif // AZTECOOPRECONDITIONER_HPP
