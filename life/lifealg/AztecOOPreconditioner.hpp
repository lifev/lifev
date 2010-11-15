//@HEADER
/*
************************************************************************

 This file is part of the LifeV Applications.
 Copyright (C) 2001-2006 EPFL, Politecnico di Milano, INRIA
               2006-2010 EPFL, Politecnico di Milano

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
 *  @file
 *  @brief AztecOO preconditioner
 *
 *  @author Cristiano Malossi <cristiano.malossi@epfl.ch>
 *  @date 17-11-2009
 */

#ifndef AZTECOOPRECONDITIONER_HPP
#define AZTECOOPRECONDITIONER_HPP 1

#include <life/lifealg/SolverTrilinos.hpp>
#include <life/lifealg/EpetraPreconditioner.hpp>

namespace LifeV {

//! AztecOOPreconditioner - The implementation of EpetraPreconditioner for AztecOO preconditioners
/*!
 *  @author Cristiano Malossi
 *
 *  This class provides the interface for using AztecOO preconditioners with SolverTrilinos.
 */
class AztecOOPreconditioner : public EpetraPreconditioner
{
public:

    //! @name Typedefs
    //@{

    typedef EpetraPreconditioner                 super;

    typedef super::operator_raw_type             operator_raw_type;
    typedef super::operator_type                 operator_type;

    typedef SolverTrilinos                       Solver_Type;

    //@}


    //! @name Constructors & Destructor
    //@{

    //! Constructor
    AztecOOPreconditioner();

    //! Destructor
    ~AztecOOPreconditioner() {}

    //@}


    //! @name Methods
    //@{

    //! Build the preconditioner
    /*!
     *  @param A the base matrix for computing the preconditioner
     */
    int buildPreconditioner( operator_type& A );

    //! Reset the preconditioner
    void precReset();

    //! AztecOO preconditioner is set?
    /*!
     *  @return true
     */
    bool set() const { return M_preconditionerCreated; }

    //@}


    //! @name Set Methods
    //@{

    //! Load parameters from a GetPot files
    /*!
     *  @param dataFile the GetPot file
     *  @param section the section inside the GetPot file containing the parameters
     */
    void setDataFromGetPot ( const GetPot& dataFile, const std::string& section, const std::string& subSection = "AztecOO" );

    //! Set the external solver (AztecOO)
    /*!
     *  @param solver reference to the AztecOO solver
     */
    void setSolver( SolverTrilinos& solver ) { M_solver = &solver; }

    //@}


    //! @name Get Methods
    //@{

    //! Compute the condition number of the preconditioner
    /*!
     *  @return Condition number of the preconditioner
     */
    Real Condest();

    //! Return the pointer to the preconditioner
    /*!
     *  @return always zero because no external precondtioner is used here!
     */
    super::prec_raw_type* getPrec();

    //! Return the shared pointer to the preconditioner
    /*!
     *  @Deprecated
     */
    super::prec_type  getPrecPtr();

    //! Return the name of the preconditioner
    /*!
     *  @return "AztecOO"
     */
    std::string precType() { return "AztecOO"; }

    //@}

private:

    Solver_Type*           M_solver;
};

inline EpetraPreconditioner* createAztecOOPreconditioner()
{
    return new AztecOOPreconditioner();
}

namespace
{
	static bool registerAztecOO = PRECFactory::instance().registerProduct( "AztecOO", &createAztecOOPreconditioner );
}

} // namespace LifeV

#endif // AZTECOOPRECONDITIONER_HPP
