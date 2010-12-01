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
    @brief A short description of the file content

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Jul 2010

    A more detailed description of the file (if necessary)
 */

#ifndef ROBININTERFACE_H
#define ROBININTERFACE_H 1

#include <life/lifecore/life.hpp>
#include <lifemc/lifesolver/BlockInterface.hpp>

#include <life/lifearray/EpetraMatrix.hpp>

namespace LifeV
{

//! robinInterface - Class for the Robin coupling of etherogeneaous problems
/*!
    @author Paolo Crosetto

    This class is a base interface for the monolithic Robin coupling of 2 problems. A coupled Dirichlet-Neumann problem
    can be transformed into a Robin one by inheriting from this class. The allocation of the memory for M_robinCoupling
    has to be done outside.

    The Robin conditions are obtained algebraically through linear combination of the lines corresponding to the Dirichlet and Neumann
    coupling conditions already present in the matrix. This is achieved by multiplying the matrix times a restriction
    matrix (M_robinCoupling) which performs the restriction to the interface and the linear combination. Then the result
    is summed to the original matrix. The update of the rhs vector must be done outside this class.
 */
class RobinInterface
{
public:


    //! @name Constructor & Destructor
    //@{

    RobinInterface():
            M_alphaf(),
            M_alphas(),
            M_robinCoupling(),
            M_robinPart(),
            M_rhsVec()
    {}

    ~RobinInterface() {}

    //@}


    //! @name Public Methods
    //@{

    //! method to set the data relative to the Robin coupling from GetPot
    /*!
      \param data: data file
      \param section: the section (usually /robin) in the GetPot file where the parameters are specified
     */
    void setRobinData(const GetPot& data, const std::string& section);

    //! method to initialize the pointer to the robin coupling part of the matrix
    /*!
      \param data: data file
      \param section: the section (usually /robin) in the GetPot file where the parameters are specified
     */
    void setRobinMatrix( BlockInterface::matrix_ptrtype& robinMatrix ) {M_robinPart=robinMatrix;}

    //! method to initialize the pointer to the robin RHS
    /*!
      \param vec: the rhs vector
     */
    void setRobinRhs( BlockInterface::vector_ptrtype& vec ) { M_rhsVec = vec; }

    //! method to apply the robin coupling to the blocks.
    /*!
      Note that the coupling matrix M_robinCoupling must be set before calling this.
      \param blockVector: the vector of blocks to couple.
     */
    void applyRobinCoupling( std::vector<BlockInterface::matrix_ptrtype> blockVector);

    //@}

protected:


    //! @name Protected Methods
    //@{

    void applyRobinCoupling( BlockInterface::matrix_ptrtype firstBlock );

    //@}


    //! @name Protected Members
    //@{

    Real                                 M_alphaf;
    Real                                 M_alphas;
    BlockInterface::matrix_ptrtype       M_robinCoupling;
    BlockInterface::matrix_ptrtype       M_robinPart;
    BlockInterface::vector_ptrtype       M_rhsVec;
    //@}

};

} // Namespace LifeV

#endif /* ROBININTERFACE_H */
