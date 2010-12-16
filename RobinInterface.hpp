/* -*- mode: c++ -*- */
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
    @brief An interface implementing the Robin coupling

    @author Paolo Crosetto <crosetto@iacspc70.epfl.ch>
    @date 09 Jul 2010

    Used for the robin coupling in a monolithic formulation of a multiphysics problem.
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

    \todo Remove this class and try to implement the same coupling otherwise.
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
    void setRobinMatrix( BlockInterface::matrixPtr_Type& robinMatrix ){M_robinPart=robinMatrix;}

    //! method to initialize the pointer to the robin RHS
    /*!
      \param vec: the rhs vector
     */
    void setRobinRhs( BlockInterface::vectorPtr_Type& vec ) { M_rhsVec = vec; }

    //! method to apply the robin coupling to the blocks.
    /*!
      Note that the coupling matrix M_robinCoupling must be set before calling this.
      \param blockVector: the vector of blocks to couple.
     */
    void applyRobinCoupling( std::vector<BlockInterface::matrixPtr_Type> blockVector);

    //@}

protected:


    //! @name Protected Methods
    //@{

    void applyRobinCoupling( BlockInterface::matrixPtr_Type firstBlock );

    //@}


    //! @name Protected Members
    //@{

    Real                                 M_alphaf;
    Real                                 M_alphas;
    BlockInterface::matrixPtr_Type       M_robinCoupling;
    BlockInterface::matrixPtr_Type       M_robinPart;
    BlockInterface::vectorPtr_Type       M_rhsVec;
    //@}

};

} // Namespace LifeV

#endif /* ROBININTERFACE_H */
