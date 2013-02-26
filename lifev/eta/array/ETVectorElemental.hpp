//@HEADER
/*
*******************************************************************************

   Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
   Copyright (C) 2010 EPFL, Politecnico di Milano, Emory UNiversity

   This file is part of the LifeV library

   LifeV is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   LifeV is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, see <http://www.gnu.org/licenses/>


*******************************************************************************
*/
//@HEADER

/*!
 *   @file
     @brief This file contains the definition of the ETVectorElemental class

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ET_VECTOR_ELEMENTAL_HPP
#define ET_VECTOR_ELEMENTAL_HPP

#include <vector>

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.hpp>

namespace LifeV
{

//! class ETVectorElemental  A class for describing an elemental vector
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class is meant to represent an ETVectorElemental. It contains mainly:

    <ol>
        <li> Simple constructors
        <li> Methods to access and modify the indexes
        <li> Methods to store and access the entries
        <li> Assembly methods to put the local entries in a global vector
    </ol>

*/
class ETVectorElemental
{

public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor with the minimal interface: number of rows has to be provided.
    explicit ETVectorElemental (const UInt& nbRow)
        :   M_rowIndices (nbRow, 0),
            M_nbRow (nbRow),
            M_rawData (new Real[nbRow])
    {
        zero();
    }

    //! Copy constructor (including deep copy of the data)
    ETVectorElemental (const ETVectorElemental& vec)
        :   M_rowIndices (vec.M_rowIndices),
            M_nbRow (vec.M_nbRow),
            M_rawData ( new Real[M_nbRow])
    {
        for (UInt iRow (0); iRow < M_nbRow; ++iRow)
        {
            M_rawData[iRow] = vec.M_rawData[iRow];
        }
    }

    //! Destructor
    ~ETVectorElemental()
    {
        delete M_rawData;
    }

    //@}


    //! @name Operators
    //@{

    //! Operator to access (read-only) the iloc-th element
    const Real& operator[] (const UInt& iloc) const
    {
        ASSERT (iloc < M_nbRow, "Access out of range");
        return M_rawData[iloc];
    }

    //! Operator to access (read-write) the iloc-th element
    Real& operator[] (const UInt& iloc)
    {
        ASSERT (iloc < M_nbRow, "Access out of range");
        return M_rawData[iloc];
    }

    //@}


    //! @name Methods
    //@{

    //! Put zero all the data stored
    void zero()
    {
        for (UInt i (0); i < M_nbRow; ++i)
        {
            M_rawData[i] = 0.0;
        }
    }

    //! Assembly procedure for a vector or a block of a vector
    /*!
    This method puts the values stored in this elemental vector into the global
    vector passed as arguement, using the positions given in the global indices
    stored in this elemental vector.
    */
    template <typename VectorType>
    void pushToGlobal (VectorType& vec)
    {
        for (UInt iRow (0); iRow < M_nbRow; ++iRow)
        {
            vec.sumIntoGlobalValues ( M_rowIndices[iRow], M_rawData[iRow] );
        }
    }

    //! Assembly procedure for a vector or a block of a vector passed in a shared_ptr
    /*!
    This method puts the values stored in this elemental vector into the global
    vector passed as argument, using the positions given in the global indices
    stored in this elemental vector.

    This is a partial specialization of the other assembly procedure of this class.
    */
    template <typename VectorType>
    void globalAssembly (boost::shared_ptr<VectorType> vec)
    {
        for (UInt iRow (0); iRow < M_nbRow; ++iRow)
        {
            vec->sumIntoGlobalValues ( M_rowIndices[iRow], M_rawData[iRow] );
        }
    }


    //! Ouput method for the sizes and the stored values
    void showMe ( std::ostream& out = std::cout ) const
    {
        out << " Elemental vector : Size " << M_nbRow << std::endl;
        for (UInt i (0); i < M_nbRow; ++i)
        {
            out << "[" << i << "] " << M_rawData[i] << " ";
        }
        out << std::endl;
    }

    //@}

    //! @name Set Methods
    //@{

    //! Setter for the value in the local position (iloc)
    Real& element (const UInt& iloc)
    {
        ASSERT (iloc < M_nbRow, "Access out of range");
        return M_rawData[iloc];
    }

    //! Setter for the global index corresponding to the iloc row of the local vector
    void setRowIndex (const UInt& iloc, const UInt& iglobal)
    {
        M_rowIndices[iloc] = iglobal;
    }

    //! Setter for all the global indices using a std::vector
    void setRowIndex (const std::vector<Int>& indicesVector)
    {
        M_rowIndices = indicesVector;
    }

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the data stored in the given elemental position
    const Real& element (const UInt& iloc) const
    {
        ASSERT (iloc < M_nbRow, "Access out of range");
        return M_rawData[iloc];
    }

    //! Getter for the full set of data
    Real* rawData() const
    {
        return M_rawData;
    }


    //! Getter for the global indices of the rows
    const std::vector<Int>& rowIndices() const
    {
        return M_rowIndices;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor, as we want at least the sizes to be defined.
    ETVectorElemental();

    //! No need for an assignement operator
    ETVectorElemental operator= (const ETVectorElemental&);

    //@}

    // Vectors storing the global indices corresponding to the entries stored
    // Here UInt would be more suitable, but it does not work with assembly functions.
    std::vector<Int> M_rowIndices;

    // Number of rows
    UInt M_nbRow;

    // Raw data
    Real* M_rawData;
};

}
#endif


