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
     @brief This file contains the definition of the ETMatrixElemental class

     @date 06/2011
     @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
 */

#ifndef ET_MATRIX_ELEMENTAL_HPP
#define ET_MATRIX_ELEMENTAL_HPP

#include <vector>

#include <boost/shared_ptr.hpp>

#include <lifev/core/LifeV.hpp>

// next needed for ROW_MAJOR ...
#include <lifev/core/array/MatrixEpetra.hpp>

namespace LifeV
{

//! class ETMatrixElemental  A class for describing an elemental matrix
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

    This class is meant to represent an ETMatrixElemental. It contains mainly:

    <ol>
        <li> Simple constructors
        <li> Methods to access and modify the indexes
        <li> Methods to store and access the entries
        <li> Assembly methods to put the local entries in a global matrix
    </ol>

*/
class ETMatrixElemental
{

public:

    //! @name Constructors & Destructor
    //@{

    //! Constructor with the minimal interface: number of columns and of rows are provided.
    ETMatrixElemental (const UInt& nbRow, const UInt& nbCol );

    //! Copy constructor (including deep copy of the data)
    ETMatrixElemental (const ETMatrixElemental& mat );

    //! Destructor
    virtual ~ETMatrixElemental();

    //@}

    //! @name Operators
    //@{

    //@}


    //! @name Methods
    //@{

    //! Put zero all the data stored
    // This method is defined in class to allow the compiler
    // to optimize it easily (used repeatedly during the assembly)
    void zero()
    {
        for (UInt i (0); i < M_nbRow; ++i)
        {
            for (UInt j (0); j < M_nbColumn; ++j)
            {
                M_rawData[i][j] = 0.0;
            }
        }
    }

    //! Assembly procedure for a matrix or a block of a matrix
    /*!
    This method puts the values stored in this elemental matrix into the global
    matrix passed as argument, using the positions given in the global indices
    stored.
    */
    // Method defined in class to allow compiler optimization
    // as this class is used repeatedly during the assembly
    template <typename MatrixType>
    void pushToGlobal (MatrixType& mat)
    {
        mat.addToCoefficients ( M_nbRow, M_nbColumn,
                                rowIndices(), columnIndices(),
                                M_rawData, Epetra_FECrsMatrix::ROW_MAJOR);
    }

    //! Assembly procedure for a matrix or a block of a matrix passed in a shared_ptr
    /*!
    This method puts the values stored in this elemental matrix into the global
    matrix passed as argument, using the positions given in the global indices
    stored.

    This is a partial specialization of the other assembly procedure of this class.
    */
    // Method defined in class to allow compiler optimization
    // as this class is used repeatedly during the assembly
    template <typename MatrixType>
    void pushToGlobal (boost::shared_ptr<MatrixType> mat)
    {
        mat->addToCoefficients ( M_nbRow, M_nbColumn,
                                 rowIndices(), columnIndices(),
                                 M_rawData, Epetra_FECrsMatrix::ROW_MAJOR);
    }

    //! Ouput method for the sizes and the stored values
    void showMe ( std::ostream& out = std::cout ) const;

    //@}

    //! @name Set Methods
    //@{

    //! Setter for the value in the local position (iloc,jloc)
    Real& element (const UInt& iloc, const UInt& jloc)
    {
        ASSERT (iloc < M_nbRow, "Try to access an element out of the elemental matrix (row)");
        ASSERT (jloc < M_nbColumn, "Try to access an element out of the elemental matrix (column)");
        return M_rawData[iloc][jloc];
    }

    //! Setter for the global index corresponding to the iloc row of the local matrix
    void setRowIndex (const UInt& iloc, const UInt& iglobal)
    {
        ASSERT (iloc < M_nbRow, "Try to set an index out of the elemental matrix (row)");
        M_rowIndices[iloc] = iglobal;
    }

    //! Setter for the global index of the rows of the local matrix
    void setRowIndex (const std::vector<Int>& indicesVector)
    {
        M_rowIndices = indicesVector;
    }

    //! Setter for the global index corresponding to the jloc column of the local matrix
    void setColumnIndex (const UInt& jloc, const UInt& jglobal)
    {
        ASSERT (jloc < M_nbColumn, "Try to set an index out of the elemental matrix (row)");
        M_columnIndices[jloc] = jglobal;
    }

    //! Setter for the global index of the columns of the local matrix
    void setColumnIndex (const std::vector<Int>& indicesVector)
    {
        M_columnIndices = indicesVector;
    }


    //@}


    //! @name Get Methods
    //@{

    //! Getter for the data stored in the given elemental position
    const Real& element (const UInt& iloc, const UInt& jloc) const
    {
        ASSERT (iloc < M_nbRow, "Try to get an element out of the elemental matrix (row)");
        ASSERT (jloc < M_nbColumn, "Try to get an element out of the elemental matrix (column)");
        return M_rawData[iloc][jloc];
    }

    //! Getter for the full set of data
    Real* const* rawData() const
    {
        return M_rawData;
    }

    //! Getter for the global indices of the rows
    const std::vector<Int>& rowIndices() const
    {
        return M_rowIndices;
    }

    //! Getter for the global indices of the columns
    const std::vector<Int>& columnIndices() const
    {
        return M_columnIndices;
    }

    //@}

private:

    //! @name Private Methods
    //@{

    //! No empty constructor, as we want at least the sizes to be defined.
    ETMatrixElemental();

    //! No need for an assignement operator
    ETMatrixElemental operator= (const ETMatrixElemental&);

    //@}

    // Vectors storing the global indices corresponding to the entries stored.
    // Here UInt would be more suitable, but it does not work with assembly functions.
    std::vector<Int> M_rowIndices;
    std::vector<Int> M_columnIndices;

    // Number of rows
    UInt M_nbRow;
    // Number of columns
    UInt M_nbColumn;

    // Raw data (stored as double pointer to avoid unnecessary conversions)
    Real** M_rawData;
};

}
#endif


