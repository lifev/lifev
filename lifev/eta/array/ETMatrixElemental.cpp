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
    @brief This file contains the implementation of the ETMatrixElemental class

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 06 2011
 */

#include <lifev/eta/array/ETMatrixElemental.hpp>

namespace LifeV
{

// ===================================================
// Constructors & Destructor
// ===================================================

ETMatrixElemental::
ETMatrixElemental (const UInt& nbRow, const UInt& nbCol )
    :   M_rowIndices (nbRow, 0),
        M_columnIndices (nbCol, 0),
        M_nbRow (nbRow),
        M_nbColumn (nbCol),
        M_rawData (new Real*[nbRow])
{
    for (UInt iRow (0); iRow < nbRow; ++iRow)
    {
        M_rawData[iRow] = new Real[nbCol];
    }

    zero();
}

ETMatrixElemental::
ETMatrixElemental (const ETMatrixElemental& mat)
    :   M_rowIndices (mat.M_rowIndices),
        M_columnIndices (mat.M_columnIndices),
        M_nbRow (mat.M_nbRow),
        M_nbColumn (mat.M_nbColumn),
        M_rawData ( new Real*[M_nbRow])
{
    for (UInt iRow (0); iRow < M_nbRow; ++iRow)
    {
        M_rawData[iRow] = new Real[M_nbColumn];
        for (UInt iCol (0); iCol < M_nbColumn; ++iCol)
        {
            M_rawData[iRow][iCol] = mat.M_rawData[iRow][iCol];
        }
    }
}

ETMatrixElemental::
~ETMatrixElemental()
{
    for (UInt iRow (0); iRow < M_nbRow; ++iRow)
    {
        delete M_rawData[iRow];
    }
    delete M_rawData;
}

// ===================================================
// Methods
// ===================================================

void
ETMatrixElemental::
showMe ( std::ostream& out ) const
{
    out << " Local matrix : " << M_nbRow << " x " << M_nbColumn << std::endl;
    for (UInt i (0); i < M_nbRow; ++i)
    {
        for (UInt j (0); j < M_nbColumn; ++j)
        {
            out << "[" << i << "][" << j << "] " << M_rawData[i][j] << " ";
        };
        out << std::endl;
    }
}


} // Namespace LifeV
