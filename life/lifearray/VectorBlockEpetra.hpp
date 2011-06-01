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
    @brief File containing the VectorBlockEpetra

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01 Jun 2011

 */

#ifndef VECTORBLOCKEPETRA_H
#define VECTORBLOCKEPETRA_H 1

#include <life/lifecore/Life.hpp>

#include <life/lifearray/VectorEpetra.hpp>
#include <life/lifearray/MapVector.hpp>

namespace LifeV {

//! VectorBlockEpetra - A block vector structure for the Epetra framework of Trilinos
/*!
    @author Samuel Quinodoz
 */

class VectorBlockEpetra
    : public VectorEpetra
{
public:

    //! @name Public Types
    //@{

    typedef Real data_type;
    typedef MapEpetra map_type;
    typedef MapVector<map_type> mapVector_type;
    typedef MapEpetraType mapType_type;
    typedef Epetra_CombineMode combine_type;

    //@}


    //! @name Constructor & Destructor
    //@{

    VectorBlockEpetra( const map_type& map, const mapType_type& mapType = Unique);

    VectorBlockEpetra( const mapVector_type& mapVector, const mapType_type& mapType = Unique);

    VectorBlockEpetra( const VectorBlockEpetra& vector);

    VectorBlockEpetra( const VectorBlockEpetra& vector, const mapType_type& mapType);

    VectorBlockEpetra( const VectorBlockEpetra& vector, const mapType_type& mapType, const combine_type& combineMode);

    ~VectorBlockEpetra(){}

    //@}


    //! @name Operators
    //@{


    //@}


    //! @name Methods
    //@{


    //@}


    //! @name Set Methods
    //@{


    //@}


    //! @name Get Methods
    //@{


    //@}

private:

    std::vector<UInt> M_blockSize;
    std::vector<UInt> M_blockFirstIndex;

};

} // Namespace LifeV

#endif /* VECTORBLOCKEPETRA_H */
