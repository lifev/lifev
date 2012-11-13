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
    @brief File containing the VectorEpetraStructured class

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 01 Jun 2011
 */

#ifndef _VECTOREPETRASTRUCTURED_HPP_
#define _VECTOREPETRASTRUCTURED_HPP_ 1

#include <lifev/core/LifeV.hpp>

#include <lifev/core/array/MapEpetra.hpp>
#include <lifev/core/array/VectorEpetra.hpp>
#include <lifev/core/array/MapVector.hpp>
#include <lifev/core/array/VectorBlockStructure.hpp>
#include <lifev/core/array/VectorEpetraStructuredView.hpp>


namespace LifeV {

//! VectorEpetraStructured - class of block vector
/*!
  @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>

  The VectorEpetraStructured class contains data related
  to block vector. It is an extension to VectorEpetra where data about blocks have
  been set. For an introduction to the block structures in LifeV, see
  \ref BlockAlgebraPage "this page".

  There are mainly two ways to define a VectorEpetraStructured:
  <ul>
  <li> Construct it using the same syntax as for LifeV::VectorEpetra and the use
  a setter for the structure.
  <li> Construct it directly with the maps of the blocks.
  </ul>
  Both ways are equivalent.

  To access the blocks, one uses then the blockView or block methods.

 */
class VectorEpetraStructured
    : public VectorEpetra
{
public:

    //! @name Public Types
    //@{

    //! Type of data stored
    typedef Real data_type;

    //! Type of the map to be used
    typedef MapEpetra map_type;

    //! Type of the MapVector to be used with this class
    typedef MapVector<map_type> mapVector_type;

    //! Type of the map (Unique/Repeated)
    typedef MapEpetraType mapType_type;

    //! Combine mode
    typedef Epetra_CombineMode combine_type;

    //! Type of the view
    typedef VectorEpetraStructuredView block_type;

    //! Pointer on the view
    typedef boost::shared_ptr<block_type> block_ptrType;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Constructor with the monolithic map
    VectorEpetraStructured( const map_type& map, const mapType_type& mapType = Unique );

    //! Construction with a vector of map
    /*!
      With this constructor, the block structure is automatically deduced from the maps in the
      vector. The monolithic map and vectors are also built by concatenating the different maps
      in the vector.
     */
    VectorEpetraStructured( const mapVector_type& mapVector, const mapType_type& mapType = Unique );

    //! Copy constructor
    VectorEpetraStructured( const VectorEpetraStructured& vector );

    //! Copy constructor with a specified map type (Repeated/Unique)
    VectorEpetraStructured( const VectorEpetraStructured& vector, const mapType_type& mapType );

    //! Copy constructor with specified map type and combine mode
    VectorEpetraStructured( const VectorEpetraStructured& vector, const mapType_type& mapType, const combine_type& combineMode );

    //! Destructor
    ~VectorEpetraStructured(){}

    //@}


    //! @name Set Methods
    //@{

    /*! Set the size of the blocks of the vector
     *  @param blockSizes Sizes of the blocks
     */
    void setBlockStructure( const std::vector<UInt>& blockSizes );

    //! Reset the block structure using the blocks of a vector of map
    /*
      The resulting block structure is symmetric (same block structure
      in the rows and in the columns).

      This method does not involve big computation overhead. Remark that
      it is not possible to change the size of the block vector
      through this method, nor its map.
      @param The MapVector containing the maps
     */
    void setBlockStructure( const mapVector_type& mapVector );

    /*! Set the block structure of the vector
     *  @param blockStructure Structure of the vector
     */
    void setBlockStructure( const VectorBlockStructure& blockStructure );

    //@}


    //! @name Get Methods
    //@{

    //! Getter for the size of the block index
    /*!
      @param index Index of the block
      @return size of the index-th block
     */
    UInt blockSize( const UInt& index ) const;

    //! Getter for the block index
    /*!
      @param index Index of the block
      @param blockView The blockView to be filled
     */
    void blockView( const UInt& index, block_type& blockView );

    //! Getter for the block index
    /*!
      @param index Index of the block
      @return The index-th block
     */
    block_ptrType block( const UInt& index );

    /*! Get the block structure of the vector
     *  @param blockStructure Structure of the vector
     */
    VectorBlockStructure blockStructure() const;

    //@}

private:

    VectorBlockStructure M_blockStructure;

};

} // Namespace LifeV

#endif /* _VECTOREPETRASTRUCTURED_HPP_ */
