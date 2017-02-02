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
    @brief A short description of the file content

    @author Samuel Quinodoz <samuel.quinodoz@epfl.ch>
    @date 07 Jun 2011

    This file contains the the VectorBlockMonolithicEpetraView implementation.
 */

#ifndef VECTOR_BLOCK_MONOLITHIC_EPETRA_VIEW_H
#define VECTOR_BLOCK_MONOLITHIC_EPETRA_VIEW_H 1

#include <lifev/core/array/VectorEpetra.hpp>

#include <lifev/core/LifeV.hpp>

#include <boost/shared_ptr.hpp>

namespace LifeV
{

//! VectorBlockMonolithicEpetraView - class representing a block in a VectorBlockMonolithicEpetra
/*!
  @author Samuel Quinodoz

  The VectorBlockMonolithicEpetraView class contains data related
  to block of a vector. It is useful to setup a clean and easy-to-use blocks management.

  For more information about the block structures in LifeV, see \ref BlockAlgebraPage "this page".
 */
class VectorBlockMonolithicEpetraView
{
public:

    //! @name Public Types
    //@{

    //! Typedef for the underlying vector type
    typedef VectorEpetra vector_Type;

    //@}


    //! @name Constructor & Destructor
    //@{

    //! Default constructor.
    VectorBlockMonolithicEpetraView();

    //! Copy constructor
    /*!
      Beware that using the copy constructor, the view is copied
      but not the underlying vector (both views still look at the
      same vector, not at different copies).
     */
    VectorBlockMonolithicEpetraView ( const VectorBlockMonolithicEpetraView& otherView );

    //! Default destructor
    ~VectorBlockMonolithicEpetraView();

    //@}


    //! @name Methods
    //@{

    //! Print the informations about the VectorBlockMonolithicEpetraView
    /*!
      @param output Stream where to print the informations
     */
    void showMe (std::ostream& output = std::cout) const;

    //! Assembly procedure
    /*!
      This procedure should always have the same behaviour and the same
      syntax than the corresponding method in MatrixEpetra.

      @param GID the global index of the element within the block viewed.
      @param value the value to be added
     */
    Int sumIntoGlobalValues ( const Int GID, const Real value ) const;

    //@}


    //! @name Set Methods
    //@{

    /*! Set all the informations relative to the block
     *  @param firstIndex First index in the block
     *  @param blockSize Number of indices in the block
     *  @param vector Vector from which the view has to be extracted
     */
    void setup ( const UInt& firstIndex,
                 const UInt& blockSize,
                 vector_Type* vector );

    //@}


    //! @name Get Methods
    //@{

    //! Returns the size of the block
    UInt blockSize() const
    {
        return M_blockSize;
    }

    //! Returns the index in the block
    UInt firstIndex() const
    {
        return M_firstIndex;
    }

    //! Returns the last index in the block
    UInt lastValidIndex() const
    {
        return M_lastValidIndex;
    }

    //! Return the shared_pointer of the Epetra_FEVector
    vector_Type* vectorPtr() const
    {
        return M_vector;
    }

    //@}


private:

    //! @name Private Methods
    //@{

    //! No assignement operator
    /*!
      The assignement operator is disabled in order to avoid confusing notation
      like vector.block(0) = vector.block(1), that would only copy the views
      and not the blocks.
     */
    VectorBlockMonolithicEpetraView operator= (const VectorBlockMonolithicEpetraView& otherView);

    //@}

    UInt M_blockSize;
    UInt M_firstIndex;
    UInt M_lastValidIndex;
    vector_Type* M_vector;
};

} // Namespace LifeV

#endif /* VECTORBLOCKVIEWEPETRA_H */
