/* -*- mode: c++ -*-

 This file is part of the LifeV Applications.

 Author(s): Gilles Fourestey <gilles.fourestey@epfl.ch>
 Simone Deparis <simone.deparis@epfl.ch>
 Date: 2006-10-04

 Copyright (C) 2009 EPFL

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 USA
 */
/*!
  \file EpetraVector.hpp

  \version 1.0
  \date 2006-10-04

  \author Simone Deparis <simone.deparis@epfl.ch>
  \author Gilles Fourestey <gilles.fourestey@epfl.ch>

  \version 1.10
  \date 15/10/2009
  \author Cristiano Malossi <cristiano.malossi@epfl.ch>

  - Added new operators: (+=, -=, *=, /=) and (+ - * /) for Element By Element and Vector-Scalar operations;
  - Added new logic operators: &&, ||, ! ;
  - Added new comparison operators: >=, <=, >, <, ==, != ;
  - Added Abs() method, and renamed scalar product to Dot() ;
  - Added ShowMe() method;
  - Added some doxygen (not finished yet) and reordered functions and methods;
*/

#ifndef _EPETRAVECTOR_HPP_
#define _EPETRAVECTOR_HPP_

#include <life/lifecore/life.hpp>
#include <life/lifealg/EpetraMap.hpp>

#include <Epetra_FEVector.h>
#include <Epetra_Export.h>

namespace LifeV {

//! EpetraVector - The Epetra Vector format Wrapper
/*!
 *  @author Gilles Fourestey, Simone Deparis, Cristiano Malossi
 *
 *  The EpetraVector class provides a general interface for the Epetra_Vector of Trilinos.
 *
 *  TODO Finish to reorder the class and to add more doxygen comments.
 */
class EpetraVector
{
public:

    typedef Epetra_FEVector vector_type;
    typedef Real data_type;

    //! @name Constructors & Destructor
    //@{

    //! Constructor - Using Maps
    EpetraVector( const EpetraMap& _map, EpetraMapType maptype = Unique );

    //! Constructor - Using Maps
    EpetraVector( const boost::shared_ptr< EpetraMap >& _map, EpetraMapType maptype = Unique );

    //! Constructor - Using Vector (without using map)
    EpetraVector( const EpetraVector& vector );

    //! Constructor - Using Vector (without using map)
    EpetraVector( const EpetraVector& vector, EpetraMapType maptype );

    //! Constructor - Using Vector (without using map)
    EpetraVector( const EpetraVector& vector,
                  EpetraMapType maptype,
                  Epetra_CombineMode combineMode );

    //! Constructor - Copies vector to FEvector that comes as Multivector
    EpetraVector( const Epetra_MultiVector& vector,
                  boost::shared_ptr< EpetraMap > _map,
                  EpetraMapType maptype );

    //! Constructor - Copies vector to a vector which resides only on the processor "reduceToProc"
    EpetraVector( const EpetraVector& vector, const int reduceToProc );

    //! Destructor
    ~EpetraVector() {}

    //@}


    //! @name Methods
    //@{

    int GlobalAssemble( Epetra_CombineMode mode = Add )
    {
        return M_epetraVector.GlobalAssemble( mode );
    }

    //! if row is mine returns the LID
    //! if row is not mine and if the numCpus > 1, returns -1
    //! if row is not mine and if the numCpus == 1, asserts
    int checkLID( const UInt row ) const;

    //! if row is mine sets this[row] = value and return true
    //! if row is not mine and if the numCpus > 1, returns false
    //! if row is not mine and if the numCpus == 1, asserts
    bool checkAndSet( const UInt row, const data_type& value, UInt offset = 0 );

    //! Set the row row of the vector to value. If it isn't on this processor,
    //! store it and send it and send it at next GlobalAssemble
    int replaceGlobalValues( std::vector< int >& rVec, std::vector< double >& datumVec );

    //! insert a global value. After insertion, you will have to call global assemble.
    int sumIntoGlobalValues( const int GID, const double value );

    /* adds to this the vector vector with an offset.
     typically to do: (u,p) += p or (u,p) += u.
     Note: the nodes to add are taken by the map of vector, hence:
     a) if this has a unique map: then vector should also (otherwise run time error)
     b) if this has a repeated map: then vector should also. (otherwise wrong)
     */
    EpetraVector& add( const EpetraVector& vector, const int offset = 0 );

    /* set this to a subset of  vector with an offset.
     typically to do: p = (u,p) or u = (u,p).
     Note: the nodes to add are taken by the map of this, hence:
     a) if vector has a unique map: then this should also (otherwise run time error)
     b) if vector has a repeated map: then this should also. (otherwise wrong)
     */
    EpetraVector& subset( const EpetraVector& vector, const UInt offset = 0 );

    //! set this to a subset of  vector with an offset.
    /*!
        similar to subset( const EpetraVector& , const int ), but with
        additional parameters:
        @param vector  vector from which to copy data
        @param map     map from which to select indeces to copy
        @param offset1 offset to apply to input vector
        @param offset2 offset to apply to this vector
    */
    EpetraVector& subset( const EpetraVector& vector,
                          const EpetraMap& map,
                          const UInt offset1,
                          const UInt offset2 );

    //! set this to a subset of  vector with an offset.
    /*!
        similar to subset( const EpetraVector& , const int ), but with
        additional parameters:
        @param vector  Epetra_MultiVector, instead of EpetraVector, from which to copy data
        @param map     map from which to select indeces to copy
        @param offset1 offset to apply to input vector
        @param offset2 offset to apply to this vector
        @param column  column of the multivector from which to extract the data
    */
    EpetraVector& subset(const Epetra_MultiVector& vector,
                         const EpetraMap&    map,
                         const UInt          offset1,
                         const UInt          offset2,
                         const UInt          column = 0);


    double Norm1   () const;
    double Norm2   () const;
    double NormInf () const;
    double MinValue() const;
    double MaxValue() const;

    void MeanValue ( double* res ) const;
    void Norm1     ( double* res ) const;
    void Norm2     ( double* res ) const;
    void NormInf   ( double* res ) const;
    void MinValue  ( double* res ) const;
    void MaxValue  ( double* res ) const;

    void Norm1     ( double& res ) const;
    void Norm2     ( double& res ) const;
    void NormInf   ( double& res ) const;
    void MinValue  ( double& res ) const;
    void MaxValue  ( double& res ) const;

    //! spy - save the values of the matrix into a file
    /*!
     * \param filename - File where to save the EpetraVector
     *
     * To read the file in Matlab type load filename;
     */
    void spy( std::string const &filename ) const;

    void ShowMe( std::ostream& output = std::cout ) const;

    //! Abs - Replace the vector with his abs.
    /*!
     * \param dataFile - Name and path of data file
     */
    void Abs( void );

    //! Abs - Compute the abs of a vector
    /*!
     * \param vector - Output vector with the abs of the EpetraVector
     */
    void Abs( EpetraVector& vector );

    //! Dot - Compute the scalar product of two vectors
    /*!
     * \param vector - Second vector for the scalar product
     */
    data_type Dot( const EpetraVector& vector ) const;

    //! Dot - Compute the scalar product of two vectors
    /*!
     * \param vector - Second vector for the scalar product
     * \param scalarProduct - results
     */
    void Dot( const EpetraVector& vector, data_type& scalarProduct );

    //@}


    //! @name Operators
    //@{

    //! Access operators
          data_type& operator[]( const UInt row );
    const data_type& operator[]( const UInt row ) const;
          data_type& operator()( const UInt row );
    const data_type& operator()( const UInt row ) const;

    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Calls Import with Add.
    EpetraVector& operator=( const EpetraVector& vector );

    //! copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
    //! try to import the values. Calls Import with Add.
    EpetraVector& operator=( const Epetra_MultiVector& vector );

    EpetraVector& operator=( data_type t );

    //! Element by Element operations (if the map is not the same, try to import values)
    EpetraVector& operator+=( const EpetraVector& vector );
    EpetraVector& operator-=( const EpetraVector& vector );
    EpetraVector& operator*=( const EpetraVector& vector );
    EpetraVector& operator/=( const EpetraVector& vector );

    //! Element by Element operations (do not modify the vector of the class)
    const EpetraVector operator+( const EpetraVector& vector ) const;
    const EpetraVector operator-( const EpetraVector& vector ) const;
    const EpetraVector operator*( const EpetraVector& vector ) const;
    const EpetraVector operator/( const EpetraVector& vector ) const;

    //! Operations with scalar values (modify the vector of the class)
    EpetraVector& operator+=( const data_type& scalar );
    EpetraVector& operator-=( const data_type& scalar );
    EpetraVector& operator*=( const data_type& scalar );
    EpetraVector& operator/=( const data_type& scalar );

    //! Operations with scalar values (do not modify the vector of the class)
    const EpetraVector operator+( const data_type& scalar ) const;
    const EpetraVector operator-( const data_type& scalar ) const;
    const EpetraVector operator*( const data_type& scalar ) const;
    const EpetraVector operator/( const data_type& scalar ) const;

    //! operator==
    /*!
     * Return a vector containing 1 where vector elements are == scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator==( const Real& scalar );

    //! operator!=
    /*!
     * Return a vector containing 1 where vector elements are != scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator!=( const Real& scalar );

    //! operator<
    /*!
     * Return a vector containing 1 where vector elements are < scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator<( const Real& scalar );

    //! operator>
    /*!
     * Return a vector containing 1 where vector elements are > scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator>( const Real& scalar );

    //! operator<=
    /*!
     * Return a vector containing 1 where vector elements are <= scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator<=( const Real& scalar );

    //! operator>=
    /*!
     * Return a vector containing 1 where vector elements are >= scalar;
     * \param scalar - Value for the comparison.
     */
    EpetraVector operator>=( const Real& scalar );

    //! Logical operator&&
    /*!
     * Return a vector containing one where both elements are != zero;
     * \param vector - Vector for the logical comparison.
     */
    EpetraVector operator&&( const EpetraVector& vector );

    //! Logical operator||
    /*!
     * Return a vector containing one where one of the elements is != zero;
     * \param vector - Vector for the logical comparison.
     */
    EpetraVector operator||( const EpetraVector& vector );

    //! Logical operator&&
    /*!
     * Return a vector containing one where the vector is equal to zero;
     */
    EpetraVector operator!( void );

    //@}


    //! @name Set Methods
    //@{

    //! Sets the combine mode for the import/export operations.
    /*!
        Most of the LifeV library is structured to use combine mode
        equal to Add. In some cases (cf test_filters) it is necessary
        to discard the data coming from other processors.
        @param combineMode combien mode to use for this vector from now on
     */
    void setCombineMode( Epetra_CombineMode combineMode );

    //! Sets the combine mode for the import/export operations to default.
    /*!
        Most of the LifeV library is structured to use combine mode
        equal to Add.
     */
    void setDefaultCombineMode( );

    //@}

    //! @name Get Methods
    //@{

    const Epetra_Comm& Comm() const
    {
        return BlockMap().Comm();
    }

    vector_type& getEpetraVector()
    {
        return M_epetraVector;
    }
    const vector_type& getEpetraVector() const
    {
        return M_epetraVector;
    }

    const Epetra_BlockMap& BlockMap() const
    {
        return M_epetraVector.Map();
    }

    EpetraMapType getMaptype() const
    {
        return M_maptype;
    }

    const EpetraMap& getMap() const
    {
        return *M_epetraMap;
    }

    const boost::shared_ptr< EpetraMap > getMap_ptr() const
    {
        return M_epetraMap;
    }

    const Epetra_Map& getEpetra_Map() const
    {
        return *( M_epetraMap->getMap( M_maptype ) );
    }

    int size() const
    {
        return M_epetraVector.GlobalLength();
    }

    //@}

private:

    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Let you decide wether to add or replace shared nodes:
    //! note:
    //!  if the original source vector vector is not repeated : use Import
    //!  if the original source vector vector is repeated : use Export
    /*
     CombineMode Valuse:
     Add         Components on the receiving processor will be added together.
     Insert 	    Off-processor components will be inserted into locations on receiving processor replacing existing values.
     InsertAdd 	Off-processor components will be inserted into locations on receiving processor and added to existing values.
     Average 	Off-processor components will be averaged with existing components on the receiving processor.
     (+ Zero and AbsMax, probably never useful)
     */
    EpetraVector& Import( const Epetra_FEVector& vector, Epetra_CombineMode combineMode );

    //! copies the value of this to a vector vector. If the map is not the same,
    //! try to import the values. Let you decide wether to add or replace shared nodes:
    //! note: tested only if the destination source vector vector is not repeated
    /*
     CombineMode Valuse:
     Add         Components on the receiving processor will be added together.
     Insert 	    Off-processor components will be inserted into locations on receiving processor replacing existing values.
     InsertAdd 	Off-processor components will be inserted into locations on receiving processor and added to existing values.
     Average 	Off-processor components will be averaged with existing components on the receiving processor.
     (+ Zero and AbsMax, probably never useful)
     */
    EpetraVector& Export( const Epetra_FEVector& vector, Epetra_CombineMode combineMode );

    boost::shared_ptr< EpetraMap > M_epetraMap;
    EpetraMapType                  M_maptype;
    vector_type                    M_epetraVector;

    Epetra_CombineMode             M_combineMode;

};

EpetraVector operator-( const EpetraVector& vector );

EpetraVector operator+( const EpetraVector::data_type& scalar, const EpetraVector& vector );
EpetraVector operator-( const EpetraVector::data_type& scalar, const EpetraVector& vector );
EpetraVector operator*( const EpetraVector::data_type& scalar, const EpetraVector& vector );

} // end namespace LifeV

#endif
