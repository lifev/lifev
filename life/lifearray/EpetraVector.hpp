/* -*- mode: c++ -*-

  This file is part of the LifeV library

  Author(s): Gilles Fourestey <gilles.fourestey@epfl.ch>
             Simone Deparis <simone.deparis@epfl.ch>
       Date: 2006-10-04

  Copyright (C) 2006 EPFL

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/**
   \file EpetraVector.hpp
   \author Gilles Fourestey <gilles.fourestey@epfl.ch>
             Simone Deparis <simone.deparis@epfl.ch>
   \date 2004-10-26
 */

#ifndef _EPETRAVECTOR_HPP_
#define _EPETRAVECTOR_HPP_

#include "Epetra_FEVector.h"
#include "life/lifecore/life.hpp"
//#include <life/lifearray/pattern.hpp>
#include "life/lifealg/EpetraMap.hpp"

#include <Epetra_Export.h>


//@@
//#define OFFSET 0

namespace LifeV
{
////////////////////////////////////////////////////////////////
//
//  Epetra Vector format Wrapper
//
///////////////////////////////////////////////////////////////



class EpetraVector
{
public:

    typedef Epetra_FEVector vector_type;
    typedef Real data_type;

//    EpetraVector(); //!< default constructor : NULL pattern
    //
    // Note that the constructors MUST be based on an existing pattern
    //
//     EpetraVector( const CSRPatt &_pattern,
//                   Epetra_Comm&   comm );
//     EpetraVector( const CSRPatt&    _pattern,
//                   const Epetra_Map& _mapRow );//,
//                  const Epetra_Map& _mapCol );

    EpetraVector( const EpetraMap& _map, EpetraMapType maptype = Unique );

    // The following constructor does not use the EpetraMap.
    //EpetraVector( const Epetra_BlockMap&    _map );

    EpetraVector( const EpetraVector& _vector);

    EpetraVector( const EpetraVector& _vector, EpetraMapType maptype);

    EpetraVector( const EpetraVector& _vector, EpetraMapType maptype,
                  Epetra_CombineMode combineMode);

    // Copies _vector to FEvector that comes as Multivector
    EpetraVector( const Epetra_MultiVector&    _vector,
                  boost::shared_ptr<EpetraMap> _map,
                  EpetraMapType                maptype );


//! Copies _vector to a vector which resides only on the processor "reduceToProc"
    EpetraVector( const EpetraVector& _vector, const int reduceToProc);

    // The following constructor does not use the EpetraMap.
    //! Copies _vector to a vector with different map (usually the repeated map from EpetraMatrix)
    //    EpetraVector( const EpetraVector&     _vector,
    //                  const Epetra_BlockMap&  _map    );

    data_type& operator[] (const UInt row);
    const data_type& operator[] (const UInt row) const;
    data_type& operator() (const UInt row);
    const data_type& operator() (const UInt row) const;


    //! if row is mine returns the LID
    //! if row is not mine and if the numCpus > 1, returns -1
    //! if row is not mine and if the numCpus == 1, asserts
    int checkLID(const UInt row) const;

    //! if row is mine sets this[row] = value and return true
    //! if row is not mine and if the numCpus > 1, returns false
    //! if row is not mine and if the numCpus == 1, asserts
    bool checkAndSet(const UInt row, const data_type& value);

    //! Set the row row of the vector to value. If it isn't on this processor,
    //! store it and send it and send it at next GlobalAssemble
    int replaceGlobalValues(const std::vector<int> rVec, const std::vector<double> datumVec);

    //! insert a global value. After insertion, you will have to call global assemble.
    int sumIntoGlobalValues (const int GID, const double value);

    double Norm2()   const;
    double NormInf() const;

    void Norm2  (double* res) const;
    void NormInf(double* res) const;

    ~EpetraVector() {};
    vector_type& getEpetraVector()             {return M_epetraVector;}
    const vector_type& getEpetraVector() const {return M_epetraVector;}

    const int  size() const { return M_epetraVector.GlobalLength(); }

    void spy ( std::string const &filename ) const;

    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Calls Import with Add.
    EpetraVector& operator=(const EpetraVector& _vector);

    //! copies the value of a Epetra_MultiVector u (assumed of width 1). If the map is not the same,
    //! try to import the values. Calls Import with Add.
    EpetraVector& operator=(const Epetra_MultiVector& _vector);

    /* adds to this the vector _vector with an offset.
       typically to do: (u,p) += p or (u,p) += u.
       Note: the nodes to add are taken by the map of _vector, hence:
       a) if this has a unique map: then _vector should also (otherwise run time error)
       b) if this has a repeated map: then _vector should also. (otherwise wrong)
     */
    EpetraVector&  add(const EpetraVector& _vector,
                       const int           offset = 0);


    /* set this to a subset of  _vector with an offset.
       typically to do: p = (u,p) or u = (u,p).
       Note: the nodes to add are taken by the map of this, hence:
       a) if _vector has a unique map: then this should also (otherwise run time error)
       b) if _vector has a repeated map: then this should also. (otherwise wrong)
     */
    EpetraVector& subset(const EpetraVector& _vector,
                                   const int           offset = 0);

    //! if the map is not the same, try to import values
    EpetraVector& operator+=(const EpetraVector& _vector);
    EpetraVector& operator-=(const EpetraVector& _vector);

    EpetraVector& operator*=(data_type t);

    data_type operator*(EpetraVector const& a) const;

    int GlobalAssemble(Epetra_CombineMode mode=Add) { return  M_epetraVector.GlobalAssemble(); }

    const Epetra_Comm& Comm() const { return BlockMap().Comm(); }

    const EpetraMapType getMaptype() const  { return  M_maptype; }

    const EpetraMap&  getMap() const  { return  *M_epetraMap; }

    const Epetra_BlockMap& BlockMap() const;




private:

    //! copies the value of a vector u. If the map is not the same,
    //! try to import the values. Let you decide wether to add or replace shared nodes:
    //! note:
    //!  if the original source vector _vector is not repeated : use Import
    //!  if the original source vector _vector is repeated : use Export
    /*
        CombineMode Valuse:
        Add         Components on the receiving processor will be added together.
        Insert 	    Off-processor components will be inserted into locations on receiving processor replacing existing values.
        InsertAdd 	Off-processor components will be inserted into locations on receiving processor and added to existing values.
        Average 	Off-processor components will be averaged with existing components on the receiving processor.
        (+ Zero and AbsMax, probably never useful)
    */
    EpetraVector& Import (const Epetra_FEVector& _vector,
                          Epetra_CombineMode combineMode);

    //! copies the value of this to a vector _vector. If the map is not the same,
    //! try to import the values. Let you decide wether to add or replace shared nodes:
    //! note: tested only if the destination source vector _vector is not repeated
    /*
        CombineMode Valuse:
        Add         Components on the receiving processor will be added together.
        Insert 	    Off-processor components will be inserted into locations on receiving processor replacing existing values.
        InsertAdd 	Off-processor components will be inserted into locations on receiving processor and added to existing values.
        Average 	Off-processor components will be averaged with existing components on the receiving processor.
        (+ Zero and AbsMax, probably never useful)
    */
    EpetraVector& Export (const Epetra_FEVector& _vector,
                          Epetra_CombineMode combineMode);


    vector_type      M_epetraVector;

    boost::shared_ptr<EpetraMap>  M_epetraMap;
    EpetraMapType    M_maptype;

};


EpetraVector operator * (EpetraVector::data_type t, const EpetraVector& _vector);


} // end namespace LifeV

//@@
//#undef OFFSET

#endif
