/*
  This file is part of the LifeV library
  Copyright (C) 2001,2002,2003,2004 EPFL, INRIA and Politechnico di Milano

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
/*----------------------------------------------------------------------*
|
| $Header: /cvsroot/lifev/lifev/life/lifearray/Attic/pattern.hpp,v 1.7 2004-09-06 11:03:49 winkelma Exp $
|
|
| #Version  0.1 Experimental   07/7/00. Luca Formaggia & Alessandro Veneziani  |
|
| #Attempt to define a pattern for vectorial problems 12/10/01. Alain Gauthier
|  Construction of a class VBRPatt which is required for the interface
|  with AZTEC.
|
| #Purposes Defines Patterns for normal and mixed matrices
|  A pattern defines the graph of a sparse matrix. The pattern class builds that
|  graph starting from a Degree of Freedom (DOF) object.
|
|  The patterns are meant to be used also as way of interrogating the local connectivities of
|  degrees of freedom. In order to avoid ambiguities the degrees of freedom ARE ALWAYS numbered
|  from 1 and they are indicated by using the type name ID (in this way we avoid ambiguity with UInt).
|  On the other hand, we may have a FORTRAN or C-like style for the RAW data storig the pattern
|  which may be CSR of MSR format. In this way we hopefully will be able to better interface external
|  Linear Algebra packages Fortran-style (like SPARSEKIT).
|
|
|
*----------------------------------------------------------------------*/
#ifndef PATTERN_OFFSET
#define PATTERN_OFFSET 0 // for the Fortran (PATTERN_OFFSET=1) vs C (PATTERN_OFFSET=0) numbering
#endif
#ifndef _PATTERN_HH
#define _PATTERN_HH
#include "lifeV.hpp"
#ifndef INDEX_T
#define INDEX_T UInt
#endif

#ifndef _VEC_UNKNOWN_HH
#include "vecUnknown.hpp"
#endif



#include<set>
#include<algorithm>
#include<string>
//#include<functional>
#include "bareItems.hpp"

namespace LifeV
{
const INDEX_T PatternOffset= PATTERN_OFFSET;

// We DEFINITIVELY need NAMESPACES (LUCA)
/*! \class PatternDefs
  This class containes some useful functions and typedefs which will be
  common to ALL patterns (base and mixed ones). I use BareEdges for the
  dynamic pattern, because I have provided the comparison function. Actually,
  just a pair<UInt,UInt> should work all the same.
*/
class PatternDefs
{
public:
  typedef INDEX_T Index_t; //!< Type for indices
  typedef vector<Index_t> Container; //!< Container for the actual (raw) pattern
  typedef Container::size_type Diff_t; //!< type for differences (offsets)
  /* Some useful converters */
  Index_t _d2i(ID const d) const; //!< From Identifier (DOF) to Index
  ID   _i2d(Index_t const i) const; //!< Form Index to Identifier
  Diff_t _i2o(Index_t const i) const; //!< From index to differences (offsets)
  Diff_t _d2o(ID const d) const; //!< from identifier to differences.
  Container & _i2d(Container & list_of_indices)const; // as Before, yet working on containers
  Container & _d2i(Container & list_of_dof)const;
  Container & _d2o(Container & list_of_dof)const;
 protected:
  typedef set<BareEdge,cmpBareItem<BareEdge> > DynPattern; //!< Container for the dynamic pattern
};


/*!\class BasePattern
  This is the building block for the pattern classes.
*/
class BasePattern :public PatternDefs
{
public:

  BasePattern();
  BasePattern(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol); //!< Here I give the n. of non zeros and matrix dimensions.
  BasePattern(class CSRPatt const &RightHandCSRP);

  inline UInt nRows() const; //!< Numer of Rows
  inline UInt nCols() const; //!< Number of Columns
  inline UInt nNz() const;   //!< Total nnz

  inline bool isEmpty() const; //!< Empty pattern

  bool diagFirst()const {return _diagfirst;}; //!< Diagonal is the first item furnished by  the raw pattern data

  void showMe(bool const verbose=false, ostream & out=cout ) const; //!< some info for the curious

 protected:
  //! It builds the dynamic pattern. It is the standard routine for a
  //! SYMMETRIC pattern associated to a single degree of freedom object.
  //! It uses the local pattern provided through the DOF object.
  //! If nbcomp > 1, it reproduces the pattern for (nbcomp X nbcomp) blocks
  template<typename DOF1>
  bool setpatt(DOF1 const & dof1, DynPattern & dynpatt,
	       UInt const nbcomp = 1);

  // Miguel 12/2003: new version handling patterns coming from IP stabilization
  /*! It builds the dynamic pattern. It uses the local pattern provided through
   *  the DOF object, augmented by the couplings needed for IP stabilization.
   *  If nbcomp > 1, it reproduces the pattern for (nbcomp X nbcomp) blocks.
   */
  template<typename DOF, typename MESH>
  bool setpatt(const DOF& dof, const MESH& mesh, DynPattern & dynpatt,
	       UInt const nbcomp);

  //! Version for mixed patterns
  //! It builds the dynamic pattern. It is the standard routine for a
  //! pattern associated to two degrees of freedom object.
  //! it uses the MixedLocalPattern class for the local patterns.
  //! See the documentation in the implementation part for more details
  template<typename DOF1, typename DOF2>
  bool setpatt(DOF1 const & dof1, DOF2 const & dof2, DynPattern & dynpatt,
	       UInt const bRows= 1, UInt const bCols= 1);

   static const UInt _defaultsize = 0; // MUST BE 0

  UInt _nnz;
  UInt _nrows;
  UInt _ncols;
  bool _filled; //!< True if the pattern has been filled!
  bool _diagfirst; //!< True if the pattern has diagfirst!
};

/*                         SPECIALISED PATTERNS                    */


///////////////////////////////////////////////////
///
///      C S R Pattern
///
//////////////////////////////////////////////////

class CSRPatt:
  public BasePattern
{
public:

  CSRPatt();
  CSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol);

  //! Constructor where I give the raw data from the exterior TODO: As it
  //! stands, This is a pretty useless contructor, unless we build a smarter
  //! version which takes in input raw data also in a other types. This
  //! indeed will help interfacing with external libraries We will need to
  //! change it to something of the type
  //
  //! template<typename Const_Iter>
  //! CSRPatt(Const_Iter &ex_ia, Const_Iter  &ex_ja )
  //! where Const_iter is an iterator to a sequence (then  a pointer as well)
  //! LUCA: I will do it soon!
  CSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja );

  CSRPatt(const CSRPatt &RightHandCSRP);
  //! Constructors for single DOF (square matrix), possibly with more than
  //! one component
  template<typename DOF> CSRPatt(DOF const  & dof, UInt const nbcomp = 1);
  template<typename DOF>
  bool buildPattern(DOF const  & dof, const UInt nbcomp);
  //! Constructors for two DOF (in general non-square matrix)
  //! This constructors are required if we want to use the pattern in a
  //! MixedPattern. It may be useful also stand alone. Look at the
  //! documentation of buildPattern<DOF1,DOF2> to have details.
  template<typename DOF1,typename DOF2>
    CSRPatt(DOF1 const  & dof1,DOF2 const  & dof2,
	    UInt const bRows= 1, UInt const bCols= 1);

  template<typename DOF1,typename DOF2>
    bool buildPattern(DOF1 const& dof1,DOF2 const  & dof2,
		      UInt const bRows= 1, UInt const bCols= 1);

  //! Version which construct two patterns: patt and its transpose one,
  //! in addition an other Container for the access to
  //! the columns of the transpose matrix is built.
  //! @author Alain, Nov. 2002
  template<typename DOF1,typename DOF2>
    friend void
    buildPattTpatt(DOF1 const& dof1, DOF2 const  & dof2,
		   CSRPatt &Patt, CSRPatt &Tpatt,
		   UInt const bRows= 1, UInt const bCols= 1);

  CSRPatt & operator= (const CSRPatt&  RhCsr);

  //! Methods that returns row data in different fashions.
  Container  ia() const {return _ia;};
  Container  ja() const {return _ja;};
  Container  jaT() const {return _jaT;};
  Index_t * giveRawCSR_ia() {return &(_ia.front());}; //!< Give ia (in a raw form)
  Index_t * giveRawCSR_ja() {return &(_ja.front());}; //!< Give ja (in a raw form)
  Index_t * giveRawCSR_jaT() {return &(_jaT.front());}; //!< Give jaT (in a raw form)
  Index_t const * giveRawCSR_ia() const {return &(_ia.front());}; //!< Give ia (in a raw form)
  Index_t const * giveRawCSR_ja() const {return &(_ja.front());}; //!< Give ja (in a raw form)
  Index_t const * giveRawCSR_jaT() const {return &(_jaT.front());}; //!< Give jaT (in a raw form)
  Container & give_ia() {return _ia;}; //!< Give ia (as container)
  Container & give_ja() {return _ja;};//!< Give ja (as container)
  Container & give_jaT() {return _jaT;};//!< Give jaT (as container)
  Container const & give_ia() const {return _ia;}; //!< Give ia (as container)
  Container const & give_ja() const {return _ja;};//!< Give ja (as container)
  Container const & give_jaT() const {return _jaT;};//!< Give jaT (as container)

  inline pair<UInt,UInt> giveMinMax() const; //!< Min and Max elements in a row.
  inline UInt nbNeighbours(ID const d)const; //!< N of neighbours of the DOF numbeered d . BEWARE d is INCLUDED!
  inline ID neighbour(ID const i,ID const d) const;//!< the i-th (start from 1) neighbour of dof d. The first is d itself
  inline void neighbours(ID const d, Container & start) const;//!< put neighbours of dof d in a list

  //! The following template function extracts a row (useful to implement A*b), returning two sequences:
  //! One with the column numbering (coldata) of one with the corresponding offsets in the vector holding
  //! the matrix value. coldata   is given as a sequence of Index_t (i.e. the content depend by the value of PATETRN_OFFSET), while
  //! the position sequence are always offsets (i.e. starting from 0)
  //! Iter is either an iterator to a container or a pointer. The function Returns the
  //! n. or row elements.
  //! IMPORTANT: for efficency reason the sequences pointed by coldata and position MUST have been
  //! dimensioned boforehand in order to have the sufficient dimension (use giveMinMax) NO CHECKS
  //! ARE MADE!!!!!!!!!
  template <typename Iter>
  inline UInt row(Diff_t const row, Iter  coldata, Iter position) const;

  //! Here the routines whcih locate the position in the vector containing the matrix value of an entry (i,j), which may be a couple
  //! of indices or of ID's (degree of freedoms identifiers)
  inline pair<Diff_t,bool> locate_dof(ID const i,ID const j) const; //!< Locate position of DOF ID couple (i,j)
  inline pair<Diff_t,bool> locate_index(Index_t const i,Index_t const j) const;//!< Locate position of index couple (i,j)

  void showMe(bool const verbose=false,ostream& c=cout) const ; //!< pattern visualization
  void spy(string const & filname="matrice.m") const; //!< pattern visualization a la Matlab

  //! column-concatenation of two CSR block patterns
  friend CSRPatt colUnify(CSRPatt const &patt1, CSRPatt const &patt2);
  //! column-concatenation of one CSR block patterns and ncolZero null columns
  friend CSRPatt colUnify(CSRPatt const &patt1, UInt const ncolZero);
  friend CSRPatt colUnify(UInt const ncolZero, CSRPatt const &patt1);
  //! row-concatenation of two CSR block patterns
  friend CSRPatt rowUnify(CSRPatt const &patt1, CSRPatt const &patt2);
  //! row-concatenation of one CSR block patterns and nrowZero null rows
  friend CSRPatt rowUnify(CSRPatt const &patt1, UInt const nrowZero);
  friend CSRPatt rowUnify(UInt const nrowZero, CSRPatt const &patt1);

  //! Return a block diagonal pattern of nblock blocks
  friend CSRPatt diagblockMatrix(CSRPatt const &patt, UInt const nblock);


protected:
  Diff_t _row_off(Index_t i)const {//!< Row offset for row index i
    return _i2o(_ia[_i2o(i)]);}
  inline bool isThere(Index_t const i,Index_t const j) const; //!< Are indices (i,j) in the pattern
  pair<Diff_t,bool> locate_pattern(Index_t const i,Index_t const j) const;//!< Are indices (i,j) in the pattern
  Container _ia; //!< point to the rows entries
  Container _ja; //!< contain col indices for each row
  Container _jaT;//!< point to the col entries of the transpose pattern
};

///////////////////////////////////////////////////
///
///      V B R Pattern
/// It holds a CSR pattern of variable size blocks
///
//////////////////////////////////////////////////

class VBRPatt: public CSRPatt
/*!\class VBRPatt
   The block pattern is given by the class CSRPatt:
   _ia : points to the location in _ja of the first block
         entry in each block row.
   _ja : contains the block column indices of the block pattern.

   BEWARE : the interpretation of _nnz, _nrows and _ncols of the class
   BasePattern is valid for the BLOCK pattern, that is
   _nnz : number of non zero BLOCKS,
   _nrows : number of BLOCK rows,
   _ncols : number of BLOCK columns.

   VBRPatt: VBR format accepts blocks of variable size.
            This pattern is more general than needed. In fact we work
	    with blocks having a fixed size.
*/
{
 public:
  VBRPatt(){}
  VBRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol);
  VBRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const
	  vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja, const
	  vector<Index_t> &ex_indx, const vector<Index_t> &ex_rpntr, const
	  vector<Index_t> &ex_cpntr);

  VBRPatt(const VBRPatt &RightHandVBRP);

  //! Constructors for single DOF (square matrix), the size of the
  //! blocks is given by GenericVecHdl.
  template<typename DOF>
    VBRPatt(DOF const  & dof, UInt const blockSize);
  bool buildPattern(UInt const blockSize);

  //! Constructors for two DOF (in general non-square matrix)
  //! as for a single DOF, the block pattern is built by CSRPatt and
  //! the size of the blocks is given by GenericVecHdl.
  template<typename DOF1,typename DOF2>
    VBRPatt(DOF1 const  & dof1,DOF2 const  & dof2, UInt const blockSize);

  VBRPatt & operator= (const VBRPatt&  RhVbr);

  //! Methods that returns row data in different fashions.
  inline Container indx() const {return _indx;}
  inline Container rpntr() const {return _rpntr;}
  inline Container cpntr() const {return _cpntr;}

  inline Index_t* giveRawVBR_ia() {return giveRawCSR_ia();} //!< Give ia (in a raw form)
  inline Index_t* giveRawVBR_ja() {return giveRawCSR_ja();} //!< Give ja (in a raw form)
  inline Index_t* giveRawVBR_indx() {return &(_indx.front());} //!< Give indx (in a raw form)
  inline Index_t* giveRawVBR_rpntr() {return &(_rpntr.front());} //!< Give rpntr (in a raw form)
  inline Index_t* giveRawVBR_cpntr() {return &(_cpntr.front());} //!< Give cpntr (in a raw form)


  inline Container & give_indx() {return _indx;}; //!< Give indx (as container)
  inline Container & give_rpntr() {return _rpntr;};//!< Give rpntr (as container)
  inline Container & give_cpntr() {return _cpntr;};//!< Give cpntr (as container)

  // No check of _rpntr size !
  inline Diff_t rbloc(Diff_t const row) const
    //!< Give the block row to which
    //!< belong row
    {
      UInt blsize=_rpntr[1]-_rpntr[0]; // size of square block
      return row/blsize+PatternOffset; // block row offset in which is row
    }

  inline Diff_t cbloc(Diff_t const col) const
    //!< Give the block column to
    //!< which belong col
    {
      UInt blsize=_rpntr[1]-_rpntr[0]; // size of square block
      return col/blsize+PatternOffset; // block row offset in which is row
    }

  inline Diff_t locr(Diff_t const row) const
    //!< Give the local row numbering
    //!< of row in the corresponding block (given by rbloc)
    {
      UInt blsize=_rpntr[1]-_rpntr[0]; // size of square block
      return row%blsize+PatternOffset;  // local row number into the block
    }

  inline Diff_t locc(Diff_t const col) const
     //!< Give the local column numbering
     //!< of row in the corresponding block (given by cbloc)
    {
      UInt blsize=_rpntr[1]-_rpntr[0]; // size of square block
      return col%blsize+PatternOffset;  // local col number into the block
    }


  //! The following template function extracts a row (useful to
  //! implement A*b), returning two sequences:
  //! coldata : column numbering (of elements and NOT blocks)
  //! position: coldata corresponding offsets in the vector holding the matrix
  //!           value.
  //! coldata is given as a sequence of Index_t (i.e. the
  //! content depends on the value of PATTERN_OFFSET), while
  //! the position sequence are always offsets (i.e. starting from 0)
  //! Iter is either an iterator to a container or a pointer. The
  //! function Returns the n. or row elements.
  //! IMPORTANT: for efficiency reason the sequences pointed by coldata
  //! and position MUST have been dimensioned before in order to
  //! have the sufficient dimension (use giveMinMax) NO CHECKS ARE MADE!!!!!
  inline UInt row(Diff_t const row, Container &  coldata, Container & position) const;

  // Here the routines which locate the position in the vector
  // containing the matrix value of an entry (i,j), which may be a
  // couple of indices or of ID's (degree of freedoms identifiers)
  // ALAIN : using locate routines member of class CSRPatt gives
  // information on the block position (element (1,1) of the block).

  //! pattern visualization
  void showMe(bool const verbose=false,ostream& c=cout) const ;
  void spy(string const & filname="matrice.m") const; //!< pattern visualization
                                                     //!< a la Matlab

 private:
  Container _indx; //!< _indx(i) points to the location in the pattern
  //!< of the (1,1) element of the i-th block entry.
  //!< _indx(_nnz+1)=1+number of non zeros elements of the pattern.
  Container _rpntr; //!< _rpntr(i)-_rpntr(0) is the row index of the
  //!< first point row in the i-th block row.
  //!< _rpntr(_nrows+1)=_rpntr(0)+number of element rows.
  Container _cpntr; //!<  the column index of the first point column in
  //!< the i-th block column.
  //!< _cpntr=_rpntr FOR SQUARE PATTERN of SQUARE BLOCKS.
  //!< _cpntr(_ncols+1)=_cpntr(0)+number of element columns.
};

//////////////////////////////////////////////////
///
///      C S R Symmetric Pattern
/// It holds only the upper triangular part
///
//////////////////////////////////////////////////

class CSRPattSymm:
  public BasePattern
{
 public:
  CSRPattSymm();
  CSRPattSymm(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol);
  CSRPattSymm(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol, const vector<Index_t> &ex_ia, const vector<Index_t> &ex_ja );
  CSRPattSymm(const CSRPattSymm &RightHandCSRP);
  template<typename DOF>CSRPattSymm(DOF const  & dof);
  CSRPattSymm & operator= (const CSRPattSymm&  RhCsr);

  template<typename DOF>
  bool buildPattern(DOF const  & dof);  // Pattern builder

  Container  ia() const {return _ia;};
  Container  ja() const {return _ja;};
  Index_t * giveRawCSR_ia() {return &(_ia.front());}; // Give ia (in a raw form)
  Index_t * giveRawCSR_ja() {return &(_ja.front());}; // Give ja (in a raw form)
  Container & give_ia() {return _ia;}; // Give ia (as container)
  Container & give_ja() {return _ja;};// Give ja (as container)
  inline pair<UInt,UInt> giveMinMax() const;

  // DO NOt use them (are very inefficient)
  UInt nbNeighbours(ID const d)const;
  ID neighbour(ID const i,ID const d)const ;
  void neighbours(ID const d, Container & start) const;
  template<typename Iter>
  inline UInt row(Diff_t const row, Iter  coldata, Iter position) const;// extracts a row (useful to implement A*b)

  inline pair<Diff_t,bool> locate_dof(ID const i,ID const j) const;

  inline pair<Diff_t,bool> locate_index(Index_t const i,Index_t const j) const;

  void showMe(bool verbose=false,ostream& c=cout) const; // pattern visualization
  void spy(string const & filname="matrice.m") const; //pattern visualization a la Matlab


  // DOF are ALWAYS numbered a la Fortran: the following definition set the correct numbering
  // as a funtion of current numbering (specified by PATTERN_OFFET)

protected:
  Diff_t _row_off(Index_t i)const {// Row offset for row index i
    return _i2o(_ia[_i2o(i)]);}
  pair<Diff_t,bool> locate_pattern(Index_t const i,Index_t const j) const;
  bool isThere(Index_t i,Index_t j) const;

 private:
  Container _ia;
  Container _ja;
};


///////////////////////////////////////////////////
///
///      M S R Pattern
///  It implements the MSR format
///
//////////////////////////////////////////////////

class MSRPatt:
  public BasePattern
{
 public:
  MSRPatt();
  MSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol);
  MSRPatt(UInt ex_nnz, UInt ex_nrow, UInt ex_ncol,const vector<Index_t> &bindx,const vector<Index_t> &ybind );
  MSRPatt(const MSRPatt &RightHandMSRP);
  MSRPatt(const CSRPatt &RightHandCSRP);
  MSRPatt & operator= (const MSRPatt&  RhMsr);
  template<typename DOF> MSRPatt(DOF const  & dof, UInt const nbcomp=1);

  // Miguel 12/2003: new version handling patterns coming from IP stabilization
  template<typename DOF, typename MESH>
    MSRPatt(const  DOF& dof, const MESH& mesh, const UInt nbcomp);

  template<typename DOF>
    bool buildPattern(DOF const  & dof, UInt const nbcomp);

  // Miguel 12/2003: new version handling patterns coming from IP stabilization
  template<typename DOF, typename MESH>
    bool buildPattern(const  DOF& dof, const MESH& mesh, const UInt nbcomp);

  Container  bindx() const {return _bindx;};
  Container  ybind() const {return _ybind;};
  Index_t * giveRaw_bindx() {return &(_bindx.front());}; // Give _bindx (in a raw form)
  Index_t * giveRaw_ybind() {return &(_ybind.front());}; // Give _ybind (in a raw form)
  Index_t const * giveRaw_bindx() const {return &(_bindx.front());}; // Give _bindx (in a raw form)
  Index_t const * giveRaw_ybind() const {return &(_ybind.front());}; // Give _ybind (in a raw form)
  Container & give_bindx() {return _bindx;};// Give _bindx (as container)
  Container & give_ybind() {return _ybind;};// Give _ybind (as container)
  Container const & give_bindx() const {return _bindx;};// Give _bindx (as container)
  Container const & give_ybind() const {return _ybind;};// Give _ybind (as container)
//
  inline UInt nbNeighbours(ID const d)const;
  inline ID neighbour(ID const i,ID const d) const;
  void neighbours(ID const d, Container & start) const;
  template <typename Iter> inline UInt row(Diff_t const row, Iter  coldata, Iter position) const;// extracts a row (useful to implement A*b)

  inline pair<Diff_t,bool> locate_dof(ID const i,ID const j) const;
//
  inline pair<UInt,UInt> giveMinMax() const;
  inline pair<Diff_t,bool> locate_index(Index_t const i,Index_t const j) const;

  void showMe(bool verbose=false, ostream& c=cout) const ; // pattern visualization
  void spy(string const & filname="matrice.m")const ; //pattern visualization a la Matlab

  // Construction of a diagonal matrix of n blocks
  friend void diagblockMatrix(MSRPatt &ans, MSRPatt const &patt, UInt const nblock);

 protected:
  Diff_t _row_off(Index_t i)const {// Row offset for row index i
    return _i2o(_bindx[_i2o(i)]);}
  pair<Diff_t,bool> locate_pattern(Index_t const i,Index_t const j) const;
  inline bool isThere(Index_t i,Index_t j)const ; // superata dalla locate

 private:
  Container _bindx;
  Container _ybind; // WARNING: AV January 2001: it helps in reading the matrix by columns ;
  // a counterpart of ybind should be added also to the CSR format
  // ybind should facilitate also the locate_pattern subroutines
  // IT IS BASED ON THE ASSUMPTION THAT THE PATTERN IS SYMMETRIC

  template<typename DOF>
  void _buildPattern(DOF const & dof, DynPattern const & dynpatt,
                     UInt const nbcomp);
};

///////////////////////////////////////////////////
///
///      Mixed Local Patterns
//////////////////////////////////////////////////
//
/* This class mixedLocalPattern<FE1, FE2> is used to provide the local
  pattern for two BasisFInite element, it corresponds to the pattern
  required to build element matrices of the type $m_{ij}=\int_K
  \phi_i^1\phi^2_j$. The number of Rows correspond to the number of local
  degrees of freedon for the finite elemnt of type 1, i.e FE1::nbNode,
  while the number of colums are FE2::nbNode. The major assumption is that
  $\forall K$ and $i,j\in K$, the support of $\phi_i^1\vert_K$ and
  $\phi_j^2\vert_K$ has non-null intersection. In other words, the local
  matrix has a  full pattern.
  The pattern follows the rule of numbering
  */
template<typename FE1, typename FE2>
class MixedLocalPattern: public PatternDefs
{
public:
  UInt nRows() const;
  UInt nCols() const;
  UInt nbPattern() const;
  UInt patternFirst(UInt const i) const; //Numbering from 0
  UInt patternSecond(UInt const i) const; //Numbering from 0
};
///////////////////////////////////////////////////
///
///      Mixed Patterns
///
//////////////////////////////////////////////////
/*
  The mixed pattern class is a very general class able to held
  multiblock matrix patterns.
  The block numbering starts from (0,0). Each block contains a pointer to a Pattern class and a couple of
  offsets, which indicate how the LOCAL block row/cols numbering has to be increased to get the
  GLOBAL numbering (i.e. the numbering associated to to the global matrix).
  It may be used in two forms
  1) as a viewer, then the local patterns are contructued externally and then "linked"
  to the mixed patter object, or
  2) by delegating to contruction of  the local patterns to the mixed pattern object.
  In the first case the destruction of the mixed pattern object will NOT imply the destruction of the
  local patterns, while in the other case everithing is destroied.
*/
template <UInt BROWS, UInt BCOLS, typename PATTERN=CSRPatt>
class
MixedPattern : public PatternDefs
{
public:

  MixedPattern();

  // Miguel 11/02: I want to construct a diagonal pattern
  // in the initialisation list, so I need a constructor which makes the pattern
  // construct and link to an external pattern (type = "diag") (Alain: type="full")
  MixedPattern(PATTERN & ex_patt, const string& type="full");

  ~MixedPattern();

  // make a diagonal pattern for vectorial problem
  void makeDiagPattern(PATTERN & ex_patt);

  inline pair<UInt,UInt> nBlocks() const; // Number of blocks (rows and columns)

  inline UInt nRows(Diff_t const m, Diff_t const n) const; // Number of rows in block (m,n)
  inline UInt nCols(Diff_t const m, Diff_t const n) const;
  inline UInt nNz(Diff_t m, Diff_t n) const; // Non zeros on Block (m,n)
  inline UInt nRows() const; // Global number of rows
  inline UInt nCols() const; // Global number of cols.
  UInt nNz() const; // Non zeros in global matrix

  // Neighbours at block level. (local ID numbering)
  UInt nbNeighbours(Diff_t const m, Diff_t const n, ID const d) const ;
  ID neighbour(Diff_t const m, Diff_t const n,ID const i, ID const d) const;
  inline void neighbours(Diff_t const m, Diff_t const n, ID const d,Container & neighs) const;
  template <typename Iter>
  inline UInt row(Diff_t const m, Diff_t const n,Diff_t const row, Iter  coldata, Iter Position) const;// extracts a row (useful to implement A*b)

  // Neighbours at global level (global ID numbering)
  UInt nbNeighbours(ID const d_g) const ;
  ID neighbour(ID const i_g, ID const d_g) const;
  void neighbours(ID const d_g, Container & neighs) const;
  template <typename Iter>
  inline UInt row(Diff_t const row, Iter  coldata, Iter position) const;// extracts a row (useful to implement A*b)

  inline       PATTERN * block_ptr(Diff_t const m, Diff_t const n); // Pointer to a a local pattern
  inline const PATTERN * block_ptr(Diff_t const m, Diff_t const n) const; // Pointer to a a local pattern

  inline pair<UInt,UInt> blockOffset(UInt const m, UInt const n) const;// The row/col offsets of the block
  pair<UInt,UInt> locateElBlock(Index_t const i_g, Index_t const j_g) const;
  //  Gives the block numbering corresponding to the  GLOBAL matrix index (i_g,j_g) Returns
  // (BROWS,BCOLS) if element not found
  pair<Diff_t,Diff_t> locateDofBlock(ID const di_g, ID const dj_g) const;
  // Give the block correponding to a THE GLOBAL DOF (di_g,dj_g)
  // Returns (BROWS,BCOLS) if dofs not found

  // local/global numbering in the block, given global numbering. It can be applied to indices and IDs.
  // I rely on implicit conversion ID->Index_t if Index_t != ID
  // LUCA: If that does not work I will do a template function pair<A,A>localNumber<T>(Diff_t,Diff_t,T,T) and the
  // necessary specialisations.
  inline pair<ID,ID> localNumber(Diff_t const m, Diff_t const n, ID const i_g, ID const j_g) const;
  inline pair<ID,ID> globalNumber(Diff_t const m, Diff_t const n, ID const i, ID const j) const;

  // Returns position in matrix  corrspondinf to a  DOF (local numbering)  in a block
  pair<Diff_t,bool> locateDof(Diff_t const m, Diff_t const n, ID const di, ID const dj) const;

  // Returns position in matrix corrspondinf to element (i,j) in a block (local numbering)
  pair<Diff_t,bool> locateIndex(Diff_t const m, Diff_t const n, Index_t const i, Index_t const j)const ;

  // I can set the size of a block without linking the block to a pattern (I may have a matrix of ZEROS!)
  // To that purpose I use the following:
  void buildZeroPattern(Diff_t const m, Diff_t const n,UInt const nrows, UInt const ncols);

  // This just links the block to an existing pattern (Viewer paradigm)
  void linkBlockToPattern(UInt const m, UInt const n,PATTERN & pattern);

  // These insetade create a pattern, by relying on the local pattern constructors
  template <typename DOF1, typename DOF2> PATTERN * buildBlock(Diff_t m, Diff_t n, DOF1 dof1, DOF2 dof2);
  template <typename DOF1> PATTERN * buildBlock(Diff_t m, Diff_t n, DOF1 dof1);

  // if two blocks have the same pattern i may avoid repeating it by linking
 // the two block. The "from" block MUST have already been set  to a local pattern.
  void linkBlocks(Diff_t const m_from, Diff_t const n_from, Diff_t const m_to, Diff_t const n_to);

  // Use with extreme care. Itdoes not delete if block is linked (viewer paradigm).
  void deleteBlock(Diff_t const m, Diff_t const n);

  // Checking and probing..
// This is true if a block is set but not linked to a pattern (zero matrix)
  inline bool isZero(Diff_t const m, Diff_t const n) const;

  // Tre is a block has been set to a  local pattern
  inline bool isSet(Diff_t const m, Diff_t const n) const;
  void showMe( bool verbose=false, ostream & c=cout) const;
  // Tests if offsets are  consistent with that of a global matrix
  bool check(bool verbose=false, ostream & c=cout) const;

protected:

private:
  void resetOffset(Diff_t const m=0, Diff_t const n=0);


  PATTERN * _blocks[BROWS][BCOLS];
  bool _linked[BROWS][BCOLS]; // indicates whether a block is actually a link to another block
  UInt _rowoff[BROWS][BCOLS]; // rows offset
  UInt _coloff[BROWS][BCOLS];// cols offset
  UInt _nrows[BROWS][BCOLS]; // rows in block
  UInt _ncols[BROWS][BCOLS]; // cols in block

};



// U T I L I T I E S ===========================================
//
// A stupid global function
//UInt addOne(UInt i){return i+1;}


// Hacking of a stl routine
template <class ForwardIterator, class T>
ForwardIterator search_binary(ForwardIterator first, ForwardIterator last, const T& value)
{
  ForwardIterator i = lower_bound(first, last, value);
  return (i != last && !(value < *i)) ? i:last;
}

// This function is used to extract pairs. TODO: Should go somewhere else!

template
<typename T1, typename T2>
inline
void
extract_pair(pair<T1,T2> const & p, T1 & t1, T2 & t2)
{
  t1=p.first;
  t2=p.second;
}

//=========================================================================================
//I M P L E M E N T A T I O N S ===========================
//=====================================================================================

inline
PatternDefs::Index_t PatternDefs::_d2i(ID const d) const {return d-1+PatternOffset;}

inline
ID   PatternDefs::_i2d(Index_t const i) const {return i+1-PatternOffset;}

inline
PatternDefs::Diff_t PatternDefs::_i2o(Index_t const i) const {return i-PatternOffset;} // From index to offsets

inline
PatternDefs::Diff_t PatternDefs::_d2o(ID const d) const {return d-1;}

inline
PatternDefs::Container & PatternDefs::_i2d(Container & list_of_indices)const{
  for (Container::iterator ip=list_of_indices.begin();ip!=list_of_indices.end();++ip) *ip=_i2d(*ip);
  return list_of_indices;
}

inline
PatternDefs::Container & PatternDefs::_d2i(Container & list_of_dof)const{
  for (Container::iterator ip=list_of_dof.begin();ip!=list_of_dof.end();++ip) *ip=_d2i(*ip);
  return list_of_dof;
}

inline
PatternDefs::Container & PatternDefs::_d2o(Container & list_of_dof)const{
  for (Container::iterator ip=list_of_dof.begin();ip!=list_of_dof.end();++ip) *ip=_d2o(*ip);
  return list_of_dof;
}


//////////////////////////////////////////////////////////////
// BasePattern
//////////////////////////////////////////
inline
UInt
BasePattern::nRows()const
{
  return _nrows;
}

inline
UInt
BasePattern::nCols()const
{
  return _ncols;
}

inline
UInt
BasePattern::nNz()const
{
  return _nnz;
}

inline
bool
BasePattern::isEmpty()const
{
  return ! _filled;
}

// Alain (nov. 2002). version for (nbcomp X nbcomp) block matrix
template<typename DOF1>
bool BasePattern::setpatt(DOF1 const & dof1, DynPattern & dynpattern,
			  UInt const nbcomp){

// The loop is on the elements rather than on the volumes:
// it should work also in 2D
  Diff_t ig,jg;

  for (UInt el=1; el<= dof1.numElements(); el++){
// We do not store the diagonal
// So inner loop starts from nbDiag:
// indeed the first nbDiag entries of patternFirst
// and patternSecond involve diagonal guys

   for (int l=dof1.fe.nbDiag(); l<dof1.fe.nbPattern(); ++l){
      ID i  = dof1.fe.patternFirst(l)+1;
      ID j  = dof1.fe.patternSecond(l)+1;

      ig = dof1.localToGlobal(el,i)-1; // I store with numering from 0 (and correct later on)
      jg = dof1.localToGlobal(el,j)-1;
	dynpattern.insert(setBareEdge(ig, jg));
    }
  }

  _nrows = nbcomp*dof1.numTotalDof();
  _ncols = _nrows; // just one DOF => square matrix
  _nnz = nbcomp*nbcomp*( dof1.numTotalDof() + 2*( dynpattern.size() ) );

  return true;
  }

//
// Miguel 12/2003
//
template<typename DOF, typename MESH>
  bool BasePattern::setpatt(const DOF& dof, const MESH& mesh, DynPattern & dynpattern,
			    UInt const nbcomp){


  Diff_t ig,jg;
  UInt idF, iElAd1, iElAd2, iFaEl1, iFaEl2;

  for (UInt el=1; el<= dof.numElements(); el++){

    // standard neighboors
    //
    for (int l=dof.fe.nbDiag(); l<dof.fe.nbPattern(); ++l){
      ID i  = dof.fe.patternFirst(l)+1;
      ID j  = dof.fe.patternSecond(l)+1;

      ig = dof.localToGlobal(el,i)-1;
      jg = dof.localToGlobal(el,j)-1;
      dynpattern.insert(setBareEdge(ig, jg));
    }
  }


  //
  // The following lines of code should be improved !!
  //
  // for each face the local numbering of the neighboors connected to this face
  //
  UInt p1[]={4,3,1,2};
  UInt p2[]={4,8,9,10,
	     3,6,7,10,
	     1,5,7,8,
	     2,5,6,9};
  UInt* a;
  UInt iop,jop,nop; // number of opposite dof

  if ( dof.fe.nbLocalDof == 4) {
    a = p1; // P1 FE
    nop = 1;
  }
  else if ( dof.fe.nbLocalDof == 10) {
    a = p2; // P2 FE
    nop = 4;
  }
  else
    ERROR_MSG("Sorry, IP stabilization only works for P1 or P2 FEM");

  //
  // non-standard neighboors
  //
  for (idF=mesh.numBFaces()+1; idF<=mesh.numFaces(); ++idF){

    iElAd1 = mesh.face(idF).ad_first();
    iFaEl1 = mesh.face(idF).pos_first();
    iElAd2 = mesh.face(idF).ad_second();
    iFaEl2 = mesh.face(idF).pos_second();
    for (iop=1; iop <= nop ; ++iop) {
      ig = dof.localToGlobal(iElAd1,a[(iFaEl1-1)*nop+iop-1])-1;
      for (jop=1; jop <= nop ; ++jop) {
	jg = dof.localToGlobal(iElAd2,a[(iFaEl2-1)*nop+jop-1])-1;
	dynpattern.insert(setBareEdge(ig, jg));
      }
    }
  }
  _nrows = nbcomp*dof.numTotalDof();
  _ncols = _nrows;
  _nnz = nbcomp*nbcomp*( dof.numTotalDof() + 2*( dynpattern.size() ) );

  return true;
}


/* THE ROUTINE THAT BUILDS THE DYN PATTERN MAY BE GLOBAL (WHY NOT) */

template<typename DOF1, typename DOF2>
bool
  BasePattern::setpatt(DOF1 const & dof1, DOF2 const & dof2,
		       DynPattern & dynpatt, UInt const bRows, UInt const bCols){

  UInt nelem=dof1.numElements();
  UInt ig,jg;
  BareEdge _be;
  for (UInt el=1; el<= nelem; el++){
    for (int l1 =0; l1 <  dof1.fe.nbLocalDof; ++l1){
      UInt i = dof1.fe.patternFirst(l1)+1; // row
      for (int l2 =0; l2 <  dof2.fe.nbLocalDof; ++l2){

	UInt j = dof2.fe.patternSecond(l2)+1; //cols

	ig=dof1.localToGlobal(el,i)-1; // I store with numering 0 always and then I correct
	jg=dof2.localToGlobal(el,j)-1;

	_be.first=ig;
	_be.second=jg;
	dynpatt.insert(_be);
      }
    }
  }
  _nrows=bRows*dof1.numTotalDof();
  _ncols=bCols*dof2.numTotalDof();
  _nnz=bRows*bCols*dynpatt.size();

  return true;
}

////////////////////////////////////////////////////////////////////////
//
// C S R Pattern
//
////////////////////////////////////////////////////////////////////////
template<typename DOF1>
CSRPatt::CSRPatt(DOF1 const & dof1, UInt const nbcomp) {
#ifdef TEST_PRE
  bool built = buildPattern(dof1, nbcomp);
  ASSERT_PRE(built,"Error in CSR Pattern construction from DOF object");
#else
  bool built = buildPattern(dof1, nbcomp);
#endif
};

// Alain (nov. 2002): version for (nbcomp X nbcomp) block matrix.
template<typename DOF1>
bool CSRPatt::buildPattern(DOF1 const & dof1, UInt const nbcomp)
{

  // This is the real constuctor.
  // It builds a FE type pattern from A SINGLE FE DOF, in CSR format;
  // I change to the actual index (i.e. I take into account of
  // PatternOffset) only at the end.

  DynPattern dynpatt;

  bool built = setpatt(dof1, dynpatt, nbcomp);

  if (! built) {
    _filled=false;
    return built;
  }
  
  Diff_t ig,jg,cur;

  // I use a modified version of CSR (compatible with the standard one), where:
  // 1) ia is dimensioned nrows+1, so that I do not have to treat the
  // last term in a special way
  // 2) The entries in ja are sorted so that the first entry for each row is the
  //    diagonal term and the other entries are in ascending order.
  //
  _ia.resize(_nrows+1,nbcomp); // "nbcomp" stands for the diagonal terms of all blocks
  // Count

  for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    jg=d->second;
    _ia[ig]+= nbcomp;
    _ia[jg]+= nbcomp;
  }

  //other components:
  Diff_t offset= dof1.numTotalDof();
  UInt icomp, jcomp;

  Container::iterator start     = _ia.begin();
  Container::iterator end       = _ia.begin()+ offset;
  for (icomp=1; icomp < nbcomp; icomp++)
    {
      Container::iterator start_copy= _ia.begin()+ icomp*offset;
      copy(start, end, start_copy);
    }

  // shift right 1 position and set _ia[0]
  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;

  // Count entries BY ACCUMULATING
  for (ig=1; ig <_nrows; ++ig) _ia[ig]+=_ia[ig-1];
  // >>>>  Now _ia[i] should point to the first available position in _ja <<<

  _ja.resize(_nnz);
  _jaT.resize(_nnz);

  // for each component
  for (icomp=0; icomp< nbcomp; icomp++){

    // We put the diagonal first
    for (ig=0; ig < offset; ++ig){

      _jaT[_ia[ig+icomp*offset]]= _ia[ig+icomp*offset];
      _ja[_ia[ig+icomp*offset]++]= ig+icomp*offset;
    }

    for (jcomp=0; jcomp< nbcomp; jcomp++){

      // >>>>  Now _ia[i] should point to the first available position in _ja <<<

      //We add the missing diagonal term for jcomp != icomp
      if (jcomp != icomp)
	for (ig=0; ig < offset; ++ig){
	  _jaT[_ia[ig+icomp*offset]]= _ia[ig+jcomp*offset];
	  _ja[_ia[ig+icomp*offset]++]=ig+jcomp*offset;
	}

      // Now the rest. The pattern is symmetric, but traditionally the CSR Pattern no (see CSRPattSymm)
      for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
	ig=d->first + icomp*offset;
	jg=d->second + jcomp*offset;
	_jaT[_ia[ig]] = _ia[jg];
	_jaT[_ia[jg]] = _ia[ig];
	_ja[_ia[ig]++] = jg;
	_ja[_ia[jg]++] = ig;
      }
    }
  }

  // shift right 1 position and set _ia[0]
  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;

  for (ig=0; ig <_nrows; ++ig){
    // We now sort the off diagonal entries
    jg=_ia[ig]+1; // I want the diagonal first!
    cur = _ia[ig+1];
    sort(_ja.begin()+jg,_ja.begin()+cur);
  }
  // Now I take account of the PatternOffset
  if ( PatternOffset != 0){
    for (Container::iterator ip=_ia.begin(); ip != _ia.end(); ++ip) *ip+=PatternOffset;
    for (Container::iterator ip=_ja.begin(); ip != _ja.end(); ++ip) *ip+=PatternOffset;
    for (Container::iterator ip=_jaT.begin(); ip != _jaT.end(); ++ip) *ip+=PatternOffset;
  }

  _filled=true;
  _diagfirst=true;

  dynpatt.clear(); // not sure if this helps...

  return built;
};

template<typename DOF1,typename DOF2>
CSRPatt::CSRPatt(DOF1 const & dof1, DOF2 const & dof2,
		 UInt const bRows, UInt const bCols) {
  bool built;
  built = buildPattern(dof1, dof2, bRows, bCols);
  ASSERT_PRE(built,"Error in CSR Pattern construction from DOF object");
};

template<typename DOF1,typename DOF2>
bool CSRPatt::buildPattern(DOF1 const & dof1, DOF2 const & dof2,
			   UInt const bRows, UInt const bCols){
  Diff_t ig,jg,cur;

  // ASSERT that the num of elements are the same in the 2 dof: we are talking of
  // the same mesh, aren't we?
  ASSERT_PRE(dof1.numElements()==dof2.numElements(), "Cannot run buildpattern on DOFs of different meshes");

  DynPattern dynpatt;

  bool built=setpatt(dof1,dof2,dynpatt,bRows,bCols);
  if(!built)return built;

  _ia.resize(_nrows+1,0);

  // The dynamic pattern now contains all (i,j) since the
  // matrix may be non-diagonal
  for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    _ia[ig]+= bCols;
  }

  //other components:
  Diff_t row_offset= dof1.numTotalDof();
  Diff_t col_offset= dof2.numTotalDof();
  UInt icomp, jcomp;

  Container::iterator start     = _ia.begin();
  Container::iterator end       = _ia.begin()+ row_offset;
  for (icomp=1; icomp < bRows; icomp++)
    {
      Container::iterator start_copy= _ia.begin()+ icomp*row_offset;
      copy(start, end, start_copy);
    }

  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;
  // ACCUMULATE
  for (ig=1; ig <_nrows; ++ig)_ia[ig]+=_ia[ig-1];
  // >>>>  Now _ia[i] should point to the first available position in _ja <<<

  _ja.resize(_nnz);

  //for each component
  for (icomp=0; icomp < bRows; icomp++){
    for (jcomp=0; jcomp < bCols; jcomp++){

      for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
	ig=d->first + icomp*row_offset;
	jg=d->second + jcomp*col_offset;
	_ja[_ia[ig]++] = jg;
      }
    }
  }

  // shift right 1 position and set _ia[0]
  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;

  for (ig=0; ig <_nrows; ++ig){
    // We now sort all entries
    jg=_ia[ig];
    cur = _ia[ig+1];
    sort(_ja.begin()+jg,_ja.begin()+cur);
  }
  // Now I take account of the PatternOffset
  if ( PatternOffset != 0){
    for (Container::iterator ip=_ia.begin(); ip != _ia.end(); ++ip) *ip+=PatternOffset;
    for (Container::iterator ip=_ja.begin(); ip != _ja.end(); ++ip) *ip+=PatternOffset;
  }
  _filled=true;
  _diagfirst=false;
  dynpatt.clear(); // not sure if this helps...
  return true;
}

inline
UInt CSRPatt::nbNeighbours(ID const d) const {
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");
  return _ia[d]-_ia[d-1];
}

inline
ID CSRPatt::neighbour(ID const i, ID const d) const {
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_BD(i>0);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");
  return _i2d(_ja[_i2o(_ia[d-1]) +_d2o(i)]);
};

inline
void CSRPatt::neighbours(ID const d, Container & neigs) const{
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");
  Container::const_iterator start1=_ja.begin()+_i2o(_ia[d-1]);
  //neigs.clear();
  neigs.resize(nbNeighbours(d),0);// copy wants a sequence!
  copy(start1,start1+nbNeighbours(d),neigs.begin());
  _i2d(neigs);// trasform to DOF ID's
};

template <typename Iter>
inline
UInt
CSRPatt::row(Diff_t const row, Iter  coldata, Iter position) const{
  // extracts a row (useful to implement A*b)
  // to make things faster no controls are made!!!
  Container::const_iterator start=_ja.begin()+_i2o(_ia[row]);
  Container::const_iterator end=_ja.begin()+_i2o(_ia[row+1]);
  copy(start,end,coldata);
  for (Index_t i=0; i< _ia[row+1]-_ia[row];++i) position[i]=_i2o(_ia[row])+i;
  return _ia[row+1]-_ia[row];
}

inline
bool CSRPatt::isThere(Index_t const i, Index_t const j) const {
  ASSERT_BD(i>=PatternOffset && i<static_cast<Index_t>(_nrows)+PatternOffset);
  ASSERT_BD(j>=PatternOffset && j<static_cast<Index_t>(_ncols)+PatternOffset);
  if (! _filled ) return false;

  Container::const_iterator start=_ja.begin()+_row_off(i);
  if (*start == j)return  true; // Diagonal term
  return binary_search(++start, _ja.begin()+_row_off(i+1), j);
}


// locate function for CSR Pattern. It returns a pair. First member
// is the position (offset) in the array correponding to (i,j), the second
// is a bool telling if that position exists (i.e. if i,j is in the
//pattern). If the boolean value is false the first member is meaningless!
inline
pair<PatternDefs::Diff_t,bool>
CSRPatt::locate_index(Index_t const i,Index_t const j) const {
  return locate_pattern(i,j);
}

inline
pair<PatternDefs::Diff_t,bool>
CSRPatt::locate_dof(ID const i,ID const j) const {
    return locate_pattern(_d2i(i),_d2i(j));
}


inline
pair<UInt,UInt>
CSRPatt::giveMinMax() const {
  Container::const_iterator current;
  UInt curr_min=_ncols;
  UInt curr_max=0;
  for (current=_ia.begin()+1;current<_ia.end();++current){
    UInt loc_col = *current - *(current-1);
    curr_min=min(curr_min,loc_col);
    curr_max=max(curr_max,loc_col);
    //if (loc_col < curr_min) curr_min=loc_col;
      //if (loc_col > curr_max) curr_max=loc_col;
  }
  return make_pair(curr_min,curr_max);
}

// Version for the construction two patterns: patt and its transpose one,
// in addition an other Container for the access to
// the columns of the transpose matrix is built.
// Alain. nov 2002.
template<typename DOF1,typename DOF2>
void buildPattTpatt(DOF1 const& dof1, DOF2 const  & dof2,
		    CSRPatt &Patt, CSRPatt &Tpatt,
		    UInt const bRows, UInt const bCols)
{
  PatternDefs::Diff_t ig,jg,cur;

  // ASSERT that the num of elements are the same in the 2 dof:
  // we are talking of the same mesh, aren't we?
  ASSERT_PRE(dof1.numElements()==dof2.numElements(), "Cannot run buildpattern on DOFs of different meshes");

  PatternDefs::DynPattern dynpatt;

  bool built= Patt.setpatt(dof1,dof2,dynpatt,bRows,bCols);
  ASSERT(built, "Cannot build pattern");
  // instead of using setpatt, Tpatt is directly given as the transpose of
  // Patt:
  Tpatt._nrows= bCols*dof2.numTotalDof();
  Tpatt._ncols= bRows*dof1.numTotalDof();
  Tpatt._nnz= Patt._nnz;

  Patt._ia.resize(Patt._nrows+1,0);
  Tpatt._ia.resize(Tpatt._nrows+1,0);

  // For Patt and Tpatt
  // cols of Patt are rows of Tpatt

  // The dynamic pattern now contains all (i,j) since the
  // matrix may be non-diagonal
  for (PatternDefs::DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    jg=d->second;
    Patt._ia[ig]+= bCols;
    Tpatt._ia[jg]+= bRows;
  }

  //other components:
  PatternDefs::Diff_t row_offset= dof1.numTotalDof();
  PatternDefs::Diff_t col_offset= dof2.numTotalDof();
  UInt icomp, jcomp;

  //Patt
  PatternDefs::Container::iterator start     = Patt._ia.begin();
  PatternDefs::Container::iterator end       = Patt._ia.begin()+ row_offset;
  for (icomp=1; icomp < bRows; icomp++)
    {
      PatternDefs::Container::iterator start_copy=
	Patt._ia.begin()+ icomp*row_offset;
      copy(start, end, start_copy);
    }
  //Tpatt
  start     = Tpatt._ia.begin();
  end       = Tpatt._ia.begin()+ col_offset;
  for (jcomp=1; jcomp < bCols; jcomp++)
    {
      PatternDefs::Container::iterator start_copy=
	Tpatt._ia.begin()+ jcomp*col_offset;
      copy(start, end, start_copy);
    }

  rotate(Patt._ia.begin(),Patt._ia.end()-1,Patt._ia.end());
  Patt._ia[0]=0;
  rotate(Tpatt._ia.begin(),Tpatt._ia.end()-1,Tpatt._ia.end());
  Tpatt._ia[0]=0;
  // ACCUMULATE
  for (ig=1; ig <Patt._nrows; ++ig)Patt._ia[ig]+=Patt._ia[ig-1];
  for (ig=1; ig <Tpatt._nrows; ++ig)Tpatt._ia[ig]+=Tpatt._ia[ig-1];
  // >>>>  Now _ia[i] should point to the first available position in _ja <<<

  Patt._ja.resize(Patt._nnz);
  Tpatt._ja.resize(Tpatt._nnz);

  Patt._jaT.resize(Patt._nnz);
  Tpatt._jaT.resize(Tpatt._nnz);

  //for each component
  for (icomp=0; icomp < bRows; icomp++){
    for (jcomp=0; jcomp < bCols; jcomp++){
      for (PatternDefs::DynPattern::iterator d=dynpatt.begin();
	   d !=dynpatt.end(); ++d){
	ig=d->first + icomp*row_offset;
	jg=d->second + jcomp*col_offset;
	Patt._jaT[Patt._ia[ig]] = Tpatt._ia[jg];
	Tpatt._jaT[Tpatt._ia[jg]] = Patt._ia[ig];
	Patt._ja[Patt._ia[ig]++] = jg;
	Tpatt._ja[Tpatt._ia[jg]++] = ig;
      }
    }
  }

  // shift right 1 position and set _ia[0]
  rotate(Patt._ia.begin(),Patt._ia.end()-1,Patt._ia.end());
  Patt._ia[0]=0;
  rotate(Tpatt._ia.begin(),Tpatt._ia.end()-1,Tpatt._ia.end());
  Tpatt._ia[0]=0;

  for (ig=0; ig <Patt._nrows; ++ig){
    // We now sort all entries
    jg=Patt._ia[ig];
    cur = Patt._ia[ig+1];
    sort(Patt._ja.begin()+jg,Patt._ja.begin()+cur);
  }

  for (ig=0; ig <Tpatt._nrows; ++ig){
    // We now sort all entries
    jg=Tpatt._ia[ig];
    cur = Tpatt._ia[ig+1];
    sort(Tpatt._ja.begin()+jg,Tpatt._ja.begin()+cur);
  }

  // Now I take account of the PatternOffset
  if ( PatternOffset != 0){
    for (PatternDefs::Container::iterator ip=Patt._ia.begin(); ip != Patt._ia.end(); ++ip) *ip+=PatternOffset;
    for (PatternDefs::Container::iterator ip=Patt._ja.begin(); ip != Patt._ja.end(); ++ip) *ip+=PatternOffset;
    for (PatternDefs::Container::iterator ip=Tpatt._ia.begin(); ip != Tpatt._ia.end(); ++ip) *ip+=PatternOffset;
    for (PatternDefs::Container::iterator ip=Tpatt._ja.begin(); ip != Tpatt._ja.end(); ++ip) *ip+=PatternOffset;
  }
  Patt._filled=true;
  Patt._diagfirst=false;
  Tpatt._filled=true;
  Tpatt._diagfirst=false;
  dynpatt.clear(); // not sure if this helps...
}


// column-concatenation of two CSR block patterns
CSRPatt colUnify(CSRPatt const &patt1, CSRPatt const &patt2);


// column-concatenation of one CSR block patterns and ncolZero null columns
// zero on the right
CSRPatt colUnify(CSRPatt const &patt1, UInt const ncolZero);


//zero on the left
CSRPatt colUnify(UInt const ncolZero, CSRPatt const &patt1);


// rows-concatenation of two CSR block patterns
CSRPatt rowUnify(CSRPatt const &patt1, CSRPatt const &patt2);

// row-concatenation of one CSR block patterns and nrowZero null rows
// zero on the below
CSRPatt rowUnify(CSRPatt const &patt1, UInt const nrowZero);

// zero on the top
CSRPatt rowUnify(UInt const nrowZero, CSRPatt const &patt1);

// Costruzione di una Matrice diagonale a blocchi
CSRPatt diagblockMatrix(CSRPatt const &patt, UInt const nblock);

////////////////////////////////////////////////////////////////////////
//
// V B R Pattern
//
////////////////////////////////////////////////////////////////////////
template<typename DOF1>
VBRPatt::VBRPatt(DOF1 const & dof1, UInt const blockSize):
  CSRPatt(dof1)
{
  bool built= buildPattern(blockSize);
  ASSERT_PRE(built,"Error in VBR Pattern construction from DOF object");
}


template<typename DOF1,typename DOF2>
VBRPatt::VBRPatt(DOF1 const  & dof1,DOF2 const  & dof2, UInt
		 const blockSize): CSRPatt(dof1,dof2)
{
  bool built = buildPattern(blockSize);
  ASSERT_PRE(built,"Error in VBR Pattern construction from DOF object");
};

inline UInt VBRPatt::row(Diff_t const row,  Container & coldata, Container & position) const
{
  // extracts a row (useful to implement A*b)
  // to make things faster no controls are made!!!
  // works on the elements and not the blocks !

  UInt blsize=_rpntr[1]-_rpntr[0]; // size of square block
  Diff_t blockr=rbloc(row);// block row offset in which is row
  Diff_t localr=locr(row);  // local row number into the block

  Diff_t start = _i2o(_ia[blockr]);
  Diff_t end = _i2o(_ia[blockr+1]);
  Container::const_iterator istart=_ja.begin() + start;//initial block number
  Container::const_iterator iend=_ja.begin() + end;//end block number(excluded)

  Container colblocdata(end-start); // column block number
  copy(istart,iend,colblocdata.begin());

  for (Diff_t ibl=start; ibl< end; ++ibl){
    for (UInt i=0; i<blsize; ++i){
      coldata.push_back(colblocdata[ibl-start]*blsize+i);
      position.push_back(_indx[ibl]+localr+i*blsize);
    };
  };

  return (end-start)*blsize;
};
////////////////////////////////////////////////////////////////////////
//
// C S R Symmetric Pattern
//
////////////////////////////////////////////////////////////////////////

template<typename DOF1>
CSRPattSymm::CSRPattSymm(DOF1 const & dof1) {
  DynPattern dynpatt;
  bool built = buildPattern(dof1);
  ASSERT_PRE(built,"Error in CSRSymm Pattern construction from DOF object");
}


template<typename DOF1>
bool CSRPattSymm::buildPattern(DOF1 const & dof1){

  DynPattern dynpatt;
  bool built = setpatt(dof1,dynpatt);
  if (! built) {
    _filled=false;
    return built;
  }

  UInt ig,jg,cur;
  UInt nnz_symm = _nrows + dynpatt.size(); // actual size of _ja when exploiting the symmetry:

  ////////////////////////////////////// Please, observe that nnz_symm = (_nnz-_nrows)/2

  // I use a modified version of CSR (compatible with the standard one), where:
  // 1) ia is dimensioned nrows+1, so that I do not have to treat the
  // last term in a special way
  // 2) The entries in ja are sorted so that the first entry for each row is the
  //    diagonal term and the other entries are in ascending order.
  //
  _ia.resize(_nrows+1,Index_t(1)); // "Index_t(1)" stands for the diagonal terms
  //  Container temp(_nrows+1,Index_t(1));
  // Count
  for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    jg=d->second;
    ++_ia[ig];
    //    ++temp[ig];
    //    ++temp[jg];
  }

  // shift right 1 position and set _ia[0]
  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;
  // Count entries BY ACCUMULATING
  for (ig=1; ig <_nrows; ++ig){
    _ia[ig]+=_ia[ig-1];
  }

  // >>>>  Now _ia[i] should point to the first available position in _ja <<<
  _ja.resize(nnz_symm);

  // We put the diagonal first
  for (ig=0; ig <_nrows; ++ig) _ja[_ia[ig]++]=ig;


  // >>>>  Now _ia[i] should point to the first available position in _ja <<<
  // Now the rest. The pattern is symmetric
  for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    jg=d->second;
    _ja[_ia[ig]++]=jg;
  }

  // shift right 1 position and set _ia[0]
  rotate(_ia.begin(),_ia.end()-1,_ia.end());
  _ia[0]=0;


  for (ig=0; ig <_nrows; ++ig){
    // We now sort the off diagonal entries
    jg=_ia[ig]; // I am storing the upper diagonal terms: it is safe starting from the diagonal
    cur = _ia[ig+1];
    sort(_ja.begin()+jg,_ja.begin()+cur);
  }

  // Now I take account of the PatternOffset
  if ( PatternOffset != 0){
    for (Container::iterator ip=_ia.begin(); ip != _ia.end(); ++ip) *ip+=PatternOffset;
    for (Container::iterator ip=_ja.begin(); ip != _ja.end(); ++ip) *ip+=PatternOffset;
  }

  _filled=true;
  _diagfirst=true;
  dynpatt.clear(); // not sure if this helps...
  return true;
};

inline
pair<PatternDefs::Diff_t,bool>
CSRPattSymm::locate_index(Index_t const i,Index_t const j) const {
  return locate_pattern(i,j);
}

inline
pair<PatternDefs::Diff_t,bool>
CSRPattSymm::locate_dof(ID const i,ID const j) const {
  return locate_pattern(_d2i(i),_d2i(j));
}

inline
pair<UInt,UInt>
CSRPattSymm::giveMinMax() const
{
  vector<UInt> accumulate(max(_nrows,_ncols),0);// just to be sure!!!
  Diff_t col;
  Diff_t row;
  for (row=0;row<_nrows;++row){
    ++accumulate[row]; // diagonal
    for (Index_t j=_ia[row]+1;j< _ia[row+1];++j){
      // off-diagonal
      col=_i2o(_ja[j]);
      ++accumulate[row];
      ++accumulate[col];
    }
  }
  UInt curr_min=*min_element(accumulate.begin(),accumulate.end());
  UInt curr_max=*max_element(accumulate.begin(),accumulate.end());
  return make_pair(curr_min,curr_max);
}

template <typename Iter>
inline
UInt
CSRPattSymm::row(Diff_t const row, Iter  coldata, Iter position) const{
  // extracts a row (useful to implement A*b)
  // to make things faster no controls are made!!!
  // REMEMBER Only upper trangular part estracted since we have a
  // symmetric matrix. TO BE CHANGED IF WE HAVE THE WHOLE Matrix
  // (with symmetric pattern)
  Container::const_iterator start=_ja.begin()+_i2o(_ia[row]);
  Container::const_iterator end=_ja.begin()+_i2o(_ia[row+1]);
  copy(start,end,coldata);
  for (unsigned int i=0; i< _ia[row+1]-_ia[row];++i) position[i]=_i2o(_ia[row])+i;
  return _ia[row+1]-_ia[row];
}


////////////////////////////////////////////////////////////////////////
//
// M S R Pattern
//
////////////////////////////////////////////////////////////////////////

template<typename DOF1>
MSRPatt::MSRPatt(DOF1 const & dof1, UInt const nbcomp)
{
  bool built;
  built = buildPattern(dof1, nbcomp);
  ASSERT_PRE(built,"Error in MSR Pattern construction from DOF object");
};

// Miguel 12/2003
//
template<typename DOF, typename MESH>
MSRPatt::MSRPatt(const DOF& dof, const MESH& mesh, const UInt nbcomp)
{
  bool built;
  built = buildPattern(dof, mesh, nbcomp);
  ASSERT_PRE(built,"Error in MSR Pattern construction from DOF object");
};


template<typename DOF1>
bool MSRPatt::buildPattern(DOF1 const & dof1, UInt const nbcomp) {
  DynPattern dynpatt;
  bool built = setpatt(dof1, dynpatt, nbcomp);
  if (built) {
      _buildPattern<DOF1>(dof1, dynpatt, nbcomp);
      dynpatt.clear(); // not sure if this helps...
  }
  return built;
}


template<typename DOF1>
void MSRPatt::_buildPattern(DOF1 const & dof1, const DynPattern & dynpatt,
                            UInt const nbcomp) {
  Index_t ig,jg,cur;

  _bindx.resize(_nnz+1,0);
  _ybind.resize(_nnz-_nrows,0);

  // Count
  for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
    ig=d->first;
    jg=d->second;
    _bindx[ig]+= nbcomp;
    _bindx[jg]+= nbcomp;
  }

  //other components
  Diff_t offset= dof1.numTotalDof();
  UInt icomp, jcomp;

  //We add the missing diagonal term for jcomp != icomp
  for (icomp=0; icomp < nbcomp; icomp++)
    for (ig=0; ig < offset; ++ig)
      _bindx[ig+icomp*offset]+=nbcomp-1;

  Container::iterator start     = _bindx.begin();
  Container::iterator end       = _bindx.begin()+ offset;
  for (icomp=1; icomp < nbcomp; icomp++)
    {
      Container::iterator start_copy= _bindx.begin()+ icomp*offset;
      copy(start, end, start_copy);
    }

  // Count entries BY ACCUMULATING
  // shift right 1 position
  rotate(_bindx.begin(),_bindx.begin()+_nrows,_bindx.begin()+_nrows+1);
  _bindx[0]=_nrows+1;

  for (ig=1; ig <static_cast<Index_t>(_nrows)+1; ++ig)
    _bindx[ig]+=_bindx[ig-1];


  // for each component
  for (icomp=0; icomp< nbcomp; icomp++){
    for (jcomp=0; jcomp< nbcomp; jcomp++){

      //We add the missing diagonal term for jcomp != icomp
      if (jcomp != icomp)
	for (ig=0; ig < offset; ++ig){
	  _ybind[_bindx[ig+icomp*offset]-_nrows-1]=_bindx[ig+jcomp*offset];
	  _bindx[_bindx[ig+icomp*offset]++]=ig+jcomp*offset;
	}

      for (DynPattern::iterator d=dynpatt.begin(); d !=dynpatt.end(); ++d){
	ig=d->first+icomp*offset;
	jg=d->second+jcomp*offset;
	_ybind[_bindx[ig]-_nrows-1]=_bindx[jg];
	_ybind[_bindx[jg]-_nrows-1]=_bindx[ig];
	_bindx[_bindx[ig]++]=jg;
	_bindx[_bindx[jg]++]=ig;
      }
    }
  }

  // shift right 1 position
  rotate(_bindx.begin(),_bindx.begin()+_nrows,_bindx.begin()+_nrows+1);
  _bindx[0]=_nrows+1;



  for (ig=0; ig <static_cast<Index_t>(_nrows); ++ig){
    // We now sort the off diagonal entries
    jg=_bindx[ig];
    cur = _bindx[ig+1];
    sort(_bindx.begin()+jg,_bindx.begin()+cur);
  }

  // do we need it 'a la Fortran?'
  if (PatternOffset != 0) for (Container::iterator ip=_bindx.begin();ip!=_bindx.end();++ip)*ip+=PatternOffset;

  _diagfirst=true;// default for MSR
  _filled=true;
}



// Miguel 12/2003
//
//
template<typename DOF, typename MESH>
bool MSRPatt::buildPattern(const DOF&  dof, const MESH& mesh,
                           const UInt nbcomp) {
  DynPattern dynpatt;
  bool built = setpatt(dof, mesh, dynpatt, nbcomp);
  if (built) {
      _buildPattern<DOF>(dof, dynpatt, nbcomp);
      dynpatt.clear();
  }
  return built;
}

inline UInt MSRPatt::nbNeighbours(ID const d)const {
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");

  return _bindx[d]-_bindx[d-1] +1; // Diagonal IS INCLUDED
}

inline  ID MSRPatt::neighbour(ID const i,ID const d) const {
  ASSERT_BD(d>0 && d<=_nrows);
  ASSERT_BD(d>0);
  ASSERT_PRE(_filled, "Cannot access an empty pattern");
  if(i==1)return d;
  return _i2d(_bindx[_i2o(_bindx[d-1])+_d2o(i)-1]);
}

template <typename Iter>
inline
UInt
MSRPatt::row(Diff_t const row, Iter  coldata, Iter position) const{
  // extracts a row (useful to implement A*b)
  // to make things faster no controls are made!!!
  Container::const_iterator start=_bindx.begin()+_i2o(_bindx[row]);
  Container::const_iterator end=_bindx.begin()+_i2o(_bindx[row+1]);
  *coldata++=row;
  copy(start,end,coldata);
  position[0]=row;
  for (UInt i=0; i< _bindx[row+1]-_bindx[row];++i) position[i+1]=_i2o(_bindx[row])+i;
  return _bindx[row+1]-_bindx[row]+1;
}

inline
bool MSRPatt::isThere(Index_t i,Index_t j)const
{
  if (isEmpty()) return false;
  ASSERT_BD(i>=PatternOffset && i<static_cast<Index_t>(_nrows)+PatternOffset);
  ASSERT_BD(j>=PatternOffset && j<static_cast<Index_t>(_ncols)+PatternOffset);
  if (i==j) return true;
  Container::const_iterator start=_bindx.begin()+_row_off(i);
  Container::const_iterator finish=_bindx.begin()+_row_off(i+1);
  return binary_search(start, finish, j);
}


inline
pair<PatternDefs::Diff_t,bool>
MSRPatt:: locate_index(Index_t const i,Index_t const j) const{
  return locate_pattern(i,j);
}

inline
pair<PatternDefs::Diff_t,bool>
MSRPatt:: locate_dof(ID const i,ID const j) const{
  return locate_pattern(_d2i(i),_d2i(j));
}

inline
pair<UInt,UInt>
MSRPatt::giveMinMax()const{
  Container::const_iterator current;
  UInt curr_min=_ncols;
  UInt curr_max=0;
  for (current=_bindx.begin()+1;current<_bindx.begin()+(_nrows+1);++current)
    {
      UInt loc_col = *current - *(current-1);
      if (loc_col < curr_min) curr_min=loc_col;
      if (loc_col > curr_max) curr_max=loc_col;
    }
  return make_pair(curr_min+1,curr_max+1); // Adding the diagonal term
}


// Construction of a diagonal block matrix. Done by A. Gilardi.
// Alain (nov. 2002), update of ybind as well !
void diagblockMatrix(MSRPatt &ans, MSRPatt const &patt, UInt const nblock);



/*****************************************************************************************/
/*                              Mixed Pattern                */
/*****************************************************************************************/

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
MixedPattern<BROWS,BCOLS,PATTERN>::
MixedPattern()
{
  for (UInt i=0; i<BROWS; i++){
    for (UInt j=0; j<BCOLS; j++){
      _rowoff[i][j]=0;
      _coloff[i][j]=0;
      _linked[i][j]=false;
      _blocks[i][j]=0;
      _nrows[i][j]=0;
      _ncols[i][j]=0;
    }
  }
}


// Miguel 11/02: I want to construct a diagonal pattern
// in the initialisation list, so I need a constructor which makes the pattern
// construct and link to an external pattern (type = "diag") (Alain: type="full")
template <UInt BROWS, UInt BCOLS, typename PATTERN>
MixedPattern<BROWS,BCOLS,PATTERN>::MixedPattern(PATTERN & ex_patt, const string& type) {

  if (type == "full" ) {
    for (UInt i=0; i<BROWS; i++){
      for (UInt j=0; j<BCOLS; j++){
	_rowoff[i][j]=0;
	_coloff[i][j]=0;
	_linked[i][j]=false;
	_blocks[i][j]=0;
	_nrows[i][j]=0;
	_ncols[i][j]=0;
      }
    }
    for (UInt j=0; j<BCOLS; j++)
      for (UInt i=0; i<BROWS; i++)
	linkBlockToPattern(i,j,ex_patt);
  }
  else if (type == "diag" ) {
    for (UInt i=0; i<BROWS; i++){
      for (UInt j=0; j<BCOLS; j++){
	_rowoff[i][j]=0;
	_coloff[i][j]=0;
	_linked[i][j]=false;
	_blocks[i][j]=0;
	_nrows[i][j]=0;
	_ncols[i][j]=0;
      }
    }
    makeDiagPattern(ex_patt);
  }
  else
    ERROR_MSG("This type of construction is not allowed");
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
MixedPattern<BROWS,BCOLS,PATTERN>::
~MixedPattern()
{
  for (UInt i=0; i<BROWS; i++){
    for (UInt j=0; j<BCOLS; j++){
      if (_blocks[i][j] != 0 && !_linked[i][j]){
	delete _blocks[i][j]; // blocks are actually deleted only if not linked
	_blocks[i][j]=0;
      }
    }
  }
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void MixedPattern<BROWS,BCOLS,PATTERN>::
makeDiagPattern(PATTERN & ex_patt)
{
  ASSERT_PRE( BROWS == BCOLS, "attempt to make a diagonal pattern for non square matrix");
  for (UInt j=0; j<BCOLS; j++)
    for (UInt i=0; i<BROWS; i++)
      if (i == j)
	linkBlockToPattern(i,j,ex_patt);
      else
	buildZeroPattern(i,j, ex_patt.nRows(), ex_patt.nRows());
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
pair<UInt,UInt>
MixedPattern<BROWS,BCOLS,PATTERN>::nBlocks() const{
  return make_pair(BROWS,BCOLS);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nRows(Diff_t m, Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _nrows[m][n];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nCols(Diff_t m, Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _ncols[m][n];
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nNz(Diff_t m, Diff_t n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");

  return _blocks[m][n] ? _blocks[m][n]->nNz():0;
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::
nRows() const{
  // I assume that the blocks are consistent
  return _rowoff[BROWS-1][BCOLS-1] + _nrows[BROWS-1][BCOLS-1];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::
nCols() const{
  // I assume that the blocks are consistent
  return _coloff[BROWS-1][BCOLS-1] + _ncols[BROWS-1][BCOLS-1];
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nNz() const{
  UInt _n=0;
  for (UInt m=0; m<BROWS; m++){
    for (UInt n=0; n<BCOLS; n++){
      _n+=nNz(m,n);
    }
  }
  return _n;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nbNeighbours(Diff_t const m,  Diff_t const n,ID const d) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return (_blocks[m][n] ? _blocks[m][n]->nbNeighbours(d):0);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
ID
MixedPattern<BROWS,BCOLS,PATTERN>::neighbour(Diff_t const m, Diff_t const n,ID const i, ID const d) const
{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  ASSERT_PRE(_blocks[m][n] !=0 , "Cannot access an empty block") ;
  return _blocks[m][n]->neighbour(i,d);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
void
MixedPattern<BROWS,BCOLS,PATTERN>::neighbours(Diff_t const m, Diff_t const n,ID const d, Container & neighs) const
{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  ASSERT_PRE(_blocks[m][n] !=0 , "Cannot access an empty block") ;
  return _blocks[m][n]->neighbours(d, neighs);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
template <typename Iter>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::row(Diff_t const m, Diff_t const n,Diff_t const nrow, Iter  coldata, Iter position) const{
  // extracts a row (useful to implement A*b)
  // to make things faster no controls are made!!!
  // REMEMBER Only upper trangular part estracted since we have a
  // symmetric matrix. TO BE CHANGED IF WE HAVE THE WHOLE Matrix
  // (with symmetric pattern)
  return _blocks[m][n]!=0?_blocks[m][n]->row(nrow,coldata,position):0;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::nbNeighbours(ID const d_g) const{
  // This is MORE complicated
  //Locate line corresponding to dof d_g in block
  pair<UInt,UInt> _b=locateDofBlock(d_g,1); // Remember that dofs are ALWAYS numbered from 1
  // line local nomber
  UInt m=_b.first;
  UInt _d=d_g-_rowoff[m][0];
  UInt _n=0;
  // Sweep columns
  for (UInt n=0; n<BCOLS; n++){
    _n+=(_blocks[m][n] ? _blocks[m][n]->nbNeighbours(_d):0);
  }
  return _n;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
ID
MixedPattern<BROWS,BCOLS,PATTERN>::neighbour(ID const i_g, ID const d_g) const //ith neighbour of d (Global!)
{
  // This is MUCH MORE complicated
  //Locate line in blocks
  pair<UInt,UInt> _b=locateDofBlock(d_g,1);
  ASSERT_PRE(_b.first < BROWS, "Invalid dof entry d_g");
  // Get local numbering of dof
  UInt m=_b.first;
  UInt _d=d_g-_rowoff[m][0];

  // In which block does i_g fall?
  // I remember that all high level interfaces use the convention of having numbering from 1
  // So i_g and d_g have always numbering from 1

  UInt _first=0;
  UInt _last=0;
  UInt _i=0;
  UInt n=0;
  // Sweep columns
  for (n=0; n<BCOLS; n++){
    _last=_first+nbNeighbours(m,n,_d);
    if(i_g >_first && i_g <= _last){
      _i=i_g-_first;// local numbering
      break;
    }
    _first=_last;
  }
  ASSERT_PRE(_i !=0, " Invalid i_g="<< i_g<< ",  i= "<< _i ) ;

  return _blocks[m][n]->neighbour(_i,_d) + _coloff[m][n];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::neighbours(ID const d_g,Container & neighs) const //ith neighbour of d (Global!)
{
  // This is MUCH MORE complicated
  //Locate line in blocks
  pair<UInt,UInt> _b=locateBlock(d_g,0);
  ASSERT_PRE(_b<BROWS, "Invalid dof entry d_g");
  // Get local numbering of dof
  UInt m=b.first;
  UInt _d=d_g-_rowoff[m][0];

  // In which block does i_g fall?
  // I remember that all high level interfaces use the convention of having numbering from 1
  // So i_g and d_g have always numbering from 1
  //neighs.clear();
  neighs.resize(nbNeighbours(d_g));
  if (neighs.empty()) return;
  Container piece;
  // Sweep columns
  for (UInt n=0; n<BCOLS; n++){
    if (_blocks[m][n] != 0 ){
      _blocks[m][n]->neighbours(_d,piece);
      for(Container::iterator ip=piece.begin();ip=!piece.end();++ip) *ip+=_coloff[m][n];
      neighs.insert(neighs.end(),piece.begin(),piece.end());
    }
  }
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
template <typename Iter>
inline
UInt
MixedPattern<BROWS,BCOLS,PATTERN>::
row(Diff_t const nrow, Iter  coldata, Iter position) const{

  Diff_t m=locateDofBlock(nrow+1,1).first;
  UInt _count=0;
  Diff_t lrow=(nrow-_rowoff[m][0]);

  UInt many=0;
  for (UInt n=0; n<BCOLS; n++){
    coldata+=many;
    position+=many;
    if(_blocks[m][n]!=0){
      many=(_blocks[m][n]->row(lrow,coldata,position));
      for (Iter itr=coldata; itr< coldata+many; ++itr)
	*itr += _coloff[m][n];
      _count+=many;
    }
  }

  return _count;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
PATTERN *
MixedPattern<BROWS,BCOLS,PATTERN>::block_ptr(Diff_t const m, Diff_t const n){
  // VERY DANGEROUS
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _blocks[m][n]; // Beware of null pointers!!
}
// const qualifyer version
// Alain 05/02.
template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
const PATTERN *
MixedPattern<BROWS,BCOLS,PATTERN>::
block_ptr(Diff_t const m, Diff_t const n) const{
  // VERY DANGEROUS
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _blocks[m][n]; // Beware of null pointers!!
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
pair<UInt,UInt>
MixedPattern<BROWS,BCOLS,PATTERN>::blockOffset(UInt const m, UInt const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return make_pair(_rowoff[m][n],_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<UInt,UInt>
MixedPattern<BROWS,BCOLS,PATTERN>::
locateElBlock(Index_t const i_g, Index_t const j_g) const
{
  UInt m=0;
  UInt n=0;
  UInt _itest=_i2o(i_g);
  UInt _jtest=_i2o(j_g);

  while(m<BROWS &&(_itest <_rowoff[m][0] || _itest >=_rowoff[m][0]+ _nrows[m][n])) ++m;
  while(n<BCOLS &&(_jtest <_coloff[0][n] || _jtest >=_coloff[0][n]+ _ncols[m][n])) ++n;

  return make_pair(m,n);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<PatternDefs::Diff_t,PatternDefs::Diff_t>
MixedPattern<BROWS,BCOLS,PATTERN>::
locateDofBlock(ID const i_g, ID const j_g) const
{
  UInt m=0;
  UInt n=0;
  UInt _itest=_d2o(i_g);
  UInt _jtest=_d2o(j_g);

  while(m<BROWS &&(_itest <_rowoff[m][0] || _itest >=_rowoff[m][0]+ _nrows[m][n])) ++m;
  while(n<BCOLS &&(_jtest <_coloff[0][n] || _jtest >=_coloff[0][n]+ _ncols[m][n])) ++n;
  return make_pair(m,n);
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
pair<ID,ID>
MixedPattern<BROWS,BCOLS,PATTERN>::localNumber(Diff_t const m, Diff_t const n,ID const i_g, ID const j_g) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");

  return make_pair(i_g-_rowoff[m][n],j_g-_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
pair<ID,ID>
MixedPattern<BROWS,BCOLS,PATTERN>::globalNumber(Diff_t const m, Diff_t const n,ID const i, ID const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return make_pair(i+_rowoff[m][n],j+_coloff[m][n]);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<PatternDefs::Diff_t,bool>
MixedPattern<BROWS,BCOLS,PATTERN>::locateIndex(Diff_t const m, Diff_t const n,Index_t const i, Index_t const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  if (_blocks[m][n]==0)return make_pair(0,false);
  return _blocks[m][n]->locateIndex(i,j);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
pair<PatternDefs::Diff_t,bool>
MixedPattern<BROWS,BCOLS,PATTERN>::locateDof(Diff_t const m, Diff_t const n,ID const i, ID const j) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  if (_blocks[m][n]==0)return make_pair(0,false);
  return _blocks[m][n]->locateDof(i,j);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
buildZeroPattern(Diff_t const m, Diff_t const n,UInt const nrows, UInt const ncols){
  _nrows[m][n]=nrows;
  _ncols[m][n]=ncols;
  resetOffset(m,n);
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
linkBlockToPattern(UInt const m, UInt const n,PATTERN & pattern){
  if (m>0)
    ASSERT_PRE( isSet(m-1,n), "pattern must be set starting from block row=0");

  if (n>0)
    ASSERT_PRE( isSet(m,n-1), "pattern must be set starting from block col=0");

  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  if(_blocks[m][n]!=0 && !_linked[m][n]) delete _blocks[m][n];

  _blocks[m][n]=&pattern;
  _nrows[m][n]=pattern.nRows();
  _ncols[m][n]=pattern.nCols();
  _linked[m][n]=true;
  resetOffset(m,n);
}

template <UInt BROWS, UInt BCOLS, typename PATTERN>
template <typename DOF1, typename DOF2>
PATTERN * MixedPattern<BROWS,BCOLS,PATTERN>::buildBlock(Diff_t m, Diff_t n, DOF1 dof1, DOF2 dof2){
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");

  if(_blocks[m][n]!=0 && !_linked[m][n]) delete _blocks[m][n];

  _blocks[m][n]=new PATTERN(dof1,dof2);
  _nrows[m][n]=_blocks[m][n]->nRows();
  _ncols[m][n]=_blocks[m][n]->nCols();
  _linked[m][n]=false;
  resetOffset(m,n);
  return _blocks[m][n];
}

template <UInt BROWS, UInt BCOLS, typename PATTERN>
template <typename DOF1>
PATTERN * MixedPattern<BROWS,BCOLS,PATTERN>::buildBlock(Diff_t m, Diff_t n, DOF1 dof1){
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");

  if(_blocks[m][n]!=0 && !_linked[m][n]) delete _blocks[m][n];

  _blocks[m][n]=new PATTERN(dof1);
  _nrows[m][n]=_blocks[m][n]->nRows();
  _ncols[m][n]=_blocks[m][n]->nCols();
  _linked[m][n]=false;
  resetOffset(m,n);
  return _blocks[m][n];
}



template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
linkBlocks(Diff_t const m_from, Diff_t const n_from, Diff_t const m_to, Diff_t const n_to){
  // ASSERT_PRE if it is ok!
  ASSERT_PRE( m_from<BROWS && m_from>=0, " Invalid block row address");
  ASSERT_PRE( n_from<BCOLS && n_from>=0, " Invalid block column address");
  ASSERT_PRE( m_to<BROWS && m_to>=0, " Invalid block row address");
  ASSERT_PRE( n_to<BCOLS && n_to>=0, " Invalid block column address");

  _blocks[m_to][n_to]=_blocks[m_from][n_from];
  _nrows[m_to][n_to]=_nrows[m_from][n_from];
  _ncols[m_to][n_to]=_ncols[m_from][n_from];
  _linked[m_to][n_to]=true;
  resetOffset(m,n);
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::
deleteBlock(Diff_t const m, Diff_t const n){
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  if(_blocks[m][n]!=0&& ! linked[m][n]) delete _blocks[m][n];
  _blocks[m][n]=0;
  _nrows[m][n]=0;
  _ncols[m][n]=0;
  _linked[m][n]=false;
  resetOffset(m,n);
}



template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
bool
MixedPattern<BROWS,BCOLS,PATTERN>::isZero(Diff_t const m, Diff_t const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _blocks[m][n]=0 && _nrows[m][n]+_ncols[m][n] != 0;
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
inline
bool
MixedPattern<BROWS,BCOLS,PATTERN>::isSet(Diff_t const m, Diff_t const n) const{
  ASSERT_PRE( m<BROWS && m>=0, " Invalid block row address");
  ASSERT_PRE( n<BCOLS && n>=0, " Invalid block column address");
  return _blocks[m][n]!=0 || _nrows[m][n]!=0 || _ncols[m][n]!=0 || _linked[m][n];
}

template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::resetOffset(Diff_t const m, Diff_t const n) // Puts offset right!
{
  UInt _m=m;
  UInt _n=n;

  while( _m<BROWS && _n <BCOLS){

    for (UInt i=_m+1; i<BROWS; i++){
      _rowoff[i][_n]=_rowoff[i-1][_n]+_nrows[i-1][_n];
    }

    for (UInt i=_n+1; i<BCOLS; i++){
      _coloff[_m][i]=_coloff[_m][i-1]+_ncols[_m][i-1];
    }
    _m++;
    _n++;
  }
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
bool
MixedPattern<BROWS,BCOLS,PATTERN>::check(bool verbose, ostream & c) const{
  bool ok=true;

  for (Diff_t m=0; m<BROWS; m++){
    for (Diff_t n=0; n<BCOLS; n++){
      ok=ok && _rowoff[m][n]==_rowoff[m][n-1] &&_coloff[m][n]==_coloff[m-1][n];
    }
  }
  for (Diff_t m=0; m<BROWS; m++){
    for (Diff_t n=0; n<BCOLS; n++){
      if(_blocks[m][n]==0&&_nrows[m][n]==0&&_ncols[m][n]==0){
	ok=false;
	c<< "Block " <<"("<<m<<", "<<n<<") unset"<<endl;
      }
    }
  }
  for (Diff_t m=0; m<BROWS; m++){
    for (Diff_t n=0; n<BCOLS; n++){
      ok=ok && _rowoff[m][n]-_rowoff[m-1][n] +_nrows[m-1][n]==0 &&
	_coloff[m][n]-_coloff[m][n-1]-_ncols[m][n-1]==0;
    }
  }

  if (verbose && ! ok){
    c<< " Error checking mixed pattern"<<endl;
    for (Diff_t m=0; m<BROWS; m++){
      for (Diff_t n=0; n<BCOLS; n++){
	if(_rowoff[m][n]!=_rowoff[m][n-1] || _coloff[m][n]!=_coloff[m-1][n]){
	  c<<"Inconsistent Offsets of Blocks:"<<endl;
	  c<<"("<<m<<", "<<n<<") Row offset= "<<_rowoff[m][n]
	   <<"("<<m<<", "<<n-1<<") Row offset= "<<_rowoff[m][n-1]<<endl;
	  c<<"("<<m<<", "<<n<<") Col offset= "<<_coloff[m][n]
	   <<"("<<m-1<<", "<<n<<") Col offset= "<<_coloff[m-1][n]<<endl;
	}
      }
    }
    for (Diff_t m=0; m<BROWS; m++){
      for (Diff_t n=0; n<BCOLS; n++){
	if(_rowoff[m][n]-_rowoff[m-1][n] +_nrows[m-1][n]!=0 ||
	   _coloff[m][n]-_coloff[m][n-1]-_ncols[m][n-1]!=0){
	  c<<"Inconsistent Number of cols/rows:"<<endl;
	  c<<"("<<m<<", "<<n<<") Row offset= "<<_rowoff[m][n];
	  c<<"("<<m-1<<", "<<n<<") Row offset= "<<_rowoff[m-1][n]<<endl;
	  c<<"("<<m-1<<", "<<n<<") N Rows= "<<_nrows[m-1][n]<<endl;
	  c<<"("<<m<<", "<<n<<") Col offset= "<<_coloff[m][n];
	  c<<"("<<m<<", "<<n-1<<") Col offset= "<<_coloff[m][n-1]<<endl;
	  c<<"("<<m<<", "<<n-1<<") N Cols= "<<_ncols[m][n-1]<<endl;
	}
      }
    }
  }
  return ok;
}


template
<UInt BROWS, UInt BCOLS, typename PATTERN>
void
MixedPattern<BROWS,BCOLS,PATTERN>::showMe(bool verbose, ostream & c) const{
  c<<endl;
  c<<" ******************* MIXED PATTERN *******************"<<endl;
  c<<nBlocks().first<<"X"<<nBlocks().second<<" Blocks"<<endl;
  c<<"Block    Row Off   Cols OFF    Nrows     Ncols     Set      Linked  "<<endl;
  for (Diff_t m=0; m<BROWS; m++){
    for (Diff_t n=0; n<BCOLS; n++){
      c<<"("<<m<<", "<<n<<") "<<_rowoff[m][n]<<" "<<_coloff[m][n]<<" "<<_nrows[m][n]<<" "<<_ncols[m][n]<<" ";
      c<<isSet(m,n)<<" "<<_linked[m][n]<<endl;
    }
  }
  if(verbose){
    c<<endl<< "Mixed pattern: single block info"<<endl;

    for (Diff_t m=0; m<BROWS; m++){
      for (Diff_t n=0; n<BCOLS; n++){
	c<<"("<<m<<", "<<n<<") ->  " <<endl;
	if (_blocks[m][n] != 0 ){
	  _blocks[m][n]->showMe(false,c);
	}
	else{
	  if(_ncols[m][n] ==0 && _nrows[m][n] ==0)
	    c<<" UNSET"<<endl;
	  else
	    c<< "Zero Pattern "<<  _nrows[m][n]<<" X "<< _ncols[m][n]<<endl;
	}
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////
// This work both for Mixed and simple pattern
//
template <typename FE1, typename FE2=FE1>
class Pattern
{
public:
  // For the laziest
  static
  UInt patternFirst(UInt m, UInt n, UInt i);
  static
  UInt patternSecond(UInt m, UInt n, UInt i);
  static
  UInt nbPattern(UInt m, UInt n);
  static
  UInt patternFirst11(const UInt i)
  {return FE1::patternFirst(i);}
  static
  UInt patternFirst22(const UInt i)
  {return FE2::patternFirst(i);}
  static
  UInt patternSecond11(const UInt i)
  {return FE1::patternSecond(i);}
  static
  UInt patternSecond22(const UInt i)
  {return FE2::patternSecond(i);}

  static
  UInt patternFirst12(const UInt i)
  {return i/FE2::nbNode;}

  static
  UInt patternSecond12(const UInt i)
  {return i % FE2::nbNode;}

  static
  UInt patternFirst21(const UInt i)
  {return i/FE1::nbNode;}
  static
  UInt patternSecond21(const UInt i)
  {return i % FE1::nbNode;}

  static UInt
  nbPattern11(){
    return FE1::nbPattern;
  }
  static UInt
  nbPattern22(){
    return FE2::nbPattern;
  }
  static UInt
  nbPattern12(){
    return FE1::nbNode*FE2::nbNode;
  }
  static UInt
  nbPattern21(){
    return FE1::nbNode*FE2::nbNode;
  }
};

/*****************************************************************************************/
/*                              Mixed Local Pattern                */
/*****************************************************************************************/

template<typename FE1, typename FE2>
UInt
MixedLocalPattern<FE1,FE2>::nRows() const
{
  return FE1::nbNode;
}


template<typename FE1, typename FE2>
UInt
MixedLocalPattern<FE1,FE2>::nCols() const
{
  return FE2::nbNode;
}

template<typename FE1, typename FE2>
UInt
MixedLocalPattern<FE1,FE2>::nbPattern() const
{
  return FE2::nbNode*FE1::nbNode;
}


template<typename FE1, typename FE2>
UInt
MixedLocalPattern<FE1,FE2>::patternFirst(UInt const i) const
{
  // Along Rows (changes slowly)
  return i/FE2::nbNode;
}

template<typename FE1, typename FE2>
UInt
MixedLocalPattern<FE1,FE2>::patternSecond(UInt const i) const
{
  // Along columns (changes fast)
  return i % FE2::nbNode;
}
}

#endif

