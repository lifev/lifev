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
/*!
  \file bcCond.cc
  \brief Implementations for bcCond.h
  \version 1.0
  \author M.A. Fernandez
  \date 07/2002

*/


#include "bcCond.hpp"

namespace LifeV
{
// ============ BCFunction_Base ================

BCFunction_Base::BCFunction_Base(const BCFunction_Base& bcf){
  _g=bcf._g;
}

//! Constructor for BCFuncion_Base from a user defined function
BCFunction_Base::BCFunction_Base(Function g){
  _g=g;
}

//! set the function after having built it.
void BCFunction_Base::setFunction(Function g){
  _g = g;
}

//! Overloading function operator by calling attribut _g
Real BCFunction_Base::operator()(const Real& t, const Real& x, const Real& y,
			    const Real& z, const ID& icomp) const {
  return _g(t,x,y,z,icomp);
}

BCFunction_Mixte::BCFunction_Mixte(const BCFunction_Mixte& bcf):BCFunction_Base(bcf._g) {
  _coef = bcf._coef;
}

BCFunction_Mixte::BCFunction_Mixte(Function g, Function coef):BCFunction_Base(g) {
  _coef = coef;
}

//! set the functions after having built it.
void BCFunction_Mixte::setFunctions_Mixte(Function g, Function coef){
  _g = g;
  _coef = coef;
}

Real BCFunction_Mixte::coef(const Real& t, const Real& x, const Real& y,
		       const Real& z, const ID& icomp) const {
  return _coef(t,x,y,z,icomp);
}


// ============ BC_Base ================
//! Constructor for BC
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		 BCFunction_Base& bcf, const vector<ID>& comp){
  _name = name;
  _flag = flag;
  _type = type;
  if (mode != Component)
    ERROR_MSG("You should use a more specific constructor for this mode");
  _mode = mode;

  switch( type ) {
  case Essential:
  case Natural:
    _bcf = new BCFunction_Base( bcf );
    break;
  case Mixte:
    // We have to hold a pointer to a mixte functor
    _bcf = new BCFunction_Mixte( static_cast<BCFunction_Mixte&>(bcf) );
    break;
  }
  _dataVector = 0;

  UInt nComp = comp.size();
  _comp.reserve(nComp);
  _comp.insert(_comp.end(),comp.begin(),comp.end());

  _finalised = 0;
}

//! Constructor for BC without components for Scalar, Tangential or Normal  mode problems
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		 BCFunction_Base& bcf) {
  _name = name;
  _flag = flag;
  _type = type;

  UInt nComp;
  switch( _mode=mode ){
  case Scalar:
    nComp = 1;
    _comp.reserve(nComp);
    _comp.push_back(1);
    switch( type ) {
    case Essential:
    case Natural:
      _bcf = new BCFunction_Base( bcf );
      break;
    case Mixte:
      // We have to hold a pointer to a mixte functor
      _bcf = new BCFunction_Mixte( static_cast<BCFunction_Mixte&>(bcf) );
      break;
    }
    break;
  case Tangential:
    nComp = nDimensions - 1;
    _comp.reserve(nComp);
    for (ID i=1; i<=nComp; ++i)
      _comp.push_back(i);
    _bcf = new BCFunction_Base( bcf );
    break;
  case Normal:
    nComp=1;
    _comp.reserve(nComp);
    _comp.push_back(nDimensions);
    _bcf = new BCFunction_Base( bcf );
    break;
  default:
    ERROR_MSG("You should use a more specific constructor for this mode");
  }
  _dataVector = 0;
  _finalised = 0;
}


//! Constructor for BC without list of components for Full mode problems
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type,
	const BCMode& mode, BCFunction_Base& bcf, const UInt& nComp){
  _name = name;
  _flag = flag;
  _type = type;

  if (mode != Full)
    ERROR_MSG("You should use a more specific constructor for this mode");
  _mode = mode;


  switch( type ) {
  case Essential:
  case Natural:
    _bcf = new BCFunction_Base( bcf );
    break;
  case Mixte:
    // We have to hold a pointer to a mixte functor
    _bcf = new BCFunction_Mixte( static_cast<BCFunction_Mixte&>(bcf) );
    break;
  }

  _comp.reserve(nComp);
  for (ID i=1; i<=nComp; ++i)
    _comp.push_back(i);

  _dataVector = 0;
  _finalised = 0;
}



//! Constructor for BC with data vector
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		 BCVector_Base& bcv, const vector<ID>& comp){
  _name = name;
  _flag = flag;
  _type = type;
  if (mode != Component)
    ERROR_MSG("You should use a more specific constructor for this mode");
  _mode = mode;

  _dataVector = 1;
  _bcv = &bcv;


  UInt nComp = comp.size();
  _comp.reserve(nComp);
  _comp.insert(_comp.end(),comp.begin(),comp.end());

  _finalised = 0;
}

//! Constructor for BC with data vector, without components for Scalar, Tangential or Normal  mode problems
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		 BCVector_Base& bcv) {
  _name = name;
  _flag = flag;
  _type = type;

  UInt nComp;
  switch( _mode=mode ){
  case Scalar:
    nComp = 1;
    _comp.reserve(nComp);
    _comp.push_back(1);

    _dataVector = 1;
    _bcv = &bcv;

    break;
  case Tangential:
    nComp = nDimensions - 1;
    _comp.reserve(nComp);
    for (ID i=1; i<=nComp; ++i)
      _comp.push_back(i);

    _dataVector = 1;
    _bcv = &bcv;

    break;
  case Normal:
    nComp=1;
    _comp.reserve(nComp);
    _comp.push_back(nDimensions);

    _dataVector = 1;
    _bcv = &bcv;

    break;
  default:
    ERROR_MSG("You should use a more specific constructor for this mode");
  }
  _finalised = 0;
}


//! Constructor for BC with data vector, without list of components for Full mode problems
BC_Base::BC_Base(const string& name, const EntityFlag& flag, const BCType& type,
	const BCMode& mode, BCVector_Base& bcv, const UInt& nComp){
  _name = name;
  _flag = flag;
  _type = type;

  if (mode != Full)
    ERROR_MSG("You should use a more specific constructor for this mode");
  _mode = mode;

  _dataVector = 1;
  _bcv = &bcv;

  _comp.reserve(nComp);
  for (ID i=1; i<=nComp; ++i)
    _comp.push_back(i);

  _finalised = 0;
}


//! Destructor (we have a vector of pointers to ID's and Functors)
BC_Base::~BC_Base() {

  // Important!!: The set member list0 is always empty at this point, it is just
  // an auxiliary container used at the moment of the boundary update (see BC_Handler::bdUpdate)

  switch( _type ) {
  case Essential:
    // Deleting identifiers
    for( IDIterator i=_idList.begin();i!=_idList.end();++i)
      delete static_cast<Identifier_Essential*>(*i);

    break;
  case Natural:
  case Mixte:
    for( IDIterator i=_idList.begin();i!=_idList.end();++i)
      // Deleting identifiers
      delete static_cast<Identifier_Natural*>(*i);
    break;
  }

  if (_dataVector) {
    _dataVector = 0;
  }
  else {
    // Deleting user defined functors
    switch( _type ) {
    case Natural: case Essential:
      delete _bcf;
      break;
    case Mixte:
      delete static_cast<BCFunction_Mixte*>(_bcf);
      break;
    }
  }
}



//! Assignment operator for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
BC_Base & BC_Base::operator=(const BC_Base& BCb) {

  _name = BCb._name;
  _flag = BCb._flag;
  _type = BCb._type;
  _mode = BCb._mode;
  _finalised = BCb._finalised;
  _dataVector = BCb._dataVector;

  switch ( _type ) {
  case Essential:
  case Natural:
    if ( _dataVector ) {
      _bcv = BCb._bcv;
    }
    else {
      _bcf = new BCFunction_Base( *(BCb._bcf) );
    }
    break;
  case Mixte:
    if ( _dataVector ) {
      _bcv = BCb._bcv;
    }
    else {
      _bcf = new BCFunction_Mixte( *( static_cast<BCFunction_Mixte*>(BCb._bcf) ));
    }
    break;
  }

  _comp = BCb._comp;

  // Important!!: The set member list0 is always empty at this point, it is just
  // an auxiliary container used at the moment of the boundary update (see BC_Handler::bdUpdate)

  // The list of ID's must be empty
  if ( !_idList.empty() || !BCb._idList.empty() )
    ERROR_MSG("The BC assigment operator does not work with lists of identifiers which are not empty");

  return *this;
}

//! Copy constructor for BC (we have a vector of pointers to ID's and a pointer to user defined functions)
BC_Base::BC_Base(const BC_Base& BCb) {
  _name = BCb._name;
  _flag = BCb._flag;
  _type = BCb._type;
  _mode = BCb._mode;
  _finalised = BCb._finalised;
  _dataVector = BCb._dataVector;

  switch ( _type ) {
  case Essential:
  case Natural:
    if ( _dataVector ) {
      _bcv = BCb._bcv;
    }
    else {
      _bcf = new BCFunction_Base( *(BCb._bcf) );
    }
    break;
  case Mixte:
    if ( _dataVector ) {
      _bcv = BCb._bcv;
    }
    else {
      _bcf = new BCFunction_Mixte( *( static_cast<BCFunction_Mixte*>(BCb._bcf) ) );
    }
    break;
  }
  _comp = BCb._comp;

  // Important!!: The set member list0 is always empty at this point, it is just
  // an auxiliary container used at the moment of the boundary update (see BC_Handler::bdUpdate)

  // The list of ID's must be empty
  if ( !_idList.empty() || !BCb._idList.empty() )
    ERROR_MSG("The BC copy constructor does not work whith list of identifiers which are not empty");
}

//! Returns the BC name
string BC_Base::name() const { return _name; }

//! Returns the BC associated flag
EntityFlag BC_Base::flag() const { return _flag; }

//! Returns the BC type
BCType BC_Base::type() const { return _type; }

//! Returns the BC mode
BCMode BC_Base::mode() const { return _mode; }

//! Returns the number of components involved in this boundary condition
UInt BC_Base::numberOfComponents() const { return _comp.size(); }

//! eturns the global i-th component involved in the boundary condition
ID BC_Base::component(const ID i) const {
  ASSERT_BD(i >= 1 && i <= _comp.size());
  return _comp[i-1];
}

//! Returns wether the list is finalised and the vector of ID's is then accessible
bool BC_Base::finalised() const { return _finalised; }

//! Overloading function operator by calling the (*_bcf)() user specified function
Real BC_Base::operator()(const Real& t, const Real& x, const Real& y,
		       const Real& z, const ID& i) const {
  return _bcf->operator()(t,x,y,z,i);
}



//! Returns a pointer  to the user defined STL functor
BCFunction_Base*  BC_Base::pointerToFunctor() const {
  return _bcf;
}

//! True is a data vector has been provided
bool BC_Base::dataVector() const {
  return _dataVector;
}



Real BC_Base::operator()(const ID& iDof, const ID& iComp) const {
  if ( _dataVector )
    return (*_bcv)(iDof,iComp);
  else
    ERROR_MSG("A data vector must be specified before calling this method");
}


//! Returns the value of the mixte coefficient (in BC Vector)
Real BC_Base::MixteCoef() const{
  if ( _dataVector )
    return (*_bcv).MixteCoef();
  else
    ERROR_MSG("A data vector must be specified before calling this method");
}

//! Returns a pointer  to the i-th elements in the (finalised) list
//! (counting from 1 ' a la FORTRAN')
Identifier_Base* BC_Base::operator()(const ID& i) const {
  return this->operator[](i-1);
}

//! Returns a pointer to the i-th elements in the (finalised) list
//! (counting from 0 ' a la C')
Identifier_Base* BC_Base::operator[](const Index_t& i) const {
  ASSERT_PRE(_finalised, "BC List should be finalised before being accessed");
  ASSERT_BD(i >= 0 && i < _idList.size());
  return _idList[i];
}


//! Add a new identifier to the preliminary list of Identifiers
void BC_Base::addIdentifier(Identifier_Base* iden) {
  list0.insert(iden);
}


//! Transfer between the list and vector containers of ID's
void BC_Base::finalise()
{
  if (! list0.empty()){
    _idList.clear();
    _idList.reserve(list0.size());
    copy(list0.begin(),list0.end(),inserter(_idList,_idList.end()));
    list0.clear();
  }
  _finalised=true;
}


//! Returns the liste size
UInt BC_Base::list_size() const { return _idList.size();}

//! Output
ostream& BC_Base::showMe(bool verbose, ostream & out) const {
  out << "********************************"<<endl;
  out << "BC Name: " << _name << endl;
  out << "Flag: " << _flag << endl;
  out << "Type: " << _type << endl;
  out << "Mode: " << _mode << endl;
  out << "Number of components: " <<_comp.size() <<endl;
  out << "List of components: ";
  for (Index_t i=0; i<_comp.size(); ++i)
    out << _comp[i] <<" ";
  out << endl;
  out << "Number of stored ID's: " << _idList.size() << endl;

  if (verbose && _finalised){
    unsigned int count(0),lines(10);
    out << "IDs in list";
    for (vector<Identifier_Base*>::const_iterator i= _idList.begin(); i != _idList.end(); i++){
      if (count++ % lines ==0){
	out << endl;
      }
      out<< (*i)->id() <<"  ";
    }
    if (count % lines !=0){ out<<endl;}
    if ( dataVector() ){
      _bcv->showMe(verbose, out);
    }
  }

  out << "********************************" << endl;
  return out;
}



// ============ BC_Handler ================

//! Constructor doing nothing (the user must call setNumber(..))
BC_Handler::BC_Handler():_nbc(0),_bdUpdateDone(0),_fullEssential(0){}


//! Set the number of BC to be stored
void BC_Handler::setNumber(const ID& nbc)
{
  _bcList.reserve(_nbc=nbc);
}

//! Constructor taking the number of BC to be stored
BC_Handler::BC_Handler(const ID& nbc,const bool& fullEssential){
  _bcList.reserve(_nbc=nbc);
  _bdUpdateDone=0;
  _fullEssential=fullEssential;
}

//! Constructor taking the number of BC to be stored
BC_Handler::BC_Handler(const ID& nbc){
  _bcList.reserve(_nbc=nbc);
  _bdUpdateDone=0;
  _fullEssential=0;
}

//! How many BC stored?
Index_t BC_Handler::size() const {
  return _bcList.size();
}

//! Is there no BC into the list?
bool BC_Handler::empty() const {
  return _bcList.empty();
}

//! Adding new BC to the list with user defined functions
void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		       BCFunction_Base& bcf, const vector<ID>& comp) {
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcf, comp) );


  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}
void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		       BCFunction_Base& bcf){
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcf) );


  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}

void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		       BCFunction_Base& bcf, const UInt& nComp){
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcf, nComp) );


  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}


//! Adding new BC to the list with data vectors
void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		       BCVector_Base& bcv, const vector<ID>& comp) {
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcv, comp) );

  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}
void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		        BCVector_Base& bcv){
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcv) );

  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}

void BC_Handler::addBC(const string& name, const EntityFlag& flag, const BCType& type, const BCMode& mode,
		       BCVector_Base& bcv, const UInt& nComp){
  if ( _nbc == _bcList.size() )
    ERROR_MSG("You cannot add another BC, the list of BC in the Handler is full");
  else if (_nbc == 0)
    ERROR_MSG("Before adding BCs, you must specify the number of BC to be stored");

  // Adding BC
  _bcList.push_back( BC_Base(name, flag, type, mode, bcv, nComp) );

  if ( _nbc == _bcList.size() )
    // Sorting list of BC. Essential BC must be treated at the end !!!!
    sort(_bcList.begin(), _bcList.end());
}


// returns true if the bdUpdate has been done before
bool BC_Handler::bdUpdateDone() const {
  return _bdUpdateDone;
}

//! returns true if all the stored BC are of Essential type
bool BC_Handler::fullEssential() const {
  return _fullEssential;
}

//! Extracting BC from the list
BC_Base& BC_Handler::operator[](const Index_t& i){
  ASSERT_PRE(_nbc == _bcList.size(), "Some BC have not been added to the list");
  return _bcList[i];
}


const BC_Base& BC_Handler::operator[](const Index_t& i) const {
  ASSERT_PRE(_nbc == _bcList.size(), "Some BC have not added to the list");
  return _bcList[i];
}


//! Ouput
ostream & BC_Handler::showMe(bool verbose, ostream & out) const
{
  out<< " Boundary Conditions Handler ====>"<<endl;
  if ( _nbc != _bcList.size() )
    out<<" Some BC have not been added to the list\n";
  out<< " Number of BC stored "<<size()<<endl;;
  out<< " List => "<<endl;
  for(UInt i=0; i<_nbc; ++i){
    _bcList[i].showMe(verbose,out);
  }
  out<< " <===========================>"<<endl;
return out;
}
}
